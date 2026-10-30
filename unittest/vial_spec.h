#pragma once

// The VIAL reference as a test suite.
//
// docs/vial_reference.md is the specification of VIAL, and its examples are the tests: every line of a ```python
// block may end in a comment that says what compiling or running it must give. This header reads the reference into
// blocks and expectations, by the grammar in its section "How to read this reference":
//
//   # -> float[433][3] Å, frame: commentary       a result: a type or a value, then stages, then commentary
//   # error E0301: commentary                     compiling fails with this code on this line
//   # warning W0601, W0405 -> float Å | missing   it compiles with these warnings, and this result
//   # warning W0405; R0101 in frames where …      warnings, then a note (here a runtime code)
//   # on T, # on P                                the structure the block runs on, from here on
//   # on T, with the names above                  continues the block before it, as one script
//   # on T, as filters in …                       each line is a filter, evaluated within a group
//
// A comment that starts with none of these is commentary and is not checked. Each phase of the implementation checks
// what it can: the parser that every block parses and gives the syntax codes it expects, the checker the types and
// compile-time codes, the interpreter the values and runtime codes. Everything the reference says is recorded here,
// so the harness can report what is still unchecked.

#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_array.h>
#include <core/md_allocator.h>

#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

#define VIAL_SPEC_MAX_CODES 4

typedef struct vial_code_t { char s[6]; } vial_code_t;   // "E0301", zero-terminated

typedef enum vial_expect_kind_t {
    VIAL_EXPECT_NONE = 0,   // commentary, or no comment
    VIAL_EXPECT_RESULT,     // -> ...
    VIAL_EXPECT_ERROR,      // error E....
    VIAL_EXPECT_WARNING,    // warning W...., with or without a result
} vial_expect_kind_t;

typedef struct vial_expect_t {
    vial_expect_kind_t kind;
    int     line;                               // 1-based line in the reference
    int     line_in_block;                      // 1-based line within the block's source
    char    structure;                          // 'T', 'P', or 0, as named by the last "on ..." comment above
    int     num_codes;
    char    codes[VIAL_SPEC_MAX_CODES][6];      // "E0301", "W0405", in the order written
    str_t   result;                             // the type or value after "->", trimmed; empty when there is none
    str_t   note;                               // what follows the result: stages, commentary
    str_t   comment;                            // the whole comment, for messages
} vial_expect_t;

typedef enum vial_block_mode_t {
    VIAL_BLOCK_SCRIPT = 0,  // a script
    VIAL_BLOCK_FILTERS,     // each line is a filter
} vial_block_mode_t;

typedef struct vial_block_t {
    int     line;               // 1-based line of the block's first source line in the reference
    str_t   section;            // the nearest heading above the block
    str_t   src;                // the block's source, a view into the reference
    char    structure;          // 'T', 'P', or 0 for none or another (such as a run with an energy file)
    bool    continues;          // "with the names above": one script with the block before it
    vial_block_mode_t mode;
    vial_expect_t* expects;     // md_array, one entry per line that has an expectation
} vial_block_t;

typedef struct vial_spec_t {
    str_t text;                 // the whole reference
    vial_block_t* blocks;       // md_array
    vial_code_t* diag_codes;    // md_array: the codes listed in the reference's diagnostics table
    int   num_unparsed;         // expectation comments that do not follow the grammar
} vial_spec_t;

static inline bool vial_spec__is_code(str_t s) {
    if (s.len < 5) return false;
    if (s.ptr[0] != 'E' && s.ptr[0] != 'W' && s.ptr[0] != 'R') return false;
    for (int i = 1; i < 5; ++i) if (!is_digit(s.ptr[i])) return false;
    return s.len == 5 || !is_digit(s.ptr[5]);
}

// The position of the comment's '#' on a line of source, outside string literals; SIZE_MAX for none
static inline size_t vial_spec__comment_pos(str_t line) {
    char quote = 0;
    for (size_t i = 0; i < line.len; ++i) {
        const char c = line.ptr[i];
        if (quote) {
            if (c == '\\') { ++i; continue; }
            if (c == quote) quote = 0;
        } else if (c == '\'' || c == '"') {
            quote = c;
        } else if (c == '#') {
            return i;
        }
    }
    return SIZE_MAX;
}

static inline bool vial_spec__starts_word(str_t s, const char* word) {
    const size_t n = strlen(word);
    return s.len >= n && strncmp(s.ptr, word, n) == 0 && (s.len == n || !is_alpha(s.ptr[n]));
}

// Splits the text after "->" into the result and what follows it. The result ends at the first ',' ':' or ';' that
// is not inside brackets, parentheses or quotes, except where ", " is followed by another type, as in
// "set[424] and set[424]" (no comma there) or "float[45] Å² and float[45][45]".
static inline void vial_spec__split_result(str_t* result, str_t* note, str_t s) {
    int depth = 0;
    char quote = 0;
    size_t end = s.len;
    for (size_t i = 0; i < s.len; ++i) {
        const char c = s.ptr[i];
        if (quote) { if (c == quote) quote = 0; continue; }
        if (c == '\'' || c == '"') { quote = c; continue; }
        if (c == '[' || c == '(') { ++depth; continue; }
        if (c == ']' || c == ')') { --depth; continue; }
        if (depth == 0 && (c == ',' || c == ':' || c == ';')) { end = i; break; }
    }
    *result = str_trim(str_substr(s, 0, end));
    *note   = end < s.len ? str_trim(str_substr(s, end + 1, SIZE_MAX)) : (str_t){0};
}

// Reads one comment. Returns false when it starts like an expectation but does not follow the grammar.
static inline bool vial_spec__parse_expect(vial_expect_t* e, str_t c) {
    c = str_trim(c);
    e->comment = c;
    if (str_begins_with(c, STR_LIT("->"))) {
        e->kind = VIAL_EXPECT_RESULT;
        vial_spec__split_result(&e->result, &e->note, str_substr(c, 2, SIZE_MAX));
        return !str_empty(e->result);
    }
    const bool err  = vial_spec__starts_word(c, "error");
    const bool warn = vial_spec__starts_word(c, "warning");
    if (!err && !warn) {
        e->kind = VIAL_EXPECT_NONE;
        return true;
    }
    e->kind = err ? VIAL_EXPECT_ERROR : VIAL_EXPECT_WARNING;
    str_t s = str_trim(str_substr(c, err ? 5 : 7, SIZE_MAX));
    // One or more codes, separated by ", "
    while (vial_spec__is_code(s)) {
        if (e->num_codes == VIAL_SPEC_MAX_CODES) return false;
        memcpy(e->codes[e->num_codes++], s.ptr, 5);
        s = str_substr(s, 5, SIZE_MAX);
        if (s.len >= 2 && s.ptr[0] == ',' && vial_spec__is_code(str_trim(str_substr(s, 1, SIZE_MAX)))) {
            s = str_trim(str_substr(s, 1, SIZE_MAX));
        } else {
            break;
        }
    }
    if (e->num_codes == 0) return false;
    s = str_trim(s);
    if (str_begins_with(s, STR_LIT("->"))) {
        if (err) return false;  // an error has no result
        vial_spec__split_result(&e->result, &e->note, str_substr(s, 2, SIZE_MAX));
        return !str_empty(e->result);
    }
    if (!str_empty(s)) {
        if (s.ptr[0] != ':' && s.ptr[0] != ';') return false;
        e->note = str_trim(str_substr(s, 1, SIZE_MAX));
    }
    return true;
}

// ### RESULTS ###
// What a result says, read for the checker: one or more types (a tuple unpacked, "set[424] and set[424]"), or a
// static value such as 442, True or 'MET1/A'; and, from the note, the stage and the pass when they are written.

#define VIAL_SPEC_MAX_TYPES 4

typedef struct vial_result_t {
    bool  is_value;                     // a value, as 442 or 'MET1/A', rather than a type
    str_t value;
    int   num_types;
    str_t types[VIAL_SPEC_MAX_TYPES];
    bool  has_stage;                    // the note names a stage
    char  stage[24];                    // "static", "frame", "trajectory" or "trajectory, per frame"
    int   pass;                         // 0 when the note names none
} vial_result_t;

static inline bool vial_spec__is_value(str_t s) {
    if (str_empty(s)) return false;
    const char c = s.ptr[0];
    return is_digit(c) || c == '-' || c == '\'' || c == '"' || str_eq(s, STR_LIT("True")) || str_eq(s, STR_LIT("False")) ||
           str_eq(s, STR_LIT("None"));
}

// Reads the stage and pass at the start of a note: "frame", "trajectory, per frame, pass 2", "static: read from ..."
static inline void vial_spec__read_stage(vial_result_t* r, str_t note) {
    static const char* stages[] = {"trajectory, per frame", "trajectory", "static", "frame"};
    for (size_t i = 0; i < 4; ++i) {
        if (vial_spec__starts_word(note, stages[i])) {
            r->has_stage = true;
            snprintf(r->stage, sizeof(r->stage), "%s", stages[i]);
            note = str_trim(str_substr(note, strlen(stages[i]), SIZE_MAX));
            if (note.len && note.ptr[0] == ',') note = str_trim(str_substr(note, 1, SIZE_MAX));
            break;
        }
    }
    if (vial_spec__starts_word(note, "pass") && note.len > 5 && is_digit(note.ptr[5])) r->pass = note.ptr[5] - '0';
}

static inline bool vial_spec_read_result(vial_result_t* r, const vial_expect_t* e) {
    memset(r, 0, sizeof(*r));
    if (str_empty(e->result)) return false;
    vial_spec__read_stage(r, e->note);
    if (vial_spec__is_value(e->result)) {
        r->is_value = true;
        r->value = e->result;
        return true;
    }
    // Types joined by " and "
    str_t s = e->result;
    while (!str_empty(s)) {
        size_t at = SIZE_MAX;
        for (size_t i = 0; i + 5 <= s.len; ++i) {
            if (memcmp(s.ptr + i, " and ", 5) == 0) { at = i; break; }
        }
        if (r->num_types == VIAL_SPEC_MAX_TYPES) return false;
        r->types[r->num_types++] = str_trim(str_substr(s, 0, at));
        s = at == SIZE_MAX ? (str_t){0} : str_substr(s, at + 5, SIZE_MAX);
    }
    // Each must start with a dtype
    static const char* dtypes[] = {"bool", "int", "float64", "float", "str", "set", "range", "drawing", "("};
    for (int k = 0; k < r->num_types; ++k) {
        bool ok = false;
        for (size_t i = 0; i < ARRAY_SIZE(dtypes) && !ok; ++i) ok = str_begins_with(r->types[k], str_from_cstr(dtypes[i]));
        if (!ok) return false;
    }
    return true;
}

// Does a printed type match the type a result writes? A letter as an extent, [n], matches any extent, and a power
// may be written Å² or Å^2.
static inline size_t vial_spec__norm_type(char* out, size_t cap, str_t s) {
    size_t n = 0;
    for (size_t i = 0; i < s.len && n + 4 < cap; ++i) {
        // ² ³ as ^2 ^3
        if ((unsigned char)s.ptr[i] == 0xC2 && i + 1 < s.len && ((unsigned char)s.ptr[i + 1] == 0xB2 || (unsigned char)s.ptr[i + 1] == 0xB3)) {
            out[n++] = '^';
            out[n++] = (unsigned char)s.ptr[i + 1] == 0xB2 ? '2' : '3';
            ++i;
            continue;
        }
        out[n++] = s.ptr[i];
    }
    out[n] = '\0';
    return n;
}

static inline bool vial_spec_type_matches(str_t expected, str_t actual) {
    char e[256], a[256];
    const size_t ne = vial_spec__norm_type(e, sizeof(e), expected);
    const size_t na = vial_spec__norm_type(a, sizeof(a), actual);
    size_t i = 0, j = 0;
    while (i < ne && j < na) {
        if (e[i] == '[' && i + 2 < ne && is_alpha(e[i + 1]) && e[i + 1] != '?') {
            // [letter...]: any extent
            size_t k = i + 1;
            while (k < ne && is_alpha(e[k])) ++k;
            if (k < ne && e[k] == ']' && a[j] == '[') {
                while (j < na && a[j] != ']') ++j;
                i = k + 1;
                j += 1;
                continue;
            }
        }
        if (e[i] != a[j]) return false;
        ++i; ++j;
    }
    return i == ne && j == na;
}

// Reads the codes of the reference's diagnostics table: rows that start with "| E0101 |"
static inline void vial_spec__read_diag_codes(vial_spec_t* spec, md_allocator_i* alloc) {
    str_t text = spec->text;
    str_t line;
    while (str_extract_line(&line, &text)) {
        if (line.len >= 9 && line.ptr[0] == '|' && line.ptr[1] == ' ') {
            str_t s = str_substr(line, 2, SIZE_MAX);
            if (vial_spec__is_code(s) && s.len > 6 && s.ptr[5] == ' ' && s.ptr[6] == '|') {
                vial_code_t code = {0};
                memcpy(code.s, s.ptr, 5);
                md_array_push(spec->diag_codes, code, alloc);
            }
        }
    }
}

static inline bool vial_spec_has_code(const vial_spec_t* spec, const char* code) {
    for (size_t i = 0; i < md_array_size(spec->diag_codes); ++i) {
        if (strncmp(spec->diag_codes[i].s, code, 5) == 0) return true;
    }
    return false;
}

// Splits the reference into blocks and expectations. Views point into text, which must outlive the spec.
static inline void vial_spec_read(vial_spec_t* spec, str_t text, md_allocator_i* alloc) {
    memset(spec, 0, sizeof(*spec));
    spec->text = text;

    str_t rest = text;
    str_t line;
    str_t section = {0};
    int line_no = 0;
    char structure = 0;

    while (str_extract_line(&line, &rest)) {
        ++line_no;
        if (line.len > 0 && line.ptr[0] == '#') {
            str_t h = line;
            while (h.len && h.ptr[0] == '#') h = str_substr(h, 1, SIZE_MAX);
            section = str_trim(h);
            continue;
        }
        if (!str_eq(str_trim_end(line), STR_LIT("```python"))) continue;

        vial_block_t block = {.line = line_no + 1, .section = section, .structure = structure};
        const char* src_beg = rest.ptr;
        const char* src_end = rest.ptr;
        int line_in_block = 0;
        bool first_comment = true;

        while (str_extract_line(&line, &rest)) {
            ++line_no;
            if (str_eq(str_trim_end(line), STR_LIT("```"))) break;
            ++line_in_block;
            src_end = line.ptr + line.len;

            const size_t pos = vial_spec__comment_pos(line);
            if (pos == SIZE_MAX) continue;
            const str_t code_part = str_trim(str_substr(line, 0, pos));
            const str_t comment   = str_trim(str_substr(line, pos + 1, SIZE_MAX));

            // A whole-line comment that names the structure: "on T", "on P", "on T, with the names above", ...
            // The first one in a block names the block's; a later one switches for the lines after it.
            if (str_empty(code_part) && str_begins_with(comment, STR_LIT("on "))) {
                const str_t what = str_substr(comment, 3, SIZE_MAX);
                structure = (what.len >= 1 && (what.ptr[0] == 'T' || what.ptr[0] == 'P') &&
                             (what.len == 1 || what.ptr[1] == ',' || what.ptr[1] == ' ')) ? what.ptr[0] : 0;
                if (first_comment) {
                    block.structure = structure;
                    size_t loc;
                    if (str_find_str(&loc, comment, STR_LIT("with the names above"))) block.continues = true;
                    if (str_find_str(&loc, comment, STR_LIT("as filters")))           block.mode = VIAL_BLOCK_FILTERS;
                }
                first_comment = false;
                continue;
            }
            first_comment = false;
            if (str_empty(code_part)) continue;  // other whole-line comments are commentary

            vial_expect_t e = {.line = line_no, .line_in_block = line_in_block, .structure = structure};
            if (!vial_spec__parse_expect(&e, comment)) {
                spec->num_unparsed += 1;
                e.kind = VIAL_EXPECT_NONE;
            }
            if (e.kind != VIAL_EXPECT_NONE) md_array_push(block.expects, e, alloc);
        }
        block.src = (str_t){src_beg, (size_t)(src_end - src_beg)};
        md_array_push(spec->blocks, block, alloc);
    }

    vial_spec__read_diag_codes(spec, alloc);
}
