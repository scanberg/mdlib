#pragma once

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#include <core/md_str.h>

struct md_allocator_i;

// JSON reader.
//
// A document is parsed once into ONE allocation from the caller's allocator - a temp arena is the
// intended home - and is then read through VALUE HANDLES: small structs passed by value, each naming
// one value of the document. Nothing in the document is converted up front; a number becomes a
// double when it is asked for, through md_parse_f64, and a string is decoded when it is copied out.
// The document is read only once parsed, so any number of threads may read it at the same time.
//
// The document refers into the text it was parsed from, which must outlive it.
//
// MISSING VALUES
//
// Every lookup that finds nothing - a key an object does not have, an index past the end, a lookup
// on something that is not a container - gives the NONE value, and every function accepts NONE and
// answers as for an empty value: count 0, type MD_JSON_TYPE_NONE, a lookup in it gives NONE again,
// and reading a number or a string out of it fails. A chain of lookups therefore needs one check, at
// the end:
//
//   md_json_val_t q = md_json_at(md_json_get(md_json_get(atom, STR_LIT("multipoles")), STR_LIT("elements")), 0);
//   double charge = 0.0;
//   if (!md_json_f64(&charge, q)) { ... the file did not have it ... }
//
// The md_json_f64 / _i64 / _bool readers leave *out untouched when they fail, so a default can be
// written into the variable first and the result ignored.
//
// ITERATION
//
//   for (md_json_val_t e = md_json_first(arr); md_json_valid(e); e = md_json_next(e)) { ... }
//
// walks the elements of an array or the member VALUES of an object, in document order; md_json_key
// gives the key of a member value, as a string value. Stepping to the next sibling is O(1), however
// large the value being stepped over. md_json_at(v, i) steps i times, so walk rather than index when
// visiting all of a large container.
//
// WHAT IS ACCEPTED
//
// RFC 8259 JSON, checked in full: structure, string escapes and the number grammar - no leading
// zeros, no '+', no bare '.' - with no comments and no trailing commas. Any value can be the root.
// Two extensions, both for files written by Python's json module:
//   NaN, Infinity and -Infinity are numbers - json.dump writes them by default, and one of them must
//   not fail a whole file
//   a UTF-8 byte order mark in front of the root is skipped
// Bytes in strings are not checked to be UTF-8; they are passed through as they are. Within an
// object a key may repeat, in which case md_json_get finds the LAST one, as Python and JavaScript do.
// The text is limited to 4 GiB.

#ifdef __cplusplus
extern "C" {
#endif

typedef enum md_json_type_t {
    MD_JSON_TYPE_NONE = 0,      // no value: what a lookup that found nothing gives
    MD_JSON_TYPE_NULL,
    MD_JSON_TYPE_BOOL,
    MD_JSON_TYPE_NUMBER,
    MD_JSON_TYPE_STRING,
    MD_JSON_TYPE_ARRAY,
    MD_JSON_TYPE_OBJECT,
} md_json_type_t;

typedef struct md_json_t md_json_t;

// A value of a document. Its fields are the reader's business; the zero struct is NONE.
typedef struct md_json_val_t {
    const md_json_t* doc;
    uint32_t idx;
} md_json_val_t;

// Where and why parsing stopped. The message is a static string.
typedef struct md_json_error_t {
    const char* message;
    size_t offset;          // byte offset into the text
    size_t line;            // 1 based
    size_t column;          // 1 based, in bytes
} md_json_error_t;

// Parses text into a document allocated from alloc. NULL when the text is not JSON or the allocation
// fails: the reason goes to err when one is given, and to the log otherwise. The allocation is
// bounded by 16 bytes per ',', ':', '[' and '{' in the text, plus a small constant.
md_json_t* md_json_parse(str_t text, struct md_allocator_i* alloc, md_json_error_t* err);

// Returns the document to the allocator it was parsed from. Not needed for an arena that is rewound.
void md_json_free(md_json_t* doc, struct md_allocator_i* alloc);

// The root value; NONE for a NULL document
md_json_val_t md_json_root(const md_json_t* doc);

md_json_type_t md_json_type(md_json_val_t v);

static inline bool md_json_valid(md_json_val_t v) { return md_json_type(v) != MD_JSON_TYPE_NONE; }

// Elements of an array, members of an object, 0 for anything else
size_t md_json_count(md_json_val_t v);

// The value of the member with this key (compared after decoding escapes); NONE when absent
md_json_val_t md_json_get(md_json_val_t obj, str_t key);

// The i:th element of an array or member value of an object; NONE past the end. O(i).
md_json_val_t md_json_at(md_json_val_t container, size_t i);

// The first element / member value, and the one after a given one; NONE at the end
md_json_val_t md_json_first(md_json_val_t container);
md_json_val_t md_json_next(md_json_val_t v);

// The key of a member value, as a string value; NONE for a value that is not an object member
md_json_val_t md_json_key(md_json_val_t member);

// Numbers. md_json_f64 gives the correctly rounded double, NaN and +-inf included. md_json_i64 takes
// an integer literal within int64, or any other number whose double is a whole number below 2^53 in
// magnitude (so 3.0 and 1e3 are integers, 2.5 and 1e300 are not). Each leaves *out untouched and
// returns false otherwise.
bool md_json_f64 (double*  out, md_json_val_t v);
bool md_json_i64 (int64_t* out, md_json_val_t v);
bool md_json_bool(bool*    out, md_json_val_t v);

// Reads the elements of an array into out, up to cap of them, stopping at the first element that is
// not a number. Returns how many were written: compare it with md_json_count to require all of them.
size_t md_json_extract_f64(double* out, size_t cap, md_json_val_t arr);
size_t md_json_extract_f32(float*  out, size_t cap, md_json_val_t arr);

// Strings.
//   _raw   the bytes between the quotes as they are written, escapes and all - the string itself
//          whenever it has no backslash in it, which is the common case for keys and names. Empty
//          for anything that is not a string.
//   _eq    whether the decoded string equals str; false for anything that is not a string
//   _copy  decodes into buf as str_copy_to_char_buf copies: at most cap - 1 bytes, zero terminated
//          when cap > 0, and returns how many bytes it wrote. A string that does not fit is cut at a
//          character boundary, never through a UTF-8 sequence
//   (none) decodes into a zero terminated copy from alloc; empty for anything that is not a string
// Decoding turns \uXXXX escapes into UTF-8; an unpaired surrogate becomes U+FFFD.
str_t  md_json_string_raw (md_json_val_t v);
bool   md_json_string_eq  (md_json_val_t v, str_t str);
size_t md_json_string_copy(char* buf, size_t cap, md_json_val_t v);
str_t  md_json_string     (md_json_val_t v, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
