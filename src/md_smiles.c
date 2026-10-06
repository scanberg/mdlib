#include <md_smiles.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_str.h>

#include <stdarg.h>
#include <stdio.h>

// OpenSMILES reader. The grammar is followed with two liberties, both unambiguous: ring bonds may follow a branch as
// well as the atom (C(C)1CCCCC1), and a '.' may open a branch ('(.C)'), as the grammar allows but most readers do not.
// See md_smiles.h for what is kept and what is not.

#define NUM_RING_NUMBERS 100

typedef struct ring_bond_t {
    int32_t  atom;          // The atom which opened it, -1 when not open
    uint8_t  order;
    uint8_t  flags;
    uint32_t offset;
} ring_bond_t;

typedef struct branch_t {
    int32_t  atom;          // The atom the branch hangs from
    uint32_t num_atoms;     // Number of atoms when the branch was opened, to tell an empty branch
    uint32_t offset;
} branch_t;

typedef struct parser_t {
    const char* str;        // The string as given, offsets are relative to it
    const char* cur;
    const char* end;

    md_array(md_smiles_atom_t) atoms;
    md_array(md_smiles_bond_t) bonds;
    md_array(branch_t)         branches;
    struct md_allocator_i*     alloc;

    // The atom the next atom, ring bond or branch attaches to. -1 at the start and after a '.'
    int32_t prev;

    // A bond symbol waiting for the atom (or ring number) which follows it
    bool     bond_pending;
    uint8_t  bond_order;
    uint8_t  bond_flags;
    uint32_t bond_offset;

    ring_bond_t ring[NUM_RING_NUMBERS];

    md_smiles_error_t* err;
    bool failed;
} parser_t;

static bool fail(parser_t* p, const char* at, const char* fmt, ...) {
    if (!p->failed) {
        p->failed = true;
        if (p->err) {
            p->err->offset = (size_t)(at - p->str);
            va_list args;
            va_start(args, fmt);
            vsnprintf(p->err->message, sizeof(p->err->message), fmt, args);
            va_end(args);
        }
    }
    return false;
}

static inline uint32_t offset_of(const parser_t* p, const char* at) {
    return (uint32_t)(at - p->str);
}

static inline char peek(const parser_t* p) {
    return p->cur < p->end ? *p->cur : '\0';
}

static inline char peek_at(const parser_t* p, size_t i) {
    return p->cur + i < p->end ? p->cur[i] : '\0';
}

// is_digit and is_whitespace come from md_str.h
static inline bool is_upper(char c) { return 'A' <= c && c <= 'Z'; }
static inline bool is_lower(char c) { return 'a' <= c && c <= 'z'; }

// Reads up to max_digits digits. Returns false if there is none.
static bool parse_uint(uint32_t* out, parser_t* p, int max_digits) {
    if (!is_digit(peek(p))) return false;
    uint32_t val = 0;
    int n = 0;
    while (n < max_digits && is_digit(peek(p))) {
        val = val * 10 + (uint32_t)(*p->cur - '0');
        p->cur++;
        n++;
    }
    *out = val;
    return true;
}

static md_atomic_number_t element_from_symbol(const char* sym, size_t len) {
    str_t s = {sym, len};
    return md_atomic_number_from_symbol(s, false);
}

// ### ATOMS ###

static bool parse_organic(parser_t* p, md_smiles_atom_t* atom) {
    const char c  = peek(p);
    const char c1 = peek_at(p, 1);
    md_atomic_number_t z = 0;
    size_t len = 1;
    uint8_t flags = 0;

    switch (c) {
    case 'B': if (c1 == 'r') { z = 35; len = 2; } else { z = 5; } break;
    case 'C': if (c1 == 'l') { z = 17; len = 2; } else { z = 6; } break;
    case 'N': z = 7;  break;
    case 'O': z = 8;  break;
    case 'P': z = 15; break;
    case 'S': z = 16; break;
    case 'F': z = 9;  break;
    case 'I': z = 53; break;
    case 'b': z = 5;  flags = MD_SMILES_ATOM_AROMATIC; break;
    case 'c': z = 6;  flags = MD_SMILES_ATOM_AROMATIC; break;
    case 'n': z = 7;  flags = MD_SMILES_ATOM_AROMATIC; break;
    case 'o': z = 8;  flags = MD_SMILES_ATOM_AROMATIC; break;
    case 'p': z = 15; flags = MD_SMILES_ATOM_AROMATIC; break;
    case 's': z = 16; flags = MD_SMILES_ATOM_AROMATIC; break;
    case '*': z = 0;  break;
    default:
        if (is_upper(c) || is_lower(c)) {
            return fail(p, p->cur, "'%c' is not an atom of the organic subset, write it in brackets: [%c...]", c, c);
        }
        return fail(p, p->cur, "Expected an atom");
    }

    atom->z = z;
    atom->flags = flags;
    p->cur += len;
    return true;
}

static bool parse_bracket(parser_t* p, md_smiles_atom_t* atom) {
    const char* open = p->cur;
    ASSERT(peek(p) == '[');
    p->cur++;

    atom->flags = MD_SMILES_ATOM_BRACKET;

    uint32_t isotope;
    if (parse_uint(&isotope, p, 5)) {
        if (isotope > UINT16_MAX) return fail(p, open + 1, "Isotope out of range");
        atom->isotope = (uint16_t)isotope;
    }

    // Symbol
    const char* sym = p->cur;
    const char c0 = peek(p);
    const char c1 = peek_at(p, 1);
    if (c0 == '*') {
        atom->z = 0;
        p->cur += 1;
    } else if (is_lower(c0)) {
        // Aromatic: se, as, te before the single letter ones
        if ((c0 == 's' && c1 == 'e') || (c0 == 'a' && c1 == 's') || (c0 == 't' && c1 == 'e')) {
            const char up[2] = {(char)(c0 - 'a' + 'A'), c1};
            atom->z = element_from_symbol(up, 2);
            p->cur += 2;
        } else if (c0 == 'b' || c0 == 'c' || c0 == 'n' || c0 == 'o' || c0 == 'p' || c0 == 's') {
            const char up = (char)(c0 - 'a' + 'A');
            atom->z = element_from_symbol(&up, 1);
            p->cur += 1;
        } else {
            return fail(p, sym, "Unknown aromatic element in brackets");
        }
        atom->flags |= MD_SMILES_ATOM_AROMATIC;
    } else if (is_upper(c0)) {
        md_atomic_number_t z = 0;
        if (is_lower(c1)) {
            z = element_from_symbol(sym, 2);
            if (z) p->cur += 2;
        }
        if (!z) {
            z = element_from_symbol(sym, 1);
            if (!z) {
                return fail(p, sym, is_lower(c1) ? "Unknown element '%c%c'" : "Unknown element '%c'", c0, c1);
            }
            p->cur += 1;
        }
        atom->z = z;
    } else {
        return fail(p, sym, "Expected an element symbol after '['");
    }

    // Chirality
    if (peek(p) == '@') {
        p->cur++;
        if (peek(p) == '@') {
            p->cur++;
            atom->flags |= MD_SMILES_ATOM_CHIRAL_CW;
        } else {
            const char a = peek(p), b = peek_at(p, 1);
            const bool ext = (a == 'T' && b == 'H') || (a == 'A' && b == 'L') || (a == 'S' && b == 'P') || (a == 'T' && b == 'B') || (a == 'O' && b == 'H');
            if (ext) {
                const char* at = p->cur;
                p->cur += 2;
                uint32_t num;
                if (!parse_uint(&num, p, 2)) return fail(p, at, "Expected a number after the chirality class");
                atom->flags |= MD_SMILES_ATOM_CHIRAL_EXT;
                atom->chiral_class = (uint8_t)num;
            } else {
                atom->flags |= MD_SMILES_ATOM_CHIRAL_CCW;
            }
        }
    }

    // Hydrogen count
    if (peek(p) == 'H') {
        p->cur++;
        uint32_t h = 1;
        parse_uint(&h, p, 2);
        if (h > UINT8_MAX) return fail(p, p->cur, "Hydrogen count out of range");
        atom->h_count = (uint8_t)h;
    }

    // Charge: +, ++, +n and the same for -
    const char sign = peek(p);
    if (sign == '+' || sign == '-') {
        const char* at = p->cur;
        p->cur++;
        int32_t magnitude = 1;
        uint32_t num;
        if (parse_uint(&num, p, 2)) {
            magnitude = (int32_t)num;
        } else {
            while (peek(p) == sign) {
                p->cur++;
                magnitude++;
            }
        }
        if (magnitude > 15) return fail(p, at, "Charge out of range");
        atom->charge = (int8_t)(sign == '+' ? magnitude : -magnitude);
    }

    // Atom class
    if (peek(p) == ':') {
        const char* at = p->cur;
        p->cur++;
        uint32_t cls;
        if (!parse_uint(&cls, p, 5) || cls > UINT16_MAX) return fail(p, at, "Expected an atom class after ':'");
        atom->atom_class = (uint16_t)cls;
    }

    if (peek(p) != ']') {
        if (p->cur >= p->end) return fail(p, open, "Unclosed '['");
        return fail(p, p->cur, "Unexpected '%c' in bracket atom", peek(p));
    }
    p->cur++;
    return true;
}

// ### BONDS ###

static bool has_bond(const parser_t* p, int32_t a, int32_t b) {
    const size_t n = md_array_size(p->bonds);
    for (size_t i = 0; i < n; ++i) {
        const md_smiles_bond_t* bond = &p->bonds[i];
        if (((int32_t)bond->a == a && (int32_t)bond->b == b) || ((int32_t)bond->a == b && (int32_t)bond->b == a)) return true;
    }
    return false;
}

static void push_bond(parser_t* p, int32_t a, int32_t b, uint8_t order, uint8_t flags) {
    md_smiles_bond_t bond = {.a = (uint32_t)a, .b = (uint32_t)b, .order = order, .flags = flags};
    md_array_push(p->bonds, bond, p->alloc);
}

static bool bond_symbol(uint8_t* order, uint8_t* flags, char c) {
    *flags = 0;
    switch (c) {
    case '-':  *order = MD_SMILES_BOND_SINGLE;    return true;
    case '=':  *order = MD_SMILES_BOND_DOUBLE;    return true;
    case '#':  *order = MD_SMILES_BOND_TRIPLE;    return true;
    case '$':  *order = MD_SMILES_BOND_QUADRUPLE; return true;
    case ':':  *order = MD_SMILES_BOND_AROMATIC;  return true;
    case '/':  *order = MD_SMILES_BOND_SINGLE; *flags = MD_SMILES_BOND_UP;   return true;
    case '\\': *order = MD_SMILES_BOND_SINGLE; *flags = MD_SMILES_BOND_DOWN; return true;
    default: return false;
    }
}

static bool parse_ring_bond(parser_t* p) {
    const char* at = p->cur;
    uint32_t num;
    if (peek(p) == '%') {
        p->cur++;
        if (!is_digit(peek(p)) || !is_digit(peek_at(p, 1))) return fail(p, at, "Expected two digits after '%%'");
        parse_uint(&num, p, 2);
    } else {
        parse_uint(&num, p, 1);
    }
    ASSERT(num < NUM_RING_NUMBERS);

    if (p->prev < 0) return fail(p, at, "Ring bond without an atom before it");

    ring_bond_t* ring = &p->ring[num];
    const uint8_t order = p->bond_pending ? p->bond_order : (uint8_t)MD_SMILES_BOND_IMPLICIT;
    const uint8_t flags = p->bond_pending ? p->bond_flags : 0;
    p->bond_pending = false;

    if (ring->atom < 0) {
        ring->atom   = p->prev;
        ring->order  = order;
        ring->flags  = flags;
        ring->offset = offset_of(p, at);
        return true;
    }

    // Closing
    if (ring->atom == p->prev) return fail(p, at, "Ring bond %u bonds an atom to itself", num);
    if (has_bond(p, ring->atom, p->prev)) return fail(p, at, "Ring bond %u duplicates a bond", num);
    if (ring->order != MD_SMILES_BOND_IMPLICIT && order != MD_SMILES_BOND_IMPLICIT && ring->order != order) {
        return fail(p, at, "The two ends of ring bond %u give different bonds", num);
    }
    const uint8_t ring_order = ring->order != MD_SMILES_BOND_IMPLICIT ? ring->order : order;
    const uint8_t ring_flags = (flags ? flags : ring->flags) | MD_SMILES_BOND_RING;
    push_bond(p, ring->atom, p->prev, ring_order, ring_flags);
    ring->atom = -1;
    return true;
}

// ### LINE ###

static bool parse_atom(parser_t* p) {
    const char* at = p->cur;
    md_smiles_atom_t atom = {0};
    if (peek(p) == '[') {
        if (!parse_bracket(p, &atom)) return false;
    } else {
        if (!parse_organic(p, &atom)) return false;
    }
    atom.offset = offset_of(p, at);

    const int32_t idx = (int32_t)md_array_size(p->atoms);
    md_array_push(p->atoms, atom, p->alloc);

    if (p->prev >= 0) {
        const uint8_t order = p->bond_pending ? p->bond_order : (uint8_t)MD_SMILES_BOND_IMPLICIT;
        const uint8_t flags = p->bond_pending ? p->bond_flags : 0;
        push_bond(p, p->prev, idx, order, flags);
    }
    p->bond_pending = false;
    p->prev = idx;
    return true;
}

static bool parse_line(parser_t* p) {
    while (p->cur < p->end) {
        const char c = peek(p);
        uint8_t order, flags;

        if (c == '[' || c == '*' || is_upper(c) || is_lower(c)) {
            if (!parse_atom(p)) return false;
        }
        else if (bond_symbol(&order, &flags, c)) {
            if (p->prev < 0) return fail(p, p->cur, "Bond without an atom before it");
            if (p->bond_pending) return fail(p, p->cur, "Two bonds in a row");
            p->bond_pending = true;
            p->bond_order   = order;
            p->bond_flags   = flags;
            p->bond_offset  = offset_of(p, p->cur);
            p->cur++;
        }
        else if (is_digit(c) || c == '%') {
            if (!parse_ring_bond(p)) return false;
        }
        else if (c == '(') {
            if (p->prev < 0) return fail(p, p->cur, "Branch without an atom before it");
            if (p->bond_pending) return fail(p, p->cur, "A bond before a branch goes inside it: C(=O)");
            branch_t branch = {.atom = p->prev, .num_atoms = (uint32_t)md_array_size(p->atoms), .offset = offset_of(p, p->cur)};
            md_array_push(p->branches, branch, p->alloc);
            p->cur++;
        }
        else if (c == ')') {
            if (md_array_size(p->branches) == 0) return fail(p, p->cur, "Unbalanced ')'");
            if (p->bond_pending) return fail(p, p->str + p->bond_offset, "Bond without an atom after it");
            const branch_t branch = md_array_back(p->branches);
            if (branch.num_atoms == md_array_size(p->atoms)) return fail(p, p->cur, "Empty branch");
            if (p->prev < 0) return fail(p, p->cur - 1, "'.' without an atom after it");
            md_array_pop(p->branches);
            p->prev = branch.atom;
            p->cur++;
        }
        else if (c == '.') {
            if (p->prev < 0) return fail(p, p->cur, "'.' without an atom before it");
            if (p->bond_pending) return fail(p, p->str + p->bond_offset, "Bond without an atom after it");
            p->prev = -1;
            p->cur++;
        }
        else if (is_whitespace(c)) {
            return fail(p, p->cur, "Unexpected whitespace");
        }
        else {
            return fail(p, p->cur, "Unexpected '%c'", c);
        }
    }

    if (p->bond_pending) return fail(p, p->str + p->bond_offset, "Bond without an atom after it");
    if (md_array_size(p->branches)) return fail(p, p->str + md_array_back(p->branches).offset, "Unclosed '('");
    if (p->prev < 0 && md_array_size(p->atoms)) return fail(p, p->cur - 1, "'.' without an atom after it");
    for (uint32_t i = 0; i < NUM_RING_NUMBERS; ++i) {
        if (p->ring[i].atom >= 0) return fail(p, p->str + p->ring[i].offset, "Ring bond %u is never closed", i);
    }
    return true;
}

static uint32_t find_root(uint32_t* parent, uint32_t i) {
    while (parent[i] != i) {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    return i;
}

static size_t count_components(size_t num_atoms, const md_smiles_bond_t* bonds, size_t num_bonds, struct md_allocator_i* alloc) {
    if (num_atoms == 0) return 0;
    uint32_t* parent = md_alloc(alloc, sizeof(uint32_t) * num_atoms);
    for (uint32_t i = 0; i < (uint32_t)num_atoms; ++i) parent[i] = i;
    size_t count = num_atoms;
    for (size_t i = 0; i < num_bonds; ++i) {
        const uint32_t a = find_root(parent, bonds[i].a);
        const uint32_t b = find_root(parent, bonds[i].b);
        if (a != b) {
            parent[a] = b;
            count--;
        }
    }
    md_free(alloc, parent, sizeof(uint32_t) * num_atoms);
    return count;
}

bool md_smiles_parse(md_smiles_t* out, str_t str, struct md_allocator_i* alloc, md_smiles_error_t* err) {
    ASSERT(out);
    ASSERT(alloc);
    MEMSET(out, 0, sizeof(md_smiles_t));
    if (err) MEMSET(err, 0, sizeof(md_smiles_error_t));

    parser_t p = {
        .str   = str.ptr,
        .cur   = str.ptr,
        .end   = str.ptr + str.len,
        .alloc = alloc,
        .prev  = -1,
        .err   = err,
    };
    for (size_t i = 0; i < NUM_RING_NUMBERS; ++i) {
        p.ring[i].atom = -1;
    }

    if (!str.ptr) {
        fail(&p, p.cur, "Empty SMILES");
        return false;
    }

    // Leading and trailing whitespace (and a terminating null) is not part of the pattern
    while (p.cur < p.end && is_whitespace(*p.cur)) p.cur++;
    while (p.end > p.cur && (is_whitespace(p.end[-1]) || p.end[-1] == '\0')) p.end--;

    bool ok = false;
    if (p.cur == p.end) {
        fail(&p, p.cur, "Empty SMILES");
    } else {
        ok = parse_line(&p);
    }

    md_array_free(p.branches, alloc);

    if (!ok) {
        md_array_free(p.atoms, alloc);
        md_array_free(p.bonds, alloc);
        return false;
    }

    out->num_atoms = md_array_size(p.atoms);
    out->atoms     = p.atoms;
    out->num_bonds = md_array_size(p.bonds);
    out->bonds     = p.bonds;
    out->num_components = count_components(out->num_atoms, out->bonds, out->num_bonds, alloc);
    out->alloc = alloc;
    return true;
}

void md_smiles_free(md_smiles_t* smiles) {
    ASSERT(smiles);
    if (smiles->alloc) {
        md_array_free(smiles->atoms, smiles->alloc);
        md_array_free(smiles->bonds, smiles->alloc);
    }
    MEMSET(smiles, 0, sizeof(md_smiles_t));
}
