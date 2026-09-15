#include <md_molden.h>

#include <md_system.h>
#include <md_types.h>
#include <md_gto.h>
#include <md_qm.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_parse.h>
#include <core/md_str.h>
#include <core/md_unit.h>
#include <core/md_vec_math.h>

#include <math.h>
#include <string.h>

#define MOLDEN_BOHR_TO_ANGSTROM 0.529177210903

// The widest contraction any published basis set uses is a long way below this; it is here so a
// shell can be read onto the stack and a corrupt count is refused rather than allocated.
#define MOLDEN_MAX_PRIMITIVES 64

// ---------------------------------------------------------------------------
// The reader's own representation
// ---------------------------------------------------------------------------
// Parse time only. Nothing here outlives md_molden_system_init_from_str, which is the whole point:
// what the file carried leaves this file as attributes on a system, never as a struct a consumer
// would then have to keep a reader alive to read.

typedef struct molden_shell_t {
    uint32_t atom_idx;
    uint32_t l;
    uint32_t primitive_offset;  // into molden_t::alpha / ::coeff
    uint32_t num_primitives;
} molden_shell_t;

typedef struct molden_orbital_t {
    double  energy;      // hartree
    double  occupation;
    str_t   symmetry;    // interned into the parse arena
    bool    beta;
    double* coefficient; // [num_ao_file], in the FILE's AO order and convention
} molden_orbital_t;

typedef struct molden_t {
    md_allocator_i* alloc;

    str_t title;
    str_t program;

    // Angular convention per angular momentum: true when the file stores the pure/spherical
    // (2l+1) set. Cartesian is the format's default and is what an unmarked file means.
    bool pure[MD_GTO_MAX_ANGULAR_MOMENTUM + 1];

    md_array(uint8_t)          atomic_number;
    md_array(dvec3_t)          coord;       // Angstrom
    md_array(str_t)            atom_label;

    md_array(molden_shell_t)   shell;
    md_array(double)           alpha;
    md_array(double)           coeff;       // as read; normalised in molden_normalise_basis

    md_array(molden_orbital_t) orbital;

    md_array(double)           frequency;     // cm^-1
    md_array(double)           ir_intensity;  // km/mol
    md_array(dvec3_t)          normal_mode;   // [num_modes * num_atoms]
    size_t                     num_modes;

    bool has_gto;
} molden_t;

static inline md_unit_t molden_unit_wavenumber(void) { return md_unit_pow(md_unit_scl(md_unit_meter(), 1.0e-2), -1); }
static inline md_unit_t molden_unit_km_per_mol(void) { return md_unit_div(md_unit_scl(md_unit_meter(), 1.0e3), md_unit_mole()); }

// ---------------------------------------------------------------------------
// Angular conventions
// ---------------------------------------------------------------------------

// The monomial exponents of each Cartesian function, in MOLDEN's order, which is not md_gto's:
//
//    6D : xx yy zz xy xz yz
//   10F : xxx yyy zzz xyy xxy xxz xzz yzz yyz xyz
//   15G : xxxx yyyy zzzz xxxy xxxz yyyx yyyz zzzx zzzy xxyy xxzz yyzz xxyz yyxz zzxy
//
// Kept as exponents rather than as a permutation table so it can be read against the format
// specification one entry at a time, and resolved to md_gto's indices by search - 15 entries, once
// per shell, against a mapping that would otherwise be two tables to keep in step.
typedef uint8_t molden_lmn_t[3];

static const molden_lmn_t molden_cart_S[1]  = {{0,0,0}};
static const molden_lmn_t molden_cart_P[3]  = {{1,0,0},{0,1,0},{0,0,1}};
static const molden_lmn_t molden_cart_D[6]  = {{2,0,0},{0,2,0},{0,0,2},{1,1,0},{1,0,1},{0,1,1}};
static const molden_lmn_t molden_cart_F[10] = {{3,0,0},{0,3,0},{0,0,3},{1,2,0},{2,1,0},{2,0,1},{1,0,2},{0,1,2},{0,2,1},{1,1,1}};
static const molden_lmn_t molden_cart_G[15] = {{4,0,0},{0,4,0},{0,0,4},{3,1,0},{3,0,1},{1,3,0},{0,3,1},{1,0,3},{0,1,3},
                                               {2,2,0},{2,0,2},{0,2,2},{2,1,1},{1,2,1},{1,1,2}};

static const molden_lmn_t* molden_cart_order(uint32_t l) {
    switch (l) {
    case 0: return molden_cart_S;
    case 1: return molden_cart_P;
    case 2: return molden_cart_D;
    case 3: return molden_cart_F;
    case 4: return molden_cart_G;
    default: return NULL;
    }
}

// md_gto's Cartesian index for a monomial, or the AO count when there is none (which cannot happen
// for a table entry of the right l, and is a loud out of range rather than a silent 0 if it does).
static uint32_t molden_cart_index(uint32_t l, const molden_lmn_t lmn) {
    const uint32_t n = md_gto_num_cart_ao(l);
    for (uint32_t c = 0; c < n; ++c) {
        int i, j, k;
        if (md_gto_cart_ijk(&i, &j, &k, l, c) && i == lmn[0] && j == lmn[1] && k == lmn[2]) {
            return c;
        }
    }
    return n;
}

// Molden's spherical order is m = 0, +1, -1, +2, -2, ... md_gto's is m ascending, -l .. +l. This is
// the file index of md_gto's function 'idx', so a caller scatters rather than gathers.
static uint32_t molden_sph_to_md_index(uint32_t l, uint32_t file_idx) {
    // file_idx 0 -> m = 0; 1 -> +1; 2 -> -1; 3 -> +2; 4 -> -2; ...
    const int m = (file_idx == 0) ? 0 : (((file_idx & 1) ? 1 : -1) * (int)((file_idx + 1) / 2));
    return (uint32_t)(m + (int)l);
}

// Whether this shell's coefficients are the pure set. s and p are the same 1 and 3 functions either
// way, and Molden orders p as x, y, z - which IS md_gto's Cartesian order - so they take the
// Cartesian path unconditionally and no conversion happens for them at all.
static bool molden_shell_is_pure(const molden_t* molden, uint32_t l) {
    return l >= 2 && l <= MD_GTO_MAX_ANGULAR_MOMENTUM && molden->pure[l];
}

static uint32_t molden_shell_num_file_ao(const molden_t* molden, uint32_t l) {
    return molden_shell_is_pure(molden, l) ? md_gto_num_sph_ao(l) : md_gto_num_cart_ao(l);
}

static size_t molden_num_file_ao(const molden_t* molden) {
    size_t n = 0;
    for (size_t i = 0; i < md_array_size(molden->shell); ++i) {
        n += molden_shell_num_file_ao(molden, molden->shell[i].l);
    }
    return n;
}

// ---------------------------------------------------------------------------
// Parsing
// ---------------------------------------------------------------------------

// A Molden number may carry a Fortran 'D' exponent (1.0D-03). parse_float does not know about it,
// and the alternative - rewriting the buffer - would mutate the caller's string.
static double molden_parse_float(str_t tok) {
    char buf[64];
    size_t n = str_copy_to_char_buf(buf, sizeof(buf), tok);
    for (size_t i = 0; i < n; ++i) {
        if (buf[i] == 'D' || buf[i] == 'd') buf[i] = 'e';
    }
    return parse_float(str_from_cstrn(buf, n));
}

// "[Atoms] (Angs)" -> section "atoms", argument "(Angs)". Returns false for anything that is not a
// section header, which is every other line in the file.
static bool molden_section_header(str_t* out_name, str_t* out_arg, str_t line) {
    line = str_trim(line);
    if (line.len < 2 || line.ptr[0] != '[') {
        return false;
    }
    size_t close;
    if (!str_find_char(&close, line, ']')) {
        return false;
    }
    *out_name = str_trim(str_substr(line, 1, close - 1));
    *out_arg  = str_trim(str_substr(line, close + 1, SIZE_MAX));
    return true;
}

static int molden_shell_type_to_l(str_t tok) {
    if (tok.len != 1) {
        return -1;
    }
    switch (tok.ptr[0]) {
    case 's': case 'S': return 0;
    case 'p': case 'P': return 1;
    case 'd': case 'D': return 2;
    case 'f': case 'F': return 3;
    case 'g': case 'G': return 4;
    default: return -1;
    }
}

static void molden_push_shell(molden_t* molden, uint32_t atom_idx, uint32_t l, const double* alpha, const double* coeff, size_t count) {
    molden_shell_t shell = {
        .atom_idx         = atom_idx,
        .l                = l,
        .primitive_offset = (uint32_t)md_array_size(molden->alpha),
        .num_primitives   = (uint32_t)count,
    };
    for (size_t i = 0; i < count; ++i) {
        md_array_push(molden->alpha, alpha[i], molden->alloc);
        md_array_push(molden->coeff, coeff[i], molden->alloc);
    }
    md_array_push(molden->shell, shell, molden->alloc);
}

static bool molden_parse_atoms(molden_t* molden, md_buffered_reader_t* reader, str_t arg) {
    // "[Atoms] Angs", "[Atoms] (AU)", "[ATOMS] AU" - the parentheses are decoration. Angstrom is
    // the default, which is what an argument-less header means.
    str_t a = arg;
    if (a.len >= 2 && a.ptr[0] == '(' && a.ptr[a.len - 1] == ')') {
        a = str_trim(str_substr(a, 1, a.len - 2));
    }
    const bool bohr = str_eq_cstr_ignore_case(a, "au") || str_eq_cstr_ignore_case(a, "bohr");

    str_t line;
    while (md_buffered_reader_peek_line(&line, reader)) {
        str_t name, sec_arg;
        if (molden_section_header(&name, &sec_arg, line)) {
            break;
        }
        md_buffered_reader_skip_line(reader);

        str_t tok[8];
        str_t rest = line;
        const size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &rest);
        if (num_tok < 6) {
            continue;
        }
        // symbol index atomic_number x y z
        md_atomic_number_t z = (md_atomic_number_t)parse_int(tok[2]);
        if (z == 0) {
            z = md_atomic_number_from_symbol(tok[0], true);
        }
        dvec3_t xyz = {
            molden_parse_float(tok[3]),
            molden_parse_float(tok[4]),
            molden_parse_float(tok[5]),
        };
        if (bohr) {
            xyz = dvec3_mul1(xyz, MOLDEN_BOHR_TO_ANGSTROM);
        }
        md_array_push(molden->atomic_number, (uint8_t)z, molden->alloc);
        md_array_push(molden->coord, xyz, molden->alloc);
        md_array_push(molden->atom_label, str_copy(tok[0], molden->alloc), molden->alloc);
    }
    return md_array_size(molden->atomic_number) > 0;
}

// [GTO] is per atom: an "<atom index> 0" line, then shells until a blank line.
//
// 'sp' is a shell with TWO coefficient columns sharing one set of exponents; it becomes an s shell
// and a p shell in that order, which is also the order its four AOs appear in, so nothing downstream
// has to know it existed.
static bool molden_parse_gto(molden_t* molden, md_buffered_reader_t* reader) {
    const size_t num_atoms = md_array_size(molden->atomic_number);
    if (num_atoms == 0) {
        MD_LOG_ERROR("MOLDEN: [GTO] appears before [Atoms]");
        return false;
    }

    double alpha[MOLDEN_MAX_PRIMITIVES];
    double coeff[MOLDEN_MAX_PRIMITIVES];
    double coeff2[MOLDEN_MAX_PRIMITIVES];

    uint32_t atom_idx = 0;
    bool have_atom = false;

    str_t line;
    while (md_buffered_reader_peek_line(&line, reader)) {
        str_t name, arg;
        if (molden_section_header(&name, &arg, line)) {
            break;
        }
        md_buffered_reader_skip_line(reader);

        str_t rest = str_trim(line);
        if (str_empty(rest)) {
            have_atom = false;
            continue;
        }

        str_t tok[4];
        str_t scan = rest;
        const size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &scan);
        if (num_tok == 0) {
            continue;
        }

        const int  l     = molden_shell_type_to_l(tok[0]);
        const bool is_sp = str_eq_cstr_ignore_case(tok[0], "sp");

        // A line is either a shell - which starts with its type letter - or the "<atom index> 0"
        // header of the next atom's block. Nothing else is legal here, and an unknown letter is a
        // shell type this library cannot evaluate rather than something to skip past quietly.
        if (l < 0 && !is_sp) {
            if (!is_int(tok[0])) {
                MD_LOG_ERROR("MOLDEN: [GTO] shell type '" STR_FMT "' is not one this library evaluates (s, p, d, f, g and sp)", STR_ARG(tok[0]));
                return false;
            }
            const int64_t idx = parse_int(tok[0]) - 1;
            if (idx < 0 || (size_t)idx >= num_atoms) {
                MD_LOG_ERROR("MOLDEN: [GTO] names atom %d, the file has %zu", (int)(idx + 1), num_atoms);
                return false;
            }
            atom_idx  = (uint32_t)idx;
            have_atom = true;
            continue;
        }

        if (!have_atom) {
            MD_LOG_ERROR("MOLDEN: [GTO] states a shell before it names an atom");
            return false;
        }

        if (num_tok < 2) {
            MD_LOG_ERROR("MOLDEN: malformed shell line '" STR_FMT "'", STR_ARG(rest));
            return false;
        }
        const size_t num_prim = (size_t)parse_int(tok[1]);
        if (num_prim == 0 || num_prim > MOLDEN_MAX_PRIMITIVES) {
            MD_LOG_ERROR("MOLDEN: shell with %zu primitives, which is outside 1..%d", num_prim, MOLDEN_MAX_PRIMITIVES);
            return false;
        }

        for (size_t p = 0; p < num_prim; ++p) {
            str_t prim_line;
            if (!md_buffered_reader_extract_line(&prim_line, reader)) {
                MD_LOG_ERROR("MOLDEN: [GTO] ends inside a contraction");
                return false;
            }
            str_t ptok[4];
            str_t pscan = prim_line;
            const size_t n = extract_tokens(ptok, ARRAY_SIZE(ptok), &pscan);
            if (n < 2) {
                MD_LOG_ERROR("MOLDEN: malformed primitive line '" STR_FMT "'", STR_ARG(str_trim(prim_line)));
                return false;
            }
            alpha[p]  = molden_parse_float(ptok[0]);
            coeff[p]  = molden_parse_float(ptok[1]);
            coeff2[p] = (n > 2) ? molden_parse_float(ptok[2]) : 0.0;
        }

        if (is_sp) {
            molden_push_shell(molden, atom_idx, 0, alpha, coeff,  num_prim);
            molden_push_shell(molden, atom_idx, 1, alpha, coeff2, num_prim);
        } else {
            molden_push_shell(molden, atom_idx, (uint32_t)l, alpha, coeff, num_prim);
        }
    }

    molden->has_gto = md_array_size(molden->shell) > 0;
    return molden->has_gto;
}

// [MO] is a run of blocks, each a few key lines followed by "<ao index> <value>" pairs. The keys
// may come in any order and any of them may be absent, so a block ends when a key line appears
// after at least one coefficient has been read - which is the only marker the format gives.
static bool molden_parse_mo(molden_t* molden, md_buffered_reader_t* reader, size_t num_ao) {
    if (num_ao == 0) {
        MD_LOG_ERROR("MOLDEN: [MO] appears before a basis is known");
        return false;
    }

    molden_orbital_t orb = {0};
    bool have_orb = false;
    size_t num_coeff = 0;

    str_t line;
    while (md_buffered_reader_peek_line(&line, reader)) {
        str_t name, arg;
        if (molden_section_header(&name, &arg, line)) {
            break;
        }
        md_buffered_reader_skip_line(reader);

        str_t s = str_trim(line);
        if (str_empty(s)) {
            continue;
        }

        size_t eq;
        const bool is_key = str_find_char(&eq, s, '=') && eq <= 8;
        if (is_key) {
            if (have_orb && num_coeff > 0) {
                md_array_push(molden->orbital, orb, molden->alloc);
                MEMSET(&orb, 0, sizeof(orb));
                have_orb  = false;
                num_coeff = 0;
            }
            if (!have_orb) {
                orb = (molden_orbital_t){0};
                orb.coefficient = (double*)md_alloc(molden->alloc, sizeof(double) * num_ao);
                MEMSET(orb.coefficient, 0, sizeof(double) * num_ao);
                have_orb = true;
            }

            str_t key = str_trim(str_substr(s, 0, eq));
            str_t val = str_trim(str_substr(s, eq + 1, SIZE_MAX));
            if (str_eq_cstr_ignore_case(key, "ene")) {
                orb.energy = molden_parse_float(val);
            } else if (str_eq_cstr_ignore_case(key, "spin")) {
                orb.beta = str_eq_cstr_ignore_case(val, "beta");
            } else if (str_eq_cstr_ignore_case(key, "occup")) {
                orb.occupation = molden_parse_float(val);
            } else if (str_eq_cstr_ignore_case(key, "sym")) {
                orb.symmetry = str_copy(val, molden->alloc);
            }
            continue;
        }

        if (!have_orb) {
            continue;
        }

        str_t tok[2];
        str_t scan = s;
        if (extract_tokens(tok, ARRAY_SIZE(tok), &scan) < 2) {
            continue;
        }
        const int64_t ao = parse_int(tok[0]) - 1;
        if (ao < 0 || (size_t)ao >= num_ao) {
            MD_LOG_ERROR("MOLDEN: [MO] names atomic orbital %d, the basis has %zu", (int)(ao + 1), num_ao);
            return false;
        }
        orb.coefficient[ao] = molden_parse_float(tok[1]);
        num_coeff += 1;
    }

    if (have_orb && num_coeff > 0) {
        md_array_push(molden->orbital, orb, molden->alloc);
    }
    return md_array_size(molden->orbital) > 0;
}

// [FR-NORM-COORD] is "vibration <n>" followed by one displacement per atom, in bohr.
static void molden_parse_normal_modes(molden_t* molden, md_buffered_reader_t* reader) {
    const size_t num_atoms = md_array_size(molden->atomic_number);
    str_t line;
    while (md_buffered_reader_peek_line(&line, reader)) {
        str_t name, arg;
        if (molden_section_header(&name, &arg, line)) {
            break;
        }
        md_buffered_reader_skip_line(reader);

        str_t s = str_trim(line);
        if (str_empty(s)) {
            continue;
        }
        if (!str_eq_cstr_n_ignore_case(s, "vibration", 9)) {
            continue;
        }
        for (size_t a = 0; a < num_atoms; ++a) {
            str_t disp_line;
            if (!md_buffered_reader_extract_line(&disp_line, reader)) {
                return;
            }
            str_t tok[4];
            str_t scan = disp_line;
            dvec3_t d = {0};
            if (extract_tokens(tok, ARRAY_SIZE(tok), &scan) >= 3) {
                d.x = molden_parse_float(tok[0]);
                d.y = molden_parse_float(tok[1]);
                d.z = molden_parse_float(tok[2]);
            }
            md_array_push(molden->normal_mode, d, molden->alloc);
        }
        molden->num_modes += 1;
    }
}

static void molden_parse_number_series(md_array(double)* out, md_buffered_reader_t* reader, md_allocator_i* alloc) {
    str_t line;
    while (md_buffered_reader_peek_line(&line, reader)) {
        str_t name, arg;
        if (molden_section_header(&name, &arg, line)) {
            break;
        }
        md_buffered_reader_skip_line(reader);
        str_t s = str_trim(line);
        if (str_empty(s)) {
            continue;
        }
        str_t tok[2];
        str_t scan = s;
        if (extract_tokens(tok, ARRAY_SIZE(tok), &scan) >= 1) {
            md_array_push(*out, molden_parse_float(tok[0]), alloc);
        }
    }
}

static void molden_skip_section(md_buffered_reader_t* reader) {
    str_t line;
    while (md_buffered_reader_peek_line(&line, reader)) {
        str_t name, arg;
        if (molden_section_header(&name, &arg, line)) {
            return;
        }
        md_buffered_reader_skip_line(reader);
    }
}

static bool molden_parse(molden_t* molden, str_t str) {
    md_buffered_reader_t reader = md_buffered_reader_from_str(str);

    // [MO] can only be read once the basis is known, and the format does not promise an order, so
    // the section positions are noted on the first pass and [MO] is read on a second.
    str_t mo_section = {0};

    bool seen_header = false;
    str_t line;
    while (md_buffered_reader_extract_line(&line, &reader)) {
        str_t name, arg;
        if (!molden_section_header(&name, &arg, line)) {
            continue;
        }

        if (str_eq_cstr_ignore_case(name, "molden format")) {
            seen_header = true;
        } else if (str_eq_cstr_ignore_case(name, "title")) {
            str_t title_line;
            if (md_buffered_reader_peek_line(&title_line, &reader)) {
                str_t n, a;
                if (!molden_section_header(&n, &a, title_line)) {
                    molden->title = str_copy(str_trim(title_line), molden->alloc);
                }
            }
            molden_skip_section(&reader);
        } else if (str_eq_cstr_ignore_case(name, "program")) {
            molden->program = str_copy(arg, molden->alloc);
            molden_skip_section(&reader);
        } else if (str_eq_cstr_ignore_case(name, "atoms")) {
            if (!molden_parse_atoms(molden, &reader, arg)) {
                MD_LOG_ERROR("MOLDEN: [Atoms] holds no atoms");
                return false;
            }
        } else if (str_eq_cstr_ignore_case(name, "gto")) {
            if (!molden_parse_gto(molden, &reader)) {
                return false;
            }
        } else if (str_eq_cstr_ignore_case(name, "5d") || str_eq_cstr_ignore_case(name, "5d7f")) {
            molden->pure[2] = true;
            molden->pure[3] = true;
        } else if (str_eq_cstr_ignore_case(name, "5d10f")) {
            molden->pure[2] = true;
            molden->pure[3] = false;
        } else if (str_eq_cstr_ignore_case(name, "7f")) {
            molden->pure[3] = true;
        } else if (str_eq_cstr_ignore_case(name, "9g")) {
            molden->pure[4] = true;
        } else if (str_eq_cstr_ignore_case(name, "6d")) {
            molden->pure[2] = false;
        } else if (str_eq_cstr_ignore_case(name, "10f")) {
            molden->pure[3] = false;
        } else if (str_eq_cstr_ignore_case(name, "15g")) {
            molden->pure[4] = false;
        } else if (str_eq_cstr_ignore_case(name, "freq")) {
            molden_parse_number_series(&molden->frequency, &reader, molden->alloc);
        } else if (str_eq_cstr_ignore_case(name, "int")) {
            molden_parse_number_series(&molden->ir_intensity, &reader, molden->alloc);
        } else if (str_eq_cstr_ignore_case(name, "fr-norm-coord")) {
            molden_parse_normal_modes(molden, &reader);
        } else if (str_eq_cstr_ignore_case(name, "mo")) {
            // Remember where it starts; the basis may not have been read yet.
            mo_section = (str_t){ reader.str.ptr, reader.str.len };
            molden_skip_section(&reader);
        } else {
            molden_skip_section(&reader);
        }
    }

    if (!seen_header) {
        MD_LOG_INFO("MOLDEN: file has no [Molden Format] header; reading it anyway");
    }
    if (md_array_size(molden->atomic_number) == 0) {
        MD_LOG_ERROR("MOLDEN: file contains no atoms");
        return false;
    }

    if (!str_empty(mo_section)) {
        if (!molden->has_gto) {
            MD_LOG_ERROR("MOLDEN: file has [MO] but no [GTO]");
            return false;
        }
        md_buffered_reader_t mo_reader = md_buffered_reader_from_str(mo_section);
        if (!molden_parse_mo(molden, &mo_reader, molden_num_file_ao(molden))) {
            return false;
        }
    }

    return true;
}

// ---------------------------------------------------------------------------
// Basis and coefficients, in md_gto's convention
// ---------------------------------------------------------------------------

// Scales the contraction so that the AO md_gto evaluates has unit norm: one factor per primitive,
// then one per shell. See md_qm.h - this is where a file's convention becomes this library's, and
// doing it per shell is what keeps a mixed spherical/Cartesian file (a [5D10F] one, say) correct.
static bool molden_normalise_basis(molden_t* molden) {
    for (size_t s = 0; s < md_array_size(molden->shell); ++s) {
        const molden_shell_t* shell = &molden->shell[s];
        double* a = molden->alpha + shell->primitive_offset;
        double* c = molden->coeff + shell->primitive_offset;

        for (uint32_t p = 0; p < shell->num_primitives; ++p) {
            if (!(a[p] > 0.0)) {
                MD_LOG_ERROR("MOLDEN: shell %zu has a non positive exponent", s);
                return false;
            }
            c[p] *= md_qm_primitive_norm_factor(shell->l, a[p]);
        }

        const double norm = md_qm_shell_norm_factor(shell->l, a, c, shell->num_primitives, molden_shell_is_pure(molden, shell->l));
        if (!(norm > 0.0)) {
            MD_LOG_ERROR("MOLDEN: shell %zu has a contraction that does not normalise", s);
            return false;
        }
        for (uint32_t p = 0; p < shell->num_primitives; ++p) {
            c[p] /= norm;
        }
    }
    return true;
}

static bool molden_gto_basis_extract(md_gto_basis_t* out, const molden_t* molden, md_allocator_i* alloc) {
    MEMSET(out, 0, sizeof(*out));
    const size_t num_shells = md_array_size(molden->shell);
    if (num_shells == 0) {
        return false;
    }

    md_array_resize(out->shells, num_shells, alloc);
    md_array_resize(out->alpha,  md_array_size(molden->alpha), alloc);
    md_array_resize(out->coeff,  md_array_size(molden->coeff), alloc);

    for (size_t s = 0; s < num_shells; ++s) {
        out->shells[s] = (md_gto_shell_t){
            .atom_idx         = molden->shell[s].atom_idx,
            .primitive_offset = molden->shell[s].primitive_offset,
            .num_primitives   = molden->shell[s].num_primitives,
            .l                = molden->shell[s].l,
        };
    }
    for (size_t p = 0; p < md_array_size(molden->alpha); ++p) {
        out->alpha[p] = (float)molden->alpha[p];
        out->coeff[p] = (float)molden->coeff[p];
    }
    out->num_shells     = (uint32_t)num_shells;
    out->num_primitives = (uint32_t)md_array_size(molden->alpha);
    return true;
}

// One orbital's coefficients, from the file's AO order and convention into md_gto's Cartesian one.
//
// Per shell, because the two conventions can be mixed within one file and because a spherical shell
// changes LENGTH on the way through: (2l+1) coefficients in, (l+1)(l+2)/2 out.
static bool molden_coefficients_to_cartesian(double* dst, size_t dst_cap, const molden_t* molden, const double* src, size_t src_count) {
    size_t si = 0, ci = 0;

    for (size_t s = 0; s < md_array_size(molden->shell); ++s) {
        const uint32_t l        = molden->shell[s].l;
        const uint32_t num_cart = md_gto_num_cart_ao(l);

        if (ci + num_cart > dst_cap) {
            return false;
        }

        if (molden_shell_is_pure(molden, l)) {
            const uint32_t num_sph = md_gto_num_sph_ao(l);
            if (si + num_sph > src_count) {
                return false;
            }
            double sph[2 * MD_GTO_MAX_ANGULAR_MOMENTUM + 1] = {0};
            for (uint32_t k = 0; k < num_sph; ++k) {
                sph[molden_sph_to_md_index(l, k)] = src[si + k];
            }
            // A one shell basis: md_gto_sph_to_cart_vector reads only 'l' off it, and going through
            // the public entry point is what keeps the expansion tables in one place.
            md_gto_shell_t shell = { .atom_idx = 0, .primitive_offset = 0, .num_primitives = 1, .l = l };
            float alpha = 1.0f, coeff = 1.0f;
            md_gto_basis_t one = { .num_shells = 1, .num_primitives = 1, .shells = &shell, .alpha = &alpha, .coeff = &coeff };
            if (md_gto_sph_to_cart_vector(dst + ci, sph, &one) != num_cart) {
                return false;
            }
            si += num_sph;
        } else {
            if (si + num_cart > src_count) {
                return false;
            }
            const molden_lmn_t* order = molden_cart_order(l);
            if (!order) {
                return false;
            }
            for (uint32_t k = 0; k < num_cart; ++k) {
                dst[ci + k] = 0.0;
            }
            for (uint32_t k = 0; k < num_cart; ++k) {
                const uint32_t c = molden_cart_index(l, order[k]);
                if (c >= num_cart) {
                    return false;
                }
                dst[ci + c] = src[si + k] * md_qm_cart_coeff_factor(l, c);
            }
            si += num_cart;
        }
        ci += num_cart;
    }

    return si == src_count && ci == dst_cap;
}

// ---------------------------------------------------------------------------
// Publishing
// ---------------------------------------------------------------------------

static void molden_publish_str(md_system_t* sys, str_t path, str_t label, str_t value) {
    if (str_empty(value)) {
        return;
    }
    md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 1, .shape = { 1 },
    };
    md_attributes_replace(&sys->attributes, &(md_attribute_desc_t){
        .path = path, .format = format, .unit = md_unit_none(), .label = label,
        .data = &value, .byte_size = sizeof(str_t),
    });
}

// A second name for an attribute already published: one datum, two paths, no copy.
static void molden_alias(md_system_t* sys, md_attribute_id_t target, str_t path) {
    if (target == MD_ATTRIBUTE_INVALID) {
        return;
    }
    const md_attribute_t* existing = md_attributes_find(&sys->attributes, path);
    if (existing) {
        md_attributes_remove(&sys->attributes, existing->id);
    }
    md_attributes_alias(&sys->attributes, target, path, (str_t){0}, (str_t){0});
}

static md_attribute_id_t molden_publish_series(md_system_t* sys, str_t path, str_t label, md_unit_t unit, const double* values, size_t count) {
    if (!values || count == 0) {
        return MD_ATTRIBUTE_INVALID;
    }
    md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)count },
    };
    return md_attributes_replace(&sys->attributes, &(md_attribute_desc_t){
        .path = path, .format = format, .unit = unit, .label = label,
        .data = values, .byte_size = count * sizeof(double),
    });
}

// "6D10F15G", "5D7F9G" and the six other spellings the markers can produce. Published because a
// consumer cannot recover it afterwards: the basis in the table is Cartesian whichever way the file
// stated it, which is the point of converting at load time.
static str_t molden_ao_convention(const molden_t* molden, char* buf, size_t cap) {
    int len = snprintf(buf, cap, "%s%s%s",
                       molden->pure[2] ? "5D"  : "6D",
                       molden->pure[3] ? "7F"  : "10F",
                       molden->pure[4] ? "9G"  : "15G");
    if (len <= 0 || (size_t)len >= cap) {
        return (str_t){0};
    }
    return str_from_cstrn(buf, (size_t)len);
}

static bool molden_publish_orbitals(md_system_t* sys, const molden_t* molden, md_allocator_i* temp) {
    const size_t num_orb = md_array_size(molden->orbital);
    if (num_orb == 0) {
        return true;
    }

    const size_t num_file_ao = molden_num_file_ao(molden);
    size_t num_cart_ao = 0;
    for (size_t s = 0; s < md_array_size(molden->shell); ++s) {
        num_cart_ao += md_gto_num_cart_ao(molden->shell[s].l);
    }

    // RESTRICTED or not. Molden tags each orbital Spin= Alpha or Beta, and a restricted calculation
    // writes only Alpha with Occup= 2 - the occupation over BOTH spins. The per channel convention
    // the orbital/ tree uses (and md_vlx.h documents) wants half of that in each channel, with beta
    // a second name for alpha, which is also what makes the total and difference densities come out
    // right. An occupation above 1 with no beta orbital in the file is the only evidence there is,
    // and it is conclusive: an unrestricted file states at most 1 per orbital.
    bool has_beta = false;
    double max_occupation = 0.0;
    for (size_t i = 0; i < num_orb; ++i) {
        has_beta = has_beta || molden->orbital[i].beta;
        max_occupation = MAX(max_occupation, molden->orbital[i].occupation);
    }
    const bool restricted = !has_beta && max_occupation > 1.0;

    md_attribute_id_t alpha_ids[4] = {0};

    for (int spin = 0; spin < 2; ++spin) {
        const bool beta = (spin == 1);

        if (beta && restricted) {
            molden_alias(sys, alpha_ids[0], STR_LIT("orbital/beta/energy"));
            molden_alias(sys, alpha_ids[1], STR_LIT("orbital/beta/occupation"));
            molden_alias(sys, alpha_ids[2], STR_LIT("orbital/beta/symmetry"));
            molden_alias(sys, alpha_ids[3], STR_LIT("orbital/beta/coefficient"));
            break;
        }

        size_t count = 0;
        for (size_t i = 0; i < num_orb; ++i) {
            count += (molden->orbital[i].beta == beta) ? 1 : 0;
        }
        if (count == 0) {
            continue;
        }

        double* energy = md_alloc(temp, sizeof(double) * count);
        double* occ    = md_alloc(temp, sizeof(double) * count);
        str_t*  sym    = md_alloc(temp, sizeof(str_t)  * count);
        double* coeff  = md_alloc(temp, sizeof(double) * count * num_cart_ao);
        if (!energy || !occ || !sym || !coeff) {
            MD_LOG_ERROR("MOLDEN: failed to allocate scratch for %zu orbitals over %zu atomic orbitals", count, num_cart_ao);
            return false;
        }

        bool have_symmetry = false;
        size_t m = 0;
        for (size_t i = 0; i < num_orb; ++i) {
            const molden_orbital_t* o = &molden->orbital[i];
            if (o->beta != beta) {
                continue;
            }
            energy[m] = o->energy;
            occ[m]    = restricted ? 0.5 * o->occupation : o->occupation;
            sym[m]    = o->symmetry;
            have_symmetry = have_symmetry || !str_empty(o->symmetry);
            if (!molden_coefficients_to_cartesian(coeff + m * num_cart_ao, num_cart_ao, molden, o->coefficient, num_file_ao)) {
                MD_LOG_ERROR("MOLDEN: orbital %zu does not span the basis", i);
                return false;
            }
            m += 1;
        }

        const str_t energy_path = beta ? STR_LIT("orbital/beta/energy")      : STR_LIT("orbital/alpha/energy");
        const str_t occ_path    = beta ? STR_LIT("orbital/beta/occupation")  : STR_LIT("orbital/alpha/occupation");
        const str_t sym_path    = beta ? STR_LIT("orbital/beta/symmetry")    : STR_LIT("orbital/alpha/symmetry");
        const str_t coeff_path  = beta ? STR_LIT("orbital/beta/coefficient") : STR_LIT("orbital/alpha/coefficient");

        const md_attribute_id_t ener_id = molden_publish_series(sys, energy_path, STR_LIT("Energy"),     md_unit_hartree(), energy, count);
        const md_attribute_id_t occ_id  = molden_publish_series(sys, occ_path,    STR_LIT("Occupation"), md_unit_none(),    occ,    count);
        md_attribute_id_t sym_id = MD_ATTRIBUTE_INVALID;

        if (have_symmetry) {
            md_attribute_format_t sym_format = {
                .type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 1, .shape = { (uint32_t)count },
            };
            sym_id = md_attributes_replace(&sys->attributes, &(md_attribute_desc_t){
                .path = sym_path, .format = sym_format, .unit = md_unit_none(), .label = STR_LIT("Symmetry"),
                .data = sym, .byte_size = count * sizeof(str_t),
            });
        }

        // f64 and not f32: md_gto takes AO coefficients as double to keep the QM code's precision at
        // the boundary, and there is no point publishing them already narrowed.
        md_attribute_format_t coeff_format = {
            .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2,
            .shape = { (uint32_t)count, (uint32_t)num_cart_ao },
        };
        const md_attribute_id_t coeff_id = md_attributes_replace(&sys->attributes, &(md_attribute_desc_t){
            .path = coeff_path, .format = coeff_format, .unit = md_unit_none(),
            .label = beta ? STR_LIT("Beta Coefficient") : STR_LIT("Alpha Coefficient"),
            .data = coeff, .byte_size = count * num_cart_ao * sizeof(double),
        });

        if (!beta) {
            alpha_ids[0] = ener_id;
            alpha_ids[1] = occ_id;
            alpha_ids[2] = sym_id;
            alpha_ids[3] = coeff_id;
        }
    }

    return true;
}

static void molden_publish_vibrations(md_system_t* sys, const molden_t* molden) {
    const size_t num_atoms = md_array_size(molden->atomic_number);
    const size_t num_freq  = md_array_size(molden->frequency);

    molden_publish_series(sys, STR_LIT("molden/vib/frequency"), STR_LIT("Frequency"), molden_unit_wavenumber(), molden->frequency, num_freq);
    molden_publish_series(sys, STR_LIT("molden/vib/ir_intensity"), STR_LIT("IR Intensity"), molden_unit_km_per_mol(),
                          molden->ir_intensity, md_array_size(molden->ir_intensity));

    if (molden->num_modes == 0 || num_atoms == 0) {
        return;
    }
    if (md_array_size(molden->normal_mode) != molden->num_modes * num_atoms) {
        MD_LOG_ERROR("MOLDEN: [FR-NORM-COORD] holds %zu displacements for %zu modes over %zu atoms",
                     md_array_size(molden->normal_mode), molden->num_modes, num_atoms);
        return;
    }

    md_attribute_format_t format = {
        .type = MD_ATTRIBUTE_TYPE_F64, .components = 3, .rank = 2,
        .shape = { (uint32_t)molden->num_modes, (uint32_t)num_atoms },
    };
    md_attributes_replace(&sys->attributes, &(md_attribute_desc_t){
        .path = STR_LIT("qm/atom/normal_mode"), .format = format, .unit = md_unit_none(),
        .label = STR_LIT("Normal Mode"),
        .data = molden->normal_mode, .byte_size = molden->num_modes * num_atoms * 3 * sizeof(double),
    });
}

// Builds the system's atoms and state, then publishes. It has to happen in this order:
// md_system_reset() clears the attribute table, so a system built after the blocks were published
// would throw away everything they put there.
static bool molden_system_begin(md_system_t* sys, md_system_state_t* state, const molden_t* molden) {
    const size_t num_atoms = md_array_size(molden->atomic_number);

    if (!sys->alloc) {
        MD_LOG_ERROR("System allocator not set");
        return false;
    }
    if (!state || !state->alloc) {
        MD_LOG_ERROR("State allocator not set");
        return false;
    }

    md_system_reset(sys);
    md_system_state_init(state, num_atoms);

    const size_t capacity = ROUND_UP(num_atoms, 16);
    md_array_resize(sys->atom.type_idx, capacity, sys->alloc);
    md_array_resize(sys->atom.flags,    capacity, sys->alloc);
    MEMSET(sys->atom.type_idx, 0, md_array_bytes(sys->atom.type_idx));
    MEMSET(sys->atom.flags,    0, md_array_bytes(sys->atom.flags));

    md_atom_type_find_or_add(&sys->atom.type, STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, sys->alloc);

    for (size_t i = 0; i < num_atoms; ++i) {
        state->x[i] = (float)molden->coord[i].x;
        state->y[i] = (float)molden->coord[i].y;
        state->z[i] = (float)molden->coord[i].z;

        const md_atomic_number_t z = molden->atomic_number[i];
        sys->atom.type_idx[i] = md_atom_type_find_or_add(&sys->atom.type, md_atomic_number_symbol(z), z,
                                                         md_atomic_number_mass(z), md_atomic_number_vdw_radius(z),
                                                         md_atomic_number_cpk_color(z), 0, sys->alloc);
    }

    sys->atom.count  = num_atoms;
    state->num_atoms = num_atoms;
    return true;
}

static bool molden_publish(md_system_t* sys, const molden_t* molden, md_allocator_i* temp) {
    if (!sys->attributes.alloc) {
        MD_LOG_ERROR("Attribute table allocator not set; the system has not been initialised");
        return false;
    }

    molden_publish_str(sys, STR_LIT("molden/title"),   STR_LIT("Title"),   molden->title);
    molden_publish_str(sys, STR_LIT("molden/program"), STR_LIT("Program"), molden->program);

    md_qm_publish_atoms(sys, molden->atomic_number, molden->coord, md_array_size(molden->atomic_number));

    if (molden->has_gto) {
        char buf[32];
        molden_publish_str(sys, STR_LIT("molden/ao_convention"), STR_LIT("AO Convention"), molden_ao_convention(molden, buf, sizeof(buf)));

        md_gto_basis_t basis = {0};
        if (molden_gto_basis_extract(&basis, molden, temp)) {
            md_qm_publish_basis(sys, &basis);
        }
        if (!molden_publish_orbitals(sys, molden, temp)) {
            return false;
        }
        md_qm_publish_overlap(sys);
        md_qm_publish_orbital_densities(sys);
    }

    molden_publish_vibrations(sys, molden);
    return true;
}

// ---------------------------------------------------------------------------
// Entry points
// ---------------------------------------------------------------------------

bool md_molden_system_init_from_str(md_system_t* sys, md_system_state_t* state, str_t str) {
    ASSERT(sys);

    md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp_scope);

    molden_t molden = { .alloc = temp_arena };
    bool result = molden_parse(&molden, str)
               && (!molden.has_gto || molden_normalise_basis(&molden))
               && molden_system_begin(sys, state, &molden)
               && molden_publish(sys, &molden, temp_arena);

    md_temp_end(temp_scope);
    return result;
}

bool md_molden_system_init_from_file(md_system_t* sys, md_system_state_t* state, str_t filename) {
    ASSERT(sys);

    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("MOLDEN: could not open file '" STR_FMT "'", STR_ARG(filename));
        return false;
    }

    const size_t size = (size_t)md_file_size(file);
    md_allocator_i* alloc = sys->alloc ? sys->alloc : md_get_heap_allocator();
    char* buf = (char*)md_alloc(alloc, size + 1);
    bool result = false;
    if (!buf) {
        MD_LOG_ERROR("MOLDEN: failed to allocate %zu bytes for '" STR_FMT "'", size, STR_ARG(filename));
    } else if (md_file_read(file, buf, size) != size) {
        MD_LOG_ERROR("MOLDEN: failed to read '" STR_FMT "'", STR_ARG(filename));
    } else {
        buf[size] = '\0';
        result = md_molden_system_init_from_str(sys, state, (str_t){ buf, size });
    }

    if (buf) {
        md_free(alloc, buf, size + 1);
    }
    md_file_close(&file);
    return result;
}

bool md_molden_file_is_molden(str_t filename) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        return false;
    }

    char buf[512];
    const size_t n = md_file_read(file, buf, sizeof(buf));
    md_file_close(&file);

    str_t str = { buf, n };
    str_t line;
    while (str_extract_line(&line, &str)) {
        line = str_trim(line);
        if (str_empty(line)) {
            continue;
        }
        str_t name, arg;
        return molden_section_header(&name, &arg, line) && str_eq_cstr_ignore_case(name, "molden format");
    }
    return false;
}
