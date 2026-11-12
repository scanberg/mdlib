#include <md_molden.h>

#include <md_system.h>
#include <md_types.h>
#include <md_gto.h>
#include <md_qm.h>
#include <md_element.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_file.h>
#include <core/md_parse.h>
#include <core/md_reader.h>
#include <core/md_str.h>
#include <core/md_format.h>
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
    md_str_t   symmetry;    // interned into the parse arena
    bool    beta;
    double* coefficient; // [num_ao_file], in the FILE's AO order and convention
} molden_orbital_t;

typedef struct molden_t {
    md_allocator_i* alloc;

    md_str_t title;
    md_str_t program;

    // Angular convention per angular momentum: true when the file stores the pure/spherical
    // (2l+1) set. Cartesian is the format's default and is what an unmarked file means.
    bool pure[MD_GTO_MAX_ANGULAR_MOMENTUM + 1];

    md_array(uint8_t)          atomic_number;
    md_array(md_f64x3_t)          coord;       // Angstrom
    md_array(md_str_t)            atom_label;

    md_array(molden_shell_t)   shell;
    md_array(double)           alpha;
    md_array(double)           coeff;       // as read; normalised in molden_normalise_basis

    md_array(molden_orbital_t) orbital;

    // [FR-COORD]: the geometry of a vibrational file, in bohr, which is all a file written for its
    // vibrations alone (CP2K, Psi4, Q-Chem) has: it has no [Atoms]
    md_array(uint8_t)          fr_atomic_number;
    md_array(md_f64x3_t)          fr_coord;      // Angstrom, converted on read
    md_array(md_str_t)            fr_label;

    md_array(double)           frequency;     // cm^-1
    md_array(double)           ir_intensity;  // km/mol
    md_array(md_f64x3_t)          normal_mode;   // [num_modes * num_atoms]
    size_t                     num_modes;

    bool has_gto;

    // How the file normalised its Cartesian functions of l >= 2 (molden_cartesian_convention): each on its
    // own (false; the format's convention, and the default) or all of a shell with its (l,0,0) monomial's
    // constant (true; Turbomole, CFOUR 2.1, PSI4 <= 1.3.2)
    bool cart_shared;
} molden_t;

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
static double molden_parse_float(md_str_t tok) {
    char buf[64];
    size_t n = md_str_copy_to_char_buf(buf, sizeof(buf), tok);
    for (size_t i = 0; i < n; ++i) {
        if (buf[i] == 'D' || buf[i] == 'd') buf[i] = 'e';
    }
    return parse_float(md_str_from_cstrn(buf, n));
}

// "[Atoms] (Angs)" -> section "atoms", argument "(Angs)". Returns false for anything that is not a
// section header, which is every other line in the file.
static bool molden_section_header(md_str_t* out_name, md_str_t* out_arg, md_str_t line) {
    line = md_str_trim(line);
    if (line.len < 2 || line.ptr[0] != '[') {
        return false;
    }
    md_str_t name, arg;
    if (!md_str_cut(&name, &arg, line, md_str_find_char(line, ']'))) {
        return false;
    }
    *out_name = md_str_trim(md_str_substr(name, 1, SIZE_MAX));   // after the '['
    *out_arg  = md_str_trim(arg);
    return true;
}

static int molden_shell_type_to_l(md_str_t tok) {
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

static bool molden_parse_atoms(molden_t* molden, md_reader_i* reader, md_str_t arg) {
    // "[Atoms] Angs", "[Atoms] (AU)", "[ATOMS] AU" - the parentheses are decoration. Angstrom is
    // the default, which is what an argument-less header means.
    md_str_t a = arg;
    if (a.len >= 2 && a.ptr[0] == '(' && a.ptr[a.len - 1] == ')') {
        a = md_str_trim(md_str_substr(a, 1, a.len - 2));
    }
    const bool bohr = md_str_eq_cstr_ignore_case(a, "au") || md_str_eq_cstr_ignore_case(a, "bohr");

    md_str_t line;
    while (md_reader_peek_line(reader, &line)) {
        md_str_t name, sec_arg;
        if (molden_section_header(&name, &sec_arg, line)) {
            break;
        }
        md_reader_skip_line(reader);

        md_str_t tok[8];
        md_str_t rest = line;
        const size_t num_tok = md_parse_next_tokens(tok, MD_COUNTOF(tok), &rest);
        if (num_tok < 6) {
            continue;
        }
        // symbol index atomic_number x y z
        md_atomic_number_t z = (md_atomic_number_t)parse_int(tok[2]);
        if (z == 0) {
            z = md_atomic_number_from_symbol(tok[0], true);
        }
        md_f64x3_t xyz = {
            molden_parse_float(tok[3]),
            molden_parse_float(tok[4]),
            molden_parse_float(tok[5]),
        };
        if (bohr) {
            xyz = md_f64x3_mul1(xyz, MOLDEN_BOHR_TO_ANGSTROM);
        }
        md_array_push(molden->atomic_number, (uint8_t)z, molden->alloc);
        md_array_push(molden->coord, xyz, molden->alloc);
        md_array_push(molden->atom_label, md_alloc_str(molden->alloc, tok[0]), molden->alloc);
    }
    return md_array_size(molden->atomic_number) > 0;
}

// [GTO] is per atom: an "<atom index> 0" line, then shells until a blank line.
//
// 'sp' is a shell with TWO coefficient columns sharing one set of exponents; it becomes an s shell
// and a p shell in that order, which is also the order its four AOs appear in, so nothing downstream
// has to know it existed.
static bool molden_parse_gto(molden_t* molden, md_reader_i* reader) {
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

    md_str_t line;
    while (md_reader_peek_line(reader, &line)) {
        md_str_t name, arg;
        if (molden_section_header(&name, &arg, line)) {
            break;
        }
        md_reader_skip_line(reader);

        md_str_t rest = md_str_trim(line);
        if (md_str_empty(rest)) {
            have_atom = false;
            continue;
        }

        md_str_t tok[4];
        md_str_t scan = rest;
        const size_t num_tok = md_parse_next_tokens(tok, MD_COUNTOF(tok), &scan);
        if (num_tok == 0) {
            continue;
        }

        const int  l     = molden_shell_type_to_l(tok[0]);
        const bool is_sp = md_str_eq_cstr_ignore_case(tok[0], "sp");

        // A line is either a shell - which starts with its type letter - or the "<atom index> 0"
        // header of the next atom's block. Nothing else is legal here, and an unknown letter is a
        // shell type this library cannot evaluate rather than something to skip past quietly.
        if (l < 0 && !is_sp) {
            if (!is_int(tok[0])) {
                MD_LOG_ERROR("MOLDEN: [GTO] shell type '" MD_STR_FMT "' is not one this library evaluates (s, p, d, f, g and sp)", MD_STR_ARG(tok[0]));
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
            MD_LOG_ERROR("MOLDEN: malformed shell line '" MD_STR_FMT "'", MD_STR_ARG(rest));
            return false;
        }
        const size_t num_prim = (size_t)parse_int(tok[1]);
        if (num_prim == 0 || num_prim > MOLDEN_MAX_PRIMITIVES) {
            MD_LOG_ERROR("MOLDEN: shell with %zu primitives, which is outside 1..%d", num_prim, MOLDEN_MAX_PRIMITIVES);
            return false;
        }

        for (size_t p = 0; p < num_prim; ++p) {
            md_str_t prim_line;
            if (!md_reader_extract_line(reader, &prim_line)) {
                MD_LOG_ERROR("MOLDEN: [GTO] ends inside a contraction");
                return false;
            }
            md_str_t ptok[4];
            md_str_t pscan = prim_line;
            const size_t n = md_parse_next_tokens(ptok, MD_COUNTOF(ptok), &pscan);
            if (n < 2) {
                MD_LOG_ERROR("MOLDEN: malformed primitive line '" MD_STR_FMT "'", MD_STR_ARG(md_str_trim(prim_line)));
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
static bool molden_parse_mo(molden_t* molden, md_reader_i* reader, size_t num_ao) {
    if (num_ao == 0) {
        MD_LOG_ERROR("MOLDEN: [MO] appears before a basis is known");
        return false;
    }

    molden_orbital_t orb = {0};
    bool have_orb = false;
    size_t num_coeff = 0;

    md_str_t line;
    while (md_reader_peek_line(reader, &line)) {
        md_str_t name, arg;
        if (molden_section_header(&name, &arg, line)) {
            break;
        }
        md_reader_skip_line(reader);

        md_str_t s = md_str_trim(line);
        if (md_str_empty(s)) {
            continue;
        }

        size_t eq;
        const bool is_key = ((eq = md_str_find_char(s, '=')) != MD_STR_NPOS) && eq <= 8;
        if (is_key) {
            if (have_orb && num_coeff > 0) {
                md_array_push(molden->orbital, orb, molden->alloc);
                md_memset(&orb, 0, sizeof(orb));
                have_orb  = false;
                num_coeff = 0;
            }
            if (!have_orb) {
                orb = (molden_orbital_t){0};
                orb.coefficient = (double*)md_alloc(molden->alloc, sizeof(double) * num_ao);
                md_memset(orb.coefficient, 0, sizeof(double) * num_ao);
                have_orb = true;
            }

            md_str_t key = md_str_trim(md_str_substr(s, 0, eq));
            md_str_t val = md_str_trim(md_str_substr(s, eq + 1, SIZE_MAX));
            if (md_str_eq_cstr_ignore_case(key, "ene")) {
                orb.energy = molden_parse_float(val);
            } else if (md_str_eq_cstr_ignore_case(key, "spin")) {
                orb.beta = md_str_eq_cstr_ignore_case(val, "beta");
            } else if (md_str_eq_cstr_ignore_case(key, "occup")) {
                orb.occupation = molden_parse_float(val);
            } else if (md_str_eq_cstr_ignore_case(key, "sym")) {
                orb.symmetry = md_alloc_str(molden->alloc, val);
            }
            continue;
        }

        if (!have_orb) {
            continue;
        }

        md_str_t tok[2];
        md_str_t scan = s;
        if (md_parse_next_tokens(tok, MD_COUNTOF(tok), &scan) < 2) {
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

// [FR-COORD] is "<symbol> x y z" per atom, in bohr.
static void molden_parse_fr_coord(molden_t* molden, md_reader_i* reader) {
    md_str_t line;
    while (md_reader_peek_line(reader, &line)) {
        md_str_t name, arg;
        if (molden_section_header(&name, &arg, line)) {
            break;
        }
        md_reader_skip_line(reader);

        md_str_t tok[5];
        md_str_t scan = line;
        if (md_parse_next_tokens(tok, MD_COUNTOF(tok), &scan) < 4) {
            continue;
        }
        const md_atomic_number_t z = md_atomic_number_from_symbol(tok[0], true);
        md_f64x3_t xyz = {
            molden_parse_float(tok[1]) * MOLDEN_BOHR_TO_ANGSTROM,
            molden_parse_float(tok[2]) * MOLDEN_BOHR_TO_ANGSTROM,
            molden_parse_float(tok[3]) * MOLDEN_BOHR_TO_ANGSTROM,
        };
        md_array_push(molden->fr_atomic_number, (uint8_t)z, molden->alloc);
        md_array_push(molden->fr_coord, xyz, molden->alloc);
        md_array_push(molden->fr_label, md_alloc_str(molden->alloc, tok[0]), molden->alloc);
    }
}

// [FR-NORM-COORD] is "vibration <n>" followed by one displacement per atom, in bohr.
static void molden_parse_normal_modes(molden_t* molden, md_reader_i* reader) {
    const size_t num_atoms = md_array_size(molden->atomic_number);
    md_str_t line;
    while (md_reader_peek_line(reader, &line)) {
        md_str_t name, arg;
        if (molden_section_header(&name, &arg, line)) {
            break;
        }
        md_reader_skip_line(reader);

        md_str_t s = md_str_trim(line);
        if (md_str_empty(s)) {
            continue;
        }
        if (!md_str_begins_with_ignore_case(s, MD_STR_LIT("vibration"))) {
            continue;
        }
        for (size_t a = 0; a < num_atoms; ++a) {
            md_str_t disp_line;
            if (!md_reader_extract_line(reader, &disp_line)) {
                return;
            }
            md_str_t tok[4];
            md_str_t scan = disp_line;
            md_f64x3_t d = {0};
            if (md_parse_next_tokens(tok, MD_COUNTOF(tok), &scan) >= 3) {
                d.x = molden_parse_float(tok[0]);
                d.y = molden_parse_float(tok[1]);
                d.z = molden_parse_float(tok[2]);
            }
            md_array_push(molden->normal_mode, d, molden->alloc);
        }
        molden->num_modes += 1;
    }
}

static void molden_parse_number_series(md_array(double)* out, md_reader_i* reader, md_allocator_i* alloc) {
    md_str_t line;
    while (md_reader_peek_line(reader, &line)) {
        md_str_t name, arg;
        if (molden_section_header(&name, &arg, line)) {
            break;
        }
        md_reader_skip_line(reader);
        md_str_t s = md_str_trim(line);
        if (md_str_empty(s)) {
            continue;
        }
        md_str_t tok[2];
        md_str_t scan = s;
        if (md_parse_next_tokens(tok, MD_COUNTOF(tok), &scan) >= 1) {
            md_array_push(*out, molden_parse_float(tok[0]), alloc);
        }
    }
}

static void molden_skip_section(md_reader_i* reader) {
    md_str_t line;
    while (md_reader_peek_line(reader, &line)) {
        md_str_t name, arg;
        if (molden_section_header(&name, &arg, line)) {
            return;
        }
        md_reader_skip_line(reader);
    }
}

static bool molden_parse(molden_t* molden, md_str_t str) {
    md_reader_i reader = md_reader_from_str(str);

    // [MO] can only be read once the basis is known, and the format does not promise an order, so
    // the section positions are noted on the first pass and [MO] is read on a second.
    md_str_t mo_section = {0};
    md_str_t fr_norm_section = {0};

    bool seen_header = false;
    md_str_t line;
    while (md_reader_extract_line(&reader, &line)) {
        md_str_t name, arg;
        if (!molden_section_header(&name, &arg, line)) {
            continue;
        }

        if (md_str_eq_cstr_ignore_case(name, "molden format")) {
            seen_header = true;
        } else if (md_str_eq_cstr_ignore_case(name, "title")) {
            md_str_t title_line;
            if (md_reader_peek_line(&reader, &title_line)) {
                md_str_t n, a;
                if (!molden_section_header(&n, &a, title_line)) {
                    molden->title = md_alloc_str(molden->alloc, md_str_trim(title_line));
                }
            }
            molden_skip_section(&reader);
        } else if (md_str_eq_cstr_ignore_case(name, "program")) {
            molden->program = md_alloc_str(molden->alloc, arg);
            molden_skip_section(&reader);
        } else if (md_str_eq_cstr_ignore_case(name, "atoms")) {
            if (!molden_parse_atoms(molden, &reader, arg)) {
                MD_LOG_ERROR("MOLDEN: [Atoms] holds no atoms");
                return false;
            }
        } else if (md_str_eq_cstr_ignore_case(name, "gto")) {
            if (!molden_parse_gto(molden, &reader)) {
                return false;
            }
        } else if (md_str_eq_cstr_ignore_case(name, "5d") || md_str_eq_cstr_ignore_case(name, "5d7f")) {
            molden->pure[2] = true;
            molden->pure[3] = true;
        } else if (md_str_eq_cstr_ignore_case(name, "5d10f")) {
            molden->pure[2] = true;
            molden->pure[3] = false;
        } else if (md_str_eq_cstr_ignore_case(name, "7f")) {
            molden->pure[3] = true;
        } else if (md_str_eq_cstr_ignore_case(name, "9g")) {
            molden->pure[4] = true;
        } else if (md_str_eq_cstr_ignore_case(name, "6d")) {
            molden->pure[2] = false;
        } else if (md_str_eq_cstr_ignore_case(name, "10f")) {
            molden->pure[3] = false;
        } else if (md_str_eq_cstr_ignore_case(name, "15g")) {
            molden->pure[4] = false;
        } else if (md_str_eq_cstr_ignore_case(name, "freq")) {
            molden_parse_number_series(&molden->frequency, &reader, molden->alloc);
        } else if (md_str_eq_cstr_ignore_case(name, "int")) {
            molden_parse_number_series(&molden->ir_intensity, &reader, molden->alloc);
        } else if (md_str_eq_cstr_ignore_case(name, "fr-coord")) {
            molden_parse_fr_coord(molden, &reader);
        } else if (md_str_eq_cstr_ignore_case(name, "fr-norm-coord")) {
            // Read once the atom count is known, which [Atoms] or [FR-COORD] may only state further down
            fr_norm_section = (md_str_t){ reader.cur, (size_t)(reader.end - reader.cur) };
            molden_skip_section(&reader);
        } else if (md_str_eq_cstr_ignore_case(name, "mo")) {
            // Remember where it starts; the basis may not have been read yet.
            mo_section = (md_str_t){ reader.cur, (size_t)(reader.end - reader.cur) };
            molden_skip_section(&reader);
        } else {
            molden_skip_section(&reader);
        }
    }

    if (!seen_header) {
        MD_LOG_INFO("MOLDEN: file has no [Molden Format] header; reading it anyway");
    }
    if (md_array_size(molden->atomic_number) == 0 && md_array_size(molden->fr_atomic_number) > 0) {
        // A file written for its vibrations alone: [FR-COORD] is its geometry
        molden->atomic_number = molden->fr_atomic_number;
        molden->coord         = molden->fr_coord;
        molden->atom_label    = molden->fr_label;
    }
    if (md_array_size(molden->atomic_number) == 0) {
        MD_LOG_ERROR("MOLDEN: file contains no atoms");
        return false;
    }

    if (!md_str_empty(fr_norm_section)) {
        md_reader_i fr_reader = md_reader_from_str(fr_norm_section);
        molden_parse_normal_modes(molden, &fr_reader);
    }

    if (!md_str_empty(mo_section)) {
        if (!molden->has_gto) {
            MD_LOG_ERROR("MOLDEN: file has [MO] but no [GTO]");
            return false;
        }
        md_reader_i mo_reader = md_reader_from_str(mo_section);
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
    md_memset(out, 0, sizeof(*out));
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
//
// A Cartesian shell's coefficients are taken as they are when the file normalised each Cartesian function on
// its own, as the format does (Molden itself, PySCF, IOData) and as md_gto does: only the order changes. A
// file written against functions that share the (l,0,0) monomial's normalisation (molden->cart_shared, see
// molden_cartesian_convention) has its coefficients scaled onto md_gto's by md_qm_cart_coeff_factor.
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
                dst[ci + c] = molden->cart_shared ? src[si + k] * md_qm_cart_coeff_factor(l, c) : src[si + k];
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

// "6D10F15G", "5D7F9G" and the six other spellings the markers can produce. Published because a
// consumer cannot recover it afterwards: the basis in the table is Cartesian whichever way the file
// stated it, which is the point of converting at load time.
static md_str_t molden_ao_convention(const molden_t* molden, char* buf, size_t cap) {
    int len = md_snprintf(buf, cap, "%s%s%s",
                       molden->pure[2] ? "5D"  : "6D",
                       molden->pure[3] ? "7F"  : "10F",
                       molden->pure[4] ? "9G"  : "15G");
    if (len <= 0 || (size_t)len >= cap) {
        return (md_str_t){0};
    }
    return md_str_from_cstrn(buf, (size_t)len);
}

static bool molden_has_cartesian_beyond_p(const molden_t* molden) {
    for (size_t s = 0; s < md_array_size(molden->shell); ++s) {
        if (molden->shell[s].l >= 2 && !molden_shell_is_pure(molden, molden->shell[s].l)) return true;
    }
    return false;
}

// HOW THE FILE NORMALISED ITS CARTESIAN FUNCTIONS. The format normalises each Cartesian function on its own -
// Molden itself, PySCF's writer and reader, IOData - and so does md_gto: then a file's coefficients are taken as
// they are. Some programs write them against functions that all share their shell's (l,0,0) normalisation
// instead, so that xy is sqrt(3) smaller than xx: Turbomole, CFOUR 2.1, PSI4 <= 1.3.2 (which IOData corrects on
// reading, as this does). The file does not say which, and the two differ only in the coefficients of the
// Cartesian functions of l >= 2 that are not powers of one coordinate. The orbitals do: a file's orbitals are
// orthonormal, and only in the convention they were written in. So where the file has Cartesian shells of
// l >= 2, the norms <psi|S|psi> of the orbitals with the most weight on those functions are compared under both,
// and the other convention is taken only when it is clearly the one they are normalised in. A file normalised in
// neither keeps the format's, with a warning.
#define MOLDEN_CONVENTION_ORBITALS 16
#define MOLDEN_CONVENTION_MAX_AO   4096    // the overlap is num_ao^2 doubles; past this the format's is assumed

static void molden_cartesian_convention(molden_t* molden, const md_gto_basis_t* basis, md_allocator_i* temp) {
    molden->cart_shared = false;
    const size_t num_ao = md_gto_basis_num_ao(basis);
    const size_t num_orb = md_array_size(molden->orbital);
    if (num_orb == 0 || num_ao == 0) return;

    // the factor each md_gto AO's coefficient takes in the shared convention; 1 where the two agree
    md_temp_scope_t scope = md_temp_begin_avoid(temp);
    double* fac = md_temp_alloc_array(scope, double, num_ao);
    bool any = false;
    size_t ao = 0;
    for (size_t s = 0; s < md_array_size(molden->shell); ++s) {
        const uint32_t l = molden->shell[s].l;
        const bool cart = !molden_shell_is_pure(molden, l);
        for (uint32_t c = 0; c < md_gto_num_cart_ao(l); ++c, ++ao) {
            fac[ao] = (cart && l >= 2) ? md_qm_cart_coeff_factor(l, c) : 1.0;
            any = any || fac[ao] != 1.0;
        }
    }
    if (!any || ao != num_ao) {
        md_temp_end(scope);
        return;
    }
    if (num_ao > MOLDEN_CONVENTION_MAX_AO) {
        MD_LOG_INFO("MOLDEN: %zu atomic orbitals; each Cartesian function is taken as normalised on its own, the format's convention, without checking", num_ao);
        md_temp_end(scope);
        return;
    }

    // every orbital in the format's convention, and the ones with the most weight where the two differ
    const size_t num_file_ao = molden_num_file_ao(molden);
    double* C = md_temp_alloc_array(scope, double, num_orb * num_ao);
    double* S = md_temp_alloc_array(scope, double, num_ao * num_ao);
    md_f64x3_t* xyz = md_temp_alloc_array(scope, md_f64x3_t, md_array_size(molden->coord));
    if (!C || !S || !xyz) {
        md_temp_end(scope);
        return;
    }
    size_t pick[MOLDEN_CONVENTION_ORBITALS];
    double weight[MOLDEN_CONVENTION_ORBITALS];
    size_t num_pick = 0;
    for (size_t i = 0; i < num_orb; ++i) {
        double* c = C + i * num_ao;
        if (!molden_coefficients_to_cartesian(c, num_ao, molden, molden->orbital[i].coefficient, num_file_ao)) {
            md_temp_end(scope);
            return;    // reported where the orbitals are published
        }
        double w = 0.0;
        for (size_t a = 0; a < num_ao; ++a) if (fac[a] != 1.0) w += c[a] * c[a];
        // keep the MOLDEN_CONVENTION_ORBITALS largest
        size_t at = num_pick;
        if (num_pick < MOLDEN_CONVENTION_ORBITALS) num_pick++;
        else if (w <= weight[num_pick - 1]) continue;
        else at = num_pick - 1;
        while (at > 0 && weight[at - 1] < w) { weight[at] = weight[at - 1]; pick[at] = pick[at - 1]; --at; }
        weight[at] = w;
        pick[at] = i;
    }
    for (size_t a = 0; a < md_array_size(molden->coord); ++a) {
        xyz[a] = (md_f64x3_t){ molden->coord[a].x * 1.8897261246257702, molden->coord[a].y * 1.8897261246257702, molden->coord[a].z * 1.8897261246257702 };
    }
    if (num_pick == 0 || weight[0] == 0.0 || md_qm_compute_overlap(S, basis, xyz, md_array_size(molden->coord)) != num_ao) {
        md_temp_end(scope);
        return;
    }

    // the largest |<psi|S|psi> - 1| of those orbitals in either convention
    double* sc = md_temp_alloc_array(scope, double, num_ao);
    double* c2 = md_temp_alloc_array(scope, double, num_ao);
    double err[2] = {0.0, 0.0};
    for (size_t k = 0; k < num_pick; ++k) {
        const double* c = C + pick[k] * num_ao;
        for (int conv = 0; conv < 2; ++conv) {
            for (size_t a = 0; a < num_ao; ++a) c2[a] = conv ? c[a] * fac[a] : c[a];
            double n = 0.0;
            for (size_t a = 0; a < num_ao; ++a) {
                const double* row = S + a * num_ao;
                double v = 0.0;
                for (size_t b = 0; b < num_ao; ++b) v += row[b] * c2[b];
                sc[a] = v;
            }
            for (size_t a = 0; a < num_ao; ++a) n += c2[a] * sc[a];
            err[conv] = fmax(err[conv], fabs(n - 1.0));
        }
    }
    // The other convention only when the orbitals are normalised in it and clearly not in the format's: the
    // basis is float, which leaves norms ~1e-7 from 1 in a file written to many digits
    molden->cart_shared = err[0] > 1.0e-5 && err[1] < 0.25 * err[0];
    const double chosen = molden->cart_shared ? err[1] : err[0];
    if (molden->cart_shared) {
        MD_LOG_INFO("MOLDEN: the Cartesian functions are normalised as a shell, not each on its own (as Turbomole, CFOUR 2.1 and PSI4 <= 1.3.2 write them): "
                    "the orbitals are normalised to %.1e that way, %.1e the format's", err[1], err[0]);
    }
    if (chosen > 1.0e-3) {
        MD_LOG_INFO("MOLDEN: the orbitals are not normalised with the Cartesian functions normalised either way (%.1e each on its own, %.1e as a shell); "
                       "read as the format has them", err[0], err[1]);
    }
    md_temp_end(scope);
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
        max_occupation = MD_MAX(max_occupation, molden->orbital[i].occupation);
    }
    const bool restricted = !has_beta && max_occupation > 1.0;

    // A restricted OPEN shell (ROHF, ROKS) states its singly occupied orbitals with Occup= 1, and those hold an
    // alpha electron, not half of one in each channel: halved, the spin density of a radical would be zero. Taken
    // as such only when every occupation is 0, 1 or 2, which is what such a calculation writes; fractional ones
    // (natural orbitals) are halved as before, which keeps their spin density at zero, as for a singlet.
    bool open_shell = false;
    if (restricted) {
        bool integral = true;
        for (size_t i = 0; i < num_orb; ++i) {
            const double o = molden->orbital[i].occupation;
            integral = integral && (o == 0.0 || o == 1.0 || o == 2.0);
            open_shell = open_shell || o == 1.0;
        }
        open_shell = open_shell && integral;
    }

    md_attribute_id_t alpha_ids[4] = {0};

    for (int spin = 0; spin < 2; ++spin) {
        const bool beta = (spin == 1);

        if (beta && restricted) {
            md_qm_alias(sys, alpha_ids[0], MD_STR_LIT("orbital/beta/energy"));
            md_qm_alias(sys, alpha_ids[2], MD_STR_LIT("orbital/beta/symmetry"));
            md_qm_alias(sys, alpha_ids[3], MD_STR_LIT("orbital/beta/coefficient"));
            if (!open_shell) {
                md_qm_alias(sys, alpha_ids[1], MD_STR_LIT("orbital/beta/occupation"));
            } else {
                // the same orbitals, with what the alpha channel does not hold
                double* occ = md_alloc(temp, sizeof(double) * num_orb);
                if (!occ) {
                    return false;
                }
                for (size_t i = 0; i < num_orb; ++i) {
                    occ[i] = molden->orbital[i].occupation - MD_MIN(molden->orbital[i].occupation, 1.0);
                }
                md_qm_publish_series(sys, MD_STR_LIT("orbital/beta/occupation"), MD_STR_LIT("Occupation"), md_unit_none(), occ, num_orb);
            }
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
        md_str_t*  sym    = md_alloc(temp, sizeof(md_str_t)  * count);
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
            occ[m]    = !restricted ? o->occupation : open_shell ? MD_MIN(o->occupation, 1.0) : 0.5 * o->occupation;
            sym[m]    = o->symmetry;
            have_symmetry = have_symmetry || !md_str_empty(o->symmetry);
            if (!molden_coefficients_to_cartesian(coeff + m * num_cart_ao, num_cart_ao, molden, o->coefficient, num_file_ao)) {
                MD_LOG_ERROR("MOLDEN: orbital %zu does not span the basis", i);
                return false;
            }
            m += 1;
        }

        const md_str_t energy_path = beta ? MD_STR_LIT("orbital/beta/energy")      : MD_STR_LIT("orbital/alpha/energy");
        const md_str_t occ_path    = beta ? MD_STR_LIT("orbital/beta/occupation")  : MD_STR_LIT("orbital/alpha/occupation");
        const md_str_t sym_path    = beta ? MD_STR_LIT("orbital/beta/symmetry")    : MD_STR_LIT("orbital/alpha/symmetry");
        const md_str_t coeff_path  = beta ? MD_STR_LIT("orbital/beta/coefficient") : MD_STR_LIT("orbital/alpha/coefficient");

        const md_attribute_id_t ener_id = md_qm_publish_series(sys, energy_path, MD_STR_LIT("Energy"),     md_unit_hartree(), energy, count);
        const md_attribute_id_t occ_id  = md_qm_publish_series(sys, occ_path,    MD_STR_LIT("Occupation"), md_unit_none(),    occ,    count);
        md_attribute_id_t sym_id = MD_ATTRIBUTE_INVALID;

        if (have_symmetry) {
            sym_id = md_qm_publish_strings(sys, sym_path, MD_STR_LIT("Symmetry"), sym, count);
        }

        // f64 and not f32: md_gto takes AO coefficients as double to keep the QM code's precision at
        // the boundary, and there is no point publishing them already narrowed.
        const md_attribute_id_t coeff_id = md_qm_publish_matrix(sys, coeff_path,
            beta ? MD_STR_LIT("Beta Coefficient") : MD_STR_LIT("Alpha Coefficient"),
            md_unit_none(), coeff, count, num_cart_ao);

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

    md_qm_publish_series(sys, MD_STR_LIT("molden/vib/frequency"), MD_STR_LIT("Frequency"), md_qm_unit_wavenumber(), molden->frequency, num_freq);
    md_qm_publish_series(sys, MD_STR_LIT("molden/vib/ir_intensity"), MD_STR_LIT("IR Intensity"), md_qm_unit_km_per_mol(),
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
    md_qm_publish(sys, MD_STR_LIT("qm/atom/normal_mode"), MD_STR_LIT("Normal Mode"), md_unit_none(),
                  format, molden->normal_mode, molden->num_modes * num_atoms * 3 * sizeof(double));
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

    const size_t capacity = MD_ROUND_UP(num_atoms, 16);
    md_array_resize(sys->atom.type_idx, capacity, sys->alloc);
    md_array_resize(sys->atom.flags,    capacity, sys->alloc);
    md_memset(sys->atom.type_idx, 0, md_array_bytes(sys->atom.type_idx));
    md_memset(sys->atom.flags,    0, md_array_bytes(sys->atom.flags));

    md_particle_type_find_or_add(&sys->atom.type, MD_STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, sys->alloc);

    for (size_t i = 0; i < num_atoms; ++i) {
        state->xyz[i] = md_f32x3_set((float)molden->coord[i].x, (float)molden->coord[i].y, (float)molden->coord[i].z);

        const md_atomic_number_t z = molden->atomic_number[i];
        sys->atom.type_idx[i] = md_particle_type_find_or_add(&sys->atom.type, md_atomic_number_symbol(z), z,
                                                         md_atomic_number_mass(z), md_atomic_number_vdw_radius(z),
                                                         md_atomic_number_cpk_color(z), 0, sys->alloc);
        // A QM calculation alone: every atom is in it (MD_ATOM_FLAG_QM)
        sys->atom.flags[i] |= MD_ATOM_FLAG_QM;
    }

    sys->atom.count  = num_atoms;
    state->num_atoms = num_atoms;
    return true;
}

static bool molden_publish(md_system_t* sys, molden_t* molden, md_allocator_i* temp) {
    if (!sys->attributes.alloc) {
        MD_LOG_ERROR("Attribute table allocator not set; the system has not been initialised");
        return false;
    }

    md_qm_publish_str(sys, MD_STR_LIT("molden/title"),   MD_STR_LIT("Title"),   molden->title);
    md_qm_publish_str(sys, MD_STR_LIT("molden/program"), MD_STR_LIT("Program"), molden->program);

    // Molden states no effective core charges ([Atoms] has one number per atom, read as the
    // element), so the nuclear charges are the atomic numbers
    md_qm_publish_atoms(sys, molden->atomic_number, NULL, molden->coord, md_array_size(molden->atomic_number));

    if (molden->has_gto) {
        char buf[32];
        md_qm_publish_str(sys, MD_STR_LIT("molden/ao_convention"), MD_STR_LIT("AO Convention"), molden_ao_convention(molden, buf, sizeof(buf)));

        md_gto_basis_t basis = {0};
        if (molden_gto_basis_extract(&basis, molden, temp)) {
            md_qm_publish_basis(sys, &basis);
            molden_cartesian_convention(molden, &basis, temp);
            if (molden->cart_shared) {
                md_qm_publish_str(sys, MD_STR_LIT("molden/cartesian_normalisation"), MD_STR_LIT("Cartesian Normalisation"), MD_STR_LIT("shell"));
            } else if (molden_has_cartesian_beyond_p(molden)) {
                md_qm_publish_str(sys, MD_STR_LIT("molden/cartesian_normalisation"), MD_STR_LIT("Cartesian Normalisation"), MD_STR_LIT("function"));
            }
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

bool md_molden_system_init_from_str(md_system_t* sys, md_system_state_t* state, md_str_t str) {
    MD_ASSERT(sys);

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

bool md_molden_system_init_from_file(md_system_t* sys, md_system_state_t* state, md_str_t filename) {
    MD_ASSERT(sys);

    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("MOLDEN: could not open file '" MD_STR_FMT "'", MD_STR_ARG(filename));
        return false;
    }

    const size_t size = (size_t)md_file_size(file);
    md_allocator_i* alloc = sys->alloc ? sys->alloc : md_get_heap_allocator();
    char* buf = (char*)md_alloc(alloc, size + 1);
    bool result = false;
    if (!buf) {
        MD_LOG_ERROR("MOLDEN: failed to allocate %zu bytes for '" MD_STR_FMT "'", size, MD_STR_ARG(filename));
    } else if (md_file_read(file, buf, size) != size) {
        MD_LOG_ERROR("MOLDEN: failed to read '" MD_STR_FMT "'", MD_STR_ARG(filename));
    } else {
        buf[size] = '\0';
        result = md_molden_system_init_from_str(sys, state, (md_str_t){ buf, size });
    }

    if (buf) {
        md_free(alloc, buf, size + 1);
    }
    md_file_close(&file);
    return result;
}

bool md_molden_file_is_molden(md_str_t filename) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        return false;
    }

    char buf[512];
    const size_t n = md_file_read(file, buf, sizeof(buf));
    md_file_close(&file);

    md_str_t str = { buf, n };
    md_str_t line;
    while (md_parse_next_line(&line, &str)) {
        line = md_str_trim(line);
        if (md_str_empty(line)) {
            continue;
        }
        md_str_t name, arg;
        return molden_section_header(&name, &arg, line) && md_str_eq_cstr_ignore_case(name, "molden format");
    }
    return false;
}
