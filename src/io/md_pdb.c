#include <md_pdb.h>
#include <md_io.h>
#include <math.h>

#include <md_system.h>
#include <md_run.h>
#include <md_infer.h>
#include <core/md_sort.h>
#include <md_element.h>

#include <core/md_common.h>
#include <core/md_array.h>
#include <core/md_str.h>
#include <core/md_parse.h>
#include <core/md_arena_allocator.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_file.h>
#include <core/md_path.h>
#include <core/md_vec_math.h>
#include <core/md_parse_simd.h>
#include <core/md_hash.h>
#include <core/md_executor.h>
#include <core/md_reader.h>

#include <hybrid_36_c.h>

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

// The PDB loader reads what a system is made of, straight from the text into the system: the atoms of the first
// model (ATOM, HETATM), where TER records end the chains, the cell of the first CRYST1 record, the bonds of the CONECT
// records and the biomolecules of REMARK 350. Every other record is passed over. A file of several models is a run as
// well (md_pdb_system_publish_run), whose index is found by a scan of the record names alone.

#define MD_PDB_CACHE_MAGIC 0x89172bab124545
#define MD_PDB_CACHE_VERSION 18     // 18: a model's atoms are those the system has (its first alternate locations)

// ### COLUMNS ###
// Columns as the format numbers them, from 1, both ends included

// The columns [beg, end], 1 based, of which a line cut short (trailing blanks trimmed) may hold only the first
static inline float extract_float(md_str_t line, size_t beg, size_t end) {
    if (line.len < beg) return 0.0f;
    return (float)md_parse_float(md_str_substr(line, beg - 1, MD_MIN(end, line.len) - beg + 1));
}

// A column the format writes with a fixed number of decimals ("%8.3f" for coordinates, "%6.2f" for
// occupancy and temperature factor): read without a scan when it is laid out that way, and parsed like
// any other column when it is not
static MD_FORCE_INLINE float extract_fixed_float(md_str_t line, size_t beg, size_t end, size_t decimals) {
    if (line.len < end) return 0.0f;
    const md_str_t field = md_str_substr(line, beg - 1, end-beg + 1);
    float value;
    return md_parse_fixed_f32(&value, field, decimals) ? value : (float)md_parse_float(field);
}

static inline char extract_char(md_str_t line, size_t idx) {
    if (line.len < idx) return 0;
    return line.ptr[idx - 1];
}

// An atom serial or a residue sequence number: a decimal number, past the width of its field hybrid-36 (A0000, then
// a0000: a letter first). 0 where the field holds neither, counted in num_bad: "*****" past the width, which many
// writers put there, or the hexadecimal VMD writes past 99999 ("186a0" for 100000), which is not hybrid-36 (a number
// which begins with a digit is decimal) and cannot be told from decimal where it has no letter (18700, 100096). Such a
// file is read by the order of its atoms; only CONECT records refer to atoms by their serials.
static inline int32_t parse_id(md_str_t line, size_t beg, size_t end, uint32_t* num_bad) {
    const uint32_t len = (uint32_t)(end - beg) + 1;
    const char* c = line.ptr + beg - 1;
    // Most are blanks and then digits (or blanks alone, 0), which hybrid-36 reads as decimal too
    uint32_t i = 0;
    while (i < len && c[i] == ' ') ++i;
    int32_t value = 0;
    uint32_t j = i;
    while (j < len && md_ascii_is_digit(c[j])) value = value * 10 + (c[j++] - '0');
    if (j == len) return value;
    // hy36decode reads hybrid-36 where the field begins with a letter and decimal otherwise, where only a '-' or a blank
    // (after digits, which it takes as a 0) is something more: anything else ('*****', VMD's hexadecimal) it rejects
    int32_t res = 0;
    if ((md_ascii_is_alpha(c[0]) || c[j] == '-' || c[j] == ' ') && !hy36decode(len, c, len, &res)) {
        return res;
    }
    if (num_bad) *num_bad += 1;
    return 0;
}

// The formal charge of the charge column (79 to 80): a digit and its sign, "1-" or "2+" (and the sign first, as some
// writers have it). NAN where there is none, which chemistry perception takes as not given (md_chem_perceive).
static MD_FORCE_INLINE float pdb_charge(md_str_t field) {
    field = md_str_trim(field);
    if (field.len == 1 && field.ptr[0] == '0') return 0.0f;
    if (field.len != 2) return NAN;
    const char a = field.ptr[0], b = field.ptr[1];
    if (md_ascii_is_digit(a) && (b == '+' || b == '-')) return (float)((b == '-' ? -1 : 1) * (a - '0'));
    if (md_ascii_is_digit(b) && (a == '+' || a == '-')) return (float)((a == '-' ? -1 : 1) * (b - '0'));
    return NAN;
}

// ### RECORDS ###

typedef enum {
    PDB_LINE_NONE,          // a record nothing reads (HEADER, REMARK, HELIX, SHEET, ...)
    PDB_LINE_ATOM,          // ATOM or HETATM, kept
    PDB_LINE_ATOM_SKIPPED,  // an alternate location other than the first
    // The records which are not of atoms, read in order once the atoms are placed
    PDB_LINE_TER,
    PDB_LINE_MODEL,
    PDB_LINE_ENDMDL,
    PDB_LINE_CRYST1,
    PDB_LINE_CONECT,
    PDB_LINE_REMARK_350,
} pdb_line_kind_t;

// The first eight characters of a line as a word, as they are in memory (little endian), and of those the first n
#define PDB_W8(a, b, c, d, e, f, g, h) ( \
      (uint64_t)(uint8_t)(a)        | ((uint64_t)(uint8_t)(b) << 8)  | ((uint64_t)(uint8_t)(c) << 16) | ((uint64_t)(uint8_t)(d) << 24) \
    | ((uint64_t)(uint8_t)(e) << 32) | ((uint64_t)(uint8_t)(f) << 40) | ((uint64_t)(uint8_t)(g) << 48) | ((uint64_t)(uint8_t)(h) << 56))
#define PDB_MASK(n) (~UINT64_C(0) >> (64 - 8 * (n)))

// The first eight characters of a line, blanks past the end of a shorter one: a record of its name alone (TER, ENDMDL,
// CRYST) is read as the same line with blanks after it. No load past the line, which may end the string.
static MD_FORCE_INLINE uint64_t pdb_line_head(md_str_t line) {
    uint64_t w;
    if (line.len >= 8) {
        md_memcpy(&w, line.ptr, 8);
    } else {
        w = 0x2020202020202020ull;
        md_memcpy(&w, line.ptr, line.len);
    }
    return w;
}

// What a line is, by its record name (columns 1 to 6): of its first eight characters, those of the name each record is
// told by, the atoms first (they are the most). TER and ATOM by their first three and four (TER may be followed by
// anything, and ATOM by a serial which runs into columns 5 and 6), CRYST1 also as CRYST with a blank where the 1 is, as
// some tools write it (the Conformational Analysis Tools of CHARMM-GUI's trajectories, with all angles 0).
static MD_FORCE_INLINE pdb_line_kind_t pdb_line_kind(md_str_t line) {
    const uint64_t w  = pdb_line_head(line);
    const uint64_t w6 = w & PDB_MASK(6);
    if (((w & PDB_MASK(4)) == PDB_W8('A','T','O','M',0,0,0,0)) | (w6 == PDB_W8('H','E','T','A','T','M',0,0))) {
        // The first alternate location only: blank or A
        const char alt_loc = extract_char(line, 17);
        return (alt_loc == ' ' || alt_loc == 'A') ? PDB_LINE_ATOM : PDB_LINE_ATOM_SKIPPED;
    }
    if ((w & PDB_MASK(3)) == PDB_W8('T','E','R',0,0,0,0,0))                         return PDB_LINE_TER;
    if ((w & PDB_MASK(5)) == PDB_W8('M','O','D','E','L',0,0,0))                     return PDB_LINE_MODEL;
    if (w6 == PDB_W8('E','N','D','M','D','L',0,0))                                  return PDB_LINE_ENDMDL;
    if (w6 == PDB_W8('C','O','N','E','C','T',0,0))                                  return PDB_LINE_CONECT;
    if (w6 == PDB_W8('C','R','Y','S','T','1',0,0) || w6 == PDB_W8('C','R','Y','S','T',' ',0,0)) return PDB_LINE_CRYST1;
    if (w == PDB_W8('R','E','M','A','R','K',' ','3') && line.len >= 10 && line.ptr[8] == '5' && line.ptr[9] == '0') return PDB_LINE_REMARK_350;
    return PDB_LINE_NONE;
}

// The cell of a CRYST1 record. No angles at all is a rectangular box whose angles were left out, not a degenerate cell.
// False for a cell of 1 x 1 x 1 at right angles, which writers put there for a structure of no cell (no periodic
// boundaries).
static bool pdb_cryst1_cell(md_unitcell_t* cell, md_str_t line) {
    float a = extract_float(line, 7, 15), b = extract_float(line, 16, 24), c = extract_float(line, 25, 33);
    float alpha = extract_float(line, 34, 40), beta = extract_float(line, 41, 47), gamma = extract_float(line, 48, 54);
    if (alpha == 0 && beta == 0 && gamma == 0) {
        alpha = beta = gamma = 90;
    }
    if (a == 1 && b == 1 && c == 1 && alpha == 90 && beta == 90 && gamma == 90) {
        return false;
    }
    // A cell of no volume (an angle missing or out of range, a length of zero) is none rather than one of infinities
    const double ca = cos(MD_DEG_TO_RAD(alpha)), cb = cos(MD_DEG_TO_RAD(beta)), cg = cos(MD_DEG_TO_RAD(gamma));
    if (!(a > 0 && b > 0 && c > 0 && alpha > 0 && alpha < 180 && beta > 0 && beta < 180 && gamma > 0 && gamma < 180) ||
        !(1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg > 1.0e-6)) {
        MD_LOG_INFO("PDB: the CRYST1 record describes no cell (%g %g %g %g %g %g)", a, b, c, alpha, beta, gamma);
        return false;
    }
    *cell = md_unitcell_from_extent_and_angles(a, b, c, alpha, beta, gamma);
    return true;
}

// A group of chains of a biomolecule of REMARK 350: the biomolecule copies the chains of the group ("APPLY THE FOLLOWING
// TO CHAINS:" and the "AND CHAINS:" lines after it) by each of the BIOMT transforms which follow the group
typedef struct pdb_group_t {
    int32_t  id;                // The number of the biomolecule, which its groups share
    uint32_t chain_offset;      // Its chains, a character each, in chains
    uint32_t chain_count;
    uint32_t transform_offset;  // Its transforms, in transforms
    uint32_t transform_count;
} pdb_group_t;

// What the records which are not of atoms give the system
typedef struct pdb_extras_t {
    bool    seen_cryst1;
    bool    has_cell;           // The first CRYST1 record gives one
    md_unitcell_t cell;
    md_array(int32_t) conect;   // The serials of the bonds of the CONECT records, two for each
    int32_t biomolecule;        // REMARK 350: the one whose groups of chains follow
    md_array(pdb_group_t) groups;
    md_array(char)   chains;
    md_array(md_float4x4_t) transforms;
    md_allocator_i*  alloc;
} pdb_extras_t;

// A serial of a CONECT record, 0 where there is none. Serials are hybrid-36 as those of the atoms are, past 99999 atoms.
static inline int32_t conect_serial(md_str_t line, size_t beg) {
    const size_t end = beg + 4;
    if (line.len < end) return 0;
    if (md_str_empty(md_str_trim(md_str_substr(line, beg - 1, 5)))) return 0;
    return parse_id(line, beg, end, NULL);
}

static void pdb_extras_record(pdb_extras_t* x, pdb_line_kind_t kind, md_str_t line) {
    md_allocator_i* alloc = x->alloc;
    switch (kind) {
    case PDB_LINE_CRYST1:
        // The first one: a file of models has one for each, the same
        if (!x->seen_cryst1) {
            x->seen_cryst1 = true;
            x->has_cell = pdb_cryst1_cell(&x->cell, md_str_trim(line));
        }
        break;
    case PDB_LINE_CONECT: {
        // The atom of columns 7 to 11 and up to four atoms bonded to it, each a bond. A record of more is continued on
        // the next line, which names the atom again.
        line = md_str_trim(line);
        const int32_t a = conect_serial(line, 7);
        if (!a) break;
        for (size_t beg = 12; beg <= 27; beg += 5) {
            const int32_t b = conect_serial(line, beg);
            if (!b) continue;
            md_array_push(x->conect, a, alloc);
            md_array_push(x->conect, b, alloc);
        }
        break;
    }
    case PDB_LINE_REMARK_350: {
        // A biomolecule is one or more groups of chains, each followed by the BIOMT transforms which copy them
        const md_str_t rest = md_str_substr(line, 11, SIZE_MAX);
        if (md_str_begins_with(rest, MD_STR_LIT("BIOMOLECULE:"))) {
            x->biomolecule = (int32_t)md_parse_int(md_str_trim(md_str_substr(rest, 12, SIZE_MAX)));
        }
        else if (md_str_begins_with(rest, MD_STR_LIT("APPLY THE FOLLOWING TO CHAINS:")) || md_str_begins_with(md_str_trim(rest), MD_STR_LIT("AND CHAINS:"))) {
            // A group begins with the chains to apply to, which AND CHAINS lines continue
            if (md_str_begins_with(rest, MD_STR_LIT("APPLY"))) {
                const pdb_group_t group = {
                    .id = x->biomolecule,
                    .chain_offset = (uint32_t)md_array_size(x->chains),
                    .transform_offset = (uint32_t)md_array_size(x->transforms),
                };
                md_array_push(x->groups, group, alloc);
            }
            pdb_group_t* group = md_array_last(x->groups);
            if (group) {
                md_str_t chains;
                md_str_cut(NULL, &chains, rest, md_str_find_char(rest, ':'));
                md_str_t chain;
                while (md_parse_next_token_delim(&chain, &chains, ',')) {
                    chain = md_str_trim(chain);
                    if (md_str_empty(chain)) continue;
                    md_array_push(x->chains, chain.ptr[0], alloc);
                    group->chain_count += 1;
                }
            } else {
                MD_LOG_DEBUG("REMARK 350: chains of no biomolecule");
            }
        }
        else if (md_str_begins_with(md_str_trim(rest), MD_STR_LIT("BIOMT")) && line.len > 18) {
            // REMARK 350   BIOMTr   n  m[r][1] m[r][2] m[r][3]  t[r], the rows of transform n of the group
            pdb_group_t* group = md_array_last(x->groups);
            const int32_t row = md_ascii_char_to_digit(line.ptr[18]) - 1;
            md_str_t fields = md_str_substr(line, 19, SIZE_MAX);
            md_str_t tok[5];
            const size_t num_tok = md_parse_next_tokens(tok, MD_COUNTOF(tok), &fields);
            if (!group || row < 0 || row > 2 || num_tok < 5) {
                MD_LOG_DEBUG("REMARK 350: a malformed BIOMT record");
            } else {
                if (row == 0) {
                    md_array_push(x->transforms, md_float4x4_ident(), alloc);
                    group->transform_count += 1;
                }
                md_float4x4_t* transform = md_array_last(x->transforms);
                for (int c = 0; c < 3; ++c) {
                    transform->elem[c][row] = (float)md_parse_float(tok[1 + c]);
                }
                transform->elem[3][row] = (float)md_parse_float(tok[4]);
            }
        }
        break;
    }
    default:
        break;
    }
}

// The assemblies of REMARK 350: the chains of a group of a biomolecule are the instances of their ids (a chain, and
// the molecules it holds: its ligands, its waters). A transform in two groups (the identity of each) is one operator,
// so that the copy it makes is one, named by its number in the first group which has it.
static void pdb_assemblies(md_system_t* sys, const pdb_extras_t* x, md_allocator_i* temp_alloc) {
    const size_t num_groups = md_array_size(x->groups);
    const size_t num_chains = md_array_size(x->chains);
    const size_t num_transforms = md_array_size(x->transforms);
    if (num_groups == 0 || num_transforms == 0) return;

    md_array(md_float4x4_t)     transform = 0;
    md_array(md_label_t) oper_id = 0;
    md_array(md_assembly_gen_t) gen = 0;

    for (size_t a = 0; a < num_groups; ++a) {
        const pdb_group_t* group = x->groups + a;
        if (group->chain_offset + group->chain_count > num_chains || group->transform_offset + group->transform_count > num_transforms) {
            MD_LOG_ERROR("REMARK 350: a group of chains out of range");
            continue;
        }

        md_array(md_instance_idx_t) inst = 0;
        for (uint32_t c = 0; c < group->chain_count; ++c) {
            const char chain = x->chains[group->chain_offset + c];
            bool found = false;
            for (size_t i = 0; i < sys->instance.count; ++i) {
                const md_str_t auth = md_system_instance_auth_id(sys, i);
                if (auth.len == 1 && auth.ptr[0] == chain) {
                    md_array_push(inst, (md_instance_idx_t)i, temp_alloc);
                    found = true;
                }
            }
            if (!found) {
                MD_LOG_DEBUG("REMARK 350: biomolecule %d names chain '%c', which has no atoms", group->id, chain);
            }
        }

        md_array(uint32_t) oper = 0;
        for (uint32_t t = 0; t < group->transform_count; ++t) {
            const md_float4x4_t M = x->transforms[group->transform_offset + t];
            size_t k = 0;
            for (; k < md_array_size(transform); ++k) {
                if (md_memcmp(&transform[k], &M, sizeof(md_float4x4_t)) == 0) break;
            }
            if (k == md_array_size(transform)) {
                char buf[16];
                const int len = snprintf(buf, sizeof(buf), "%u", t + 1);
                md_array_push(transform, M, temp_alloc);
                md_array_push(oper_id, md_label_make((md_str_t){buf, (size_t)MD_MAX(len, 0)}), temp_alloc);
            }
            md_array_push(oper, (uint32_t)k, temp_alloc);
        }

        char buf[16];
        const int len = snprintf(buf, sizeof(buf), "%d", group->id);
        const md_assembly_gen_t row = {
            .assembly_id = md_label_make((md_str_t){buf, (size_t)MD_MAX(len, 0)}),
            .inst        = inst,
            .inst_count  = md_array_size(inst),
            .oper        = oper,
            .oper_count  = md_array_size(oper),
        };
        md_array_push(gen, row, temp_alloc);
    }

    md_system_assemblies_init(sys, gen, md_array_size(gen), transform, oper_id, md_array_size(transform));
}

// ---------------------------------------------------------------------------------------------------------------------
// THE SYSTEM
//
// What a load gives is the core: the atoms, where they are, their types, the residues, the units of the file (its
// chains), the bonds the file gives (CONECT), the cell and the assemblies. Nothing is derived from the coordinates here,
// neither the bonds nor what the residues are by them: that is md_system_infer's, which a caller runs when (and
// where, on what thread) it wants.

// The atoms of serials, by the serials of the atoms (0 where a serial is not a number of the format, which no record
// refers to). Where the serials count up by one from the first, as most writers number them, an atom is its serial less
// the first; otherwise the serials are sorted, the atoms of a serial in their order, and a serial is the first of them.
typedef struct pdb_serial_index_t {
    uint32_t        num_atoms;
    uint32_t        first;      // The serial of atom 0, where the serials count up by one: else 0
    const int32_t*  serial;
    uint32_t*       order;      // The atoms by their serials, where they do not
} pdb_serial_index_t;

static void pdb_serial_index_init(pdb_serial_index_t* idx, const int32_t* serial, size_t num_atoms, md_allocator_i* alloc) {
    *idx = (pdb_serial_index_t){ .num_atoms = (uint32_t)num_atoms, .serial = serial };
    if (num_atoms == 0) return;
    bool counts_up = serial[0] > 0;
    for (size_t i = 1; i < num_atoms && counts_up; ++i) counts_up = serial[i] == serial[0] + (int32_t)i;
    if (counts_up) {
        idx->first = (uint32_t)serial[0];
        return;
    }
    idx->order = md_alloc(alloc, num_atoms * sizeof(uint32_t));
    md_argsort_radix_u32(idx->order, (const uint32_t*)serial, num_atoms);
}

// The atom of a serial, or UINT32_MAX for none
static inline uint32_t pdb_serial_index_find(const pdb_serial_index_t* idx, int32_t s) {
    if (s == 0) return UINT32_MAX;
    if (idx->first) {
        const uint32_t i = (uint32_t)s - idx->first;
        return i < idx->num_atoms ? i : UINT32_MAX;
    }
    if (!idx->order) return UINT32_MAX;
    // The first of the atoms whose serial is not less than s
    const uint32_t key = (uint32_t)s;
    size_t lo = 0, hi = idx->num_atoms;
    while (lo < hi) {
        const size_t mid = (lo + hi) / 2;
        if ((uint32_t)idx->serial[idx->order[mid]] < key) lo = mid + 1;
        else hi = mid;
    }
    return (lo < idx->num_atoms && (uint32_t)idx->serial[idx->order[lo]] == key) ? idx->order[lo] : UINT32_MAX;
}

// The bonds of the CONECT records (pairs of serials), of origin MD_BOND_ORIGIN_FILE, in the order the records first
// give them. Each bond is in the records of both its atoms, and some writers repeat it for its bond order: the pairs
// are sorted to find the first of each, which is linear where a map of pairs of neighbouring atoms is not.
static void pdb_conect_bonds(md_system_t* sys, const int32_t* atom_serial, const int32_t* conect, size_t num_conect, md_allocator_i* temp_arena) {
    pdb_serial_index_t idx;
    pdb_serial_index_init(&idx, atom_serial, sys->atom.count, temp_arena);

    uint64_t* key = md_alloc(temp_arena, num_conect * sizeof(uint64_t));
    size_t num_keys = 0;
    for (size_t i = 0; i < num_conect; ++i) {
        const uint32_t a = pdb_serial_index_find(&idx, conect[2 * i + 0]);
        const uint32_t b = pdb_serial_index_find(&idx, conect[2 * i + 1]);
        if (a == UINT32_MAX || b == UINT32_MAX || a == b) continue;
        key[num_keys++] = ((uint64_t)MD_MIN(a, b) << 32) | MD_MAX(a, b);
    }
    if (num_keys == 0) return;

    // Stable: of the keys of one pair, the first record's comes first
    uint32_t* order = md_alloc(temp_arena, num_keys * sizeof(uint32_t));
    md_argsort_radix_u64(order, key, num_keys);
    uint8_t* keep = md_alloc(temp_arena, num_keys);
    md_memset(keep, 0, num_keys);
    size_t num_bonds = 0;
    for (size_t i = 0; i < num_keys; ++i) {
        if (i == 0 || key[order[i]] != key[order[i - 1]]) {
            keep[order[i]] = 1;
            num_bonds += 1;
        }
    }

    const md_bond_flags_t file = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_FILE);
    md_array_ensure(sys->bond.pairs, sys->bond.count + num_bonds, sys->alloc);
    md_array_ensure(sys->bond.flags, sys->bond.count + num_bonds, sys->alloc);
    for (size_t i = 0; i < num_keys; ++i) {
        if (!keep[i]) continue;
        const md_atom_pair_t pair = {{ (md_atom_idx_t)(key[i] >> 32), (md_atom_idx_t)(key[i] & 0xFFFFFFFF) }};
        md_array_push_no_grow(sys->bond.pairs, pair);
        md_array_push_no_grow(sys->bond.flags, file);
    }
    sys->bond.count += num_bonds;
    md_bond_build_connectivity(&sys->bond, sys->atom.count, sys->alloc);
}

// What follows the atoms and residues: the cell, the bonds of the CONECT records (by the serials of the atoms loaded),
// the units of the file (the runs of residues of one chain id) and the assemblies of REMARK 350
static void pdb_system_finish(md_system_t* sys, md_system_state_t* state, const pdb_extras_t* x, const int32_t* atom_serial, const md_str_t* res_chain_ids, const bool* res_after_ter, md_allocator_i* temp_arena) {
    if (x->has_cell) {
        state->unitcell = x->cell;
    }

    // The bonds of the CONECT records, of origin MD_BOND_ORIGIN_FILE: inference keeps them and adds the bonds they leave
    // out (CONECT records usually cover the ligands only, in a coarse grained structure from martinize2 all of it).
    // Their serials name atoms of the first model; alternate locations other than the first are not loaded.
    const size_t num_conect = md_array_size(x->conect) / 2;
    if (num_conect > 0 && atom_serial) {
        pdb_conect_bonds(sys, atom_serial, x->conect, num_conect, temp_arena);
    }

    // The units of the file are its chains: the runs of residues of one chain id, which a TER record ends (two chains of
    // one id written one after the other are two: the segments of an EM structure)
    md_system_instances_from_chain_ids(sys, res_chain_ids, res_after_ter);

    // The assemblies name the chains by the ids they have in the file, which are those of the instances
    pdb_assemblies(sys, x, temp_arena);
}

// ---------------------------------------------------------------------------------------------------------------------
// THE LOAD
//
// Straight from the text, a block at a time (read over the threads of the executor where the source reads at offsets,
// a file), each block in passes and each pass in pieces at line boundaries, a few for each thread (md_reader.h: blocks
// and pieces). Three loops over the pieces of a pass, in parallel:
//
//   count  what a piece holds: its atoms (an alternate location other than the first is none), the residues which
//          begin at its atoms after its first, its other records of interest, the residue fields of its first and last
//          atom and whether a TER record is before the first or after the last
//   parse  each of its atoms where it goes in the arrays of the load, its residues where they go in the residues, its
//          records where they go in the records of the pass, and each atom's type in a table of the piece's own
//   type   each atom the system's type of its table's
//
// Between them the calling thread works with the pieces, never with their atoms: it places them (where their atoms,
// residues and records go, and whether the first atom of each begins a residue, which the last atom before it decides),
// makes the types of their tables the system's, piece by piece, which is the order the atoms first have them, and reads
// the few records of interest in order.
//
// Of a file of several models the atoms of the first are loaded (all of them with MD_LOAD_FLAG_ALL_MODELS). Once
// it has ended the load goes to the end of the file, past the last MODEL or ENDMDL record, for the CONECT records.

#define PDB_RES_FIELDS_BEG  17      // The residue name, the chain id, the sequence number and the insertion code:
#define PDB_RES_FIELDS_END  27      // columns 18 to 27

// The residue fields of an atom line, its ten characters as two words (blanks past the end of a shorter line): the chain
// id, as it is in the file, the fifth
typedef struct pdb_res_fields_t {
    uint64_t lo;
    uint16_t hi;
} pdb_res_fields_t;

static MD_FORCE_INLINE pdb_res_fields_t pdb_res_fields(md_str_t line) {
    pdb_res_fields_t f;
    if (line.len >= PDB_RES_FIELDS_END) {
        md_memcpy(&f.lo, line.ptr + PDB_RES_FIELDS_BEG, 8);
        md_memcpy(&f.hi, line.ptr + PDB_RES_FIELDS_BEG + 8, 2);
    } else {
        char c[PDB_RES_FIELDS_END - PDB_RES_FIELDS_BEG];
        md_memset(c, ' ', sizeof(c));
        if (line.len > PDB_RES_FIELDS_BEG) md_memcpy(c, line.ptr + PDB_RES_FIELDS_BEG, line.len - PDB_RES_FIELDS_BEG);
        md_memcpy(&f.lo, c, 8);
        md_memcpy(&f.hi, c + 8, 2);
    }
    return f;
}

static MD_FORCE_INLINE bool pdb_res_fields_eq(const pdb_res_fields_t* a, const pdb_res_fields_t* b) {
    return a->lo == b->lo && a->hi == b->hi;
}

static MD_FORCE_INLINE char pdb_res_fields_chain(const pdb_res_fields_t* f) {
    return (char)(f->lo >> (8 * (21 - PDB_RES_FIELDS_BEG)));
}

// A field of up to 4 columns from offset (from 0) as a label, without the blanks around it: what is at or below ' ' (a
// tab, a carriage return, a nul) at its ends, which no name has. One word, no loop over its characters.
static MD_FORCE_INLINE md_label_t pdb_label(md_str_t line, size_t offset, size_t width) {
    uint32_t w = 0x20202020u;
    if (line.len >= offset + width) {
        if (width == 4) md_memcpy(&w, line.ptr + offset, 4);
        else            { uint16_t h; md_memcpy(&h, line.ptr + offset, 2); w = 0x20200000u | h; }
    } else if (line.len > offset) {
        md_memcpy(&w, line.ptr + offset, line.len - offset);
    }
    // 0x80 in every byte above ' ': its low 7 bits are, or its high bit is set (no borrow crosses a byte)
    const uint32_t mark = (((w | 0x80808080u) - 0x21212121u) | w) & 0x80808080u;
    md_label_t lbl = {0};
    if (mark) {
        const uint32_t lo = md_ctz32(mark) >> 3;
        const uint32_t hi = (31 - md_clz32(mark)) >> 3;
        const uint64_t v = ((uint64_t)w >> (8 * lo)) & ((UINT64_C(1) << (8 * (hi - lo + 1))) - 1);
        md_memcpy(lbl.buf, &v, sizeof(v));
    }
    return lbl;
}

// The three coordinates of an atom line, "%8.3f" each from column 31: at once (md_parse_fixed_float3_w8), and one by one
// (md_parse_fixed_f32, else md_parse_float) where any is laid out otherwise or the line is too short for them
static inline void pdb_coords_scalar(float out[3], md_str_t line) {
    out[0] = extract_fixed_float(line, 31, 38, 3);
    out[1] = extract_fixed_float(line, 39, 46, 3);
    out[2] = extract_fixed_float(line, 47, 54, 3);
}

static MD_FORCE_INLINE void pdb_coords(float out[3], md_str_t line) {
    if (line.len < 54 || !md_parse_fixed_float3_w8(out, line.ptr + 30, 3)) pdb_coords_scalar(out, line);
}

typedef struct pdb_residue_start_t {
    uint32_t         atom;      // of the load
    md_sequence_id_t seq_id;    // INT_MAX where the field is not a number of the format
    md_label_t       name;
    char             chain_id;  // As the file has it: ' ' for none
    bool             after_ter; // A TER record before it: it begins a unit of the file
} pdb_residue_start_t;

typedef struct pdb_record_t {
    size_t          offset;     // of the line in the pass
    size_t          len;
    size_t          atom;       // the atoms of the load before it
    pdb_line_kind_t kind;
} pdb_record_t;

// The types of the atoms of a piece, by their names and elements, numbered from 1 in the order the atoms first have
// them: what the parse gives an atom until the system's types replace them. A piece of more types than the table holds
// (no file of atoms has as many atom names in 256 KB) gives the atoms past them PDB_TYPE_OVERFLOW, typed one by one.
#define PDB_PIECE_TYPES     1024
#define PDB_TYPE_SLOTS      2048
#define PDB_TYPE_OVERFLOW   ((md_particle_type_idx_t)-1)

typedef struct pdb_types_t {
    uint32_t    count;
    bool        overflow;
    uint16_t    slot[PDB_TYPE_SLOTS];       // entry + 1, 0 for none
    uint64_t    key[PDB_PIECE_TYPES];       // md_label_key of the name
    md_label_t  name[PDB_PIECE_TYPES];
    uint8_t     z[PDB_PIECE_TYPES];
    md_particle_type_idx_t type[PDB_PIECE_TYPES];  // The system's, once made
} pdb_types_t;

static inline md_particle_type_idx_t pdb_types_find_or_add(pdb_types_t* t, md_label_t name, uint8_t z) {
    const uint64_t key = md_label_key(name);
    uint32_t h = (uint32_t)(((key ^ ((uint64_t)z << 56)) * 0x9E3779B97F4A7C15ull) >> 53) & (PDB_TYPE_SLOTS - 1);
    for (;;) {
        const uint32_t s = t->slot[h];
        if (s == 0) {
            if (t->count == PDB_PIECE_TYPES) {
                t->overflow = true;
                return PDB_TYPE_OVERFLOW;
            }
            const uint32_t e = t->count++;
            t->key[e] = key;
            t->name[e] = name;
            t->z[e] = z;
            t->slot[h] = (uint16_t)(e + 1);
            return (md_particle_type_idx_t)(e + 1);
        }
        if (t->key[s - 1] == key && t->z[s - 1] == z) return (md_particle_type_idx_t)s;
        h = (h + 1) & (PDB_TYPE_SLOTS - 1);
    }
}

enum { PDB_COL_OCCUPANCY, PDB_COL_B_FACTOR, PDB_COL_CHARGE, PDB_NUM_COLUMNS };

// A line of a pass the parse reads, as the count found it: an atom kept, or a record (TER to REMARK 350). The lines are
// split and told apart once, by the count: the parse takes them from these. Of an atom, whether it begins a residue and
// whether a TER record is between it and the atom before (in its piece: its first atom's are the place's to tell).
#define PDB_LINE_MAX_LEN ((1u << 26) - 1)       // a longer line (no record or atom is) is cut at that

typedef struct pdb_line_t {
    uint32_t offset;            // of the pass
    uint32_t len        : 26;
    uint32_t kind       : 4;    // pdb_line_kind_t
    uint32_t begins_res : 1;
    uint32_t after_ter  : 1;
} pdb_line_t;

// No line of these is shorter than 3 characters and its newline: of len bytes, a piece has at most len / 4 + 1
#define PDB_LINES_OF(len) ((len) / 4 + 1)

typedef struct pdb_piece_t {
    md_text_piece_t range;
    pdb_types_t* types;
    pdb_line_t* lines;          // its lines of the index, room for PDB_LINES_OF its bytes
    // The count
    uint32_t num_lines;
    uint32_t num_atoms;         // kept
    uint32_t num_records;
    uint32_t num_inner_res;     // residues which begin at an atom of the piece after its first
    bool     ter_first;         // a TER record before its first atom; of a piece of no atoms, anywhere
    bool     ter_last;          // a TER record after its last atom
    pdb_res_fields_t first, last;   // the residue fields of its first and last atom
    // Where it goes
    size_t   first_atom;        // of the load
    size_t   first_record;      // of the pass
    size_t   first_res;         // of the load
    bool     first_begins_res;  // its first atom begins a residue
    bool     first_after_ter;   // and a TER record is between it and the atom before
    // The parse
    uint32_t num_bad_ids;
    uint32_t num_bad_seq;
    bool     unassigned;        // an atom without an element
    uint32_t col_first[PDB_NUM_COLUMNS];    // the bits of each column's value of the first atom
    bool     col_varies[PDB_NUM_COLUMNS];   // another atom's has other bits
} pdb_piece_t;

static inline uint32_t pdb_f32_bits(float x) {
    uint32_t u;
    md_memcpy(&u, &x, sizeof(u));
    return u;
}

// The atoms of a piece as it left the columns: written where a column varies in the piece, else all of the one value
// of its first atom, which pdb_fill_columns writes once the column is known to vary over the load
typedef struct pdb_col_run_t {
    size_t   first_atom;
    uint32_t num_atoms;
    uint32_t col_first[PDB_NUM_COLUMNS];
    bool     col_varies[PDB_NUM_COLUMNS];
} pdb_col_run_t;

// Where the time of a load goes, in the debug log: of the passes (the reads, the three loops and what the calling
// thread does between them) and of making the system
enum { PDB_PHASE_READ, PDB_PHASE_COUNT, PDB_PHASE_PLACE, PDB_PHASE_PARSE, PDB_PHASE_GATHER, PDB_PHASE_TYPE, PDB_PHASE_RECORDS,
       PDB_PHASE_ATOMS, PDB_PHASE_RESIDUES, PDB_PHASE_TYPING, PDB_PHASE_FINISH, PDB_NUM_PHASES };

// What a load has: the pass it reads, and what it has read
typedef struct pdb_load_t {
    md_system_t*        sys;
    md_system_state_t*  state;
    md_allocator_i*     temp;
    bool                concat_models;

    // The pass
    const char*         text;
    md_array(md_text_piece_t) ranges;
    md_array(pdb_piece_t) pieces;
    md_array(pdb_record_t) records;
    pdb_types_t*        types;          // A table for each piece
    pdb_line_t*         lines;          // The index of the lines of a pass, the pieces' after one another
    size_t              lines_cap;
    bool                atoms;          // false once the first model ended: the atom lines are passed over

    // Of the load, by atom (the coordinates, types and flags where the state and the system keep them)
    md_array(int32_t)    serial;        // what CONECT records refer to atoms by
    md_array(md_label_t) atom_name;     // of the atoms without an element (pdb_fill_names gives the others theirs)
    md_array(float)      column[PDB_NUM_COLUMNS];  // occupancy, b factor, formal charge: of the pieces where they vary
    md_array(pdb_col_run_t) col_runs;   // what each piece left of the columns (pdb_fill_columns writes the rest)
    md_array(pdb_residue_start_t) res_starts;
    size_t              num_atoms;
    size_t              num_res;

    // The last atom before the pass: whether there is one, its residue fields and whether a TER record is after it
    bool                has_prev;
    bool                prev_ter;
    pdb_res_fields_t    prev;

    size_t              num_models;
    size_t              model_beg;      // The first model: [model_beg, num_atoms) once its MODEL record is seen
    uint32_t            col_first[PDB_NUM_COLUMNS];
    bool                col_varies[PDB_NUM_COLUMNS];    // A column of one value is not published
    uint32_t            num_bad_ids;
    uint32_t            num_bad_seq;
    bool                unassigned;
    md_hashmap32_t      type_map;
    pdb_extras_t        extras;
    double              ms[PDB_NUM_PHASES];     // where the time went, for the debug log
} pdb_load_t;

static void pdb_count_lines(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    const pdb_load_t* L = (const pdb_load_t*)arg;
    for (uint32_t p = beg; p < end; ++p) {
        pdb_piece_t* piece = &L->pieces[p];
        // Few counts in the loop, which has more to keep than registers: the lines from where the index ends, the records
        // from the lines and the atoms
        const char* text = L->text;
        pdb_line_t* out = piece->lines;
        uint32_t num_atoms = 0, num_inner_res = 0;
        bool ter = false, ter_first = false;
        pdb_res_fields_t prev = {0};
        md_parse_lines_t lines = md_parse_lines(L->text + piece->range.beg, L->text + piece->range.end);
        md_str_t line;
        while (md_parse_lines_next(&lines, &line)) {
            const pdb_line_kind_t kind = pdb_line_kind(line);
            if (kind == PDB_LINE_ATOM) {
                if (!L->atoms) continue;
                // A residue is its name, its sequence number and its insertion code in its chain: 100, 100A and 100B
                // are three residues of the same number (the insertions of Kabat numbered antibodies), also when they
                // share a name. A TER record ends the residue before it, whatever the residue after it is called, and
                // the chain: what follows is a unit of its own (md_system_instances_from_chain_ids).
                const pdb_res_fields_t fields = pdb_res_fields(line);
                bool begins = false;
                if (num_atoms == 0) {
                    piece->first = fields;
                    ter_first = ter;
                } else if (ter || !pdb_res_fields_eq(&fields, &prev)) {
                    num_inner_res += 1;
                    begins = true;
                }
                *out++ = (pdb_line_t){ .offset = (uint32_t)(line.ptr - text), .len = (uint32_t)MD_MIN(line.len, PDB_LINE_MAX_LEN),
                                       .kind = PDB_LINE_ATOM, .begins_res = begins, .after_ter = ter };
                prev = fields;
                ter = false;
                num_atoms += 1;
            } else if (kind >= PDB_LINE_TER) {
                if (kind == PDB_LINE_TER) ter = true;
                *out++ = (pdb_line_t){ .offset = (uint32_t)(line.ptr - text), .len = (uint32_t)MD_MIN(line.len, PDB_LINE_MAX_LEN), .kind = kind };
            }
        }
        const uint32_t num_lines = (uint32_t)(out - piece->lines);
        const uint32_t num_records = num_lines - num_atoms;
        MD_ASSERT(num_lines <= PDB_LINES_OF(piece->range.end - piece->range.beg));
        piece->num_lines = num_lines;
        piece->num_atoms = num_atoms;
        piece->num_records = num_records;
        piece->num_inner_res = num_inner_res;
        piece->ter_first = num_atoms ? ter_first : ter;
        piece->ter_last = num_atoms ? ter : false;
        piece->last = prev;
    }
}

// A column of the atom lines as the last atom had it: its characters (up to 8, in the low bytes of a word), and the value
// they are. UINT64_MAX is none: no field of up to 6 characters is that.
typedef struct pdb_field_cache_t {
    uint64_t raw;
    float    value;
} pdb_field_cache_t;

// Columns beg to end (from 1, both included) of the line: the occupancy and the b factor ("%6.2f"), or the charge
static MD_FORCE_INLINE float pdb_field_cached(pdb_field_cache_t* cache, md_str_t line, size_t beg, size_t end) {
    const size_t width = end - beg + 1;
    uint64_t raw = 0;
    if (line.len >= end) {
        md_memcpy(&raw, line.ptr + beg - 1, width);
        raw |= (uint64_t)width << 56;                               // and how much of it there is
    } else {
        if (line.len >= beg) md_memcpy(&raw, line.ptr + beg - 1, line.len - beg + 1);     // The end of a shorter line
        raw |= (uint64_t)(line.len - MD_MIN(line.len, beg - 1)) << 56;
    }
    if (raw == cache->raw) return cache->value;
    cache->raw = raw;
    cache->value = width == 2 ? pdb_charge(md_str_substr(line, beg - 1, 2)) : extract_fixed_float(line, beg, end, 2);
    return cache->value;
}

static void pdb_parse_lines(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    const pdb_load_t* L = (const pdb_load_t*)arg;
    md_float3_t* xyz = L->state->xyz;
    md_particle_type_idx_t* type_idx = L->sys->atom.type_idx;
    md_atom_flags_t* flags = L->sys->atom.flags;
    for (uint32_t p = beg; p < end; ++p) {
        pdb_piece_t* piece = &L->pieces[p];
        pdb_types_t* types = piece->types;
        piece->num_bad_ids = 0;
        piece->num_bad_seq = 0;
        piece->unassigned = false;
        types->count = 0;
        types->overflow = false;
        if (L->atoms && piece->num_atoms) md_memset(types->slot, 0, sizeof(types->slot));

        // The type of the atom before, which the next often shares (a water's hydrogens, an ion after another)
        uint64_t prev_key = 0;
        uint8_t  prev_z = 0;
        md_particle_type_idx_t prev_type = 0;
        // The occupancy, b factor and charge of most atoms are those of the atom before: the same characters, which are
        // not parsed again
        pdb_field_cache_t cache[PDB_NUM_COLUMNS] = { {UINT64_MAX}, {UINT64_MAX}, {UINT64_MAX} };
        // Of the columns: the bits of the first atom's value, and whether another atom's differ (kept here, not in the
        // piece, which the stores of the atoms' values would make the compiler read again for every atom)
        uint32_t col_first[PDB_NUM_COLUMNS] = {0};
        bool col_varies[PDB_NUM_COLUMNS] = {false};
        bool any_varies = false;
        // The characters of the columns of the atom before, as two words: an atom of the same has its values, and while
        // no column varies, nothing is written. UINT64_MAX in the second is none (a line too short for the words).
        uint64_t cols_lo = 0, cols_hi = UINT64_MAX;
        float col[PDB_NUM_COLUMNS] = {0};

        pdb_residue_start_t* res = L->res_starts + piece->first_res;
        pdb_record_t* rec = L->records + piece->first_record;
        size_t atom = piece->first_atom;
        // The lines of the piece as the count found them: its atoms (none once the first model ended) and records
        const pdb_line_t* lines = piece->lines;
        for (uint32_t i = 0; i < piece->num_lines; ++i) {
            const pdb_line_t ln = lines[i];
            const md_str_t line = { L->text + ln.offset, ln.len };
            if (ln.kind != PDB_LINE_ATOM) {
                rec->offset = ln.offset;
                rec->len    = ln.len;
                rec->atom   = atom;
                rec->kind   = (pdb_line_kind_t)ln.kind;
                ++rec;
                continue;
            }

            // Where a residue begins, as the count found it, and the first atom's as the place found it
            const bool first = atom == piece->first_atom;
            if (first ? piece->first_begins_res : ln.begins_res) {
                uint32_t bad = 0;
                const int32_t seq_id = line.len >= 26 ? parse_id(line, 23, 26, &bad) : 0;
                piece->num_bad_seq += bad;
                res->atom      = (uint32_t)atom;
                res->seq_id    = bad ? INT_MAX : seq_id;
                res->name      = pdb_label(line, 17, 4);
                const pdb_res_fields_t fields = pdb_res_fields(line);
                res->chain_id  = pdb_res_fields_chain(&fields);
                res->after_ter = first ? piece->first_after_ter : ln.after_ter;
                ++res;
            }

            pdb_coords(xyz[atom].elem, line);
            flags[atom] = MD_ATOM_FLAG_NONE;
            const md_label_t name = pdb_label(line, 12, 4);
            const md_label_t element = pdb_label(line, 76, 2);
            if (element.buf[0] == 0) {
                // Typed by its name with the others of its residue, once all are read
                L->atom_name[atom] = name;
                type_idx[atom] = 0;
                piece->unassigned = true;
            } else {
                const uint8_t z = (uint8_t)md_atomic_number_from_symbol(md_label_str(&element), true);
                const uint64_t key = md_label_key(name);
                if (key != prev_key || z != prev_z || prev_type == 0) {
                    prev_type = pdb_types_find_or_add(types, name, z);
                    prev_key = key;
                    prev_z = z;
                }
                type_idx[atom] = prev_type;
            }
            // Columns 55-66 (the occupancy and the b factor), and of 79-80 (the charge) what the line has of them and how
            // much: at most 0x2FFFF, so the second word is never UINT64_MAX
            uint64_t lo = 0, hi = UINT64_MAX;
            if (line.len >= 66) {
                uint32_t mid, charge = 0;
                md_memcpy(&lo, line.ptr + 54, 8);
                md_memcpy(&mid, line.ptr + 62, 4);
                if (line.len >= 80) {
                    uint16_t h;
                    md_memcpy(&h, line.ptr + 78, 2);
                    charge = 0x20000u | h;
                } else if (line.len == 79) {
                    charge = 0x10000u | (uint8_t)line.ptr[78];
                }
                hi = (uint64_t)charge << 32 | mid;
            }
            const bool same = hi != UINT64_MAX && lo == cols_lo && hi == cols_hi;
            cols_lo = lo;
            cols_hi = hi;
            if (!same) {
                col[PDB_COL_OCCUPANCY] = pdb_field_cached(&cache[PDB_COL_OCCUPANCY], line, 55, 60);
                col[PDB_COL_B_FACTOR]  = pdb_field_cached(&cache[PDB_COL_B_FACTOR],  line, 61, 66);
                col[PDB_COL_CHARGE]    = pdb_field_cached(&cache[PDB_COL_CHARGE],    line, 79, 80);
            }
            // A column of one value (the occupancy of 1 and b factor of 0 of a simulation, no charges) is neither written
            // nor published: told by bits, as publish tells it. Once one varies in the piece, the atoms before have the
            // first's value.
            if (atom == piece->first_atom) {
                for (int k = 0; k < PDB_NUM_COLUMNS; ++k) col_first[k] = pdb_f32_bits(col[k]);
            }
            // The same values as the atom before, while all are the first's: nothing to tell or write
            for (int k = 0; k < PDB_NUM_COLUMNS && (!same || any_varies); ++k) {
                if (!col_varies[k]) {
                    if (pdb_f32_bits(col[k]) == col_first[k]) continue;
                    col_varies[k] = true;
                    any_varies = true;
                    float first_value;
                    md_memcpy(&first_value, &col_first[k], sizeof(first_value));
                    for (size_t i = piece->first_atom; i < atom; ++i) L->column[k][i] = first_value;
                }
                L->column[k][atom] = col[k];
            }
            L->serial[atom] = line.len >= 11 ? parse_id(line, 7, 11, &piece->num_bad_ids) : 0;
            ++atom;
        }
        for (int k = 0; k < PDB_NUM_COLUMNS; ++k) {
            piece->col_first[k] = col_first[k];
            piece->col_varies[k] = col_varies[k];
        }
        MD_ASSERT(atom == piece->first_atom + (L->atoms ? piece->num_atoms : 0));
    }
}

// The atoms of each piece their types in the system: of the piece's table, made the system's
static void pdb_type_lines(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    const pdb_load_t* L = (const pdb_load_t*)arg;
    for (uint32_t p = beg; p < end; ++p) {
        const pdb_piece_t* piece = &L->pieces[p];
        const pdb_types_t* types = piece->types;
        md_particle_type_idx_t* t = L->sys->atom.type_idx + piece->first_atom;
        for (uint32_t i = 0; types->count && i < piece->num_atoms; ++i) {
            if (t[i] > 0) t[i] = types->type[t[i] - 1];
        }
    }
}

// The atoms of a piece whose types are past its table: their names and elements read again, and typed one by one (no
// file of atoms has as many atom names in a piece; this is to be right, not fast)
static void pdb_type_overflow(pdb_load_t* L, const pdb_piece_t* piece) {
    md_system_t* sys = L->sys;
    size_t atom = piece->first_atom;
    for (uint32_t i = 0; i < piece->num_lines; ++i) {
        const pdb_line_t ln = piece->lines[i];
        if (ln.kind != PDB_LINE_ATOM) continue;
        const md_str_t line = { L->text + ln.offset, ln.len };
        if (sys->atom.type_idx[atom] == PDB_TYPE_OVERFLOW) {
            const md_str_t name = md_str_trim(md_str_substr(line, 12, 4));
            const md_atomic_number_t z = md_atomic_number_from_symbol(md_str_trim(md_str_substr(line, 76, 2)), true);
            sys->atom.type_idx[atom] = md_system_particle_type_find_or_add(sys, &L->type_map, name, z);
        }
        ++atom;
    }
}

// Room for count elements in an array of the load, which grows by doubling past what was reserved
#define PDB_ROOM(arr, count, alloc) do { \
        if ((count) > md_array_capacity(arr)) md_array_ensure((arr), MD_MAX((size_t)(count), 2 * md_array_capacity(arr)), (alloc)); \
        md_array_resize((arr), (count), (alloc)); \
    } while (0)

// Where the atoms, residues and records of the pieces of a pass go, and whether the first atom of each begins a residue:
// as the last atom before it, which the pieces before it give, has it. Room for them.
static void pdb_place_pieces(pdb_load_t* L) {
    size_t num_atoms = L->num_atoms, num_res = L->num_res, num_records = 0;
    for (size_t p = 0; p < md_array_size(L->pieces); ++p) {
        pdb_piece_t* piece = &L->pieces[p];
        piece->first_atom = num_atoms;
        piece->first_record = num_records;
        piece->first_res = num_res;
        num_atoms += piece->num_atoms;
        num_records += piece->num_records;
        if (piece->num_atoms) {
            const bool ter = L->prev_ter || piece->ter_first;
            piece->first_begins_res = !L->has_prev || ter || !pdb_res_fields_eq(&piece->first, &L->prev);
            piece->first_after_ter = L->has_prev && ter;
            num_res += piece->first_begins_res + piece->num_inner_res;
            L->has_prev = true;
            L->prev = piece->last;
            L->prev_ter = piece->ter_last;
        } else {
            L->prev_ter |= piece->ter_first;
        }
    }

    md_allocator_i* temp = L->temp;
    PDB_ROOM(L->state->xyz, num_atoms, L->state->alloc);
    PDB_ROOM(L->sys->atom.type_idx, num_atoms, L->sys->alloc);
    PDB_ROOM(L->sys->atom.flags, num_atoms, L->sys->alloc);
    PDB_ROOM(L->serial, num_atoms, temp);
    PDB_ROOM(L->atom_name, num_atoms, temp);
    for (int c = 0; c < PDB_NUM_COLUMNS; ++c) PDB_ROOM(L->column[c], num_atoms, temp);
    PDB_ROOM(L->res_starts, num_res, temp);
    PDB_ROOM(L->records, num_records, temp);
}

// What the parse of a pass gives that is not of its atoms: the types of the pieces' tables, the system's in the order
// of the pieces (the order the atoms first have them), and what the pieces say of the columns and of the numbers
static void pdb_gather_pieces(pdb_load_t* L) {
    const size_t num_pieces = md_array_size(L->pieces);
    for (size_t p = 0; p < num_pieces; ++p) {
        const pdb_piece_t* piece = &L->pieces[p];
        pdb_types_t* t = piece->types;
        for (uint32_t e = 0; e < t->count; ++e) {
            t->type[e] = md_system_particle_type_find_or_add(L->sys, &L->type_map, md_label_str(&t->name[e]), (md_atomic_number_t)t->z[e]);
        }
        if (piece->num_atoms) {
            const bool first = piece->first_atom == 0;
            pdb_col_run_t run = { .first_atom = piece->first_atom, .num_atoms = piece->num_atoms };
            for (int c = 0; c < PDB_NUM_COLUMNS; ++c) {
                if (first) L->col_first[c] = piece->col_first[c];
                L->col_varies[c] |= piece->col_varies[c] || piece->col_first[c] != L->col_first[c];
                run.col_first[c] = piece->col_first[c];
                run.col_varies[c] = piece->col_varies[c];
            }
            md_array_push(L->col_runs, run, L->temp);
        }
        L->num_bad_ids += piece->num_bad_ids;
        L->num_bad_seq += piece->num_bad_seq;
        L->unassigned  |= piece->unassigned;
    }
}

// The records of a pass, in order: the models end the atoms, the others are the extras
static void pdb_read_records(pdb_load_t* L) {
    const size_t num_records = md_array_size(L->records);
    for (size_t r = 0; r < num_records; ++r) {
        const pdb_record_t* rec = &L->records[r];
        switch (rec->kind) {
        case PDB_LINE_MODEL:
        case PDB_LINE_ENDMDL:
            if (!L->concat_models && L->atoms) {
                if (rec->kind == PDB_LINE_MODEL && L->num_models == 0) {
                    L->model_beg = rec->atom;
                } else if (L->num_models > 0) {
                    // The end of the first model: its ENDMDL, or the MODEL record of a second without one
                    L->num_atoms = rec->atom;
                    while (L->num_res > 0 && L->res_starts[L->num_res - 1].atom >= L->num_atoms) --L->num_res;
                    L->atoms = false;
                }
            }
            if (rec->kind == PDB_LINE_MODEL) L->num_models += 1;
            break;
        case PDB_LINE_CRYST1:
        case PDB_LINE_CONECT:
        case PDB_LINE_REMARK_350:
            pdb_extras_record(&L->extras, rec->kind, (md_str_t){ L->text + rec->offset, rec->len });
            break;
        default:
            break;
        }
    }
}

// The time since *t to phase, and *t now
static inline void pdb_lap(pdb_load_t* L, int phase, md_tick_t* t) {
    const md_tick_t now = md_tick_now();
    L->ms[phase] += md_tick_to_milliseconds(now - *t);
    *t = now;
}

// A pass: false if the executor was cancelled
static bool pdb_read_pass(pdb_load_t* L, const char* text, size_t len, const md_executor_t* ex) {
    const size_t num_pieces = md_text_num_pieces(len, md_executor_num_threads(ex));
    md_array_resize(L->ranges, num_pieces, L->temp);
    md_array_resize(L->pieces, num_pieces, L->temp);
    md_text_split(L->ranges, num_pieces, text, len);
    // The index: of each piece PDB_LINES_OF its bytes, after one another. More than the load reserved only for a pass of
    // a line longer than a pass (which has room for a pass of short lines).
    if (PDB_LINES_OF(len) + num_pieces > L->lines_cap) {
        L->lines_cap = PDB_LINES_OF(len) + num_pieces;
        L->lines = md_alloc(L->temp, L->lines_cap * sizeof(pdb_line_t));
    }
    pdb_line_t* lines = L->lines;
    for (size_t p = 0; p < num_pieces; ++p) {
        L->pieces[p] = (pdb_piece_t){ .range = L->ranges[p], .types = &L->types[p], .lines = lines };
        lines += PDB_LINES_OF(L->ranges[p].end - L->ranges[p].beg);
    }
    L->text = text;

    md_tick_t t = md_tick_now();
    if (!md_parallel_for(ex, (uint32_t)num_pieces, 1, pdb_count_lines, L)) return false;
    pdb_lap(L, PDB_PHASE_COUNT, &t);
    pdb_place_pieces(L);
    pdb_lap(L, PDB_PHASE_PLACE, &t);
    if (!md_parallel_for(ex, (uint32_t)num_pieces, 1, pdb_parse_lines, L)) return false;
    pdb_lap(L, PDB_PHASE_PARSE, &t);
    pdb_gather_pieces(L);
    pdb_lap(L, PDB_PHASE_GATHER, &t);
    if (!md_parallel_for(ex, (uint32_t)num_pieces, 1, pdb_type_lines, L)) return false;
    for (size_t p = 0; p < num_pieces; ++p) {
        if (L->pieces[p].types->overflow) pdb_type_overflow(L, &L->pieces[p]);
    }
    pdb_lap(L, PDB_PHASE_TYPE, &t);
    L->num_atoms = md_array_size(L->serial);
    L->num_res = md_array_size(L->res_starts);
    pdb_read_records(L);
    pdb_lap(L, PDB_PHASE_RECORDS, &t);
    return true;
}

// Where the end of a file of several models is, from the first model's end at from: the line after the last MODEL or
// ENDMDL record, where the CONECT records follow. Looked for in the last bytes of the source, more of them until one
// is found; from when there is none in the last PDB_TAIL_MAX bytes (a file of one model, whose CONECT records follow its
// end, or of models larger than that), or the source cannot go there.
#define PDB_TAIL_MIN MD_MEGABYTES(1)
#define PDB_TAIL_MAX MD_MEGABYTES(64)

static int64_t pdb_tail(md_reader_i* src, int64_t from) {
    const int64_t size = md_reader_size(src);
    if (size < 0 || size - from <= PDB_TAIL_MIN) return from;
    for (int64_t window = PDB_TAIL_MIN; ; window *= 4) {
        const int64_t start = MD_MAX(from, size - window);
        // From the line at start or after it: the line through the byte before
        if (!md_reader_seek(src, start - 1)) return from;
        md_reader_skip_line(src);
        int64_t found = -1;
        md_str_t line;
        while (md_reader_extract_line(src, &line)) {
            const pdb_line_kind_t kind = pdb_line_kind(line);
            if (kind == PDB_LINE_MODEL || kind == PDB_LINE_ENDMDL) found = md_reader_tell(src);
        }
        if (found >= 0) return found;
        if (start == from || window >= PDB_TAIL_MAX) return from;
    }
}

// The atoms of the first model, where atoms come before its MODEL record (a file whose first model is not the first of
// its atoms): those before it are not
static void pdb_first_model(pdb_load_t* L) {
    const size_t beg = L->model_beg;
    const size_t n = L->num_atoms - beg;
    md_memmove(L->state->xyz,         L->state->xyz + beg,         n * sizeof(md_float3_t));
    md_memmove(L->sys->atom.type_idx, L->sys->atom.type_idx + beg, n * sizeof(md_particle_type_idx_t));
    md_memmove(L->sys->atom.flags,    L->sys->atom.flags + beg,    n * sizeof(md_atom_flags_t));
    md_memmove(L->serial,             L->serial + beg,             n * sizeof(int32_t));
    if (L->unassigned) md_memmove(L->atom_name, L->atom_name + beg, n * sizeof(md_label_t));
    for (int c = 0; c < PDB_NUM_COLUMNS; ++c) {
        if (L->col_varies[c]) md_memmove(L->column[c], L->column[c] + beg, n * sizeof(float));
    }
    // The residue of the first atom of the model, which may have begun before it
    pdb_residue_start_t* rs = L->res_starts;
    pdb_residue_start_t first = {0};
    for (size_t r = 0; r < L->num_res && rs[r].atom <= beg; ++r) first = rs[r];
    size_t k = 0;
    for (size_t r = 0; r < L->num_res; ++r) {
        if (rs[r].atom < beg || rs[r].atom >= L->num_atoms) continue;
        rs[k] = rs[r];
        rs[k++].atom -= (uint32_t)beg;
    }
    if (n > 0 && (k == 0 || rs[0].atom != 0)) {
        first.atom = 0;
        md_array_push(L->res_starts, first, L->temp);     // Room for one more
        rs = L->res_starts;
        md_memmove(rs + 1, rs, k * sizeof(pdb_residue_start_t));
        rs[0] = first;
        k += 1;
    }
    L->num_res = k;
    L->num_atoms = n;
}

// The columns which vary over the load, where they did not in a piece: the one value of the piece. Of the atoms read
// (a run may go past the end of the first model).
static void pdb_fill_columns(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    const pdb_load_t* L = (const pdb_load_t*)arg;
    for (uint32_t r = beg; r < end; ++r) {
        const pdb_col_run_t* run = &L->col_runs[r];
        const size_t run_end = MD_MIN(run->first_atom + run->num_atoms, L->num_atoms);
        for (int c = 0; c < PDB_NUM_COLUMNS; ++c) {
            if (!L->col_varies[c] || run->col_varies[c]) continue;
            float v;
            md_memcpy(&v, &run->col_first[c], sizeof(v));
            for (size_t i = run->first_atom; i < run_end; ++i) L->column[c][i] = v;
        }
    }
}

// The names of the atoms which have types (an element), where others have none and are typed by their names with the
// rest of their residue: the names of their types
#define PDB_NAME_CHUNK 65536

static void pdb_fill_names(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    const pdb_load_t* L = (const pdb_load_t*)arg;
    const md_particle_type_data_t* type = &L->sys->atom.type;
    const md_particle_type_idx_t* type_idx = L->sys->atom.type_idx;
    const size_t i_end = MD_MIN((size_t)end * PDB_NAME_CHUNK, L->num_atoms);
    for (size_t i = (size_t)beg * PDB_NAME_CHUNK; i < i_end; ++i) {
        if (type_idx[i] > 0) L->atom_name[i] = type->name[type_idx[i]];
    }
}

// The residues of the system from where they begin, and their chain ids as the file has them (a character, which is an
// id whatever it is; a blank is none), a chunk of PDB_RES_CHUNK at a time
#define PDB_RES_CHUNK 16384

typedef struct pdb_residues_t {
    md_system_t* sys;
    const pdb_residue_start_t* starts;
    size_t num_res;
    char*  chain_chars;
    md_str_t* chain_ids;
    bool*  after_ter;
} pdb_residues_t;

static void pdb_build_residues(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    pdb_residues_t* rb = (pdb_residues_t*)arg;
    md_residue_data_t* res = &rb->sys->residue;
    const size_t r_end = MD_MIN((size_t)end * PDB_RES_CHUNK, rb->num_res);
    for (size_t r = (size_t)beg * PDB_RES_CHUNK; r < r_end; ++r) {
        const pdb_residue_start_t* rs = &rb->starts[r];
        res->atom_offset[r] = rs->atom;
        res->name[r]   = rs->name;
        res->seq_id[r] = rs->seq_id;
        res->flags[r]  = MD_RESIDUE_FLAG_NONE;
        rb->chain_chars[r] = rs->chain_id;
        rb->chain_ids[r] = md_str_trim((md_str_t){ &rb->chain_chars[r], 1 });
        rb->after_ter[r] = rs->after_ter;
    }
}

// The system of what the load read: false if the executor was cancelled
static bool pdb_make_system(pdb_load_t* L, const md_executor_t* ex) {
    md_system_t* sys = L->sys;
    md_system_state_t* state = L->state;
    md_allocator_i* temp = L->temp;
    md_tick_t t = md_tick_now();
    bool any_column = false;
    for (int c = 0; c < PDB_NUM_COLUMNS; ++c) any_column |= L->col_varies[c];
    if (any_column && !md_parallel_for(ex, (uint32_t)md_array_size(L->col_runs), 1, pdb_fill_columns, L)) return false;
    if (L->model_beg > 0) pdb_first_model(L);
    const size_t num_atoms = L->num_atoms, num_res = L->num_res;

    if (L->num_bad_ids + L->num_bad_seq > 0) {
        MD_LOG_DEBUG("%u atom serials or residue numbers are not numbers of the PDB format ('*****' past the width of the field, "
            "or hexadecimal as VMD writes past 99999): the atoms are taken in their order", L->num_bad_ids + L->num_bad_seq);
    }

    // The atoms
    md_array_resize(state->xyz, num_atoms, state->alloc);
    state->num_atoms = num_atoms;
    sys->atom.count = num_atoms;
    // The atoms' types and flags are written as they are read (pdb_parse_lines), the padding past them here
    const size_t capacity = MD_ROUND_UP(num_atoms, 16);
    PDB_ROOM(sys->atom.type_idx, capacity, sys->alloc);
    PDB_ROOM(sys->atom.flags, capacity, sys->alloc);
    if (capacity) {
        md_memset(sys->atom.type_idx + num_atoms, 0, (capacity - num_atoms) * sizeof(md_particle_type_idx_t));
        md_memset(sys->atom.flags + num_atoms, 0, (capacity - num_atoms) * sizeof(md_atom_flags_t));
    }
    pdb_lap(L, PDB_PHASE_ATOMS, &t);

    // The residues
    md_array_resize(sys->residue.atom_offset, num_res + 1, sys->alloc);
    md_array_resize(sys->residue.name,   num_res, sys->alloc);
    md_array_resize(sys->residue.seq_id, num_res, sys->alloc);
    md_array_resize(sys->residue.flags,  num_res, sys->alloc);
    sys->residue.count = num_res;
    pdb_residues_t rb = {
        .sys = sys,
        .starts = L->res_starts,
        .num_res = num_res,
        .chain_chars = md_alloc(temp, MD_MAX(num_res, 1)),
        .chain_ids = md_alloc(temp, MD_MAX(num_res, 1) * sizeof(md_str_t)),
        .after_ter = md_alloc(temp, MD_MAX(num_res, 1) * sizeof(bool)),
    };
    if (!md_parallel_for(ex, (uint32_t)MD_DIV_UP(num_res, PDB_RES_CHUNK), 1, pdb_build_residues, &rb)) return false;
    if (L->num_bad_seq > 0) {
        for (size_t r = 0; r < num_res; ++r) {
            // Not a number of the format: the one after the residue before
            if (sys->residue.seq_id[r] == INT_MAX) sys->residue.seq_id[r] = r > 0 ? sys->residue.seq_id[r - 1] + 1 : 1;
        }
    }
    sys->residue.atom_offset[num_res] = (uint32_t)num_atoms;    // Final sentinel
    pdb_lap(L, PDB_PHASE_RESIDUES, &t);

    if (L->unassigned) {
        if (!md_parallel_for(ex, (uint32_t)MD_DIV_UP(num_atoms, PDB_NAME_CHUNK), 1, pdb_fill_names, L)) return false;
        md_system_infer_particle_types_from_labels(sys, L->atom_name, ex);
    }
    pdb_lap(L, PDB_PHASE_TYPING, &t);

    // The columns which hold more than one value: md_mmcif publishes the same, and they have to agree
    const md_str_t col_path[PDB_NUM_COLUMNS] = { MD_STR_LIT("atom/occupancy"), MD_STR_LIT("atom/b_factor"), MD_STR_LIT("atom/formal_charge") };
    const md_unit_t col_unit[PDB_NUM_COLUMNS] = { md_unit_none(), md_unit_pow(md_unit_angstrom(), 2), md_unit_none() };
    for (int c = 0; c < PDB_NUM_COLUMNS; ++c) {
        // Chemistry perception starts from the formal charges the file gives
        if (L->col_varies[c]) md_attributes_publish_atom_column(&sys->attributes, col_path[c], col_unit[c], 1, L->column[c], num_atoms);
    }

    pdb_system_finish(sys, state, &L->extras, L->serial, rb.chain_ids, rb.after_ter, temp);
    pdb_lap(L, PDB_PHASE_FINISH, &t);

    const double* ms = L->ms;
    MD_LOG_DEBUG("PDB: %zu atoms, %zu residues, in ms: read %.1f | count %.1f, place %.1f, parse %.1f, gather %.1f, type %.1f, records %.1f | "
        "atoms %.1f, residues %.1f, typing %.1f, finish %.1f", num_atoms, num_res, ms[PDB_PHASE_READ], ms[PDB_PHASE_COUNT], ms[PDB_PHASE_PLACE],
        ms[PDB_PHASE_PARSE], ms[PDB_PHASE_GATHER], ms[PDB_PHASE_TYPE], ms[PDB_PHASE_RECORDS], ms[PDB_PHASE_ATOMS], ms[PDB_PHASE_RESIDUES],
        ms[PDB_PHASE_TYPING], ms[PDB_PHASE_FINISH]);
    return true;
}

static bool pdb_load(md_system_t* sys, md_system_state_t* state, md_reader_i* src, const md_load_desc_t* desc) {
    const md_executor_t* ex = desc ? desc->executor : NULL;
    MD_ASSERT(sys);
    if (!sys->alloc) {
        MD_LOG_ERROR("System allocator not set");
        return false;
    }
    if (!state || !state->alloc) {
        MD_LOG_ERROR("State allocator not set");
        return false;
    }
    md_system_reset(sys);
    md_system_state_init(state, 0);

    // Everything of the load but what the system keeps is in one temp scope, rewound at the end. The tasks of the
    // executor take nothing of it (they allocate nothing): a scope of their thread's own, the calling thread's too, would
    // not collide with it. Arrays are reserved for as many atoms as there can be, so that none grows (and copies) in the
    // arena while another is above it: the reservation is address space, and only what is written becomes memory.
    md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp = md_temp_allocator(temp_scope);

    const uint32_t num_threads = md_executor_num_threads(ex);
    const size_t pass_size = md_text_pass_size(num_threads);
    const size_t max_pieces = md_text_num_pieces(pass_size, num_threads);

    pdb_load_t L = {
        .sys = sys,
        .state = state,
        .temp = temp,
        .concat_models = desc && (desc->flags & MD_LOAD_FLAG_ALL_MODELS),
        .atoms = true,
        .type_map = { .allocator = temp },
        .extras = { .alloc = temp },
    };

    // As many atoms as there can be, as an atom line of coordinates has 54 characters at least
    const int64_t src_size = md_reader_size(src);
    const size_t max_atoms = src_size >= 0 ? (size_t)src_size / 55 + 16 : 0;
    md_array_ensure(state->xyz, max_atoms, state->alloc);
    md_array_ensure(sys->atom.type_idx, MD_ROUND_UP(max_atoms, 16), sys->alloc);
    md_array_ensure(sys->atom.flags, MD_ROUND_UP(max_atoms, 16), sys->alloc);
    md_array_ensure(L.res_starts, max_atoms, temp);
    md_array_ensure(L.serial, max_atoms, temp);
    md_array_ensure(L.atom_name, max_atoms, temp);
    for (int c = 0; c < PDB_NUM_COLUMNS; ++c) md_array_ensure(L.column[c], max_atoms, temp);
    // A run for each piece of a pass: no more than a piece of every half of MD_TEXT_MIN_PIECE_SIZE
    md_array_ensure(L.col_runs, (src_size >= 0 ? (size_t)src_size / (MD_TEXT_MIN_PIECE_SIZE / 2) : 0) + 64, temp);
    md_array_ensure(L.pieces, max_pieces, temp);
    md_array_ensure(L.ranges, max_pieces, temp);
    L.types = md_alloc(temp, max_pieces * sizeof(pdb_types_t));
    // The index of the lines of a pass: address space for a pass of the shortest lines, of which those of atoms are written
    L.lines_cap = PDB_LINES_OF(pass_size) + max_pieces;
    L.lines = md_alloc(temp, L.lines_cap * sizeof(pdb_line_t));

    // Add default unknown particle type at index 0
    md_particle_type_find_or_add(&sys->atom.type, MD_STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, sys->alloc);

    // hy36decode fills its tables of digits on its first call, unguarded (a second thread would read them half filled):
    // once here, before the tasks call it
    {
        int unused;
        hy36decode(5, "    1", 5, &unused);
    }

    bool result = false;
    bool tail = false;      // gone to the end of the file, for its CONECT records
    for (;;) {
        md_tick_t t = md_tick_now();
        const md_str_t block = md_reader_block_parallel(src, ex);
        pdb_lap(&L, PDB_PHASE_READ, &t);
        if (block.len == 0) break;
        // The block a pass at a time: what the lines are, then what they hold, of text the first left in the caches
        size_t taken = 0;
        while (taken < block.len && (L.atoms || tail || L.concat_models)) {
            size_t len = md_text_pass_len_of(block.ptr + taken, block.len - taken, pass_size);
            // A pass whose first line is longer than a pass is all of the block: of more than 32 bits of offsets (a
            // string of one), its first line alone
            if (len > UINT32_MAX) len = MD_MIN(block.len - taken, (size_t)(md_parse_find_byte(block.ptr + taken, block.ptr + block.len, '\n') - (block.ptr + taken)) + 1);
            if (!pdb_read_pass(&L, block.ptr + taken, len, ex)) goto done;     // Cancelled
            taken += len;
        }
        md_reader_advance(src, taken);

        // The first model ended: on to the end of the file, for its CONECT records
        if (!L.atoms && !tail) {
            tail = true;
            if (!md_reader_seek(src, pdb_tail(src, md_reader_tell(src)))) {
                MD_LOG_ERROR("PDB: the source cannot be read from where the first model ends");
                goto done;
            }
        }
    }
    result = pdb_make_system(&L, ex);

done:
    // The scratch goes with the scope
    md_temp_end(temp_scope);
    if (!result) {
        md_system_reset(sys);
        md_system_state_init(state, 0);
    }
    return result;
}

#undef PDB_ROOM

bool md_pdb_system_init_from_str(md_system_t* sys, md_system_state_t* state, md_str_t str, const md_load_desc_t* desc) {
    md_reader_i src = md_reader_from_str(str);
    return pdb_load(sys, state, &src, desc);
}

bool md_pdb_system_init_from_file(md_system_t* sys, md_system_state_t* state, md_str_t filename, const md_load_desc_t* desc) {
    // The reader's buffer in a temp scope, below the one of the load (pdb_load): a few passes of it
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_file_reader_t fr;
    const size_t cap = MD_MAX((size_t)MD_TEXT_BLOCK_SIZE, 2 * md_text_pass_size(md_executor_num_threads(desc ? desc->executor : NULL)));
    if (!md_file_reader_open(&fr, filename, cap, md_temp_allocator(temp))) {
        MD_LOG_ERROR("Could not open file '" MD_STR_FMT "'", MD_STR_ARG(filename));
        md_temp_end(temp);
        return false;
    }
    // Several models are a trajectory as well: md_pdb_system_publish_run makes them a run.
    const bool result = pdb_load(sys, state, &fr.r, desc);
    md_file_reader_close(&fr);
    md_temp_end(temp);
    return result;
}

typedef struct pdb_cache_t {
    md_run_cache_header_t header;
    md_unitcell_t cell;
    int64_t* frame_offsets;
} pdb_cache_t;

// The cache beside path when it was made from the file as it is now: header, cell, num_frames + 1
// model offsets.
static bool try_read_cache(pdb_cache_t* cache, md_str_t path, md_allocator_i* alloc) {
    MD_ASSERT(cache);
    MD_ASSERT(alloc);
    md_file_t file = {0};
    if (!md_run_cache_open(&file, &cache->header, path, MD_PDB_CACHE_MAGIC, MD_PDB_CACHE_VERSION)) {
        return false;
    }
    const size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
    cache->frame_offsets = md_alloc(alloc, offset_bytes);
    const bool ok =
        md_file_read(file, &cache->cell, sizeof(cache->cell)) == sizeof(cache->cell) &&
        md_file_read(file, cache->frame_offsets, offset_bytes) == offset_bytes &&
        md_file_tell(file) == (int64_t)md_file_size(file);
    if (!ok) {
        MD_LOG_ERROR("The PDB cache beside '" MD_STR_FMT "' is incomplete", MD_STR_ARG(path));
        md_free(alloc, cache->frame_offsets, offset_bytes);
        cache->frame_offsets = NULL;
    }
    md_file_close(&file);
    return ok;
}

static bool write_cache(const pdb_cache_t* cache, md_str_t path, const md_file_info_t* scanned) {
    md_file_t file = {0};
    if (!md_run_cache_create(&file, path, scanned, MD_PDB_CACHE_MAGIC, MD_PDB_CACHE_VERSION, cache->header.num_atoms, cache->header.num_frames)) {
        return false;
    }
    const size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
    const bool ok =
        md_file_write(file, &cache->cell, sizeof(cache->cell)) == sizeof(cache->cell) &&
        md_file_write(file, cache->frame_offsets, offset_bytes) == offset_bytes;
    if (!ok) {
        MD_LOG_ERROR("Failed to write the PDB cache beside '" MD_STR_FMT "'", MD_STR_ARG(path));
    }
    md_file_close(&file);
    return ok;
}

// Where each model is and the file's cell, from the cache beside the file when that is current and from a scan of the
// file otherwise (writing the cache unless told not to). The scan reads the names of the records, and of the lines
// nothing else but the first CRYST1's: a model begins at its MODEL record, and its atoms are those the system has of
// it (the first alternate locations, pdb_line_kind), up to its ENDMDL record or the next MODEL record.
static bool pdb_index_load(pdb_cache_t* cache, md_str_t filename, md_run_flags_t flags, md_allocator_i* alloc) {
    md_memset(cache, 0, sizeof(*cache));

    if (try_read_cache(cache, filename, alloc)) {
        return true;
    }
    md_memset(cache, 0, sizeof(*cache));

    // The file as it is before the scan is what the cache is stamped with
    md_file_info_t scanned = { 0 };
    if (!md_file_info_extract_from_path(filename, &scanned)) {
        MD_LOG_ERROR("Failed to extract file info from path '" MD_STR_FMT "'", MD_STR_ARG(filename));
        return false;
    }

    md_allocator_i* heap = md_get_heap_allocator();
    md_file_reader_t fr;
    if (!md_file_reader_open(&fr, filename, MD_READER_DEFAULT_CAPACITY, heap)) {
        MD_LOG_ERROR("Could not open file '" MD_STR_FMT "'", MD_STR_ARG(filename));
        return false;
    }

    md_array(int64_t)  model_offset = 0;
    md_array(uint32_t) model_atoms = 0;
    pdb_extras_t extras = { .alloc = heap };    // The cell
    bool in_model = false;
    md_str_t line;
    for (;;) {
        const int64_t at = md_reader_tell(&fr.r);   // The offset of the line
        if (!md_reader_extract_line(&fr.r, &line)) break;
        const pdb_line_kind_t kind = pdb_line_kind(line);
        switch (kind) {
        case PDB_LINE_ATOM:
            if (in_model) *md_array_last(model_atoms) += 1;
            break;
        case PDB_LINE_MODEL:
            md_array_push(model_offset, at, heap);
            md_array_push(model_atoms, 0, heap);
            in_model = true;
            break;
        case PDB_LINE_ENDMDL:
            in_model = false;
            break;
        case PDB_LINE_CRYST1:
            pdb_extras_record(&extras, kind, line);
            break;
        default:
            break;
        }
    }
    md_file_reader_close(&fr);

    bool success = false;
    const size_t num_models = md_array_size(model_offset);
    if (num_models <= 1) {
        MD_LOG_INFO("The PDB file does not contain multiple model entries and cannot be read as a trajectory");
        goto done;
    }
    // Models of the same atoms (an empty one aside)
    const size_t num_atoms = model_atoms[0];
    if (!num_atoms) {
        MD_LOG_ERROR("The PDB file models are empty and cannot be read as a trajectory");
        goto done;
    }
    for (size_t i = 1; i < num_models; ++i) {
        if (model_atoms[i] && model_atoms[i] != num_atoms) {
            MD_LOG_ERROR("The PDB file models are not of equal length and cannot be read as a trajectory");
            goto done;
        }
    }

    cache->header.num_atoms = num_atoms;
    cache->header.num_frames = num_models;
    cache->frame_offsets = md_alloc(alloc, (num_models + 1) * sizeof(int64_t));
    md_memcpy(cache->frame_offsets, model_offset, num_models * sizeof(int64_t));
    cache->frame_offsets[num_models] = (int64_t)scanned.size;
    if (extras.has_cell) {
        cache->cell = extras.cell;
    }

    if (!(flags & MD_RUN_FLAG_DISABLE_CACHE_WRITE)) {
        // A cache that cannot be written only costs the next load a scan
        write_cache(cache, filename, &scanned);
    }
    success = true;

done:
    md_array_free(model_offset, heap);
    md_array_free(model_atoms, heap);
    return success;
}

// ### RUN ###

// The coordinates of atoms [first, first + count) of one model, packed xyz, straight off its ATOM and HETATM records:
// those the system has (the first alternate locations), read as the load reads them. Returns how many were written.
static size_t pdb_parse_model_coords(float* xyz, md_str_t str, size_t first, size_t count) {
    md_str_t line;
    size_t i = 0;
    size_t written = 0;
    while (written < count && md_parse_next_line(&line, &str)) {
        if (pdb_line_kind(line) != PDB_LINE_ATOM) continue;
        if (i >= first) {
            xyz[written * 3 + 0] = extract_fixed_float(line, 31, 38, 3);
            xyz[written * 3 + 1] = extract_fixed_float(line, 39, 46, 3);
            xyz[written * 3 + 2] = extract_fixed_float(line, 47, 54, 3);
            written += 1;
        }
        i += 1;
    }
    return written;
}

// <run>/atom/position: one read of the model's text, parsed into the caller's buffer. One atom asks
// for the text up to that atom only.
static size_t pdb_position_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    const md_system_t* sys = (const md_system_t*)user_data;
    MD_ASSERT(sys);
    if (!slice || slice->num_idx == 0 || slice->num_idx > 2) return 0;

    md_run_source_t src;
    if (!md_run_source(&src, &sys->attributes, attr, MD_STR_LIT("atom/position"))) return 0;

    const uint32_t frame = slice->idx[0];
    const size_t N = attr->format.shape[1];
    size_t first = 0, count = N;
    if (slice->num_idx == 2) {
        if (slice->idx[1] >= N) return 0;
        first = slice->idx[1];
        count = 1;
    }
    if (cap != count * 3) return 0;

    const size_t frame_size = (size_t)src.size[frame];
    md_temp_scope_t temp = md_temp_begin();
    size_t written = 0;
    char* text = md_temp_alloc(temp, MD_MAX(frame_size, 1));
    if (text && md_attribute_io_read_at(io, src.path, src.offset[frame], text, frame_size) == frame_size) {
        if (pdb_parse_model_coords((float*)dst, (md_str_t){ text, frame_size }, first, count) == count) {
            written = cap;
        } else {
            MD_LOG_ERROR("PDB: model %u of '" MD_STR_FMT "' has fewer than %zu atoms", frame, MD_STR_ARG(src.path), first + count);
        }
    } else {
        MD_LOG_ERROR("PDB: Failed to read model %u from '" MD_STR_FMT "'", frame, MD_STR_ARG(src.path));
    }
    md_temp_end(temp);
    return written;
}

bool md_pdb_system_publish_run(md_system_t* sys, md_str_t filename, md_str_t run, uint32_t flags) {
    MD_ASSERT(sys);
    char path_buf[4096];
    const size_t path_len = md_path_write_canonical(path_buf, sizeof(path_buf), filename);
    const md_str_t path = { path_buf, path_len };

    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    bool result = false;

    pdb_cache_t index;
    if (path_len == 0 || !pdb_index_load(&index, path, flags, arena)) {
        goto done;
    }

    const size_t F = index.header.num_frames;
    const size_t N = index.header.num_atoms;
    double*  times = md_alloc(arena, F * sizeof(double));
    int64_t* sizes = md_alloc(arena, F * sizeof(int64_t));
    float*   boxes = md_alloc(arena, F * 9 * sizeof(float));

    // One CRYST1 for the file, so one cell for every model; zero without one
    float A[3][3] = {0};
    if (index.cell.flags) {
        md_unitcell_A_extract_float(A, &index.cell);
    }
    for (size_t i = 0; i < F; ++i) {
        times[i] = (double)i;
        sizes[i] = index.frame_offsets[i + 1] - index.frame_offsets[i];
        md_memcpy(boxes + i * 9, A, sizeof(A));
    }

    const md_attribute_virtual_t virt = { .provider = pdb_position_provider, .user_data = sys };
    const md_run_desc_t desc = {
        .num_frames    = F,
        .num_atoms     = N,
        .time          = times,
        .time_unit     = md_unit_none(),    // a model has no time: ordinals
        .unitcell      = boxes,
        .source_path   = path,
        .source_offset = index.frame_offsets,
        .source_size   = sizes,
        .position_virt = &virt,
    };
    result = md_run_publish(sys, run, &desc);

done:
    md_arena_allocator_destroy(arena);
    return result;
}

bool md_pdb_system_write_state_to_file(md_file_t file, const struct md_system_t* sys, const md_system_state_t* state, const int32_t* atom_indices, size_t num_atoms, int model_num) {
    if (model_num > 0) {
        md_file_printf(file, "MODEL     %4d\n", model_num);
    }

	md_part_idx_t prev_part_idx = INT_MIN;
    for (size_t i = 0; i < num_atoms; ++i) {
        int idx = atom_indices ? atom_indices[i] : (int)i;

        char atom_id[6] = "";
        {
            const char* errmsg = hy36encode(5, idx + 1, atom_id);
            if (errmsg) {
                MD_LOG_DEBUG("Failed to encode atom index '%d', error: %s. This atom will be written with a default atom id of 1", idx + 1, errmsg);
            }
        }
        
        // Get element symbol
		md_atomic_number_t atomic_nr = md_atom_atomic_number(&sys->atom, idx);

		char element[3] = "";
        md_str_copy_to_char_buf(element, sizeof(element), md_atomic_number_symbol(atomic_nr));
        
        // Get atom name (use element symbol if not available)
		char name[5] = "";
		md_str_copy_to_char_buf(name, sizeof(name), md_atom_name(&sys->atom, idx));
        
        // Get residue name
		md_residue_idx_t res_idx = md_system_residue_find_by_atom_idx(sys, idx);
		char resname[5] = "";
		md_str_copy_to_char_buf(resname, sizeof(resname), md_system_residue_name(sys, res_idx));

		md_sequence_id_t res_seq = md_system_residue_seq_id(sys, res_idx);
        
        // The chain id, the label the atom is filed under (blank for one in no chain)
        char chain_id[4] = " ";
        const md_chain_idx_t chain_idx = md_system_chain_find_by_atom_idx(sys, idx);
        if (chain_idx >= 0) md_str_copy_to_char_buf(chain_id, sizeof(chain_id), md_system_chain_id(sys, (size_t)chain_idx));
        const md_part_idx_t part_idx = md_system_part_find_by_atom_idx(sys, idx);

        char res_seq_str[5] = "";
        {
            const char* errmsg = hy36encode(4, res_seq, res_seq_str);
            if (errmsg) {
                MD_LOG_DEBUG("Failed to encode residue sequence id '%d' for atom index %d, error: %s.", res_seq, idx, errmsg);
            }
        }

        // A TER ends each part (a strand, a ligand, a water) before the next one begins
        if (prev_part_idx != INT_MIN && part_idx != prev_part_idx) {
            md_file_printf(file, "TER\n");
        }
        prev_part_idx = part_idx;

        // PDB format:
        // ATOM serial name altLoc resName chainID resSeq iCode x y z occupancy tempFactor element charge
        md_file_printf(file,
            "%-6s%5s %-4s%1c%3s %1c%4s%1c   "
            "%8.3f%8.3f%8.3f"
            "%6.2f%6.2f          "
            "%2s%2s\n",
            "ATOM",
            atom_id,            // serial number
            name,               // atom name
            ' ',                // altLoc
            resname,            // residue name
            chain_id[0],        // chain ID
            res_seq_str,        // residue sequence number
            ' ',                // iCode
            state->xyz[idx].x,      // x coordinate
            state->xyz[idx].y,      // y coordinate
            state->xyz[idx].z,      // z coordinate
            1.0,                // occupancy
            0.0,                // bfactor
            element,            // element symbol
            ""                  // charge
        );


    }
    
    if (model_num > 0) {
        md_file_printf(file, "ENDMDL\n");
    } else {
        md_file_printf(file, "END\n");
    }

    return true;
}

#ifdef __cplusplus
}
#endif
