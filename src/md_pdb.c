#include <md_pdb.h>

#include <md_system.h>
#include <md_util.h>

#include <core/md_common.h>
#include <core/md_array.h>
#include <core/md_str.h>
#include <core/md_parse.h>
#include <core/md_arena_allocator.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_vec_math.h>
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

static inline float extract_float(str_t line, size_t beg, size_t end) {
    if (line.len < end) return 0.0f;
    return (float)parse_float(str_substr(line, beg - 1, end-beg + 1));
}

// A column the format writes with a fixed number of decimals ("%8.3f" for coordinates, "%6.2f" for
// occupancy and temperature factor): read without a scan when it is laid out that way, and parsed like
// any other column when it is not
static inline float extract_fixed_float(str_t line, size_t beg, size_t end, size_t decimals) {
    if (line.len < end) return 0.0f;
    const str_t field = str_substr(line, beg - 1, end-beg + 1);
    float value;
    return md_parse_fixed_f32(&value, field, decimals) ? value : (float)parse_float(field);
}

static inline char extract_char(str_t line, size_t idx) {
    if (line.len < idx) return 0;
    return line.ptr[idx - 1];
}

// An atom serial or a residue sequence number: a decimal number, past the width of its field hybrid-36 (A0000, then
// a0000: a letter first). 0 where the field holds neither, counted in num_bad: "*****" past the width, which many
// writers put there, or the hexadecimal VMD writes past 99999 ("186a0" for 100000), which is not hybrid-36 (a number
// which begins with a digit is decimal) and cannot be told from decimal where it has no letter (18700, 100096). Such a
// file is read by the order of its atoms; only CONECT records refer to atoms by their serials.
static inline int32_t parse_id(str_t line, size_t beg, size_t end, uint32_t* num_bad) {
	uint32_t len = (uint32_t)(end - beg) + 1;
    int32_t  res = 0;
    if (hy36decode(len, line.ptr + beg - 1, len, &res)) {
        if (num_bad) *num_bad += 1;
        return 0;
    }
    return res;
}

// The formal charge of the charge column (79 to 80): a digit and its sign, "1-" or "2+" (and the sign first, as some
// writers have it). NAN where there is none, which chemistry perception takes as not given (md_chem_perceive).
static inline float pdb_charge(str_t field) {
    field = str_trim(field);
    if (field.len == 1 && field.ptr[0] == '0') return 0.0f;
    if (field.len != 2) return NAN;
    const char a = field.ptr[0], b = field.ptr[1];
    if (is_digit(a) && (b == '+' || b == '-')) return (float)((b == '-' ? -1 : 1) * (a - '0'));
    if (is_digit(b) && (a == '+' || a == '-')) return (float)((a == '-' ? -1 : 1) * (b - '0'));
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

// The record of the cell: CRYST1, and CRYST with a blank where the 1 is, as some tools write it (the Conformational
// Analysis Tools of CHARMM-GUI's trajectories, with all angles 0)
static inline bool is_cryst1_record(str_t line) {
    return str_eq_cstr_n(line, "CRYST", 5) && (line.len == 5 || line.ptr[5] == '1' || line.ptr[5] == ' ');
}

static inline pdb_line_kind_t pdb_line_kind(str_t line) {
    // A TER record may be the three letters alone (GROMACS, and the EM structures of the PDB written by some tools)
    if (str_eq_cstr_n(line, "TER", 3)) return PDB_LINE_TER;
    if (line.len < 6) return PDB_LINE_NONE;
    if (str_eq_cstr_n(line, "ATOM", 4) || str_eq_cstr_n(line, "HETATM", 6)) {
        // The first alternate location only: blank or A
        const char alt_loc = extract_char(line, 17);
        return (alt_loc == ' ' || alt_loc == 'A') ? PDB_LINE_ATOM : PDB_LINE_ATOM_SKIPPED;
    }
    if (str_eq_cstr_n(line, "MODEL", 5))        return PDB_LINE_MODEL;
    if (str_eq_cstr_n(line, "ENDMDL", 6))       return PDB_LINE_ENDMDL;
    if (str_eq_cstr_n(line, "CONECT", 6))       return PDB_LINE_CONECT;
    if (is_cryst1_record(line))                 return PDB_LINE_CRYST1;
    if (str_eq_cstr_n(line, "REMARK 350", 10))  return PDB_LINE_REMARK_350;
    return PDB_LINE_NONE;
}

// The cell of a CRYST1 record. No angles at all is a rectangular box whose angles were left out, not a degenerate cell.
// False for a cell of 1 x 1 x 1 at right angles, which writers put there for a structure of no cell (no periodic
// boundaries).
static bool pdb_cryst1_cell(md_unitcell_t* cell, str_t line) {
    float a = extract_float(line, 7, 15), b = extract_float(line, 16, 24), c = extract_float(line, 25, 33);
    float alpha = extract_float(line, 34, 40), beta = extract_float(line, 41, 47), gamma = extract_float(line, 48, 54);
    if (alpha == 0 && beta == 0 && gamma == 0) {
        alpha = beta = gamma = 90;
    }
    if (a == 1 && b == 1 && c == 1 && alpha == 90 && beta == 90 && gamma == 90) {
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
    md_array(mat4_t) transforms;
    md_allocator_i*  alloc;
} pdb_extras_t;

// A serial of a CONECT record, 0 where there is none. Serials are hybrid-36 as those of the atoms are, past 99999 atoms.
static inline int32_t conect_serial(str_t line, size_t beg) {
    const size_t end = beg + 4;
    if (line.len < end) return 0;
    if (str_empty(str_trim(str_substr(line, beg - 1, 5)))) return 0;
    return parse_id(line, beg, end, NULL);
}

static void pdb_extras_record(pdb_extras_t* x, pdb_line_kind_t kind, str_t line) {
    md_allocator_i* alloc = x->alloc;
    switch (kind) {
    case PDB_LINE_CRYST1:
        // The first one: a file of models has one for each, the same
        if (!x->seen_cryst1) {
            x->seen_cryst1 = true;
            x->has_cell = pdb_cryst1_cell(&x->cell, str_trim(line));
        }
        break;
    case PDB_LINE_CONECT: {
        // The atom of columns 7 to 11 and up to four atoms bonded to it, each a bond. A record of more is continued on
        // the next line, which names the atom again.
        line = str_trim(line);
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
        const str_t rest = str_substr(line, 11, SIZE_MAX);
        if (str_eq_cstr_n(rest, "BIOMOLECULE:", 12)) {
            x->biomolecule = (int32_t)parse_int(str_trim(str_substr(rest, 12, SIZE_MAX)));
        }
        else if (str_eq_cstr_n(rest, "APPLY THE FOLLOWING TO CHAINS:", 30) || str_eq_cstr_n(str_trim(rest), "AND CHAINS:", 11)) {
            // A group begins with the chains to apply to, which AND CHAINS lines continue
            if (str_eq_cstr_n(rest, "APPLY", 5)) {
                const pdb_group_t group = {
                    .id = x->biomolecule,
                    .chain_offset = (uint32_t)md_array_size(x->chains),
                    .transform_offset = (uint32_t)md_array_size(x->transforms),
                };
                md_array_push(x->groups, group, alloc);
            }
            pdb_group_t* group = md_array_last(x->groups);
            if (group) {
                size_t colon = 0;
                str_find_char(&colon, rest, ':');
                str_t chains = str_substr(rest, colon + 1, SIZE_MAX);
                str_t chain;
                while (extract_token_delim(&chain, &chains, ',')) {
                    chain = str_trim(chain);
                    if (str_empty(chain)) continue;
                    md_array_push(x->chains, chain.ptr[0], alloc);
                    group->chain_count += 1;
                }
            } else {
                MD_LOG_DEBUG("REMARK 350: chains of no biomolecule");
            }
        }
        else if (str_eq_cstr_n(str_trim(rest), "BIOMT", 5) && line.len > 18) {
            // REMARK 350   BIOMTr   n  m[r][1] m[r][2] m[r][3]  t[r], the rows of transform n of the group
            pdb_group_t* group = md_array_last(x->groups);
            const int32_t row = char_to_digit(line.ptr[18]) - 1;
            str_t fields = str_substr(line, 19, SIZE_MAX);
            str_t tok[5];
            const size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &fields);
            if (!group || row < 0 || row > 2 || num_tok < 5) {
                MD_LOG_DEBUG("REMARK 350: a malformed BIOMT record");
            } else {
                if (row == 0) {
                    md_array_push(x->transforms, mat4_ident(), alloc);
                    group->transform_count += 1;
                }
                mat4_t* transform = md_array_last(x->transforms);
                for (int c = 0; c < 3; ++c) {
                    transform->elem[c][row] = (float)parse_float(tok[1 + c]);
                }
                transform->elem[3][row] = (float)parse_float(tok[4]);
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

    md_array(mat4_t)     transform = 0;
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
                const str_t auth = md_system_instance_auth_id(sys, i);
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
            const mat4_t M = x->transforms[group->transform_offset + t];
            size_t k = 0;
            for (; k < md_array_size(transform); ++k) {
                if (MEMCMP(&transform[k], &M, sizeof(mat4_t)) == 0) break;
            }
            if (k == md_array_size(transform)) {
                char buf[16];
                const int len = snprintf(buf, sizeof(buf), "%u", t + 1);
                md_array_push(transform, M, temp_alloc);
                md_array_push(oper_id, make_label((str_t){buf, (size_t)MAX(len, 0)}), temp_alloc);
            }
            md_array_push(oper, (uint32_t)k, temp_alloc);
        }

        char buf[16];
        const int len = snprintf(buf, sizeof(buf), "%d", group->id);
        const md_assembly_gen_t row = {
            .assembly_id = make_label((str_t){buf, (size_t)MAX(len, 0)}),
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
// neither the bonds nor what the residues are by them: that is md_util_system_infer's, which a caller runs when (and
// where, on what thread) it wants.

// What follows the atoms and residues: the cell, the bonds of the CONECT records (by the serials of the atoms loaded),
// the units of the file (the runs of residues of one chain id) and the assemblies of REMARK 350
static void pdb_system_finish(md_system_t* sys, md_system_state_t* state, const pdb_extras_t* x, const int32_t* atom_serial, const str_t* res_chain_ids, const bool* res_after_ter, md_allocator_i* temp_arena) {
    if (x->has_cell) {
        state->unitcell = x->cell;
    }

    // The bonds of the CONECT records, of origin MD_BOND_ORIGIN_FILE: inference keeps them and adds the bonds they leave
    // out (CONECT records usually cover the ligands only, in a coarse grained structure from martinize2 all of it).
    // Their serials name atoms of the first model; alternate locations other than the first are not loaded.
    const size_t num_conect = md_array_size(x->conect) / 2;
    if (num_conect > 0 && !atom_serial) {
        MD_LOG_INFO("PDB: CONECT records which are not at the end of the file, where the format has them, are not read");
    }
    if (num_conect > 0 && atom_serial) {
        md_hashmap32_t serial_to_atom = { .allocator = temp_arena };
        md_hashmap_reserve(&serial_to_atom, sys->atom.count);
        for (size_t i = 0; i < sys->atom.count; ++i) {
            const uint64_t key = (uint64_t)(uint32_t)atom_serial[i];
            // 0: a serial which is not a number of the format (parse_id), which no record can refer to
            if (key && !md_hashmap_get(&serial_to_atom, key)) md_hashmap_add(&serial_to_atom, key, (uint32_t)i);
        }

        // Each bond is in the records of both its atoms, and some writers repeat it for its bond order
        md_hashmap32_t seen = { .allocator = temp_arena };
        md_hashmap_reserve(&seen, num_conect);
        const md_bond_flags_t file = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_FILE);
        for (size_t i = 0; i < num_conect; ++i) {
            const uint32_t* a = md_hashmap_get(&serial_to_atom, (uint64_t)(uint32_t)x->conect[2 * i + 0]);
            const uint32_t* b = md_hashmap_get(&serial_to_atom, (uint64_t)(uint32_t)x->conect[2 * i + 1]);
            if (!a || !b || *a == *b) continue;
            const uint32_t lo = MIN(*a, *b), hi = MAX(*a, *b);
            const uint64_t pair_key = ((uint64_t)lo << 32) | hi;
            if (md_hashmap_get(&seen, pair_key)) continue;
            md_hashmap_add(&seen, pair_key, 1);
            md_atom_pair_t pair = {{ (md_atom_idx_t)lo, (md_atom_idx_t)hi }};
            md_array_push(sys->bond.pairs, pair, sys->alloc);
            md_array_push(sys->bond.flags, file, sys->alloc);
            sys->bond.count += 1;
        }
        if (sys->bond.count) {
            md_bond_build_connectivity(&sys->bond, sys->atom.count, sys->alloc);
        }
    }

    // The units of the file are its chains: the runs of residues of one chain id, which a TER record ends (two chains of
    // one id written one after the other are two: the segments of an EM structure)
    md_util_system_instances_from_chain_ids(sys, res_chain_ids, res_after_ter);

    // The assemblies name the chains by the ids they have in the file, which are those of the instances
    pdb_assemblies(sys, x, temp_arena);
}

// ---------------------------------------------------------------------------------------------------------------------
// Straight from the text, a block at a time, each block in pieces at line boundaries on as many threads as the executor
// has (md_reader.h: blocks and pieces). A piece counts the atoms it keeps (an alternate location other than the first is not) and the
// other records of interest first, which places them; then it reads its atoms into arrays where they go, and the
// residues where the residue fields of an atom differ from those of the atom before (or a TER record is between them).
// The other records are few: they are read in order, on the calling thread, once the pieces are done.
//
// Of a file of several models the atoms of the first are loaded (all of them, MD_PDB_OPTION_CONCAT_MODELS): once it
// ends, the rest of the file is only looked through for the records which are not of atoms (CONECT, which come last).

#define PDB_RES_FIELDS_BEG  17      // The residue name, the chain id, the sequence number and the insertion code:
#define PDB_RES_FIELDS_END  27      // columns 18 to 27
#define PDB_NO_ELEMENT      0xFF

// The residue fields of an atom line, the chain id as the loader takes it (a letter, or none)
typedef struct pdb_res_fields_t {
    char c[PDB_RES_FIELDS_END - PDB_RES_FIELDS_BEG];
} pdb_res_fields_t;

static inline pdb_res_fields_t pdb_res_fields(str_t line) {
    pdb_res_fields_t f;
    MEMSET(f.c, ' ', sizeof(f.c));
    if (line.len > PDB_RES_FIELDS_BEG) {
        MEMCPY(f.c, line.ptr + PDB_RES_FIELDS_BEG, MIN(line.len - PDB_RES_FIELDS_BEG, sizeof(f.c)));
    }
    char* chain = &f.c[21 - PDB_RES_FIELDS_BEG];
    if (!is_alpha(*chain)) *chain = ' ';
    return f;
}

typedef struct pdb_residue_start_t {
    uint32_t         atom;      // of the load
    md_sequence_id_t seq_id;
    md_label_t       name;
    char             chain_id;
    bool             after_ter;     // A TER record before it: it begins a unit of the file
} pdb_residue_start_t;

typedef struct pdb_record_t {
    size_t          offset;     // of the line in the block
    size_t          len;
    size_t          atom;       // the atoms of the load before it
    pdb_line_kind_t kind;
} pdb_record_t;

typedef struct pdb_piece_t {
    md_text_piece_t range;
    size_t   num_atoms;         // kept
    size_t   num_records;
    size_t   first_atom;        // of the block
    size_t   first_record;
    size_t   num_residues;      // found, at res_starts + first_atom
    uint32_t num_bad_ids;
    bool     unassigned;        // an atom without an element
} pdb_piece_t;

enum { PDB_COL_OCCUPANCY, PDB_COL_B_FACTOR, PDB_COL_CHARGE, PDB_NUM_COLUMNS };

typedef struct pdb_block_t {
    const char*           text;
    size_t                len;
    pdb_piece_t*          pieces;
    size_t                atom_beg;     // of the load, the block's first atom
    bool                  atoms;        // false once the first model ended: the atom lines are passed over
    bool                  has_prev;     // an atom before the block
    pdb_res_fields_t      prev;         // its residue fields
    bool                  prev_ter;     // a TER record after it, before the block
    pdb_residue_start_t*  res_starts;   // [atoms of the block]
    pdb_record_t*         records;      // [records of the block]
    vec3_t*               xyz;          // Of the load, by atom: the state's
    int32_t*              serial;       // Of the load, by atom, where CONECT records need them (pdb_has_conect): else NULL
    // Of the block, by its atom: what is resolved once the block is read
    md_label_t*           atom_name;
    uint8_t*              element;      // PDB_NO_ELEMENT when the line has none
    float*                column[PDB_NUM_COLUMNS];  // occupancy, b factor, formal charge
} pdb_block_t;

static void pdb_count_lines(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    pdb_block_t* b = (pdb_block_t*)arg;
    for (uint32_t p = beg; p < end; ++p) {
        pdb_piece_t* piece = &b->pieces[p];
        size_t num_atoms = 0, num_records = 0;
        const char* c = b->text + piece->range.beg;
        const char* e = b->text + piece->range.end;
        while (c < e) {
            const pdb_line_kind_t kind = pdb_line_kind(md_text_next_line(&c, e));
            if (kind == PDB_LINE_ATOM) num_atoms += b->atoms ? 1 : 0;
            else if (kind >= PDB_LINE_TER) num_records += 1;
        }
        piece->num_atoms = num_atoms;
        piece->num_records = num_records;
    }
}

static void pdb_parse_lines(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    pdb_block_t* b = (pdb_block_t*)arg;
    for (uint32_t p = beg; p < end; ++p) {
        pdb_piece_t* piece = &b->pieces[p];
        piece->num_residues = 0;
        piece->num_bad_ids = 0;
        piece->unassigned = false;

        // The atom before the piece, and whether a TER record is between: in the block before it, or the block before
        bool has_prev = false, ter = false;
        pdb_res_fields_t prev = {0};
        if (b->atoms && piece->num_atoms > 0) {
            size_t at = piece->range.beg;
            while (at > 0) {
                const str_t line = md_text_line_before(b->text, at);
                at = (size_t)(line.ptr - b->text);
                const pdb_line_kind_t kind = pdb_line_kind(line);
                if (kind == PDB_LINE_TER) ter = true;
                else if (kind == PDB_LINE_ATOM) {
                    has_prev = true;
                    prev = pdb_res_fields(line);
                    break;
                }
            }
            if (!has_prev && b->has_prev) {
                has_prev = true;
                prev = b->prev;
                ter |= b->prev_ter;
            }
        }

        pdb_residue_start_t* res = b->res_starts + piece->first_atom;
        pdb_record_t* rec = b->records + piece->first_record;
        size_t atom = b->atom_beg + piece->first_atom;
        const char* c = b->text + piece->range.beg;
        const char* e = b->text + piece->range.end;
        while (c < e) {
            const str_t line = md_text_next_line(&c, e);
            const pdb_line_kind_t kind = pdb_line_kind(line);
            if (kind == PDB_LINE_ATOM && b->atoms) {
                // A residue is its name, its sequence number and its insertion code in its chain: 100, 100A and 100B
                // are three residues of the same number (the insertions of Kabat numbered antibodies), also when they
                // share a name. A TER record ends the residue before it, whatever the residue after it is called, and
                // the chain: what follows is a unit of its own (md_util_system_instances_from_chain_ids).
                const pdb_res_fields_t fields = pdb_res_fields(line);
                if (!has_prev || ter || MEMCMP(fields.c, prev.c, sizeof(fields.c)) != 0) {
                    res->atom     = (uint32_t)atom;
                    res->seq_id   = line.len >= 26 ? parse_id(line, 23, 26, &piece->num_bad_ids) : 0;
                    res->name     = make_label(str_trim(str_substr(line, 17, 4)));
                    res->chain_id = fields.c[21 - PDB_RES_FIELDS_BEG];
                    res->after_ter = has_prev && ter;
                    ++res;
                }
                has_prev = true;
                ter = false;
                prev = fields;

                b->xyz[atom] = vec3_set(extract_fixed_float(line, 31, 38, 3), extract_fixed_float(line, 39, 46, 3), extract_fixed_float(line, 47, 54, 3));
                const size_t k = atom - b->atom_beg;
                b->atom_name[k] = make_label(str_trim(str_substr(line, 12, 4)));
                const str_t element = str_trim(str_substr(line, 76, 2));
                if (element.len) {
                    b->element[k] = (uint8_t)md_atomic_number_from_symbol(element, true);
                } else {
                    b->element[k] = PDB_NO_ELEMENT;
                    piece->unassigned = true;
                }
                b->column[PDB_COL_OCCUPANCY][k] = extract_fixed_float(line, 55, 60, 2);
                b->column[PDB_COL_B_FACTOR][k]  = extract_fixed_float(line, 61, 66, 2);
                b->column[PDB_COL_CHARGE][k]    = pdb_charge(str_substr(line, 78, 2));
                if (b->serial) {
                    b->serial[atom] = line.len >= 11 ? parse_id(line, 7, 11, &piece->num_bad_ids) : 0;
                }
                ++atom;
            } else if (kind >= PDB_LINE_TER) {
                if (kind == PDB_LINE_TER) ter = true;
                rec->offset = (size_t)(line.ptr - b->text);
                rec->len    = line.len;
                rec->atom   = atom;
                rec->kind   = kind;
                ++rec;
            }
        }
        piece->num_residues = (size_t)(res - (b->res_starts + piece->first_atom));
    }
}

// A column of one value for every atom (the occupancy of 1 and the b factor of 0 of a file from a simulation, no formal
// charge at all) is that value: an array is made once an atom differs, holding the value for the atoms before it. Bits
// are compared, as md_attributes_publish_atom_column does to publish no column of one value (NAN is none).
typedef struct pdb_column_t {
    bool     seen;
    uint32_t first;             // The bits of the first atom's value
    md_array(float) data;       // NULL while every value is the first
    size_t   reserve;           // Of the array, once made: for as many atoms as there can be
} pdb_column_t;

static inline uint32_t pdb_f32_bits(float x) {
    uint32_t u;
    MEMCPY(&u, &x, sizeof(u));
    return u;
}

// The values of count atoms after the beg before them
static void pdb_column_append(pdb_column_t* col, const float* values, size_t beg, size_t count, md_allocator_i* alloc) {
    if (count == 0) return;
    if (!col->data) {
        if (!col->seen) {
            col->seen = true;
            col->first = pdb_f32_bits(values[0]);
        }
        size_t i = 0;
        while (i < count && pdb_f32_bits(values[i]) == col->first) ++i;
        if (i == count) return;
        float first;
        MEMCPY(&first, &col->first, sizeof(first));
        md_array_ensure(col->data, MAX(beg + count, col->reserve), alloc);
        md_array_resize(col->data, beg, alloc);
        for (size_t j = 0; j < beg; ++j) col->data[j] = first;
    }
    md_array_push_array(col->data, values, count, alloc);
}

// Room for count elements in an array of the load, which grows by doubling past what was reserved
#define PDB_ROOM(arr, count, alloc) do { \
        if ((count) > md_array_capacity(arr)) md_array_ensure((arr), MAX((size_t)(count), 2 * md_array_capacity(arr)), (alloc)); \
        md_array_resize((arr), (count), (alloc)); \
    } while (0)

// The CONECT records of a file are at its end (after the coordinates, before MASTER and END): when the end has none,
// the serials of the atoms, which only they refer to, are not kept. A source which cannot be looked at the end of (a
// stream) is taken to have them.
#define PDB_CONECT_PROBE KILOBYTES(256)

static bool pdb_has_conect(md_reader_i* src) {
    const int64_t size = md_reader_size(src);
    if (size < 0) return true;
    const int64_t from = MAX(0, size - (int64_t)PDB_CONECT_PROBE);
    if (!md_reader_seek(src, from)) return true;
    bool found = false;
    str_t line;
    if (from > 0) md_reader_skip_line(src);     // The rest of a line which began before
    while (!found && md_reader_extract_line(src, &line)) {
        found = str_eq_cstr_n(line, "CONECT", 6);
    }
    if (!md_reader_rewind(src)) {
        MD_LOG_ERROR("PDB: the source cannot be read from its start again");
    }
    return found;
}

static bool pdb_load(md_system_t* sys, md_system_state_t* state, md_reader_i* src, md_pdb_options_t options, const md_executor_t* ex) {
    ASSERT(sys);
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
    // executor take nothing of it (pdb_count_lines, pdb_parse_lines allocate nothing): a scope of their thread's own,
    // the calling thread's too, would not collide with it. Arrays are reserved for as many atoms as there can be, so
    // that none grows (and copies) in the arena while another is above it: the reservation is address space, and only
    // what is written becomes memory.
    md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp_scope);
    md_allocator_i* scratch = temp_arena;

    const bool concat_models = (options & MD_PDB_OPTION_CONCAT_MODELS) != 0;

    pdb_block_t b = { .atoms = true };
    md_array(pdb_piece_t)         pieces = 0;
    md_array(md_text_piece_t)     ranges = 0;
    md_array(pdb_residue_start_t) res_starts = 0;   // of the load
    md_array(pdb_residue_start_t) block_res = 0;    // of a block, at the place of the atom of each piece
    md_array(pdb_record_t)        records = 0;      // of a block
    // Of a block
    md_array(md_label_t)          atom_name = 0;
    md_array(uint8_t)             element = 0;
    md_array(float)               column[PDB_NUM_COLUMNS] = {0};
    // Of the load, where they are needed
    md_array(int32_t)             serial = 0;       // Where the file has CONECT records
    md_array(md_label_t)          labels = 0;       // Where atoms have no element: their names, to be typed by
    pdb_column_t                  col[PDB_NUM_COLUMNS] = {0};
    pdb_extras_t extras = { .alloc = temp_arena };
    uint32_t num_bad_ids = 0;
    bool result = false;

    // The first model: [model_beg, model_end) of the atoms read, once its MODEL and ENDMDL records are seen
    size_t num_atoms = 0;
    size_t num_models = 0;
    size_t model_beg = 0;

    // Straight into the state and the system: as many atoms as there can be, as an atom line of coordinates has 54
    // characters at least. The room beyond the atoms is address space, not memory, until written.
    const bool keep_serial = pdb_has_conect(src);
    const int64_t src_size = md_reader_size(src);
    const size_t max_atoms = src_size >= 0 ? (size_t)src_size / 55 + 16 : 0;
    {
        md_array_ensure(state->xyz, max_atoms, state->alloc);
        md_array_ensure(sys->atom.type_idx, ROUND_UP(max_atoms, 16), sys->alloc);
        md_array_ensure(res_starts, max_atoms, scratch);
        if (keep_serial) md_array_ensure(serial, max_atoms, scratch);
        // Of a block: of as many atoms as a block of the reader holds (a string is one block)
        const size_t block_cap = src->cap ? src->cap : (size_t)MAX(0, src_size);
        const size_t block_atoms = block_cap / 55 + 16;
        md_array_ensure(atom_name, block_atoms, scratch);
        md_array_ensure(element,   block_atoms, scratch);
        md_array_ensure(block_res, block_atoms, scratch);
        for (int c = 0; c < PDB_NUM_COLUMNS; ++c) md_array_ensure(column[c], block_atoms, scratch);
        for (int c = 0; c < PDB_NUM_COLUMNS; ++c) col[c].reserve = max_atoms;
        md_array_ensure(pieces, md_text_num_pieces(block_cap, md_executor_num_threads(ex)), scratch);
        md_array_ensure(ranges, md_text_num_pieces(block_cap, md_executor_num_threads(ex)), scratch);
    }

    // The types, as the atoms are read: of the name and the element where the line gives one
	// Add default unknown particle type at index 0
    md_particle_type_find_or_add(&sys->atom.type, STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, sys->alloc);
    md_hashmap32_t type_map = { .allocator = temp_arena };
    // The atom before's type, which the next often shares (a water's hydrogens, an ion after another)
    uint64_t prev_name = 0;
    uint8_t  prev_element = PDB_NO_ELEMENT;
    md_particle_type_idx_t prev_type = 0;

    const uint32_t num_threads = md_executor_num_threads(ex);
    for (;;) {
        const str_t block = md_reader_block(src);
        if (block.len == 0) break;

        const size_t num_pieces = md_text_num_pieces(block.len, num_threads);
        md_array_resize(pieces, num_pieces, scratch);
        md_array_resize(ranges, num_pieces, scratch);
        md_text_split(ranges, num_pieces, block.ptr, block.len);
        for (size_t p = 0; p < num_pieces; ++p) pieces[p] = (pdb_piece_t){ .range = ranges[p] };

        b.text = block.ptr;
        b.len = block.len;
        b.pieces = pieces;
        b.atom_beg = num_atoms;

        if (!md_parallel_for(ex, (uint32_t)num_pieces, 1, pdb_count_lines, &b)) goto done;   // Cancelled
        size_t block_atoms = 0, block_records = 0;
        for (size_t p = 0; p < num_pieces; ++p) {
            pieces[p].first_atom = block_atoms;
            pieces[p].first_record = block_records;
            block_atoms += pieces[p].num_atoms;
            block_records += pieces[p].num_records;
        }

        // Room for the atoms of the block, where they go
        const size_t total = num_atoms + block_atoms;
        PDB_ROOM(state->xyz, total, state->alloc);
        PDB_ROOM(sys->atom.type_idx, total, sys->alloc);
        if (keep_serial) PDB_ROOM(serial, total, scratch);
        PDB_ROOM(atom_name, block_atoms, scratch);
        PDB_ROOM(element,   block_atoms, scratch);
        for (int c = 0; c < PDB_NUM_COLUMNS; ++c) PDB_ROOM(column[c], block_atoms, scratch);
        PDB_ROOM(block_res, block_atoms, scratch);
        PDB_ROOM(records, block_records, scratch);
        b.xyz = state->xyz;
        b.serial = keep_serial ? serial : NULL;
        b.atom_name = atom_name;
        b.element = element;
        for (int c = 0; c < PDB_NUM_COLUMNS; ++c) b.column[c] = column[c];
        b.records = records;
        // Each piece writes its residues at the place of its first atom, put in order below
        b.res_starts = block_res;

        if (!md_parallel_for(ex, (uint32_t)num_pieces, 1, pdb_parse_lines, &b)) goto done;   // Cancelled

        // The residues in order, after those before the block
        bool unassigned = false;
        for (size_t p = 0; p < num_pieces; ++p) {
            const pdb_piece_t* piece = &pieces[p];
            md_array_push_array(res_starts, block_res + piece->first_atom, piece->num_residues, scratch);
            num_bad_ids += piece->num_bad_ids;
            unassigned |= piece->unassigned;
        }
        sys->residue.count = md_array_size(res_starts);   // Counted here, made below

        // The types of the block's atoms. Atoms without an element are typed by their names once all are read: the
        // names are kept from the first such atom on, and those of the atoms before it are those of their types.
        for (size_t k = 0; k < block_atoms; ++k) {
            if (element[k] == PDB_NO_ELEMENT) {
                sys->atom.type_idx[num_atoms + k] = 0;
                continue;
            }
            const uint64_t name = md_label_key(atom_name[k]);
            if (name != prev_name || element[k] != prev_element) {
                prev_type = md_util_particle_type_find_or_add(sys, &type_map, md_label_str(&atom_name[k]), (md_atomic_number_t)element[k]);
                prev_name = name;
                prev_element = element[k];
            }
            sys->atom.type_idx[num_atoms + k] = prev_type;
        }
        if (unassigned && !labels) {
            md_array_ensure(labels, MAX(total, max_atoms), scratch);
            md_array_resize(labels, num_atoms, scratch);
            for (size_t i = 0; i < num_atoms; ++i) labels[i] = sys->atom.type.name[sys->atom.type_idx[i]];
        }
        if (labels) md_array_push_array(labels, atom_name, block_atoms, scratch);
        for (int c = 0; c < PDB_NUM_COLUMNS; ++c) pdb_column_append(&col[c], column[c], num_atoms, block_atoms, scratch);

        // The residue fields of the last atom, and a TER after it, for the next block
        if (block_atoms > 0) {
            for (size_t at = block.len; at > 0;) {
                const str_t line = md_text_line_before(block.ptr, at);
                at = (size_t)(line.ptr - block.ptr);
                if (pdb_line_kind(line) == PDB_LINE_ATOM) {
                    b.has_prev = true;
                    b.prev = pdb_res_fields(line);
                    break;
                }
            }
            b.prev_ter = false;
        }
        num_atoms = total;

        // The other records, in order: the models end the atoms
        for (size_t r = 0; r < block_records; ++r) {
            const pdb_record_t* rec = &records[r];
            const str_t line = { block.ptr + rec->offset, rec->len };
            switch (rec->kind) {
            case PDB_LINE_TER:
                if (rec->atom == num_atoms) b.prev_ter = true;
                break;
            case PDB_LINE_MODEL:
                if (!concat_models && b.atoms) {
                    if (num_models == 0) {
                        model_beg = rec->atom;
                    } else {
                        // A second model, without an ENDMDL of the first
                        num_atoms = rec->atom;
                        b.atoms = false;
                    }
                }
                num_models += 1;
                break;
            case PDB_LINE_ENDMDL:
                if (!concat_models && b.atoms && num_models > 0) {
                    num_atoms = rec->atom;
                    b.atoms = false;
                }
                break;
            case PDB_LINE_CRYST1:
            case PDB_LINE_CONECT:
            case PDB_LINE_REMARK_350:
                pdb_extras_record(&extras, rec->kind, line);
                break;
            default:
                break;
            }
        }
        md_reader_advance(src, block.len);
    }

    // The atoms of the first model: those before its MODEL record and after its end are not
    {
        const size_t n = num_atoms - model_beg;
        if (model_beg > 0) {
            MEMMOVE(state->xyz,         state->xyz + model_beg,         n * sizeof(vec3_t));
            MEMMOVE(sys->atom.type_idx, sys->atom.type_idx + model_beg, n * sizeof(md_particle_type_idx_t));
            if (serial) MEMMOVE(serial, serial + model_beg, n * sizeof(int32_t));
            if (labels) MEMMOVE(labels, labels + model_beg, n * sizeof(md_label_t));
            for (int c = 0; c < PDB_NUM_COLUMNS; ++c) {
                if (col[c].data) MEMMOVE(col[c].data, col[c].data + model_beg, n * sizeof(float));
            }
        }
        // The residue of the first atom of the model, which may have begun before it
        pdb_residue_start_t first = {0};
        for (size_t r = 0; r < sys->residue.count && res_starts[r].atom <= model_beg; ++r) first = res_starts[r];
        size_t num_res = 0;
        for (size_t r = 0; r < sys->residue.count; ++r) {
            pdb_residue_start_t rs = res_starts[r];
            if (rs.atom < model_beg || rs.atom >= num_atoms) continue;
            rs.atom -= (uint32_t)model_beg;
            res_starts[num_res++] = rs;
        }
        if (n > 0 && (num_res == 0 || res_starts[0].atom != 0)) {
            first.atom = 0;
            md_array_push(res_starts, first, scratch);     // Room for one more
            MEMMOVE(res_starts + 1, res_starts, num_res * sizeof(pdb_residue_start_t));
            res_starts[0] = first;
            num_res += 1;
        }
        sys->residue.count = num_res;
        num_atoms = n;
    }

    if (num_bad_ids > 0) {
        MD_LOG_DEBUG("%u atom serials or residue numbers are not numbers of the PDB format ('*****' past the width of the field, "
            "or hexadecimal as VMD writes past 99999): the atoms are taken in their order", num_bad_ids);
    }

    // The atoms
    md_array_resize(state->xyz, num_atoms, state->alloc);
    state->num_atoms = num_atoms;
    sys->atom.count = num_atoms;
    {
        const size_t capacity = ROUND_UP(num_atoms, 16);
        PDB_ROOM(sys->atom.type_idx, capacity, sys->alloc);
        MEMSET(sys->atom.type_idx + num_atoms, 0, (capacity - num_atoms) * sizeof(md_particle_type_idx_t));
        md_array_resize(sys->atom.flags, capacity, sys->alloc);
        MEMSET(sys->atom.flags, 0, capacity * sizeof(md_atom_flags_t));
    }

    // The residues, and their chain ids
    const size_t num_res = sys->residue.count;
    md_array(str_t) res_chain_ids = 0;
    md_array_resize(res_chain_ids, num_res, temp_arena);
    md_array(bool) res_after_ter = 0;
    md_array_resize(res_after_ter, num_res, temp_arena);
    md_array_resize(sys->residue.atom_offset, num_res + 1, sys->alloc);
    md_array_resize(sys->residue.name,   num_res, sys->alloc);
    md_array_resize(sys->residue.seq_id, num_res, sys->alloc);
    md_array_resize(sys->residue.flags,  num_res, sys->alloc);
    char* chain_chars = md_alloc(temp_arena, MAX(num_res, 1));
    for (size_t r = 0; r < num_res; ++r) {
        const pdb_residue_start_t* rs = &res_starts[r];
        md_sequence_id_t seq_id = rs->seq_id;
        if (seq_id == INT_MAX) {
			// Missing residue sequence id, assign one based on previous residue
			seq_id = r > 0 ? sys->residue.seq_id[r - 1] + 1 : 1;
		}
        sys->residue.atom_offset[r] = rs->atom;
        sys->residue.name[r]   = rs->name;
        sys->residue.seq_id[r] = seq_id;
        sys->residue.flags[r]  = MD_RESIDUE_FLAG_NONE;
        chain_chars[r] = rs->chain_id;
        res_chain_ids[r] = str_trim((str_t){ &chain_chars[r], 1 });
        res_after_ter[r] = rs->after_ter;
    }
    sys->residue.atom_offset[num_res] = (uint32_t)num_atoms;    // Final sentinel

    if (labels) {
        md_util_system_infer_particle_types_from_labels(sys, labels);
    }

    // The columns which hold more than one value: md_mmcif publishes the same, and they have to agree
    if (col[PDB_COL_OCCUPANCY].data) {
        md_attributes_publish_atom_column(&sys->attributes, STR_LIT("atom/occupancy"), md_unit_none(), 1, col[PDB_COL_OCCUPANCY].data, num_atoms);
    }
    if (col[PDB_COL_B_FACTOR].data) {
        md_attributes_publish_atom_column(&sys->attributes, STR_LIT("atom/b_factor"), md_unit_pow(md_unit_angstrom(), 2), 1, col[PDB_COL_B_FACTOR].data, num_atoms);
    }
    if (col[PDB_COL_CHARGE].data) {
        // Chemistry perception starts from what the file gives
        md_attributes_publish_atom_column(&sys->attributes, STR_LIT("atom/formal_charge"), md_unit_none(), 1, col[PDB_COL_CHARGE].data, num_atoms);
    }
    pdb_system_finish(sys, state, &extras, serial, res_chain_ids, res_after_ter, temp_arena);
    result = true;

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

bool md_pdb_system_init_from_str(md_system_t* sys, md_system_state_t* state, str_t str, md_pdb_options_t options, const md_executor_t* executor) {
    md_reader_i src = md_reader_from_str(str);
    return pdb_load(sys, state, &src, options, executor);
}

bool md_pdb_system_init_from_file(md_system_t* sys, md_system_state_t* state, str_t filename, md_pdb_options_t options, const md_executor_t* executor) {
    // The reader's buffer in a temp scope, below the one of the load (pdb_load)
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_file_reader_t fr;
    if (!md_file_reader_open(&fr, filename, MD_TEXT_BLOCK_SIZE, md_temp_allocator(temp))) {
        MD_LOG_ERROR("Could not open file '" STR_FMT "'", STR_ARG(filename));
        md_temp_end(temp);
        return false;
    }
    // Several models are a trajectory as well: md_pdb_system_publish_run makes them a run.
    const bool result = pdb_load(sys, state, &fr.r, options, executor);
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
static bool try_read_cache(pdb_cache_t* cache, str_t path, md_allocator_i* alloc) {
    ASSERT(cache);
    ASSERT(alloc);
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
        MD_LOG_ERROR("The PDB cache beside '" STR_FMT "' is incomplete", STR_ARG(path));
        md_free(alloc, cache->frame_offsets, offset_bytes);
        cache->frame_offsets = NULL;
    }
    md_file_close(&file);
    return ok;
}

static bool write_cache(const pdb_cache_t* cache, str_t path, const md_file_info_t* scanned) {
    md_file_t file = {0};
    if (!md_run_cache_create(&file, path, scanned, MD_PDB_CACHE_MAGIC, MD_PDB_CACHE_VERSION, cache->header.num_atoms, cache->header.num_frames)) {
        return false;
    }
    const size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
    const bool ok =
        md_file_write(file, &cache->cell, sizeof(cache->cell)) == sizeof(cache->cell) &&
        md_file_write(file, cache->frame_offsets, offset_bytes) == offset_bytes;
    if (!ok) {
        MD_LOG_ERROR("Failed to write the PDB cache beside '" STR_FMT "'", STR_ARG(path));
    }
    md_file_close(&file);
    return ok;
}

// Where each model is and the file's cell, from the cache beside the file when that is current and from a scan of the
// file otherwise (writing the cache unless told not to). The scan reads the names of the records, and of the lines
// nothing else but the first CRYST1's: a model begins at its MODEL record, and its atoms are those the system has of
// it (the first alternate locations, pdb_line_kind), up to its ENDMDL record or the next MODEL record.
static bool pdb_index_load(pdb_cache_t* cache, str_t filename, md_run_flags_t flags, md_allocator_i* alloc) {
    MEMSET(cache, 0, sizeof(*cache));

    if (try_read_cache(cache, filename, alloc)) {
        return true;
    }
    MEMSET(cache, 0, sizeof(*cache));

    // The file as it is before the scan is what the cache is stamped with
    md_file_info_t scanned = { 0 };
    if (!md_file_info_extract_from_path(filename, &scanned)) {
        MD_LOG_ERROR("Failed to extract file info from path '" STR_FMT "'", STR_ARG(filename));
        return false;
    }

    md_allocator_i* heap = md_get_heap_allocator();
    md_file_reader_t fr;
    if (!md_file_reader_open(&fr, filename, MD_READER_DEFAULT_CAPACITY, heap)) {
        MD_LOG_ERROR("Could not open file '" STR_FMT "'", STR_ARG(filename));
        return false;
    }

    md_array(int64_t)  model_offset = 0;
    md_array(uint32_t) model_atoms = 0;
    pdb_extras_t extras = { .alloc = heap };    // The cell
    bool in_model = false;
    str_t line;
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
    MEMCPY(cache->frame_offsets, model_offset, num_models * sizeof(int64_t));
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
static size_t pdb_parse_model_coords(float* xyz, str_t str, size_t first, size_t count) {
    str_t line;
    size_t i = 0;
    size_t written = 0;
    while (written < count && str_extract_line(&line, &str)) {
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
    ASSERT(sys);
    if (!slice || slice->num_idx == 0 || slice->num_idx > 2) return 0;

    md_run_source_t src;
    if (!md_run_source(&src, &sys->attributes, attr, STR_LIT("atom/position"))) return 0;

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
    char* text = md_temp_alloc(temp, MAX(frame_size, 1));
    if (text && md_attribute_io_read_at(io, src.path, src.offset[frame], text, frame_size) == frame_size) {
        if (pdb_parse_model_coords((float*)dst, (str_t){ text, frame_size }, first, count) == count) {
            written = cap;
        } else {
            MD_LOG_ERROR("PDB: model %u of '" STR_FMT "' has fewer than %zu atoms", frame, STR_ARG(src.path), first + count);
        }
    } else {
        MD_LOG_ERROR("PDB: Failed to read model %u from '" STR_FMT "'", frame, STR_ARG(src.path));
    }
    md_temp_end(temp);
    return written;
}

bool md_pdb_system_publish_run(md_system_t* sys, str_t filename, str_t run, uint32_t flags) {
    ASSERT(sys);
    char path_buf[4096];
    const size_t path_len = md_path_write_canonical(path_buf, sizeof(path_buf), filename);
    const str_t path = { path_buf, path_len };

    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
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
        MEMCPY(boxes + i * 9, A, sizeof(A));
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
        str_copy_to_char_buf(element, sizeof(element), md_atomic_number_symbol(atomic_nr));
        
        // Get atom name (use element symbol if not available)
		char name[5] = "";
		str_copy_to_char_buf(name, sizeof(name), md_atom_name(&sys->atom, idx));
        
        // Get residue name
		md_residue_idx_t res_idx = md_system_residue_find_by_atom_idx(sys, idx);
		char resname[5] = "";
		str_copy_to_char_buf(resname, sizeof(resname), md_system_residue_name(sys, res_idx));

		md_sequence_id_t res_seq = md_system_residue_seq_id(sys, res_idx);
        
        // The chain id, the label the atom is filed under (blank for one in no chain)
        char chain_id[4] = " ";
        const md_chain_idx_t chain_idx = md_system_chain_find_by_atom_idx(sys, idx);
        if (chain_idx >= 0) str_copy_to_char_buf(chain_id, sizeof(chain_id), md_system_chain_id(sys, (size_t)chain_idx));
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
