#include "md_infer_internal.h"
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <md_dssp.h>
#include <core/md_geom.h>

// Polymer inference: the roles of the atoms of amino acids and nucleotides, the classification of residues, the
// protein and nucleic acid backbones and their chain breaks, and the backbone of a state (md_util.h).

bool md_backbone_angles_compute(md_backbone_angles_t backbone_angles[], size_t capacity, const md_float3_t* xyz, const md_unitcell_t* cell, const md_protein_backbone_data_t* backbone) {
    if (!backbone_angles) return false;
    if (capacity == 0) return false;

    md_memset(backbone_angles, 0, sizeof(md_backbone_angles_t) * capacity);

    if (!backbone) return false;
    if (!xyz) return false;
    if (!backbone->segment.atoms) return false;

    for (size_t bb_idx = 0; bb_idx < backbone->range.count; ++bb_idx) {
        const md_urange_t range = {backbone->range.offset[bb_idx], backbone->range.offset[bb_idx + 1]};
        MD_ASSERT(range.end <= capacity);

        // A trace has one bead per residue, and so no phi and psi
        if (range.end - range.beg < 4 || md_protein_backbone_range_is_trace(backbone, bb_idx)) {
            continue;
        }

        for (size_t i = range.beg + 1; i + 1 < range.end; ++i) {
			int c_prev_idx = backbone->segment.atoms[i - 1].c;
			int ca_idx     = backbone->segment.atoms[i].ca;
			int  c_idx     = backbone->segment.atoms[i].c;
			int  n_idx     = backbone->segment.atoms[i].n;
			int n_next_idx = backbone->segment.atoms[i + 1].n;

			md_float3_t c_prev = xyz[c_prev_idx];
            md_float3_t n      = xyz[n_idx];
            md_float3_t ca     = xyz[ca_idx];
            md_float3_t c      = xyz[c_idx];
            md_float3_t n_next = xyz[n_next_idx];

            md_float3_t d[4] = {
                md_float3_sub(n, c_prev),
                md_float3_sub(ca, n),
                md_float3_sub(c, ca),
                md_float3_sub(n_next, c),
            };

            md_geom_min_image((float*)d, 0, MD_COUNTOF(d), cell);

            backbone_angles[i].phi = md_float3_dihedral_angle(d[0], d[1], d[2]);
            backbone_angles[i].psi = md_float3_dihedral_angle(d[1], d[2], d[3]);
        }
    }

    return true;
}

// ### THE BACKBONE OF A STATE ###

MD_STATIC_ASSERT(sizeof(md_secondary_structure_t) == sizeof(int32_t), "The secondary structure is stored as I32");
MD_STATIC_ASSERT(sizeof(md_backbone_angles_t) == 2 * sizeof(float), "The backbone angles are stored as F32 x 2");

static const md_attribute_t* state_backbone_attr(const md_system_state_t* state, md_str_t path, md_attribute_type_t type, uint32_t components, size_t num_segments) {
    if (!state || num_segments == 0) return NULL;
    const md_attribute_t* attr = md_attributes_find(&state->attributes, path);
    if (!md_attribute_view(attr, type, components, 1) || attr->format.shape[0] != num_segments) return NULL;
    return attr;
}

static void* state_backbone_write(md_system_state_t* state, const md_system_t* sys, md_str_t path, md_attribute_type_t type, uint32_t components, md_unit_t unit, md_str_t label) {
    if (!state || !sys) return NULL;
    const size_t num_segments = sys->protein_backbone.segment.count;
    if (num_segments == 0) return NULL;
    if (!state->attributes.alloc) {
        // The table of a state allocates with the state's allocator (see md_system_state_t); a view has neither
        if (!state->alloc) return NULL;
        state->attributes.alloc = state->alloc;
    }

    const md_attribute_t* attr = state_backbone_attr(state, path, type, components, num_segments);
    md_attribute_id_t id = attr ? attr->id : MD_ATTRIBUTE_INVALID;
    if (id == MD_ATTRIBUTE_INVALID) {
        const md_attribute_desc_t desc = {
            .path   = path,
            .format = { .type = type, .components = components, .rank = 1, .shape = { (uint32_t)num_segments } },
            .unit   = unit,
            .label  = label,
        };
        id = md_attributes_replace(&state->attributes, &desc);
        if (id == MD_ATTRIBUTE_INVALID) return NULL;
    }
    md_attributes_touch(&state->attributes, id);
    return md_attributes_data(&state->attributes, id, type);
}

const md_backbone_angles_t* md_system_state_backbone_angles(const md_system_state_t* state, const md_system_t* sys) {
    if (!sys) return NULL;
    const md_attribute_t* attr = state_backbone_attr(state, MD_STR_LIT(MD_BACKBONE_ANGLE_PATH), MD_ATTRIBUTE_TYPE_F32, 2, sys->protein_backbone.segment.count);
    return attr ? (const md_backbone_angles_t*)attr->data : NULL;
}

const md_secondary_structure_t* md_system_state_secondary_structure(const md_system_state_t* state, const md_system_t* sys) {
    if (!sys) return NULL;
    const md_attribute_t* attr = state_backbone_attr(state, MD_STR_LIT(MD_BACKBONE_SECONDARY_STRUCTURE_PATH), MD_ATTRIBUTE_TYPE_I32, 1, sys->protein_backbone.segment.count);
    return attr ? (const md_secondary_structure_t*)attr->data : NULL;
}

md_backbone_angles_t* md_system_state_backbone_angles_write(md_system_state_t* state, const md_system_t* sys) {
    return (md_backbone_angles_t*)state_backbone_write(state, sys, MD_STR_LIT(MD_BACKBONE_ANGLE_PATH), MD_ATTRIBUTE_TYPE_F32, 2, md_unit_radian(), MD_STR_LIT("Backbone Angles"));
}

md_secondary_structure_t* md_system_state_secondary_structure_write(md_system_state_t* state, const md_system_t* sys) {
    return (md_secondary_structure_t*)state_backbone_write(state, sys, MD_STR_LIT(MD_BACKBONE_SECONDARY_STRUCTURE_PATH), MD_ATTRIBUTE_TYPE_I32, 1, md_unit_none(), MD_STR_LIT("Secondary Structure"));
}

bool md_system_state_backbone_compute(md_system_state_t* state, const md_system_t* sys) {
    if (!state || !sys || !md_system_state_has_coords(state)) return false;
    const size_t num_segments = sys->protein_backbone.segment.count;
    if (num_segments == 0) return false;
    md_backbone_angles_t*     angle = md_system_state_backbone_angles_write(state, sys);
    md_secondary_structure_t* ss    = md_system_state_secondary_structure_write(state, sys);
    if (!angle || !ss) return false;
    const bool ok_angle = md_backbone_angles_compute(angle, num_segments, state->xyz, &state->unitcell, &sys->protein_backbone);
    const bool ok_ss    = md_dssp_compute(ss, num_segments, state->xyz, &state->unitcell, &sys->protein_backbone);
    return ok_angle && ok_ss;
}

// Of segment i of an atomistic backbone: 'G' a glycine, 'P' a proline, 0 any other residue. By the name of its
// residue where that is the name of an amino acid, else (a segment of no residue, or of one the data does not name as
// one: "MOL" of a LAMMPS data file) by its atoms: a glycine has no CB, and the N of a
// proline closes the ring N-CA-CB-CG-CD: it is bonded to a CD (besides its CA and the C of the segment before it, if
// it has_prev) which is bonded to a CG which is bonded to the CB. The carbonyl of a cap (ACE) on an N terminus closes
// none.
static char ramachandran_residue(const md_system_t* sys, const md_protein_backbone_data_t* backbone, size_t i, bool has_prev) {
    const md_residue_idx_t res_idx = backbone->segment.res_idx[i];
    if (res_idx >= 0 && md_resname_amino_acid(md_system_residue_name(sys, res_idx))) {
        const md_str_t name = md_system_residue_name(sys, res_idx);
        return md_str_eq(name, MD_STR_LIT("GLY")) ? 'G' : md_str_eq(name, MD_STR_LIT("PRO")) ? 'P' : 0;
    }
    const md_amino_acid_atoms_t* atoms = &backbone->segment.atoms[i];
    const md_atom_idx_t prev_c = has_prev ? backbone->segment.atoms[i - 1].c : -1;
    if (atoms->cb < 0) return 'G';
    md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)atoms->n);
    for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
        const md_atom_idx_t cd = md_bond_iter_atom_index(&it);
        if (cd == atoms->ca || cd == prev_c || md_atom_atomic_number(&sys->atom, (size_t)cd) != MD_Z_C) continue;
        md_bond_iter_t it_cd = md_bond_iter(&sys->bond, (size_t)cd);
        for (; md_bond_iter_has_next(&it_cd); md_bond_iter_next(&it_cd)) {
            const md_atom_idx_t cg = md_bond_iter_atom_index(&it_cd);
            if (cg != atoms->n && md_bond_find(&sys->bond, cg, atoms->cb) >= 0) return 'P';
        }
    }
    return 0;
}

bool md_backbone_ramachandran_classify(md_ramachandran_type_t ramachandran_types[], size_t capacity, const md_system_t* sys) {
    MD_ASSERT(ramachandran_types);
    md_memset(ramachandran_types, MD_RAMACHANDRAN_TYPE_UNKNOWN, sizeof(md_ramachandran_type_t) * capacity);

    if (capacity == 0) return false;

    const md_protein_backbone_data_t* backbone = &sys->protein_backbone;
    if (backbone->segment.count == 0) return false;
    if (backbone->range.count == 0 || !backbone->range.offset) return false;

    MD_ASSERT(backbone->segment.res_idx);

    // The types are MolProbity's (mmtbx.validation.ramalyze), decided in its order: glycine, proline, then pre-proline
    // for any other residue followed by a proline. A glycine or a proline before a proline keeps its own type.
    //
    // Walked range by range, because a residue is pre-proline only when the proline follows it in the same chain. The
    // last segment of a range has no successor; segment i + 1 is then the first residue of another chain.
    for (size_t range_idx = 0; range_idx < backbone->range.count; ++range_idx) {
        // A trace has no phi and psi to plot: its segments stay UNKNOWN
        if (md_protein_backbone_range_is_trace(backbone, range_idx)) continue;
        const size_t range_beg = backbone->range.offset[range_idx];
        const size_t range_end = backbone->range.offset[range_idx + 1];
        // Written up to the capacity, but the successor is read from the system, so the last written type is still
        // decided by the segment after it.
        const size_t end = MD_MIN(range_end, capacity);

        for (size_t i = range_beg; i < end; ++i) {
            const char res = ramachandran_residue(sys, backbone, i, i > range_beg);
            if (res == 'G') {
                ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_GLYCINE;
            } else if (res == 'P') {
                ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_PROLINE;
            } else if (i + 1 < range_end && ramachandran_residue(sys, backbone, i + 1, true) == 'P') {
                ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_PREPROL;
            } else {
                ramachandran_types[i] = MD_RAMACHANDRAN_TYPE_GENERAL;
            }
        }
    }

    return true;
}

#define MIN_POLYMER_RES_LEN 2

// @TODO: Convert to a table
static bool monatomic_ion_element(md_atomic_number_t elem) {
    switch (elem) {
    case 1:   // H+ / H-
    case 3:   // Li+
    case 9:   // F-
    case 11:  // Na+
    case 17:  // Cl-
    case 19:  // K+
    case 35:  // Br-
    case 37:  // Rb+
    case 47:  // Ag+
    case 53:  // I-
    case 55:  // Cs+
    case 4:   // Be2+
    case 8:   // O2-
    case 12:  // Mg2+
    case 20:  // Ca2+
    case 16:  // S2-
    case 30:  // Zn2+
    case 34:  // Se2-
    case 38:  // Sr2+
    case 48:  // Cd2+
    case 56:  // Ba2+
    case 66:  // Hg2+
    case 7:   // N3-
    case 13:  // Al3+
    case 15:  // P3-
    case 29:  // Cu1+ / Cu2+
    case 24:  // Cr2+ / Cr3+
    case 25:  // Mn2+ / Mn3+
    case 26:  // Fe2+ / Fe3+
    case 27:  // Co2+ / Co3+
    case 50:  // Sn2+ / Sn4+
    case 68:  // Pb2+ / Pb4+
        return true;
    default:
        return false;
    }
}

// helper function to get atom indices within a 'tree branch' of the connectivity graph of a system.
static inline bool get_branch_atom_indices(md_array(int)* out_indices, const md_bond_data_t* bond, int start_atom, int exclude_atom, md_allocator_i* alloc) {
    md_array_push(*out_indices, start_atom, alloc);
    md_bond_iter_t iter = md_bond_iter(bond, start_atom);
    while (md_bond_iter_has_next(&iter)) {
        int neighbor = md_bond_iter_atom_index(&iter);
        if (neighbor != exclude_atom) {
            if (!get_branch_atom_indices(out_indices, bond, neighbor, start_atom, alloc)) {
                return false;
            }
        }
        md_bond_iter_next(&iter);
    }
    return true;
}

// ### BACKBONES FROM THE BONDS ###
//
// The backbones of an atomistic system are found from its covalent bonds alone (backbone_paths_ncco, _sugar),
// whatever its atoms and residues are called: an XYZ file or the geometry of a QM calculation has no residues, a LAMMPS
// data file names its atoms by their types and its residues "MOL". The names play no part. What is found is what a
// cartoon and the secondary structure need: the atoms of each segment, in order along each chain.
//
//   AMINO ACIDS  N-CA-C-O, the bond orders not tested: they are md_chem's perception of bonds which were themselves
//                inferred from distances, and one wrong bond (a clash between two chains taken for a bond) moves the
//                double bond of a carbonyl onto its C-N and loses the residue. Without orders N-CA-CB-OG of a serine or a
//                threonine matches too, as does N-CD-CG-OD1 of a hydroxyproline. The matches are taken best first, each
//                where its N, CA and C are not taken yet: one whose C is bonded to the N of another match (a peptide
//                bond), then one whose C has two oxygens (the carboxylate of a C terminus), then one whose C-O is as short
//                as a double bond (which tells the carbonyl of a C terminus without its OXT from a hydroxyl). CB is the
//                heavy atom on CA besides N and C (none for a glycine), HN the hydrogen of an N which has one.
//   NUCLEOTIDES  The sugar with its backbone, O5'-C5'-C4'(-O4'-C1'-C2'-C3'-)C3'-O3', and the P bonded to the O5'. A
//                sugar matches in one way only.
//   CHAINS       A residue follows the one whose C (O3') is bonded to its N (P). A chain runs from a residue which follows
//                none; what is left are rings (a cyclic peptide), cut before their first residue. A chain shorter than
//                min_length is no backbone: that leaves out the free amino acids and nucleotides (ATP, GTP) and the
//                groups which match without being residues (N7-C5-C6=O6 of a guanine, the ethanolamine of a lipid).
//   RESIDUES     A segment is of the residue which holds its CA (C4'), where that residue holds the CA of no other
//                segment; else of none (-1): a system without residues, or a residue which is a whole chain (a LAMMPS
//                molecule). A range is within one part, the part of its residues, and is cut where that changes; one of
//                residues of none is of no part (-1). Two consecutive residues of a part whose segments are not linked
//                are a CHAIN BREAK (md_util.h), unless the chain ends there.
//
// The ranges come in the order of their first atom (CA or C4'), with the traces of the coarse grained residues.

// Ångström below which a bond from C to O is double: C=O is 1.23 long, C-O 1.43
#define CARBONYL_MAX_DIST 1.33f

static inline bool atoms_closer_than(const md_system_state_t* state, md_atom_idx_t a, md_atom_idx_t b, float max_dist) {
    if (!md_system_state_has_coords(state) || (size_t)a >= state->num_atoms || (size_t)b >= state->num_atoms) return false;
    md_float3_t d = md_float3_sub(md_state_coord(state, (size_t)b), md_state_coord(state, (size_t)a));
    md_geom_min_image(d.elem, 0, 1, &state->unitcell);
    return md_float3_dot(d, d) < max_dist * max_dist;
}

// The first neighbour of atom a of element z, which is neither x nor y (z 0: any element but hydrogen)
static md_atom_idx_t bonded_atom_other(const md_system_t* sys, md_atom_idx_t a, md_atomic_number_t z, md_atom_idx_t x, md_atom_idx_t y) {
    md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)a);
    for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
        const md_atom_idx_t b = md_bond_iter_atom_index(&it);
        if (b == x || b == y) continue;
        const md_atomic_number_t zb = md_atom_atomic_number(&sys->atom, (size_t)b);
        if (z ? zb == z : zb != MD_Z_H) return b;
    }
    return -1;
}

// The paths of the backbone atoms in the covalent bond graph (the bonds md_bond_iter walks: no elastic or metallic
// bonds), by the elements of the atoms alone: bond orders, hydrogens and charges are not tested. A coordination bond
// takes no part, as it is between a metal and an atom of none of these paths. atom_idx holds count paths of width
// atoms each, in the order of their first atom.
typedef struct backbone_paths_t {
    int32_t* atom_idx;
    size_t   count;
} backbone_paths_t;

#define FOR_BONDED(x, a, z)                                                                                      \
    for (md_bond_iter_t it_##x = md_bond_iter(&sys->bond, (size_t)(a)); md_bond_iter_has_next(&it_##x); md_bond_iter_next(&it_##x)) \
        for (md_atom_idx_t x = md_bond_iter_atom_index(&it_##x); x >= 0 && md_atom_atomic_number(&sys->atom, (size_t)x) == (z); x = -1)

// N-CA-C-O: an N, a C, a C and an O bonded in a row
static bool backbone_paths_ncco(backbone_paths_t* out, const md_system_t* sys, md_allocator_i* alloc) {
    md_array(int32_t) idx = 0;
    for (size_t i = 0; i < sys->atom.count; ++i) {
        const md_atom_idx_t n = (md_atom_idx_t)i;
        if (md_atom_atomic_number(&sys->atom, i) != MD_Z_N) continue;
        FOR_BONDED(ca, n, MD_Z_C) {
            FOR_BONDED(c, ca, MD_Z_C) {
                FOR_BONDED(o, c, MD_Z_O) {
                    const int32_t path[4] = { n, ca, c, o };
                    md_array_push_array(idx, path, 4, alloc);
                }
            }
        }
    }
    out->atom_idx = idx;
    out->count = md_array_size(idx) / 4;
    return out->count > 0;
}

// The sugar of a nucleotide with its backbone, O5'-C5'-C4'(-O4'-C1'-C2'-C3'-)C3'-O3': the ring C4' O4' C1' C2' C3' of
// an oxygen and four carbons, C5' a carbon on C4' with the O5' on it, O3' an oxygen on C3'. Eight atoms, all different.
static bool backbone_paths_sugar(backbone_paths_t* out, const md_system_t* sys, md_allocator_i* alloc) {
    md_array(int32_t) idx = 0;
    for (size_t i = 0; i < sys->atom.count; ++i) {
        const md_atom_idx_t c4 = (md_atom_idx_t)i;
        if (md_atom_atomic_number(&sys->atom, i) != MD_Z_C) continue;
        FOR_BONDED(c5, c4, MD_Z_C) {
            FOR_BONDED(o5, c5, MD_Z_O) {
                FOR_BONDED(o4, c4, MD_Z_O) {
                    FOR_BONDED(c1, o4, MD_Z_C) {
                        if (c1 == c4 || c1 == c5) continue;
                        FOR_BONDED(c2, c1, MD_Z_C) {
                            if (c2 == c4 || c2 == c5) continue;
                            FOR_BONDED(c3, c2, MD_Z_C) {
                                if (c3 == c1 || c3 == c4 || c3 == c5 || md_bond_find(&sys->bond, c3, c4) < 0) continue;
                                FOR_BONDED(o3, c3, MD_Z_O) {
                                    if (o3 == o4 || o3 == o5) continue;
                                    const int32_t path[8] = { o5, c5, c4, o4, c1, c2, c3, o3 };
                                    md_array_push_array(idx, path, 8, alloc);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    out->atom_idx = idx;
    out->count = md_array_size(idx) / 8;
    return out->count > 0;
}

#undef FOR_BONDED

static int cmp_u64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a;
    const uint64_t y = *(const uint64_t*)b;
    return (x > y) - (x < y);
}

// The order of the chains of num residues, from the residue each follows (prev, -1 for none) and is followed by (next):
// out_order lists the residues chain by chain, a chain beginning at each residue of out_begin (num + 1 entries each).
// Returns the number of chains.
static size_t backbone_chains(uint32_t* out_order, uint32_t* out_begin, const int32_t* prev, const int32_t* next, size_t num, md_temp_scope_t temp) {
    uint8_t* done = md_temp_alloc_zero_array(temp, uint8_t, num);
    size_t num_order = 0;
    size_t num_chains = 0;
    // Open chains from their first residue, then the rings that are left
    for (int pass = 0; pass < 2; ++pass) {
        for (size_t u = 0; u < num; ++u) {
            if (done[u] || (pass == 0 && prev[u] >= 0)) continue;
            out_begin[num_chains++] = (uint32_t)num_order;
            for (int32_t v = (int32_t)u; v >= 0 && !done[v]; v = next[v]) {
                done[v] = 1;
                out_order[num_order++] = (uint32_t)v;
            }
        }
    }
    out_begin[num_chains] = (uint32_t)num_order;
    return num_chains;
}

// The residues of the polymer chains found from the bonds: the backbone atoms of each (aa of a protein, na of a nucleic
// acid), the residue each follows (prev) and is followed by (next), -1 for none, and the chains they form, every one
// whatever its length: chain c is the residues order[begin[c] .. begin[c + 1]).
typedef struct backbone_chains_t {
    size_t num;
    md_amino_acid_atoms_t*   aa;
    md_nucleic_acid_atoms_t* na;
    int32_t*  prev;                 // [num]
    int32_t*  next;                 // [num]
    uint32_t* order;                // [num]
    uint32_t* begin;                // [num_chains + 1]
    size_t num_chains;
} backbone_chains_t;

static bool protein_chains_from_bonds(backbone_chains_t* out, const md_system_t* sys, const md_system_state_t* state, md_temp_scope_t temp) {
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    backbone_paths_t res = {0};
    // Columns: N CA C O
    if (!backbone_paths_ncco(&res, sys, temp_alloc)) return false;

    const size_t num_atoms = sys->atom.count;
    const size_t num_matches = res.count;

    // The score of each match, its key: score descending, then the order found
    uint8_t* is_n = md_temp_alloc_zero_array(temp, uint8_t, num_atoms);
    for (size_t i = 0; i < num_matches; ++i) {
        is_n[res.atom_idx[i * 4 + 0]] = 1;
    }
    uint64_t* key = md_temp_alloc_array(temp, uint64_t, num_matches);
    for (size_t i = 0; i < num_matches; ++i) {
        const int32_t* m = res.atom_idx + i * 4;
        uint32_t score = 0;
        uint32_t num_o = 0;
        md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)m[2]);
        for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
            const md_atom_idx_t b = md_bond_iter_atom_index(&it);
            const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, (size_t)b);
            if (z == MD_Z_N && b != m[0] && is_n[b]) score |= 4;
            if (z == MD_Z_O) num_o += 1;
        }
        if (num_o >= 2) score |= 2;
        if (atoms_closer_than(state, m[2], m[3], CARBONYL_MAX_DIST)) score |= 1;
        key[i] = (uint64_t)(7 - score) << 32 | i;
    }
    qsort(key, num_matches, sizeof(uint64_t), cmp_u64);

    // Taken best first where N, CA and C are free; the residue whose N, CA or C an atom is, -1 for none
    int32_t* residue_of = md_temp_alloc_array(temp, int32_t, num_atoms);
    md_memset(residue_of, -1, sizeof(int32_t) * num_atoms);
    uint64_t* by_ca = md_temp_alloc_array(temp, uint64_t, num_matches);
    size_t num = 0;
    for (size_t k = 0; k < num_matches; ++k) {
        const uint32_t i = (uint32_t)key[k];
        const int32_t* m = res.atom_idx + i * 4;
        if (residue_of[m[0]] >= 0 || residue_of[m[1]] >= 0 || residue_of[m[2]] >= 0) continue;
        residue_of[m[0]] = residue_of[m[1]] = residue_of[m[2]] = 0;
        by_ca[num++] = (uint64_t)m[1] << 32 | i;
    }
    qsort(by_ca, num, sizeof(uint64_t), cmp_u64);

    md_amino_acid_atoms_t* res_atoms = md_temp_alloc_array(temp, md_amino_acid_atoms_t, num);
    for (size_t u = 0; u < num; ++u) {
        const int32_t* m = res.atom_idx + (uint32_t)by_ca[u] * 4;
        md_amino_acid_atoms_t* a = &res_atoms[u];
        a->n  = m[0];
        a->ca = m[1];
        a->c  = m[2];
        a->o  = m[3];
        a->cb = bonded_atom_other(sys, a->ca, 0, a->n, a->c);
        a->hn = bonded_atom_other(sys, a->n, MD_Z_H, -1, -1);
        if (a->hn >= 0 && bonded_atom_other(sys, a->n, MD_Z_H, a->hn, -1) >= 0) a->hn = -1;
        residue_of[a->n] = residue_of[a->ca] = residue_of[a->c] = (int32_t)u;
    }

    // Residue v follows u when the C of u is bonded to the N of v
    int32_t* prev = md_temp_alloc_array(temp, int32_t, num);
    int32_t* next = md_temp_alloc_array(temp, int32_t, num);
    md_memset(prev, -1, sizeof(int32_t) * num);
    md_memset(next, -1, sizeof(int32_t) * num);
    for (size_t u = 0; u < num; ++u) {
        md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)res_atoms[u].c);
        for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
            const md_atom_idx_t b = md_bond_iter_atom_index(&it);
            const int32_t v = residue_of[b];
            if (v >= 0 && v != (int32_t)u && res_atoms[v].n == b && prev[v] < 0) {
                next[u] = v;
                prev[v] = (int32_t)u;
                break;
            }
        }
    }

    uint32_t* order = md_temp_alloc_array(temp, uint32_t, num);
    uint32_t* begin = md_temp_alloc_array(temp, uint32_t, num + 1);
    const size_t num_chains = backbone_chains(order, begin, prev, next, num, temp);

    out->num = num;
    out->aa = res_atoms;
    out->prev = prev;
    out->next = next;
    out->order = order;
    out->begin = begin;
    out->num_chains = num_chains;
    return num > 0;
}

static bool nucleic_chains_from_bonds(backbone_chains_t* out, const md_system_t* sys, md_temp_scope_t temp) {
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    backbone_paths_t res = {0};
    // Columns: O5' C5' C4' O4' C1' C2' C3' O3'
    if (!backbone_paths_sugar(&res, sys, temp_alloc)) return false;

    const size_t num_atoms = sys->atom.count;

    // A sugar matches in one way only; the sugars by their C4'
    uint64_t* by_c4 = md_temp_alloc_array(temp, uint64_t, res.count);
    for (size_t i = 0; i < res.count; ++i) {
        by_c4[i] = (uint64_t)res.atom_idx[i * 8 + 2] << 32 | i;
    }
    qsort(by_c4, res.count, sizeof(uint64_t), cmp_u64);

    // The residue whose P an atom is, -1 for none
    int32_t* residue_of_p = md_temp_alloc_array(temp, int32_t, num_atoms);
    md_memset(residue_of_p, -1, sizeof(int32_t) * num_atoms);
    md_nucleic_acid_atoms_t* res_atoms = md_temp_alloc_array(temp, md_nucleic_acid_atoms_t, res.count);
    size_t num = 0;
    for (size_t k = 0; k < res.count; ++k) {
        if (k > 0 && by_c4[k] >> 32 == by_c4[k - 1] >> 32) continue;
        const int32_t* m = res.atom_idx + (uint32_t)by_c4[k] * 8;
        md_nucleic_acid_atoms_t* a = &res_atoms[num];
        a->o5 = m[0];
        a->c5 = m[1];
        a->c4 = m[2];
        a->o4 = m[3];
        a->c1 = m[4];
        a->c2 = m[5];
        a->c3 = m[6];
        a->o3 = m[7];
        a->p  = bonded_atom_other(sys, a->o5, MD_Z_P, -1, -1);
        if (a->p >= 0 && residue_of_p[a->p] < 0) residue_of_p[a->p] = (int32_t)num;
        num += 1;
    }

    // Residue v follows u when the O3' of u is bonded to the P of v
    int32_t* prev = md_temp_alloc_array(temp, int32_t, num);
    int32_t* next = md_temp_alloc_array(temp, int32_t, num);
    md_memset(prev, -1, sizeof(int32_t) * num);
    md_memset(next, -1, sizeof(int32_t) * num);
    for (size_t u = 0; u < num; ++u) {
        md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)res_atoms[u].o3);
        for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
            const int32_t v = residue_of_p[md_bond_iter_atom_index(&it)];
            if (v >= 0 && v != (int32_t)u && prev[v] < 0) {
                next[u] = v;
                prev[v] = (int32_t)u;
                break;
            }
        }
    }

    uint32_t* order = md_temp_alloc_array(temp, uint32_t, num);
    uint32_t* begin = md_temp_alloc_array(temp, uint32_t, num + 1);
    const size_t num_chains = backbone_chains(order, begin, prev, next, num, temp);

    out->num = num;
    out->na = res_atoms;
    out->prev = prev;
    out->next = next;
    out->order = order;
    out->begin = begin;
    out->num_chains = num_chains;
    return num > 0;
}

// Gives the flag to start_atom and to what its bonds reach from it within the residue (range): not past exclude_atom,
// nor into an atom with any of stop_flags (the backbone, for a side chain). A bond leaving the residue is not followed:
// the peptide bond (which a proline side chain reaches through its N), a disulfide, a metal coordinated to it. Nor is a
// coordination bond, where the range is not a residue's (POLYMER ATOMS).
static inline void topo_floodfill_flag(md_system_t* sys, int start_atom, int exclude_atom, md_atom_flags_t flag, md_urange_t range, md_atom_flags_t stop_flags, bool covalent_only) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);
    md_array(int) stack = 0;
    md_array_push(stack, start_atom, alloc);
    while (md_array_size(stack) > 0) {
        int atom_idx = stack[md_array_size(stack) - 1];
        md_array_pop(stack);
        if ((sys->atom.flags[atom_idx] & flag) == 0) {
            sys->atom.flags[atom_idx] |= flag;
            md_bond_iter_t iter = md_bond_iter(&sys->bond, atom_idx);
            while (md_bond_iter_has_next(&iter)) {
                int nbr = md_bond_iter_atom_index(&iter);
                const bool follow = !covalent_only || md_bond_kind((md_bond_flags_t)md_bond_iter_bond_flags(&iter)) == MD_BOND_KIND_COVALENT;
                if (follow && nbr != exclude_atom && range.beg <= (uint32_t)nbr && (uint32_t)nbr < range.end &&
                    (sys->atom.flags[nbr] & flag) == 0 && (sys->atom.flags[nbr] & stop_flags) == 0) {
                    md_array_push(stack, nbr, alloc);
                }
                md_bond_iter_next(&iter);
            }
        }
    }
    md_temp_end(temp);
}

// Whether the atoms of the range are all coarse grained beads (or their virtual sites)
static bool atoms_are_beads(const md_system_t* sys, md_urange_t range) {
    if (range.beg >= range.end) return false;
    for (uint32_t i = range.beg; i < range.end; ++i) {
        if (md_atom_particle_kind(&sys->atom, i) == MD_PARTICLE_ATOM || md_atom_atomic_number(&sys->atom, i) != 0) return false;
    }
    return true;
}

// Whether the atom has a covalent bond (not a coordination bond to a metal)
static bool atom_bonded_covalently(const md_system_t* sys, md_atom_idx_t atom) {
    md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)atom);
    for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
        if (md_bond_kind((md_bond_flags_t)md_bond_iter_bond_flags(&it)) == MD_BOND_KIND_COVALENT) return true;
    }
    return false;
}

// True when every neighbor of the atom other than 'except' is of element allowed_z: the atom ends its chain
static bool atom_is_chain_end(const md_system_t* sys, int atom, int except, md_atomic_number_t allowed_z) {
    md_bond_iter_t it = md_bond_iter(&sys->bond, atom);
    while (md_bond_iter_has_next(&it)) {
        const int nbr = md_bond_iter_atom_index(&it);
        if (nbr != except && md_atom_atomic_number(&sys->atom, nbr) != allowed_z) {
            return false;
        }
        md_bond_iter_next(&it);
    }
    return true;
}

// Whether the atom is bonded to an atom of element z outside of the range (another residue)
static bool atom_bonded_outside(const md_system_t* sys, int atom, md_urange_t range, md_atomic_number_t z) {
    md_bond_iter_t it = md_bond_iter(&sys->bond, atom);
    while (md_bond_iter_has_next(&it)) {
        const int nbr = md_bond_iter_atom_index(&it);
        if ((nbr < (int)range.beg || nbr >= (int)range.end) && md_atom_atomic_number(&sys->atom, nbr) == z) {
            return true;
        }
        md_bond_iter_next(&it);
    }
    return false;
}

// ### POLYMER ATOMS ###
//
// The amino acids and nucleotides of an atomistic system are found once, from its covalent bonds (BACKBONES FROM THE
// BONDS), and their atoms given their roles:
//   MD_ATOM_FLAG_BACKBONE      the N, CA and C of an amino acid; the P, O5', C5', C4', C3' and O3' of a nucleotide
//   MD_ATOM_FLAG_SIDE_CHAIN    what hangs off the CA, from the CB on
//   MD_ATOM_FLAG_NUCLEOSIDE    the sugar (C1', C2', C3', C4', O4') and the base
//   MD_ATOM_FLAG_NUCLEOBASE    the base
//   MD_ATOM_FLAG_TERMINAL_BEG  a free amino group, a free 5' hydroxyl
//   MD_ATOM_FLAG_TERMINAL_END  a free carboxyl group, a free 3' hydroxyl
// Not every match is a residue of a chain. One which is all of a residue of its own (a residue which holds no other)
// is one where the system names it a standard residue, or where it is linked to another residue by its backbone (a
// peptide bond from its N or C, a phosphodiester bond from its P or O3'): that keeps out the ligands which hold the
// atoms of one (GTP, ATP, NAD; an amide next to an alpha carbon) while modified residues in a chain (MSE, PSU) and caps
// (ACE, NME) come in. Any other (a system without residues, a residue which is a whole chain as a LAMMPS molecule is)
// is one where it is linked to another match: the N7-C5-C6=O6 of a guanine is not. A side chain, a base or a terminal
// group reaches as far as the residue where the match is all of one, else as far as the covalent bonds go up to a
// backbone: without residues the glycan on an asparagine is of its side chain.
//
// The roles are what everything after reads the atoms of a residue of a chain by, rather than their names: the kinds of
// the residues (md_system_infer_residue_flags), the ends of a chain (mdi_polymer_chain_ends_between), the links over a
// gap (chain_gap_link), the backbones (md_backbone_infer). Beads keep the roles their model gave them.

static inline bool atom_has_flags(const md_system_t* sys, md_atom_idx_t a, md_atom_flags_t flags) {
    return ((uint32_t)md_atom_flags(&sys->atom, (size_t)a) & (uint32_t)flags) == (uint32_t)flags;
}

// The first neighbour of atom a of element z (0: any element but hydrogen) with all of the flags with and none of
// without, which is neither x nor y; -1 for none
md_atom_idx_t mdi_bonded_atom_with(const md_system_t* sys, md_atom_idx_t a, md_atomic_number_t z, md_atom_flags_t with, md_atom_flags_t without, md_atom_idx_t x, md_atom_idx_t y) {
    if (a < 0) return -1;
    md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)a);
    for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
        const md_atom_idx_t b = md_bond_iter_atom_index(&it);
        if (b == x || b == y) continue;
        const md_atomic_number_t zb = md_atom_atomic_number(&sys->atom, (size_t)b);
        if (z ? zb != z : zb == MD_Z_H) continue;
        const uint32_t f = (uint32_t)md_atom_flags(&sys->atom, (size_t)b);
        if ((f & (uint32_t)with) == (uint32_t)with && (f & (uint32_t)without) == 0) return b;
    }
    return -1;
}

// The oxygen of the carbonyl of C: the one without a hydrogen where there are two (a carboxylic acid), else the first
static md_atom_idx_t carbonyl_oxygen(const md_system_t* sys, md_atom_idx_t c) {
    md_atom_idx_t first = -1;
    md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)c);
    for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
        const md_atom_idx_t o = md_bond_iter_atom_index(&it);
        if (md_atom_atomic_number(&sys->atom, (size_t)o) != MD_Z_O) continue;
        if (mdi_bonded_atom_with(sys, o, MD_Z_H, 0, 0, -1, -1) < 0) return o;
        if (first < 0) first = o;
    }
    return first;
}

// The atoms of the amino acid whose alpha carbon is ca, by their roles: CA is the backbone carbon bonded to a backbone
// N and to the backbone C of a carbonyl (which tells it from the C, whose backbone carbon is the CA). CB is the heavy
// atom on CA besides N and C (none for a glycine), HN the hydrogen of an N which has one. False where ca is no CA.
static bool amino_acid_atoms_of_ca(md_amino_acid_atoms_t* out, const md_system_t* sys, md_atom_idx_t ca) {
    if (md_atom_atomic_number(&sys->atom, (size_t)ca) != MD_Z_C || !atom_has_flags(sys, ca, MD_ATOM_FLAG_BACKBONE)) return false;
    const md_atom_idx_t n = mdi_bonded_atom_with(sys, ca, MD_Z_N, MD_ATOM_FLAG_BACKBONE, 0, -1, -1);
    if (n < 0) return false;
    md_atom_idx_t c = -1;
    md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)ca);
    for (; md_bond_iter_has_next(&it) && c < 0; md_bond_iter_next(&it)) {
        const md_atom_idx_t b = md_bond_iter_atom_index(&it);
        if (md_atom_atomic_number(&sys->atom, (size_t)b) == MD_Z_C && atom_has_flags(sys, b, MD_ATOM_FLAG_BACKBONE) &&
            mdi_bonded_atom_with(sys, b, MD_Z_O, 0, 0, -1, -1) >= 0) {
            c = b;
        }
    }
    if (c < 0) return false;
    md_amino_acid_atoms_t a;
    a.n  = n;
    a.ca = ca;
    a.c  = c;
    a.o  = carbonyl_oxygen(sys, c);
    a.cb = mdi_bonded_atom_with(sys, ca, 0, 0, 0, n, c);
    a.hn = mdi_bonded_atom_with(sys, n, MD_Z_H, 0, 0, -1, -1);
    if (a.hn >= 0 && mdi_bonded_atom_with(sys, n, MD_Z_H, 0, 0, a.hn, -1) >= 0) a.hn = -1;
    if (out) *out = a;
    return true;
}

// The atoms of the nucleotide whose C4' is c4, by their roles: C4' is the carbon of the backbone and of the sugar bonded
// to a carbon of the backbone alone (C5'). P is -1 where the O5' has none. False where c4 is no C4'.
static bool nucleotide_atoms_of_c4(md_nucleic_acid_atoms_t* out, const md_system_t* sys, md_atom_idx_t c4) {
    const md_atom_flags_t BB = MD_ATOM_FLAG_BACKBONE, NS = MD_ATOM_FLAG_NUCLEOSIDE, NB = MD_ATOM_FLAG_NUCLEOBASE;
    if (md_atom_atomic_number(&sys->atom, (size_t)c4) != MD_Z_C || !atom_has_flags(sys, c4, BB | NS)) return false;
    md_nucleic_acid_atoms_t a;
    a.c4 = c4;
    a.c5 = mdi_bonded_atom_with(sys, c4, MD_Z_C, BB, NS, -1, -1);
    a.c3 = mdi_bonded_atom_with(sys, c4, MD_Z_C, BB | NS, 0, -1, -1);
    a.o4 = mdi_bonded_atom_with(sys, c4, MD_Z_O, NS, BB, -1, -1);
    a.o5 = mdi_bonded_atom_with(sys, a.c5, MD_Z_O, BB, 0, -1, -1);
    a.o3 = mdi_bonded_atom_with(sys, a.c3, MD_Z_O, BB, 0, -1, -1);
    a.c1 = mdi_bonded_atom_with(sys, a.o4, MD_Z_C, NS, BB, c4, -1);
    a.c2 = mdi_bonded_atom_with(sys, a.c3, MD_Z_C, NS, BB | NB, c4, -1);
    a.p  = mdi_bonded_atom_with(sys, a.o5, MD_Z_P, BB, 0, -1, -1);
    if (a.c5 < 0 || a.c3 < 0 || a.o4 < 0 || a.o5 < 0 || a.o3 < 0 || a.c1 < 0 || a.c2 < 0) return false;
    if (out) *out = a;
    return true;
}

static inline bool atom_in_range(md_atom_idx_t a, md_urange_t range) {
    return a >= 0 && range.beg <= (uint32_t)a && (uint32_t)a < range.end;
}

// The atoms of the one amino acid of a chain the residue (its atoms range) is: false where it holds the backbone of
// none, or of more than one (a residue which is a whole chain)
bool mdi_residue_amino_acid_atoms(md_amino_acid_atoms_t* out, const md_system_t* sys, md_urange_t range) {
    if (out) md_memset(out, -1, sizeof(*out));
    size_t count = 0;
    for (uint32_t i = range.beg; i < range.end; ++i) {
        md_amino_acid_atoms_t a;
        if (amino_acid_atoms_of_ca(&a, sys, (md_atom_idx_t)i) && atom_in_range(a.n, range) && atom_in_range(a.c, range)) {
            if (count++ == 0 && out) *out = a;
        }
    }
    return count == 1;
}

// As mdi_residue_amino_acid_atoms, of a nucleotide
bool mdi_residue_nucleotide_atoms(md_nucleic_acid_atoms_t* out, const md_system_t* sys, md_urange_t range) {
    if (out) md_memset(out, -1, sizeof(*out));
    size_t count = 0;
    for (uint32_t i = range.beg; i < range.end; ++i) {
        md_nucleic_acid_atoms_t a;
        if (nucleotide_atoms_of_c4(&a, sys, (md_atom_idx_t)i) && atom_in_range(a.o5, range) && atom_in_range(a.o3, range)) {
            if (count++ == 0 && out) *out = a;
        }
    }
    return count == 1;
}

// The residue each match of the chains is all of, where it is the only match the residue holds; -1 for none
static int32_t* chain_units_own_residue(const md_system_t* sys, const backbone_chains_t* ch, md_temp_scope_t temp) {
    int32_t* res = md_temp_alloc_array(temp, int32_t, MD_MAX(1, ch->num));
    md_memset(res, -1, sizeof(int32_t) * MD_MAX(1, ch->num));
    if (sys->residue.count == 0) return res;
    uint32_t* count = md_temp_alloc_zero_array(temp, uint32_t, sys->residue.count);
    for (size_t u = 0; u < ch->num; ++u) {
        const md_atom_idx_t key = ch->aa ? ch->aa[u].ca : ch->na[u].c4;
        const md_residue_idx_t r = md_system_residue_find_by_atom_idx(sys, (size_t)key);
        if (r < 0) continue;
        const md_urange_t range = md_system_residue_atom_range(sys, (size_t)r);
        const bool whole = ch->aa ? (atom_in_range(ch->aa[u].n, range) && atom_in_range(ch->aa[u].c, range))
                                  : (atom_in_range(ch->na[u].o5, range) && atom_in_range(ch->na[u].o3, range));
        if (!whole) continue;
        res[u] = r;
        count[r] += 1;
    }
    for (size_t u = 0; u < ch->num; ++u) {
        if (res[u] >= 0 && count[res[u]] != 1) res[u] = -1;
    }
    return res;
}

static void perceive_polymer_atoms(md_system_t* sys, const md_system_state_t* state) {
    const size_t num_atoms = sys->atom.count;
    if (num_atoms == 0 || !sys->atom.type_idx) return;
    if (md_array_size(sys->atom.flags) < num_atoms) {
        const size_t size = md_array_size(sys->atom.flags);
        md_array_resize(sys->atom.flags, num_atoms, sys->alloc);
        md_memset(sys->atom.flags + size, 0, (num_atoms - size) * sizeof(md_atom_flags_t));
    }

    // Found anew: the roles of atoms are this step's alone (those of beads their model's)
    const uint32_t roles = MD_ATOM_FLAG_BACKBONE | MD_ATOM_FLAG_SIDE_CHAIN | MD_ATOM_FLAG_NUCLEOSIDE | MD_ATOM_FLAG_NUCLEOBASE | MD_ATOM_FLAG_TERMINAL_BEG | MD_ATOM_FLAG_TERMINAL_END;
    for (size_t i = 0; i < num_atoms; ++i) {
        if (md_atom_particle_kind(&sys->atom, i) == MD_PARTICLE_ATOM) {
            sys->atom.flags[i] = (md_atom_flags_t)((uint32_t)sys->atom.flags[i] & ~roles);
        }
    }
    if (sys->bond.count == 0) return;

    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    const md_urange_t all = { 0, (uint32_t)num_atoms };

    backbone_chains_t prot = {0};
    if (protein_chains_from_bonds(&prot, sys, state, temp)) {
        const int32_t* own = chain_units_own_residue(sys, &prot, temp);
        uint8_t* take = md_temp_alloc_zero_array(temp, uint8_t, prot.num);
        for (size_t u = 0; u < prot.num; ++u) {
            const md_amino_acid_atoms_t* a = &prot.aa[u];
            if (own[u] >= 0) {
                const md_urange_t range = md_system_residue_atom_range(sys, (size_t)own[u]);
                take[u] = md_resname_amino_acid(md_system_residue_name(sys, (size_t)own[u])) ||
                          atom_bonded_outside(sys, a->n, range, MD_Z_C) || atom_bonded_outside(sys, a->c, range, MD_Z_N);
            } else {
                take[u] = prot.prev[u] >= 0 || prot.next[u] >= 0;
            }
        }
        // The backbones of all first, so no side chain reaches into another's backbone over a disulfide
        for (size_t u = 0; u < prot.num; ++u) {
            if (!take[u]) continue;
            sys->atom.flags[prot.aa[u].n]  |= MD_ATOM_FLAG_BACKBONE;
            sys->atom.flags[prot.aa[u].ca] |= MD_ATOM_FLAG_BACKBONE;
            sys->atom.flags[prot.aa[u].c]  |= MD_ATOM_FLAG_BACKBONE;
        }
        for (size_t u = 0; u < prot.num; ++u) {
            if (!take[u]) continue;
            const md_amino_acid_atoms_t* a = &prot.aa[u];
            const bool in_residue = own[u] >= 0;
            const md_urange_t range = in_residue ? md_system_residue_atom_range(sys, (size_t)own[u]) : all;
            if (a->cb >= 0) {
                topo_floodfill_flag(sys, a->cb, a->ca, MD_ATOM_FLAG_SIDE_CHAIN, range, MD_ATOM_FLAG_BACKBONE, !in_residue);
            }
            // A free amino group (only hydrogens besides CA) begins the chain, a free carboxyl group (only oxygens
            // besides CA) ends it
            if (atom_is_chain_end(sys, a->n, a->ca, MD_Z_H)) {
                topo_floodfill_flag(sys, a->n, a->ca, MD_ATOM_FLAG_TERMINAL_BEG, range, MD_ATOM_FLAG_BACKBONE | MD_ATOM_FLAG_SIDE_CHAIN, !in_residue);
            }
            if (atom_is_chain_end(sys, a->c, a->ca, MD_Z_O)) {
                topo_floodfill_flag(sys, a->c, a->ca, MD_ATOM_FLAG_TERMINAL_END, range, MD_ATOM_FLAG_BACKBONE | MD_ATOM_FLAG_SIDE_CHAIN, !in_residue);
            }
        }
    }

    backbone_chains_t nucl = {0};
    if (nucleic_chains_from_bonds(&nucl, sys, temp)) {
        const int32_t* own = chain_units_own_residue(sys, &nucl, temp);
        uint8_t* take = md_temp_alloc_zero_array(temp, uint8_t, nucl.num);
        for (size_t u = 0; u < nucl.num; ++u) {
            const md_nucleic_acid_atoms_t* a = &nucl.na[u];
            if (own[u] >= 0) {
                const md_urange_t range = md_system_residue_atom_range(sys, (size_t)own[u]);
                take[u] = md_resname_nucleotide(md_system_residue_name(sys, (size_t)own[u])) ||
                          (a->p >= 0 && atom_bonded_outside(sys, a->p, range, MD_Z_O)) || atom_bonded_outside(sys, a->o3, range, MD_Z_P);
            } else {
                take[u] = nucl.prev[u] >= 0 || nucl.next[u] >= 0;
            }
        }
        for (size_t u = 0; u < nucl.num; ++u) {
            if (!take[u]) continue;
            const md_nucleic_acid_atoms_t* a = &nucl.na[u];
            if (a->p >= 0) sys->atom.flags[a->p] |= MD_ATOM_FLAG_BACKBONE;
            sys->atom.flags[a->o5] |= MD_ATOM_FLAG_BACKBONE;
            sys->atom.flags[a->c5] |= MD_ATOM_FLAG_BACKBONE;
            sys->atom.flags[a->c4] |= MD_ATOM_FLAG_BACKBONE | MD_ATOM_FLAG_NUCLEOSIDE;
            sys->atom.flags[a->c3] |= MD_ATOM_FLAG_BACKBONE | MD_ATOM_FLAG_NUCLEOSIDE;
            sys->atom.flags[a->o3] |= MD_ATOM_FLAG_BACKBONE;
            sys->atom.flags[a->c1] |= MD_ATOM_FLAG_NUCLEOSIDE;
            sys->atom.flags[a->c2] |= MD_ATOM_FLAG_NUCLEOSIDE;
            sys->atom.flags[a->o4] |= MD_ATOM_FLAG_NUCLEOSIDE;
        }
        for (size_t u = 0; u < nucl.num; ++u) {
            if (!take[u]) continue;
            const md_nucleic_acid_atoms_t* a = &nucl.na[u];
            const bool in_residue = own[u] >= 0;
            const md_urange_t range = in_residue ? md_system_residue_atom_range(sys, (size_t)own[u]) : all;
            // The base, from the atom of it bonded to the C1'
            const md_atom_idx_t base = mdi_bonded_atom_with(sys, a->c1, 0, 0, 0, a->c2, a->o4);
            if (base >= 0) {
                topo_floodfill_flag(sys, base, a->c1, MD_ATOM_FLAG_NUCLEOBASE | MD_ATOM_FLAG_NUCLEOSIDE, range, MD_ATOM_FLAG_BACKBONE, !in_residue);
            }
            // A free 3' hydroxyl ends the chain, a free 5' hydroxyl begins it
            if (atom_is_chain_end(sys, a->o3, a->c3, MD_Z_H)) {
                topo_floodfill_flag(sys, a->o3, a->c3, MD_ATOM_FLAG_TERMINAL_END, range, 0, !in_residue);
            }
            if (atom_is_chain_end(sys, a->o5, a->c5, MD_Z_H)) {
                topo_floodfill_flag(sys, a->o5, a->c5, MD_ATOM_FLAG_TERMINAL_BEG, range, 0, !in_residue);
            }
        }
    }

    md_temp_end(temp);
}

// The residues of the chains by the roles of their atoms (POLYMER ATOMS), as protein_chains_from_bonds: in the order of
// their CA, a residue following the one whose C is bonded to its N
static bool protein_chains_from_roles(backbone_chains_t* out, const md_system_t* sys, md_temp_scope_t temp) {
    const size_t num_atoms = sys->atom.count;
    if (!sys->atom.flags || num_atoms == 0) return false;
    md_amino_acid_atoms_t* units = NULL;
    size_t num = 0;
    for (size_t i = 0; i < num_atoms; ++i) {
        num += amino_acid_atoms_of_ca(NULL, sys, (md_atom_idx_t)i);
    }
    if (num == 0) return false;
    units = md_temp_alloc_array(temp, md_amino_acid_atoms_t, num);
    int32_t* unit_of_n = md_temp_alloc_array(temp, int32_t, num_atoms);
    md_memset(unit_of_n, -1, sizeof(int32_t) * num_atoms);
    size_t k = 0;
    for (size_t i = 0; i < num_atoms && k < num; ++i) {
        if (amino_acid_atoms_of_ca(&units[k], sys, (md_atom_idx_t)i)) {
            unit_of_n[units[k].n] = (int32_t)k;
            k += 1;
        }
    }
    int32_t* prev = md_temp_alloc_array(temp, int32_t, num);
    int32_t* next = md_temp_alloc_array(temp, int32_t, num);
    md_memset(prev, -1, sizeof(int32_t) * num);
    md_memset(next, -1, sizeof(int32_t) * num);
    for (size_t u = 0; u < num; ++u) {
        md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)units[u].c);
        for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
            const int32_t v = unit_of_n[md_bond_iter_atom_index(&it)];
            if (v >= 0 && v != (int32_t)u && prev[v] < 0) {
                next[u] = v;
                prev[v] = (int32_t)u;
                break;
            }
        }
    }
    uint32_t* order = md_temp_alloc_array(temp, uint32_t, num);
    uint32_t* begin = md_temp_alloc_array(temp, uint32_t, num + 1);
    out->num = num;
    out->aa = units;
    out->prev = prev;
    out->next = next;
    out->order = order;
    out->begin = begin;
    out->num_chains = backbone_chains(order, begin, prev, next, num, temp);
    return true;
}

// As protein_chains_from_roles, of the nucleotides: in the order of their C4', one following that whose O3' is bonded to
// its P
static bool nucleic_chains_from_roles(backbone_chains_t* out, const md_system_t* sys, md_temp_scope_t temp) {
    const size_t num_atoms = sys->atom.count;
    if (!sys->atom.flags || num_atoms == 0) return false;
    size_t num = 0;
    for (size_t i = 0; i < num_atoms; ++i) {
        num += nucleotide_atoms_of_c4(NULL, sys, (md_atom_idx_t)i);
    }
    if (num == 0) return false;
    md_nucleic_acid_atoms_t* units = md_temp_alloc_array(temp, md_nucleic_acid_atoms_t, num);
    int32_t* unit_of_p = md_temp_alloc_array(temp, int32_t, num_atoms);
    md_memset(unit_of_p, -1, sizeof(int32_t) * num_atoms);
    size_t k = 0;
    for (size_t i = 0; i < num_atoms && k < num; ++i) {
        if (nucleotide_atoms_of_c4(&units[k], sys, (md_atom_idx_t)i)) {
            if (units[k].p >= 0) unit_of_p[units[k].p] = (int32_t)k;
            k += 1;
        }
    }
    int32_t* prev = md_temp_alloc_array(temp, int32_t, num);
    int32_t* next = md_temp_alloc_array(temp, int32_t, num);
    md_memset(prev, -1, sizeof(int32_t) * num);
    md_memset(next, -1, sizeof(int32_t) * num);
    for (size_t u = 0; u < num; ++u) {
        md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)units[u].o3);
        for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
            const int32_t v = unit_of_p[md_bond_iter_atom_index(&it)];
            if (v >= 0 && v != (int32_t)u && prev[v] < 0) {
                next[u] = v;
                prev[v] = (int32_t)u;
                break;
            }
        }
    }
    uint32_t* order = md_temp_alloc_array(temp, uint32_t, num);
    uint32_t* begin = md_temp_alloc_array(temp, uint32_t, num + 1);
    out->num = num;
    out->na = units;
    out->prev = prev;
    out->next = next;
    out->order = order;
    out->begin = begin;
    out->num_chains = backbone_chains(order, begin, prev, next, num, temp);
    return true;
}

// The ends of coarse grained chains: an amino acid or a nucleotide of beads whose backbone is not bonded to that of the
// residue before (after) it begins (ends) its chain. Martini 3 gives the backbone bead at either end of a protein chain
// a regular type, charged (Q5) or neutral (P6, P1, P2) alike, where Gly, Ala, Val and Pro otherwise have a small one: a
// small backbone bead (54) at an end gets the mass and size of a regular one.
static void infer_coarse_grained_chain_ends(md_system_t* sys) {
    if (!md_system_is_coarse_grained(sys) || !sys->residue.flags || sys->bond.conn.offset == NULL) return;
    for (size_t ci = 0; ci < sys->residue.count; ++ci) {
        const md_residue_kind_t kind = md_system_residue_kind(sys, ci);
        if (kind != MD_RESIDUE_KIND_AMINO_ACID && kind != MD_RESIDUE_KIND_NUCLEOTIDE) continue;
        if (sys->residue.flags[ci] & MD_RESIDUE_FLAG_RESOLVED) continue;
        const md_urange_t range = md_system_residue_atom_range(sys, ci);
        if (range.beg >= range.end || md_atom_particle_kind(&sys->atom, range.beg) != MD_PARTICLE_BEAD) continue;
        const md_urange_t prev = ci > 0 ? md_system_residue_atom_range(sys, ci - 1) : (md_urange_t){0, 0};
        const md_urange_t next = ci + 1 < sys->residue.count ? md_system_residue_atom_range(sys, ci + 1) : (md_urange_t){0, 0};

        bool to_prev = false, to_next = false;
        for (uint32_t i = range.beg; i < range.end; ++i) {
            if (!(md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_BACKBONE)) continue;
            md_bond_iter_t it = md_bond_iter(&sys->bond, i);
            while (md_bond_iter_has_next(&it)) {
                const uint32_t j = (uint32_t)md_bond_iter_atom_index(&it);
                if (md_atom_flags(&sys->atom, j) & MD_ATOM_FLAG_BACKBONE) {
                    if (prev.beg <= j && j < prev.end) to_prev = true;
                    if (next.beg <= j && j < next.end) to_next = true;
                }
                md_bond_iter_next(&it);
            }
        }
        if (!to_prev) sys->residue.flags[ci] |= MD_RESIDUE_FLAG_TERMINAL_BEG;
        if (!to_next) sys->residue.flags[ci] |= MD_RESIDUE_FLAG_TERMINAL_END;

        if (kind == MD_RESIDUE_KIND_AMINO_ACID && (!to_prev || !to_next)) {
            for (uint32_t i = range.beg; i < range.end; ++i) {
                const md_particle_type_idx_t t = sys->atom.type_idx[i];
                if (md_particle_type_mass(&sys->atom.type, t) != 54.0f || !md_str_eq_cstr(md_particle_type_name(&sys->atom.type, t), "BB")) continue;
                sys->atom.type_idx[i] = mdi_find_or_add_particle_type(&sys->atom.type, md_particle_type_name(&sys->atom.type, t), 0, 72.0f, MARTINI_R,
                    md_particle_type_color(&sys->atom.type, t), md_particle_type_flags(&sys->atom.type, t), sys->alloc);
            }
        }
    }
}

// Classifies the residues (amino acid, nucleotide, water, ion, lipid) from their names, atoms and bonds, and gives the
// atoms of the amino acids and nucleotides whose backbone could be verified their roles.
// What residue res_idx is (md_system_infer_residue_flags), fallback where nothing tells. The flags it has besides
// (RESOLVED, the terminal flags) go to out_flags; ambiguous is 1 (2) for an amino acid (nucleotide) by name alone.
static md_residue_kind_t classify_residue(const md_system_t* sys, size_t res_idx, md_residue_kind_t fallback, md_residue_flags_t* out_flags, int* ambiguous) {
    const md_str_t res_name = md_system_residue_name(sys, res_idx);
    const md_urange_t res_range = md_system_residue_atom_range(sys, res_idx);
    const size_t len = (size_t)(res_range.end - res_range.beg);

    // A molecule of the Martini force fields (a lipid, water, an ion, a sugar, ...): the kind of its topology. Its
    // name alone may say something else: PCA, a lipid of Martini 2, is pyroglutamic acid in the PDB.
    if (len <= 255 && atoms_are_beads(sys, res_range)) {
        md_str_t labels[256];
        uint8_t map[256];
        for (size_t i = 0; i < len; ++i) labels[i] = md_atom_name(&sys->atom, res_range.beg + i);
        const martini_db_residue_t* mol = mdi_martini_db_match_molecule(MARTINI_ANY, res_name, labels, len, map);
        // Renamed (CHL for CHOL): by its beads, unless it is named as a residue of a chain
        if (!mol && !md_resname_amino_acid(res_name) && !md_resname_nucleotide(res_name)) {
            mol = mdi_martini_db_match_beads(MARTINI_ANY, labels, NULL, len, map);
        }
        if (mol) return mol->kind;
    }

    // Amino acids and nucleotides are classified whatever their kind, as a predefined coarse grained type may have
    // given the kind but the atoms of an atomistic one still have their roles: a residue which holds the backbone of
    // one amino acid (nucleotide) of a chain is one (POLYMER ATOMS)
    md_amino_acid_atoms_t prot_atoms;
    if (mdi_residue_amino_acid_atoms(&prot_atoms, sys, res_range)) {
        *out_flags |= MD_RESIDUE_FLAG_RESOLVED;
        if (md_atom_flags(&sys->atom, (size_t)prot_atoms.n) & MD_ATOM_FLAG_TERMINAL_BEG) *out_flags |= MD_RESIDUE_FLAG_TERMINAL_BEG;
        if (md_atom_flags(&sys->atom, (size_t)prot_atoms.c) & MD_ATOM_FLAG_TERMINAL_END) *out_flags |= MD_RESIDUE_FLAG_TERMINAL_END;
        return MD_RESIDUE_KIND_AMINO_ACID;
    }
    if (md_resname_amino_acid(res_name)) {
        *ambiguous = 1;
        return MD_RESIDUE_KIND_AMINO_ACID;
    }

    md_nucleic_acid_atoms_t nucl_atoms;
    if (mdi_residue_nucleotide_atoms(&nucl_atoms, sys, res_range)) {
        *out_flags |= MD_RESIDUE_FLAG_RESOLVED;
        if (md_atom_flags(&sys->atom, (size_t)nucl_atoms.o5) & MD_ATOM_FLAG_TERMINAL_BEG) *out_flags |= MD_RESIDUE_FLAG_TERMINAL_BEG;
        if (md_atom_flags(&sys->atom, (size_t)nucl_atoms.o3) & MD_ATOM_FLAG_TERMINAL_END) *out_flags |= MD_RESIDUE_FLAG_TERMINAL_END;
        return MD_RESIDUE_KIND_NUCLEOTIDE;
    }
    if (md_resname_nucleotide(res_name)) {
        *ambiguous = 2;
        return MD_RESIDUE_KIND_NUCLEOTIDE;
    }

    // Water by its atoms, or by name when its hydrogens are left out (united atom, coarse grained)
    if (atoms_are_water(sys, res_range) ||
       ((len == 1 || len == 3) && md_resname_water(res_name)) ||
       (len == 1 && md_resname_water(md_atom_name(&sys->atom, res_range.beg)))) {
        return MD_RESIDUE_KIND_WATER;
    }
    // An ion by its name, or by its element where it is bonded to nothing but by coordination: a lone O, S or N
    // covalently bonded (the O between two Si of a silica network, a residue to an atom) is no ion
    if (len == 1 && (md_resname_ion(res_name) ||
        (monatomic_ion_element(md_atom_atomic_number(&sys->atom, res_range.beg)) && !atom_bonded_covalently(sys, (md_atom_idx_t)res_range.beg)))) {
        return MD_RESIDUE_KIND_ION;
    }
    if (md_resname_lipid(res_name)) {
        return MD_RESIDUE_KIND_LIPID;
    }
    return fallback;
}

bool md_system_infer_residue_flags(md_system_t* sys, const md_system_state_t* state) {
    if (!sys) {
        MD_LOG_ERROR("Missing system");
        return false;
    }
    // The atoms of the chains first: what follows reads them (POLYMER ATOMS)
    if (sys->atom.count) {
        MD_ASSERT(sys->alloc);
        perceive_polymer_atoms(sys, state);
    }
    if (sys->residue.count == 0) {
        return true;
    }
    MD_ASSERT(sys->alloc);
    md_system_topology_changed(sys);

    if (md_array_size(sys->residue.flags) < sys->residue.count) {
        const size_t size = md_array_size(sys->residue.flags);
        md_array_resize(sys->residue.flags, sys->residue.count, sys->alloc);
        md_memset(sys->residue.flags + size, 0, (sys->residue.count - size) * sizeof(md_residue_flags_t));
    }
    if (md_array_size(sys->atom.flags) < sys->atom.count) {
        const size_t size = md_array_size(sys->atom.flags);
        md_array_resize(sys->atom.flags, sys->atom.count, sys->alloc);
        md_memset(sys->atom.flags + size, 0, (sys->atom.count - size) * sizeof(md_atom_flags_t));
    }

    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    md_array(int) ambigous_amino_acid_indices = 0;
    md_array(int) ambigous_nucleotide_indices = 0;

    // What the classification finds is found anew: the flags it sets are cleared first (the breaks are the backbones')
    for (size_t i = 0; i < sys->residue.count; ++i) {
        sys->residue.flags[i] = (md_residue_flags_t)((uint32_t)sys->residue.flags[i] & ~(uint32_t)(MD_RESIDUE_FLAG_RESOLVED | MD_RESIDUE_FLAG_TERMINAL_BEG | MD_RESIDUE_FLAG_TERMINAL_END));
    }

    for (size_t residue_idx = 0; residue_idx < sys->residue.count; ++residue_idx) {
        const md_urange_t res_range = md_system_residue_atom_range(sys, residue_idx);
        const bool beads = atoms_are_beads(sys, res_range);

        // A residue of beads a predefined coarse grained type classified keeps its kind (what the model says it is,
        // whatever it is named). One of atoms is classified anew each time from what the system is now, and is OTHER
        // where nothing tells what it is: a kind found before (an ion before the bonds joined it to a network, an
        // amino acid before its bonds were edited) is not kept.
        const md_residue_kind_t kind = md_system_residue_kind(sys, residue_idx);
        if (beads && (kind == MD_RESIDUE_KIND_WATER || kind == MD_RESIDUE_KIND_ION || kind == MD_RESIDUE_KIND_LIPID)) {
            continue;
        }
        md_residue_flags_t flags = MD_RESIDUE_FLAG_NONE;
        int ambiguous = 0;
        const md_residue_kind_t found = classify_residue(sys, residue_idx, beads ? kind : MD_RESIDUE_KIND_OTHER, &flags, &ambiguous);
        set_residue_kind(sys, residue_idx, found);
        sys->residue.flags[residue_idx] |= flags;
        if (ambiguous == 1) md_array_push(ambigous_amino_acid_indices, (int)residue_idx, temp_alloc);
        if (ambiguous == 2) md_array_push(ambigous_nucleotide_indices, (int)residue_idx, temp_alloc);
    }

    infer_coarse_grained_chain_ends(sys);

    const size_t max_warnings = 10;
    
    // Print warnings about unsafe assignments
    size_t num_ambigous_amino_acids = md_array_size(ambigous_amino_acid_indices);
    if (num_ambigous_amino_acids > 0) {
        MD_LOG_DEBUG("Structure contains %zu amino acid residue names, but their backbone cannot be mapped from atom labels:", num_ambigous_amino_acids);
        for (size_t i = 0; i < MD_MIN(num_ambigous_amino_acids, max_warnings); ++i) {
            int res_idx = ambigous_amino_acid_indices[i];
            MD_LOG_DEBUG(" - [%d] " MD_STR_FMT " (%d)", res_idx + 1, MD_STR_ARG(md_system_residue_name(sys, res_idx)), md_system_residue_seq_id(sys, res_idx));
        }
        if (num_ambigous_amino_acids > max_warnings) {
            MD_LOG_DEBUG(" - And %zu more...", num_ambigous_amino_acids - max_warnings);
        }
    }

    size_t num_ambigous_nucleotides = md_array_size(ambigous_nucleotide_indices);
    if (num_ambigous_nucleotides > 0) {
        MD_LOG_DEBUG("Structure contains %zu nucleotide residue names, but their backbone cannot be mapped from atom labels:", num_ambigous_nucleotides);
        for (size_t i = 0; i < MD_MIN(num_ambigous_nucleotides, max_warnings); ++i) {
            int res_idx = ambigous_nucleotide_indices[i];
            MD_LOG_DEBUG(" - [%d] " MD_STR_FMT " (%d)", res_idx + 1, MD_STR_ARG(md_system_residue_name(sys, res_idx)), md_system_residue_seq_id(sys, res_idx));
        }
        if (num_ambigous_nucleotides > max_warnings) {
            MD_LOG_DEBUG(" - And %zu more...", num_ambigous_nucleotides - max_warnings);
        }
    }
    
    md_temp_end(temp);
    return true;
}


// The segments of one range: their residues res_idx, or consecutive residues from res_base without it (-1: none)
static inline void commit_protein_backbone(const md_amino_acid_atoms_t* bb_atoms, size_t bb_length, md_residue_idx_t res_base, const md_residue_idx_t* res_idx_of, md_part_idx_t part_idx, md_protein_backbone_data_t* bb_data, md_allocator_i* alloc) {
    uint32_t offset = (uint32_t)md_array_size(bb_data->segment.atoms);
    md_array_push_array(bb_data->segment.atoms, bb_atoms, bb_length, alloc);
    md_array_push(bb_data->range.offset, offset, alloc);
    md_array_push(bb_data->range.part_idx, part_idx, alloc);

    for (size_t i = 0; i < bb_length; ++i) {
        const md_residue_idx_t res_idx = res_idx_of ? res_idx_of[i] : res_base < 0 ? -1 : res_base + (md_residue_idx_t)i;
        md_array_push(bb_data->segment.res_idx, res_idx, alloc);
    }
}

// The bead which stands for the backbone of an amino acid of a coarse grained model (the BB bead of Martini): its one
// bead with the backbone role. -1 for a residue which is no amino acid, or has none or several.
static md_atom_idx_t amino_acid_backbone_bead(const md_system_t* sys, size_t res_idx) {
    if (md_system_residue_kind(sys, res_idx) != MD_RESIDUE_KIND_AMINO_ACID) return -1;
    const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
    md_atom_idx_t bead = -1;
    for (uint32_t i = range.beg; i < range.end; ++i) {
        if ((md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_BACKBONE) && md_atom_particle_kind(&sys->atom, i) == MD_PARTICLE_BEAD) {
            if (bead >= 0) return -1;
            bead = (md_atom_idx_t)i;
        }
    }
    return bead;
}

// The backbone beads of a nucleotide of a coarse grained model, in the fields of the atoms they hold (Martini 2 DNA and
// RNA, Uusitalo et al. 2015, 2017; Martini 3 RNA of reForge): the phosphate bead BB1 in p, the sugar bead BB2 which
// holds the C4' in c4 and the sugar bead BB3 which holds the C1' in c1, the other fields -1. A 5' end has no BB1. False
// for a residue which is no nucleotide, or lacks BB2 or BB3. Their o5 is -1, which that of atoms never is
// (nucleotide_atoms_of_c4 requires it): that tells beads from atoms (nucleic_atoms_are_beads).
static bool nucleotide_backbone_beads(md_nucleic_acid_atoms_t* out, const md_system_t* sys, size_t res_idx) {
    if (md_system_residue_kind(sys, res_idx) != MD_RESIDUE_KIND_NUCLEOTIDE) return false;
    md_nucleic_acid_atoms_t na;
    md_memset(&na, -1, sizeof(na));
    const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
    for (uint32_t i = range.beg; i < range.end; ++i) {
        if (!(md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_BACKBONE) || md_atom_particle_kind(&sys->atom, i) != MD_PARTICLE_BEAD) continue;
        const md_str_t name = md_atom_name(&sys->atom, i);
        if      (md_str_eq_cstr(name, "BB1")) na.p  = (md_atom_idx_t)i;
        else if (md_str_eq_cstr(name, "BB2")) na.c4 = (md_atom_idx_t)i;
        else if (md_str_eq_cstr(name, "BB3")) na.c1 = (md_atom_idx_t)i;
    }
    if (out) *out = na;
    return na.c4 >= 0 && na.c1 >= 0;
}

static inline bool nucleic_atoms_are_beads(const md_nucleic_acid_atoms_t* na) {
    return na->o5 < 0;
}

// ### CHAIN BREAKS ### (see md_util.h)

// Whether the beads a and b of consecutive residues of a coarse grained chain link their backbones. The bonds decide when
// the system has bonds for both; beads without bonds are linked when their sequence numbers are consecutive (seq_next),
// not by their distance: beads are bonded farther apart than atoms are. A bead which is not there (-1) links nothing.
// The backbones of atoms are linked by their bonds alone (BACKBONES FROM THE BONDS).
static bool chain_beads_linked(const md_system_t* sys, md_atom_idx_t a, md_atom_idx_t b, bool seq_next) {
    if (a < 0 || b < 0) return false;
    const md_bond_iter_t it_a = md_bond_iter(&sys->bond, (size_t)a);
    const md_bond_iter_t it_b = md_bond_iter(&sys->bond, (size_t)b);
    if (md_bond_iter_has_next(&it_a) && md_bond_iter_has_next(&it_b)) {
        return md_bond_find(&sys->bond, a, b) >= 0;
    }
    return seq_next;
}

// Whether the consecutive amino acids of a trace with the backbone beads prev and next are linked: by a bond between them
static bool protein_trace_linked(const md_system_t* sys, const md_amino_acid_atoms_t* prev, const md_amino_acid_atoms_t* next, bool seq_next) {
    return chain_beads_linked(sys, prev->ca, next->ca, seq_next);
}

// Whether the consecutive nucleotides of beads prev and next are linked: by a bond from a sugar bead of prev to the
// phosphate bead of next (BB3 to BB1 in Martini 2, BB2 and BB3 to BB1 in Martini 3). Without its phosphate bead (a 5'
// end) next is linked to nothing.
static bool nucleic_beads_linked(const md_system_t* sys, const md_nucleic_acid_atoms_t* prev, const md_nucleic_acid_atoms_t* next, bool seq_next) {
    return chain_beads_linked(sys, prev->c1, next->p, seq_next) || chain_beads_linked(sys, prev->c4, next->p, seq_next);
}

// The atoms of element z bonded to the atom, and whether it has bonds at all
static size_t bonded_count_of_element(const md_system_t* sys, md_atom_idx_t atom, md_atomic_number_t z, bool* has_bonds) {
    size_t count = 0;
    *has_bonds = false;
    if (atom < 0) return 0;
    md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)atom);
    while (md_bond_iter_has_next(&it)) {
        *has_bonds = true;
        count += md_atom_atomic_number(&sys->atom, md_bond_iter_atom_index(&it)) == z;
        md_bond_iter_next(&it);
    }
    return count;
}

// Whether the carboxyl group of an amino acid (its C) is free: its C bonded to two oxygens, O and OXT. Without bonds, by
// the name of the second oxygen.
static bool amino_acid_carboxyl_free(const md_system_t* sys, size_t res_idx, md_atom_idx_t c) {
    (void)res_idx;
    bool has_bonds;
    return bonded_count_of_element(sys, c, MD_Z_O, &has_bonds) >= 2;
}

// Whether the chain ends between the consecutive residues a and b of a protein chain whose backbones are not linked,
// rather than going on over a break. Only what is unambiguous counts, since a crystal structure rarely has the
// hydrogens and often leaves out the OXT: a free carboxyl group on a (O and OXT), a free amino group on b (two or
// three hydrogens on its N), or b a lone amino acid with a free carboxyl group of its own: a free amino acid (the
// arginine bound to 1LAF), which a PDB file can give the chain id of the chain and the next number. The beads of a
// trace (coarse grained) tell none of this.
static bool protein_chain_ends_between(const md_system_t* sys, size_t a, size_t b, const md_amino_acid_atoms_t* atoms_a, const md_amino_acid_atoms_t* atoms_b) {
    if (atoms_a->n < 0) return false;
    if (amino_acid_carboxyl_free(sys, a, atoms_a->c)) return true;
    if (amino_acid_carboxyl_free(sys, b, atoms_b->c)) return true;
    bool has_bonds;
    return bonded_count_of_element(sys, atoms_b->n, MD_Z_H, &has_bonds) >= 2;
}

// Whether the chain ends between the consecutive residues of a nucleic acid chain whose backbones are not linked: b
// has no P (a 5' end, which a crystal structure gives the first residue of a strand), or a hydrogen on the O3' of a or
// the O5' of b makes a free 3' or 5' hydroxyl. Beads tell only the first: a b without its phosphate bead.
static bool nucleic_chain_ends_between(const md_system_t* sys, const md_nucleic_acid_atoms_t* atoms_a, const md_nucleic_acid_atoms_t* atoms_b) {
    if (nucleic_atoms_are_beads(atoms_a)) return atoms_b->p < 0;
    if (atoms_a->o3 < 0 || atoms_b->p < 0) return true;
    bool has_bonds;
    return bonded_count_of_element(sys, atoms_a->o3, MD_Z_H, &has_bonds) > 0 || bonded_count_of_element(sys, atoms_b->o5, MD_Z_H, &has_bonds) > 0;
}

// Whether a polymer chain ends between the consecutive residues a and b = a + 1, rather than going on over a break, by
// what their atoms or beads tell (CHAIN BREAKS in md_backbone.h): the amino acids by protein_chain_ends_between, the
// nucleotides of atoms or of beads by nucleic_chain_ends_between. Whether their backbones are linked is not asked.
bool mdi_polymer_chain_ends_between(const md_system_t* sys, size_t a, size_t b) {
    const md_urange_t range_a = md_system_residue_atom_range(sys, a);
    const md_urange_t range_b = md_system_residue_atom_range(sys, b);
    md_amino_acid_atoms_t aa_a, aa_b;
    md_nucleic_acid_atoms_t na_a, na_b;
    if (mdi_residue_amino_acid_atoms(&aa_a, sys, range_a) && mdi_residue_amino_acid_atoms(&aa_b, sys, range_b)) {
        return protein_chain_ends_between(sys, a, b, &aa_a, &aa_b);
    }
    if (mdi_residue_nucleotide_atoms(&na_a, sys, range_a) && mdi_residue_nucleotide_atoms(&na_b, sys, range_b)) {
        return nucleic_chain_ends_between(sys, &na_a, &na_b);
    }
    if (nucleotide_backbone_beads(&na_a, sys, a) && nucleotide_backbone_beads(&na_b, sys, b)) {
        return nucleic_chain_ends_between(sys, &na_a, &na_b);
    }
    return false;
}

// Marks the break between the consecutive residues a and b = a + 1 of a chain. The chain goes on past the break, so
// neither residue ends it: what the classification took for a terminus there (the carboxyl group of a residue whose
// successor is missing, the amino group of one whose predecessor is, a hydroxyl on either side of a nucleic acid
// break) is not one, and loses its terminal flags.
// A break of the chain between the residues a and b = a + 1 (CHAIN BREAKS in md_backbone.h)
static int cmp_residue_idx(const void* a, const void* b) {
    const md_residue_idx_t x = *(const md_residue_idx_t*)a, y = *(const md_residue_idx_t*)b;
    return (x > y) - (x < y);
}

static void mark_chain_break(md_backbone_data_t* out, size_t a, size_t b, md_allocator_i* alloc) {
    MD_ASSERT(b == a + 1);
    (void)b;
    md_array_push(out->breaks, (md_residue_idx_t)a, alloc);
}

size_t md_system_protein_segment_atoms(md_atom_idx_t out[], size_t cap, const md_system_t* sys, size_t seg_idx) {
    MD_ASSERT(sys);
    const md_protein_backbone_data_t* bb = &sys->protein_backbone;
    if (seg_idx >= bb->segment.count || !bb->segment.atoms) return 0;
    const md_residue_idx_t res_idx = bb->segment.res_idx ? bb->segment.res_idx[seg_idx] : -1;
    if (res_idx >= 0 && (size_t)res_idx < sys->residue.count) {
        const md_urange_t range = md_system_residue_atom_range(sys, (size_t)res_idx);
        for (uint32_t i = range.beg; i < range.end && i - range.beg < cap; ++i) out[i - range.beg] = (md_atom_idx_t)i;
        return range.end - range.beg;
    }

    const md_amino_acid_atoms_t* a = &bb->segment.atoms[seg_idx];
    if (a->ca < 0) return 0;
    if (a->n < 0) {
        // A bead of a trace stands for itself
        if (cap > 0) out[0] = a->ca;
        return 1;
    }
    enum { MAX_ATOMS = 512 };
    md_atom_idx_t found[MAX_ATOMS];
    size_t num = 0;
    found[num++] = a->ca;
    for (size_t k = 0; k < num; ++k) {
        const md_atom_idx_t x = found[k];
        const bool x_sulfur = md_atom_atomic_number(&sys->atom, (size_t)x) == MD_Z_S;
        md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)x);
        for (; md_bond_iter_has_next(&it) && num < (size_t)MAX_ATOMS; md_bond_iter_next(&it)) {
            const md_atom_idx_t y = md_bond_iter_atom_index(&it);
            if (md_bond_kind((md_bond_flags_t)md_bond_iter_bond_flags(&it)) != MD_BOND_KIND_COVALENT) continue;
            if ((md_atom_flags(&sys->atom, (size_t)y) & MD_ATOM_FLAG_BACKBONE) && y != a->n && y != a->c) continue;
            if (x_sulfur && md_atom_atomic_number(&sys->atom, (size_t)y) == MD_Z_S) continue;
            bool seen = false;
            for (size_t j = 0; j < num && !seen; ++j) seen = found[j] == y;
            if (!seen) found[num++] = y;
        }
    }
    for (size_t k = 0; k < num && k < cap; ++k) out[k] = found[k];
    return num;
}

bool md_system_chain_break_atoms(md_atom_idx_t out_atoms[2], const md_system_t* sys, size_t res_idx) {
    MD_ASSERT(out_atoms);
    out_atoms[0] = -1;
    out_atoms[1] = -1;
    if (!sys || res_idx + 1 >= sys->residue.count) return false;
    if (!(md_system_residue_flags(sys, res_idx) & MD_RESIDUE_FLAG_BREAK_AFTER)) return false;

    // The atoms of the link the backbones were found by (md_backbone_infer): of the segments of the two
    // residues, the C and N of amino acids (the beads of a trace), else the O3' and P of nucleotides (the BB3 and the
    // phosphate bead of beads)
    const size_t next = res_idx + 1;
    md_atom_idx_t a = -1;
    md_atom_idx_t b = -1;
    const md_protein_backbone_data_t* pb = &sys->protein_backbone;
    int32_t sa = -1, sb = -1;
    for (size_t s = 0; pb->segment.res_idx && s < pb->segment.count; ++s) {
        if (pb->segment.res_idx[s] == (md_residue_idx_t)res_idx) sa = (int32_t)s;
        if (pb->segment.res_idx[s] == (md_residue_idx_t)next)    sb = (int32_t)s;
    }
    if (sa >= 0 && sb >= 0) {
        const md_amino_acid_atoms_t* x = &pb->segment.atoms[sa];
        const md_amino_acid_atoms_t* y = &pb->segment.atoms[sb];
        a = x->n < 0 ? x->ca : x->c;
        b = y->n < 0 ? y->ca : y->n;
    } else {
        const md_nucleic_backbone_data_t* nb = &sys->nucleic_backbone;
        sa = sb = -1;
        for (size_t s = 0; nb->segment.res_idx && s < nb->segment.count; ++s) {
            if (nb->segment.res_idx[s] == (md_residue_idx_t)res_idx) sa = (int32_t)s;
            if (nb->segment.res_idx[s] == (md_residue_idx_t)next)    sb = (int32_t)s;
        }
        if (sa >= 0 && sb >= 0) {
            const md_nucleic_acid_atoms_t* x = &nb->segment.atoms[sa];
            const md_nucleic_acid_atoms_t* y = &nb->segment.atoms[sb];
            a = nucleic_atoms_are_beads(x) ? x->c1 : x->o3;
            b = y->p;
        }
    }
    if (a < 0 || b < 0) return false;
    out_atoms[0] = a;
    out_atoms[1] = b;
    return true;
}

static inline void flush_protein_backbone_segment(md_amino_acid_atoms_t* backbone_atoms, md_residue_idx_t* res_base, size_t part_idx, md_protein_backbone_data_t* bb_data, md_allocator_i* alloc, size_t min_backbone_length) {
    size_t backbone_length = md_array_size(backbone_atoms);
    if (backbone_length >= min_backbone_length && *res_base != -1) {
        commit_protein_backbone(backbone_atoms, backbone_length, *res_base, NULL, (md_part_idx_t)part_idx, bb_data, alloc);
    }
    md_array_shrink(backbone_atoms, 0);
    *res_base = -1;
}

// As commit_protein_backbone
static inline void commit_nucleic_backbone(const md_nucleic_acid_atoms_t* bb_atoms, size_t bb_length, md_residue_idx_t res_base, const md_residue_idx_t* res_idx_of, md_part_idx_t part_idx, md_nucleic_backbone_data_t* bb_data, md_allocator_i* alloc) {
    uint32_t offset = (uint32_t)md_array_size(bb_data->segment.atoms);
    md_array_push_array(bb_data->segment.atoms, bb_atoms, bb_length, alloc);
    md_array_push(bb_data->range.offset, offset, alloc);
    md_array_push(bb_data->range.part_idx, part_idx, alloc);

    for (size_t i = 0; i < bb_length; ++i) {
        const md_residue_idx_t res_idx = res_idx_of ? res_idx_of[i] : res_base < 0 ? -1 : res_base + (md_residue_idx_t)i;
        md_array_push(bb_data->segment.res_idx, res_idx, alloc);
    }
}

static inline void flush_nucleic_backbone_segment(md_nucleic_acid_atoms_t* backbone_atoms, md_residue_idx_t* res_base, size_t part_idx, md_nucleic_backbone_data_t* bb_data, md_allocator_i* alloc, size_t min_backbone_length) {
    size_t backbone_length = md_array_size(backbone_atoms);
    if (backbone_length >= min_backbone_length && *res_base != -1) {
        commit_nucleic_backbone(backbone_atoms, backbone_length, *res_base, NULL, (md_part_idx_t)part_idx, bb_data, alloc);
    }
    md_array_shrink(backbone_atoms, 0);
    *res_base = -1;
}

// The residue of each residue of the chains (RESIDUES above): of its CA (C4'), -1 where that residue holds the CA (C4')
// of another as well, or the system has none
static int32_t* chain_units_residue(const md_system_t* sys, const backbone_chains_t* ch, md_temp_scope_t temp) {
    int32_t* res = md_temp_alloc_array(temp, int32_t, MD_MAX(1, ch->num));
    md_memset(res, -1, sizeof(int32_t) * MD_MAX(1, ch->num));
    if (sys->residue.count == 0) return res;
    uint32_t* count = md_temp_alloc_zero_array(temp, uint32_t, sys->residue.count);
    for (size_t u = 0; u < ch->num; ++u) {
        const md_atom_idx_t key = ch->aa ? ch->aa[u].ca : ch->na[u].c4;
        res[u] = md_system_residue_find_by_atom_idx(sys, (size_t)key);
        if (res[u] >= 0) count[res[u]] += 1;
    }
    for (size_t u = 0; u < ch->num; ++u) {
        if (res[u] >= 0 && count[res[u]] != 1) res[u] = -1;
    }
    return res;
}

static inline md_part_idx_t residue_part(const md_system_t* sys, md_residue_idx_t res_idx) {
    return (res_idx >= 0 && sys->part.count) ? md_system_part_find_by_residue_idx(sys, (size_t)res_idx) : -1;
}

// Commits the chains of at least min_length residues as ranges, cut where the part of their residues changes
static void commit_backbone_chains(md_backbone_data_t* out, const md_system_t* sys, const backbone_chains_t* ch, const int32_t* unit_res, size_t min_length, md_temp_scope_t temp, md_allocator_i* alloc) {
    md_amino_acid_atoms_t*   aa  = ch->aa ? md_temp_alloc_array(temp, md_amino_acid_atoms_t, ch->num) : NULL;
    md_nucleic_acid_atoms_t* na  = ch->na ? md_temp_alloc_array(temp, md_nucleic_acid_atoms_t, ch->num) : NULL;
    md_residue_idx_t*        res = md_temp_alloc_array(temp, md_residue_idx_t, MD_MAX(1, ch->num));
    for (size_t c = 0; c < ch->num_chains; ++c) {
        size_t len = 0;
        md_part_idx_t part = -1;
        for (uint32_t j = ch->begin[c]; j <= ch->begin[c + 1]; ++j) {
            const bool end = j == ch->begin[c + 1];
            const uint32_t u = end ? 0 : ch->order[j];
            const md_part_idx_t p = end ? -1 : residue_part(sys, unit_res[u]);
            if (len > 0 && (end || p != part)) {
                if (len >= min_length) {
                    if (aa) commit_protein_backbone(aa, len, -1, res, part, &out->protein, alloc);
                    else    commit_nucleic_backbone(na, len, -1, res, part, &out->nucleic, alloc);
                }
                len = 0;
            }
            if (end) break;
            part = p;
            if (aa) aa[len] = ch->aa[u];
            else    na[len] = ch->na[u];
            res[len] = unit_res[u];
            len += 1;
        }
    }
}

// Marks the chain breaks between consecutive residues of a part which are not linked, unless the chain ends there
// (CHAIN BREAKS in md_backbone.h). Both are residues of the chains which are of the polymer: in a part of its kind (a
// peptide, a nucleic acid), or on a chain of at least min_length residues where the residues are not classified (a
// LAMMPS molecule of residues named "MOL"). The groups which match without being residues (N7-C5-C6=O6 of a guanine)
// are neither.
static void mark_backbone_chain_breaks(md_backbone_data_t* out, const md_system_t* sys, const backbone_chains_t* ch, const int32_t* unit_res, size_t min_length, md_temp_scope_t temp, md_allocator_i* alloc) {
    if (sys->residue.count < 2 || sys->part.count == 0) return;
    uint32_t* chain_len = md_temp_alloc_array(temp, uint32_t, MD_MAX(1, ch->num));
    for (size_t c = 0; c < ch->num_chains; ++c) {
        for (uint32_t j = ch->begin[c]; j < ch->begin[c + 1]; ++j) chain_len[ch->order[j]] = ch->begin[c + 1] - ch->begin[c];
    }
    int32_t* res_unit = md_temp_alloc_array(temp, int32_t, sys->residue.count);
    md_memset(res_unit, -1, sizeof(int32_t) * sys->residue.count);
    for (size_t u = 0; u < ch->num; ++u) {
        const md_residue_idx_t r = unit_res[u];
        if (r < 0) continue;
        const md_part_idx_t part = residue_part(sys, r);
        const md_part_kind_t kind = part >= 0 ? md_system_part_kind(sys, (size_t)part) : MD_PART_KIND_UNKNOWN;
        const bool polymer = ch->aa ? kind == MD_PART_KIND_PEPTIDE : md_part_kind_is_nucleic_acid(kind);
        if (polymer || chain_len[u] >= min_length) res_unit[r] = (int32_t)u;
    }
    for (size_t r = 0; r + 1 < sys->residue.count; ++r) {
        const int32_t u = res_unit[r];
        const int32_t v = res_unit[r + 1];
        if (u < 0 || v < 0 || ch->next[u] == v) continue;
        const md_part_idx_t part = residue_part(sys, (md_residue_idx_t)r);
        if (part < 0 || part != residue_part(sys, (md_residue_idx_t)(r + 1))) continue;
        const bool ends = ch->aa ? protein_chain_ends_between(sys, r, r + 1, &ch->aa[u], &ch->aa[v]) : nucleic_chain_ends_between(sys, &ch->na[u], &ch->na[v]);
        if (!ends) mark_chain_break(out, r, r + 1, alloc);
    }
}

// Orders the ranges of the backbone by the first atom of each (CA or C4'): those found from the bonds and the traces
// of the beads, which come from two passes, as the residues are along each chain
#define SORT_BACKBONE_RANGES(NAME, BB_T, ATOM_T, KEY)                                                                       \
static void NAME(BB_T* bb, md_temp_scope_t temp) {                                                                         \
    const size_t num_ranges = md_array_size(bb->range.offset);                                                             \
    const size_t num_segs = md_array_size(bb->segment.atoms);                                                              \
    if (num_ranges < 2) return;                                                                                             \
    uint64_t* key = md_temp_alloc_array(temp, uint64_t, num_ranges);                                                        \
    for (size_t r = 0; r < num_ranges; ++r) key[r] = (uint64_t)(uint32_t)bb->segment.atoms[bb->range.offset[r]].KEY << 32 | r; \
    qsort(key, num_ranges, sizeof(uint64_t), cmp_u64);                                                                      \
    uint32_t* off = md_temp_alloc_array(temp, uint32_t, num_ranges);                                                        \
    md_part_idx_t* part = md_temp_alloc_array(temp, md_part_idx_t, num_ranges);                                             \
    ATOM_T* atoms = md_temp_alloc_array(temp, ATOM_T, MD_MAX(1, num_segs));                                                 \
    md_residue_idx_t* res = md_temp_alloc_array(temp, md_residue_idx_t, MD_MAX(1, num_segs));                               \
    md_memcpy(off, bb->range.offset, sizeof(uint32_t) * num_ranges);                                                        \
    md_memcpy(part, bb->range.part_idx, sizeof(md_part_idx_t) * num_ranges);                                                \
    md_memcpy(atoms, bb->segment.atoms, sizeof(ATOM_T) * num_segs);                                                         \
    md_memcpy(res, bb->segment.res_idx, sizeof(md_residue_idx_t) * num_segs);                                               \
    uint32_t w = 0;                                                                                                         \
    for (size_t k = 0; k < num_ranges; ++k) {                                                                               \
        const uint32_t r = (uint32_t)key[k];                                                                                \
        const uint32_t beg = off[r];                                                                                        \
        const uint32_t end = r + 1 < num_ranges ? off[r + 1] : (uint32_t)num_segs;                                          \
        bb->range.offset[k] = w;                                                                                            \
        bb->range.part_idx[k] = part[r];                                                                                    \
        md_memcpy(bb->segment.atoms + w, atoms + beg, sizeof(ATOM_T) * (end - beg));                                        \
        md_memcpy(bb->segment.res_idx + w, res + beg, sizeof(md_residue_idx_t) * (end - beg));                              \
        w += end - beg;                                                                                                     \
    }                                                                                                                       \
}

SORT_BACKBONE_RANGES(sort_protein_backbone_ranges, md_protein_backbone_data_t, md_amino_acid_atoms_t, ca)
SORT_BACKBONE_RANGES(sort_nucleic_backbone_ranges, md_nucleic_backbone_data_t, md_nucleic_acid_atoms_t, c4)

#undef SORT_BACKBONE_RANGES

bool md_backbone_infer(md_backbone_data_t* out, const md_system_t* sys, md_allocator_i* alloc) {
    MD_ASSERT(out);
    MD_ASSERT(sys);
    MD_ASSERT(alloc);
    MD_ASSERT(&out->protein != &out->protein && &out->nucleic != &out->nucleic);
    md_backbone_data_free(out, alloc);

    md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp_scope);

    static const size_t MIN_BACKBONE_LENGTH = 3;

    if (sys->bond.count && sys->atom.type_idx) {
        // Atoms: by the roles the bonds gave them, whatever the atoms and residues are called (POLYMER ATOMS)
        const size_t temp_pos = md_vm_arena_get_pos(temp_arena);
        backbone_chains_t prot = {0};
        if (protein_chains_from_roles(&prot, sys, temp_scope)) {
            const int32_t* unit_res = chain_units_residue(sys, &prot, temp_scope);
            commit_backbone_chains(out, sys, &prot, unit_res, MIN_BACKBONE_LENGTH, temp_scope, alloc);
            mark_backbone_chain_breaks(out, sys, &prot, unit_res, MIN_BACKBONE_LENGTH, temp_scope, alloc);
        }
        backbone_chains_t nucl = {0};
        if (nucleic_chains_from_roles(&nucl, sys, temp_scope)) {
            const int32_t* unit_res = chain_units_residue(sys, &nucl, temp_scope);
            commit_backbone_chains(out, sys, &nucl, unit_res, MIN_BACKBONE_LENGTH, temp_scope, alloc);
            mark_backbone_chain_breaks(out, sys, &nucl, unit_res, MIN_BACKBONE_LENGTH, temp_scope, alloc);
        }
        md_vm_arena_set_pos_back(temp_arena, temp_pos);
    }

    if (sys->part.count && sys->atom.type_idx) {
        // Beads: the traces of a coarse grained model, by the backbone beads its residues were given (the BB bead of
        // Martini, BB1 to BB3 of its nucleotides), residue by residue along each strand. The bonds of beads say nothing
        // of what a bead is.
        static const size_t MAX_BACKBONE_LENGTH = 1 << 15;
        size_t temp_pos = md_vm_arena_get_pos(temp_arena);

        {
            md_array(md_amino_acid_atoms_t) backbone_atoms = 0;
            md_array_ensure(backbone_atoms, MAX_BACKBONE_LENGTH, temp_arena);
            md_residue_idx_t res_base = -1;
            md_sequence_id_t prev_seq_id = -1;

            for (size_t part_idx = 0; part_idx < sys->part.count; ++part_idx) {
                if (md_system_part_kind(sys, part_idx) != MD_PART_KIND_PEPTIDE) {
                    continue;
                }

                md_urange_t range = md_system_part_residue_range(sys, part_idx);
                for (size_t res_idx = range.beg; res_idx < range.end; ++res_idx) {
                    const md_atom_idx_t bead = amino_acid_backbone_bead(sys, res_idx);
                    if (bead >= 0) {
                        md_amino_acid_atoms_t atoms;
                        md_memset(&atoms, -1, sizeof(atoms));
                        atoms.ca = bead;
                        md_sequence_id_t seq_id = md_system_residue_seq_id(sys, res_idx);

                        if (md_array_size(backbone_atoms) > 0) {
                            // The residue before it in the chain: a range ends where the chain breaks (CHAIN BREAKS in
                            // md_util.h), where the beads are not bonded, whatever the sequence numbers say
                            const md_amino_acid_atoms_t prev = *md_array_last(backbone_atoms);
                            if (!protein_trace_linked(sys, &prev, &atoms, seq_id == prev_seq_id + 1)) {
                                mark_chain_break(out, res_idx - 1, res_idx, alloc);
                                flush_protein_backbone_segment(backbone_atoms, &res_base, part_idx, &out->protein, alloc, MIN_BACKBONE_LENGTH);
                            }
                        }

                        if (res_base == -1) {
                            res_base = (md_residue_idx_t)res_idx;
                        }
                        md_array_push(backbone_atoms, atoms, temp_arena);
                        prev_seq_id = seq_id;
                        continue;
                    }

                    flush_protein_backbone_segment(backbone_atoms, &res_base, part_idx, &out->protein, alloc, MIN_BACKBONE_LENGTH);
                    prev_seq_id = -1;
                }
                flush_protein_backbone_segment(backbone_atoms, &res_base, part_idx, &out->protein, alloc, MIN_BACKBONE_LENGTH);
                prev_seq_id = -1;
            }
        }
        md_vm_arena_set_pos_back(temp_arena, temp_pos);

        {
            md_array(md_nucleic_acid_atoms_t) backbone_atoms = 0;
            md_array_ensure(backbone_atoms, MAX_BACKBONE_LENGTH, temp_arena);
            md_residue_idx_t res_base = -1;
            md_sequence_id_t prev_seq_id = -1;

            for (size_t part_idx = 0; part_idx < sys->part.count; ++part_idx) {
                if (!md_part_kind_is_nucleic_acid(md_system_part_kind(sys, part_idx))) {
                    continue;
                }

                md_urange_t range = md_system_part_residue_range(sys, part_idx);
                for (md_residue_idx_t res_idx = (int)range.beg; res_idx < (int)range.end; ++res_idx) {
                    md_nucleic_acid_atoms_t atoms;
                    if (nucleotide_backbone_beads(&atoms, sys, (size_t)res_idx)) {
                        md_sequence_id_t seq_id = md_system_residue_seq_id(sys, res_idx);

                        if (md_array_size(backbone_atoms) > 0) {
                            // The residue before it in the chain: a range ends where the chain breaks (CHAIN BREAKS in
                            // md_util.h); a residue without its phosphate bead, which is linked to nothing before it,
                            // begins one
                            const md_nucleic_acid_atoms_t prev = *md_array_last(backbone_atoms);
                            if (!nucleic_beads_linked(sys, &prev, &atoms, seq_id == prev_seq_id + 1)) {
                                if (!nucleic_chain_ends_between(sys, &prev, &atoms)) {
                                    mark_chain_break(out, (size_t)res_idx - 1, (size_t)res_idx, alloc);
                                }
                                flush_nucleic_backbone_segment(backbone_atoms, &res_base, part_idx, &out->nucleic, alloc, MIN_BACKBONE_LENGTH);
                            }
                        }

                        if (res_base == -1) {
                            res_base = res_idx;
                        }
                        md_array_push(backbone_atoms, atoms, temp_arena);
                        prev_seq_id = seq_id;
                        continue;
                    }

                    flush_nucleic_backbone_segment(backbone_atoms, &res_base, part_idx, &out->nucleic, alloc, MIN_BACKBONE_LENGTH);
                    prev_seq_id = -1;
                }
                flush_nucleic_backbone_segment(backbone_atoms, &res_base, part_idx, &out->nucleic, alloc, MIN_BACKBONE_LENGTH);
                prev_seq_id = -1;
            }
        }
        md_vm_arena_set_pos_back(temp_arena, temp_pos);
    }

    sort_protein_backbone_ranges(&out->protein, temp_scope);
    sort_nucleic_backbone_ranges(&out->nucleic, temp_scope);

    out->protein.range.count = md_array_size(out->protein.range.offset);
    if (out->protein.range.count) {
        md_array_push(out->protein.range.offset, (uint32_t)md_array_size(out->protein.segment.atoms), alloc);
    }
    out->protein.segment.count = md_array_size(out->protein.segment.atoms);

    out->nucleic.range.count = md_array_size(out->nucleic.range.offset);
    if (out->nucleic.range.count) {
        md_array_push(out->nucleic.range.offset, (uint32_t)md_array_size(out->nucleic.segment.atoms), alloc);
    }
    out->nucleic.segment.count = md_array_size(out->nucleic.segment.atoms);

    // The breaks, each once and in order of the residues
    const size_t num_breaks = md_array_size(out->breaks);
    if (num_breaks > 1) {
        qsort(out->breaks, num_breaks, sizeof(md_residue_idx_t), cmp_residue_idx);
        size_t w = 1;
        for (size_t k = 1; k < num_breaks; ++k) {
            if (out->breaks[k] != out->breaks[w - 1]) out->breaks[w++] = out->breaks[k];
        }
        md_array_shrink(out->breaks, w);
    }
    out->num_breaks = md_array_size(out->breaks);

    md_temp_end(temp_scope);
    return true;
}

void md_backbone_data_free(md_backbone_data_t* bb, md_allocator_i* alloc) {
    MD_ASSERT(bb);
    MD_ASSERT(alloc);
    md_array_free(bb->protein.range.offset, alloc);
    md_array_free(bb->protein.range.part_idx, alloc);
    md_array_free(bb->protein.segment.atoms, alloc);
    md_array_free(bb->protein.segment.res_idx, alloc);
    md_array_free(bb->nucleic.range.offset, alloc);
    md_array_free(bb->nucleic.range.part_idx, alloc);
    md_array_free(bb->nucleic.segment.atoms, alloc);
    md_array_free(bb->nucleic.segment.res_idx, alloc);
    md_array_free(bb->breaks, alloc);
    md_memset(bb, 0, sizeof(*bb));
}

void md_system_set_chain_breaks(md_system_t* sys, const md_residue_idx_t breaks[], size_t num_breaks) {
    MD_ASSERT(sys);
    if (sys->residue.count == 0) return;
    if (md_array_size(sys->residue.flags) < sys->residue.count) {
        const size_t size = md_array_size(sys->residue.flags);
        md_array_resize(sys->residue.flags, sys->residue.count, sys->alloc);
        md_memset(sys->residue.flags + size, 0, (sys->residue.count - size) * sizeof(md_residue_flags_t));
    }
    for (size_t i = 0; i < sys->residue.count; ++i) {
        sys->residue.flags[i] = (md_residue_flags_t)((uint32_t)sys->residue.flags[i] & ~(uint32_t)(MD_RESIDUE_FLAG_BREAK_BEFORE | MD_RESIDUE_FLAG_BREAK_AFTER));
    }
    for (size_t k = 0; k < num_breaks; ++k) {
        const size_t a = (size_t)breaks[k];
        const size_t b = a + 1;
        if (breaks[k] < 0 || b >= sys->residue.count) continue;
        // A break is no end: the residues on either side lose the terminal flags their classification gave them
        sys->residue.flags[a] = (md_residue_flags_t)(((uint32_t)sys->residue.flags[a] | MD_RESIDUE_FLAG_BREAK_AFTER)  & ~(uint32_t)MD_RESIDUE_FLAG_TERMINAL_END);
        sys->residue.flags[b] = (md_residue_flags_t)(((uint32_t)sys->residue.flags[b] | MD_RESIDUE_FLAG_BREAK_BEFORE) & ~(uint32_t)MD_RESIDUE_FLAG_TERMINAL_BEG);
        if (sys->atom.flags) {
            const md_urange_t range_a = md_system_residue_atom_range(sys, a);
            const md_urange_t range_b = md_system_residue_atom_range(sys, b);
            for (uint32_t i = range_a.beg; i < range_a.end; ++i) {
                sys->atom.flags[i] = (md_atom_flags_t)((uint32_t)sys->atom.flags[i] & ~(uint32_t)MD_ATOM_FLAG_TERMINAL_END);
            }
            for (uint32_t i = range_b.beg; i < range_b.end; ++i) {
                sys->atom.flags[i] = (md_atom_flags_t)((uint32_t)sys->atom.flags[i] & ~(uint32_t)MD_ATOM_FLAG_TERMINAL_BEG);
            }
        }
    }
    md_system_topology_changed(sys);
}
