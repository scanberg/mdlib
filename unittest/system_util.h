#pragma once

// What the tests of systems share: inferring a part of a system on its own, building small systems by hand, and the
// reference systems several of them check the counts of.

#include "utest.h"

#include <md_system.h>
#include <md_infer.h>
#include <md_pdb.h>
#include <md_gro.h>
#include <md_run.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_arena_allocator.h>
#include <core/md_vec_math.h>

#include <string.h>

// What md_system_infer gives the system of its molecules
static inline bool test_infer_molecules(md_system_t* sys) {
    md_molecule_data_t molecules = {0};
    // Of the derived: in the system's derived arena, which they are freed with
    const bool ok = md_molecule_infer(&molecules, sys, md_system_derived_alloc(sys));
    sys->molecule = molecules;
    return ok;
}

static inline bool test_infer_chains(md_system_t* sys) {
    md_chain_data_t chains = {0};
    const bool ok = md_chain_infer(&chains, sys, md_system_derived_alloc(sys));
    sys->chain = chains;
    return ok;
}

// What a structure file's loader gave before it gave only the core (md_gro, md_pdb, md_mmcif): the bonds, what the
// residues are by them, and the chains. The tests below were written against that.
static inline bool test_infer_as_loaded(md_system_t* sys, md_system_state_t* state) {
    return md_system_infer(sys, state, NULL) && test_infer_chains(sys);
}

// The tests lay out coordinates per axis, which reads well; a state holds them packed.
static inline md_float3_t* pack_xyz(md_float3_t* dst, const float* x, const float* y, const float* z, size_t n) {
    for (size_t i = 0; i < n; ++i) dst[i] = md_float3_set(x[i], y[i], z[i]);
    return dst;
}

// Builds a system of 'num_chains' disconnected linear chains of 'chain_len' atoms each.
// Atom c*chain_len + i is bonded to its two neighbours within the chain and nothing else.
static inline void make_chain_system(md_system_t* sys, md_allocator_i* alloc, int num_chains, int chain_len) {
    const int n = num_chains * chain_len;
    sys->alloc = alloc;
    sys->atom.count = (size_t)n;
    // The bonds are of the derived, in its arena (md_system_t, ALLOCATORS), which md_system_free gives back
    md_allocator_i* dalloc = md_system_derived_alloc(sys);
    for (int c = 0; c < num_chains; ++c) {
        for (int i = 0; i < chain_len - 1; ++i) {
            md_atom_pair_t pair = { c * chain_len + i, c * chain_len + i + 1 };
            md_array_push(sys->bond.pairs, pair, dalloc);
            sys->bond.count += 1;
        }
    }
    md_bond_build_connectivity(&sys->bond, (size_t)n, dalloc);
}

// Appends a residue of coarse grained beads. The BB bead of an amino acid residue is its backbone.
static inline void cg_add_residue(md_system_t* sys, const char* res_name, int seq_id, md_residue_kind_t res_kind, const char* atom_names[], size_t num_atoms, md_allocator_i* alloc) {
    // The flags are of the derived, in its arena where the system has an allocator (md_system_t, ALLOCATORS)
    md_allocator_i* dalloc = sys->alloc ? md_system_derived_alloc(sys) : alloc;
    if (sys->atom.type.count == 0) {
        // Type 0 is the "not found" sentinel
        md_particle_type_find_or_add(&sys->atom.type, MD_STR_LIT("?"), 0, 0, 0, 0, 0, alloc);
    }
    if (sys->residue.count == 0) {
        md_array_push(sys->residue.atom_offset, 0, alloc);
    }
    for (size_t i = 0; i < num_atoms; ++i) {
        const md_particle_type_flags_t type_flags = md_particle_type_flags_set_particle_kind(MD_PARTICLE_TYPE_FLAG_NONE, MD_PARTICLE_BEAD);
        const md_atom_flags_t flags = (res_kind == MD_RESIDUE_KIND_AMINO_ACID && strcmp(atom_names[i], "BB") == 0) ? MD_ATOM_FLAG_BACKBONE : MD_ATOM_FLAG_NONE;
        md_particle_type_idx_t type = md_particle_type_find_or_add(&sys->atom.type, md_str_from_cstr(atom_names[i]), 0, 50.0f, 2.35f, 0xFFFFFFFF, type_flags, alloc);
        md_array_push(sys->atom.type_idx, type, alloc);
        md_array_push(sys->atom.flags, flags, dalloc);
        sys->atom.count += 1;
    }
    md_array_push(sys->residue.name, md_label_make(md_str_from_cstr(res_name)), alloc);
    md_array_push(sys->residue.seq_id, seq_id, alloc);
    md_array_push(sys->residue.flags, MD_RESIDUE_FLAG_NONE, dalloc);
    md_array_push(sys->residue.atom_offset, (uint32_t)sys->atom.count, alloc);
    sys->residue.count += 1;
    md_system_set_residue_kind(sys, sys->residue.count - 1, res_kind);
}

static inline void check_molecule_invariants(int* utest_result, const md_system_t* sys) {
    size_t total = 0;
    for (size_t s = 0; s < md_system_molecule_count(sys); ++s) {
        const md_molecule_t mol = md_system_molecule(sys, s);
        total += mol.count;
        EXPECT_EQ(mol.parent_idx[0], mol.atom_idx[0]);
        for (size_t k = 0; k < mol.count; ++k) {
            const int32_t atom   = mol.atom_idx[k];
            const int32_t parent = mol.parent_idx[k];
            EXPECT_EQ(md_system_molecule_atom_slot(sys, atom), (int32_t)(sys->molecule.offset[s] + k));
            if (k > 0) {
                EXPECT_NE(parent, atom);
                EXPECT_LT(md_system_molecule_atom_slot(sys, parent), md_system_molecule_atom_slot(sys, atom));
            }
        }
    }
    EXPECT_EQ(total, sys->atom.count);
}

// ### STRUCTURES ###
//
// Whether molecule m is molecule r slot by slot, as the structures of the molecules (md_molecule_data_t) are told,
// checked here on their own terms: the same particle types, the same residues relative to the first slot's, each atom
// with as many bonds, to the same slots within the molecule or to atoms of the same types out of it, of the same kinds,
// origins and orders (but those perceived). A system may lack connectivity (no bonds), flags (pairs alone), types and residues.
static inline bool structure_test_res(const md_system_t* sys, int32_t a, int32_t a0, int64_t* out) {
    if (!sys->residue.count || !sys->residue.atom_offset) { *out = 0; return true; }
    const int r = md_system_residue_find_by_atom_idx(sys, (size_t)a), r0 = md_system_residue_find_by_atom_idx(sys, (size_t)a0);
    *out = r < 0 ? INT64_MIN : (int64_t)r - r0;
    return true;
}

static inline bool structure_test_same(const md_system_t* sys, size_t m, size_t r) {
    const md_molecule_data_t* mol = &sys->molecule;
    const uint32_t mb = mol->offset[m], me = mol->offset[m + 1], rb = mol->offset[r], re = mol->offset[r + 1];
    if (me - mb != re - rb) return false;
    const md_bond_data_t* bd = &sys->bond;
    const int32_t m0 = mol->atom_idx[mb], r0 = mol->atom_idx[rb];
    for (uint32_t k = 0; k < me - mb; ++k) {
        const int32_t a = mol->atom_idx[mb + k], b = mol->atom_idx[rb + k];
        if (sys->atom.type_idx && sys->atom.type_idx[a] != sys->atom.type_idx[b]) return false;
        int64_t ra, rr;
        structure_test_res(sys, a, m0, &ra);
        structure_test_res(sys, b, r0, &rr);
        if (ra != rr) return false;
        const uint32_t na = bd->conn.offset ? bd->conn.offset[a + 1] - bd->conn.offset[a] : 0;
        const uint32_t nb = bd->conn.offset ? bd->conn.offset[b + 1] - bd->conn.offset[b] : 0;
        if (na != nb) return false;
        for (uint32_t j = 0; j < na; ++j) {
            const int32_t oa = bd->conn.atom_idx[bd->conn.offset[a] + j], ob = bd->conn.atom_idx[bd->conn.offset[b] + j];
            const uint32_t sa = (uint32_t)mol->atom_slot[oa], sb = (uint32_t)mol->atom_slot[ob];
            const bool ia = sa >= mb && sa < me, ib = sb >= rb && sb < re;
            if (ia != ib) return false;
            if (ia && sa - mb != sb - rb) return false;
            if (!ia && sys->atom.type_idx && sys->atom.type_idx[oa] != sys->atom.type_idx[ob]) return false;
            if (bd->flags) {
                const md_bond_flags_t fa = bd->flags[bd->conn.bond_idx[bd->conn.offset[a] + j]];
                const md_bond_flags_t fb = bd->flags[bd->conn.bond_idx[bd->conn.offset[b] + j]];
                if (md_bond_kind(fa) != md_bond_kind(fb) || md_bond_origin(fa) != md_bond_origin(fb)) return false;
                // An order perceived is the chemistry's, of the structure's atoms, not part of it
                const bool perceived = (fa & MD_BOND_FLAG_ORDER_PERCEIVED) || (fb & MD_BOND_FLAG_ORDER_PERCEIVED);
                if (!perceived && md_bond_order(fa) != md_bond_order(fb)) return false;
            }
        }
    }
    return true;
}

// The structures of the molecules: each its reference's own, every copy its reference slot by slot, no two references
// alike (up to max_refs of them compared pairwise), and the structures valid for the topology as it is
static inline void check_structure_invariants(int* utest_result, md_system_t* sys, size_t max_refs) {
    const md_molecule_data_t* mol = &sys->molecule;
    if (mol->count == 0) return;
    ASSERT_TRUE(mol->structure_idx != NULL);
    ASSERT_GT(mol->structure_count, (size_t)0);
    ASSERT_LE(mol->structure_count, mol->count);
    for (size_t s = 0; s < mol->structure_count; ++s) {
        const int32_t r = mol->structure_ref[s];
        ASSERT_TRUE(r >= 0 && (size_t)r < mol->count);
        EXPECT_EQ(mol->structure_idx[r], (int32_t)s);
    }
    for (size_t m = 0; m < mol->count; ++m) {
        const int32_t s = mol->structure_idx[m];
        ASSERT_TRUE(s >= 0 && (size_t)s < mol->structure_count);
        const int32_t r = mol->structure_ref[s];
        if (r != (int32_t)m) EXPECT_TRUE_MSG(structure_test_same(sys, m, (size_t)r), "a copy is not its reference");
    }
    if (mol->structure_count <= max_refs) {
        for (size_t s = 0; s < mol->structure_count; ++s) {
            for (size_t t = s + 1; t < mol->structure_count; ++t) {
                EXPECT_FALSE_MSG(structure_test_same(sys, (size_t)mol->structure_ref[s], (size_t)mol->structure_ref[t]), "two structures alike");
            }
        }
    }
    EXPECT_TRUE(md_system_molecule_structures_valid(sys));
}

// The systems several tests check counts of, loaded and inferred as an application loads them
typedef struct reference_systems_t {
    md_allocator_i* alloc;
    md_system_t mol_ala;            // 1ALA-560ns.pdb, its frames published as run/ala
    md_system_t mol_pftaa;          // pftaa.gro
    md_system_t mol_nucleotides;    // nucleotides.gro
    md_system_t mol_centered;       // centered.gro
} reference_systems_t;

// Into *ref: four systems and their derived data, which there is no reason to copy out of a function
static inline void reference_systems_load(reference_systems_t* ref) {
    *ref = (reference_systems_t){ .alloc = md_vm_arena_create(MD_GIGABYTES(1)) };
    md_allocator_i* alloc = ref->alloc;

    ref->mol_ala.alloc = alloc;
    md_system_state_t mol_ala_state = { .alloc = alloc };
    md_pdb_system_init_from_file(&ref->mol_ala, &mol_ala_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), NULL);
    md_pdb_system_publish_run(&ref->mol_ala, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), MD_STR_LIT("run/ala"), MD_RUN_FLAG_DISABLE_CACHE_WRITE);
    md_system_infer(&ref->mol_ala, &mol_ala_state, NULL);

    ref->mol_pftaa.alloc = alloc;
    md_system_state_t mol_pftaa_state = { .alloc = alloc };
    md_gro_system_init_from_file(&ref->mol_pftaa, &mol_pftaa_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/pftaa.gro"), NULL);
    md_system_infer(&ref->mol_pftaa, &mol_pftaa_state, NULL);

    ref->mol_nucleotides.alloc = alloc;
    md_system_state_t mol_nucleotides_state = { .alloc = alloc };
    md_gro_system_init_from_file(&ref->mol_nucleotides, &mol_nucleotides_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/nucleotides.gro"), NULL);
    md_system_infer(&ref->mol_nucleotides, &mol_nucleotides_state, NULL);

    ref->mol_centered.alloc = alloc;
    md_system_state_t mol_centered_state = { .alloc = alloc };
    md_gro_system_init_from_file(&ref->mol_centered, &mol_centered_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro"), NULL);
    md_system_infer(&ref->mol_centered, &mol_centered_state, NULL);
}

static inline void reference_systems_free(reference_systems_t* ref) {
    md_vm_arena_destroy(ref->alloc);
}
