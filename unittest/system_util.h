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

// What md_system_infer gives the system of its molecules (MD_SYSTEM_INFER_MOLECULE_BIT)
static inline bool test_infer_molecules(md_system_t* sys) {
    md_molecule_data_t molecules = {0};
    const bool ok = md_molecule_infer(&molecules, sys, sys->alloc);
    md_molecule_data_free(&sys->molecule, sys->alloc);
    sys->molecule = molecules;
    return ok;
}

static inline bool test_infer_chains(md_system_t* sys) {
    md_chain_data_t chains = {0};
    const bool ok = md_chain_infer(&chains, sys, sys->alloc);
    md_chain_data_free(&sys->chain, sys->alloc);
    sys->chain = chains;
    return ok;
}

// What a structure file's loader gave before it gave only the core (md_gro, md_pdb, md_mmcif): the bonds, what the
// residues are by them, and the chains. The tests below were written against that.
static inline bool test_infer_as_loaded(md_system_t* sys, md_system_state_t* state) {
    return md_system_infer(sys, state, MD_SYSTEM_INFER_BOND_BIT | MD_SYSTEM_INFER_RESIDUE_BIT, NULL) && test_infer_chains(sys);
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
    for (int c = 0; c < num_chains; ++c) {
        for (int i = 0; i < chain_len - 1; ++i) {
            md_atom_pair_t pair = { c * chain_len + i, c * chain_len + i + 1 };
            md_array_push(sys->bond.pairs, pair, alloc);
            sys->bond.count += 1;
        }
    }
    md_bond_build_connectivity(&sys->bond, (size_t)n, alloc);
}

// Appends a residue of coarse grained beads. The BB bead of an amino acid residue is its backbone.
static inline void cg_add_residue(md_system_t* sys, const char* res_name, int seq_id, md_residue_kind_t res_kind, const char* atom_names[], size_t num_atoms, md_allocator_i* alloc) {
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
        md_array_push(sys->atom.flags, flags, alloc);
        sys->atom.count += 1;
    }
    md_array_push(sys->residue.name, md_label_make(md_str_from_cstr(res_name)), alloc);
    md_array_push(sys->residue.seq_id, seq_id, alloc);
    md_array_push(sys->residue.flags, MD_RESIDUE_FLAG_NONE, alloc);
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

// The systems several tests check counts of, loaded and inferred as an application loads them
typedef struct reference_systems_t {
    md_allocator_i* alloc;
    md_system_t mol_ala;            // 1ALA-560ns.pdb, its frames published as run/ala
    md_system_t mol_pftaa;          // pftaa.gro
    md_system_t mol_nucleotides;    // nucleotides.gro
    md_system_t mol_centered;       // centered.gro
} reference_systems_t;

static inline reference_systems_t reference_systems_load(void) {
    reference_systems_t ref = { .alloc = md_vm_arena_create(MD_GIGABYTES(1)) };
    md_allocator_i* alloc = ref.alloc;

    ref.mol_ala.alloc = alloc;
    md_system_state_t mol_ala_state = { .alloc = alloc };
    md_pdb_system_init_from_file(&ref.mol_ala, &mol_ala_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), NULL);
    md_pdb_system_publish_run(&ref.mol_ala, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), MD_STR_LIT("run/ala"), MD_RUN_FLAG_DISABLE_CACHE_WRITE);
    md_system_infer(&ref.mol_ala, &mol_ala_state, MD_SYSTEM_INFER_ALL, NULL);

    ref.mol_pftaa.alloc = alloc;
    md_system_state_t mol_pftaa_state = { .alloc = alloc };
    md_gro_system_init_from_file(&ref.mol_pftaa, &mol_pftaa_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/pftaa.gro"), NULL);
    md_system_infer(&ref.mol_pftaa, &mol_pftaa_state, MD_SYSTEM_INFER_ALL, NULL);

    ref.mol_nucleotides.alloc = alloc;
    md_system_state_t mol_nucleotides_state = { .alloc = alloc };
    md_gro_system_init_from_file(&ref.mol_nucleotides, &mol_nucleotides_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/nucleotides.gro"), NULL);
    md_system_infer(&ref.mol_nucleotides, &mol_nucleotides_state, MD_SYSTEM_INFER_ALL, NULL);

    ref.mol_centered.alloc = alloc;
    md_system_state_t mol_centered_state = { .alloc = alloc };
    md_gro_system_init_from_file(&ref.mol_centered, &mol_centered_state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro"), NULL);
    md_system_infer(&ref.mol_centered, &mol_centered_state, MD_SYSTEM_INFER_ALL, NULL);

    return ref;
}

static inline void reference_systems_free(reference_systems_t* ref) {
    md_vm_arena_destroy(ref->alloc);
}
