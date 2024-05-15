#include "utest.h"
#include <string.h>

#include <md_pdb.h>
#include <md_gro.h>
#include <md_xyz.h>
#include <md_trajectory.h>
#include <md_molecule.h>
#include <md_util.h>
#include <md_smiles.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_str_builder.h>
#include <core/md_bitfield.h>

#include "rmsd.h"

struct util {
    md_allocator_i* alloc;
    md_molecule_t mol_ala;
    md_molecule_t mol_pftaa;
    md_molecule_t mol_nucleotides;
    md_molecule_t mol_centered;
    md_molecule_t mol_dna;
    md_molecule_t mol_trp;
    md_molecule_t mol_aspirine;

    md_trajectory_i* traj_ala;
};

UTEST_F_SETUP(util) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    utest_fixture->alloc = alloc;

    md_pdb_molecule_api()->init_from_file(&utest_fixture->mol_ala, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), NULL, alloc);
    md_util_molecule_postprocess(&utest_fixture->mol_ala, alloc, MD_UTIL_POSTPROCESS_ALL);

    md_gro_molecule_api()->init_from_file(&utest_fixture->mol_pftaa, STR_LIT(MD_UNITTEST_DATA_DIR "/pftaa.gro"), NULL, alloc);
    md_util_molecule_postprocess(&utest_fixture->mol_pftaa, alloc, MD_UTIL_POSTPROCESS_ALL);

    md_gro_molecule_api()->init_from_file(&utest_fixture->mol_nucleotides, STR_LIT(MD_UNITTEST_DATA_DIR "/nucleotides.gro"), NULL, alloc);
    md_util_molecule_postprocess(&utest_fixture->mol_nucleotides, alloc, MD_UTIL_POSTPROCESS_ALL);

    md_gro_molecule_api()->init_from_file(&utest_fixture->mol_centered, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro"), NULL, alloc);
    md_util_molecule_postprocess(&utest_fixture->mol_centered, alloc, MD_UTIL_POSTPROCESS_ALL);

    md_gro_molecule_api()->init_from_file(&utest_fixture->mol_dna, STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro"), NULL, alloc);
    md_util_molecule_postprocess(&utest_fixture->mol_dna, alloc, MD_UTIL_POSTPROCESS_ALL);

    md_gro_molecule_api()->init_from_file(&utest_fixture->mol_trp, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.gro"), NULL, alloc);
    md_util_molecule_postprocess(&utest_fixture->mol_trp, alloc, MD_UTIL_POSTPROCESS_ALL);

    md_gro_molecule_api()->init_from_file(&utest_fixture->mol_aspirine, STR_LIT(MD_UNITTEST_DATA_DIR "/inside-md-pullout.gro"), NULL, alloc);
    md_util_molecule_postprocess(&utest_fixture->mol_aspirine, alloc, MD_UTIL_POSTPROCESS_ALL);

    utest_fixture->traj_ala = md_pdb_trajectory_create(STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), alloc, MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE);
}

UTEST_F_TEARDOWN(util) {
    md_vm_arena_destroy(utest_fixture->alloc);
}

UTEST_F(util, bonds) {
    EXPECT_EQ(152,    utest_fixture->mol_ala.bond.count);
    EXPECT_EQ(55,     utest_fixture->mol_pftaa.bond.count);
    EXPECT_EQ(40,     utest_fixture->mol_nucleotides.bond.count);
    EXPECT_EQ(163504, utest_fixture->mol_centered.bond.count);
}

UTEST_F(util, chain) {
    EXPECT_EQ(1,   utest_fixture->mol_ala.chain.count);
	EXPECT_EQ(0,   utest_fixture->mol_pftaa.chain.count);
	EXPECT_EQ(0,   utest_fixture->mol_nucleotides.chain.count);
	EXPECT_EQ(253, utest_fixture->mol_centered.chain.count);

    md_molecule_t pdb_mol;
    md_pdb_molecule_api()->init_from_file(&pdb_mol, STR_LIT(MD_UNITTEST_DATA_DIR "/dppc64_chain_test.pdb"), NULL, utest_fixture->alloc);
    md_util_molecule_postprocess(&pdb_mol, utest_fixture->alloc, MD_UTIL_POSTPROCESS_ALL);

    const md_molecule_t* mol = &utest_fixture->mol_centered;
    if (mol->chain.count == 0) return;

    size_t ref_size = md_chain_atom_count(mol->chain, 0);
    for (size_t i = 1; i < mol->chain.count; ++i) {
        size_t size = md_chain_atom_count(mol->chain, i);
        EXPECT_EQ(ref_size, size);
    }
}

UTEST_F(util, backbone) {
	EXPECT_EQ(0,        utest_fixture->mol_pftaa.protein_backbone.range.count);
	EXPECT_EQ(0,        utest_fixture->mol_nucleotides.protein_backbone.range.count);
    EXPECT_EQ(1,        utest_fixture->mol_ala.protein_backbone.range.count);
    EXPECT_EQ(15,       utest_fixture->mol_ala.protein_backbone.count);
	EXPECT_EQ(253,      utest_fixture->mol_centered.protein_backbone.range.count);
    EXPECT_EQ(10626,    utest_fixture->mol_centered.protein_backbone.count); // Should be equal to the total count of residues in chains
}

UTEST_F(util, structure) {
    size_t num_structures_pftaa = md_index_data_count(utest_fixture->mol_pftaa.structure);
    EXPECT_EQ(1, num_structures_pftaa);
    size_t num_structures_nucleotides = md_index_data_count(utest_fixture->mol_nucleotides.structure);
    EXPECT_EQ(2, num_structures_nucleotides);
    size_t num_structures_ala = md_index_data_count(utest_fixture->mol_ala.structure);
	EXPECT_EQ(1, num_structures_ala);
	size_t num_structures_centered = md_index_data_count(utest_fixture->mol_centered.structure);
	EXPECT_EQ(253+61, num_structures_centered);
}

UTEST_F(util, rmsd) {
    md_allocator_i* alloc = utest_fixture->alloc;
    md_molecule_t* mol = &utest_fixture->mol_ala;
    md_trajectory_i* traj = utest_fixture->traj_ala;
    ASSERT_TRUE(mol);
    ASSERT_TRUE(traj);

    const int64_t stride = mol->atom.count;
    const int64_t mem_size = stride * 6 * sizeof(float) + stride * 6 * sizeof(double);
    float* mem = md_alloc(alloc, mem_size);
    float* x[2] = {
        mem + stride * 0,
        mem + stride * 1,
    };
    float* y[2] = {
        mem + stride * 2,
        mem + stride * 3,
    };
    float* z[2] = {
        mem + stride * 4,
        mem + stride * 5,
    };

    const float* w = mol->atom.mass;

    double* xyz0 = (double*)(mem + stride * 6);
    double* xyz1 = (double*)(mem + stride * 6) + stride * 3;

    md_trajectory_load_frame(traj, 0, NULL, x[0], y[0], z[0]);
    md_trajectory_load_frame(traj, 1, NULL, x[1], y[1], z[1]);

    for (int64_t i = 0; i < mol->atom.count; ++i) {
        xyz0[i * 3 + 0] = x[0][i];
        xyz0[i * 3 + 1] = y[0][i];
        xyz0[i * 3 + 2] = z[0][i];

        xyz1[i * 3 + 0] = x[1][i];
        xyz1[i * 3 + 1] = y[1][i];
        xyz1[i * 3 + 2] = z[1][i];
    }

    // Reference
    double ref_rmsd;
    fast_rmsd((double(*)[3])xyz0, (double(*)[3])xyz1, (int)mol->atom.count, &ref_rmsd);

    // Our implementation
    const float* const cx[2] = { x[0], x[1] };
    const float* const cy[2] = { y[0], y[1] };
    const float* const cz[2] = { z[0], z[1] };

    vec3_t com[2] = {
        md_util_com_compute(x[0], y[0], z[0], w, 0, mol->atom.count, 0),
        md_util_com_compute(x[1], y[1], z[1], w, 0, mol->atom.count, 0),
    };
    double rmsd = md_util_rmsd_compute(cx, cy, cz, w, 0, mol->atom.count, com);
    
    EXPECT_LE(fabs(ref_rmsd - rmsd), 0.1);

    md_free(alloc, mem, mem_size);
}

UTEST(util, com) {
    /* DISCLAIMER
        Computing the center of mass is a bit tricky, because of the periodic boundary conditions.
        The problem occurs when we have structures which have an extent which covers more than half the box.
        There are many different variations of how to handle this and it seems none of them are perfect.
        Unless you handpick your algorithm based on some external knowledge of the structure you are dealing with.

        In this case, we have settled on the trigonometric approach presented in:
        Bai, Linge, and David Breen. "Calculating center of mass in an unbounded 2D environment." Journal of Graphics Tools 13.4 (2008): 53-60.
    */
    vec3_t pbc_ext = {5,0,0};
    md_unit_cell_t unit_cell = md_util_unit_cell_from_extent(5,0,0);
    {
        const vec4_t xyzw[] = {
            {1,0,0,1},
            {2,0,0,1},
            {3,0,0,1},
            {4,0,0,1},
        };

        vec3_t com = md_util_com_compute_vec4(xyzw, 0, ARRAY_SIZE(xyzw), &unit_cell);
        EXPECT_NEAR(2.5f, com.x, 1.0E-5F);
        EXPECT_EQ(0, com.y);
        EXPECT_EQ(0, com.z);
    }
    
    {
        const vec4_t xyzw[] = {
            {0,0,0,1},
            {5,0,0,1},
        };

        vec3_t com = md_util_com_compute_vec4(xyzw, 0, ARRAY_SIZE(xyzw), &unit_cell);
		com = vec3_deperiodize(com, (vec3_t){ 0,0,0 }, pbc_ext);
        EXPECT_NEAR(0, com.x, 1.0E-5F);
        EXPECT_EQ(0, com.y);
        EXPECT_EQ(0, com.z);
    }

    {
        const vec4_t xyzw[] = {
            {0,0,0,1},
            {4,0,0,1},
        };

        /*
        Here we expect the 4 to wrap around to -1,
        then added to the 0, producing a center of mass of -0.5.
        which is then placed within the period to 4.5.
        */

        vec3_t com = md_util_com_compute_vec4(xyzw, 0, ARRAY_SIZE(xyzw), &unit_cell);
        com = vec3_deperiodize(com, vec3_mul_f(pbc_ext, 0.5f), pbc_ext);
        EXPECT_NEAR(4.5f, com.x, 1.0E-5F);
        EXPECT_EQ(0, com.y);
        EXPECT_EQ(0, com.z);
    }

    {
        const vec3_t pbc_ext = { 5,0,0 };

        const vec4_t pos0[] = {
            {4,0,0, 1},
            {5,0,0, 1},
            {6,0,0, 1},
            {7,0,0, 1},
        };

        const vec4_t pos1[] = {
            {4,0,0, 1},
            {0,0,0, 1},
            {1,0,0, 1},
            {2,0,0, 1},
        };

        const vec4_t pos2[] = {
            {-1,0, 0, 1},
            {0 ,0, 0, 1},
            {1 ,0, 0, 1},
            {2 ,0, 0, 1},
        };

        vec3_t com0 = md_util_com_compute_vec4(pos0, 0, ARRAY_SIZE(pos0), &unit_cell);
        vec3_t com1 = md_util_com_compute_vec4(pos1, 0, ARRAY_SIZE(pos1), &unit_cell);
        vec3_t com2 = md_util_com_compute_vec4(pos2, 0, ARRAY_SIZE(pos2), &unit_cell);

        com0 = vec3_deperiodize(com0, (vec3_t){ 0,0,0 }, pbc_ext);
        com1 = vec3_deperiodize(com1, (vec3_t){ 0,0,0 }, pbc_ext);
        com2 = vec3_deperiodize(com2, (vec3_t){ 0,0,0 }, pbc_ext);
        
        EXPECT_NEAR(0.5f, com0.x, 1.0E-5F);
        EXPECT_NEAR(0.5f, com1.x, 1.0E-5F);
        EXPECT_NEAR(0.5f, com2.x, 1.0E-5F);
    }
}

UTEST_F(util, structures) {
    size_t num_structures = 0;

    num_structures = md_index_data_count(utest_fixture->mol_nucleotides.structure);
    EXPECT_EQ(num_structures, 2);

    num_structures = md_index_data_count(utest_fixture->mol_ala.structure);
	EXPECT_EQ(num_structures, 1);

	num_structures = md_index_data_count(utest_fixture->mol_pftaa.structure);
	EXPECT_EQ(num_structures, 1);

	num_structures = md_index_data_count(utest_fixture->mol_centered.structure);
    EXPECT_EQ(num_structures, 253 + 61); // Chains + PFTAA
}

UTEST_F(util, rings_common) {
    int64_t num_rings = 0;

    num_rings = md_index_data_count(utest_fixture->mol_nucleotides.ring);
    EXPECT_EQ(num_rings, 4);
    
    num_rings = md_index_data_count(utest_fixture->mol_ala.ring);
    EXPECT_EQ(num_rings, 0);
   
    num_rings = md_index_data_count(utest_fixture->mol_pftaa.ring);
    EXPECT_EQ(num_rings, 5);

    num_rings = md_index_data_count(utest_fixture->mol_centered.ring);
    EXPECT_EQ(num_rings, 2076);
}

UTEST(util, rings_c60) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
	md_molecule_t mol = {0};
	md_pdb_molecule_api()->init_from_file(&mol, STR_LIT(MD_UNITTEST_DATA_DIR "/c60.pdb"), NULL, alloc);
	md_util_molecule_postprocess(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

	EXPECT_EQ(mol.atom.count, 60);
	EXPECT_EQ(mol.bond.count, 90);

    const size_t num_rings = md_index_data_count(mol.ring);
    EXPECT_EQ(num_rings, 32);

    const size_t num_structures = md_index_data_count(mol.structure);
    EXPECT_EQ(num_structures, 1);

    md_vm_arena_destroy(alloc);
}

UTEST(util, rings_c720) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_molecule_t mol = {0};
    md_xyz_molecule_api()->init_from_file(&mol, STR_LIT(MD_UNITTEST_DATA_DIR "/c720.xyz"), NULL, alloc);
    md_util_molecule_postprocess(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    EXPECT_EQ(mol.atom.count, 720);
    EXPECT_EQ(mol.bond.count, 1080);

    const size_t num_rings = md_index_data_count(mol.ring);
    EXPECT_EQ(num_rings, 362);

    const size_t num_structures = md_index_data_count(mol.structure);
    EXPECT_EQ(num_structures, 1);

    md_vm_arena_destroy(alloc);
}

UTEST(util, rings_14kr) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_molecule_t mol = {0};
    md_pdb_molecule_api()->init_from_file(&mol, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), NULL, alloc);
    md_util_molecule_postprocess(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    const size_t num_rings = md_index_data_count(mol.ring);
    EXPECT_EQ(num_rings, 207);

    md_vm_arena_destroy(alloc);
}

UTEST(util, rings_trytophan_pdb) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_molecule_t mol = {0};
    md_pdb_molecule_api()->init_from_file(&mol, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan.pdb"), NULL, alloc);
    md_util_molecule_postprocess(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    const size_t num_rings = md_index_data_count(mol.ring);
    EXPECT_EQ(num_rings, 2);

    const size_t num_structures = md_index_data_count(mol.structure);
    EXPECT_EQ(num_structures, 1);

    md_vm_arena_destroy(alloc);
}

UTEST(util, rings_trytophan_xyz) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_molecule_t mol = {0};
    md_xyz_molecule_api()->init_from_file(&mol, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan.xyz"), NULL, alloc);
    md_util_molecule_postprocess(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    const size_t num_rings = md_index_data_count(mol.ring);
    EXPECT_EQ(num_rings, 2);

    const size_t num_structures = md_index_data_count(mol.structure);
    EXPECT_EQ(num_structures, 1);

    md_vm_arena_destroy(alloc);
}

UTEST(util, rings_full) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_molecule_t mol = {0};
    md_xyz_molecule_api()->init_from_file(&mol, STR_LIT(MD_UNITTEST_DATA_DIR "/full.xyz"), NULL, alloc);
    md_util_molecule_postprocess(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    const size_t num_rings = md_index_data_count(mol.ring);
    EXPECT_EQ(num_rings, 195);

    const size_t num_structures = md_index_data_count(mol.structure);
    EXPECT_EQ(num_structures, 1);

    md_vm_arena_destroy(alloc);
}

UTEST(util, rings_ciprofloxacin) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_molecule_t mol = {0};
    md_pdb_molecule_api()->init_from_file(&mol, STR_LIT(MD_UNITTEST_DATA_DIR "/ciprofloxacin.pdb"), NULL, alloc);
    md_util_molecule_postprocess(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    const int64_t num_rings = md_index_data_count(mol.ring);
    ASSERT_EQ(num_rings, 4);
    EXPECT_EQ(md_index_range_size(mol.ring, 0), 3);
    EXPECT_EQ(md_index_range_size(mol.ring, 1), 6);
    EXPECT_EQ(md_index_range_size(mol.ring, 2), 6);
    EXPECT_EQ(md_index_range_size(mol.ring, 3), 6);

    md_vm_arena_destroy(alloc);
}

UTEST_F(util, structure_matching_amyloid_chain) {
    md_allocator_i* alloc = utest_fixture->alloc;
    md_molecule_t* mol = &utest_fixture->mol_centered;

#if 1
    {
        // Test for the chains
        const int ref_structure_idx = 0;
        int* ref_idx = md_index_range_beg(mol->structure, ref_structure_idx);
        size_t ref_len = md_index_range_size(mol->structure, ref_structure_idx);

        // Prune Hydrogen
        if (true) {
            md_array(int) new_idx = 0;
            for (size_t i = 0; i < ref_len; ++i) {
                if (mol->atom.element[i] != 1) {
                    md_array_push(new_idx, ref_idx[i], alloc);
                }
            }
            ref_idx = new_idx;
            ref_len = md_array_size(new_idx);
        }

        md_timestamp_t t0 = md_time_current();
        md_index_data_t result = md_util_match_by_element(ref_idx, ref_len, MD_UTIL_MATCH_MODE_FIRST, MD_UTIL_MATCH_LEVEL_CHAIN, mol, alloc);
        md_timestamp_t t1 = md_time_current();
        printf("time: %f ms\n", md_time_as_milliseconds(t1-t0));
        size_t result_count = md_index_data_count(result);
        printf("result count: %zu\n", result_count);
        EXPECT_EQ(result_count, 253);
    }
#endif
}

UTEST_F(util, structure_matching_PFTAA) {
    md_allocator_i* alloc = utest_fixture->alloc;
    md_molecule_t* mol = &utest_fixture->mol_centered;

#if 1
    {
        // Test for the PFTAAs
        const int ref_structure_idx = 253;
        const int* ref_idx = md_index_range_beg(mol->structure, ref_structure_idx);
        const size_t ref_size = md_index_range_size(mol->structure, ref_structure_idx);

        md_timestamp_t t0 = md_time_current();
        md_index_data_t result = md_util_match_by_element(ref_idx, ref_size, MD_UTIL_MATCH_MODE_FIRST, MD_UTIL_MATCH_LEVEL_RESIDUE, mol, alloc);
        md_timestamp_t t1 = md_time_current();
        printf("time: %f ms\n", md_time_as_milliseconds(t1-t0));
        size_t result_count = md_index_data_count(result);
        printf("result count: %zu\n", result_count);
        EXPECT_EQ(result_count, 61);
    }
#endif
}

UTEST_F(util, structure_matching_PFTAA_ring) {
    md_allocator_i* alloc = utest_fixture->alloc;
    md_molecule_t* mol = &utest_fixture->mol_pftaa;
    {
#if 1
        // Rings, represents a ring within the molecule
        const int ref_idx[] = {19,20,21,22,24};
        const size_t ref_size = ARRAY_SIZE(ref_idx);

        md_index_data_t result = md_util_match_by_element(ref_idx, ref_size, MD_UTIL_MATCH_MODE_UNIQUE, MD_UTIL_MATCH_LEVEL_STRUCTURE, mol, alloc);
        const size_t result_count = md_index_data_count(result);
        EXPECT_EQ(result_count, 5);
#endif
#if 0
        for (size_t i = 0; i < result_count; ++i) {
            printf("[");
            for (int* it = md_index_range_beg(result, i); it != md_index_range_end(result, i); ++it) {
				printf("%d", *it + 1);
                if (it < md_index_range_end(result, i) - 1) printf(", ");
			}
            printf("]\n");
        }

        md_index_data_free(&result, alloc);
#endif
    }
    {
#if 1
        // Symmetry, represents half molecule minus the ring which holds the two symmetric parts together.
        const int ref_idx[] = {0,1,2,3,4,5,6,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,46,47,48};
        const int64_t ref_len = ARRAY_SIZE(ref_idx);
        md_timestamp_t t0 = md_time_current();
        md_index_data_t result = md_util_match_by_element(ref_idx, ref_len, MD_UTIL_MATCH_MODE_UNIQUE, MD_UTIL_MATCH_LEVEL_STRUCTURE, mol, alloc);
        md_timestamp_t t1 = md_time_current();
        printf("time: %f ms\n", md_time_as_milliseconds(t1-t0));
        size_t result_count = md_index_data_count(result);
        EXPECT_EQ(result_count, 2);

        md_index_data_free(&result, alloc);
#endif
    }
}

UTEST_F(util, structure_matching_smiles) {
    md_allocator_i* alloc = utest_fixture->alloc;
    
    // There are natural variations in how the alanine molecule is structured depending on its environment. (terminal residue etc)
    //const char ALANINE1[] = "N[C@@H]([CH3])C(=O)"; 
    //const char ALANINE2[] = "N[C@@H]([CH3])C(-O)(-O)";

    const str_t ALANINE         = STR_LIT("[NH][C@@H]([CH3])CO");
    const str_t ARGININE        = STR_LIT("[NH][C@@H](CO)[CH2][CH2][CH2][NH]C([NH2])[NH2]");
    const str_t GLYCINE         = STR_LIT("[NH][CH2]CO");
    const str_t ASPARAGINE      = STR_LIT("[NH][C@@H](CO)[CH2]C(O)[NH2]");
    const str_t ASPARTATE       = STR_LIT("[NH][C@@H](CO)[CH2]C(O)O");
    const str_t CYSTEINE        = STR_LIT("[NH][C@@H](CO)[CH2][SH]");
    const str_t GLUTAMIC_ACID   = STR_LIT("[NH][C@@H](CO)[CH2][CH2]C(O)O");
    const str_t GLUTAMINE       = STR_LIT("[NH][C@@H](CO)[CH2][CH2]C(O)[NH2]");
    const str_t HISTIDINE       = STR_LIT("[NH][C@@H](CO)[CH2]C1:N:[CH]:[NH]:[CH]:1");
    const str_t ISOLEUCINE      = STR_LIT("[NH][C@@H](CO)[CH]([CH3])[CH2][CH3]");
    const str_t LEUCINE         = STR_LIT("[NH][C@@H](CO)[CH2][CH]([CH3])[CH3]");
    const str_t LYSINE          = STR_LIT("[NH][C@@H](CO)[CH2][CH2][CH2][CH2][NH3]");
    const str_t METHIONINE      = STR_LIT("[NH][C@@H](CO)[CH2][CH2]S[CH3]");
    const str_t PHENYLALANINE   = STR_LIT("[NH][C@@H](CO)[CH2]C1:[CH]:[CH]:[CH]:[CH]:[CH]:1");
    const str_t PROLINE         = STR_LIT("N1[C@@H](CO)[CH2][CH2][CH2]1");
    const str_t SERINE          = STR_LIT("[NH][C@@H](CO)[CH2][OH]");
    const str_t THREONINE       = STR_LIT("[NH][C@@H](CO)[CH]([CH3])[OH]");
    const str_t TRYPTOPHAN      = STR_LIT("[NH][C@@H](CO)[CH2]C1:[CH]:[NH]:C2:[CH]:[CH]:[CH]:[CH]:C:1:2");
    const str_t TYROSINE        = STR_LIT("[NH][C@@H](CO)[CH2]C1:[CH]:[CH]:C([OH]):[CH]:[CH]:1");
    const str_t VALINE          = STR_LIT("[NH][C@@H](CO)[CH]([CH3])[CH3]");
    
    const str_t SELENOCYSTEINE  = STR_LIT("[NH][C@@H](CO)[CH2][SeH]");
    const str_t PYRROLYSINE     = STR_LIT("[NH][C@@H](CO)[CH2][CH2][CH2][CH2][NH]C(=O)C1=N[CH][CH2][CH]1[CH3]");

    const str_t DA              = STR_LIT("P(O)(O)O[CH2][CH]1[CH](O)[CH2][CH](O1)N2CNC3C(N)NCNC23");
    const str_t DC              = STR_LIT("P(O)(O)O[CH2][CH]1[CH](O)[CH2][CH](O1)N2C(O)NC(N)CC2");
    const str_t DG              = STR_LIT("P(O)(O)O[CH2][CH]1[CH](O)[CH2][CH](O1)N2C3NC(N)NC(O)C3NC2");
    const str_t DT              = STR_LIT("P(O)(O)O[CH2][CH]1[CH](O)[CH2][CH](O1)N2C(O)NC(O)C(C2)C");
    const str_t DU              = STR_LIT("P(O)(O)O[CH2][CH]1[CH](O)[CH2][CH](O1)N2C(O)NC(O)CC2");

    typedef struct {
        str_t name;
        str_t smiles;
    } res_t;

    typedef struct {
        str_t name;
        const md_molecule_t* mol;
    } mol_t;

    res_t residues[] = {
        {STR_LIT("ALA"), ALANINE},
        {STR_LIT("ARG"), ARGININE},
        {STR_LIT("ASN"), ASPARAGINE},
        {STR_LIT("ASP"), ASPARTATE},
        {STR_LIT("CYS"), CYSTEINE},
        {STR_LIT("GLY"), GLYCINE},
        {STR_LIT("GLU"), GLUTAMIC_ACID},
        {STR_LIT("GLN"), GLUTAMINE},
        {STR_LIT("HIS"), HISTIDINE},
        {STR_LIT("ILE"), ISOLEUCINE},
        {STR_LIT("LEU"), LEUCINE},
        {STR_LIT("LYS"), LYSINE},
        {STR_LIT("MET"), METHIONINE},
        {STR_LIT("PHE"), PHENYLALANINE},
        {STR_LIT("PRO"), PROLINE},
        {STR_LIT("SER"), SERINE},
        {STR_LIT("THR"), THREONINE},
        {STR_LIT("TRP"), TRYPTOPHAN},
        {STR_LIT("TYR"), TYROSINE},
        {STR_LIT("VAL"), VALINE},

        {STR_LIT("SEC"), SELENOCYSTEINE},
        {STR_LIT("PYR"), PYRROLYSINE},

        {STR_LIT("DA"),  DA},
        {STR_LIT("DC"),  DC},
        {STR_LIT("DG"),  DG},
        {STR_LIT("DT"),  DT},
    };

    mol_t molecules[] = {
        {STR_LIT("ALA"), &utest_fixture->mol_ala},
        {STR_LIT("AMYLOID PFTAA"), &utest_fixture->mol_centered},
        {STR_LIT("NUCLEOTIDES"), &utest_fixture->mol_nucleotides},
        {STR_LIT("DNA"), &utest_fixture->mol_dna},
        {STR_LIT("TRP"), &utest_fixture->mol_trp},
        //{STR_LIT("ASPIRINE"), &utest_fixture->mol_aspirine},
    };

    md_timestamp_t t0 = md_time_current();
    for (const mol_t* mol = molecules; mol != molecules + ARRAY_SIZE(molecules); ++mol) {
        md_vm_arena_temp_t temp = md_vm_arena_temp_begin(alloc);
        for (const res_t* res = residues; res != residues + ARRAY_SIZE(residues); ++res) {
            
            md_bitfield_t bf_ref = md_bitfield_create(alloc);
            for (size_t i = 0; i < mol->mol->residue.count; ++i) {
                if (str_eq(LBL_TO_STR(mol->mol->residue.name[i]), res->name)) {
                    md_bitfield_set_bit(&bf_ref, i);
                }
            }
            size_t resname_count = md_bitfield_popcount(&bf_ref);

            md_index_data_t match_result = md_util_match_smiles(res->smiles, MD_UTIL_MATCH_MODE_FIRST, MD_UTIL_MATCH_LEVEL_RESIDUE, mol->mol, alloc);
            size_t match_count = md_index_data_count(match_result);
            EXPECT_EQ(resname_count, match_count);
            if (resname_count != match_count) {
                md_bitfield_t bf = md_bitfield_create(alloc);
                for (size_t i = 0; i < match_count; ++i) {
                    int* atom_idx = md_index_range_beg(match_result, i);
                    md_residue_idx_t res_idx = mol->mol->atom.res_idx[atom_idx[0]];
                    md_bitfield_set_bit(&bf, res_idx);
                }
                md_bitfield_xor_inplace(&bf, &bf_ref);
                printf("Mismatch in dataset '" STR_FMT "' for residues with name '" STR_FMT "': resname_count: %zu, match_count: %zu\n", STR_ARG(mol->name), STR_ARG(res->name), resname_count, match_count);
                printf("The residues indices are:\n");
                md_bitfield_iter_t it = md_bitfield_iter_create(&bf);
                while (md_bitfield_iter_next(&it)) {
                    int idx = md_bitfield_iter_idx(&it);
                    printf("%i ", idx + 1);
                }
                printf("\n");
            }
        }
        md_vm_arena_temp_end(temp);
    }
    md_timestamp_t t1 = md_time_current();
    printf("time: %f ms\n", md_time_as_milliseconds(t1-t0));

}

UTEST(util, parse_smiles) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    {
        const char input[] = "C1=CC=CC=C1";
        md_smiles_node_t smiles[sizeof(input)];
        int64_t len = md_smiles_parse(smiles, ARRAY_SIZE(smiles), input, sizeof(input));

        ASSERT_EQ(len, 11);
        EXPECT_EQ(smiles[0].type, MD_SMILES_NODE_ATOM);
        EXPECT_EQ(smiles[0].atom.element, 6);
        
        EXPECT_EQ(smiles[1].type, MD_SMILES_NODE_BRIDGE);
        EXPECT_EQ(smiles[1].bridge.index, 1);
        
        EXPECT_EQ(smiles[2].type, MD_SMILES_NODE_BOND);
        EXPECT_EQ(smiles[2].bond.symbol, '=');

        EXPECT_EQ(smiles[3].type, MD_SMILES_NODE_ATOM);
		EXPECT_EQ(smiles[3].atom.element, 6);
        
        EXPECT_EQ(smiles[4].type, MD_SMILES_NODE_ATOM);
		EXPECT_EQ(smiles[4].atom.element, 6);
        
        EXPECT_EQ(smiles[5].type, MD_SMILES_NODE_BOND);
		EXPECT_EQ(smiles[5].bond.symbol, '=');

        EXPECT_EQ(smiles[6].type, MD_SMILES_NODE_ATOM);
		EXPECT_EQ(smiles[6].atom.element, 6);

        EXPECT_EQ(smiles[7].type, MD_SMILES_NODE_ATOM);
		EXPECT_EQ(smiles[7].atom.element, 6);

        EXPECT_EQ(smiles[8].type, MD_SMILES_NODE_BOND);
        EXPECT_EQ(smiles[8].bond.symbol, '=');

        EXPECT_EQ(smiles[9].type, MD_SMILES_NODE_ATOM);
        EXPECT_EQ(smiles[9].atom.element, 6);

		EXPECT_EQ(smiles[10].type, MD_SMILES_NODE_BRIDGE);
        EXPECT_EQ(smiles[10].bridge.index, 1);
    }

    {
        const char input[] = "C1=CC=C(CO[CH2:2])C=C1";
        md_smiles_node_t smiles[sizeof(input)];
        int64_t len = md_smiles_parse(smiles, ARRAY_SIZE(smiles), input, sizeof(input));
        
        EXPECT_EQ(len, 16);
        EXPECT_EQ(smiles[0].type, MD_SMILES_NODE_ATOM);
        EXPECT_EQ(smiles[0].atom.element, 6);
    }

    md_vm_arena_destroy(alloc);
}

UTEST(util, radix_sort) {
    uint32_t arr[] = { 1, 278, 128312745, 4, 5, 0, 12382, 26, 12, 12, 7 };
    size_t len = ARRAY_SIZE(arr);
    md_util_sort_radix_inplace_uint32(arr, len);

    for (size_t i = 0; i < len - 1; ++i) {
    	EXPECT_LE(arr[i], arr[i+1]);
    }
}
