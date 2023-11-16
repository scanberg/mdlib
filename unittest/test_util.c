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

#include "rmsd.h"

struct util {
    md_allocator_i* alloc;
    md_molecule_t mol_ala;
    md_molecule_t mol_pftaa;
    md_molecule_t mol_nucleotides;
    md_molecule_t mol_centered;

    md_trajectory_i* traj_ala;
};

UTEST_F_SETUP(util) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(1));
    utest_fixture->alloc = alloc;

    md_pdb_molecule_api()->init_from_file(&utest_fixture->mol_ala, STR(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), alloc);
    md_util_postprocess_molecule(&utest_fixture->mol_ala, alloc, MD_UTIL_POSTPROCESS_ALL);

    md_gro_molecule_api()->init_from_file(&utest_fixture->mol_pftaa, STR(MD_UNITTEST_DATA_DIR "/pftaa.gro"), alloc);
    md_util_postprocess_molecule(&utest_fixture->mol_pftaa, alloc, MD_UTIL_POSTPROCESS_ALL);

    md_gro_molecule_api()->init_from_file(&utest_fixture->mol_nucleotides, STR(MD_UNITTEST_DATA_DIR "/nucleotides.gro"), alloc);
    md_util_postprocess_molecule(&utest_fixture->mol_nucleotides, alloc, MD_UTIL_POSTPROCESS_ALL);

    md_gro_molecule_api()->init_from_file(&utest_fixture->mol_centered, STR(MD_UNITTEST_DATA_DIR "/centered.gro"), alloc);
    md_util_postprocess_molecule(&utest_fixture->mol_centered, alloc, MD_UTIL_POSTPROCESS_ALL);

    utest_fixture->traj_ala = md_pdb_trajectory_create(STR(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), alloc);
}

UTEST_F_TEARDOWN(util) {
    md_arena_allocator_destroy(utest_fixture->alloc);
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
    md_vec3_soa_t coord[2] = {
        {
            mem + stride * 0,
            mem + stride * 1,
            mem + stride * 2,
        },
        {
            mem + stride * 3,
            mem + stride * 4,
            mem + stride * 5,
        },
    };

    double* xyz0 = (double*)(mem + stride * 6);
    double* xyz1 = (double*)(mem + stride * 6) + stride * 3;

    md_trajectory_load_frame(traj, 0, NULL, coord[0].x, coord[0].y, coord[0].z);
    md_trajectory_load_frame(traj, 1, NULL, coord[1].x, coord[1].y, coord[1].z);

    for (int64_t i = 0; i < mol->atom.count; ++i) {
        xyz0[i * 3 + 0] = coord[0].x[i];
        xyz0[i * 3 + 1] = coord[0].y[i];
        xyz0[i * 3 + 2] = coord[0].z[i];

        xyz1[i * 3 + 0] = coord[1].x[i];
        xyz1[i * 3 + 1] = coord[1].y[i];
        xyz1[i * 3 + 2] = coord[1].z[i];
    }

    // Reference
    double ref_rmsd;
    fast_rmsd((double(*)[3])xyz0, (double(*)[3])xyz1, (int)mol->atom.count, &ref_rmsd);

    // Our implementation
    vec3_t com[2] = {
        md_util_compute_com(coord[0].x, coord[0].y, coord[0].z, mol->atom.mass, 0, mol->atom.count),
        md_util_compute_com(coord[1].x, coord[1].y, coord[1].z, mol->atom.mass, 0, mol->atom.count),
    };
    double rmsd = md_util_compute_rmsd(coord, com, mol->atom.mass, mol->atom.count);
    
    EXPECT_LE(fabs(ref_rmsd - rmsd), 0.1);

    md_free(alloc, mem, mem_size);
}

UTEST(util, com) {
    /* DISCLAIMER
        Computing the center of mass is a bit tricky, because of the periodic boundary conditions.
        The problem occurs when we have structures which have an extent which covers more than half the box.
        There are many different variations of how to handle this and it seems none of them are perfect.
        Unless you handpick your algorithm based on some external knowledge of the structure you are dealing with.

        In this case, we have settled on an algorithm which iteratively accumulates the center of mass by deperiodizing
        new points with respect to the previous point. This works in most cases when dealing with long structures, given that
        the indices or subsequent points are not too far apart, which is almost always the case in real world datasets.
    */
    {
        const vec4_t xyzw[] = {
            {1,0,0,1},
            {2,0,0,1},
            {3,0,0,1},
            {4,0,0,1},
        };
        const vec3_t pbc_ext = { 5,0,0 };

        vec3_t com = md_util_compute_com_vec4_ortho(xyzw, 0, ARRAY_SIZE(xyzw), pbc_ext);
        EXPECT_NEAR(2.5f, com.x, 1.0E-5F);
        EXPECT_EQ(0, com.y);
        EXPECT_EQ(0, com.z);
    }
    
    {
        const vec4_t xyzw[] = {
            {0,0,0,1},
            {5,0,0,1},
        };
        const vec3_t pbc_ext = { 5,0,0 };

        vec3_t com = md_util_compute_com_vec4_ortho(xyzw, 0, ARRAY_SIZE(xyzw), pbc_ext);
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
        const vec3_t pbc_ext = { 5,0,0 };

        /*
        Here we expect the 4 to wrap around to -1,
        then added to the 0, producing a center of mass of -0.5.
        which is then placed within the period to 4.5.
        */

        vec3_t com = md_util_compute_com_vec4_ortho(xyzw, 0, ARRAY_SIZE(xyzw), pbc_ext);
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

        vec3_t com0 = md_util_compute_com_vec4_ortho(pos0, 0, ARRAY_SIZE(pos0), pbc_ext);
        vec3_t com1 = md_util_compute_com_vec4_ortho(pos1, 0, ARRAY_SIZE(pos1), pbc_ext);
        vec3_t com2 = md_util_compute_com_vec4_ortho(pos2, 0, ARRAY_SIZE(pos2), pbc_ext);

        com0 = vec3_deperiodize(com0, (vec3_t){ 0,0,0 }, pbc_ext);
        com1 = vec3_deperiodize(com1, (vec3_t){ 0,0,0 }, pbc_ext);
        com2 = vec3_deperiodize(com2, (vec3_t){ 0,0,0 }, pbc_ext);
        
        EXPECT_NEAR(0.5f, com0.x, 1.0E-5F);
        EXPECT_NEAR(0.5f, com1.x, 1.0E-5F);
        EXPECT_NEAR(0.5f, com2.x, 1.0E-5F);
    }
}

UTEST_F(util, structures) {
    int64_t num_structures = 0;

    num_structures = md_index_data_count(utest_fixture->mol_nucleotides.structures);
    EXPECT_EQ(num_structures, 2);

    num_structures = md_index_data_count(utest_fixture->mol_ala.structures);
	EXPECT_EQ(num_structures, 1);

	num_structures = md_index_data_count(utest_fixture->mol_pftaa.structures);
	EXPECT_EQ(num_structures, 1);

	num_structures = md_index_data_count(utest_fixture->mol_centered.structures);
    const int64_t expected_count = utest_fixture->mol_centered.chain.count + 61;
    EXPECT_EQ(num_structures, expected_count);
}

UTEST_F(util, rings_common) {
    int64_t num_rings = 0;

    num_rings = md_index_data_count(utest_fixture->mol_nucleotides.rings);
    EXPECT_EQ(num_rings, 4);
    
    num_rings = md_index_data_count(utest_fixture->mol_ala.rings);
    EXPECT_EQ(num_rings, 0);
   
    num_rings = md_index_data_count(utest_fixture->mol_pftaa.rings);
    EXPECT_EQ(num_rings, 5);

    num_rings = md_index_data_count(utest_fixture->mol_centered.rings);
    EXPECT_EQ(num_rings, 2076);
}

UTEST(util, rings_c60) {
	md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(1));
	md_molecule_t mol = {0};
	md_pdb_molecule_api()->init_from_file(&mol, STR(MD_UNITTEST_DATA_DIR "/c60.pdb"), alloc);
	md_util_postprocess_molecule(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

	EXPECT_EQ(mol.atom.count, 60);
	EXPECT_EQ(mol.bond.count, 90);

    const int64_t num_rings = md_index_data_count(mol.rings);
    EXPECT_EQ(num_rings, 32);

    const int64_t num_structures = md_index_data_count(mol.structures);
    EXPECT_EQ(num_structures, 1);

	md_arena_allocator_destroy(alloc);
}

UTEST(util, rings_14kr) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(1));
    md_molecule_t mol = {0};
    md_pdb_molecule_api()->init_from_file(&mol, STR(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), alloc);
    md_util_postprocess_molecule(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    const int64_t num_rings = md_index_data_count(mol.rings);
    EXPECT_EQ(num_rings, 207);

    md_arena_allocator_destroy(alloc);
}

UTEST(util, rings_trytophan_pdb) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(1));
    md_molecule_t mol = {0};
    md_pdb_molecule_api()->init_from_file(&mol, STR(MD_UNITTEST_DATA_DIR "/tryptophan.pdb"), alloc);
    md_util_postprocess_molecule(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    const int64_t num_rings = md_index_data_count(mol.rings);
    EXPECT_EQ(num_rings, 2);

    const int64_t num_structures = md_index_data_count(mol.structures);
    EXPECT_EQ(num_structures, 1);

    md_arena_allocator_destroy(alloc);
}

UTEST(util, rings_trytophan_xyz) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(1));
    md_molecule_t mol = {0};
    md_xyz_molecule_api()->init_from_file(&mol, STR(MD_UNITTEST_DATA_DIR "/tryptophan.xyz"), alloc);
    md_util_postprocess_molecule(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    const int64_t num_rings = md_index_data_count(mol.rings);
    EXPECT_EQ(num_rings, 2);

    const int64_t num_structures = md_index_data_count(mol.structures);
    EXPECT_EQ(num_structures, 1);

    md_arena_allocator_destroy(alloc);
}

UTEST(util, rings_full) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(1));
    md_molecule_t mol = {0};
    md_xyz_molecule_api()->init_from_file(&mol, STR(MD_UNITTEST_DATA_DIR "/full.xyz"), alloc);
    md_util_postprocess_molecule(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    const int64_t num_rings = md_index_data_count(mol.rings);
    EXPECT_EQ(num_rings, 195);

    const int64_t num_structures = md_index_data_count(mol.structures);
    EXPECT_EQ(num_structures, 1);

    md_arena_allocator_destroy(alloc);
}

UTEST(util, rings_ciprofloxacin) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(1));
    md_molecule_t mol = {0};
    md_pdb_molecule_api()->init_from_file(&mol, STR(MD_UNITTEST_DATA_DIR "/ciprofloxacin.pdb"), alloc);
    md_util_postprocess_molecule(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    const int64_t num_rings = md_index_data_count(mol.rings);
    ASSERT_EQ(num_rings, 4);
    EXPECT_EQ(md_index_range_size(mol.rings, 0), 3);
    EXPECT_EQ(md_index_range_size(mol.rings, 1), 6);
    EXPECT_EQ(md_index_range_size(mol.rings, 2), 6);
    EXPECT_EQ(md_index_range_size(mol.rings, 3), 6);

    md_arena_allocator_destroy(alloc);
}

UTEST_F(util, structure_matching_amyloid) {
    md_allocator_i* alloc = utest_fixture->alloc;
    md_molecule_t* mol = &utest_fixture->mol_centered;

#if 0
    {
        // Test for the chains
        const int ref_structure_idx = 1;
        const int* ref_idx = md_index_range_beg(mol->structures, ref_structure_idx);
        const int64_t ref_size = md_index_range_size(mol->structures, ref_structure_idx);
    
        md_index_data_t result = md_util_structure_find_equivalent(ref_idx, ref_size, mol, alloc);
        const int result_count = (int)md_index_data_count(result);
        EXPECT_EQ(result_count, 253);
    }
#endif
#if 0
    {
        // Test for the PFTAAs
        const int ref_structure_idx = 253;
        const int* ref_idx = md_index_range_beg(mol->structures, ref_structure_idx);
        const int64_t ref_size = md_index_range_size(mol->structures, ref_structure_idx);

        md_timestamp_t t0 = md_time_current();
        md_index_data_t result = md_util_structure_find_equivalent(ref_idx, ref_size, mol, alloc);
        md_timestamp_t t1 = md_time_current();
        printf("time: %f ms\n", md_time_as_milliseconds(t1-t0));
        const int64_t result_count = md_index_data_count(result);
        EXPECT_EQ(result_count, 61);
    }
#endif
}

UTEST_F(util, structure_matching_PFTAA) {
    md_allocator_i* alloc = utest_fixture->alloc;
    md_molecule_t* mol = &utest_fixture->mol_pftaa;
    {
#if 0
        // Rings, represents a ring within the molecule
        const int ref_idx[] = {19,20,21,22,24};
        const int64_t ref_size = ARRAY_SIZE(ref_idx);

        md_index_data_t result = md_util_structure_find_equivalent(ref_idx, ref_size, mol, alloc);
        const int64_t result_count = md_index_data_count(result);
        EXPECT_EQ(result_count, 5);
#endif
#if 0
        for (int64_t i = 0; i < result_count; ++i) {
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
#if 0
        // Symmetry, represents half molecule minus the ring which holds the two symmetric parts together.
        const int ref_idx[] = {0,1,2,3,4,5,6,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,46,47,48};
        const int64_t ref_len = ARRAY_SIZE(ref_idx);

        md_index_data_t result = md_util_structure_find_equivalent(ref_idx, ref_len, mol, alloc);
        const int64_t result_count = md_index_data_count(result);
        EXPECT_EQ(result_count, 2);

        md_index_data_free(&result, alloc);
#endif
    }
}

UTEST_F(util, structure_matching_smiles_ala) {
    md_allocator_i* alloc = utest_fixture->alloc;
    
    const char ALANINE1[] = "N[C@@H]([CH3])C(=O)";
    const char ALANINE2[] = "N[C@@H]([CH3])C(-O)(-O)";

    const char ARGININE[] = "N[C@@H](CCCNC(=N)N)C(=O)";
    const char ASPARAGINE[] = "C(C(C(=O)O)N)C(=O)N";
    const char ASPARTATE[] = "C(C(C(=O)[O-])N)C(=O)[O-]";
    const char CYSTEINE[] = "C(C(C(=O)O)N)S";
    const char GLUTAMATE[] = "C(CC(=O)[O-])C(C(=O)[O-])N";
    const char GLUTAMINE[] = "C(CC(=O)[O-])C(C(=O)N)N";
    const char GLYCINE[] = "C(C(=O)O)N";
    const char HISTIDINE[] = "C(CC1=CNC=N1)C(C(=O)O)N";
    const char ISOLEUCINE[] = "CC(C)C(C(=O)O)N";
    const char LEUCINE[] = "CC(C)CC(C(=O)O)N";
    const char LYSINE[] = "C(CCN)CC(C(=O)O)N";
    const char METHIONINE[] = "CSCC(C(=O)O)N";
    const char PHENYLALANINE[] = "C1=CC=C(C=C1)CC(C(=O)O)N";
    const char PROLINE[] = "C1CC(NC1)C(=O)O";
    const char SERINE[] = "C(C(C(=O)O)O)N";
    const char THREONINE[] = "C(C(C(C(=O)O)O)O)N";
    const char TRYPTOPHAN[] = "C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N";
    const char TYROSINE[] = "C1=CC=C(C=C1)C(C(=O)O)N";
    const char VALINE[] = "CC(C)C(C(=O)O)N";

    const char SELENOCYSTINE[] = "C(C(C(=O)O)N)[Se]";
    const char PYRROLYSINE[] = "CC1CC=NC1C(=O)NCCCCC(C(=O)O)N";

    {
#if 0
        const md_molecule_t* mol = &utest_fixture->mol_ala;
        md_index_data_t res1 = md_util_structure_find_equivalent_smiles((str_t){ALANINE1, sizeof(ALANINE1)}, mol, alloc);
        md_index_data_t res2 = md_util_structure_find_equivalent_smiles((str_t){ALANINE2, sizeof(ALANINE2)}, mol, alloc);
        const int64_t res_count1 = md_index_data_count(res1);
        const int64_t res_count2 = md_index_data_count(res2);
        const int64_t count = res_count1 + res_count2;
        EXPECT_EQ(count, 15);
#endif
    }
    {
#if 0
        const md_molecule_t* mol = &utest_fixture->mol_centered;
        md_index_data_t res1 = md_util_structure_find_equivalent_smiles((str_t){ALANINE1, sizeof(ALANINE1)}, mol, alloc);
        md_index_data_t res2 = md_util_structure_find_equivalent_smiles((str_t){ALANINE2, sizeof(ALANINE2)}, mol, alloc);
        const int64_t res_count1 = md_index_data_count(res1);
        const int64_t res_count2 = md_index_data_count(res2);
        const int64_t count = res_count1 + res_count2;
        EXPECT_EQ(count, 1012);
        
        for (int64_t i = 0; i < res_count1; ++i) {
            int* it = md_index_range_beg(res1, i);
            printf("[%d] resname: %s\n", utest_fixture->mol_centered.atom.res_idx[*it], utest_fixture->mol_centered.atom.resname[*it].buf);
        }

        for (int64_t i = 0; i < res_count2; ++i) {
            int* it = md_index_range_beg(res2, i);
            printf("[%d] resname: %s\n", utest_fixture->mol_centered.atom.res_idx[*it], utest_fixture->mol_centered.atom.resname[*it].buf);
        }
#endif
    }
}

UTEST(util, parse_smiles) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(1));

    {
        const char input[] = "C1=CC=CC=C1";
        md_smiles_node_t smiles[sizeof(input)];
        int64_t len = md_smiles_parse(smiles, ARRAY_SIZE(smiles), input, sizeof(input));

        ASSERT_EQ(len, 11);
        EXPECT_EQ(smiles[0].type, MD_SMILES_NODE_ATOM);
        EXPECT_EQ(smiles[0].atom.symbol[0], 'C');
        
        EXPECT_EQ(smiles[1].type, MD_SMILES_NODE_RING_CLOSURE);
        EXPECT_EQ(smiles[1].ring.index, 1);
        
        EXPECT_EQ(smiles[2].type, MD_SMILES_NODE_BOND);
        EXPECT_EQ(smiles[2].bond.symbol, '=');

        EXPECT_EQ(smiles[3].type, MD_SMILES_NODE_ATOM);
		EXPECT_EQ(smiles[3].atom.symbol[0], 'C');
        
        EXPECT_EQ(smiles[4].type, MD_SMILES_NODE_ATOM);
		EXPECT_EQ(smiles[4].atom.symbol[0], 'C');
        
        EXPECT_EQ(smiles[5].type, MD_SMILES_NODE_BOND);
		EXPECT_EQ(smiles[5].bond.symbol, '=');

        EXPECT_EQ(smiles[6].type, MD_SMILES_NODE_ATOM);
		EXPECT_EQ(smiles[6].atom.symbol[0], 'C');

        EXPECT_EQ(smiles[7].type, MD_SMILES_NODE_ATOM);
		EXPECT_EQ(smiles[7].atom.symbol[0], 'C');

        EXPECT_EQ(smiles[8].type, MD_SMILES_NODE_BOND);
        EXPECT_EQ(smiles[8].bond.symbol, '=');

        EXPECT_EQ(smiles[9].type, MD_SMILES_NODE_ATOM);
        EXPECT_EQ(smiles[9].atom.symbol[0], 'C');

		EXPECT_EQ(smiles[10].type, MD_SMILES_NODE_RING_CLOSURE);
        EXPECT_EQ(smiles[10].ring.index, 1);

        md_array_free(smiles, alloc);
    }

    {
        const char input[] = "C1=CC=C(CO[CH2:2])C=C1";
        md_smiles_node_t smiles[sizeof(input)];
        int64_t len = md_smiles_parse(smiles, ARRAY_SIZE(smiles), input, sizeof(input));
        
        EXPECT_EQ(len, 16);
        EXPECT_EQ(smiles[0].type, MD_SMILES_NODE_ATOM);
        EXPECT_EQ(smiles[0].atom.symbol[0], 'C');



        md_array_free(smiles, alloc);
    }

    md_arena_allocator_destroy(alloc);
}