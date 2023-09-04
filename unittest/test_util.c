#include "utest.h"
#include <string.h>

#include <md_pdb.h>
#include <md_gro.h>
#include <md_trajectory.h>
#include <md_molecule.h>
#include <md_util.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>

#include "rmsd.h"

UTEST(util, rmsd) {
    md_allocator_i* alloc = md_heap_allocator;
    const str_t path = STR(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");

    md_pdb_data_t pdb_data = {0};
    ASSERT_TRUE(md_pdb_data_parse_file(&pdb_data, path, alloc));

    md_molecule_t mol = {0};
    EXPECT_TRUE(md_pdb_molecule_init(&mol, &pdb_data, alloc));

    md_trajectory_i* traj = md_pdb_trajectory_create(path, alloc);
    EXPECT_TRUE(traj);

    const int64_t stride = mol.atom.count;
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

    for (int64_t i = 0; i < mol.atom.count; ++i) {
        xyz0[i * 3 + 0] = coord[0].x[i];
        xyz0[i * 3 + 1] = coord[0].y[i];
        xyz0[i * 3 + 2] = coord[0].z[i];

        xyz1[i * 3 + 0] = coord[1].x[i];
        xyz1[i * 3 + 1] = coord[1].y[i];
        xyz1[i * 3 + 2] = coord[1].z[i];
    }

    // Reference
    double ref_rmsd;
    fast_rmsd((double(*)[3])xyz0, (double(*)[3])xyz1, (int)mol.atom.count, &ref_rmsd);

    // Our implementation
    vec3_t com[2] = {
        md_util_compute_com_soa(coord[0].x, coord[0].y, coord[0].z, mol.atom.mass, mol.atom.count),
        md_util_compute_com_soa(coord[1].x, coord[1].y, coord[1].z, mol.atom.mass, mol.atom.count),
    };
    double rmsd = md_util_compute_rmsd(coord, com, mol.atom.mass, mol.atom.count);
    
    EXPECT_LE(fabs(ref_rmsd - rmsd), 0.1);

    md_molecule_free(&mol, alloc);
    md_pdb_trajectory_free(traj);
    md_pdb_data_free(&pdb_data, alloc);
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
        const vec3_t pos[] = {
            {1,0,0},
            {2,0,0},
            {3,0,0},
            {4,0,0},
        };
        const vec3_t pbc_ext = { 5,0,0 };

        vec3_t com = md_util_compute_com_ortho(pos, 0, ARRAY_SIZE(pos), pbc_ext);
        EXPECT_NEAR(2.5f, com.x, 1.0E-5F);
        EXPECT_EQ(0, com.y);
        EXPECT_EQ(0, com.z);
    }
    
    {
        const vec3_t pos[] = {
            {0,0,0},
            {5,0,0},
        };
        const vec3_t pbc_ext = { 5,0,0 };

        vec3_t com = md_util_compute_com_ortho(pos, 0, ARRAY_SIZE(pos), pbc_ext);
		com = vec3_deperiodize(com, (vec3_t){ 0,0,0 }, pbc_ext);
        EXPECT_NEAR(0, com.x, 1.0E-5F);
        EXPECT_EQ(0, com.y);
        EXPECT_EQ(0, com.z);
    }

    {
        const vec3_t pos[] = {
            {0,0,0},
            {4,0,0},
        };
        const vec3_t pbc_ext = { 5,0,0 };

        /*
        Here we expect the 4 to wrap around to -1,
        then added to the 0, producing a center of mass of -0.5.
        which is then placed within the period to 4.5.
        */

        vec3_t com = md_util_compute_com_ortho(pos, 0, ARRAY_SIZE(pos), pbc_ext);
        com = vec3_deperiodize(com, vec3_mul_f(pbc_ext, 0.5f), pbc_ext);
        EXPECT_NEAR(4.5f, com.x, 1.0E-5F);
        EXPECT_EQ(0, com.y);
        EXPECT_EQ(0, com.z);
    }

    {
        const vec3_t pbc_ext = { 5,0,0 };

        const vec3_t pos0[] = {
            {4,0,0},
            {5,0,0},
            {6,0,0},
            {7,0,0},
            {8,0,0},
        };

        const vec3_t pos1[] = {
            {4,0,0},
            {0,0,0},
            {1,0,0},
            {2,0,0},
            {3,0,0},
        };

        const vec3_t pos2[] = {
            {-1,0,0},
            {0,0,0},
            {1,0,0},
            {2,0,0},
            {3,0,0},
        };

        // All of the given positions, are the same under PBC. We can also see that the mid point in the sequence is 1.
        // But the sequence itself spans the entire extent of the periodic bounds.
        // What is the expected center of mass?
        // The result will differ depending on the underlying technique used.
        // If the trigonometric version is used, the points will be evenly spread across the domain, resulting in a atan2(0,0) which is not defined.
        // In such case, we fall back to picking the center of the domain.

        vec3_t com0 = md_util_compute_com_ortho(pos0, NULL, ARRAY_SIZE(pos0), pbc_ext);
        vec3_t com1 = md_util_compute_com_ortho(pos1, NULL, ARRAY_SIZE(pos1), pbc_ext);
        vec3_t com2 = md_util_compute_com_ortho(pos2, NULL, ARRAY_SIZE(pos2), pbc_ext);

        com0 = vec3_deperiodize(com0, (vec3_t){ 0,0,0 }, pbc_ext);
        com1 = vec3_deperiodize(com1, (vec3_t){ 0,0,0 }, pbc_ext);
        com2 = vec3_deperiodize(com2, (vec3_t){ 0,0,0 }, pbc_ext);
        
        EXPECT_NEAR(1.0f, com0.x, 1.0E-5F);
        EXPECT_NEAR(1.0f, com1.x, 1.0E-5F);
        EXPECT_NEAR(1.0f, com2.x, 1.0E-5F);
    }
}

UTEST(util, structure) {
    md_allocator_i* arena = md_arena_allocator_create(md_heap_allocator, MD_ARENA_ALLOCATOR_DEFAULT_PAGE_SIZE);
    md_molecule_t mol = {0};
    int64_t num_rings = 0;
    int64_t num_structures = 0;

    ASSERT_TRUE(md_gro_molecule_api()->init_from_file(&mol, STR(MD_UNITTEST_DATA_DIR "/nucleotides.gro"), arena));
    md_util_postprocess_molecule(&mol, arena, MD_UTIL_POSTPROCESS_ELEMENT_BIT | MD_UTIL_POSTPROCESS_BOND_BIT | MD_UTIL_POSTPROCESS_CONNECTIVITY_BIT);

    num_structures = md_index_data_count(mol.structures);
    EXPECT_EQ(num_structures, 2);

    num_rings = md_index_data_count(mol.rings);
    EXPECT_EQ(num_rings, 4);

    

    ASSERT_TRUE(md_pdb_molecule_api()->init_from_file(&mol, STR(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), arena));
    md_util_postprocess_molecule(&mol, arena, MD_UTIL_POSTPROCESS_ELEMENT_BIT | MD_UTIL_POSTPROCESS_BOND_BIT | MD_UTIL_POSTPROCESS_CONNECTIVITY_BIT);

    num_structures = md_index_data_count(mol.structures);
	EXPECT_EQ(num_structures, 1);
    
    num_rings = md_index_data_count(mol.rings);
    EXPECT_EQ(num_rings, 0);
    
    md_arena_allocator_reset(arena);
    MEMSET(&mol, 0, sizeof(md_molecule_t));

    
    
    ASSERT_TRUE(md_gro_molecule_api()->init_from_file(&mol, STR(MD_UNITTEST_DATA_DIR "/pftaa.gro"), arena));
    md_util_postprocess_molecule(&mol, arena, MD_UTIL_POSTPROCESS_ELEMENT_BIT | MD_UTIL_POSTPROCESS_BOND_BIT | MD_UTIL_POSTPROCESS_CONNECTIVITY_BIT);
    
    num_structures = md_index_data_count(mol.structures);
    EXPECT_EQ(num_structures, 1);

    num_rings = md_index_data_count(mol.rings);
    EXPECT_EQ(num_rings, 5);

    md_arena_allocator_reset(arena);
    MEMSET(&mol, 0, sizeof(md_molecule_t));
    
    ASSERT_TRUE(md_gro_molecule_api()->init_from_file(&mol, STR(MD_UNITTEST_DATA_DIR "/centered.gro"), arena));
    md_util_postprocess_molecule(&mol, arena, MD_UTIL_POSTPROCESS_ELEMENT_BIT | MD_UTIL_POSTPROCESS_BOND_BIT | MD_UTIL_POSTPROCESS_CONNECTIVITY_BIT | MD_UTIL_POSTPROCESS_CHAINS_BIT);
    
    num_structures = md_index_data_count(mol.structures);
    const int64_t expected_count = mol.chain.count + 61;
    EXPECT_EQ(num_structures, expected_count);

    num_rings = md_index_data_count(mol.rings);
    EXPECT_EQ(num_rings, 2076);
    
    md_arena_allocator_destroy(arena);
}