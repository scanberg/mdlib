#include "utest.h"
#include <string.h>

#include <md_pdb.h>
#include <md_gro.h>
#include <md_xyz.h>
#include <md_mmcif.h>
#include <md_trajectory.h>
#include <md_molecule.h>
#include <md_util.h>
#include <md_smiles.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_str_builder.h>
#include <core/md_bitfield.h>
#include <core/md_hash.h>

#include "rmsd.h"

struct util {
    md_allocator_i* alloc;
    md_system_t mol_ala;
    md_system_t mol_pftaa;
    md_system_t mol_nucleotides;
    md_system_t mol_centered;
    md_system_t mol_dna;
    md_system_t mol_trp;
    md_system_t mol_aspirine;

    md_system_t mol_1fez;
    md_system_t mol_2or2;
    md_system_t mol_1k4r;
    md_system_t mol_8g7u;

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

    md_mmcif_molecule_api()->init_from_file(&utest_fixture->mol_1fez, STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif"), NULL, alloc);
    md_util_molecule_postprocess(&utest_fixture->mol_1fez, alloc, MD_UTIL_POSTPROCESS_ALL);

    md_mmcif_molecule_api()->init_from_file(&utest_fixture->mol_2or2, STR_LIT(MD_UNITTEST_DATA_DIR "/2or2.cif"), NULL, alloc);
    md_util_molecule_postprocess(&utest_fixture->mol_2or2, alloc, MD_UTIL_POSTPROCESS_ALL);

    md_pdb_molecule_api()->init_from_file(&utest_fixture->mol_1k4r, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), NULL, alloc);
    md_util_molecule_postprocess(&utest_fixture->mol_1k4r, alloc, MD_UTIL_POSTPROCESS_ALL);

    md_mmcif_molecule_api()->init_from_file(&utest_fixture->mol_8g7u, STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif"), NULL, alloc);
    md_util_molecule_postprocess(&utest_fixture->mol_8g7u, alloc, MD_UTIL_POSTPROCESS_ALL);

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

UTEST_F(util, inst) {
    EXPECT_EQ(1,   utest_fixture->mol_ala.inst.count);
	EXPECT_EQ(1,   utest_fixture->mol_pftaa.inst.count);
	EXPECT_EQ(2,   utest_fixture->mol_nucleotides.inst.count);
	EXPECT_EQ(253 + 61, utest_fixture->mol_centered.inst.count);

    const md_system_t* sys = &utest_fixture->mol_centered;
    ASSERT(sys->inst.count > 253);

    size_t ref_size = md_system_inst_atom_count(sys, 0);
    for (size_t i = 1; i < 253; ++i) {
        size_t size = md_system_inst_atom_count(sys, i);
        EXPECT_EQ(ref_size, size);
    }
}

UTEST_F(util, backbone) {
	EXPECT_EQ(0,        utest_fixture->mol_pftaa.protein_backbone.range.count);
	EXPECT_EQ(0,        utest_fixture->mol_nucleotides.protein_backbone.range.count);
    EXPECT_EQ(1,        utest_fixture->mol_ala.protein_backbone.range.count);
    EXPECT_EQ(15,       utest_fixture->mol_ala.protein_backbone.segment.count);
	EXPECT_EQ(253,      utest_fixture->mol_centered.protein_backbone.range.count);
    EXPECT_EQ(10626,    utest_fixture->mol_centered.protein_backbone.segment.count); // Should be equal to the total count of residues in chains
}

UTEST_F(util, structure) {
    size_t num_structures_pftaa = md_index_data_num_ranges(&utest_fixture->mol_pftaa.structure);
    EXPECT_EQ(1, num_structures_pftaa);
    size_t num_structures_nucleotides = md_index_data_num_ranges(&utest_fixture->mol_nucleotides.structure);
    EXPECT_EQ(2, num_structures_nucleotides);
    size_t num_structures_ala = md_index_data_num_ranges(&utest_fixture->mol_ala.structure);
	EXPECT_EQ(1, num_structures_ala);
	size_t num_structures_centered = md_index_data_num_ranges(&utest_fixture->mol_centered.structure);
	EXPECT_EQ(253+61, num_structures_centered);
}

UTEST_F(util, rmsd) {
    md_allocator_i* alloc = utest_fixture->alloc;
    md_system_t* mol = &utest_fixture->mol_ala;
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

    float* w = md_alloc(alloc, sizeof(float) * stride);
	md_atom_extract_masses(w, 0, mol->atom.count, &mol->atom);

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
    const float* const cw[2] = { w, w };

    vec3_t com[2] = {
        md_util_com_compute(x[0], y[0], z[0], w, 0, mol->atom.count, 0),
        md_util_com_compute(x[1], y[1], z[1], w, 0, mol->atom.count, 0),
    };
    double rmsd = md_util_rmsd_compute(cx, cy, cz, cw, 0, mol->atom.count, com);
    
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
		com = vec3_deperiodize_ortho(com, (vec3_t){ 0,0,0 }, pbc_ext);
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
        com = vec3_deperiodize_ortho(com, vec3_mul_f(pbc_ext, 0.5f), pbc_ext);
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

        com0 = vec3_deperiodize_ortho(com0, (vec3_t){ 0,0,0 }, pbc_ext);
        com1 = vec3_deperiodize_ortho(com1, (vec3_t){ 0,0,0 }, pbc_ext);
        com2 = vec3_deperiodize_ortho(com2, (vec3_t){ 0,0,0 }, pbc_ext);
        
        EXPECT_NEAR(0.5f, com0.x, 1.0E-5F);
        EXPECT_NEAR(0.5f, com1.x, 1.0E-5F);
        EXPECT_NEAR(0.5f, com2.x, 1.0E-5F);
    }
}

UTEST_F(util, structures) {
    size_t num_structures = 0;

    num_structures = md_index_data_num_ranges(&utest_fixture->mol_nucleotides.structure);
    EXPECT_EQ(num_structures, 2);

    num_structures = md_index_data_num_ranges(&utest_fixture->mol_ala.structure);
	EXPECT_EQ(num_structures, 1);

	num_structures = md_index_data_num_ranges(&utest_fixture->mol_pftaa.structure);
	EXPECT_EQ(num_structures, 1);

	num_structures = md_index_data_num_ranges(&utest_fixture->mol_centered.structure);
    EXPECT_EQ(num_structures, 253 + 61); // Chains + PFTAA
}

UTEST_F(util, rings_common) {
    int64_t num_rings = 0;

    num_rings = md_index_data_num_ranges(&utest_fixture->mol_nucleotides.ring);
    EXPECT_EQ(num_rings, 4);
    
    num_rings = md_index_data_num_ranges(&utest_fixture->mol_ala.ring);
    EXPECT_EQ(num_rings, 0);
   
    num_rings = md_index_data_num_ranges(&utest_fixture->mol_pftaa.ring);
    EXPECT_EQ(num_rings, 5);

    num_rings = md_index_data_num_ranges(&utest_fixture->mol_centered.ring);
    EXPECT_EQ(num_rings, 2076);
}

UTEST(util, rings_c60) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
	md_system_t mol = {0};
	md_pdb_molecule_api()->init_from_file(&mol, STR_LIT(MD_UNITTEST_DATA_DIR "/c60.pdb"), NULL, alloc);
	md_util_molecule_postprocess(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

	EXPECT_EQ(mol.atom.count, 60);
	EXPECT_EQ(mol.bond.count, 90);

    const size_t num_rings = md_index_data_num_ranges(&mol.ring);
    EXPECT_EQ(num_rings, 32);

    const size_t num_structures = md_index_data_num_ranges(&mol.structure);
    EXPECT_EQ(num_structures, 1);

    md_vm_arena_destroy(alloc);
}

UTEST(util, rings_c720) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t mol = {0};
    md_xyz_molecule_api()->init_from_file(&mol, STR_LIT(MD_UNITTEST_DATA_DIR "/c720.xyz"), NULL, alloc);
    md_util_molecule_postprocess(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    EXPECT_EQ(mol.atom.count, 720);
    EXPECT_EQ(mol.bond.count, 1080);

    const size_t num_rings = md_index_data_num_ranges(&mol.ring);
    EXPECT_EQ(num_rings, 362);

    const size_t num_structures = md_index_data_num_ranges(&mol.structure);
    EXPECT_EQ(num_structures, 1);

    md_vm_arena_destroy(alloc);
}

UTEST(util, rings_14kr) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t mol = {0};
    md_pdb_molecule_api()->init_from_file(&mol, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), NULL, alloc);
    md_util_molecule_postprocess(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    const size_t num_rings = md_index_data_num_ranges(&mol.ring);
    EXPECT_EQ(num_rings, 207);

    md_vm_arena_destroy(alloc);
}

UTEST(util, rings_trytophan_pdb) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t mol = {0};
    md_pdb_molecule_api()->init_from_file(&mol, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan.pdb"), NULL, alloc);
    md_util_molecule_postprocess(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    const size_t num_rings = md_index_data_num_ranges(&mol.ring);
    EXPECT_EQ(num_rings, 2);

    const size_t num_structures = md_index_data_num_ranges(&mol.structure);
    EXPECT_EQ(num_structures, 1);

    md_vm_arena_destroy(alloc);
}

UTEST(util, rings_trytophan_xyz) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t mol = {0};
    md_xyz_molecule_api()->init_from_file(&mol, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan.xyz"), NULL, alloc);
    md_util_molecule_postprocess(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    const size_t num_rings = md_index_data_num_ranges(&mol.ring);
    EXPECT_EQ(num_rings, 2);

    const size_t num_structures = md_index_data_num_ranges(&mol.structure);
    EXPECT_EQ(num_structures, 1);

    md_vm_arena_destroy(alloc);
}

UTEST(util, rings_full) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t mol = {0};
    md_xyz_molecule_api()->init_from_file(&mol, STR_LIT(MD_UNITTEST_DATA_DIR "/full.xyz"), NULL, alloc);
    md_util_molecule_postprocess(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    const size_t num_rings = md_index_data_num_ranges(&mol.ring);
    EXPECT_EQ(num_rings, 195);

    const size_t num_structures = md_index_data_num_ranges(&mol.structure);
    EXPECT_EQ(num_structures, 1);

    md_vm_arena_destroy(alloc);
}

UTEST(util, rings_ciprofloxacin) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t mol = {0};
    md_pdb_molecule_api()->init_from_file(&mol, STR_LIT(MD_UNITTEST_DATA_DIR "/ciprofloxacin.pdb"), NULL, alloc);
    md_util_molecule_postprocess(&mol, alloc, MD_UTIL_POSTPROCESS_ALL);

    const int64_t num_rings = md_index_data_num_ranges(&mol.ring);
    ASSERT_EQ(num_rings, 4);
    EXPECT_EQ(md_index_range_size(&mol.ring, 0), 3);
    EXPECT_EQ(md_index_range_size(&mol.ring, 1), 6);
    EXPECT_EQ(md_index_range_size(&mol.ring, 2), 6);
    EXPECT_EQ(md_index_range_size(&mol.ring, 3), 6);

    md_vm_arena_destroy(alloc);
}

UTEST_F(util, structure_matching_amyloid_chain) {
    md_allocator_i* alloc = utest_fixture->alloc;
    md_system_t* mol = &utest_fixture->mol_centered;

#if 1
    {
        // Test for the chains
        const int ref_structure_idx = 0;
        int*   ref_idx = md_index_range_beg(&mol->structure,  ref_structure_idx);
        size_t ref_len = md_index_range_size(&mol->structure, ref_structure_idx);

        // Prune Hydrogen
        if (true) {
            md_array(int) new_idx = 0;
            for (size_t i = 0; i < ref_len; ++i) {
                if (md_atom_atomic_number(&mol->atom, ref_idx[i]) > 1) {
                    md_array_push(new_idx, ref_idx[i], alloc);
                }
            }
            ref_idx = new_idx;
            ref_len = md_array_size(new_idx);
        }

        md_timestamp_t t0 = md_time_current();
        md_index_data_t result = md_util_match_by_element(ref_idx, ref_len, MD_UTIL_MATCH_MODE_FIRST, MD_UTIL_MATCH_LEVEL_INSTANCE, mol, alloc);
        md_timestamp_t t1 = md_time_current();
        printf("time: %f ms\n", md_time_as_milliseconds(t1 - t0));
        size_t result_count = md_index_data_num_ranges(&result);
        printf("result count: %zu\n", result_count);
        EXPECT_EQ(result_count, 253);
    }
#endif
}

UTEST_F(util, structure_matching_PFTAA) {
    md_allocator_i* alloc = utest_fixture->alloc;
    md_system_t* mol = &utest_fixture->mol_centered;

#if 1
    {
        // Test for the PFTAAs
        const int ref_structure_idx = 253;
        const int*   ref_idx  = md_index_range_beg (&mol->structure, ref_structure_idx);
        const size_t ref_size = md_index_range_size(&mol->structure, ref_structure_idx);

        md_timestamp_t t0 = md_time_current();
        md_index_data_t result = md_util_match_by_element(ref_idx, ref_size, MD_UTIL_MATCH_MODE_FIRST, MD_UTIL_MATCH_LEVEL_COMPONENT, mol, alloc);
        md_timestamp_t t1 = md_time_current();
        printf("time: %f ms\n", md_time_as_milliseconds(t1 - t0));
        size_t result_count = md_index_data_num_ranges(&result);
        printf("result count: %zu\n", result_count);
        EXPECT_EQ(result_count, 61);
    }
#endif
}

UTEST_F(util, structure_matching_PFTAA_ring) {
    md_allocator_i* alloc = utest_fixture->alloc;
    md_system_t* mol = &utest_fixture->mol_pftaa;
    {
#if 1
        // Rings, represents a ring within the molecule
        const int ref_idx[] = {19,20,21,22,24};
        const size_t ref_size = ARRAY_SIZE(ref_idx);

        md_index_data_t result = md_util_match_by_element(ref_idx, ref_size, MD_UTIL_MATCH_MODE_UNIQUE, MD_UTIL_MATCH_LEVEL_STRUCTURE, mol, alloc);
        const size_t result_count = md_index_data_num_ranges(&result);
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
        size_t result_count = md_index_data_num_ranges(&result);
        EXPECT_EQ(result_count, 2);

        md_index_data_free(&result);
#endif
    }
}


UTEST_F(util, structure_matching_smiles) {
    md_allocator_i* alloc = utest_fixture->alloc;
    
    // These are smiles patterns for the common amino acids and nucleotides

// These are the common parts of the amino acids (with the exception of Glycine and Proline)
#define AA_INT "[NH][C@@H](CO)"
#define AA_NT  "[NH3][C@@H](CO)"
#define AA_CT  "[NH][C@@H](C(O)O)"

    const str_t ALANINE         = STR_LIT(AA_INT "[CH3]");
    const str_t ALANINE_NT      = STR_LIT(AA_NT  "[CH3]");
    const str_t ALANINE_CT      = STR_LIT(AA_CT  "[CH3]");

    const str_t ARGININE        = STR_LIT(AA_INT "[CH2][CH2][CH2][NH]C([NH2])[NH2]");
    const str_t ARGININE_NT     = STR_LIT(AA_NT  "[CH2][CH2][CH2][NH]C([NH2])[NH2]");
    const str_t ARGININE_CT     = STR_LIT(AA_CT  "[CH2][CH2][CH2][NH]C([NH2])[NH2]");

    const str_t GLYCINE         = STR_LIT("[NH][CH2]CO");
    const str_t GLYCINE_NT      = STR_LIT("[NH3][CH2]CO");
    const str_t GLYCINE_CT      = STR_LIT("[NH][CH2]C(O)O");

    const str_t ASPARAGINE      = STR_LIT(AA_INT "[CH2]C(O)[NH2]");
    const str_t ASPARAGINE_NT   = STR_LIT(AA_NT  "[CH2]C(O)[NH2]");
    const str_t ASPARAGINE_CT   = STR_LIT(AA_CT  "[CH2]C(O)[NH2]");

    const str_t ASPARTATE       = STR_LIT(AA_INT "[CH2]C(O)O");
    const str_t ASPARTATE_NT    = STR_LIT(AA_NT  "[CH2]C(O)O");
    const str_t ASPARTATE_CT    = STR_LIT(AA_CT  "[CH2]C(O)O");

    const str_t CYSTEINE        = STR_LIT(AA_INT "[CH2][SH]");
    const str_t CYSTEINE_NT     = STR_LIT(AA_NT  "[CH2][SH]");
    const str_t CYSTEINE_CT     = STR_LIT(AA_CT  "[CH2][SH]");

    const str_t GLUTAMIC_ACID    = STR_LIT(AA_INT "[CH2][CH2]C(O)O");
    const str_t GLUTAMIC_ACID_NT = STR_LIT(AA_NT  "[CH2][CH2]C(O)O");
    const str_t GLUTAMIC_ACID_CT = STR_LIT(AA_CT  "[CH2][CH2]C(O)O");

    const str_t GLUTAMINE       = STR_LIT(AA_INT "[CH2][CH2]C(O)[NH2]");
    const str_t GLUTAMINE_NT    = STR_LIT(AA_NT  "[CH2][CH2]C(O)[NH2]");
    const str_t GLUTAMINE_CT    = STR_LIT(AA_CT  "[CH2][CH2]C(O)[NH2]");

    const str_t HISTIDINE       = STR_LIT(AA_INT "[CH2]C1:N:[CH]:[NH]:[CH]:1");
    const str_t HISTIDINE_NT    = STR_LIT(AA_NT  "[CH2]C1:N:[CH]:[NH]:[CH]:1");
    const str_t HISTIDINE_CT    = STR_LIT(AA_CT  "[CH2]C1:N:[CH]:[NH]:[CH]:1");

    const str_t ISOLEUCINE      = STR_LIT(AA_INT "[CH]([CH3])[CH2][CH3]");
    const str_t ISOLEUCINE_NT   = STR_LIT(AA_NT  "[CH]([CH3])[CH2][CH3]");
    const str_t ISOLEUCINE_CT   = STR_LIT(AA_CT "[CH]([CH3])[CH2][CH3]");

    const str_t LEUCINE         = STR_LIT(AA_INT "[CH2][CH]([CH3])[CH3]");
    const str_t LEUCINE_NT      = STR_LIT(AA_NT  "[CH2][CH]([CH3])[CH3]");
    const str_t LEUCINE_CT      = STR_LIT(AA_CT  "[CH2][CH]([CH3])[CH3]");

    const str_t LYSINE          = STR_LIT(AA_INT "[CH2][CH2][CH2][CH2][NH3]");
    const str_t METHIONINE      = STR_LIT(AA_INT "[CH2][CH2]S[CH3]");
    const str_t PHENYLALANINE   = STR_LIT(AA_INT "[CH2]C1:[CH]:[CH]:[CH]:[CH]:[CH]:1");

    const str_t PROLINE         = STR_LIT("N1[C@@H](CO)[CH2][CH2][CH2]1");
    const str_t PROLINE_NT      = STR_LIT("N1[C@H](CO)[CH2][CH2][CH2]1");
    const str_t PROLINE_CT      = STR_LIT("N1[C@@H](C(O)O)[CH2][CH2][CH2]1");

    const str_t SERINE          = STR_LIT(AA_INT "[CH2][OH]");
    const str_t SERINE_NT       = STR_LIT(AA_NT  "[CH2][OH]");
    const str_t SERINE_CT       = STR_LIT(AA_CT  "[CH2][OH]");

    const str_t THREONINE       = STR_LIT(AA_INT "[CH]([CH3])[OH]");
    const str_t THREONINE_NT    = STR_LIT(AA_NT  "[CH]([CH3])[OH]");
    const str_t THREONINE_CT    = STR_LIT(AA_CT  "[CH]([CH3])[OH]");

    const str_t TRYPTOPHAN      = STR_LIT(AA_INT "[CH2]C1:[CH]:[NH]:C2:[CH]:[CH]:[CH]:[CH]:C:1:2");
    const str_t TRYPTOPHAN_NT   = STR_LIT(AA_NT "[CH2]C1:[CH]:[NH]:C2:[CH]:[CH]:[CH]:[CH]:C:1:2");
    const str_t TRYPTOPHAN_CT   = STR_LIT(AA_CT "[CH2]C1:[CH]:[NH]:C2:[CH]:[CH]:[CH]:[CH]:C:1:2");

    const str_t TYROSINE        = STR_LIT(AA_INT "[CH2]C1:[CH]:[CH]:C([OH]):[CH]:[CH]:1");
    const str_t TYROSINE_NT     = STR_LIT(AA_NT  "[CH2]C1:[CH]:[CH]:C([OH]):[CH]:[CH]:1");
    const str_t TYROSINE_CT     = STR_LIT(AA_CT  "[CH2]C1:[CH]:[CH]:C([OH]):[CH]:[CH]:1");

    const str_t VALINE          = STR_LIT(AA_INT "[CH]([CH3])[CH3]");
    const str_t VALINE_NT       = STR_LIT(AA_NT "[CH]([CH3])[CH3]");
    const str_t VALINE_CT       = STR_LIT(AA_CT "[CH]([CH3])[CH3]");
    
    const str_t SELENOCYSTEINE    = STR_LIT(AA_INT "[CH2][SeH]");
    const str_t SELENOCYSTEINE_NT = STR_LIT(AA_NT  "[CH2][SeH]");
    const str_t SELENOCYSTEINE_CT = STR_LIT(AA_CT  "[CH2][SeH]");

    const str_t PYRROLYSINE     = STR_LIT(AA_INT "[CH2][CH2][CH2][CH2][NH]C(=O)C1=N[CH][CH2][CH]1[CH3]");
    const str_t PYRROLYSINE_NT  = STR_LIT(AA_NT "[CH2][CH2][CH2][CH2][NH]C(=O)C1=N[CH][CH2][CH]1[CH3]");
    const str_t PYRROLYSINE_CT  = STR_LIT(AA_CT "[CH2][CH2][CH2][CH2][NH]C(=O)C1=N[CH][CH2][CH]1[CH3]");

#define NUCL_INT "P(O)(O)O[CH2][CH]1[CH](O)[CH2][CH](O1)N"
#define NUCL_5T  "[OH][CH2][CH]1[CH](O)[CH2][CH](O1)N"

    const str_t DA              = STR_LIT(NUCL_INT "2CNC3C(N)NCNC23");
    const str_t DC              = STR_LIT(NUCL_INT "2C(O)NC(N)CC2");
    const str_t DG              = STR_LIT(NUCL_INT "2C3NC(N)NC(O)C3NC2");
    const str_t DT              = STR_LIT(NUCL_INT "2C(O)NC(O)C(C2)C");
    const str_t DU              = STR_LIT(NUCL_INT "2C(O)NC(O)CC2");

    // Terminal variation of DNA molecules
    const str_t DA_alt          = STR_LIT(NUCL_5T "2CNC3C(N)NCNC23");
    const str_t DC_alt          = STR_LIT(NUCL_5T "2C(O)NC(N)CC2");
    const str_t DG_alt          = STR_LIT(NUCL_5T "2C3NC(N)NC(O)C3NC2");
    const str_t DT_alt          = STR_LIT(NUCL_5T "2C(O)NC(O)C(C2)C");
    const str_t DU_alt          = STR_LIT(NUCL_5T "2C(O)NC(O)CC2");

    typedef struct {
        str_t name;
        struct {
            str_t smiles;
            md_util_match_flags_t flags;
        } pattern[4];
    } res_t;

    typedef struct {
        str_t name;
        const md_system_t* sys;
    } test_sys_t;

#define PROT_FLAGS 0
#define NUCL_FLAGS MD_UTIL_MATCH_FLAGS_NO_CH

    res_t residues[] = {
        {STR_LIT("ALA"), {
            {ALANINE,                       PROT_FLAGS},
            {ALANINE_NT,                    PROT_FLAGS},
            {ALANINE_CT,                    PROT_FLAGS}},
        },
        {STR_LIT("ARG"), {
            {ARGININE,                      PROT_FLAGS},
            {ARGININE_NT,                   PROT_FLAGS},
            {ARGININE_CT,                   PROT_FLAGS}},
        },
        {STR_LIT("ASN"), {
            {ASPARAGINE,                    PROT_FLAGS},
            {ASPARAGINE_NT,                 PROT_FLAGS},
            {ASPARAGINE_CT,                 PROT_FLAGS}},
        },
        {STR_LIT("ASP"), {
            {ASPARTATE,                     PROT_FLAGS},
            {ASPARTATE_NT,                  PROT_FLAGS},
            {ASPARTATE_CT,                  PROT_FLAGS}},
        },
        {STR_LIT("CYS"), {
            {CYSTEINE,                      PROT_FLAGS},
            {CYSTEINE_NT,                   PROT_FLAGS},
            {CYSTEINE_CT,                   PROT_FLAGS}},
        },
        // Glycine is a b*tch. It has no sidechain. Therefore it will essentially match against every amino acid pattern
        // Thus, it has to be handled with extra care to avoid false positives.
        {STR_LIT("GLY"), {
            {GLYCINE,                       PROT_FLAGS | MD_UTIL_MATCH_FLAGS_STRICT_EDGE_COUNT},
            {GLYCINE_NT,                    PROT_FLAGS | MD_UTIL_MATCH_FLAGS_STRICT_EDGE_COUNT},
            {GLYCINE_CT,                    PROT_FLAGS | MD_UTIL_MATCH_FLAGS_STRICT_EDGE_COUNT}}
        },
        {STR_LIT("GLU"), {
            {GLUTAMIC_ACID,                 PROT_FLAGS},
            {GLUTAMIC_ACID_NT,              PROT_FLAGS},
            {GLUTAMIC_ACID_CT,              PROT_FLAGS}},
        },
        {STR_LIT("GLN"), {
            {GLUTAMINE,                     PROT_FLAGS},
            {GLUTAMINE_NT,                  PROT_FLAGS},
            {GLUTAMINE_CT,                  PROT_FLAGS}},
        },
        {STR_LIT("HIS"), {
            {HISTIDINE,                     PROT_FLAGS},
            {HISTIDINE_NT,                  PROT_FLAGS},
            {HISTIDINE_CT,                  PROT_FLAGS}},
        },
        {STR_LIT("ILE"), {
            {ISOLEUCINE,                    PROT_FLAGS},
            {ISOLEUCINE_NT,                 PROT_FLAGS},
            {ISOLEUCINE_CT,                 PROT_FLAGS}},
        },
        {STR_LIT("LEU"), {LEUCINE,          PROT_FLAGS}},
        {STR_LIT("LYS"), {LYSINE,           PROT_FLAGS}},
        {STR_LIT("MET"), {METHIONINE,       PROT_FLAGS}},
        {STR_LIT("PHE"), {PHENYLALANINE,    PROT_FLAGS}},
        {STR_LIT("PRO"), {PROLINE,          PROT_FLAGS}},
        {STR_LIT("SER"), {SERINE,           PROT_FLAGS}},
        {STR_LIT("THR"), {THREONINE,        PROT_FLAGS}},
        {STR_LIT("TRP"), {TRYPTOPHAN,       PROT_FLAGS}},
        {STR_LIT("TYR"), {TYROSINE,         PROT_FLAGS}},
        {STR_LIT("VAL"), {VALINE,           PROT_FLAGS}},

        {STR_LIT("SEC"), {SELENOCYSTEINE,   PROT_FLAGS}},
        {STR_LIT("PYR"), {PYRROLYSINE,      PROT_FLAGS}},

        {STR_LIT("DA"),  {{DA, NUCL_FLAGS}, {DA_alt, NUCL_FLAGS}}},
        {STR_LIT("DC"),  {{DC, NUCL_FLAGS}, {DC_alt, NUCL_FLAGS}}},
        {STR_LIT("DG"),  {{DG, NUCL_FLAGS}, {DG_alt, NUCL_FLAGS}}},
        {STR_LIT("DT"),  {{DT, NUCL_FLAGS}, {DT_alt, NUCL_FLAGS}}},
    };

    test_sys_t test_mols[] = {
        {STR_LIT("ALA"), &utest_fixture->mol_ala},
        {STR_LIT("AMYLOID PFTAA"), &utest_fixture->mol_centered},
        {STR_LIT("NUCLEOTIDES"), &utest_fixture->mol_nucleotides},
        {STR_LIT("DNA"), &utest_fixture->mol_dna},
        {STR_LIT("TRP"), &utest_fixture->mol_trp},
        {STR_LIT("1K4R"), &utest_fixture->mol_1k4r},
        {STR_LIT("2OR2"), &utest_fixture->mol_2or2},
        {STR_LIT("1FEZ"), &utest_fixture->mol_1fez},
        //{STR_LIT("ASPIRINE"), &utest_fixture->mol_aspirine},
    };

    md_timestamp_t t0 = md_time_current();
    for (const test_sys_t* test = test_mols; test != test_mols + ARRAY_SIZE(test_mols); ++test) {
        md_vm_arena_temp_t temp = md_vm_arena_temp_begin(alloc);
        for (const res_t* res = residues; res != residues + ARRAY_SIZE(residues); ++res) {
            const md_system_t* sys = test->sys;
            md_array(md_comp_idx_t) ref_list = 0;
            
            for (size_t i = 0; i < sys->comp.count; ++i) {
                if (str_eq(LBL_TO_STR(sys->comp.name[i]), res->name)) {
                    md_array_push(ref_list, i, alloc);
                }
            }
            size_t ref_count = md_array_size(ref_list);

            if (ref_count == 0) {
                continue;
            }

            bool has_h = false;
            bool has_ch = false;

            for (size_t ref_idx = 0; ref_idx < ref_count; ++ref_idx) {
                md_urange_t atom_range = md_comp_atom_range(&sys->comp, ref_list[ref_idx]);
                for (size_t i = atom_range.beg; i < atom_range.end; ++i) {
					md_atomic_number_t z_i = md_atom_atomic_number(&sys->atom, i);
                    if (z_i == MD_Z_H) {
                        has_h = true;
                    }

                    if (z_i == MD_Z_C) {
						const md_atomic_number_t z_j[] = { MD_Z_H };
                        if (md_atom_is_connected_to_atomic_numbers(&sys->atom, &sys->bond, i, z_j, ARRAY_SIZE(z_j))) {
                            has_ch = true;
                            break;
						}
                    }
                }
                if (has_h && has_ch) {
                    break;
                }
            }

            md_util_match_flags_t extra_flags = 0;
            if (!has_h) {
                extra_flags |= MD_UTIL_MATCH_FLAGS_NO_H;
            }
            if (!has_ch) {
                extra_flags |= MD_UTIL_MATCH_FLAGS_NO_CH | MD_UTIL_MATCH_FLAGS_STRICT_EDGE_COUNT;
            }

            md_index_data_t match_result = { .alloc = alloc };

            for (size_t p_idx = 0; p_idx < ARRAY_SIZE(res->pattern); ++p_idx) {
                str_t smiles = res->pattern[p_idx].smiles;
                md_util_match_flags_t flags = res->pattern[p_idx].flags | extra_flags;

                if (str_empty(smiles)) {
                    break;
                }

                md_util_match_smiles(&match_result, smiles, MD_UTIL_MATCH_MODE_FIRST, MD_UTIL_MATCH_LEVEL_COMPONENT, flags, sys, alloc);
                
                size_t tot_count = md_index_data_num_ranges(&match_result);
                if (tot_count >= ref_count) {
                    break;
                }
            }
            size_t match_count = md_index_data_num_ranges(&match_result);

            EXPECT_EQ(ref_count, match_count);
            if (ref_count != match_count) {

                md_hashset_t match_set = {.allocator = alloc};
                for (size_t i = 0; i < match_count; ++i) {
                    int* atom_idx = md_index_range_beg(&match_result, i);
                    md_comp_idx_t res_idx = md_comp_find_by_atom_idx(&sys->comp, atom_idx[0]);
                    md_hashset_add(&match_set, res_idx);
                }
                
                printf("Mismatch in dataset '" STR_FMT "' for residues with name '" STR_FMT "': resname_count: %zu, match_count: %zu\n", STR_ARG(test->name), STR_ARG(res->name), ref_count, match_count);
                printf("The mismatcing residues indices are:\n");
                for (size_t i = 0; i < ref_count; ++i) {
                    if (!md_hashset_get(&match_set, ref_list[i])) {
                        printf("%i ", ref_list[i] + 1);
                    }
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


UTEST(integration, multiple_formats_tryptophan) {
    md_allocator_i* alloc = md_get_heap_allocator();
    
    // Load same tryptophan molecule from different formats
    md_system_t mol_pdb = {0};
    md_system_t mol_gro = {0};
    md_system_t mol_xyz = {0};
    
    bool result_pdb = md_pdb_molecule_api()->init_from_file(&mol_pdb, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan.pdb"), NULL, alloc);
    bool result_gro = md_gro_molecule_api()->init_from_file(&mol_gro, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.gro"), NULL, alloc);
    bool result_xyz = md_xyz_molecule_api()->init_from_file(&mol_xyz, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan.xyz"), NULL, alloc);
    
    ASSERT_TRUE(result_pdb);
    ASSERT_TRUE(result_gro);
    ASSERT_TRUE(result_xyz);
    
    // Test that all files loaded successfully and have reasonable atom counts
    EXPECT_EQ(mol_pdb.atom.count, 28);
    EXPECT_EQ(mol_xyz.atom.count, 27);  // XYZ has 27 atoms (different representation)
    // GRO might have a different structure (trajectory frame), so we don't compare directly
    EXPECT_GT(mol_gro.atom.count, 0);
    
    // All should have some atoms representing tryptophan
    EXPECT_GT(mol_pdb.atom.count, 20);
    EXPECT_GT(mol_xyz.atom.count, 20);
    EXPECT_GT(mol_gro.atom.count, 20);
    
    md_system_free(&mol_pdb, alloc);
    md_system_free(&mol_gro, alloc);
    md_system_free(&mol_xyz, alloc);
}
