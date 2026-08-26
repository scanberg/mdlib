#include "utest.h"
#include <string.h>
#include <math.h>

#include <md_pdb.h>
#include <md_gro.h>
#include <md_xyz.h>
#include <md_mmcif.h>
#include <md_trajectory.h>
#include <md_system.h>
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
};

UTEST_F_SETUP(util) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    utest_fixture->alloc = alloc;

    utest_fixture->mol_ala.alloc = alloc;
    md_system_state_t mol_ala_state = { .alloc = alloc };
    md_pdb_system_init_from_file(&utest_fixture->mol_ala, &mol_ala_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE);
    md_util_system_infer(&utest_fixture->mol_ala, &mol_ala_state, MD_UTIL_INFER_ALL);

    utest_fixture->mol_pftaa.alloc = alloc;
    md_system_state_t mol_pftaa_state = { .alloc = alloc };
    md_gro_system_init_from_file(&utest_fixture->mol_pftaa, &mol_pftaa_state, STR_LIT(MD_UNITTEST_DATA_DIR "/pftaa.gro"));
    md_util_system_infer(&utest_fixture->mol_pftaa, &mol_pftaa_state, MD_UTIL_INFER_ALL);

    utest_fixture->mol_nucleotides.alloc = alloc;
    md_system_state_t mol_nucleotides_state = { .alloc = alloc };
    md_gro_system_init_from_file(&utest_fixture->mol_nucleotides, &mol_nucleotides_state, STR_LIT(MD_UNITTEST_DATA_DIR "/nucleotides.gro"));
    md_util_system_infer(&utest_fixture->mol_nucleotides, &mol_nucleotides_state, MD_UTIL_INFER_ALL);

    utest_fixture->mol_centered.alloc = alloc;
    md_system_state_t mol_centered_state = { .alloc = alloc };
    md_gro_system_init_from_file(&utest_fixture->mol_centered, &mol_centered_state, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro"));
    md_util_system_infer(&utest_fixture->mol_centered, &mol_centered_state, MD_UTIL_INFER_ALL);

    utest_fixture->mol_dna.alloc = alloc;
    md_system_state_t mol_dna_state = { .alloc = alloc };
    md_gro_system_init_from_file(&utest_fixture->mol_dna, &mol_dna_state, STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro"));
    md_util_system_infer(&utest_fixture->mol_dna, &mol_dna_state, MD_UTIL_INFER_ALL);

    utest_fixture->mol_trp.alloc = alloc;
    md_system_state_t mol_trp_state = { .alloc = alloc };
    md_gro_system_init_from_file(&utest_fixture->mol_trp, &mol_trp_state, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.gro"));
    md_util_system_infer(&utest_fixture->mol_trp, &mol_trp_state, MD_UTIL_INFER_ALL);

    utest_fixture->mol_aspirine.alloc = alloc;
    md_system_state_t mol_aspirine_state = { .alloc = alloc };
    md_gro_system_init_from_file(&utest_fixture->mol_aspirine, &mol_aspirine_state, STR_LIT(MD_UNITTEST_DATA_DIR "/inside-md-pullout.gro"));
    md_util_system_infer(&utest_fixture->mol_aspirine, &mol_aspirine_state, MD_UTIL_INFER_ALL);

    utest_fixture->mol_1fez.alloc = alloc;
    md_system_state_t mol_1fez_state = { .alloc = alloc };
    md_mmcif_system_init_from_file(&utest_fixture->mol_1fez, &mol_1fez_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif"));
    md_util_system_infer(&utest_fixture->mol_1fez, &mol_1fez_state, MD_UTIL_INFER_ALL);

    utest_fixture->mol_2or2.alloc = alloc;
    md_system_state_t mol_2or2_state = { .alloc = alloc };
    md_mmcif_system_init_from_file(&utest_fixture->mol_2or2, &mol_2or2_state, STR_LIT(MD_UNITTEST_DATA_DIR "/2or2.cif"));
    md_util_system_infer(&utest_fixture->mol_2or2, &mol_2or2_state, MD_UTIL_INFER_ALL);

    utest_fixture->mol_1k4r.alloc = alloc;
    md_system_state_t mol_1k4r_state = { .alloc = alloc };
    md_pdb_system_init_from_file(&utest_fixture->mol_1k4r, &mol_1k4r_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), MD_PDB_OPTION_NONE);
    md_util_system_infer(&utest_fixture->mol_1k4r, &mol_1k4r_state, MD_UTIL_INFER_ALL);

    utest_fixture->mol_8g7u.alloc = alloc;
    md_system_state_t mol_8g7u_state = { .alloc = alloc };
    md_mmcif_system_init_from_file(&utest_fixture->mol_8g7u, &mol_8g7u_state, STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif"));
    md_util_system_infer(&utest_fixture->mol_8g7u, &mol_8g7u_state, MD_UTIL_INFER_ALL);
}

UTEST_F_TEARDOWN(util) {
    md_vm_arena_destroy(utest_fixture->alloc);
}

UTEST(util, hbonds) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* arena = md_temp_allocator(temp);

    md_system_t sys = { .alloc = arena };
    md_system_state_t sys_state = { .alloc = arena };
    md_gro_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro"));
    md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);

    md_hydrogen_bond_data_t hbond_data = {0};
    md_util_hydrogen_bond_init(&hbond_data, &sys, arena);

    EXPECT_LT(0, hbond_data.candidate.donor.count);
    EXPECT_LT(0, hbond_data.candidate.acceptor.count);

    md_util_hydrogen_bond_infer(&hbond_data, sys_state.x, sys_state.y, sys_state.z, &sys_state.unitcell, 3.0, 150.0);
    EXPECT_LT(0, hbond_data.num_bonds);

    md_temp_end(temp);
} 

UTEST_F(util, bonds) {
    EXPECT_EQ(152,      utest_fixture->mol_ala.bond.count);
    EXPECT_EQ(55,       utest_fixture->mol_pftaa.bond.count);
    EXPECT_EQ(40,       utest_fixture->mol_nucleotides.bond.count);
    EXPECT_EQ(163504,   utest_fixture->mol_centered.bond.count);
}

UTEST_F(util, inst) {
    EXPECT_EQ(1,        utest_fixture->mol_ala.instance.count);
	EXPECT_EQ(1,        utest_fixture->mol_pftaa.instance.count);
	EXPECT_EQ(2,        utest_fixture->mol_nucleotides.instance.count);
	EXPECT_EQ(253 + 61, utest_fixture->mol_centered.instance.count);

    const md_system_t* sys = &utest_fixture->mol_centered;
    ASSERT(sys->instance.count > 253);

    size_t ref_size = md_system_instance_atom_count(sys, 0);
    for (size_t i = 1; i < 253; ++i) {
        size_t size = md_system_instance_atom_count(sys, i);
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
    size_t num_structures_pftaa = md_structure_count(&utest_fixture->mol_pftaa.structure);
    EXPECT_EQ(1, num_structures_pftaa);
    size_t num_structures_nucleotides = md_structure_count(&utest_fixture->mol_nucleotides.structure);
    EXPECT_EQ(2, num_structures_nucleotides);
    size_t num_structures_ala = md_structure_count(&utest_fixture->mol_ala.structure);
	EXPECT_EQ(1, num_structures_ala);
	size_t num_structures_centered = md_structure_count(&utest_fixture->mol_centered.structure);
	EXPECT_EQ(253+61, num_structures_centered);
}

UTEST_F(util, rmsd) {
    md_allocator_i* alloc = utest_fixture->alloc;
    md_system_t* mol = &utest_fixture->mol_ala;
    md_trajectory_i* traj = utest_fixture->mol_ala.trajectory;
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

    md_trajectory_load_frame(traj, 0, &(md_system_state_t){0, x[0], y[0], z[0], {0}});
    md_trajectory_load_frame(traj, 1, &(md_system_state_t){0, x[1], y[1], z[1], {0}});

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
    md_unitcell_t cell = md_unitcell_from_extent(5,0,0);
    {
        const vec4_t xyzw[] = {
            {1,0,0,1},
            {2,0,0,1},
            {3,0,0,1},
            {4,0,0,1},
        };

        vec3_t com = md_util_com_compute_vec4(xyzw, 0, ARRAY_SIZE(xyzw), &cell);
        EXPECT_NEAR(2.5f, com.x, 1.0E-5F);
        EXPECT_EQ(0, com.y);
        EXPECT_EQ(0, com.z);
    }
    
    {
        const vec4_t xyzw[] = {
            {0,0,0,1},
            {5,0,0,1},
        };

        vec3_t com = md_util_com_compute_vec4(xyzw, 0, ARRAY_SIZE(xyzw), &cell);
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

        vec3_t com = md_util_com_compute_vec4(xyzw, 0, ARRAY_SIZE(xyzw), &cell);
        com = vec3_deperiodize_ortho(com, vec3_mul1(pbc_ext, 0.5f), pbc_ext);
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

        vec3_t com0 = md_util_com_compute_vec4(pos0, 0, ARRAY_SIZE(pos0), &cell);
        vec3_t com1 = md_util_com_compute_vec4(pos1, 0, ARRAY_SIZE(pos1), &cell);
        vec3_t com2 = md_util_com_compute_vec4(pos2, 0, ARRAY_SIZE(pos2), &cell);

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

    num_structures = md_structure_count(&utest_fixture->mol_nucleotides.structure);
    EXPECT_EQ(num_structures, 2);

    num_structures = md_structure_count(&utest_fixture->mol_ala.structure);
	EXPECT_EQ(num_structures, 1);

	num_structures = md_structure_count(&utest_fixture->mol_pftaa.structure);
	EXPECT_EQ(num_structures, 1);

	num_structures = md_structure_count(&utest_fixture->mol_centered.structure);
    EXPECT_EQ(num_structures, 253 + 61); // Chains + PFTAA
}

// Builds a system of 'num_chains' disconnected linear chains of 'chain_len' atoms each.
// Atom c*chain_len + i is bonded to its two neighbours within the chain and nothing else.
static void make_chain_system(md_system_t* sys, md_allocator_i* alloc, int num_chains, int chain_len) {
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

// The structure hierarchy must be rooted at the graph center and every atom must appear after the
// atom it was reached from. md_util_unwrap_structure depends on
// that ordering, and atom_slot must be a consistent inverse of atom_idx.
UTEST(util, structure_hierarchy) {
    md_allocator_i* alloc = md_get_heap_allocator();

    {
        // A chain of nine has a single center, the middle atom
        md_system_t sys = {0};
        make_chain_system(&sys, alloc, 1, 9);
        ASSERT_TRUE(md_util_system_infer_structures(&sys));

        ASSERT_EQ(md_structure_count(&sys.structure), 1);
        EXPECT_EQ(sys.structure.atom_idx[0], 4);
        EXPECT_EQ(sys.structure.parent_idx[0], sys.structure.atom_idx[0]); // the root is its own parent

        for (size_t slot = 0; slot < 9; ++slot) {
            const int32_t atom   = sys.structure.atom_idx[slot];
            const int32_t parent = sys.structure.parent_idx[slot];
            EXPECT_EQ(md_structure_atom_slot(&sys.structure, atom), (int32_t)slot);
            EXPECT_EQ(md_structure_atom_parent(&sys.structure, atom), parent);
            if (parent != atom) {
                // slot(parent) < slot(child), the invariant a single forward pass relies on
                EXPECT_LT(md_structure_atom_slot(&sys.structure, parent), (int32_t)slot);
            }
        }

        md_system_free(&sys);
    }

    {
        // Separate components are rooted independently and stored back to back
        md_system_t sys = {0};
        make_chain_system(&sys, alloc, 2, 5);
        ASSERT_TRUE(md_util_system_infer_structures(&sys));

        ASSERT_EQ(md_structure_count(&sys.structure), 2);
        EXPECT_EQ(sys.structure.offset[0], 0u);
        EXPECT_EQ(sys.structure.offset[1], 5u);
        EXPECT_EQ(sys.structure.offset[2], 10u);
        EXPECT_EQ(sys.structure.atom_idx[0], 2);
        EXPECT_EQ(sys.structure.atom_idx[5], 7);

        md_system_free(&sys);
    }

    {
        // A lone atom and a bonded pair must not trip the double sweep used to locate the center
        md_system_t sys = {0};
        sys.alloc = alloc;
        sys.atom.count = 3;
        md_atom_pair_t pair = {1, 2};
        md_array_push(sys.bond.pairs, pair, alloc);
        sys.bond.count = 1;
        md_bond_build_connectivity(&sys.bond, 3, alloc);
        ASSERT_TRUE(md_util_system_infer_structures(&sys));

        ASSERT_EQ(md_structure_count(&sys.structure), 2);
        for (size_t slot = 0; slot < 3; ++slot) {
            EXPECT_EQ(md_structure_atom_slot(&sys.structure, sys.structure.atom_idx[slot]), (int32_t)slot);
        }

        md_system_free(&sys);
    }
}

// A chain of nine atoms spaced 2A along x inside a 10A box, so the wrapped input folds back on
// itself twice. md_util_unwrap_structure must recover a straight line.
//
// @NOTE: this is the regression test for two bugs that were live together here. min_image_ortho was
// being fed a half extent where it wants a reciprocal one, and the result was written as x + dx
// instead of parent + dx, which is wrong for every non root atom even with no image correction at
// all. Either one alone destroys the spacing checked below.
UTEST(util, unwrap_structure_ortho) {
    md_allocator_i* alloc = md_get_heap_allocator();

    md_system_t sys = {0};
    make_chain_system(&sys, alloc, 1, 9);
    ASSERT_TRUE(md_util_system_infer_structures(&sys));

    float x[9], y[9] = {0}, z[9] = {0};
    for (int i = 0; i < 9; ++i) {
        x[i] = fmodf(1.0f + 2.0f * i, 10.0f);
    }
    md_system_state_t state = { .num_atoms = 9, .x = x, .y = y, .z = z, .unitcell = md_unitcell_from_extent(10, 10, 10) };

    md_structure_t structure = {0};
    ASSERT_TRUE(md_structure_extract(&structure, &sys.structure, 0));
    md_util_unwrap_structure(&state, &structure);

    for (int i = 0; i < 8; ++i) {
        EXPECT_NEAR(x[i + 1] - x[i], 2.0f, 1.0e-4f);
        EXPECT_NEAR(y[i], 0.0f, 1.0e-4f);
        EXPECT_NEAR(z[i], 0.0f, 1.0e-4f);
    }

    md_system_free(&sys);
}

// The same chain folded back by whole lattice vectors of a sheared cell. Unwrapping must recover the
// original straight line up to a single rigid lattice translation of the whole structure.
UTEST(util, unwrap_structure_triclinic) {
    md_allocator_i* alloc = md_get_heap_allocator();

    md_system_t sys = {0};
    make_chain_system(&sys, alloc, 1, 9);
    ASSERT_TRUE(md_util_system_infer_structures(&sys));

    md_unitcell_t cell = md_unitcell_from_basis_parameters(10, 10, 10, 2, 1, 3);
    float A[3][3] = {0};
    md_unitcell_A_extract_float(A, &cell);

    float tx[9], ty[9], tz[9];
    float x[9], y[9], z[9];
    for (int i = 0; i < 9; ++i) {
        tx[i] = 1.0f + 2.00f * i;
        ty[i] = 1.0f + 0.50f * i;
        tz[i] = 1.0f + 0.25f * i;

        const int na = (i >= 5) ? -1 : 0;   // displace the tail by one lattice vector along a
        const int nb = (i >= 7) ? -1 : 0;   // and the last two also along b

        x[i] = tx[i] + na * A[0][0] + nb * A[1][0];
        y[i] = ty[i] + na * A[0][1] + nb * A[1][1];
        z[i] = tz[i] + na * A[0][2] + nb * A[1][2];
    }
    md_system_state_t state = { .num_atoms = 9, .x = x, .y = y, .z = z, .unitcell = cell };

    md_structure_t structure = {0};
    ASSERT_TRUE(md_structure_extract(&structure, &sys.structure, 0));
    md_util_unwrap_structure(&state, &structure);

    const float ox = x[0] - tx[0];
    const float oy = y[0] - ty[0];
    const float oz = z[0] - tz[0];

    for (int i = 0; i < 9; ++i) {
        EXPECT_NEAR(x[i], tx[i] + ox, 1.0e-3f);
        EXPECT_NEAR(y[i], ty[i] + oy, 1.0e-3f);
        EXPECT_NEAR(z[i], tz[i] + oz, 1.0e-3f);
    }

    md_system_free(&sys);
}

// md_util_min_image_vec3 / _vec4 reduce a separation vector to its shortest periodic image, so no
// component may exceed half the box and the result must differ from the input by whole box lengths.
UTEST(util, min_image) {
    const float e = 10.0f;
    md_unitcell_t cell = md_unitcell_from_extent(e, e, e);

    vec3_t dx3[4] = {
        { 1.0f, -2.0f,  3.0f},      // already minimal, must be left alone
        {12.0f,  0.0f,  0.0f},
        {-27.0f, 41.0f, -8.0f},
        { 4.9f, -4.9f,  4.9f},      // just inside the half box
    };
    const vec3_t in3[4] = { dx3[0], dx3[1], dx3[2], dx3[3] };

    md_util_min_image_vec3(dx3, 4, &cell);

    for (int i = 0; i < 4; ++i) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_LE(fabsf(dx3[i].elem[c]), 0.5f * e + 1.0e-4f);
            const float shift = in3[i].elem[c] - dx3[i].elem[c];
            EXPECT_NEAR(shift, e * nearbyintf(shift / e), 1.0e-3f);
        }
    }

    EXPECT_NEAR(dx3[0].x,  1.0f, 1.0e-4f);
    EXPECT_NEAR(dx3[0].y, -2.0f, 1.0e-4f);
    EXPECT_NEAR(dx3[0].z,  3.0f, 1.0e-4f);
    EXPECT_NEAR(dx3[1].x,  2.0f, 1.0e-4f);
    EXPECT_NEAR(dx3[3].x,  4.9f, 1.0e-4f);

    vec4_t dx4[2] = { {12.0f, 0.0f, 0.0f, 7.0f}, {-27.0f, 41.0f, -8.0f, 8.0f} };
    md_util_min_image_vec4(dx4, 2, &cell);
    EXPECT_NEAR(dx4[0].x, 2.0f, 1.0e-4f);
    EXPECT_NEAR(dx4[1].x, 3.0f, 1.0e-4f);
    EXPECT_NEAR(dx4[1].y, 1.0f, 1.0e-4f);
    EXPECT_NEAR(dx4[1].z, 2.0f, 1.0e-4f);
    EXPECT_NEAR(dx4[0].w, 7.0f, 1.0e-4f);   // w must be left untouched
    EXPECT_NEAR(dx4[1].w, 8.0f, 1.0e-4f);
}


// ---- periodic image selection -------------------------------------------------------------

#define PBC_NP 12

// A lumpy, chiral-ish shape so no two points are interchangeable and the fit is well determined.
static void pbc_make_shape(vec4_t out[PBC_NP], float radius) {
    for (int i = 0; i < PBC_NP; ++i) {
        const float t = (float)i / PBC_NP * 6.2831853f;
        const float u = (float)(i % 5) / 5.0f;
        out[i] = vec4_set(radius * cosf(t) * (0.4f + u),
                          radius * sinf(t) * (0.5f + 0.5f * u),
                          radius * (u - 0.5f),
                          1.0f + 0.05f * i);
    }
}

static mat3_t pbc_rot(float angle) {
    const float c = cosf(angle), s = sinf(angle);
    const vec3_t k = vec3_normalize(vec3_set(0.3f, 0.6f, 0.74f));
    mat3_t m;
    m.elem[0][0]=c+k.x*k.x*(1-c);     m.elem[1][0]=k.x*k.y*(1-c)-k.z*s; m.elem[2][0]=k.x*k.z*(1-c)+k.y*s;
    m.elem[0][1]=k.y*k.x*(1-c)+k.z*s; m.elem[1][1]=c+k.y*k.y*(1-c);     m.elem[2][1]=k.y*k.z*(1-c)-k.x*s;
    m.elem[0][2]=k.z*k.x*(1-c)-k.y*s; m.elem[1][2]=k.z*k.y*(1-c)+k.x*s; m.elem[2][2]=c+k.z*k.z*(1-c);
    return m;
}

static float pbc_wrap(float v, float L) { v = fmodf(v, L); return v < 0.0f ? v + L : v; }

// Angle between R and the inverse of R_true, in degrees. Zero when R undoes R_true.
static float pbc_angle_err(mat3_t R, mat3_t R_true) {
    const mat3_t E = mat3_mul(R, R_true);
    const float tr = E.elem[0][0] + E.elem[1][1] + E.elem[2][2];
    return acosf(CLAMP((tr - 1.0f) * 0.5f, -1.0f, 1.0f)) * 180.0f / 3.14159265f;
}

// A set broken across a boundary must come back whole, with a centre that does not depend on which
// images the input happened to arrive in.
UTEST(util, pbc_deperiodize_self) {
    const float L = 20.0f;
    md_unitcell_t cell = md_unitcell_from_extent(L, L, L);

    vec4_t shape[PBC_NP];
    pbc_make_shape(shape, 3.0f);

    // Place it straddling the origin corner, then wrap every point independently
    vec4_t pts[PBC_NP];
    for (int i = 0; i < PBC_NP; ++i) {
        const vec3_t v = vec3_from_vec4(shape[i]);
        pts[i] = vec4_set(pbc_wrap(v.x, L), pbc_wrap(v.y, L), pbc_wrap(v.z, L), shape[i].w);
    }

    vec3_t com = {0};
    ASSERT_TRUE(md_util_deperiodize_self_vec4(pts, PBC_NP, &cell, &com));

    // Every pairwise separation must match the original shape: the set is whole again
    for (int i = 0; i < PBC_NP; ++i) {
        for (int j = i + 1; j < PBC_NP; ++j) {
            const float got = vec3_distance(vec3_from_vec4(pts[i]), vec3_from_vec4(pts[j]));
            const float ref = vec3_distance(vec3_from_vec4(shape[i]), vec3_from_vec4(shape[j]));
            EXPECT_NEAR(got, ref, 1.0e-3f);
        }
    }

    // The centre must sit inside the reassembled set, not at some average of scattered images
    const vec3_t shape_com = md_util_com_compute_vec4(shape, 0, PBC_NP, 0);
    for (int i = 0; i < PBC_NP; ++i) {
        const vec3_t d = vec3_sub(vec3_from_vec4(pts[i]), com);
        const vec3_t r = vec3_sub(vec3_from_vec4(shape[i]), shape_com);
        EXPECT_NEAR(d.x, r.x, 1.0e-3f);
        EXPECT_NEAR(d.y, r.y, 1.0e-3f);
        EXPECT_NEAR(d.z, r.z, 1.0e-3f);
    }

    // Already consistent input must be left alone
    vec4_t again[PBC_NP];
    MEMCPY(again, pts, sizeof(pts));
    vec3_t com2 = {0};
    ASSERT_TRUE(md_util_deperiodize_self_vec4(again, PBC_NP, &cell, &com2));
    for (int i = 0; i < PBC_NP; ++i) {
        EXPECT_NEAR(again[i].x, pts[i].x, 1.0e-4f);
        EXPECT_NEAR(again[i].w, pts[i].w, 1.0e-6f);  // weights carried through untouched
    }
}

// A set that arrives whole, in an image other than the reference one, must be left in THAT image.
//
// The alternation seeds from the circular mean, which always lands in the reference cell, so
// without an explicit correction the whole set is quietly translated into image (0,0,0). That is
// invisible to a caller reading relative geometry and a whole cell vector wrong for one that
// combines the reported centre with coordinates it never passed in - which is how viamd's recenter
// came to place a structure in the middle of the WRONG image on nojump trajectories.
UTEST(util, pbc_deperiodize_self_keeps_image) {
    const float L = 20.0f;
    md_unitcell_t cell = md_unitcell_from_extent(L, L, L);

    vec4_t shape[PBC_NP];
    pbc_make_shape(shape, 3.0f);
    const vec3_t shape_com = md_util_com_compute_vec4(shape, 0, PBC_NP, 0);

    for (int nx = -2; nx <= 2; ++nx) {
        for (int nz = -1; nz <= 1; ++nz) {
            const vec3_t off = vec3_set(nx * L, 10.0f, 10.0f + nz * L);

            vec4_t pts[PBC_NP];
            for (int i = 0; i < PBC_NP; ++i) {
                pts[i] = vec4_from_vec3(vec3_add(vec3_from_vec4(shape[i]), off), shape[i].w);
            }

            vec3_t com = {0};
            ASSERT_TRUE(md_util_deperiodize_self_vec4(pts, PBC_NP, &cell, &com));

            // Nothing was broken to begin with, so nothing may move
            for (int i = 0; i < PBC_NP; ++i) {
                const vec3_t want = vec3_add(vec3_from_vec4(shape[i]), off);
                EXPECT_NEAR(pts[i].x, want.x, 1.0e-3f);
                EXPECT_NEAR(pts[i].y, want.y, 1.0e-3f);
                EXPECT_NEAR(pts[i].z, want.z, 1.0e-3f);
            }

            // and the centre must be reported in the image the points actually occupy
            const vec3_t want_com = vec3_add(shape_com, off);
            EXPECT_NEAR(com.x, want_com.x, 1.0e-3f);
            EXPECT_NEAR(com.y, want_com.y, 1.0e-3f);
            EXPECT_NEAR(com.z, want_com.z, 1.0e-3f);
        }
    }
}

// Same contract for the joint fit: out_com and the placed points come back in the target's image,
// so a caller can build a transform against coordinates it did not hand over. With a rotation in
// play the old behaviour was not merely one cell vector off - it was off by R times a cell vector,
// which is not a periodic image at all.
UTEST(util, pbc_optimal_rotation_keeps_image) {
    const float L = 20.0f;
    md_unitcell_t cell = md_unitcell_from_extent(L, L, L);

    vec4_t ref[PBC_NP];
    pbc_make_shape(ref, 3.0f);
    const vec3_t ref_com = md_util_com_compute_vec4(ref, 0, PBC_NP, 0);

    const mat3_t R_true = pbc_rot(40.0f * 3.14159265f / 180.0f);

    for (int nx = -1; nx <= 2; ++nx) {
        const vec3_t c_true = vec3_set(9.0f + nx * L, 11.0f, 8.0f - nx * L);

        vec4_t trg[PBC_NP];
        for (int i = 0; i < PBC_NP; ++i) {
            const vec3_t v = vec3_add(mat3_mul_vec3(R_true, vec3_sub(vec3_from_vec4(ref[i]), ref_com)), c_true);
            trg[i] = vec4_from_vec3(v, ref[i].w);
        }

        mat3_t R = {0};
        vec3_t com = {0};
        vec4_t placed[PBC_NP];
        const float residual = md_util_optimal_rotation_pbc_vec4(&R, &com, placed, ref, ref_com, trg, PBC_NP, &cell);

        EXPECT_LT(residual, 1.0e-2f);
        EXPECT_NEAR(com.x, c_true.x, 1.0e-2f);
        EXPECT_NEAR(com.y, c_true.y, 1.0e-2f);
        EXPECT_NEAR(com.z, c_true.z, 1.0e-2f);

        // Untouched input, so the placed points must be exactly where they came in
        for (int i = 0; i < PBC_NP; ++i) {
            EXPECT_NEAR(placed[i].x, trg[i].x, 1.0e-2f);
            EXPECT_NEAR(placed[i].y, trg[i].y, 1.0e-2f);
            EXPECT_NEAR(placed[i].z, trg[i].z, 1.0e-2f);
        }

        // The recentering transform viamd builds must land the set on the cell centre, in image 0.
        const vec3_t centre = vec3_set(0.5f * L, 0.5f * L, 0.5f * L);
        vec3_t sum = vec3_zero();
        float  wsum = 0.0f;
        for (int i = 0; i < PBC_NP; ++i) {
            const vec3_t q = vec3_add(centre, mat3_mul_vec3(R, vec3_sub(vec3_from_vec4(trg[i]), com)));
            sum = vec3_add(sum, vec3_mul1(q, trg[i].w));
            wsum += trg[i].w;
        }
        sum = vec3_div1(sum, wsum);
        EXPECT_NEAR(sum.x, centre.x, 1.0e-2f);
        EXPECT_NEAR(sum.y, centre.y, 1.0e-2f);
        EXPECT_NEAR(sum.z, centre.z, 1.0e-2f);
    }
}

// The image the set is left in has to survive a basis with off diagonal terms, where the lattice
// vector being undone is not axis aligned.
UTEST(util, pbc_deperiodize_self_keeps_image_triclinic) {
    const float a = 47.8497f;
    const float c = 33.8349f;
    md_unitcell_t cell = md_unitcell_from_basis_parameters(a, a, c, 0.0, a * 0.5f, a * 0.5f);
    ASSERT_TRUE(md_unitcell_is_triclinic(&cell));

    mat3_t A = {0};
    md_unitcell_A_extract_float(A.elem, &cell);

    vec4_t shape[PBC_NP];
    pbc_make_shape(shape, 5.0f);
    const vec3_t shape_com = md_util_com_compute_vec4(shape, 0, PBC_NP, 0);
    const vec3_t base = mat3_mul_vec3(A, vec3_set(0.5f, 0.5f, 0.5f));

    for (int nx = -1; nx <= 1; ++nx) {
        for (int ny = -1; ny <= 1; ++ny) {
            const vec3_t off = vec3_add(base, mat3_mul_vec3(A, vec3_set((float)nx, (float)ny, 0.0f)));

            vec4_t pts[PBC_NP];
            for (int i = 0; i < PBC_NP; ++i) {
                pts[i] = vec4_from_vec3(vec3_add(vec3_from_vec4(shape[i]), off), shape[i].w);
            }

            vec3_t com = {0};
            ASSERT_TRUE(md_util_deperiodize_self_vec4(pts, PBC_NP, &cell, &com));

            const vec3_t want_com = vec3_add(shape_com, off);
            EXPECT_NEAR(com.x, want_com.x, 1.0e-2f);
            EXPECT_NEAR(com.y, want_com.y, 1.0e-2f);
            EXPECT_NEAR(com.z, want_com.z, 1.0e-2f);
        }
    }
}

// Every other test in this block uses a cube, which hides anything that only goes wrong when the
// basis has off diagonal terms. This one uses a rhombic dodecahedron - the shape GROMACS writes for
// -bt dodecahedron, and the shape that exposed the bug this test now guards.
//
// The circular mean has to be taken in FRACTIONAL space and carried back through the full basis.
// Handling each Cartesian axis on its own is only valid when the basis is diagonal; do it on a
// triclinic cell and the centre comes back displaced by a lattice vector, which then drags the
// whole deperiodize / align chain into the wrong periodic image.
UTEST(util, pbc_com_triclinic) {
    // a = b, c = a/sqrt(2), third vector leaning by half a cell in x and y
    const float a = 47.8497f;
    const float c = 33.8349f;
    md_unitcell_t cell = md_unitcell_from_basis_parameters(a, a, c, 0.0, a * 0.5f, a * 0.5f);
    ASSERT_TRUE(md_unitcell_is_triclinic(&cell));

    mat3_t A = {0};
    md_unitcell_A_extract_float(A.elem, &cell);
    mat3_t I = {0};
    md_unitcell_I_extract_float(I.elem, &cell);

    // A compact blob sitting well inside the cell, nowhere near a boundary.
    vec4_t pts[PBC_NP];
    const vec3_t centre = mat3_mul_vec3(A, vec3_set(0.41f, 0.63f, 0.48f));
    for (int i = 0; i < PBC_NP; ++i) {
        const float t = (float)i / PBC_NP * 6.2831853f;
        pts[i] = vec4_set(centre.x + 6.0f * cosf(t), centre.y + 6.0f * sinf(t), centre.z + 3.0f * cosf(t * 2.0f), 1.0f + (float)(i % 3));
    }

    // Nothing is wrapped, so the periodic mean must agree with the plain mean. It is the same set
    // of points either way.
    const vec3_t com_plain = md_util_com_compute_vec4(pts, 0, PBC_NP, 0);
    const vec3_t com_pbc   = md_util_com_compute_vec4(pts, 0, PBC_NP, &cell);
    EXPECT_NEAR(com_pbc.x, com_plain.x, 0.2f);
    EXPECT_NEAR(com_pbc.y, com_plain.y, 0.2f);
    EXPECT_NEAR(com_pbc.z, com_plain.z, 0.2f);

    // Same blob, every point folded into the primary cell the way a trajectory stores it. The
    // recovered centre must be the same one, in the same image - not a lattice vector away.
    vec4_t wrapped[PBC_NP];
    for (int i = 0; i < PBC_NP; ++i) {
        vec3_t f = mat3_mul_vec3(I, vec3_from_vec4(pts[i]));
        f.x -= floorf(f.x);
        f.y -= floorf(f.y);
        f.z -= floorf(f.z);
        wrapped[i] = vec4_from_vec3(mat3_mul_vec3(A, f), pts[i].w);
    }

    vec3_t com_wrapped = {0};
    ASSERT_TRUE(md_util_deperiodize_self_vec4(wrapped, PBC_NP, &cell, &com_wrapped));
    EXPECT_NEAR(com_wrapped.x, com_plain.x, 0.05f);
    EXPECT_NEAR(com_wrapped.y, com_plain.y, 0.05f);
    EXPECT_NEAR(com_wrapped.z, com_plain.z, 0.05f);

    // and the blob itself must be back in one piece, in that same image
    for (int i = 0; i < PBC_NP; ++i) {
        EXPECT_NEAR(wrapped[i].x, pts[i].x, 0.05f);
        EXPECT_NEAR(wrapped[i].y, pts[i].y, 0.05f);
        EXPECT_NEAR(wrapped[i].z, pts[i].z, 0.05f);
    }
}

// ------------------------------------------------------------------------------------------------
// Periodic centre of mass.
//
// The circular mean has to be taken in FRACTIONAL space and carried back out through the full
// basis. Resolving each Cartesian axis against itself alone is only valid when the basis is
// diagonal, so an ORTHORHOMBIC cell hides that mistake completely and a triclinic one does not.
// Everything below therefore runs against both a cube and a rhombic dodecahedron - the shape
// GROMACS writes for -bt dodecahedron, and the shape that surfaced this.
//
// md_util_com_compute dispatches to one of four internal variants depending on whether weights and
// indices are present, and each of those has an AVX512, an AVX2, an SSE2 and a scalar remainder
// path. The counts below straddle the block sizes (4 / 8 / 16) so that on any given build both the
// vector body and the remainder loop are exercised, and every variant is called for each count.

#define COM_MAX_PTS 253

static md_unitcell_t com_cell_ortho(void) {
    return md_unitcell_from_extent(47.8497, 47.8497, 33.8349);
}

static md_unitcell_t com_cell_triclinic(void) {
    const double a = 47.8497;   // a = b, c = a/sqrt(2), third vector leaning half a cell in x and y
    const double c = 33.8349;
    return md_unitcell_from_basis_parameters(a, a, c, 0.0, a * 0.5, a * 0.5);
}

// A compact blob sitting well inside the cell, nowhere near a boundary, so the periodic mean has an
// unambiguous answer: the ordinary arithmetic mean of the very same points.
static void com_make_blob(float* x, float* y, float* z, float* w, size_t count, const md_unitcell_t* cell) {
    mat3_t A = {0};
    md_unitcell_A_extract_float(A.elem, cell);
    const vec3_t centre = mat3_mul_vec3(A, vec3_set(0.41f, 0.63f, 0.48f));
    for (size_t i = 0; i < count; ++i) {
        const float t = (float)i * 0.7913f;
        x[i] = centre.x + 3.0f * cosf(t);
        y[i] = centre.y + 3.0f * sinf(t * 1.3f);
        z[i] = centre.z + 2.0f * cosf(t * 0.7f);
        w[i] = 1.0f + (float)(i % 7);
    }
}

static const size_t COM_COUNTS[] = { 1, 2, 3, 5, 7, 8, 15, 16, 17, 33, 61, 64, 253 };

// With nothing wrapped, the periodic mean is just the mean. All four variants have to say so, for
// both cell shapes, at every count - and the vec4 entry point has to agree with the float one.
UTEST(util, com_pbc_matches_plain_mean) {
    md_unitcell_t cells[2];
    cells[0] = com_cell_ortho();
    cells[1] = com_cell_triclinic();
    ASSERT_TRUE(md_unitcell_is_orthorhombic(&cells[0]));
    ASSERT_TRUE(md_unitcell_is_triclinic(&cells[1]));

    for (int ci = 0; ci < 2; ++ci) {
        for (size_t k = 0; k < ARRAY_SIZE(COM_COUNTS); ++k) {
            const size_t n = COM_COUNTS[k];
            float x[COM_MAX_PTS], y[COM_MAX_PTS], z[COM_MAX_PTS], w[COM_MAX_PTS];
            int32_t idx[COM_MAX_PTS];
            vec4_t xyzw[COM_MAX_PTS];
            com_make_blob(x, y, z, w, n, &cells[ci]);
            for (size_t i = 0; i < n; ++i) {
                idx[i]  = (int32_t)i;
                xyzw[i] = vec4_set(x[i], y[i], z[i], w[i]);
            }

            const vec3_t plain_u = md_util_com_compute(x, y, z, NULL, NULL, n, NULL);
            const vec3_t plain_w = md_util_com_compute(x, y, z, w,    NULL, n, NULL);

            // the four internal variants, in order: _com_pbc, _com_pbc_w, _com_pbc_i, _com_pbc_iw
            const vec3_t pbc    = md_util_com_compute(x, y, z, NULL, NULL, n, &cells[ci]);
            const vec3_t pbc_w  = md_util_com_compute(x, y, z, w,    NULL, n, &cells[ci]);
            const vec3_t pbc_i  = md_util_com_compute(x, y, z, NULL, idx,  n, &cells[ci]);
            const vec3_t pbc_iw = md_util_com_compute(x, y, z, w,    idx,  n, &cells[ci]);

            // and the vec4 entry point, which carries its weight in w
            const vec3_t pbc_v4 = md_util_com_compute_vec4(xyzw, NULL, n, &cells[ci]);
            const vec3_t pbc_v4i = md_util_com_compute_vec4(xyzw, idx, n, &cells[ci]);

            for (int e = 0; e < 3; ++e) {
                EXPECT_NEAR(pbc.elem[e],     plain_u.elem[e], 0.1f);
                EXPECT_NEAR(pbc_i.elem[e],   plain_u.elem[e], 0.1f);
                EXPECT_NEAR(pbc_w.elem[e],   plain_w.elem[e], 0.1f);
                EXPECT_NEAR(pbc_iw.elem[e],  plain_w.elem[e], 0.1f);
                EXPECT_NEAR(pbc_v4.elem[e],  plain_w.elem[e], 0.1f);
                EXPECT_NEAR(pbc_v4i.elem[e], plain_w.elem[e], 0.1f);
            }

            // indexed and contiguous walk the same points, so they must agree exactly, not merely
            // to within the tolerance above
            for (int e = 0; e < 3; ++e) {
                EXPECT_NEAR(pbc_i.elem[e],  pbc.elem[e],   1.0e-4f);
                EXPECT_NEAR(pbc_iw.elem[e], pbc_w.elem[e], 1.0e-4f);
            }
        }
    }
}

// The defining property, and the one that needs no tolerance argument: moving individual points by
// whole lattice vectors does not change the configuration, so it must not move the centre. A
// per axis circular mean is periodic in the cell EXTENT rather than in the lattice, which is the
// same thing for a cube and is not the same thing at all for a triclinic cell.
UTEST(util, com_pbc_invariant_under_lattice_shift) {
    md_unitcell_t cells[2];
    cells[0] = com_cell_ortho();
    cells[1] = com_cell_triclinic();

    for (int ci = 0; ci < 2; ++ci) {
        mat3_t A = {0};
        md_unitcell_A_extract_float(A.elem, &cells[ci]);

        for (size_t k = 0; k < ARRAY_SIZE(COM_COUNTS); ++k) {
            const size_t n = COM_COUNTS[k];
            float x[COM_MAX_PTS], y[COM_MAX_PTS], z[COM_MAX_PTS], w[COM_MAX_PTS];
            int32_t idx[COM_MAX_PTS];
            com_make_blob(x, y, z, w, n, &cells[ci]);
            for (size_t i = 0; i < n; ++i) idx[i] = (int32_t)i;

            const vec3_t before    = md_util_com_compute(x, y, z, NULL, NULL, n, &cells[ci]);
            const vec3_t before_w  = md_util_com_compute(x, y, z, w,    NULL, n, &cells[ci]);
            const vec3_t before_i  = md_util_com_compute(x, y, z, NULL, idx,  n, &cells[ci]);
            const vec3_t before_iw = md_util_com_compute(x, y, z, w,    idx,  n, &cells[ci]);

            // scatter the points across images - a different lattice vector for every third one
            for (size_t i = 0; i < n; ++i) {
                if (i % 3 != 0) continue;
                // NOTE: cast before subtracting - i is size_t and the subtraction would wrap
                const vec3_t nvec = vec3_set((float)((int)(i % 5) - 2), (float)((int)(i % 3) - 1), (float)((int)(i % 7) - 3));
                const vec3_t shift = mat3_mul_vec3(A, nvec);
                x[i] += shift.x;
                y[i] += shift.y;
                z[i] += shift.z;
            }

            const vec3_t after    = md_util_com_compute(x, y, z, NULL, NULL, n, &cells[ci]);
            const vec3_t after_w  = md_util_com_compute(x, y, z, w,    NULL, n, &cells[ci]);
            const vec3_t after_i  = md_util_com_compute(x, y, z, NULL, idx,  n, &cells[ci]);
            const vec3_t after_iw = md_util_com_compute(x, y, z, w,    idx,  n, &cells[ci]);

            for (int e = 0; e < 3; ++e) {
                EXPECT_NEAR(after.elem[e],    before.elem[e],    0.02f);
                EXPECT_NEAR(after_w.elem[e],  before_w.elem[e],  0.02f);
                EXPECT_NEAR(after_i.elem[e],  before_i.elem[e],  0.02f);
                EXPECT_NEAR(after_iw.elem[e], before_iw.elem[e], 0.02f);
            }
        }
    }
}

// The indexed variants must read through the index list and nothing else. Bury the real points in a
// larger array whose unselected slots hold coordinates far away, and the answer must not move.
UTEST(util, com_pbc_indices_select) {
    md_unitcell_t cells[2];
    cells[0] = com_cell_ortho();
    cells[1] = com_cell_triclinic();

    for (int ci = 0; ci < 2; ++ci) {
        const size_t n = 61;
        float xs[COM_MAX_PTS], ys[COM_MAX_PTS], zs[COM_MAX_PTS], ws[COM_MAX_PTS];
        com_make_blob(xs, ys, zs, ws, n, &cells[ci]);

        const vec3_t want   = md_util_com_compute(xs, ys, zs, NULL, NULL, n, &cells[ci]);
        const vec3_t want_w = md_util_com_compute(xs, ys, zs, ws,   NULL, n, &cells[ci]);

        // stride the real points through a 3x larger array, junk in between
        float bx[COM_MAX_PTS * 3], by[COM_MAX_PTS * 3], bz[COM_MAX_PTS * 3], bw[COM_MAX_PTS * 3];
        int32_t idx[COM_MAX_PTS];
        for (size_t i = 0; i < n * 3; ++i) {
            bx[i] = -931.0f + (float)i;
            by[i] =  757.0f - (float)i;
            bz[i] =  613.0f + (float)(i * 2);
            bw[i] =  99.0f;
        }
        for (size_t i = 0; i < n; ++i) {
            const size_t slot = i * 3 + 2;
            bx[slot] = xs[i]; by[slot] = ys[i]; bz[slot] = zs[i]; bw[slot] = ws[i];
            idx[i] = (int32_t)slot;
        }

        const vec3_t got   = md_util_com_compute(bx, by, bz, NULL, idx, n, &cells[ci]);
        const vec3_t got_w = md_util_com_compute(bx, by, bz, bw,   idx, n, &cells[ci]);

        for (int e = 0; e < 3; ++e) {
            EXPECT_NEAR(got.elem[e],   want.elem[e],   1.0e-4f);
            EXPECT_NEAR(got_w.elem[e], want_w.elem[e], 1.0e-4f);
        }
    }
}

// Weights have to actually weight. Two points, one ten times heavier, in a triclinic cell: the
// centre belongs near the heavy one, and nowhere near the midpoint.
UTEST(util, com_pbc_weights_apply) {
    md_unitcell_t cell = com_cell_triclinic();
    mat3_t A = {0};
    md_unitcell_A_extract_float(A.elem, &cell);

    const vec3_t p0 = mat3_mul_vec3(A, vec3_set(0.30f, 0.30f, 0.30f));
    const vec3_t p1 = mat3_mul_vec3(A, vec3_set(0.34f, 0.36f, 0.38f));

    float x[2] = { p0.x, p1.x };
    float y[2] = { p0.y, p1.y };
    float z[2] = { p0.z, p1.z };
    float w[2] = { 10.0f, 1.0f };

    const vec3_t expect = vec3_div1(vec3_add(vec3_mul1(p0, 10.0f), p1), 11.0f);
    const vec3_t got    = md_util_com_compute(x, y, z, w, NULL, 2, &cell);

    EXPECT_NEAR(got.x, expect.x, 0.05f);
    EXPECT_NEAR(got.y, expect.y, 0.05f);
    EXPECT_NEAR(got.z, expect.z, 0.05f);
}

// The rotation must be recovered from a target whose points arrive scattered across images, for any
// rotation, without ever consulting topology.
UTEST(util, pbc_optimal_rotation) {
    const float L = 20.0f;
    md_unitcell_t cell = md_unitcell_from_extent(L, L, L);

    for (int deg = 0; deg <= 180; deg += 20) {
        vec4_t ref[PBC_NP];
        pbc_make_shape(ref, 3.0f);
        const vec3_t ref_com = md_util_com_compute_vec4(ref, 0, PBC_NP, 0);

        const mat3_t R_true = pbc_rot(deg * 3.14159265f / 180.0f);
        const vec3_t c_true = { 3.0f, 17.0f, 9.0f };   // deliberately near two boundaries

        vec4_t trg[PBC_NP];
        for (int i = 0; i < PBC_NP; ++i) {
            const vec3_t v = vec3_add(mat3_mul_vec3(R_true, vec3_sub(vec3_from_vec4(ref[i]), ref_com)), c_true);
            trg[i] = vec4_set(pbc_wrap(v.x, L), pbc_wrap(v.y, L), pbc_wrap(v.z, L), ref[i].w);
        }

        mat3_t R = {0};
        vec3_t com = {0};
        vec4_t placed[PBC_NP];
        const float residual = md_util_optimal_rotation_pbc_vec4(&R, &com, placed, ref, ref_com, trg, PBC_NP, &cell);

        EXPECT_LT(residual, 1.0e-2f);
        EXPECT_LT(pbc_angle_err(R, R_true), 0.5f);

        // The placed points must satisfy the transform the call reported
        for (int i = 0; i < PBC_NP; ++i) {
            const vec3_t p = vec3_sub(vec3_from_vec4(ref[i]), ref_com);
            const vec3_t q = mat3_mul_vec3(R, vec3_sub(vec3_from_vec4(placed[i]), com));
            EXPECT_NEAR(q.x, p.x, 1.0e-2f);
            EXPECT_NEAR(q.y, p.y, 1.0e-2f);
            EXPECT_NEAR(q.z, p.z, 1.0e-2f);
            EXPECT_NEAR(placed[i].w, ref[i].w, 1.0e-6f);
        }
    }
}

// Beyond the regime where any per point image choice can work, the fit must SAY so rather than
// return a confident wrong answer. A set of radius 8.5 in a 20A cell flipped end for end is past it.
UTEST(util, pbc_optimal_rotation_reports_ambiguity) {
    const float L = 20.0f;
    md_unitcell_t cell = md_unitcell_from_extent(L, L, L);

    vec4_t ref[PBC_NP];
    pbc_make_shape(ref, 8.5f);
    const vec3_t ref_com = md_util_com_compute_vec4(ref, 0, PBC_NP, 0);

    const mat3_t R_true = pbc_rot(3.14159265f);
    const vec3_t c_true = { 3.0f, 17.0f, 9.0f };

    vec4_t trg[PBC_NP];
    for (int i = 0; i < PBC_NP; ++i) {
        const vec3_t v = vec3_add(mat3_mul_vec3(R_true, vec3_sub(vec3_from_vec4(ref[i]), ref_com)), c_true);
        trg[i] = vec4_set(pbc_wrap(v.x, L), pbc_wrap(v.y, L), pbc_wrap(v.z, L), ref[i].w);
    }

    mat3_t R = {0};
    vec3_t com = {0};
    const float residual = md_util_optimal_rotation_pbc_vec4(&R, &com, NULL, ref, ref_com, trg, PBC_NP, &cell);

    // The certificate is the whole point: a residual on the order of half a cell means do not trust it
    EXPECT_GT(residual, 0.25f * L);
}

// Degenerate inputs must not misbehave
UTEST(util, pbc_optimal_rotation_degenerate) {
    md_unitcell_t cell = md_unitcell_from_extent(20, 20, 20);
    vec4_t ref[PBC_NP];
    pbc_make_shape(ref, 3.0f);
    const vec3_t ref_com = md_util_com_compute_vec4(ref, 0, PBC_NP, 0);

    mat3_t R = {0};
    vec3_t com = {0};

    // Empty set
    EXPECT_NEAR(md_util_optimal_rotation_pbc_vec4(&R, &com, NULL, ref, ref_com, ref, 0, &cell), 0.0f, 1.0e-6f);

    // No cell at all: still a plain Kabsch fit, identity here since target == reference
    const float residual = md_util_optimal_rotation_pbc_vec4(&R, &com, NULL, ref, ref_com, ref, PBC_NP, NULL);
    EXPECT_LT(residual, 1.0e-3f);
    EXPECT_LT(pbc_angle_err(R, mat3_ident()), 0.5f);

    // A set with no coordinates to speak of
    vec3_t c = {0};
    EXPECT_TRUE(md_util_deperiodize_self_vec4(ref, 0, &cell, &c));
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
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp_scope);
	md_system_t sys = { .alloc = alloc };
	md_system_state_t sys_state = { .alloc = alloc };
	md_pdb_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/c60.pdb"), MD_PDB_OPTION_NONE);
	md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);

	EXPECT_EQ(sys.atom.count, 60);
	EXPECT_EQ(sys.bond.count, 90);

    const size_t num_rings = md_index_data_num_ranges(&sys.ring);
    EXPECT_EQ(num_rings, 32);

    const size_t num_structures = md_structure_count(&sys.structure);
    EXPECT_EQ(num_structures, 1);

    md_temp_end(temp_scope);
}

UTEST(util, rings_c720) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp_scope);
    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    md_xyz_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/c720.xyz"), MD_XYZ_OPTION_NONE);
    md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);

    EXPECT_EQ(sys.atom.count, 720);
    EXPECT_EQ(sys.bond.count, 1080);

    const size_t num_rings = md_index_data_num_ranges(&sys.ring);
    EXPECT_EQ(num_rings, 362);

    const size_t num_structures = md_structure_count(&sys.structure);
    EXPECT_EQ(num_structures, 1);

    md_temp_end(temp_scope);
}

UTEST(util, rings_14kr) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp_scope);
    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    md_pdb_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), MD_PDB_OPTION_NONE);
    md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);

    const size_t num_rings = md_index_data_num_ranges(&sys.ring);
    EXPECT_EQ(num_rings, 207);

    md_temp_end(temp_scope);
}

UTEST(util, rings_trytophan_pdb) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp_scope);
    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    md_pdb_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan.pdb"), MD_PDB_OPTION_NONE);
    md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);

    const size_t num_rings = md_index_data_num_ranges(&sys.ring);
    EXPECT_EQ(num_rings, 2);

    const size_t num_structures = md_structure_count(&sys.structure);
    EXPECT_EQ(num_structures, 1);

    md_temp_end(temp_scope);
}

UTEST(util, rings_trytophan_xyz) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp_scope);
    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    md_xyz_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan.xyz"), MD_XYZ_OPTION_NONE);
    md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);

    const size_t num_rings = md_index_data_num_ranges(&sys.ring);
    EXPECT_EQ(num_rings, 2);

    const size_t num_structures = md_structure_count(&sys.structure);
    EXPECT_EQ(num_structures, 1);

    md_temp_end(temp_scope);
}

UTEST(util, rings_full) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp_scope);
    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    md_xyz_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/full.xyz"), MD_XYZ_OPTION_NONE);
    md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);

    const size_t num_rings = md_index_data_num_ranges(&sys.ring);
    EXPECT_EQ(num_rings, 195);

    const size_t num_structures = md_structure_count(&sys.structure);
    EXPECT_EQ(num_structures, 1);

    md_temp_end(temp_scope);
}

UTEST(util, rings_ciprofloxacin) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp_scope);
    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    md_pdb_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/ciprofloxacin.pdb"), MD_PDB_OPTION_NONE);
    md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);

    const int64_t num_rings = md_index_data_num_ranges(&sys.ring);
    ASSERT_EQ(num_rings, 4);
    EXPECT_EQ(md_index_range_size(&sys.ring, 0), 3);
    EXPECT_EQ(md_index_range_size(&sys.ring, 1), 6);
    EXPECT_EQ(md_index_range_size(&sys.ring, 2), 6);
    EXPECT_EQ(md_index_range_size(&sys.ring, 3), 6);

    md_temp_end(temp_scope);
}

UTEST_F(util, structure_matching_amyloid_chain) {
    md_allocator_i* alloc = utest_fixture->alloc;
    md_system_t* sys = &utest_fixture->mol_centered;

#if 1
    {
        // Test for the chains
        md_structure_t ref_structure = {0};
        md_structure_extract(&ref_structure, &sys->structure, 0);

        const int32_t* ref_idx = ref_structure.atom_idx;
        size_t ref_len = ref_structure.count;

        // Prune Hydrogen
        if (true) {
            md_array(int) new_idx = 0;
            for (size_t i = 0; i < ref_structure.count; ++i) {
                if (md_atom_atomic_number(&sys->atom, ref_idx[i]) > 1) {
                    md_array_push(new_idx, ref_idx[i], alloc);
                }
            }
            ref_idx = new_idx;
            ref_len = md_array_size(new_idx);
        }

        md_timestamp_t t0 = md_time_now();
        md_index_data_t result = md_util_match_by_element(ref_idx, ref_len, MD_UTIL_MATCH_MODE_FIRST, MD_UTIL_MATCH_LEVEL_INSTANCE, sys, alloc);
        md_timestamp_t t1 = md_time_now();
        printf("time: %f ms\n", md_time_as_milliseconds(t1 - t0));
        size_t result_count = md_index_data_num_ranges(&result);
        printf("result count: %zu\n", result_count);
        EXPECT_EQ(result_count, 253);
    }
#endif
}

UTEST_F(util, structure_matching_PFTAA) {
    md_allocator_i* alloc = utest_fixture->alloc;
    md_system_t* sys = &utest_fixture->mol_centered;

#if 1
    {
        // Test for the PFTAAs
        const int ref_structure_idx = 253;
        md_structure_t ref_structure = { 0 };
        md_structure_extract(&ref_structure, &sys->structure, ref_structure_idx);
        const int*   ref_idx  = ref_structure.atom_idx;
        const size_t ref_size = ref_structure.count;

        md_timestamp_t t0 = md_time_now();
        md_index_data_t result = md_util_match_by_element(ref_idx, ref_size, MD_UTIL_MATCH_MODE_FIRST, MD_UTIL_MATCH_LEVEL_COMPONENT, sys, alloc);
        md_timestamp_t t1 = md_time_now();
        printf("time: %f ms\n", md_time_as_milliseconds(t1 - t0));
        size_t result_count = md_index_data_num_ranges(&result);
        printf("result count: %zu\n", result_count);
        EXPECT_EQ(result_count, 61);
    }
#endif
}

UTEST_F(util, structure_matching_PFTAA_ring) {
    md_allocator_i* alloc = utest_fixture->alloc;
    md_system_t* sys = &utest_fixture->mol_pftaa;
    {
#if 1
        // Rings, represents a ring within the molecule
        const int ref_idx[] = {19,20,21,22,24};
        const size_t ref_size = ARRAY_SIZE(ref_idx);

        md_index_data_t result = md_util_match_by_element(ref_idx, ref_size, MD_UTIL_MATCH_MODE_UNIQUE, MD_UTIL_MATCH_LEVEL_STRUCTURE, sys, alloc);
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
        md_timestamp_t t0 = md_time_now();
        md_index_data_t result = md_util_match_by_element(ref_idx, ref_len, MD_UTIL_MATCH_MODE_UNIQUE, MD_UTIL_MATCH_LEVEL_STRUCTURE, sys, alloc);
        md_timestamp_t t1 = md_time_now();
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
    const str_t ISOLEUCINE_CT   = STR_LIT(AA_CT  "[CH]([CH3])[CH2][CH3]");

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
        //{STR_LIT("2OR2"), &utest_fixture->mol_2or2},
        //{STR_LIT("1FEZ"), &utest_fixture->mol_1fez},
        //{STR_LIT("ASPIRINE"), &utest_fixture->mol_aspirine},
    };

    md_timestamp_t t0 = md_time_now();
    for (const test_sys_t* test = test_mols; test != test_mols + ARRAY_SIZE(test_mols); ++test) {
        md_temp_scope_t temp = md_temp_begin_in(alloc);
        for (const res_t* res = residues; res != residues + ARRAY_SIZE(residues); ++res) {
            const md_system_t* sys = test->sys;
            md_array(md_component_idx_t) ref_list = 0;
            
            for (size_t i = 0; i < sys->component.count; ++i) {
                if (str_eq(LBL_TO_STR(sys->component.name[i]), res->name)) {
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
                md_urange_t atom_range = md_component_atom_range(&sys->component, ref_list[ref_idx]);
                for (size_t i = atom_range.beg; i < atom_range.end; ++i) {
					md_atomic_number_t z_i = md_atom_atomic_number(&sys->atom, i);
                    if (z_i == MD_Z_H) {
                        has_h = true;
                    }

                    if (z_i == MD_Z_C) {
						const md_atomic_number_t z_j[] = { MD_Z_H };
                        if (md_atom_is_connected_to_atomic_numbers(&sys->atom, &sys->bond, i, z_j, ARRAY_SIZE(z_j))) {
                            has_ch = true;
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

			// Store all matching components in here, use a set to avoid duplicates
            md_hashset_t match_set = {.allocator = alloc};

            for (size_t p_idx = 0; p_idx < ARRAY_SIZE(res->pattern); ++p_idx) {
                str_t smiles = res->pattern[p_idx].smiles;
                md_util_match_flags_t flags = res->pattern[p_idx].flags | extra_flags;

                if (str_empty(smiles)) {
                    break;
                }

				md_index_data_t pattern_result = { .alloc = alloc };
                size_t matches = md_util_match_smiles(&pattern_result, smiles, MD_UTIL_MATCH_MODE_FIRST, MD_UTIL_MATCH_LEVEL_COMPONENT, flags, sys, alloc);

				size_t num_pattern_matches = md_index_data_num_ranges(&pattern_result);
                for (size_t i = 0; i < num_pattern_matches; ++i) {
					const int* atom_idx = md_index_range_beg(&pattern_result, i);
					md_component_idx_t res_idx = md_component_find_by_atom_idx(&sys->component, atom_idx[0]);
					md_hashset_add(&match_set, res_idx);
                }
            }
            size_t match_count = match_set.num_used;

            EXPECT_EQ(ref_count, match_count);
            if (ref_count != match_count) {
                
                printf("Mismatch in dataset '" STR_FMT "' for residues with name '" STR_FMT "': resname_count: %zu, match_count: %zu\n", STR_ARG(test->name), STR_ARG(res->name), ref_count, match_count);
                printf("Residues missed:\n");
                for (size_t i = 0; i < ref_count; ++i) {
                    if (!md_hashset_get(&match_set, ref_list[i])) {
                        printf("%i ", ref_list[i] + 1);
                    }
                }
            }
        }
        md_temp_end(temp);
    }
    md_timestamp_t t1 = md_time_now();
    printf("time: %f ms\n", md_time_as_milliseconds(t1-t0));

}

UTEST(util, parse_smiles) {
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
}

UTEST(util, radix_sort) {
    uint32_t arr[] = { 1, 278, 128312745, 4, 5, 0, 12382, 26, 12, 14, 7 };
    size_t len = ARRAY_SIZE(arr);

    uint32_t idx[16];

    md_util_sort_radix_uint32(idx, arr, len);
    for (size_t i = 0; i < len - 1; ++i) {
        EXPECT_LE(arr[idx[i]], arr[idx[i+1]]);
    }

    md_util_sort_radix_inplace_uint32(arr, len);

    for (size_t i = 0; i < len - 1; ++i) {
    	EXPECT_LE(arr[i], arr[i+1]);
    }

    md_temp_scope_t temp = md_temp_begin();

    size_t N = 10000000;
    uint32_t* values   = md_temp_alloc_array(temp, uint32_t, N);
    uint32_t* indices  = md_temp_alloc_array(temp, uint32_t, N);

    for (size_t i = 0; i < N; ++i) {
        values[i] = (uint32_t)rand() * rand();
    }

    md_timestamp_t t0, t1;

    t0 = md_time_now();
    md_util_sort_radix_uint32(indices, values, N);
    t1 = md_time_now();

    printf("Time for radix index sort: %.4f ms\n", md_time_as_milliseconds(t1 - t0));

    t0 = md_time_now();
    md_util_sort_radix_inplace_uint32(values, N);
    t1 = md_time_now();

    printf("Time for radix inplace sort: %.4f ms\n", md_time_as_milliseconds(t1 - t0));
    md_temp_end(temp);
}

static inline bool init_system(md_system_t* sys, md_system_state_t* sys_state, str_t path) {
    str_t ext;
    if (!extract_ext(&ext, path)) {
        return false;
    }

    if (str_eq_ignore_case(ext, STR_LIT("pdb"))) {
        return md_pdb_system_init_from_file(sys, sys_state, path, MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE);
    } else
    if (str_eq_ignore_case(ext, STR_LIT("gro"))) {
        return md_gro_system_init_from_file(sys, sys_state, path);
    } else
    if (str_eq_ignore_case(ext, STR_LIT("cif"))) {
        return md_mmcif_system_init_from_file(sys, sys_state, path);
    }

    return false;
}

UTEST(util, entity_instance) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp_scope);
    ASSERT(alloc);

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb")));
        EXPECT_GT(md_system_atom_count(&sys),   0);
        ASSERT_EQ(md_system_entity_count(&sys), 1);
        EXPECT_EQ(md_system_entity_flags(&sys, 0), MD_FLAG_POLYMER | MD_FLAG_POLYPEPTIDE);
        
        ASSERT_EQ(md_system_instance_count(&sys), 1);
        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 0), STR_LIT("A")));
        EXPECT_TRUE(str_empty(md_system_instance_auth_id(&sys, 0)));
    }

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb")));
        EXPECT_GT(md_system_atom_count(&sys),   0);
        ASSERT_EQ(md_system_entity_count(&sys), 1);
        EXPECT_EQ(md_system_entity_flags(&sys, 0), MD_FLAG_POLYMER | MD_FLAG_POLYPEPTIDE);
        
        ASSERT_EQ(md_system_instance_count(&sys), 3);
        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 0), STR_LIT("A")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 0), STR_LIT("A")));

        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 1), STR_LIT("B")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 1), STR_LIT("B")));

        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 2), STR_LIT("C")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 2), STR_LIT("C")));
    }

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1LAF.pdb")));
        EXPECT_GT(md_system_atom_count(&sys),   0);
        ASSERT_EQ(md_system_entity_count(&sys), 3);
        EXPECT_EQ(md_system_entity_flags(&sys, 0), MD_FLAG_POLYMER | MD_FLAG_POLYPEPTIDE);
        EXPECT_EQ(md_system_entity_flags(&sys, 1), MD_FLAG_HETERO);
        EXPECT_EQ(md_system_entity_flags(&sys, 2), MD_FLAG_HETERO | MD_FLAG_WATER);
        
        ASSERT_EQ(md_system_instance_count(&sys), 3);
        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 0), STR_LIT("A")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 0), STR_LIT("E")));

        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 1), STR_LIT("B")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 1), STR_LIT("E")));

        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 2), STR_LIT("C")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 2), STR_LIT("E")));
    }

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb")));
        EXPECT_GT(md_system_atom_count(&sys),   0);

        ASSERT_EQ(md_system_entity_count(&sys), 8);
        EXPECT_EQ(md_system_entity_flags(&sys, 0), MD_FLAG_POLYMER | MD_FLAG_POLYPEPTIDE);
        EXPECT_EQ(md_system_entity_flags(&sys, 1), MD_FLAG_POLYMER | MD_FLAG_POLYPEPTIDE);
        EXPECT_EQ(md_system_entity_flags(&sys, 2), MD_FLAG_HETERO);
        EXPECT_EQ(md_system_entity_flags(&sys, 3), MD_FLAG_HETERO | MD_FLAG_ION);
        EXPECT_EQ(md_system_entity_flags(&sys, 4), MD_FLAG_HETERO);
        EXPECT_EQ(md_system_entity_flags(&sys, 5), MD_FLAG_HETERO | MD_FLAG_WATER);
        EXPECT_EQ(md_system_entity_flags(&sys, 6), MD_FLAG_HETERO);
        EXPECT_EQ(md_system_entity_flags(&sys, 7), MD_FLAG_HETERO);

        ASSERT_EQ(md_system_instance_count(&sys), 11);
        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 0),      STR_LIT("A")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 0), STR_LIT("A")));
        EXPECT_EQ(md_system_instance_entity_idx(&sys, 0),       0);

        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 1),      STR_LIT("B")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 1), STR_LIT("B")));
        EXPECT_EQ(md_system_instance_entity_idx(&sys, 1),       1);

        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 2),      STR_LIT("C")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 2), STR_LIT("A")));
        EXPECT_EQ(md_system_instance_entity_idx(&sys, 2),       2);

        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 3),      STR_LIT("D")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 3), STR_LIT("A")));
        EXPECT_EQ(md_system_instance_entity_idx(&sys, 3),       3);

        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 4),      STR_LIT("E")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 4), STR_LIT("A")));
        EXPECT_EQ(md_system_instance_entity_idx(&sys, 4),       4);

        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 5),      STR_LIT("F")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 5), STR_LIT("A")));
        EXPECT_EQ(md_system_instance_entity_idx(&sys, 5),       4);

        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 6),      STR_LIT("G")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 6), STR_LIT("A")));
        EXPECT_EQ(md_system_instance_entity_idx(&sys, 6),       5);

        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 7),      STR_LIT("H")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 7), STR_LIT("B")));
        EXPECT_EQ(md_system_instance_entity_idx(&sys, 7),       6);

        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 8),      STR_LIT("I")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 8), STR_LIT("B")));
        EXPECT_EQ(md_system_instance_entity_idx(&sys, 8),       4);

        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 9),      STR_LIT("J")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 9), STR_LIT("B")));
        EXPECT_EQ(md_system_instance_entity_idx(&sys, 9),       5);

        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 10),      STR_LIT("K")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 10), STR_LIT("C")));
        EXPECT_EQ(md_system_instance_entity_idx(&sys, 10),       7);
    }

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/dppc64.pdb")));
        EXPECT_GT(md_system_atom_count(&sys),   0);

        ASSERT_EQ(md_system_entity_count(&sys), 2);
        EXPECT_EQ(md_system_entity_flags(&sys, 0), 0);
        EXPECT_EQ(md_system_entity_flags(&sys, 1), MD_FLAG_WATER);

        ASSERT_EQ(md_system_instance_count(&sys), 65);
        for (size_t i = 0; i < 64; ++i) {
            EXPECT_EQ(md_system_instance_entity_idx(&sys, i), 0);
        }

        EXPECT_EQ(md_system_instance_entity_idx(&sys, 64), 1);
    }

    md_temp_end(temp_scope);
}