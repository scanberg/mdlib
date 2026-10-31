#include "utest.h"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <md_pdb.h>
#include <md_gro.h>
#include <md_xyz.h>
#include <md_mmcif.h>
#include <md_tpr.h>
#include <md_lammps.h>
#include <md_system.h>
#include <md_util.h>
#include <md_dssp.h>
#include <md_itp.h>
#include <core/md_writer.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_bitset.h>
#include <core/md_hash.h>

#include "rmsd.h"
#include "system_invariants.h"

// The tests below lay out coordinates per axis, which reads well; a state holds them packed.
static vec3_t* pack_xyz(vec3_t* dst, const float* x, const float* y, const float* z, size_t n) {
    for (size_t i = 0; i < n; ++i) dst[i] = vec3_set(x[i], y[i], z[i]);
    return dst;
}

// md_util_com_compute over planar inputs: packed first. With indices the source is as long as the
// largest index reaches.
static vec3_t com_planar(const float* x, const float* y, const float* z, const float* w, const int32_t* idx, size_t n, const md_unitcell_t* cell) {
    size_t len = n;
    if (idx) {
        len = 0;
        for (size_t i = 0; i < n; ++i) len = MAX(len, (size_t)idx[i] + 1);
    }
    md_temp_scope_t temp = md_temp_begin();
    vec3_t* xyz = md_temp_alloc_array(temp, vec3_t, ALIGN_TO(len, 16));
    pack_xyz(xyz, x, y, z, len);
    const vec3_t com = md_util_com_compute(xyz, w, idx, n, cell);
    md_temp_end(temp);
    return com;
}

static void unpack_xyz(float* x, float* y, float* z, const vec3_t* src, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        x[i] = src[i].x;
        y[i] = src[i].y;
        z[i] = src[i].z;
    }
}

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
    md_pdb_system_publish_run(&utest_fixture->mol_ala, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), STR_LIT("run/ala"), MD_RUN_FLAG_DISABLE_CACHE_WRITE);
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

UTEST_F(util, bonds) {
    EXPECT_EQ(152,      utest_fixture->mol_ala.bond.count);
    EXPECT_EQ(55,       utest_fixture->mol_pftaa.bond.count);
    EXPECT_EQ(40,       utest_fixture->mol_nucleotides.bond.count);
    EXPECT_EQ(163504,   utest_fixture->mol_centered.bond.count);
}

UTEST_F(util, part) {
    EXPECT_EQ(1,        utest_fixture->mol_ala.part.count);
	EXPECT_EQ(1,        utest_fixture->mol_pftaa.part.count);
	EXPECT_EQ(2,        utest_fixture->mol_nucleotides.part.count);
	EXPECT_EQ(253 + 61, utest_fixture->mol_centered.part.count);

    const md_system_t* sys = &utest_fixture->mol_centered;
    ASSERT(sys->part.count > 253);

    const md_urange_t ref = md_system_part_atom_range(sys, 0);
    for (size_t i = 1; i < 253; ++i) {
        const md_urange_t r = md_system_part_atom_range(sys, i);
        EXPECT_EQ(ref.end - ref.beg, r.end - r.beg);
        EXPECT_EQ(md_system_part_type_idx(sys, 0), md_system_part_type_idx(sys, i));
    }
    // The 253 chains are made alike: one molecule type
    EXPECT_EQ(md_system_molecule_type_idx(sys, 0), md_system_molecule_type_idx(sys, 252));
}

UTEST_F(util, backbone) {
	EXPECT_EQ(0,        utest_fixture->mol_pftaa.protein_backbone.range.count);
	EXPECT_EQ(0,        utest_fixture->mol_nucleotides.protein_backbone.range.count);
    EXPECT_EQ(1,        utest_fixture->mol_ala.protein_backbone.range.count);
    EXPECT_EQ(15,       utest_fixture->mol_ala.protein_backbone.segment.count);
	EXPECT_EQ(253,      utest_fixture->mol_centered.protein_backbone.range.count);
    EXPECT_EQ(10626,    utest_fixture->mol_centered.protein_backbone.segment.count); // Should be equal to the total count of residues in chains
}

UTEST_F(util, molecule_count) {
    size_t num_molecules_pftaa = md_system_molecule_count(&utest_fixture->mol_pftaa);
    EXPECT_EQ(1, num_molecules_pftaa);
    size_t num_molecules_nucleotides = md_system_molecule_count(&utest_fixture->mol_nucleotides);
    EXPECT_EQ(2, num_molecules_nucleotides);
    size_t num_molecules_ala = md_system_molecule_count(&utest_fixture->mol_ala);
	EXPECT_EQ(1, num_molecules_ala);
	size_t num_molecules_centered = md_system_molecule_count(&utest_fixture->mol_centered);
	EXPECT_EQ(253+61, num_molecules_centered);
}

UTEST_F(util, rmsd) {
    md_allocator_i* alloc = utest_fixture->alloc;
    md_system_t* mol = &utest_fixture->mol_ala;
    ASSERT_TRUE(mol);

    const size_t N = mol->atom.count;
    const size_t cap = ALIGN_TO(N, 16);
    vec3_t* xyz[2] = {
        md_alloc(alloc, cap * sizeof(vec3_t)),
        md_alloc(alloc, cap * sizeof(vec3_t)),
    };
    double* xyz0 = md_alloc(alloc, N * 3 * sizeof(double));
    double* xyz1 = md_alloc(alloc, N * 3 * sizeof(double));

    float* w = md_alloc(alloc, sizeof(float) * N);
	md_atom_extract_masses(w, 0, mol->atom.count, &mol->atom);

    {
        const str_t paths[] = { STR_INIT("atom/position") };
        md_system_extract_t* ex = md_system_extract_begin(mol, STR_LIT("run/ala"), paths, 1, md_get_heap_allocator());
        ASSERT_TRUE(ex != NULL);
        ASSERT_TRUE(md_system_extract_frame(ex, 0, &(md_system_state_t){ .xyz = xyz[0] }));
        ASSERT_TRUE(md_system_extract_frame(ex, 1, &(md_system_state_t){ .xyz = xyz[1] }));
        md_system_extract_end(ex);
    }

    for (size_t i = 0; i < N; ++i) {
        xyz0[i * 3 + 0] = xyz[0][i].x;
        xyz0[i * 3 + 1] = xyz[0][i].y;
        xyz0[i * 3 + 2] = xyz[0][i].z;

        xyz1[i * 3 + 0] = xyz[1][i].x;
        xyz1[i * 3 + 1] = xyz[1][i].y;
        xyz1[i * 3 + 2] = xyz[1][i].z;
    }

    // Reference
    double ref_rmsd;
    fast_rmsd((double(*)[3])xyz0, (double(*)[3])xyz1, (int)N, &ref_rmsd);

    // Our implementation
    const vec3_t* const cxyz[2] = { xyz[0], xyz[1] };
    const float* const cw[2] = { w, w };

    vec3_t com[2] = {
        md_util_com_compute(xyz[0], w, 0, N, 0),
        md_util_com_compute(xyz[1], w, 0, N, 0),
    };
    double rmsd = md_util_rmsd_compute(cxyz, cw, 0, N, com);

    EXPECT_LE(fabs(ref_rmsd - rmsd), 0.1);

    md_free(alloc, xyz[0], cap * sizeof(vec3_t));
    md_free(alloc, xyz[1], cap * sizeof(vec3_t));
    md_free(alloc, xyz0, N * 3 * sizeof(double));
    md_free(alloc, xyz1, N * 3 * sizeof(double));
    md_free(alloc, w, N * sizeof(float));
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

UTEST_F(util, molecules) {
    size_t num_molecules = 0;

    num_molecules = md_system_molecule_count(&utest_fixture->mol_nucleotides);
    EXPECT_EQ(num_molecules, 2);

    num_molecules = md_system_molecule_count(&utest_fixture->mol_ala);
	EXPECT_EQ(num_molecules, 1);

	num_molecules = md_system_molecule_count(&utest_fixture->mol_pftaa);
	EXPECT_EQ(num_molecules, 1);

	num_molecules = md_system_molecule_count(&utest_fixture->mol_centered);
    EXPECT_EQ(num_molecules, 253 + 61); // Chains + PFTAA
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

// The molecule hierarchy must be rooted at the graph center and every atom must appear after the
// atom it was reached from. md_util_unwrap_molecule depends on
// that ordering, and atom_slot must be a consistent inverse of atom_idx.
UTEST(util, molecule_hierarchy) {
    md_allocator_i* alloc = md_get_heap_allocator();

    {
        // A chain of nine has a single center, the middle atom
        md_system_t sys = {0};
        make_chain_system(&sys, alloc, 1, 9);
        ASSERT_TRUE(md_util_system_infer_molecules(&sys));

        ASSERT_EQ(md_system_molecule_count(&sys), 1);
        EXPECT_EQ(sys.molecule.atom_idx[0], 4);
        EXPECT_EQ(sys.molecule.parent_idx[0], sys.molecule.atom_idx[0]); // the root is its own parent

        for (size_t slot = 0; slot < 9; ++slot) {
            const int32_t atom   = sys.molecule.atom_idx[slot];
            const int32_t parent = sys.molecule.parent_idx[slot];
            EXPECT_EQ(md_system_molecule_atom_slot(&sys, atom), (int32_t)slot);
            EXPECT_EQ(md_system_molecule_atom_parent(&sys, atom), parent);
            if (parent != atom) {
                // slot(parent) < slot(child), the invariant a single forward pass relies on
                EXPECT_LT(md_system_molecule_atom_slot(&sys, parent), (int32_t)slot);
            }
        }

        md_system_free(&sys);
    }

    {
        // Separate residues are rooted independently and stored back to back
        md_system_t sys = {0};
        make_chain_system(&sys, alloc, 2, 5);
        ASSERT_TRUE(md_util_system_infer_molecules(&sys));

        ASSERT_EQ(md_system_molecule_count(&sys), 2);
        EXPECT_EQ(sys.molecule.offset[0], 0u);
        EXPECT_EQ(sys.molecule.offset[1], 5u);
        EXPECT_EQ(sys.molecule.offset[2], 10u);
        EXPECT_EQ(sys.molecule.atom_idx[0], 2);
        EXPECT_EQ(sys.molecule.atom_idx[5], 7);

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
        ASSERT_TRUE(md_util_system_infer_molecules(&sys));

        ASSERT_EQ(md_system_molecule_count(&sys), 2);
        for (size_t slot = 0; slot < 3; ++slot) {
            EXPECT_EQ(md_system_molecule_atom_slot(&sys, sys.molecule.atom_idx[slot]), (int32_t)slot);
        }

        md_system_free(&sys);
    }
}

// Appends a residue of coarse grained beads. The BB bead of an amino acid residue is its backbone.
static void cg_add_residue(md_system_t* sys, const char* res_name, int seq_id, md_residue_kind_t res_kind, const char* atom_names[], size_t num_atoms, md_allocator_i* alloc) {
    if (sys->atom.type.count == 0) {
        // Type 0 is the "not found" sentinel
        md_particle_type_find_or_add(&sys->atom.type, STR_LIT("?"), 0, 0, 0, 0, 0, alloc);
    }
    if (sys->residue.count == 0) {
        md_array_push(sys->residue.atom_offset, 0, alloc);
    }
    for (size_t i = 0; i < num_atoms; ++i) {
        const md_particle_type_flags_t type_flags = md_particle_type_flags_set_particle_kind(MD_PARTICLE_TYPE_FLAG_NONE, MD_PARTICLE_BEAD);
        const md_atom_flags_t flags = (res_kind == MD_RESIDUE_KIND_AMINO_ACID && strcmp(atom_names[i], "BB") == 0) ? MD_ATOM_FLAG_BACKBONE : MD_ATOM_FLAG_NONE;
        md_particle_type_idx_t type = md_particle_type_find_or_add(&sys->atom.type, str_from_cstr(atom_names[i]), 0, 50.0f, 2.35f, 0xFFFFFFFF, type_flags, alloc);
        md_array_push(sys->atom.type_idx, type, alloc);
        md_array_push(sys->atom.flags, flags, alloc);
        sys->atom.count += 1;
    }
    md_array_push(sys->residue.name, make_label(str_from_cstr(res_name)), alloc);
    md_array_push(sys->residue.seq_id, seq_id, alloc);
    md_array_push(sys->residue.flags, md_residue_flags_set_kind(MD_RESIDUE_FLAG_NONE, res_kind), alloc);
    md_array_push(sys->residue.atom_offset, (uint32_t)sys->atom.count, alloc);
    sys->residue.count += 1;
}

static void check_molecule_invariants(int* utest_result, const md_system_t* sys) {
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

// Coarse grained systems frequently carry no bonds. The molecule hierarchy must then come from the
// residue hierarchy: beads hang off their residue's anchor, consecutive polymer residues are
// chained, non polymer residues stay separate, and nothing is written to sys->bond.
UTEST(util, molecule_hierarchy_coarse_grained_without_bonds) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    md_system_t sys = {0};
    sys.alloc = alloc;

    const char* ala[] = {"BB"};
    const char* lys[] = {"BB", "SC1", "SC2"};
    const char* phe[] = {"BB", "SC1", "SC2", "SC3"};
    const char* popc[] = {"NC3", "PO4", "GL1", "GL2", "C1A", "C1B"};
    const char* w[] = {"W"};

    cg_add_residue(&sys, "LYS", 1, MD_RESIDUE_KIND_AMINO_ACID, lys, 3, alloc); // atoms 0-2
    cg_add_residue(&sys, "ALA", 2, MD_RESIDUE_KIND_AMINO_ACID, ala, 1, alloc); // atom      3
    cg_add_residue(&sys, "PHE", 3, MD_RESIDUE_KIND_AMINO_ACID, phe, 4, alloc); // atoms 4-7
    cg_add_residue(&sys, "ALA", 7, MD_RESIDUE_KIND_AMINO_ACID, ala, 1, alloc); // atom      8, sequence gap: new chain
    cg_add_residue(&sys, "ALA", 8, MD_RESIDUE_KIND_AMINO_ACID, ala, 1, alloc); // atom      9
    cg_add_residue(&sys, "POPC", 9, 0, popc, 6, alloc);                  // atoms 10-15
    cg_add_residue(&sys, "POPC", 10, 0, popc, 6, alloc);                 // atoms 16-21, must not join the previous lipid
    cg_add_residue(&sys, "W", 11, MD_RESIDUE_KIND_WATER, w, 1, alloc);             // atom  22

    ASSERT_TRUE(md_util_system_infer_molecules(&sys));
    EXPECT_EQ(sys.bond.count, 0u);

    ASSERT_EQ(md_system_molecule_count(&sys), 5u);
    check_molecule_invariants(utest_result, &sys);

    // The first chain is a path of three anchors, so the middle one is the root
    EXPECT_EQ(sys.molecule.offset[1] - sys.molecule.offset[0], 8u);
    EXPECT_EQ(sys.molecule.atom_idx[0], 3);
    EXPECT_EQ(md_system_molecule_atom_parent(&sys, 0), 3);
    EXPECT_EQ(md_system_molecule_atom_parent(&sys, 1), 0);
    EXPECT_EQ(md_system_molecule_atom_parent(&sys, 2), 0);
    EXPECT_EQ(md_system_molecule_atom_parent(&sys, 4), 3);
    EXPECT_EQ(md_system_molecule_atom_parent(&sys, 7), 4);

    // Second chain, broken off by the sequence gap
    EXPECT_EQ(sys.molecule.offset[2] - sys.molecule.offset[1], 2u);

    // Each lipid on its own, rooted at its first bead, and the water bead alone
    EXPECT_EQ(sys.molecule.offset[3] - sys.molecule.offset[2], 6u);
    EXPECT_EQ(sys.molecule.offset[4] - sys.molecule.offset[3], 6u);
    EXPECT_EQ(sys.molecule.offset[5] - sys.molecule.offset[4], 1u);
    EXPECT_EQ(md_system_molecule_atom_parent(&sys, 15), 10);
    EXPECT_EQ(md_system_molecule_atom_parent(&sys, 21), 16);

    md_vm_arena_destroy(alloc);
}

// Where bonds exist they are used, and links only fill in what they leave disconnected. A system that is
// not coarse grained gets no hierarchy: only what makes a molecule of a residue (synthesize_molecule_links).
UTEST(util, molecule_hierarchy_coarse_grained_partial_bonds) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    {
        md_system_t sys = {0};
        sys.alloc = alloc;
        const char* lys[] = {"BB", "SC1", "SC2"};
        cg_add_residue(&sys, "LYS", 1, MD_RESIDUE_KIND_AMINO_ACID, lys, 3, alloc);
        cg_add_residue(&sys, "LYS", 2, MD_RESIDUE_KIND_AMINO_ACID, lys, 3, alloc);

        // SC1-SC2 bonded within each residue, nothing else
        md_atom_pair_t p0 = {{1, 2}};
        md_atom_pair_t p1 = {{4, 5}};
        md_array_push(sys.bond.pairs, p0, alloc);
        md_array_push(sys.bond.pairs, p1, alloc);
        sys.bond.count = 2;
        md_bond_build_connectivity(&sys.bond, sys.atom.count, alloc);

        ASSERT_TRUE(md_util_system_infer_molecules(&sys));
        EXPECT_EQ(sys.bond.count, 2u);
        ASSERT_EQ(md_system_molecule_count(&sys), 1u);
        check_molecule_invariants(utest_result, &sys);

        // SC2 is reached through its bond to SC1, not through a redundant link to BB
        EXPECT_EQ(md_system_molecule_atom_parent(&sys, 2), 1);
        EXPECT_EQ(md_system_molecule_atom_parent(&sys, 5), 4);
    }

    {
        md_system_t sys = {0};
        sys.alloc = alloc;
        const char* lys[] = {"BB", "SC1", "SC2"};
        cg_add_residue(&sys, "LYS", 1, MD_RESIDUE_KIND_AMINO_ACID, lys, 3, alloc);
        cg_add_residue(&sys, "LYS", 2, MD_RESIDUE_KIND_AMINO_ACID, lys, 3, alloc);
        for (size_t i = 0; i < sys.atom.type.count; ++i) {
            sys.atom.type.flags[i] = md_particle_type_flags_set_particle_kind(sys.atom.type.flags[i], MD_PARTICLE_ATOM);
        }

        // Atoms without bonds: each residue a molecule, its atoms joined to its first. Without instances there are no
        // chains to join the residues.
        ASSERT_TRUE(md_util_system_infer_molecules(&sys));
        EXPECT_EQ(md_system_molecule_count(&sys), 2u);
        check_molecule_invariants(utest_result, &sys);
    }

    md_vm_arena_destroy(alloc);
}

// A metal coordinated by residues is a molecule of its own, and so are the waters of its shell: the coordination bonds
// between residues hold no molecule together. Unless a topology gives them, as bonds within one of its molecules.
static size_t molecule_size_of_atom(const md_system_t* sys, int32_t atom);

UTEST(util, coordinated_metals_are_molecules_of_their_own) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    {
        // The Mg of tubulin, between the phosphates of its GTP and four waters
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb"), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
        ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
        size_t num_mg = 0;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            if (md_atom_atomic_number(&sys.atom, i) != MD_Z_Mg) continue;
            num_mg += 1;
            EXPECT_EQ(1u, molecule_size_of_atom(&sys, (int32_t)i));
        }
        EXPECT_GT(num_mg, 0u);
        for (size_t c = 0; c < sys.residue.count; ++c) {
            if (md_system_residue_kind(&sys, c) != MD_RESIDUE_KIND_WATER) continue;
            const md_urange_t r = md_system_residue_atom_range(&sys, c);
            EXPECT_EQ((size_t)(r.end - r.beg), molecule_size_of_atom(&sys, (int32_t)r.beg));
        }
        check_molecule_invariants(utest_result, &sys);
    }
    {
        // The Mg of 1fez, coordinated by its protein; and then the same bonds as a topology would give them
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif")));
        ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
        int32_t mg = -1;
        for (size_t i = 0; i < sys.atom.count && mg < 0; ++i) if (md_atom_atomic_number(&sys.atom, i) == MD_Z_Mg) mg = (int32_t)i;
        ASSERT_GE(mg, 0);
        EXPECT_EQ(1u, molecule_size_of_atom(&sys, mg));
        size_t num_coord = 0;
        for (size_t b = 0; b < sys.bond.count; ++b) {
            if (md_bond_kind(sys.bond.flags[b]) != MD_BOND_KIND_COORDINATE) continue;
            sys.bond.flags[b] = md_bond_flags_set_origin(sys.bond.flags[b], MD_BOND_ORIGIN_TOPOLOGY);
            num_coord += 1;
        }
        EXPECT_GT(num_coord, 0u);
        ASSERT_TRUE(md_util_system_infer_molecules(&sys));
        EXPECT_LT(1000u, molecule_size_of_atom(&sys, mg));
        check_molecule_invariants(utest_result, &sys);
    }
    md_vm_arena_destroy(alloc);
}

static size_t molecule_size_of_atom(const md_system_t* sys, int32_t atom) {
    return md_system_molecule_atom_count(sys, (size_t)md_system_molecule_find_by_atom_idx(sys, (size_t)atom));
}

// A molecule holds its residues and its polymer chains whole, whatever its bonds say: a chain with residues missing
// from a crystal structure is one molecule, and the waters of a TIP4P model are one each with their virtual sites
UTEST(util, molecules_hold_residues_and_chains) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    // The molecule of an atom, by the slots of the structures
    #define MOLECULE_OF(sys, atom) md_system_molecule_find_by_atom_idx(&(sys), (size_t)(atom))

    {
        // 8g7u: chains with residues missing, each one molecule over its gaps
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif")));
        ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
        size_t breaks = 0;
        for (size_t i = 0; i < sys.residue.count; ++i) breaks += (md_system_residue_flags(&sys, i) & MD_RESIDUE_FLAG_BREAK_AFTER) != 0;
        EXPECT_GT(breaks, 0u);
        size_t polymers = 0;
        for (size_t i = 0; i < md_system_part_count(&sys); ++i) {
            if (!md_part_kind_is_polymer(md_system_part_kind(&sys, i))) continue;
            polymers += 1;
            const md_urange_t r = md_system_part_atom_range(&sys, i);
            const size_t mol = MOLECULE_OF(sys, r.beg);
            for (uint32_t a = r.beg; a < r.end; ++a) EXPECT_EQ(mol, MOLECULE_OF(sys, a));
        }
        EXPECT_GT(polymers, 1u);
        EXPECT_LE(md_system_molecule_count(&sys), md_system_part_count(&sys));
        check_molecule_invariants(utest_result, &sys);
    }

    {
        // Two DNA strands under one chain id, numbered on: the 5' end of the second has no phosphate, so the chain ends
        // there rather than going on over a gap. Two molecules, one chain.
        static const char pdb[] =
            "ATOM      1  O5'  DA A   1      24.130 139.750  27.200  1.00  0.00           O\n"
            "ATOM      2  C5'  DA A   1      23.020 138.940  27.630  1.00  0.00           C\n"
            "ATOM      3  C4'  DA A   1      21.700 139.500  27.160  1.00  0.00           C\n"
            "ATOM      4  O4'  DA A   1      21.830 139.890  25.780  1.00  0.00           O\n"
            "ATOM      5  C1'  DA A   1      21.020 139.070  24.950  1.00  0.00           C\n"
            "ATOM      6  C3'  DA A   1      20.540 138.510  27.210  1.00  0.00           C\n"
            "ATOM      7  C2'  DA A   1      20.540 137.920  25.820  1.00  0.00           C\n"
            "ATOM      8  O3'  DA A   1      19.300 139.190  27.440  1.00  0.00           O\n"
            "ATOM      9  P    DT A   2      17.930 138.350  27.500  1.00  0.00           P\n"
            "ATOM     10  O5'  DT A   2      17.220 138.660  26.110  1.00  0.00           O\n"
            "ATOM     11  C5'  DT A   2      16.740 139.970  25.830  1.00  0.00           C\n"
            "ATOM     12  C4'  DT A   2      15.860 139.970  24.600  1.00  0.00           C\n"
            "ATOM     13  O4'  DT A   2      16.660 139.700  23.420  1.00  0.00           O\n"
            "ATOM     14  C1'  DT A   2      16.150 138.560  22.730  1.00  0.00           C\n"
            "ATOM     15  C3'  DT A   2      14.730 138.940  24.610  1.00  0.00           C\n"
            "ATOM     16  C2'  DT A   2      15.240 137.830  23.710  1.00  0.00           C\n"
            "ATOM     17  O3'  DT A   2      13.540 139.550  24.110  1.00  0.00           O\n"
            "ATOM     18  O5'  DA A   3      76.380 128.690 -24.150  1.00  0.00           O\n"
            "ATOM     19  C5'  DA A   3      76.750 129.990 -24.620  1.00  0.00           C\n"
            "ATOM     20  C4'  DA A   3      78.030 130.440 -23.960  1.00  0.00           C\n"
            "ATOM     21  O4'  DA A   3      77.720 131.090 -22.700  1.00  0.00           O\n"
            "ATOM     22  C1'  DA A   3      78.420 130.480 -21.630  1.00  0.00           C\n"
            "ATOM     23  C3'  DA A   3      78.960 129.270 -23.630  1.00  0.00           C\n"
            "ATOM     24  C2'  DA A   3      78.820 129.100 -22.130  1.00  0.00           C\n"
            "ATOM     25  O3'  DA A   3      80.320 129.520 -23.980  1.00  0.00           O\n"
            "ATOM     26  P    DT A   4      81.440 128.440 -23.560  1.00  0.00           P\n"
            "ATOM     27  O5'  DT A   4      81.940 128.990 -22.150  1.00  0.00           O\n"
            "ATOM     28  C5'  DT A   4      82.550 130.270 -22.050  1.00  0.00           C\n"
            "ATOM     29  C4'  DT A   4      83.510 130.280 -20.890  1.00  0.00           C\n"
            "ATOM     30  O4'  DT A   4      82.760 130.310 -19.650  1.00  0.00           O\n"
            "ATOM     31  C1'  DT A   4      83.000 129.130 -18.890  1.00  0.00           C\n"
            "ATOM     32  C3'  DT A   4      84.360 129.010 -20.840  1.00  0.00           C\n"
            "ATOM     33  C2'  DT A   4      83.640 128.130 -19.840  1.00  0.00           C\n"
            "ATOM     34  O3'  DT A   4      85.690 129.290 -20.420  1.00  0.00           O\n"
            "END\n";
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, str_from_cstr(pdb), MD_PDB_OPTION_NONE));
        ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
        ASSERT_EQ(sys.atom.count, 34u);
        EXPECT_EQ(md_system_instance_count(&sys), 1u);
        EXPECT_EQ(md_system_part_count(&sys), 2u);
        EXPECT_EQ(md_system_chain_count(&sys), 1u);
        EXPECT_EQ(md_system_chain_atom_count(&sys, 0), 34u);
        EXPECT_EQ(md_system_molecule_count(&sys), 2u);
        EXPECT_NE(MOLECULE_OF(sys, 0), MOLECULE_OF(sys, 17));
        EXPECT_EQ(sys.molecule.offset[sys.molecule.count], 34u);
        for (size_t s = 0; s < md_array_size(sys.molecule.atom_idx); ++s) {
            EXPECT_TRUE(0 <= sys.molecule.atom_idx[s] && sys.molecule.atom_idx[s] < 34);
        }
        check_molecule_invariants(utest_result, &sys);
    }

    {
        // Ions and a water given without any bonds: a molecule each, as their residues are
        static const char pdb[] =
            "HETATM    1 NA    NA A   1       0.000   0.000   0.000  1.00  0.00          NA\n"
            "HETATM    2 CL    CL A   2       5.000   0.000   0.000  1.00  0.00          CL\n"
            "HETATM    3  O   HOH A   3      10.000   0.000   0.000  1.00  0.00           O\n"
            "END\n";
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, str_from_cstr(pdb), MD_PDB_OPTION_NONE));
        ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
        EXPECT_EQ(md_system_bond_count(&sys), 0u);
        EXPECT_EQ(md_system_molecule_count(&sys), 3u);
        check_molecule_invariants(utest_result, &sys);
    }

    #undef MOLECULE_OF
    md_vm_arena_destroy(alloc);
}

// A bondless coarse grained chain wrapped into a small box must come out whole. The BB beads sit 3.8A
// apart along x in a 10A box and each carries a side chain bead 3A up in y, across the cell boundary.
UTEST(util, unwrap_molecule_coarse_grained_without_bonds) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    md_system_t sys = {0};
    sys.alloc = alloc;
    const char* res[] = {"BB", "SC1"};
    enum { N = 7 };
    for (int i = 0; i < N; ++i) {
        cg_add_residue(&sys, "LEU", i + 1, MD_RESIDUE_KIND_AMINO_ACID, res, 2, alloc);
    }
    ASSERT_TRUE(md_util_system_infer_molecules(&sys));
    ASSERT_EQ(md_system_molecule_count(&sys), 1u);

    float x[2 * N], y[2 * N], z[2 * N];
    for (int i = 0; i < N; ++i) {
        x[2 * i]     = fmodf(1.0f + 3.8f * i, 10.0f);
        y[2 * i]     = 8.5f;
        z[2 * i]     = 5.0f;
        x[2 * i + 1] = x[2 * i];
        y[2 * i + 1] = fmodf(8.5f + 3.0f, 10.0f); // 3A up, wrapped to the bottom of the cell
        z[2 * i + 1] = 5.0f;
    }
    vec3_t xyz[2 * N];
    md_system_state_t state = { .num_atoms = 2 * N, .xyz = pack_xyz(xyz, x, y, z, 2 * N), .unitcell = md_unitcell_from_extent(10, 10, 10) };

    md_util_unwrap_system(&state, &sys);
    unpack_xyz(x, y, z, xyz, state.num_atoms);

    for (int i = 0; i + 1 < N; ++i) {
        EXPECT_NEAR(x[2 * (i + 1)] - x[2 * i], 3.8f, 1.0e-4f);
    }
    for (int i = 0; i < N; ++i) {
        EXPECT_NEAR(x[2 * i + 1] - x[2 * i], 0.0f, 1.0e-4f);
        EXPECT_NEAR(y[2 * i + 1] - y[2 * i], 3.0f, 1.0e-4f);
    }

    md_vm_arena_destroy(alloc);
}

// The coarse grained branch of covalent bond inference used to skip the first residue, test every bead
// against itself (always under the cutoff, so a self bond) and emit every other pair twice.
UTEST(util, infer_bonds_coarse_grained) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    md_system_t sys = {0};
    sys.alloc = alloc;
    const char* lys[] = {"BB", "SC1", "SC2"};
    cg_add_residue(&sys, "LYS", 1, MD_RESIDUE_KIND_AMINO_ACID, lys, 3, alloc);
    cg_add_residue(&sys, "LYS", 2, MD_RESIDUE_KIND_AMINO_ACID, lys, 3, alloc);

    // BB beads 3.8A apart, side chain beads 3A above their BB and 3A from each other
    float x[] = {0, 0, 0,   3.8f, 3.8f, 3.8f};
    float y[] = {0, 3, 6,   0,    3,    6};
    float z[] = {0, 0, 0,   0,    0,    0};
    vec3_t xyz[6];
    md_system_state_t state = { .num_atoms = 6, .xyz = pack_xyz(xyz, x, y, z, 6) };

    md_bond_data_t bond = {0};
    md_util_infer_covalent_bonds(&bond, &state, &sys, alloc);

    // (0,1) (1,2) (3,4) (4,5) within residues, (0,3) between backbones. 0-2 is 6A, over the cutoff.
    EXPECT_EQ(bond.count, 5u);
    bool has[6][6] = {0};
    for (size_t i = 0; i < bond.count; ++i) {
        const int a = bond.pairs[i].idx[0];
        const int b = bond.pairs[i].idx[1];
        EXPECT_NE(a, b);
        EXPECT_FALSE(has[a][b] || has[b][a]);
        has[a][b] = true;
    }
    EXPECT_TRUE(has[0][1] && has[1][2] && has[3][4] && has[4][5] && has[0][3]);
    for (size_t i = 0; i < bond.count; ++i) {
        EXPECT_EQ((int)bond.flags[i], (int)md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_INFERRED));
    }

    md_vm_arena_destroy(alloc);
}

// martinize2 (vermouth 0.15, -ff martini3001) on residues 384-392 of 1k4r (YVGELSHQW), with the His renamed HSP as
// martinize2 names a charged one. Martini 3 has Tyr SC1-SC4 and Trp SC1-SC5. The file carries no bonds, so whether
// the chain holds together rests on the backbone beads.
static const char martini3_peptide_gro[] =
    "Martini 3 peptide (martinize2, YVGELSHQW, HIS as HSP)\n"
    "   26\n"
    "    1TYR     BB    1   1.878   1.149   1.615\n"
    "    1TYR    SC1    2   1.634   1.202   1.502\n"
    "    1TYR    SC2    3   1.450   1.367   1.498\n"
    "    1TYR    SC3    4   1.552   1.300   1.291\n"
    "    1TYR    SC4    5   1.363   1.462   1.288\n"
    "    2VAL     BB    6   1.793   1.092   1.968\n"
    "    2VAL    SC1    7   2.047   1.095   2.033\n"
    "    3GLY     BB    8   1.530   1.000   1.928\n"
    "    4GLU     BB    9   1.365   1.301   2.000\n"
    "    4GLU    SC1   10   1.000   1.167   2.124\n"
    "    5LEU     BB   11   1.625   1.416   1.973\n"
    "    5LEU    SC1   12   1.822   1.530   2.206\n"
    "    6SER     BB   13   1.746   1.709   1.799\n"
    "    6SER    SC1   14   1.552   1.747   1.679\n"
    "    7HSP     BB   15   1.960   1.818   1.504\n"
    "    7HSP    SC1   16   2.135   1.902   1.685\n"
    "    7HSP    SC2   17   2.220   2.095   1.803\n"
    "    7HSP    SC3   18   2.253   2.095   1.631\n"
    "    8GLN     BB   19   2.101   2.029   1.314\n"
    "    8GLN    SC1   20   1.781   2.141   1.166\n"
    "    9TRP     BB   21   2.304   1.959   1.000\n"
    "    9TRP    SC1   22   2.470   1.827   1.121\n"
    "    9TRP    SC2   23   2.679   1.768   1.005\n"
    "    9TRP    SC3   24   2.708   1.885   1.138\n"
    "    9TRP    SC4   25   2.709   2.018   1.298\n"
    "    9TRP    SC5   26   2.890   1.973   1.194\n"
    "   5.00000   5.00000   5.00000\n";

// The same as martinize2 wrote it, without its CONECT records
static const char martini3_peptide_pdb[] =
    "ATOM      1 BB   TYR A 384     -61.526 -17.356 229.581  1.00  0.00\n"
    "ATOM      2 SC1  TYR A 384     -63.958 -16.821 228.450  1.00  0.00\n"
    "ATOM      3 SC2  TYR A 384     -65.803 -15.175 228.406  1.00  0.00\n"
    "ATOM      4 SC3  TYR A 384     -64.786 -15.840 226.333  1.00  0.00\n"
    "ATOM      5 SC4  TYR A 384     -66.677 -14.229 226.304  1.00  0.00\n"
    "ATOM      6 BB   VAL A 385     -62.376 -17.927 233.108  1.00  0.00\n"
    "ATOM      7 SC1  VAL A 385     -59.828 -17.893 233.758  1.00  0.00\n"
    "ATOM      8 BB   GLY A 386     -65.004 -18.845 232.710  1.00  0.00\n"
    "ATOM      9 BB   GLU A 387     -66.648 -15.831 233.428  1.00  0.00\n"
    "ATOM     10 SC1  GLU A 387     -70.302 -17.177 234.666  1.00  0.00\n"
    "ATOM     11 BB   LEU A 388     -64.049 -14.684 233.156  1.00  0.00\n"
    "ATOM     12 SC1  LEU A 388     -62.078 -13.547 235.491  1.00  0.00\n"
    "ATOM     13 BB   SER A 389     -62.839 -11.758 231.415  1.00  0.00\n"
    "ATOM     14 SC1  SER A 389     -64.783 -11.373 230.216  1.00  0.00\n"
    "ATOM     15 BB   HSP A 390     -60.700 -10.669 228.464  1.00  0.00\n"
    "ATOM     16 SC1  HSP A 390     -58.952  -9.827 230.275  1.00  0.00\n"
    "ATOM     17 SC2  HSP A 390     -58.102  -7.893 231.462  1.00  0.00\n"
    "ATOM     18 SC3  HSP A 390     -57.772  -7.898 229.740  1.00  0.00\n"
    "ATOM     19 BB   GLN A 391     -59.295  -8.559 226.571  1.00  0.00\n"
    "ATOM     20 SC1  GLN A 391     -62.491  -7.433 225.082  1.00  0.00\n"
    "ATOM     21 BB   TRP A 392     -57.267  -9.252 223.427  1.00  0.00\n"
    "ATOM     22 SC1  TRP A 392     -55.603 -10.572 224.636  1.00  0.00\n"
    "ATOM     23 SC2  TRP A 392     -53.517 -11.161 223.475  1.00  0.00\n"
    "ATOM     24 SC3  TRP A 392     -53.223  -9.997 224.803  1.00  0.00\n"
    "ATOM     25 SC4  TRP A 392     -53.209  -8.667 226.404  1.00  0.00\n"
    "ATOM     26 SC5  TRP A 392     -51.401  -9.117 225.366  1.00  0.00\n"
    "TER      27      TRP A 392\n"
    "END\n";

// Every particle a bead without an element, every residue an amino acid whose BB bead is its one backbone bead and
// whose other beads are side chain, and the whole a single peptide
static void check_martini_peptide(int* utest_result, md_system_t* sys, md_system_state_t* state, size_t num_residues) {
    for (size_t t = 1; t < sys->atom.type.count; ++t) {
        EXPECT_EQ(0, sys->atom.type.z[t]);
        EXPECT_EQ(MD_PARTICLE_BEAD, md_particle_type_particle_kind(&sys->atom.type, t));
    }
    ASSERT_EQ(num_residues, md_system_residue_count(sys));
    for (size_t c = 0; c < num_residues; ++c) {
        EXPECT_EQ(MD_RESIDUE_KIND_AMINO_ACID, md_system_residue_kind(sys, c));
        const md_urange_t range = md_system_residue_atom_range(sys, c);
        size_t num_backbone = 0;
        for (uint32_t i = range.beg; i < range.end; ++i) {
            const md_atom_flags_t flags = md_atom_flags(&sys->atom, i);
            const bool bb = str_eq(md_atom_name(&sys->atom, i), STR_LIT("BB"));
            EXPECT_EQ((int)bb,  (int)((flags & MD_ATOM_FLAG_BACKBONE)   != 0));
            EXPECT_EQ((int)!bb, (int)((flags & MD_ATOM_FLAG_SIDE_CHAIN) != 0));
            num_backbone += (flags & MD_ATOM_FLAG_BACKBONE) != 0;
        }
        EXPECT_EQ((size_t)1, num_backbone);
    }

    ASSERT_TRUE(md_util_system_infer(sys, state, MD_UTIL_INFER_ALL));
    ASSERT_EQ((size_t)1, md_system_part_count(sys));
    EXPECT_EQ(num_residues, md_system_part_residue_range(sys, 0).end - md_system_part_residue_range(sys, 0).beg);
    EXPECT_EQ(MD_PART_KIND_PEPTIDE, md_system_part_kind(sys, 0));
}

UTEST(util, martini3_peptide_gro) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, str_from_cstr(martini3_peptide_gro)));
    check_martini_peptide(utest_result, &sys, &state, 9);
    md_vm_arena_destroy(alloc);
}

UTEST(util, martini3_peptide_pdb) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, str_from_cstr(martini3_peptide_pdb), MD_PDB_OPTION_NONE));
    check_martini_peptide(utest_result, &sys, &state, 9);
    md_vm_arena_destroy(alloc);
}

// The residue names beyond the standard ones that martinize2 writes, in a chain of one each with the beads it gives
// them: the protonation states (CHARMM, Amber, GROMOS and Martini names), the lipidated residues of Martini 3 with their
// tails, Ala, which has a side chain bead in Martini 3 and none in Martini 2, and the charge dummies (SCP, SCN) of
// polarizable Martini 2. BB beads 3.5A apart, and the other beads of a residue 3A apart off its BB.
UTEST(util, martini_residue_names) {
    static const struct { const char* name; const char* beads[11]; } res[] = {
        {"HIS", {"BB", "SC1", "SC2", "SC3"}},
        {"HSD", {"BB", "SC1", "SC2", "SC3"}},
        {"HSE", {"BB", "SC1", "SC2", "SC3"}},
        {"HSP", {"BB", "SC1", "SC2", "SC3"}},
        {"HID", {"BB", "SC1", "SC2", "SC3"}},
        {"HIE", {"BB", "SC1", "SC2", "SC3"}},
        {"HIP", {"BB", "SC1", "SC2", "SC3"}},
        {"HISE", {"BB", "SC1", "SC2", "SC3"}},
        {"HISH", {"BB", "SC1", "SC2", "SC3"}},
        {"HIH", {"BB", "SC1", "SC2", "SC3"}},
        {"LYS", {"BB", "SC1", "SC2"}},
        {"LSN", {"BB", "SC1", "SC2"}},
        {"LYN", {"BB", "SC1", "SC2"}},
        {"LYS0", {"BB", "SC1", "SC2"}},
        {"ASP", {"BB", "SC1"}},
        {"ASH", {"BB", "SC1"}},
        {"ASPP", {"BB", "SC1"}},
        {"ASP0", {"BB", "SC1"}},
        {"GLU", {"BB", "SC1"}},
        {"GLH", {"BB", "SC1"}},
        {"GLUP", {"BB", "SC1"}},
        {"GLU0", {"BB", "SC1"}},
        {"ARG", {"BB", "SC1", "SC2"}},
        {"ARG0", {"BB", "SC1", "SC2"}},
        {"CYS", {"BB", "SC1"}},
        {"CYX", {"BB", "SC1"}},
        {"CYSF", {"BB", "SC1", "C1", "C2", "C3", "C4", "C5", "C6"}},
        {"CYSG", {"BB", "SC1", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"}},
        {"CYSP", {"BB", "SC1", "C1", "C2", "C3", "C4"}},
        {"GLY", {"BB"}},
        {"GLYM", {"BB", "C1", "C2", "C3", "C4"}},
        {"ALA", {"BB", "SC1"}},
        {"PRO", {"BB", "SC1"}},
        {"ASN", {"BB", "SC1", "SCP", "SCN"}},
        {"LYS", {"BB", "SC1", "SC2", "SCP"}},
    };
    const int n = (int)ARRAY_SIZE(res);

    int num_beads = 0;
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < (int)ARRAY_SIZE(res[i].beads) && res[i].beads[k]; ++k) num_beads += 1;
    }
    char text[16384];
    int len = snprintf(text, sizeof(text), "Martini residue names\n%5d\n", num_beads);
    int atom_nr = 0;
    for (int i = 0; i < n; ++i) {
        const float x = 1.0f + 0.35f * (float)i;
        for (int k = 0; k < (int)ARRAY_SIZE(res[i].beads) && res[i].beads[k]; ++k) {
            len += snprintf(text + len, sizeof(text) - len, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", i + 1, res[i].name, res[i].beads[k], ++atom_nr, x, 1.0f + 0.3f * (float)k, 1.0f);
        }
    }
    len += snprintf(text + len, sizeof(text) - len, "  20.00000  20.00000  20.00000\n");
    ASSERT_LT(len, (int)sizeof(text));

    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, (str_t){text, (size_t)len}));
    ASSERT_EQ((size_t)n, md_system_residue_count(&sys));
    check_martini_peptide(utest_result, &sys, &state, (size_t)n);
    md_vm_arena_destroy(alloc);
}

// The bonds of the martinize2 topology of the peptide above (its CONECT records), as atom index pairs
static const int martini3_peptide_topology[][2] = {
    {0, 1}, {0, 5}, {1, 2}, {1, 3}, {2, 3}, {2, 4}, {3, 4}, {5, 6}, {5, 7}, {7, 8}, {8, 9}, {8, 10}, {10, 11}, {10, 12},
    {12, 13}, {12, 14}, {14, 15}, {14, 18}, {15, 16}, {15, 17}, {16, 17}, {18, 19}, {18, 20}, {20, 21}, {21, 22}, {21, 24},
    {22, 24}, {22, 25}, {24, 25},
};

static bool pair_in(const int (*pairs)[2], size_t num_pairs, int a, int b) {
    for (size_t i = 0; i < num_pairs; ++i) {
        if ((pairs[i][0] == a && pairs[i][1] == b) || (pairs[i][0] == b && pairs[i][1] == a)) return true;
    }
    return false;
}

// Every inferred bond is one of the topology's, and the topology's bonds it leaves out are those named
static void check_bonds_against_topology(int* utest_result, const md_system_t* sys, const int (*topology)[2], size_t num_topology, const int (*missing)[2], size_t num_missing) {
    for (size_t i = 0; i < sys->bond.count; ++i) {
        const int a = sys->bond.pairs[i].idx[0];
        const int b = sys->bond.pairs[i].idx[1];
        EXPECT_TRUE_MSG(pair_in(topology, num_topology, a, b), "a bond the topology does not have");
    }
    for (size_t i = 0; i < num_topology; ++i) {
        const bool found = md_bond_find(&sys->bond, topology[i][0], topology[i][1]) != -1;
        EXPECT_EQ((int)!pair_in(missing, num_missing, topology[i][0], topology[i][1]), (int)found);
    }
}

// Inferred from the beads alone, the bonds of the peptide are the topology's: those of each residue (the 4.7 A SC2-SC4
// of the Trp ring too) and the link from each BB to the next. Trp's SC3, a virtual site, is bonded to SC5, the first
// bead it is constructed from, as md_tpr bonds it (martinize2 writes no CONECT record for it). With the CONECT records,
// the bonds are theirs, and inference adds that one.
UTEST(util, martini3_peptide_bonds) {
    int topology[ARRAY_SIZE(martini3_peptide_topology) + 1][2];
    MEMCPY(topology, martini3_peptide_topology, sizeof(martini3_peptide_topology));
    topology[ARRAY_SIZE(martini3_peptide_topology)][0] = 23;
    topology[ARRAY_SIZE(martini3_peptide_topology)][1] = 25;

    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    {
        md_system_t sys = { .alloc = alloc };
        md_system_state_t state = { .alloc = alloc };
        ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, str_from_cstr(martini3_peptide_pdb), MD_PDB_OPTION_NONE));
        check_bonds_against_topology(utest_result, &sys, (const int (*)[2])topology, ARRAY_SIZE(topology), NULL, 0);
        EXPECT_EQ(1u, md_bond_conn_count(&sys.bond, 23));
    }
    {
        // The same with its CONECT records, as martinize2 wrote it
        char text[8192];
        int len = snprintf(text, sizeof(text), "%s", martini3_peptide_pdb);
        len -= (int)strlen("END\n");
        for (size_t i = 0; i < ARRAY_SIZE(martini3_peptide_topology); ++i) {
            len += snprintf(text + len, sizeof(text) - len, "CONECT%5d%5d\n", martini3_peptide_topology[i][0] + 1, martini3_peptide_topology[i][1] + 1);
        }
        len += snprintf(text + len, sizeof(text) - len, "END\n");
        md_system_t sys = { .alloc = alloc };
        md_system_state_t state = { .alloc = alloc };
        ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, (str_t){text, (size_t)len}, MD_PDB_OPTION_NONE));
        check_bonds_against_topology(utest_result, &sys, (const int (*)[2])topology, ARRAY_SIZE(topology), NULL, 0);
        EXPECT_EQ(ARRAY_SIZE(topology), sys.bond.count);
        for (size_t i = 0; i < sys.bond.count; ++i) {
            const bool vsite = pair_in((const int (*)[2])topology + ARRAY_SIZE(martini3_peptide_topology), 1, sys.bond.pairs[i].idx[0], sys.bond.pairs[i].idx[1]);
            EXPECT_EQ(vsite ? MD_BOND_ORIGIN_INFERRED : MD_BOND_ORIGIN_FILE, md_bond_origin(sys.bond.flags[i]));
        }
    }
    md_vm_arena_destroy(alloc);
}

// Polarizable Martini 2 (martinize2 -ff martini22p on villin, residues 53-57): the polar SC1 of Thr and Ser is a
// massless virtual site between two charge dummies, Arg's SC2 carries one. The bonds are the topology's as GROMACS has
// them, all of them: the virtual site is bonded to SCP, the first bead it is constructed from, and the dummies to each
// other (martinize2's CONECT records bond SCN to SC1 as well).
static const char martini22p_peptide_pdb[] =
    "ATOM      1 BB   MET A  53       9.910  -7.189  67.361  1.00  0.00\n"
    "ATOM      2 SC1  MET A  53       7.345  -7.793  70.918  1.00  0.00\n"
    "ATOM      3 BB   THR A  54       7.578  -5.503  65.552  1.00  0.00\n"
    "ATOM      4 SC1  THR A  54       6.937  -7.801  64.213  1.00  0.00\n"
    "ATOM      5 SCP  THR A  54       6.589  -7.804  63.606  1.00  0.00\n"
    "ATOM      6 SCN  THR A  54       7.285  -7.798  64.820  1.00  0.00\n"
    "ATOM      7 BB   ARG A  55       6.195  -2.964  64.028  1.00  0.00\n"
    "ATOM      8 SC1  ARG A  55       7.904  -0.987  62.207  1.00  0.00\n"
    "ATOM      9 SC2  ARG A  55       8.217  -1.441  59.751  1.00  0.00\n"
    "ATOM     10 SCP  ARG A  55       9.317  -1.441  59.751  1.00  0.00\n"
    "ATOM     11 BB   SER A  56       3.400  -3.710  64.222  1.00  0.00\n"
    "ATOM     12 SC1  SER A  56       2.364  -4.957  62.032  1.00  0.00\n"
    "ATOM     13 SCP  SER A  56       1.924  -5.331  61.636  1.00  0.00\n"
    "ATOM     14 SCN  SER A  56       2.803  -4.583  62.428  1.00  0.00\n"
    "ATOM     15 BB   ALA A  57       3.595  -4.775  67.236  1.00  0.00\n"
    "END\n";

UTEST(util, martini22p_peptide) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, str_from_cstr(martini22p_peptide_pdb), MD_PDB_OPTION_NONE));

    static const int topology[][2] = {
        {0, 1}, {0, 2}, {2, 3}, {2, 6}, {3, 4}, {4, 5}, {6, 7}, {6, 10}, {7, 8}, {8, 9}, {10, 11}, {10, 14}, {11, 12}, {12, 13},
    };
    check_bonds_against_topology(utest_result, &sys, topology, ARRAY_SIZE(topology), NULL, 0);

    // The masses of polarizable Martini 2
    EXPECT_EQ(72.0f, md_atom_mass(&sys.atom, 1));     // Met SC1
    EXPECT_EQ( 0.0f, md_atom_mass(&sys.atom, 3));     // Thr SC1, a virtual site
    EXPECT_EQ(36.0f, md_atom_mass(&sys.atom, 4));     // A charge dummy
    EXPECT_EQ(36.0f, md_atom_mass(&sys.atom, 8));     // Arg SC2, charged
    EXPECT_NEAR(1.0f, md_atom_radius(&sys.atom, 4), 1e-4f);
    for (size_t i = 0; i < sys.atom.count; ++i) {
        EXPECT_EQ(MD_PARTICLE_BEAD, md_atom_particle_kind(&sys.atom, i));
    }
    check_martini_peptide(utest_result, &sys, &state, 5);
    md_vm_arena_destroy(alloc);
}

// The masses and sizes of Martini beads are those of the version of Martini the system is from: the same names stand
// for regular, small and tiny beads in Martini 3, for beads of 72 in Martini 2
UTEST(util, martini_masses_by_version) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    {
        // Martini 3 (Tyr has SC4 and Trp SC5)
        md_system_t sys = { .alloc = alloc };
        md_system_state_t state = { .alloc = alloc };
        ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, str_from_cstr(martini3_peptide_gro)));
        EXPECT_EQ(72.0f, md_atom_mass(&sys.atom, 0));     // Tyr BB, regular
        EXPECT_EQ(36.0f, md_atom_mass(&sys.atom, 4));     // Tyr SC4, tiny
        EXPECT_EQ(54.0f, md_atom_mass(&sys.atom, 5));     // Val BB, small
        EXPECT_EQ(54.0f, md_atom_mass(&sys.atom, 7));     // Gly BB
        EXPECT_EQ(54.0f, md_atom_mass(&sys.atom, 11));    // Leu SC1
        EXPECT_EQ(72.0f, md_atom_mass(&sys.atom, 9));     // Glu SC1
        EXPECT_EQ(36.0f, md_atom_mass(&sys.atom, 16));    // His SC2
        EXPECT_EQ( 0.0f, md_atom_mass(&sys.atom, 23));    // Trp SC3, a virtual site
        EXPECT_EQ(36.0f, md_atom_mass(&sys.atom, 25));    // Trp SC5
        EXPECT_NEAR(2.64f, md_atom_radius(&sys.atom, 0), 0.01f);
        EXPECT_NEAR(2.30f, md_atom_radius(&sys.atom, 11), 0.01f);
        EXPECT_NEAR(1.91f, md_atom_radius(&sys.atom, 23), 0.01f);
        // The SC1 beads of Tyr (tiny) and Leu (small) are different types
        EXPECT_NE(sys.atom.type_idx[1], sys.atom.type_idx[11]);
    }
    {
        // Martini 2 (Trp has SC1-SC4, Ala no side chain bead): 72, but 45 for the small beads of the ring
        const char* text =
            "Martini 2\n"
            "    9\n"
            "    1ALA     BB    1   1.000   1.000   1.000\n"
            "    2TRP     BB    2   1.350   1.000   1.000\n"
            "    2TRP    SC1    3   1.350   1.300   1.000\n"
            "    2TRP    SC2    4   1.350   1.600   1.000\n"
            "    2TRP    SC3    5   1.650   1.600   1.000\n"
            "    2TRP    SC4    6   1.650   1.300   1.000\n"
            "    3GLY     BB    7   1.700   1.000   1.000\n"
            "    4LEU     BB    8   2.050   1.000   1.000\n"
            "    4LEU    SC1    9   2.050   1.300   1.000\n"
            "   5.00000   5.00000   5.00000\n";
        md_system_t sys = { .alloc = alloc };
        md_system_state_t state = { .alloc = alloc };
        ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, str_from_cstr(text)));
        for (size_t i = 0; i < sys.atom.count; ++i) {
            const bool ring = 2 <= i && i <= 5;
            EXPECT_EQ(ring ? 45.0f : 72.0f, md_atom_mass(&sys.atom, i));
        }
        EXPECT_NEAR(2.64f, md_atom_radius(&sys.atom, 6), 0.01f);   // Gly BB
        EXPECT_NEAR(2.41f, md_atom_radius(&sys.atom, 4), 0.01f);   // Trp SC3
        EXPECT_NEAR(2.64f, md_atom_radius(&sys.atom, 8), 0.01f);   // Leu SC1
    }
    md_vm_arena_destroy(alloc);
}

// A branched molecule of no polymer (a glycan) is one instance: its third residue is bonded to the first, not to the
// second before it. A second one, bonded to nothing of the first, is an instance of its own.
UTEST(util, branched_molecule_is_one_instance) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    static const char pdb[] =
        "HETATM    1  C1  SGR A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  O1  SGR A   1       1.400   0.000   0.000  1.00  0.00           O\n"
        "HETATM    3  C1  SGR A   2       2.800   0.000   0.000  1.00  0.00           C\n"
        "HETATM    4  O2  SGR A   2       4.200   0.000   0.000  1.00  0.00           O\n"
        "HETATM    5  C1  SGR A   3       0.700   1.212   0.000  1.00  0.00           C\n"
        "HETATM    6  O3  SGR A   3       0.700   2.612   0.000  1.00  0.00           O\n"
        "HETATM    7  C1  SGR A   4      20.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    8  O1  SGR A   4      21.400   0.000   0.000  1.00  0.00           O\n"
        "HETATM    9  C1  SGR A   5      22.800   0.000   0.000  1.00  0.00           C\n"
        "HETATM   10  O2  SGR A   5      24.200   0.000   0.000  1.00  0.00           O\n"
        "END\n";
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, str_from_cstr(pdb), MD_PDB_OPTION_NONE));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    ASSERT_EQ(5u, md_system_residue_count(&sys));
    // The branch: the third bonded to the first
    EXPECT_TRUE(md_bond_find(&sys.bond, 0, 4) >= 0 || md_bond_find(&sys.bond, 1, 4) >= 0);
    EXPECT_LT(md_bond_find(&sys.bond, 2, 4), 0);
    ASSERT_EQ(2u, md_system_part_count(&sys));
    EXPECT_EQ(3u, md_system_part_residue_range(&sys, 0).end - md_system_part_residue_range(&sys, 0).beg);
    EXPECT_EQ(2u, md_system_part_residue_range(&sys, 1).end - md_system_part_residue_range(&sys, 1).beg);
    EXPECT_EQ(2u, md_system_molecule_count(&sys));
    md_vm_arena_destroy(alloc);
}

// The caps of a peptide (ACE, NME) are of its strand: bonded into it by its backbone. The alanine dipeptide of tleap.
UTEST(util, caps_are_of_their_strand) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    static const char pdb[] =
        "ATOM      1 HH31 ACE     1       2.000   1.000  -0.000  1.00  0.00\n"
        "ATOM      2  CH3 ACE     1       2.000   2.090   0.000  1.00  0.00\n"
        "ATOM      3 HH32 ACE     1       1.486   2.454   0.890  1.00  0.00\n"
        "ATOM      4 HH33 ACE     1       1.486   2.454  -0.890  1.00  0.00\n"
        "ATOM      5  C   ACE     1       3.427   2.641  -0.000  1.00  0.00\n"
        "ATOM      6  O   ACE     1       4.391   1.877  -0.000  1.00  0.00\n"
        "ATOM      7  N   ALA     2       3.555   3.970  -0.000  1.00  0.00\n"
        "ATOM      8  H   ALA     2       2.733   4.556  -0.000  1.00  0.00\n"
        "ATOM      9  CA  ALA     2       4.853   4.614  -0.000  1.00  0.00\n"
        "ATOM     10  HA  ALA     2       5.408   4.316   0.890  1.00  0.00\n"
        "ATOM     11  CB  ALA     2       5.661   4.221  -1.232  1.00  0.00\n"
        "ATOM     12  HB1 ALA     2       5.123   4.521  -2.131  1.00  0.00\n"
        "ATOM     13  HB2 ALA     2       6.630   4.719  -1.206  1.00  0.00\n"
        "ATOM     14  HB3 ALA     2       5.809   3.141  -1.241  1.00  0.00\n"
        "ATOM     15  C   ALA     2       4.713   6.129   0.000  1.00  0.00\n"
        "ATOM     16  O   ALA     2       3.601   6.653   0.000  1.00  0.00\n"
        "ATOM     17  N   NME     3       5.846   6.835   0.000  1.00  0.00\n"
        "ATOM     18  H   NME     3       6.737   6.359  -0.000  1.00  0.00\n"
        "ATOM     19  CH3 NME     3       5.846   8.284   0.000  1.00  0.00\n"
        "ATOM     20 HH31 NME     3       4.819   8.648   0.000  1.00  0.00\n"
        "ATOM     21 HH32 NME     3       6.360   8.648   0.890  1.00  0.00\n"
        "ATOM     22 HH33 NME     3       6.360   8.648  -0.890  1.00  0.00\n"
        "END\n";
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, str_from_cstr(pdb), MD_PDB_OPTION_NONE));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    ASSERT_EQ(3u, md_system_residue_count(&sys));
    EXPECT_EQ((int)MD_RESIDUE_KIND_AMINO_ACID, (int)md_system_residue_kind(&sys, 1));
    EXPECT_EQ(1u, md_system_part_count(&sys));
    EXPECT_EQ(1u, md_system_molecule_count(&sys));
    md_vm_arena_destroy(alloc);
}

// A modified residue larger than any standard one: the palmitoylated cysteine of a CHARMM-GUI spike protein (6VSB),
// named CYS, 58 atoms with its hydrogens. Its backbone is found by its atom names whatever its size, as a nucleotide's is.
static const char palmitoyl_cys_pdb[] =
    "ATOM  19115  N   CYS A1236     120.160 136.980 189.740  0.00  0.00           N\n"
    "ATOM  19116  HN  CYS A1236     120.310 136.990 190.720  0.00  0.00           H\n"
    "ATOM  19117  CA  CYS A1236     121.240 136.510 188.940  0.00  0.00           C\n"
    "ATOM  19118  HA  CYS A1236     121.650 137.310 188.340  0.00  0.00           H\n"
    "ATOM  19119  C   CYS A1236     120.710 135.530 187.950  0.00  0.00           C\n"
    "ATOM  19120  O   CYS A1236     120.880 135.590 186.730  0.00  0.00           O\n"
    "ATOM  19121  CB  CYS A1236     122.260 135.850 189.880  0.00  0.00           C\n"
    "ATOM  19122  HB1 CYS A1236     121.720 135.070 190.470  0.00  0.00           H\n"
    "ATOM  19123  HB2 CYS A1236     123.050 135.500 189.190  0.00  0.00           H\n"
    "ATOM  19124  SG  CYS A1236     123.200 136.950 191.060  0.00  0.00           S\n"
    "ATOM  19125  C1  CYS A1236     124.530 135.950 191.790  0.00  0.00           C\n"
    "ATOM  19126  O1  CYS A1236     124.750 134.930 191.150  0.00  0.00           O\n"
    "ATOM  19127  C2  CYS A1236     125.130 136.410 193.130  0.00  0.00           C\n"
    "ATOM  19128  H2A CYS A1236     124.670 137.340 193.530  0.00  0.00           H\n"
    "ATOM  19129  H2B CYS A1236     126.240 136.490 193.130  0.00  0.00           H\n"
    "ATOM  19130  C3  CYS A1236     124.810 135.380 194.180  0.00  0.00           C\n"
    "ATOM  19131  H3A CYS A1236     125.220 134.400 193.880  0.00  0.00           H\n"
    "ATOM  19132  H3B CYS A1236     123.720 135.180 194.310  0.00  0.00           H\n"
    "ATOM  19133  C4  CYS A1236     125.540 135.780 195.480  0.00  0.00           C\n"
    "ATOM  19134  H4A CYS A1236     125.330 136.760 195.940  0.00  0.00           H\n"
    "ATOM  19135  H4B CYS A1236     126.620 135.670 195.230  0.00  0.00           H\n"
    "ATOM  19136  C5  CYS A1236     125.150 134.810 196.570  0.00  0.00           C\n"
    "ATOM  19137  H5A CYS A1236     125.320 133.760 196.230  0.00  0.00           H\n"
    "ATOM  19138  H5B CYS A1236     124.110 134.990 196.910  0.00  0.00           H\n"
    "ATOM  19139  C6  CYS A1236     126.030 135.050 197.800  0.00  0.00           C\n"
    "ATOM  19140  H6A CYS A1236     125.820 136.010 198.320  0.00  0.00           H\n"
    "ATOM  19141  H6B CYS A1236     127.100 134.980 197.490  0.00  0.00           H\n"
    "ATOM  19142  C7  CYS A1236     125.820 133.990 198.900  0.00  0.00           C\n"
    "ATOM  19143  H7A CYS A1236     126.090 132.980 198.520  0.00  0.00           H\n"
    "ATOM  19144  H7B CYS A1236     124.760 133.900 199.240  0.00  0.00           H\n"
    "ATOM  19145  C8  CYS A1236     126.660 134.290 200.100  0.00  0.00           C\n"
    "ATOM  19146  H8A CYS A1236     126.210 135.060 200.780  0.00  0.00           H\n"
    "ATOM  19147  H8B CYS A1236     127.610 134.760 199.770  0.00  0.00           H\n"
    "ATOM  19148  C9  CYS A1236     127.010 133.060 200.990  0.00  0.00           C\n"
    "ATOM  19149  H9A CYS A1236     127.370 132.230 200.350  0.00  0.00           H\n"
    "ATOM  19150  H9B CYS A1236     125.990 132.790 201.330  0.00  0.00           H\n"
    "ATOM  19151  C10 CYS A1236     128.130 133.330 202.020  0.00  0.00           C\n"
    "ATOM  19152 H10A CYS A1236     128.000 134.290 202.560  0.00  0.00           H\n"
    "ATOM  19153 H10B CYS A1236     129.160 133.250 201.630  0.00  0.00           H\n"
    "ATOM  19154  C11 CYS A1236     127.990 132.280 203.120  0.00  0.00           C\n"
    "ATOM  19155 H11A CYS A1236     128.250 131.250 202.800  0.00  0.00           H\n"
    "ATOM  19156 H11B CYS A1236     126.970 132.430 203.520  0.00  0.00           H\n"
    "ATOM  19157  C12 CYS A1236     129.020 132.650 204.310  0.00  0.00           C\n"
    "ATOM  19158 H12A CYS A1236     128.830 133.620 204.800  0.00  0.00           H\n"
    "ATOM  19159 H12B CYS A1236     129.970 132.610 203.740  0.00  0.00           H\n"
    "ATOM  19160  C13 CYS A1236     128.870 131.530 205.400  0.00  0.00           C\n"
    "ATOM  19161 H13A CYS A1236     128.760 130.550 204.900  0.00  0.00           H\n"
    "ATOM  19162 H13B CYS A1236     127.910 131.720 205.930  0.00  0.00           H\n"
    "ATOM  19163  C14 CYS A1236     130.040 131.660 206.440  0.00  0.00           C\n"
    "ATOM  19164 H14A CYS A1236     131.060 131.450 206.040  0.00  0.00           H\n"
    "ATOM  19165 H14B CYS A1236     129.820 130.980 207.280  0.00  0.00           H\n"
    "ATOM  19166  C15 CYS A1236     130.200 133.070 207.240  0.00  0.00           C\n"
    "ATOM  19167 H15A CYS A1236     129.130 133.390 207.220  0.00  0.00           H\n"
    "ATOM  19168 H15B CYS A1236     130.700 133.870 206.650  0.00  0.00           H\n"
    "ATOM  19169  C16 CYS A1236     130.740 133.110 208.620  0.00  0.00           C\n"
    "ATOM  19170 H16A CYS A1236     130.110 132.440 209.250  0.00  0.00           H\n"
    "ATOM  19171 H16B CYS A1236     131.740 132.620 208.620  0.00  0.00           H\n"
    "ATOM  19172 H16C CYS A1236     130.750 134.090 209.140  0.00  0.00           H\n"
    "END\n";

UTEST(util, amino_acid_of_any_size) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, str_from_cstr(palmitoyl_cys_pdb), MD_PDB_OPTION_NONE));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    ASSERT_EQ(1u, md_system_residue_count(&sys));
    ASSERT_EQ(58u, (unsigned)sys.atom.count);
    EXPECT_EQ((int)MD_RESIDUE_KIND_AMINO_ACID, (int)md_system_residue_kind(&sys, 0));
    size_t backbone = 0, side = 0;
    for (size_t i = 0; i < sys.atom.count; ++i) {
        const md_atom_flags_t f = md_atom_flags(&sys.atom, i);
        backbone += (f & MD_ATOM_FLAG_BACKBONE) != 0;
        side += (f & MD_ATOM_FLAG_SIDE_CHAIN) != 0;
        if (str_eq_cstr(md_atom_name(&sys.atom, i), "C16")) EXPECT_TRUE((f & MD_ATOM_FLAG_SIDE_CHAIN) != 0);
    }
    EXPECT_EQ(3u, backbone);        // N, CA, C
    EXPECT_GT(side, 40u);           // CB, SG and the palmitoyl, with their hydrogens
    md_vm_arena_destroy(alloc);
}

// The roles of the atoms of a residue stay within it: a side chain does not run on through the N of a proline into the
// backbone before it, nor through a coordination bond into a metal, and a base is found from C1' whatever atom is
// listed first among its bonds.
UTEST(util, atom_roles_stay_in_their_residue) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/1a64.pdb"), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
        ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
        size_t num_pro = 0, num_with_side_chain = 0, num_aa = 0;
        for (size_t c = 0; c < sys.residue.count; ++c) {
            if (md_system_residue_kind(&sys, c) != MD_RESIDUE_KIND_AMINO_ACID) continue;
            num_aa += 1;
            const str_t resname = md_system_residue_name(&sys, c);
            const bool pro = str_eq_cstr(resname, "PRO");
            num_pro += pro;
            size_t side = 0;
            const md_urange_t r = md_system_residue_atom_range(&sys, c);
            for (uint32_t i = r.beg; i < r.end; ++i) {
                const md_atom_flags_t f = md_atom_flags(&sys.atom, i);
                EXPECT_FALSE((f & MD_ATOM_FLAG_SIDE_CHAIN) && (f & MD_ATOM_FLAG_BACKBONE));
                if (!(f & MD_ATOM_FLAG_SIDE_CHAIN)) continue;
                side += 1;
                const str_t name = md_atom_name(&sys.atom, i);
                if (pro) EXPECT_TRUE(name.len >= 2 && (name.ptr[0] == 'C' || name.ptr[0] == 'H') && strchr("BGD", name.ptr[1]));
            }
            num_with_side_chain += side > 0;
        }
        EXPECT_GT(num_pro, 0u);
        EXPECT_GT(num_with_side_chain, num_aa / 2);
    }
    {
        // The Mg of 1fez is coordinated by an aspartate
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif")));
        ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
        size_t num_mg = 0;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            if (md_atom_atomic_number(&sys.atom, i) != MD_Z_Mg) continue;
            num_mg += 1;
            EXPECT_EQ(0, (int)(md_atom_flags(&sys.atom, i) & (MD_ATOM_FLAG_SIDE_CHAIN | MD_ATOM_FLAG_BACKBONE)));
        }
        EXPECT_EQ(2u, num_mg);
    }
    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        ASSERT_TRUE(md_gro_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro")));
        ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
        size_t num_nuc = 0, num_with_base = 0;
        for (size_t c = 0; c < sys.residue.count; ++c) {
            if (md_system_residue_kind(&sys, c) != MD_RESIDUE_KIND_NUCLEOTIDE) continue;
            num_nuc += 1;
            size_t base = 0;
            const md_urange_t r = md_system_residue_atom_range(&sys, c);
            for (uint32_t i = r.beg; i < r.end; ++i) {
                const md_atom_flags_t f = md_atom_flags(&sys.atom, i);
                EXPECT_FALSE((f & MD_ATOM_FLAG_NUCLEOBASE) && (f & MD_ATOM_FLAG_BACKBONE));
                base += (f & MD_ATOM_FLAG_NUCLEOBASE) != 0;
            }
            num_with_base += base > 0;
        }
        EXPECT_GT(num_nuc, 0u);
        EXPECT_EQ(num_nuc, num_with_base);
        // Nor does a base reach another residue
        for (size_t i = 0; i < sys.atom.count; ++i) {
            const md_residue_idx_t c = md_system_residue_find_by_atom_idx(&sys, i);
            if (md_atom_flags(&sys.atom, i) & MD_ATOM_FLAG_NUCLEOBASE) EXPECT_EQ((int)MD_RESIDUE_KIND_NUCLEOTIDE, (int)md_system_residue_kind(&sys, c));
        }
    }
    md_vm_arena_destroy(alloc);
}

// Martini 2 DNA (Uusitalo et al. 2015), mapped from four nucleotides of a strand of the nucleosome in nucl-dna.gro
// (TGCA, -60 to -57): the phosphate BB1 and sugar BB2, BB3 are the backbone, the rest the base. Each bead is the centre
// of its atoms: BB1 P OP1 OP2 O5', BB2 C5' O4' C4', BB3 C3' C2' C1'; SC1 N9 C4, SC2 C2 N3 (and N2 of G), SC3 C6 N1 and
// N6 (O6 of G), SC4 C8 N7 C5 of a purine; SC1 N1 C6, SC2 N3 C2 O2, SC3 C5 C4 and N4 (O4 C7 of T) of a pyrimidine.
// Within a residue their distances are the bond and constraint lengths of the topology within 0.03 nm. The bonds are
// those of the topology (polyply's DNA_M2), BB3 to the next BB1 included.
static const char martini2_dna_pdb[] =
    "ATOM      1 BB1  DT  A -60      -3.350 101.668  23.203  1.00  0.00\n"
    "ATOM      2 BB2  DT  A -60      -2.020 101.387  20.003  1.00  0.00\n"
    "ATOM      3 BB3  DT  A -60      -1.460  99.483  19.917  1.00  0.00\n"
    "ATOM      4 SC1  DT  A -60       0.850  99.875  21.345  1.00  0.00\n"
    "ATOM      5 SC2  DT  A -60       2.520  99.690  19.910  1.00  0.00\n"
    "ATOM      6 SC3  DT  A -60       2.758  99.835  23.080  1.00  0.00\n"
    "ATOM      7 BB1  DG  A -59      -3.655  97.578  17.750  1.00  0.00\n"
    "ATOM      8 BB2  DG  A -59      -1.057  97.917  15.617  1.00  0.00\n"
    "ATOM      9 BB3  DG  A -59      -0.357  96.060  15.630  1.00  0.00\n"
    "ATOM     10 SC1  DG  A -59       1.545  96.355  18.015  1.00  0.00\n"
    "ATOM     11 SC2  DG  A -59       4.450  96.427  17.973  1.00  0.00\n"
    "ATOM     12 SC3  DG  A -59       3.530  96.343  20.680  1.00  0.00\n"
    "ATOM     13 SC4  DG  A -59       0.887  96.270  19.463  1.00  0.00\n"
    "ATOM     14 BB1  DC  A -58      -0.695  94.367  12.558  1.00  0.00\n"
    "ATOM     15 BB2  DC  A -58       2.600  93.820  13.640  1.00  0.00\n"
    "ATOM     16 BB3  DC  A -58       2.597  91.947  14.360  1.00  0.00\n"
    "ATOM     17 SC1  DC  A -58       2.015  92.455  17.020  1.00  0.00\n"
    "ATOM     18 SC2  DC  A -58       3.943  92.957  18.023  1.00  0.00\n"
    "ATOM     19 SC3  DC  A -58       1.660  92.757  19.563  1.00  0.00\n"
    "ATOM     20 BB1  DA  A -57       3.948  90.553  10.992  1.00  0.00\n"
    "ATOM     21 BB2  DA  A -57       6.963  91.600  12.027  1.00  0.00\n"
    "ATOM     22 BB3  DA  A -57       7.717  89.830  12.417  1.00  0.00\n"
    "ATOM     23 SC1  DA  A -57       6.215  89.725  15.085  1.00  0.00\n"
    "ATOM     24 SC2  DA  A -57       6.710  90.325  17.270  1.00  0.00\n"
    "ATOM     25 SC3  DA  A -57       4.267  89.427  17.783  1.00  0.00\n"
    "ATOM     26 SC4  DA  A -57       4.760  89.077  14.940  1.00  0.00\n"
    "END\n";

UTEST(util, martini2_dna) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, str_from_cstr(martini2_dna_pdb), MD_PDB_OPTION_NONE));
    ASSERT_EQ(4u, md_system_residue_count(&sys));

    for (size_t c = 0; c < 4; ++c) {
        EXPECT_EQ(MD_RESIDUE_KIND_NUCLEOTIDE, md_system_residue_kind(&sys, c));
        const md_urange_t range = md_system_residue_atom_range(&sys, c);
        for (uint32_t i = range.beg; i < range.end; ++i) {
            const str_t name = md_atom_name(&sys.atom, i);
            const md_atom_flags_t flags = md_atom_flags(&sys.atom, i);
            EXPECT_EQ(MD_PARTICLE_BEAD, md_atom_particle_kind(&sys.atom, i));
            EXPECT_EQ(0, md_atom_atomic_number(&sys.atom, i));
            if (str_eq(name, STR_LIT("BB1"))) {
                EXPECT_EQ((int)MD_ATOM_FLAG_BACKBONE, (int)flags);
                EXPECT_EQ(72.0f, md_atom_mass(&sys.atom, i));
            } else if (str_eq(name, STR_LIT("BB2")) || str_eq(name, STR_LIT("BB3"))) {
                EXPECT_EQ((int)(MD_ATOM_FLAG_BACKBONE | MD_ATOM_FLAG_NUCLEOSIDE), (int)flags);
                EXPECT_EQ(45.0f, md_atom_mass(&sys.atom, i));
            } else {
                EXPECT_EQ((int)(MD_ATOM_FLAG_NUCLEOBASE | MD_ATOM_FLAG_NUCLEOSIDE), (int)flags);
                EXPECT_EQ(45.0f, md_atom_mass(&sys.atom, i));
            }
        }
    }

    static const int topology[][2] = {
        {0, 1}, {1, 2}, {2, 3}, {2, 6}, {3, 4}, {3, 5}, {4, 5}, {6, 7}, {7, 8}, {8, 9}, {8, 13}, {9, 10}, {9, 12}, {10, 11},
        {10, 12}, {11, 12}, {13, 14}, {14, 15}, {15, 16}, {15, 19}, {16, 17}, {16, 18}, {17, 18}, {19, 20}, {20, 21},
        {21, 22}, {22, 23}, {22, 25}, {23, 24}, {23, 25}, {24, 25},
    };
    check_bonds_against_topology(utest_result, &sys, topology, ARRAY_SIZE(topology), NULL, 0);

    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    ASSERT_EQ(1u, md_system_part_count(&sys));
    EXPECT_EQ(MD_PART_KIND_DNA, md_system_part_kind(&sys, 0));
    md_vm_arena_destroy(alloc);
}

// Martini 2 DNA likewise: a second strand under the chain id of the first, numbered on, begins without its phosphate
// bead. The chain ends there: two molecules of one chain, each strand whole.
static const char martini2_dna_second_strand_pdb[] =
    "ATOM     27 BB2  DT  A -56      37.980 101.387  20.003  1.00  0.00\n"
    "ATOM     28 BB3  DT  A -56      38.540  99.483  19.917  1.00  0.00\n"
    "ATOM     29 SC1  DT  A -56      40.850  99.875  21.345  1.00  0.00\n"
    "ATOM     30 SC2  DT  A -56      42.520  99.690  19.910  1.00  0.00\n"
    "ATOM     31 SC3  DT  A -56      42.758  99.835  23.080  1.00  0.00\n"
    "ATOM     32 BB1  DG  A -55      36.345  97.578  17.750  1.00  0.00\n"
    "ATOM     33 BB2  DG  A -55      38.943  97.917  15.617  1.00  0.00\n"
    "ATOM     34 BB3  DG  A -55      39.643  96.060  15.630  1.00  0.00\n"
    "ATOM     35 SC1  DG  A -55      41.545  96.355  18.015  1.00  0.00\n"
    "ATOM     36 SC2  DG  A -55      44.450  96.427  17.973  1.00  0.00\n"
    "ATOM     37 SC3  DG  A -55      43.530  96.343  20.680  1.00  0.00\n"
    "ATOM     38 SC4  DG  A -55      40.887  96.270  19.463  1.00  0.00\n"
    "ATOM     39 BB1  DC  A -54      39.305  94.367  12.558  1.00  0.00\n"
    "ATOM     40 BB2  DC  A -54      42.600  93.820  13.640  1.00  0.00\n"
    "ATOM     41 BB3  DC  A -54      42.597  91.947  14.360  1.00  0.00\n"
    "ATOM     42 SC1  DC  A -54      42.015  92.455  17.020  1.00  0.00\n"
    "ATOM     43 SC2  DC  A -54      43.943  92.957  18.023  1.00  0.00\n"
    "ATOM     44 SC3  DC  A -54      41.660  92.757  19.563  1.00  0.00\n"
    "ATOM     45 BB1  DA  A -53      43.948  90.553  10.992  1.00  0.00\n"
    "ATOM     46 BB2  DA  A -53      46.963  91.600  12.027  1.00  0.00\n"
    "ATOM     47 BB3  DA  A -53      47.717  89.830  12.417  1.00  0.00\n"
    "ATOM     48 SC1  DA  A -53      46.215  89.725  15.085  1.00  0.00\n"
    "ATOM     49 SC2  DA  A -53      46.710  90.325  17.270  1.00  0.00\n"
    "ATOM     50 SC3  DA  A -53      44.267  89.427  17.783  1.00  0.00\n"
    "ATOM     51 SC4  DA  A -53      44.760  89.077  14.940  1.00  0.00\n"
    "END\n";

UTEST(util, martini2_dna_strands_under_one_chain_id) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    char text[sizeof(martini2_dna_pdb) + sizeof(martini2_dna_second_strand_pdb)];
    const size_t first_len = sizeof(martini2_dna_pdb) - 1 - 4;     // without its END line
    memcpy(text, martini2_dna_pdb, first_len);
    memcpy(text + first_len, martini2_dna_second_strand_pdb, sizeof(martini2_dna_second_strand_pdb));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, str_from_cstr(text), MD_PDB_OPTION_NONE));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    ASSERT_EQ(8u, md_system_residue_count(&sys));
    ASSERT_EQ(51u, (unsigned)sys.atom.count);
    EXPECT_TRUE(md_system_is_coarse_grained(&sys));
    EXPECT_EQ(2u, md_system_part_count(&sys));
    EXPECT_EQ(1u, md_system_chain_count(&sys));
    EXPECT_EQ(51u, md_system_chain_atom_count(&sys, 0));
    EXPECT_EQ(2u, md_system_molecule_count(&sys));
    for (size_t i = 0; i < md_array_size(sys.molecule.atom_idx); ++i) {
        EXPECT_TRUE(0 <= sys.molecule.atom_idx[i] && sys.molecule.atom_idx[i] < 51);
    }
    check_molecule_invariants(utest_result, &sys);
    md_vm_arena_destroy(alloc);
}

// Between residues: the numbering may jump forward over a bond (an engineered insert) but starting over is a new chain,
// however close its first residue is. Two cysteines of different chains 2.4 A apart are bridged.
UTEST(util, martini_chain_links) {
    const char* text =
        "Chains\n"
        "   10\n"
        "    1ALA     BB    1   1.000   1.000   1.000\n"
        "    1ALA    SC1    2   1.000   0.750   1.000\n"
        "    2CYS     BB    3   1.350   1.000   1.000\n"
        "    2CYS    SC1    4   1.350   1.100   1.200\n"
        "  921ALA     BB    5   1.700   1.000   1.000\n"
        "  921ALA    SC1    6   1.700   0.750   1.000\n"
        "    1ALA     BB    7   1.700   1.350   1.000\n"
        "    1ALA    SC1    8   1.950   1.350   1.000\n"
        "    2CYS     BB    9   1.350   1.350   1.000\n"
        "    2CYS    SC1   10   1.350   1.250   1.380\n"
        "   5.00000   5.00000   5.00000\n";
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, str_from_cstr(text)));
    static const int expected[][2] = {
        {0, 1}, {2, 3}, {4, 5}, {6, 7}, {8, 9},     // BB-SC1
        {0, 2}, {2, 4},                             // 1-2, and 2-921 over the jump
        {6, 8},                                     // the second chain
        {3, 9},                                     // the disulfide
    };
    check_bonds_against_topology(utest_result, &sys, expected, ARRAY_SIZE(expected), NULL, 0);
    EXPECT_EQ(-1, md_bond_find(&sys.bond, 4, 6));  // 3.5 A, but the numbering starts over
    md_vm_arena_destroy(alloc);
}

// A residue whose beads are not those of a topology (here Trp of Martini 3 written without its virtual site SC3) is bonded
// by its geometry: the tree of its beads and the ring bonds besides, but for the 4.7 A SC2-SC4 of the Trp ring. Between
// residues numbered one after the other a link may be stretched to 5.5 A, as in a simulation; over a jump in the
// numbering it is one up to 4.3 A.
UTEST(util, martini_bonds_by_geometry) {
    const char* text =
        "Geometry\n"
        "   11\n"
        "    9TRP     BB    1   2.304   1.959   1.000\n"
        "    9TRP    SC1    2   2.470   1.827   1.121\n"
        "    9TRP    SC2    3   2.679   1.768   1.005\n"
        "    9TRP    SC4    4   2.709   2.018   1.298\n"
        "    9TRP    SC5    5   2.890   1.973   1.194\n"
        "   10ALA     BB    6   2.304   1.959   1.500\n"
        "   10ALA    SC1    7   2.104   1.959   1.500\n"
        "   11ALA     BB    8   2.304   1.959   1.950\n"
        "   11ALA    SC1    9   2.104   1.959   1.950\n"
        "   20ALA     BB   10   2.304   1.959   2.410\n"
        "   20ALA    SC1   11   2.104   1.959   2.410\n"
        "   5.00000   5.00000   5.00000\n";
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, str_from_cstr(text)));
    static const int expected[][2] = {
        {0, 1}, {1, 2}, {2, 4}, {3, 4}, {1, 3},     // Trp, without SC2-SC4 (2-3)
        {5, 6}, {7, 8}, {9, 10},                    // Ala BB-SC1
        {0, 5},                                     // 5.0 A, numbered next
        {5, 7},                                     // 4.5 A, numbered next
    };
    check_bonds_against_topology(utest_result, &sys, expected, ARRAY_SIZE(expected), NULL, 0);
    EXPECT_EQ(-1, md_bond_find(&sys.bond, 7, 9));  // 4.6 A over a jump in the numbering
    md_vm_arena_destroy(alloc);
}

// Martini systems simulated with GROMACS (test_data/martini, see its README): a frame of each as .gro, and its topology
// as .tpr. Read from the frame alone, every bead gets the mass, size, particle kind and role the topology gives it,
// every residue its kind, and the bonds are those of the topology, no more and no fewer. The frames are of a running
// simulation: bonds stretched, molecules broken across the box.
UTEST(util, martini_reference_systems) {
    static const char* names[] = {
        "martini3_membrane",            // POPC, POPE, POPS, CHOL (virtual sites), W, NA, CL
        "martini3_lipids",              // 24 lipids of the Martini 3 lipidome: phospholipids, sphingolipids, glycerides, ...
        "martini3_peptides",            // Charged termini (Pro, Ala), neutral termini (Gly), a disulfide
        "martini3_rna",                 // An RNA duplex (reForge), virtual sites in the bases
        "martini3_small_molecules",     // Small molecules (rings of virtual sites), solvents, sugars, ions, SW, TW, W
        "martini3_polymer",             // PEO with hydroxyl ends and polystyrene (polyply): residues linked by the force field
        "martini2_membrane_protein",    // Ubiquitin on DPPC, DOPC, CHOL; W, WF, NA+, CL-
        "martini2_lipids",              // 24 lipids of the Martini 2 lipidome, surfactants, sterols
        "martini2p_protein",            // Polarizable Martini 2: villin, PW, charge dummies, virtual sites
        "martini2_dna",                 // A DNA duplex (martinize-dna), its elastic network left out
        "martini2_rna",                 // An RNA duplex (martinize-nucleotide)
        "martini2_small_molecules",     // The solvents, sugars and ions of Martini 2
    };
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    for (size_t n = 0; n < ARRAY_SIZE(names); ++n) {
        char tpr_path[256], gro_path[256];
        snprintf(tpr_path, sizeof(tpr_path), MD_UNITTEST_DATA_DIR "/martini/%s.tpr", names[n]);
        snprintf(gro_path, sizeof(gro_path), MD_UNITTEST_DATA_DIR "/martini/%s.gro", names[n]);
        md_system_t top = { .alloc = alloc }, gro = { .alloc = alloc };
        md_system_state_t top_state = { .alloc = alloc }, gro_state = { .alloc = alloc };
        ASSERT_TRUE(md_tpr_system_init_from_file(&top, &top_state, str_from_cstr(tpr_path)));
        ASSERT_TRUE(md_gro_system_init_from_file(&gro, &gro_state, str_from_cstr(gro_path)));
        ASSERT_EQ(top.atom.count, gro.atom.count);
        ASSERT_EQ(top.residue.count, gro.residue.count);

        size_t bad_atoms = 0;
        for (size_t i = 0; i < top.atom.count; ++i) {
            const bool same = fabsf(md_atom_mass(&top.atom, i) - md_atom_mass(&gro.atom, i)) < 1e-3f &&
                fabsf(md_atom_radius(&top.atom, i) - md_atom_radius(&gro.atom, i)) < 0.02f &&
                md_atom_particle_kind(&top.atom, i) == md_atom_particle_kind(&gro.atom, i) &&
                md_atom_flags(&top.atom, i) == md_atom_flags(&gro.atom, i) &&
                md_atom_atomic_number(&gro.atom, i) == 0;
            bad_atoms += !same;
        }
        EXPECT_EQ_MSG((size_t)0, bad_atoms, names[n]);

        size_t bad_residues = 0;
        for (size_t c = 0; c < top.residue.count; ++c) {
            bad_residues += md_system_residue_kind(&top, c) != md_system_residue_kind(&gro, c);
        }
        EXPECT_EQ_MSG((size_t)0, bad_residues, names[n]);

        size_t missing = 0, extra = 0;
        for (size_t b = 0; b < top.bond.count; ++b) {
            missing += md_bond_find(&gro.bond, top.bond.pairs[b].idx[0], top.bond.pairs[b].idx[1]) == -1;
        }
        for (size_t b = 0; b < gro.bond.count; ++b) {
            extra += md_bond_find(&top.bond, gro.bond.pairs[b].idx[0], gro.bond.pairs[b].idx[1]) == -1;
        }
        EXPECT_EQ_MSG((size_t)0, missing, names[n]);
        EXPECT_EQ_MSG((size_t)0, extra, names[n]);
        EXPECT_GT(top.bond.count, (size_t)0);
    }
    md_vm_arena_destroy(alloc);
}

// The nucleic acids of Martini (test_data/martini): each strand a coarse grained nucleic backbone of all its
// nucleotides, read from the frame alone or with its topology. A segment holds the phosphate bead BB1 in p (none at the
// 5' end, which starts a strand) and the sugar beads BB2 and BB3 in c4 and c1, the C4' and C1' they hold; no breaks.
UTEST(util, martini_nucleic_backbones) {
    static const struct {
        const char* name;
        size_t      strand;     // Nucleotides of each of its two strands
        md_part_kind_t kind;
    } systems[] = {
        { "martini2_dna", 12, MD_PART_KIND_DNA },
        { "martini2_rna", 14, MD_PART_KIND_RNA },
        { "martini3_rna", 22, MD_PART_KIND_RNA },
    };
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    for (size_t n = 0; n < ARRAY_SIZE(systems); ++n) {
        for (int with_tpr = 0; with_tpr < 2; ++with_tpr) {
            char path[256];
            snprintf(path, sizeof(path), MD_UNITTEST_DATA_DIR "/martini/%s.%s", systems[n].name, with_tpr ? "tpr" : "gro");
            md_system_t sys = { .alloc = alloc };
            md_system_state_t state = { .alloc = alloc };
            ASSERT_TRUE(with_tpr ? md_tpr_system_init_from_file(&sys, &state, str_from_cstr(path)) : md_gro_system_init_from_file(&sys, &state, str_from_cstr(path)));
            ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));

            const md_nucleic_backbone_data_t* bb = &sys.nucleic_backbone;
            ASSERT_EQ_MSG((size_t)2, bb->range.count, path);
            ASSERT_EQ_MSG(2 * systems[n].strand, bb->segment.count, path);
            for (size_t r = 0; r < bb->range.count; ++r) {
                EXPECT_TRUE(md_nucleic_backbone_range_is_coarse_grained(bb, r));
                EXPECT_EQ(systems[n].strand, (size_t)(bb->range.offset[r + 1] - bb->range.offset[r]));
                EXPECT_TRUE(md_part_kind_is_nucleic_acid(md_system_part_kind(&sys, bb->range.part_idx[r])));
                EXPECT_EQ((int)systems[n].kind, (int)md_system_part_kind(&sys, bb->range.part_idx[r]));
                for (uint32_t s = bb->range.offset[r]; s < bb->range.offset[r + 1]; ++s) {
                    const md_nucleic_acid_atoms_t a = bb->segment.atoms[s];
                    const md_residue_idx_t res = bb->segment.res_idx[s];
                    EXPECT_EQ(-1, a.o5);
                    EXPECT_EQ(-1, a.o3);
                    EXPECT_TRUE(str_eq(md_atom_name(&sys.atom, a.c4), STR_LIT("BB2")));
                    EXPECT_TRUE(str_eq(md_atom_name(&sys.atom, a.c1), STR_LIT("BB3")));
                    EXPECT_EQ(res, (md_residue_idx_t)md_system_residue_find_by_atom_idx(&sys, (size_t)a.c4));
                    if (s == bb->range.offset[r]) {
                        EXPECT_EQ(-1, a.p);
                    } else {
                        ASSERT_GE(a.p, 0);
                        EXPECT_TRUE(str_eq(md_atom_name(&sys.atom, a.p), STR_LIT("BB1")));
                    }
                    EXPECT_EQ(0, (int)(md_system_residue_flags(&sys, (size_t)res) & (MD_RESIDUE_FLAG_BREAK_BEFORE | MD_RESIDUE_FLAG_BREAK_AFTER)));
                }
            }
        }
    }
    md_vm_arena_destroy(alloc);
}

// A strand of Martini 2 DNA (one chain, the first seven nucleotides of martini/martini2_dna but the fourth) is broken
// where the fourth is missing: no bead of the third is bonded to the BB1 of the fifth, 8.8 A away. Two ranges, the
// residues on either side of the break marked as such, and the link over it from the BB3 to the BB1. A nucleotide with
// its BB1 does not start a strand.
static const char martini2_dna_gap_pdb[] =
    "ATOM      1 BB2   DC A   1      33.620  44.890  45.110  1.00  0.00\n"
    "ATOM      2 BB3   DC A   1      33.350  43.570  43.660  1.00  0.00\n"
    "ATOM      3 SC1   DC A   1      31.610  41.910  45.080  1.00  0.00\n"
    "ATOM      4 SC2   DC A   1      32.610  40.180  45.990  1.00  0.00\n"
    "ATOM      5 SC3   DC A   1      29.800  40.440  46.390  1.00  0.00\n"
    "ATOM      6 BB1   DG A   2      36.200  43.640  41.430  1.00  0.00\n"
    "ATOM      7 BB2   DG A   2      36.750  40.590  42.920  1.00  0.00\n"
    "ATOM      8 BB3   DG A   2      36.400  39.240  41.580  1.00  0.00\n"
    "ATOM      9 SC1   DG A   2      33.990  38.330  42.990  1.00  0.00\n"
    "ATOM     10 SC2   DG A   2      33.420  35.860  44.500  1.00  0.00\n"
    "ATOM     11 SC3   DG A   2      30.930  37.410  44.160  1.00  0.00\n"
    "ATOM     12 SC4   DG A   2      32.640  39.140  42.680  1.00  0.00\n"
    "ATOM     13 BB1   DC A   3      39.410  37.770  40.920  1.00  0.00\n"
    "ATOM     14 BB2   DC A   3      37.700  34.940  42.220  1.00  0.00\n"
    "ATOM     15 BB3   DC A   3      36.840  33.840  40.970  1.00  0.00\n"
    "ATOM     16 SC1   DC A   3      34.700  35.320  40.920  1.00  0.00\n"
    "ATOM     17 SC2   DC A   3      33.060  34.200  41.880  1.00  0.00\n"
    "ATOM     18 SC3   DC A   3      32.370  36.510  40.360  1.00  0.00\n"
    "ATOM     19 BB1   DA A   5      35.490  25.450  38.650  1.00  0.00\n"
    "ATOM     20 BB2   DA A   5      32.050  25.760  39.770  1.00  0.00\n"
    "ATOM     21 BB3   DA A   5      30.960  25.540  38.160  1.00  0.00\n"
    "ATOM     22 SC1   DA A   5      30.910  28.390  37.490  1.00  0.00\n"
    "ATOM     23 SC2   DA A   5      29.070  29.740  37.340  1.00  0.00\n"
    "ATOM     24 SC3   DA A   5      30.630  31.410  35.990  1.00  0.00\n"
    "ATOM     25 SC4   DA A   5      32.100  29.000  36.570  1.00  0.00\n"
    "ATOM     26 BB1   DA A   6      29.490  22.800  38.560  1.00  0.00\n"
    "ATOM     27 BB2   DA A   6      26.770  24.910  38.020  1.00  0.00\n"
    "ATOM     28 BB3   DA A   6      26.250  24.770  36.130  1.00  0.00\n"
    "ATOM     29 SC1   DA A   6      27.410  27.350  35.400  1.00  0.00\n"
    "ATOM     30 SC2   DA A   6      26.480  29.330  34.710  1.00  0.00\n"
    "ATOM     31 SC3   DA A   6      28.770  29.740  33.420  1.00  0.00\n"
    "ATOM     32 SC4   DA A   6      28.780  27.100  34.580  1.00  0.00\n"
    "ATOM     33 BB1   DT A   7      23.110  23.480  35.870  1.00  0.00\n"
    "ATOM     34 BB2   DT A   7      22.770  26.750  34.350  1.00  0.00\n"
    "ATOM     35 BB3   DT A   7      22.660  26.780  32.400  1.00  0.00\n"
    "ATOM     36 SC1   DT A   7      25.320  26.770  32.060  1.00  0.00\n"
    "ATOM     37 SC2   DT A   7      25.960  28.690  31.290  1.00  0.00\n"
    "ATOM     38 SC3   DT A   7      27.780  26.040  31.410  1.00  0.00\n"
    "END\n";

UTEST(util, martini_nucleic_break) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, str_from_cstr(martini2_dna_gap_pdb), MD_PDB_OPTION_NONE));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    ASSERT_EQ((size_t)6, md_system_residue_count(&sys));
    EXPECT_EQ(-1, md_bond_find(&sys.bond, 14, 18));    // BB3 of 3 and BB1 of 5
    EXPECT_NE(-1, md_bond_find(&sys.bond, 7, 12));     // BB3 of 2 and BB1 of 3

    const md_nucleic_backbone_data_t* bb = &sys.nucleic_backbone;
    ASSERT_EQ((size_t)2, bb->range.count);
    EXPECT_EQ((size_t)6, bb->segment.count);
    EXPECT_TRUE(md_nucleic_backbone_range_is_coarse_grained(bb, 0));
    EXPECT_TRUE(md_nucleic_backbone_range_is_coarse_grained(bb, 1));
    EXPECT_EQ(3u, bb->range.offset[1]);
    EXPECT_EQ(-1, bb->segment.atoms[0].p);
    EXPECT_EQ(18, bb->segment.atoms[3].p);
    EXPECT_TRUE(md_system_residue_flags(&sys, 2) & MD_RESIDUE_FLAG_BREAK_AFTER);
    EXPECT_TRUE(md_system_residue_flags(&sys, 3) & MD_RESIDUE_FLAG_BREAK_BEFORE);
    md_atom_idx_t link[2];
    ASSERT_TRUE(md_util_chain_break_atoms(link, &sys, 2));
    EXPECT_EQ(14, link[0]);
    EXPECT_EQ(18, link[1]);
    EXPECT_FALSE(md_util_chain_break_atoms(link, &sys, 1));
    md_vm_arena_destroy(alloc);
}

// Bond inference marks what it produces as of origin MD_BOND_ORIGIN_INFERRED, and re-inference replaces only those.
// A bond of another origin (read from a file, user defined, from a topology) is kept wherever it sits in the arrays,
// is never duplicated by an inferred bond naming the same pair, and survives an inference that cannot run.
UTEST(util, infer_bonds_replaces_only_inferred) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    // Four carbons 1.5A apart along x: inference finds 0-1, 1-2, 2-3
    md_system_t sys = { .alloc = alloc };
    md_particle_type_find_or_add(&sys.atom.type, STR_LIT("?"), 0, 0, 0, 0, 0, alloc);
    md_particle_type_idx_t c = md_particle_type_find_or_add(&sys.atom.type, STR_LIT("C"), 6, 12.011f, 0.76f, 0, 0, alloc);
    for (int i = 0; i < 5; ++i) {
        md_array_push(sys.atom.type_idx, c, alloc);
        md_array_push(sys.atom.flags, MD_ATOM_FLAG_NONE, alloc);
        sys.atom.count += 1;
    }
    float x[] = {0.0f, 1.5f, 3.0f, 4.5f, 20.0f};
    float y[] = {0, 0, 0, 0, 0};
    float z[] = {0, 0, 0, 0, 0};
    vec3_t xyz[5];
    md_system_state_t state = { .num_atoms = 5, .xyz = pack_xyz(xyz, x, y, z, 5) };

    md_util_infer_covalent_bonds(&sys.bond, &state, &sys, alloc);
    ASSERT_EQ(sys.bond.count, 3u);
    for (size_t i = 0; i < sys.bond.count; ++i) {
        EXPECT_EQ(MD_BOND_ORIGIN_INFERRED, md_bond_origin(sys.bond.flags[i]));
    }

    // Kept bonds placed FIRST, so nothing may rely on them trailing: a file bond duplicating an inferred pair,
    // a file bond inference would never find, a user bond and a topology bond.
    md_bond_flags_t* flags = 0;
    md_atom_pair_t*  pairs = 0;
    const md_atom_pair_t kept_pairs[] = {{{1, 2}}, {{3, 4}}, {{0, 4}}, {{2, 4}}};
    const md_bond_flags_t kept_flags[] = {MD_BOND_FLAG_NONE, MD_BOND_FLAG_NONE, md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_USER), md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_TOPOLOGY)};
    for (size_t i = 0; i < 4; ++i) {
        md_array_push(pairs, kept_pairs[i], alloc);
        md_array_push(flags, kept_flags[i], alloc);
    }
    for (size_t i = 0; i < sys.bond.count; ++i) {
        md_array_push(pairs, sys.bond.pairs[i], alloc);
        md_array_push(flags, sys.bond.flags[i], alloc);
    }
    md_bond_data_clear(&sys.bond);
    for (size_t i = 0; i < md_array_size(pairs); ++i) {
        md_array_push(sys.bond.pairs, pairs[i], alloc);
        md_array_push(sys.bond.flags, flags[i], alloc);
    }
    sys.bond.count = md_array_size(pairs);

    md_util_infer_covalent_bonds(&sys.bond, &state, &sys, alloc);

    // Inferred 0-1 and 2-3 (1-2 is already kept), then the four kept bonds in their original order
    ASSERT_EQ(sys.bond.count, 6u);
    size_t num_inferred = 0;
    for (size_t i = 0; i < sys.bond.count; ++i) {
        if (md_bond_origin(sys.bond.flags[i]) == MD_BOND_ORIGIN_INFERRED) num_inferred += 1;
    }
    EXPECT_EQ(num_inferred, 2u);
    for (size_t i = 0; i < 4; ++i) {
        EXPECT_EQ(sys.bond.pairs[2 + i].idx[0], kept_pairs[i].idx[0]);
        EXPECT_EQ(sys.bond.pairs[2 + i].idx[1], kept_pairs[i].idx[1]);
        EXPECT_EQ((int)sys.bond.flags[2 + i], (int)kept_flags[i]);
    }

    // No coordinates: nothing is inferred, the kept bonds are still there
    md_system_state_t empty = { .num_atoms = 5 };
    md_util_infer_covalent_bonds(&sys.bond, &empty, &sys, alloc);
    EXPECT_EQ(sys.bond.count, 4u);

    md_vm_arena_destroy(alloc);
}

// A chain of nine atoms spaced 2A along x inside a 10A box, so the wrapped input folds back on
// itself twice. md_util_unwrap_molecule must recover a straight line.
//
// @NOTE: this is the regression test for two bugs that were live together here. min_image_ortho was
// being fed a half extent where it wants a reciprocal one, and the result was written as x + dx
// instead of parent + dx, which is wrong for every non root atom even with no image correction at
// all. Either one alone destroys the spacing checked below.
UTEST(util, unwrap_molecule_ortho) {
    md_allocator_i* alloc = md_get_heap_allocator();

    md_system_t sys = {0};
    make_chain_system(&sys, alloc, 1, 9);
    ASSERT_TRUE(md_util_system_infer_molecules(&sys));

    float x[9], y[9] = {0}, z[9] = {0};
    for (int i = 0; i < 9; ++i) {
        x[i] = fmodf(1.0f + 2.0f * i, 10.0f);
    }
    vec3_t xyz[9];
    md_system_state_t state = { .num_atoms = 9, .xyz = pack_xyz(xyz, x, y, z, 9), .unitcell = md_unitcell_from_extent(10, 10, 10) };

    const md_molecule_t mol = md_system_molecule(&sys, 0);
    ASSERT_LT(0u, (unsigned)mol.count);
    md_util_unwrap_molecule(&state, &mol);
    unpack_xyz(x, y, z, xyz, state.num_atoms);

    for (int i = 0; i < 8; ++i) {
        EXPECT_NEAR(x[i + 1] - x[i], 2.0f, 1.0e-4f);
        EXPECT_NEAR(y[i], 0.0f, 1.0e-4f);
        EXPECT_NEAR(z[i], 0.0f, 1.0e-4f);
    }

    md_system_free(&sys);
}

// The same chain folded back by whole lattice vectors of a sheared cell. Unwrapping must recover the
// original straight line up to a single rigid lattice translation of the whole structure.
UTEST(util, unwrap_molecule_triclinic) {
    md_allocator_i* alloc = md_get_heap_allocator();

    md_system_t sys = {0};
    make_chain_system(&sys, alloc, 1, 9);
    ASSERT_TRUE(md_util_system_infer_molecules(&sys));

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
    vec3_t xyz[9];
    md_system_state_t state = { .num_atoms = 9, .xyz = pack_xyz(xyz, x, y, z, 9), .unitcell = cell };

    const md_molecule_t mol = md_system_molecule(&sys, 0);
    ASSERT_LT(0u, (unsigned)mol.count);
    md_util_unwrap_molecule(&state, &mol);
    unpack_xyz(x, y, z, xyz, state.num_atoms);

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
// residue may exceed half the box and the result must differ from the input by whole box lengths.
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
        const float residual = md_util_optimal_rotation_pbc_vec4_iter(&R, &com, placed, ref, ref_com, trg, PBC_NP, &cell, 8, 1.0e-6f);

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

            const vec3_t plain_u = com_planar(x, y, z, NULL, NULL, n, NULL);
            const vec3_t plain_w = com_planar(x, y, z, w,    NULL, n, NULL);

            // the four internal variants, in order: _com_pbc, _com_pbc_w, _com_pbc_i, _com_pbc_iw
            const vec3_t pbc    = com_planar(x, y, z, NULL, NULL, n, &cells[ci]);
            const vec3_t pbc_w  = com_planar(x, y, z, w,    NULL, n, &cells[ci]);
            const vec3_t pbc_i  = com_planar(x, y, z, NULL, idx,  n, &cells[ci]);
            const vec3_t pbc_iw = com_planar(x, y, z, w,    idx,  n, &cells[ci]);

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

            const vec3_t before    = com_planar(x, y, z, NULL, NULL, n, &cells[ci]);
            const vec3_t before_w  = com_planar(x, y, z, w,    NULL, n, &cells[ci]);
            const vec3_t before_i  = com_planar(x, y, z, NULL, idx,  n, &cells[ci]);
            const vec3_t before_iw = com_planar(x, y, z, w,    idx,  n, &cells[ci]);

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

            const vec3_t after    = com_planar(x, y, z, NULL, NULL, n, &cells[ci]);
            const vec3_t after_w  = com_planar(x, y, z, w,    NULL, n, &cells[ci]);
            const vec3_t after_i  = com_planar(x, y, z, NULL, idx,  n, &cells[ci]);
            const vec3_t after_iw = com_planar(x, y, z, w,    idx,  n, &cells[ci]);

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

        const vec3_t want   = com_planar(xs, ys, zs, NULL, NULL, n, &cells[ci]);
        const vec3_t want_w = com_planar(xs, ys, zs, ws,   NULL, n, &cells[ci]);

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

        const vec3_t got   = com_planar(bx, by, bz, NULL, idx, n, &cells[ci]);
        const vec3_t got_w = com_planar(bx, by, bz, bw,   idx, n, &cells[ci]);

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
    const vec3_t got    = com_planar(x, y, z, w, NULL, 2, &cell);

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
        const float residual = md_util_optimal_rotation_pbc_vec4_iter(&R, &com, placed, ref, ref_com, trg, PBC_NP, &cell, 8, 1.0e-6f);

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
    const float residual = md_util_optimal_rotation_pbc_vec4_iter(&R, &com, NULL, ref, ref_com, trg, PBC_NP, &cell, 8, 1.0e-6f);

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
    EXPECT_NEAR(md_util_optimal_rotation_pbc_vec4_iter(&R, &com, NULL, ref, ref_com, ref, 0, &cell, 8, 1.0e-6f), 0.0f, 1.0e-6f);

    // No cell at all: still a plain Kabsch fit, identity here since target == reference
    const float residual = md_util_optimal_rotation_pbc_vec4_iter(&R, &com, NULL, ref, ref_com, ref, PBC_NP, NULL, 8, 1.0e-6f);
    EXPECT_LT(residual, 1.0e-3f);
    EXPECT_LT(pbc_angle_err(R, mat3_ident()), 0.5f);

    // A set with no coordinates to speak of
    vec3_t c = {0};
    EXPECT_TRUE(md_util_deperiodize_self_vec4(ref, 0, &cell, &c));
}

UTEST_F(util, rings_common) {
    int64_t num_rings = 0;

    num_rings = md_ring_data_count(&utest_fixture->mol_nucleotides.ring);
    EXPECT_EQ(num_rings, 4);
    
    num_rings = md_ring_data_count(&utest_fixture->mol_ala.ring);
    EXPECT_EQ(num_rings, 0);
   
    num_rings = md_ring_data_count(&utest_fixture->mol_pftaa.ring);
    EXPECT_EQ(num_rings, 5);

    num_rings = md_ring_data_count(&utest_fixture->mol_centered.ring);
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

    const size_t num_rings = md_ring_data_count(&sys.ring);
    EXPECT_EQ(num_rings, 32);

    const size_t num_molecules = md_system_molecule_count(&sys);
    EXPECT_EQ(num_molecules, 1);

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

    const size_t num_rings = md_ring_data_count(&sys.ring);
    EXPECT_EQ(num_rings, 362);

    const size_t num_molecules = md_system_molecule_count(&sys);
    EXPECT_EQ(num_molecules, 1);

    md_temp_end(temp_scope);
}

UTEST(util, rings_14kr) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp_scope);
    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    md_pdb_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), MD_PDB_OPTION_NONE);
    md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);

    const size_t num_rings = md_ring_data_count(&sys.ring);
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

    const size_t num_rings = md_ring_data_count(&sys.ring);
    EXPECT_EQ(num_rings, 2);

    const size_t num_molecules = md_system_molecule_count(&sys);
    EXPECT_EQ(num_molecules, 1);

    md_temp_end(temp_scope);
}

UTEST(util, rings_trytophan_xyz) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp_scope);
    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    md_xyz_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan.xyz"), MD_XYZ_OPTION_NONE);
    md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);

    const size_t num_rings = md_ring_data_count(&sys.ring);
    EXPECT_EQ(num_rings, 2);

    const size_t num_molecules = md_system_molecule_count(&sys);
    EXPECT_EQ(num_molecules, 1);

    md_temp_end(temp_scope);
}

UTEST(util, rings_full) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp_scope);
    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    md_xyz_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/full.xyz"), MD_XYZ_OPTION_NONE);
    md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);

    const size_t num_rings = md_ring_data_count(&sys.ring);
    EXPECT_EQ(num_rings, 195);

    const size_t num_molecules = md_system_molecule_count(&sys);
    EXPECT_EQ(num_molecules, 1);

    md_temp_end(temp_scope);
}

UTEST(util, rings_ciprofloxacin) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp_scope);
    md_system_t sys = { .alloc = alloc };
    md_system_state_t sys_state = { .alloc = alloc };
    md_pdb_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/ciprofloxacin.pdb"), MD_PDB_OPTION_NONE);
    md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);

    const int64_t num_rings = md_ring_data_count(&sys.ring);
    ASSERT_EQ(num_rings, 4);
    EXPECT_EQ(md_ring_size(&sys.ring, 0), 3);
    EXPECT_EQ(md_ring_size(&sys.ring, 1), 6);
    EXPECT_EQ(md_ring_size(&sys.ring, 2), 6);
    EXPECT_EQ(md_ring_size(&sys.ring, 3), 6);

    md_temp_end(temp_scope);
}

// Beads make no rings. Their cycles are the rigid triangles Martini closes where the atoms had a ring (the side chains
// of PHE, TYR, TRP, HIS, the small molecules) and, in a protein, those its secondary structure bonds and elastic
// network close: thousands for an elastic network, all of which the ring search and the aromaticity of
// md_chem_perceive went through.
UTEST(util, rings_coarse_grained) {
    static const char* files[] = { "martini3_small_molecules.tpr", "martini2_membrane_protein.tpr", "martini3_peptides.tpr" };
    for (size_t f = 0; f < sizeof(files) / sizeof(files[0]); ++f) {
        md_temp_scope_t temp_scope = md_temp_begin();
        md_allocator_i* alloc = md_temp_allocator(temp_scope);
        md_system_t sys = { .alloc = alloc };
        md_system_state_t sys_state = { .alloc = alloc };
        char path[512];
        snprintf(path, sizeof(path), MD_UNITTEST_DATA_DIR "/martini/%s", files[f]);
        ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &sys_state, str_from_cstr(path)));
        ASSERT_TRUE(md_system_is_coarse_grained(&sys));
        // The bonds of the topology close cycles: more bonds than a forest of the molecules has
        EXPECT_GT(sys.bond.count + sys.instance.count, sys.atom.count);

        md_infer_flags_t flags = MD_UTIL_INFER_ALL & ~MD_UTIL_INFER_BOND_BIT;
        EXPECT_TRUE(md_util_system_infer(&sys, &sys_state, flags));
        EXPECT_EQ(0, (int)md_ring_data_count(&sys.ring));
        md_temp_end(temp_scope);
    }
}

// splitmix64: the same keys on every platform (unlike rand), with all 64 bits random
static uint64_t radix_test_rand(uint64_t* state) {
    uint64_t z = (*state += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

// Sizes around the insertion sort of small arrays (64), the 256 buckets of a pass, and one past the 16 bits of an index
static const size_t radix_test_sizes[] = { 0, 1, 2, 3, 16, 63, 64, 65, 255, 256, 257, 1000, 65537 };
#define RADIX_TEST_MAX_SIZE 65537

// The distributions of keys: the passes the sorts skip (an odd number of passes left leaves the keys in the scratch),
// long runs of equal keys for the stability of the indices, keys in order (left as they are) and in order but the last,
// and the two ends of the range
enum { RADIX_TEST_DIST_COUNT = 11 };
static const char* radix_test_dist_name[RADIX_TEST_DIST_COUNT] = {
    "random", "few distinct", "all equal", "ascending", "descending", "top digit only", "three low digits",
    "packed index pairs", "first and last digit", "both ends of the range", "ascending but the last",
};

static uint64_t radix_test_key64(int dist, size_t i, size_t n, uint64_t* rng) {
    const uint64_t r = radix_test_rand(rng);
    switch (dist) {
    case 0: return r;
    case 1: return r % 7;
    case 2: return 0xDEADBEEFCAFEF00DULL;
    case 3: return (uint64_t)i;
    case 4: return UINT64_MAX - (uint64_t)i;
    case 5: return (r << 56) | 0x0000123456789ABCULL;
    case 6: return (r & 0xFFFFFF) | 0xAB00000000000000ULL;
    case 7: return ((r % 5000) << 32) | ((r >> 32) % 5000);
    case 8: return (r & 0xFF000000000000FFULL) | 0x0011223344556600ULL;
    case 9: return (r & 1) ? UINT64_MAX - (r >> 1) % 1000 : (r >> 1) % 1000;
    case 10: return i + 1 == n ? 0 : (uint64_t)i + 1;
    default: return 0;
    }
}

static uint32_t radix_test_key32(int dist, size_t i, size_t n, uint64_t* rng) {
    const uint64_t r = radix_test_rand(rng);
    switch (dist) {
    case 0: return (uint32_t)r;
    case 1: return (uint32_t)(r % 7);
    case 2: return 0xCAFEF00Du;
    case 3: return (uint32_t)i;
    case 4: return UINT32_MAX - (uint32_t)i;
    case 5: return (uint32_t)(r << 24) | 0x123456u;
    case 6: return (uint32_t)(r & 0xFFFFFF) | 0xAB000000u;
    case 7: return (uint32_t)(((r % 5000) << 16) | ((r >> 32) % 5000));
    case 8: return (uint32_t)(r & 0xFF0000FFu) | 0x00112200u;
    case 9: return (r & 1) ? UINT32_MAX - (uint32_t)((r >> 1) % 1000) : (uint32_t)((r >> 1) % 1000);
    case 10: return i + 1 == n ? 0 : (uint32_t)i + 1;
    default: return 0;
    }
}

typedef struct radix_test_ref_t {
    uint64_t key;
    uint32_t idx;
} radix_test_ref_t;

static int radix_test_ref_cmp(const void* a, const void* b) {
    const radix_test_ref_t* x = (const radix_test_ref_t*)a;
    const radix_test_ref_t* y = (const radix_test_ref_t*)b;
    if (x->key != y->key) return x->key < y->key ? -1 : 1;
    if (x->idx != y->idx) return x->idx < y->idx ? -1 : 1;
    return 0;
}

// The reference, from qsort: the keys in order, and the indices of a stable sort (equal keys by index)
static void radix_test_reference(uint64_t* out_sorted, uint32_t* out_idx, radix_test_ref_t* scratch, const uint64_t* keys, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        scratch[i] = (radix_test_ref_t){ keys[i], (uint32_t)i };
    }
    qsort(scratch, n, sizeof(radix_test_ref_t), radix_test_ref_cmp);
    for (size_t i = 0; i < n; ++i) {
        out_sorted[i] = scratch[i].key;
        out_idx[i] = scratch[i].idx;
    }
}

#define RADIX_TEST_GUARD 0xA5A5A5A5u

// Every size and distribution against the reference: the keys sorted in place, and the indices exactly those of a stable
// sort. A guard past the end of each output catches a write beyond count.
UTEST(util, radix_sort_u32) {
    md_temp_scope_t temp = md_temp_begin();
    const size_t cap = RADIX_TEST_MAX_SIZE + 1;
    uint32_t* keys       = md_temp_alloc_array(temp, uint32_t, cap);
    uint32_t* sorted     = md_temp_alloc_array(temp, uint32_t, cap);
    uint32_t* idx        = md_temp_alloc_array(temp, uint32_t, cap);
    uint64_t* keys_wide  = md_temp_alloc_array(temp, uint64_t, cap);
    uint64_t* ref_sorted = md_temp_alloc_array(temp, uint64_t, cap);
    uint32_t* ref_idx    = md_temp_alloc_array(temp, uint32_t, cap);
    radix_test_ref_t* scratch = md_temp_alloc_array(temp, radix_test_ref_t, cap);

    uint64_t rng = 0x5EED32;
    for (size_t s = 0; s < ARRAY_SIZE(radix_test_sizes); ++s) {
        const size_t n = radix_test_sizes[s];
        for (int d = 0; d < RADIX_TEST_DIST_COUNT; ++d) {
            char msg[128];
            snprintf(msg, sizeof(msg), "%s, n = %zu", radix_test_dist_name[d], n);

            for (size_t i = 0; i < n; ++i) {
                keys[i] = radix_test_key32(d, i, n, &rng);
                keys_wide[i] = keys[i];
            }
            radix_test_reference(ref_sorted, ref_idx, scratch, keys_wide, n);

            memcpy(sorted, keys, n * sizeof(uint32_t));
            sorted[n] = RADIX_TEST_GUARD;
            md_util_sort_radix_u32(sorted, n);
            bool same = true;
            for (size_t i = 0; i < n; ++i) {
                same &= sorted[i] == (uint32_t)ref_sorted[i];
            }
            EXPECT_TRUE_MSG(same, msg);
            EXPECT_EQ(RADIX_TEST_GUARD, sorted[n]);

            idx[n] = RADIX_TEST_GUARD;
            md_util_argsort_radix_u32(idx, keys, n);
            EXPECT_TRUE_MSG(n == 0 || memcmp(idx, ref_idx, n * sizeof(uint32_t)) == 0, msg);
            EXPECT_EQ(RADIX_TEST_GUARD, idx[n]);
        }
    }

    md_temp_end(temp);
}

UTEST(util, radix_sort_u64) {
    md_temp_scope_t temp = md_temp_begin();
    const size_t cap = RADIX_TEST_MAX_SIZE + 1;
    uint64_t* keys       = md_temp_alloc_array(temp, uint64_t, cap);
    uint64_t* sorted     = md_temp_alloc_array(temp, uint64_t, cap);
    uint32_t* idx        = md_temp_alloc_array(temp, uint32_t, cap);
    uint64_t* ref_sorted = md_temp_alloc_array(temp, uint64_t, cap);
    uint32_t* ref_idx    = md_temp_alloc_array(temp, uint32_t, cap);
    radix_test_ref_t* scratch = md_temp_alloc_array(temp, radix_test_ref_t, cap);

    uint64_t rng = 0x5EED64;
    for (size_t s = 0; s < ARRAY_SIZE(radix_test_sizes); ++s) {
        const size_t n = radix_test_sizes[s];
        for (int d = 0; d < RADIX_TEST_DIST_COUNT; ++d) {
            char msg[128];
            snprintf(msg, sizeof(msg), "%s, n = %zu", radix_test_dist_name[d], n);

            for (size_t i = 0; i < n; ++i) {
                keys[i] = radix_test_key64(d, i, n, &rng);
            }
            radix_test_reference(ref_sorted, ref_idx, scratch, keys, n);

            memcpy(sorted, keys, n * sizeof(uint64_t));
            sorted[n] = RADIX_TEST_GUARD;
            md_util_sort_radix_u64(sorted, n);
            EXPECT_TRUE_MSG(n == 0 || memcmp(sorted, ref_sorted, n * sizeof(uint64_t)) == 0, msg);
            EXPECT_EQ((uint64_t)RADIX_TEST_GUARD, sorted[n]);

            idx[n] = RADIX_TEST_GUARD;
            md_util_argsort_radix_u64(idx, keys, n);
            EXPECT_TRUE_MSG(n == 0 || memcmp(idx, ref_idx, n * sizeof(uint32_t)) == 0, msg);
            EXPECT_EQ(RADIX_TEST_GUARD, idx[n]);
        }
    }

    md_temp_end(temp);
}

// Timings on large arrays of random keys, each result checked for order
UTEST(util, radix_sort) {
    uint32_t arr[] = { 1, 278, 128312745, 4, 5, 0, 12382, 26, 12, 14, 7 };
    size_t len = ARRAY_SIZE(arr);

    uint32_t idx[16];

    md_util_argsort_radix_u32(idx, arr, len);
    for (size_t i = 0; i < len - 1; ++i) {
        EXPECT_LE(arr[idx[i]], arr[idx[i+1]]);
    }

    md_util_sort_radix_u32(arr, len);

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

    md_tick_t t0, t1;
    bool ordered;

    t0 = md_tick_now();
    md_util_argsort_radix_u32(indices, values, N);
    t1 = md_tick_now();

    printf("Time for radix argsort (u32): %.4f ms\n", md_tick_to_milliseconds(t1 - t0));
    ordered = true;
    for (size_t i = 1; i < N; ++i) ordered &= values[indices[i - 1]] <= values[indices[i]];
    EXPECT_TRUE(ordered);

    t0 = md_tick_now();
    md_util_sort_radix_u32(values, N);
    t1 = md_tick_now();

    printf("Time for radix sort (u32): %.4f ms\n", md_tick_to_milliseconds(t1 - t0));
    ordered = true;
    for (size_t i = 1; i < N; ++i) ordered &= values[i - 1] <= values[i];
    EXPECT_TRUE(ordered);

    uint64_t* values64 = md_temp_alloc_array(temp, uint64_t, N);
    uint64_t rng = 0x5EED;
    for (size_t i = 0; i < N; ++i) {
        values64[i] = radix_test_rand(&rng);
    }

    t0 = md_tick_now();
    md_util_argsort_radix_u64(indices, values64, N);
    t1 = md_tick_now();

    printf("Time for radix argsort (u64): %.4f ms\n", md_tick_to_milliseconds(t1 - t0));
    ordered = true;
    for (size_t i = 1; i < N; ++i) ordered &= values64[indices[i - 1]] <= values64[indices[i]];
    EXPECT_TRUE(ordered);

    t0 = md_tick_now();
    md_util_sort_radix_u64(values64, N);
    t1 = md_tick_now();

    printf("Time for radix sort (u64): %.4f ms\n", md_tick_to_milliseconds(t1 - t0));
    ordered = true;
    for (size_t i = 1; i < N; ++i) ordered &= values64[i - 1] <= values64[i];
    EXPECT_TRUE(ordered);

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

static md_part_type_flags_t inferred_kind(md_part_kind_t kind) {
    return md_part_type_flags_set_kind(MD_PART_TYPE_FLAG_INFERRED, kind);
}

// The chain of the first residue of a part, "" for one in no chain
static str_t part_chain_id(const md_system_t* sys, size_t part_idx) {
    const md_chain_idx_t c = md_system_chain_find_by_residue_idx(sys, md_system_part_residue_range(sys, part_idx).beg);
    return c >= 0 ? md_system_chain_id(sys, (size_t)c) : STR_LIT("");
}

// Units only where the data gives them (the chains of a PDB file), parts and part types always inferred from them
UTEST(util, part_and_part_type) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp_scope);
    ASSERT(alloc);

    {
        // No chain ids: no units, and the parts are found over the whole system
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb")));
        EXPECT_GT(md_system_atom_count(&sys),   0);
        EXPECT_EQ(md_system_instance_count(&sys), 0);
        EXPECT_EQ(md_system_chain_count(&sys), 0);
        ASSERT_TRUE(md_util_system_infer_parts(&sys));
        ASSERT_EQ(md_system_part_type_count(&sys), 1);
        EXPECT_EQ(md_system_part_type_flags(&sys, 0), inferred_kind(MD_PART_KIND_PEPTIDE));
        EXPECT_TRUE(str_eq(md_system_part_type_description(&sys, 0), STR_LIT("peptide")));
        ASSERT_EQ(md_system_part_count(&sys), 1);
        EXPECT_EQ(md_system_part_type_idx(&sys, 0), 0);
    }

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb")));
        EXPECT_GT(md_system_atom_count(&sys),   0);

        // The units are the chains of the file, named by them and by no made up label
        ASSERT_EQ(md_system_instance_count(&sys), 3);
        static const char* ids[] = {"A", "B", "C"};
        for (size_t i = 0; i < 3; ++i) {
            EXPECT_TRUE(str_empty(md_system_instance_id(&sys, i)));
            EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, i), str_from_cstr(ids[i])));
        }
        ASSERT_TRUE(md_util_system_infer_parts(&sys));
        ASSERT_EQ(md_system_part_type_count(&sys), 1);
        EXPECT_EQ(md_system_part_type_flags(&sys, 0), inferred_kind(MD_PART_KIND_PEPTIDE));
        ASSERT_EQ(md_system_part_count(&sys), 3);
        for (size_t p = 0; p < 3; ++p) {
            EXPECT_EQ(md_system_part_type_idx(&sys, p), 0);
            EXPECT_TRUE(str_eq(part_chain_id(&sys, p), str_from_cstr(ids[p])));
        }
    }

    {
        // One chain, the arginine bound to it and its 106 waters, each water a part of its own. The arginine is numbered
        // on after the chain under its id, but is a free amino acid with a carboxyl group of its own: the chain ends
        // before it (CHAIN BREAKS in md_util.h)
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1LAF.pdb")));
        EXPECT_GT(md_system_atom_count(&sys),   0);
        ASSERT_TRUE(md_util_system_infer_parts(&sys));
        ASSERT_EQ(md_system_part_type_count(&sys), 3);
        EXPECT_EQ(md_system_part_type_flags(&sys, 0), inferred_kind(MD_PART_KIND_PEPTIDE));
        EXPECT_EQ(md_system_part_type_flags(&sys, 2), inferred_kind(MD_PART_KIND_WATER));
        EXPECT_TRUE(str_eq(md_system_part_type_description(&sys, 1), STR_LIT("ARG")));

        ASSERT_EQ(md_system_part_count(&sys), 1 + 1 + 106);
        EXPECT_EQ(md_system_part_residue_range(&sys, 0).end - md_system_part_residue_range(&sys, 0).beg, 238);
        EXPECT_EQ(md_system_part_residue_range(&sys, 1).end - md_system_part_residue_range(&sys, 1).beg, 1);
        for (size_t p = 0; p < md_system_part_count(&sys); ++p) {
            EXPECT_TRUE(str_eq(part_chain_id(&sys, p), STR_LIT("E")));
            if (p < 2) continue;
            EXPECT_EQ(md_system_part_type_idx(&sys, p), 2);
            EXPECT_EQ(md_system_part_residue_range(&sys, p).end - md_system_part_residue_range(&sys, p).beg, 1);
        }
    }

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb")));
        EXPECT_GT(md_system_atom_count(&sys),   0);

        // GTP and GDP are named like nucleotides atom by atom, but are not linked into a chain: ligands
        size_t num_nucleotides = 0;
        for (size_t i = 0; i < md_system_residue_count(&sys); ++i) {
            num_nucleotides += md_system_residue_kind(&sys, i) == MD_RESIDUE_KIND_NUCLEOTIDE;
        }
        EXPECT_EQ(num_nucleotides, 0);

        ASSERT_TRUE(md_util_system_infer_parts(&sys));
        static const struct { md_part_kind_t kind; const char* desc; } types[] = {
            {MD_PART_KIND_PEPTIDE, "peptide"}, {MD_PART_KIND_PEPTIDE, "peptide"}, {MD_PART_KIND_NON_POLYMER, "GTP"}, {MD_PART_KIND_NON_POLYMER, "MG"},
            {MD_PART_KIND_NON_POLYMER, "SO4"}, {MD_PART_KIND_WATER, "water"}, {MD_PART_KIND_NON_POLYMER, "GDP"}, {MD_PART_KIND_NON_POLYMER, "VLB"},
        };
        ASSERT_EQ(md_system_part_type_count(&sys), ARRAY_SIZE(types));
        for (size_t i = 0; i < ARRAY_SIZE(types); ++i) {
            EXPECT_EQ(md_system_part_type_flags(&sys, i), inferred_kind(types[i].kind));
            EXPECT_TRUE(str_eq(md_system_part_type_description(&sys, i), str_from_cstr(types[i].desc)));
        }

        // The strands, then the parts of chain A, of chain B and of chain C: the two SO4 and the five waters of chain A
        // are parts of their own
        static const struct { const char* chain; int type; } parts[] = {
            {"A", 0}, {"B", 1},
            {"A", 2}, {"A", 3}, {"A", 4}, {"A", 4}, {"A", 5}, {"A", 5}, {"A", 5}, {"A", 5}, {"A", 5},
            {"B", 6}, {"B", 4}, {"B", 5}, {"B", 5}, {"B", 5},
            {"C", 7},
        };
        ASSERT_EQ(md_system_part_count(&sys), ARRAY_SIZE(parts));
        for (size_t i = 0; i < ARRAY_SIZE(parts); ++i) {
            EXPECT_TRUE(str_eq(part_chain_id(&sys, i), str_from_cstr(parts[i].chain)));
            EXPECT_EQ(md_system_part_type_idx(&sys, i), parts[i].type);
        }
    }

    {
        // 64 lipids, the first 26 with chain ids of their own, and the waters: a unit for each id, and one for the rest,
        // which the file names nothing
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/dppc64.pdb")));
        EXPECT_GT(md_system_atom_count(&sys),   0);
        EXPECT_EQ(md_system_instance_count(&sys), 27);
        EXPECT_EQ(md_system_chain_count(&sys), 26);
        EXPECT_TRUE(str_empty(md_system_instance_auth_id(&sys, 26)));

        ASSERT_TRUE(md_util_system_infer_parts(&sys));
        ASSERT_EQ(md_system_part_type_count(&sys), 2);
        EXPECT_EQ(md_system_part_type_flags(&sys, 0), inferred_kind(MD_PART_KIND_NON_POLYMER));
        EXPECT_EQ(md_system_part_type_flags(&sys, 1), inferred_kind(MD_PART_KIND_WATER));
        EXPECT_TRUE(str_eq(md_system_part_type_description(&sys, 0), STR_LIT("DPP")));

        ASSERT_EQ(md_system_part_count(&sys), md_system_residue_count(&sys));
        for (size_t i = 0; i < md_system_part_count(&sys); ++i) {
            EXPECT_EQ(md_system_part_type_idx(&sys, i), i < 64 ? 0 : 1);
        }
    }

    md_temp_end(temp_scope);
}

// Part types the file does not define are inferred and flagged as such, file defined ones are not
UTEST(util, part_type_derived) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp_scope);

    {
        // Nucleosome: two DNA strands and eight histones, a .gro carries no units and no types
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_gro_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro")));
        md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);
        EXPECT_EQ(md_system_instance_count(&sys), 0);

        ASSERT_EQ(md_system_part_type_count(&sys), 6);
        // The DNA strands are DNA (by their residue names), which the nucleic backbone extraction relies upon
        EXPECT_EQ(md_system_part_type_flags(&sys, 0), inferred_kind(MD_PART_KIND_DNA));
        EXPECT_EQ(md_system_part_type_flags(&sys, 1), inferred_kind(MD_PART_KIND_DNA));
        EXPECT_TRUE(str_eq(md_system_part_type_description(&sys, 0), STR_LIT("DNA")));
        for (size_t i = 2; i < 6; ++i) {
            EXPECT_EQ(md_system_part_type_flags(&sys, i), inferred_kind(MD_PART_KIND_PEPTIDE));
        }
        EXPECT_EQ(md_system_part_count(&sys), 10);
        EXPECT_EQ(sys.nucleic_backbone.range.count, 2);
        EXPECT_EQ(sys.protein_backbone.range.count, 8);

        // Deriving again starts over rather than appending
        ASSERT_TRUE(md_util_system_infer_parts(&sys));
        EXPECT_EQ(md_system_part_type_count(&sys), 6);
        EXPECT_EQ(md_system_part_count(&sys), 10);
    }

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif")));
        md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);

        // Each asym is one part, of the entity of the file
        ASSERT_GT(md_system_part_type_count(&sys), 0);
        for (size_t i = 0; i < md_system_part_type_count(&sys); ++i) {
            EXPECT_FALSE(md_system_part_type_flags(&sys, i) & MD_PART_TYPE_FLAG_INFERRED);
        }
        ASSERT_EQ(md_system_part_count(&sys), md_system_instance_count(&sys));
        for (size_t p = 0; p < md_system_part_count(&sys); ++p) {
            EXPECT_EQ(md_system_part_type_idx(&sys, p), md_system_instance_part_type_idx(&sys, p));
        }
    }

    md_temp_end(temp_scope);
}

// mmCIF: the part types are the file's entities, and the molecules of a non-polymer asym (the waters of a chain) are
// instances of their own which share the asym's id
UTEST(util, part_type_instance_mmcif) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif")));

        static const struct { md_part_kind_t kind; const char* desc; } types[] = {
            {MD_PART_KIND_PEPTIDE, "PHOSPHONOACETALDEHYDE HYDROLASE"}, {MD_PART_KIND_NON_POLYMER, "MAGNESIUM ION"},
            {MD_PART_KIND_NON_POLYMER, "TUNGSTATE(VI)ION"}, {MD_PART_KIND_WATER, "water"},
        };
        ASSERT_EQ(md_system_part_type_count(&sys), ARRAY_SIZE(types));
        for (size_t i = 0; i < ARRAY_SIZE(types); ++i) {
            EXPECT_EQ(md_system_part_type_kind(&sys, i), types[i].kind);
            EXPECT_TRUE(str_eq(md_system_part_type_description(&sys, i), str_from_cstr(types[i].desc)));
        }

        static const struct { const char* id; const char* auth; int type; } instances[] = {
            {"A", "A", 0}, {"B", "B", 0}, {"E", "A", 1}, {"F", "B", 2}, {"G", "B", 1},
            {"K", "A", 3}, {"K", "A", 3}, {"K", "A", 3}, {"K", "A", 3},
            {"L", "B", 3}, {"L", "B", 3}, {"L", "B", 3}, {"L", "B", 3},
        };
        ASSERT_EQ(md_system_instance_count(&sys), ARRAY_SIZE(instances));
        for (size_t i = 0; i < ARRAY_SIZE(instances); ++i) {
            EXPECT_TRUE(str_eq(md_system_instance_id(&sys, i),      str_from_cstr(instances[i].id)));
            EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, i), str_from_cstr(instances[i].auth)));
            EXPECT_EQ(md_system_instance_part_type_idx(&sys, i),    instances[i].type);
            if (i >= 2) EXPECT_EQ(md_system_instance_residue_count(&sys, i), 1);
        }
        // MG is an ion as a residue, whatever its entity is
        EXPECT_EQ(md_system_residue_kind(&sys, md_system_instance_residue_range(&sys, 2).beg), MD_RESIDUE_KIND_ION);
    }

    {
        // RNA strands which begin with a GTP (the 5' triphosphate): linked into the chain, so a nucleotide
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif")));
        md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);
        ASSERT_GE(md_system_part_type_count(&sys), 5);
        EXPECT_EQ(md_system_part_type_kind(&sys, 2), MD_PART_KIND_RNA);
        EXPECT_EQ(md_system_part_type_kind(&sys, 3), MD_PART_KIND_RNA);
        size_t resolved = 0;
        for (size_t p = 0; p < md_system_part_count(&sys); ++p) {
            if (!md_part_kind_is_nucleic_acid(md_system_part_kind(&sys, p))) continue;
            const md_urange_t range = md_system_part_residue_range(&sys, p);
            EXPECT_TRUE(str_eq(md_system_residue_name(&sys, range.beg), STR_LIT("GTP")) || str_eq(md_system_residue_name(&sys, range.beg), STR_LIT("UTP")));
            for (uint32_t c = range.beg; c < range.end; ++c) {
                EXPECT_EQ(md_system_residue_kind(&sys, c), MD_RESIDUE_KIND_NUCLEOTIDE);
                resolved += (md_system_residue_flags(&sys, c) & MD_RESIDUE_FLAG_RESOLVED) != 0;
            }
        }
        EXPECT_EQ(resolved, 48);
        EXPECT_EQ(sys.nucleic_backbone.range.count, 2);
    }

    md_vm_arena_destroy(alloc);
}

// A hypermodified nucleotide of a tRNA: wybutosine (YYG) at 37 of the anticodon loop, between A36 and A38, its heavy
// atoms only (39 of them, more than any standard nucleotide has with its hydrogens). Its name is no standard one, but it
// is linked into the chain by its phosphodiester bonds: a nucleotide, on the backbone with its neighbours, its base the
// base. From the tRNA of rendering_all_atom.gro, a membrane complex of a deposited structure.
UTEST(util, hypermodified_nucleotide) {
    const char* text =
        "tRNA A36 YYG37 A38\n"
        "   83\n"
        "   36A        P 4250  15.017  13.816   1.786\n"
        "   36A      OP1 4251  14.910  13.913   1.824\n"
        "   36A      OP2 4252  15.149  13.867   1.738\n"
        "   36A      O5' 4253  15.043  13.715   1.907\n"
        "   36A      C5' 4254  14.936  13.633   1.952\n"
        "   36A      C4' 4255  14.986  13.521   2.041\n"
        "   36A      O4' 4256  15.062  13.421   1.969\n"
        "   36A      C3' 4257  15.077  13.551   2.159\n"
        "   36A      O3' 4258  15.009  13.616   2.268\n"
        "   36A      C2' 4259  15.122  13.410   2.199\n"
        "   36A      O2' 4260  15.033  13.344   2.287\n"
        "   36A      C1' 4261  15.120  13.334   2.065\n"
        "   36A       N9 4262  15.254  13.293   2.026\n"
        "   36A       C8 4263  15.349  13.351   1.945\n"
        "   36A       N7 4264  15.464  13.287   1.944\n"
        "   36A       C5 4265  15.442  13.178   2.028\n"
        "   36A       C6 4266  15.525  13.070   2.070\n"
        "   36A       N6 4267  15.652  13.053   2.030\n"
        "   36A       N1 4268  15.471  12.979   2.155\n"
        "   36A       C2 4269  15.344  12.994   2.192\n"
        "   36A       N3 4270  15.257  13.089   2.159\n"
        "   36A       C4 4271  15.313  13.179   2.076\n"
        "   37YYG      P 4272  15.093  13.697   2.379\n"
        "   37YYG    OP1 4273  14.999  13.793   2.445\n"
        "   37YYG    OP2 4274  15.218  13.749   2.316\n"
        "   37YYG    O5' 4275  15.137  13.584   2.483\n"
        "   37YYG    C5' 4276  15.038  13.502   2.550\n"
        "   37YYG    C4' 4277  15.103  13.379   2.608\n"
        "   37YYG    O4' 4278  15.148  13.291   2.502\n"
        "   37YYG    C3' 4279  15.229  13.406   2.688\n"
        "   37YYG    O3' 4280  15.207  13.449   2.822\n"
        "   37YYG    C2' 4281  15.304  13.273   2.681\n"
        "   37YYG    O2' 4282  15.260  13.178   2.775\n"
        "   37YYG    C1' 4283  15.270  13.227   2.539\n"
        "   37YYG     N9 4284  15.375  13.266   2.444\n"
        "   37YYG     C8 4285  15.370  13.366   2.351\n"
        "   37YYG     N7 4286  15.479  13.380   2.282\n"
        "   37YYG     C5 4287  15.563  13.281   2.332\n"
        "   37YYG     C6 4288  15.697  13.247   2.295\n"
        "   37YYG     O6 4289  15.767  13.299   2.207\n"
        "   37YYG     N1 4290  15.746  13.142   2.373\n"
        "   37YYG     C2 4291  15.678  13.080   2.473\n"
        "   37YYG     N2 4292  15.750  12.987   2.538\n"
        "   37YYG     N3 4293  15.553  13.107   2.506\n"
        "   37YYG     C3 4294  15.481  13.032   2.611\n"
        "   37YYG     C4 4295  15.501  13.209   2.432\n"
        "   37YYG    C10 4296  15.995  12.939   2.545\n"
        "   37YYG    C11 4297  15.877  13.004   2.490\n"
        "   37YYG    C12 4298  15.876  13.098   2.377\n"
        "   37YYG    C13 4299  15.995  13.119   2.280\n"
        "   37YYG    C14 4300  15.988  13.001   2.179\n"
        "   37YYG    C15 4301  16.023  13.027   2.031\n"
        "   37YYG    C16 4302  16.005  12.889   1.971\n"
        "   37YYG    O17 4303  15.889  12.843   1.960\n"
        "   37YYG    O18 4304  16.103  12.820   1.932\n"
        "   37YYG    C19 4305  16.182  12.868   1.821\n"
        "   37YYG    N20 4306  15.929  13.130   1.971\n"
        "   37YYG    C21 4307  15.923  13.166   1.837\n"
        "   37YYG    O22 4308  15.817  13.160   1.771\n"
        "   37YYG    O23 4309  16.028  13.206   1.781\n"
        "   37YYG    C24 4310  16.092  13.322   1.834\n"
        "   38A        P 4311  15.313  13.550   2.889\n"
        "   38A      OP1 4312  15.266  13.584   3.026\n"
        "   38A      OP2 4313  15.340  13.658   2.792\n"
        "   38A      O5' 4314  15.445  13.462   2.903\n"
        "   38A      C5' 4315  15.447  13.344   2.984\n"
        "   38A      C4' 4316  15.584  13.282   2.978\n"
        "   38A      O4' 4317  15.607  13.230   2.845\n"
        "   38A      C3' 4318  15.693  13.385   2.998\n"
        "   38A      O3' 4319  15.716  13.413   3.135\n"
        "   38A      C2' 4320  15.810  13.328   2.919\n"
        "   38A      O2' 4321  15.885  13.229   2.988\n"
        "   38A      C1' 4322  15.736  13.269   2.799\n"
        "   38A       N9 4323  15.716  13.370   2.695\n"
        "   38A       C8 4324  15.597  13.432   2.664\n"
        "   38A       N7 4325  15.606  13.519   2.566\n"
        "   38A       C5 4326  15.740  13.513   2.529\n"
        "   38A       C6 4327  15.814  13.580   2.429\n"
        "   38A       N6 4328  15.759  13.669   2.345\n"
        "   38A       N1 4329  15.946  13.553   2.419\n"
        "   38A       C2 4330  15.999  13.462   2.503\n"
        "   38A       N3 4331  15.940  13.391   2.600\n"
        "   38A       C4 4332  15.809  13.422   2.608\n"
        "122.00000 52.00000 36.00000\n"
        ;
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, str_from_cstr(text)));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    ASSERT_EQ((size_t)3, md_system_residue_count(&sys));
    EXPECT_TRUE(str_eq(md_system_residue_name(&sys, 1), STR_LIT("YYG")));
    for (size_t c = 0; c < 3; ++c) {
        EXPECT_EQ(MD_RESIDUE_KIND_NUCLEOTIDE, md_system_residue_kind(&sys, c));
        EXPECT_TRUE(md_system_residue_flags(&sys, c) & MD_RESIDUE_FLAG_RESOLVED);
    }
    ASSERT_EQ((size_t)1, md_system_part_count(&sys));
    EXPECT_EQ(MD_PART_KIND_RNA, md_system_part_kind(&sys, 0));
    ASSERT_EQ((size_t)1, sys.nucleic_backbone.range.count);
    EXPECT_EQ((size_t)3, sys.nucleic_backbone.segment.count);

    // The backbone of the YYG is its P O5' C5' C4' C3' O3', its base all of it but the phosphate and the sugar: the
    // tricyclic base and its side chain
    const md_urange_t range = md_system_residue_atom_range(&sys, 1);
    size_t base = 0, backbone = 0;
    for (uint32_t i = range.beg; i < range.end; ++i) {
        base     += (md_atom_flags(&sys.atom, i) & MD_ATOM_FLAG_NUCLEOBASE) != 0;
        backbone += (md_atom_flags(&sys.atom, i) & MD_ATOM_FLAG_BACKBONE) != 0;
    }
    EXPECT_EQ((size_t)27, base);
    EXPECT_EQ((size_t)6, backbone);
    md_vm_arena_destroy(alloc);
}

// The chains are the chain ids the data gives (the author's), whatever carries them, and nothing where it gives none
static size_t chain_atoms(const md_system_t* sys, size_t c) {
    return md_system_chain_atom_count(sys, c);
}

UTEST(util, chains) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    {
        // One chain of a polymer, its arginine and its 106 waters: they carry its id, they are of it
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1LAF.pdb")));
        ASSERT_EQ(1, (int)md_system_chain_count(&sys));
        EXPECT_TRUE(str_eq(md_system_chain_id(&sys, 0), STR_LIT("E")));
        EXPECT_EQ(md_system_atom_count(&sys), chain_atoms(&sys, 0));
        // The protein, the arginine and the waters are consecutive residues: one run
        EXPECT_EQ(1, (int)(md_system_chain_run_end(&sys, 0) - md_system_chain_run_beg(&sys, 0)));
        for (size_t c = 0; c < md_system_residue_count(&sys); ++c) EXPECT_EQ(0, md_system_chain_find_by_residue_idx(&sys, c));
        EXPECT_EQ(0, md_system_chain_find_by_atom_idx(&sys, md_system_atom_count(&sys) - 1));
    }

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb")));
        ASSERT_EQ(3, (int)md_system_chain_count(&sys));
        static const char* ids[] = {"A", "B", "C"};
        size_t total = 0;
        for (size_t c = 0; c < 3; ++c) {
            EXPECT_TRUE(str_eq(md_system_chain_id(&sys, c), str_from_cstr(ids[c])));
            total += chain_atoms(&sys, c);
        }
        EXPECT_EQ(md_system_atom_count(&sys), total);
    }

    {
        // mmCIF: the chains are the author's ids, not the asyms. Asyms E, F, G, K and L hold ligands, ions and waters
        // the author filed under the chains A and B.
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif")));
        ASSERT_EQ(2, (int)md_system_chain_count(&sys));
        EXPECT_TRUE(str_eq(md_system_chain_id(&sys, 0), STR_LIT("A")));
        EXPECT_TRUE(str_eq(md_system_chain_id(&sys, 1), STR_LIT("B")));
        EXPECT_EQ(md_system_atom_count(&sys), chain_atoms(&sys, 0) + chain_atoms(&sys, 1));
        for (size_t i = 0; i < md_system_instance_count(&sys); ++i) {
            const md_urange_t residues = md_system_instance_residue_range(&sys, i);
            for (uint32_t k = residues.beg; k < residues.end; ++k) {
                const md_chain_idx_t c = md_system_chain_find_by_residue_idx(&sys, k);
                ASSERT_GE(c, 0);
                EXPECT_TRUE(str_eq(md_system_chain_id(&sys, (size_t)c), md_system_instance_auth_id(&sys, i)));
            }
        }
        // Chain A is its protein; after the protein of B, its Mg; after the WO4 and Mg of B, its waters: three runs
        EXPECT_EQ(3, (int)(md_system_chain_run_end(&sys, 0) - md_system_chain_run_beg(&sys, 0)));
    }

    {
        // No ids in the data: no chains, and the strands of the polymers are parts
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro")));
        ASSERT_TRUE(md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL));
        EXPECT_EQ(0, (int)md_system_chain_count(&sys));
        size_t polymers = 0;
        for (size_t p = 0; p < md_system_part_count(&sys); ++p) {
            polymers += md_part_kind_is_polymer(md_system_part_kind(&sys, p));
            EXPECT_EQ(-1, md_system_chain_find_by_residue_idx(&sys, md_system_part_residue_range(&sys, p).beg));
        }
        EXPECT_EQ(253, (int)polymers);
    }

    {
        // A PDB file without chain ids (as simulations write them) is the same
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb")));
        EXPECT_EQ(0, (int)md_system_chain_count(&sys));
    }

    {
        // A topology names no chains either: the peptide is a polymer part, the waters and ions parts of their own
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/tpr/peptide_tip4p.tpr")));
        ASSERT_TRUE(md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL & ~MD_UTIL_INFER_BOND_BIT));
        EXPECT_EQ(0, (int)md_system_chain_count(&sys));
        EXPECT_EQ(MD_PART_KIND_PEPTIDE, md_system_part_kind(&sys, 0));
        EXPECT_EQ(MD_PART_KIND_WATER, md_system_part_kind(&sys, 1));
    }

    {
        // martinize-dna writes the two strands of a duplex as one molecule type: one molecule, as GROMACS has it, of
        // two parts, a strand each, as the .gro of the same system gives them
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini2_dna.tpr")));
        ASSERT_TRUE(md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL & ~MD_UTIL_INFER_BOND_BIT));
        EXPECT_EQ(0, (int)md_system_chain_count(&sys));
        const md_urange_t duplex = md_system_instance_residue_range(&sys, 0);
        EXPECT_EQ(24u, duplex.end - duplex.beg);
        EXPECT_EQ(0, md_system_part_find_by_residue_idx(&sys, duplex.beg));
        EXPECT_EQ(1, md_system_part_find_by_residue_idx(&sys, duplex.end - 1));
        EXPECT_EQ(12u, md_system_part_residue_range(&sys, 0).end - md_system_part_residue_range(&sys, 0).beg);
        const md_urange_t atoms = md_system_instance_atom_range(&sys, 0);
        const md_molecule_idx_t mol = md_system_molecule_find_by_atom_idx(&sys, atoms.beg);
        EXPECT_EQ(atoms.end - atoms.beg, md_system_molecule_atom_count(&sys, (size_t)mol));
        EXPECT_EQ(mol, md_system_molecule_find_by_atom_idx(&sys, atoms.end - 1));

        md_system_t gro = {.alloc = alloc};
        md_system_state_t gro_state = { .alloc = alloc };
        ASSERT_TRUE(md_gro_system_init_from_file(&gro, &gro_state, STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini2_dna.gro")));
        ASSERT_TRUE(md_util_system_infer(&gro, &gro_state, MD_UTIL_INFER_ALL));
        for (size_t p = 0; p < 2; ++p) {
            const md_urange_t a = md_system_part_atom_range(&sys, p), b = md_system_part_atom_range(&gro, p);
            EXPECT_EQ(a.end - a.beg, b.end - b.beg);
            EXPECT_TRUE(md_part_kind_is_nucleic_acid(md_system_part_kind(&gro, p)));
        }
    }

    md_vm_arena_destroy(alloc);
}

// A molecule of a topology is one molecule whatever its bonds, as GROMACS has it: here two waters made one molecule of
// the topology by hand, which nothing else would join
UTEST(util, topology_molecules_are_whole) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t sys_state = { .alloc = alloc };
    ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/tpr/peptide_tip3p.tpr")));
    ASSERT_TRUE(md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL & ~MD_UTIL_INFER_BOND_BIT));
    const size_t num_inst = md_system_instance_count(&sys);
    EXPECT_EQ(num_inst, md_system_molecule_count(&sys));
    const md_molecule_type_idx_t water = md_system_instance_molecule_type_idx(&sys, 1);
    ASSERT_GE(water, 0);
    ASSERT_EQ(water, md_system_instance_molecule_type_idx(&sys, 2));

    // Instances 1 and 2 as one molecule of the topology
    for (size_t i = 2; i < num_inst; ++i) {
        sys.instance.res_offset[i] = sys.instance.res_offset[i + 1];
        sys.instance.id[i] = sys.instance.id[i + 1 < num_inst ? i + 1 : i];
        sys.instance.auth_id[i] = sys.instance.auth_id[i + 1 < num_inst ? i + 1 : i];
        sys.instance.part_type_idx[i] = sys.instance.part_type_idx[i + 1 < num_inst ? i + 1 : i];
        sys.instance.mol_type_idx[i] = sys.instance.mol_type_idx[i + 1 < num_inst ? i + 1 : i];
    }
    sys.instance.count -= 1;
    ASSERT_TRUE(md_util_system_infer_molecules(&sys));
    EXPECT_EQ(num_inst - 1, md_system_molecule_count(&sys));
    check_molecule_invariants(utest_result, &sys);

    // Not a topology's: two molecules
    sys.instance.mol_type_idx[1] = -1;
    ASSERT_TRUE(md_util_system_infer_molecules(&sys));
    EXPECT_EQ(num_inst, md_system_molecule_count(&sys));
    md_vm_arena_destroy(alloc);
}

// Topologies: each molecule type is a molecule type of the system, each molecule an instance of it, and what the parts
// of the molecules are is inferred from the residues
UTEST(util, molecule_type_topology) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/tpr/peptide_tip4p.tpr")));

        static const struct { md_part_kind_t kind; const char* name; size_t count; } types[] = {
            {MD_PART_KIND_PEPTIDE, "Protein_chain_U", 1}, {MD_PART_KIND_WATER, "SOL", 506}, {MD_PART_KIND_NON_POLYMER, "NA", 2}, {MD_PART_KIND_NON_POLYMER, "CL", 2},
        };
        ASSERT_EQ(md_system_molecule_type_count(&sys), ARRAY_SIZE(types));
        size_t inst = 0;
        for (size_t t = 0; t < ARRAY_SIZE(types); ++t) {
            EXPECT_EQ(md_system_molecule_type_flags(&sys, t), MD_MOLECULE_TYPE_FLAG_NONE);
            EXPECT_TRUE(str_eq(md_system_molecule_type_name(&sys, t), str_from_cstr(types[t].name)));
            for (size_t k = 0; k < types[t].count; ++k, ++inst) {
                EXPECT_EQ(md_system_instance_molecule_type_idx(&sys, inst), (int)t);
                EXPECT_TRUE(str_empty(md_system_instance_id(&sys, inst)));
            }
        }
        ASSERT_EQ(md_system_instance_count(&sys), inst);
        EXPECT_EQ(md_system_instance_residue_count(&sys, 0), 2);

        // Inference leaves the topology's own alone: every molecule is of its type, and the parts are classified
        ASSERT_TRUE(md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL & ~MD_UTIL_INFER_BOND_BIT));
        EXPECT_EQ(md_system_molecule_type_count(&sys), ARRAY_SIZE(types));
        EXPECT_EQ(md_system_instance_count(&sys), inst);
        ASSERT_EQ(md_system_molecule_count(&sys), inst);
        inst = 0;
        for (size_t t = 0; t < ARRAY_SIZE(types); ++t) {
            for (size_t k = 0; k < types[t].count; ++k, ++inst) {
                const md_urange_t atoms = md_system_instance_atom_range(&sys, inst);
                const md_molecule_idx_t mol = md_system_molecule_find_by_atom_idx(&sys, atoms.beg);
                EXPECT_EQ(md_system_molecule_type_idx(&sys, (size_t)mol), (int)t);
                const md_part_idx_t part = md_system_part_find_by_atom_idx(&sys, atoms.beg);
                EXPECT_EQ(md_system_part_kind(&sys, (size_t)part), types[t].kind);
            }
        }
    }

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/tpr/martini3.tpr")));
        // The peptide, two lipids, four ions and 583 waters, as the structures are
        ASSERT_EQ(md_system_molecule_type_count(&sys), 5);
        EXPECT_EQ(md_system_instance_count(&sys), 1 + 2 + 4 + 583);
        ASSERT_TRUE(md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL & ~MD_UTIL_INFER_BOND_BIT));
        EXPECT_EQ(md_system_part_kind(&sys, (size_t)md_system_part_find_by_residue_idx(&sys, md_system_instance_residue_range(&sys, 0).beg)), MD_PART_KIND_PEPTIDE);
        EXPECT_EQ(md_system_part_kind(&sys, (size_t)md_system_part_find_by_residue_idx(&sys, md_system_instance_residue_range(&sys, 1).beg)), MD_PART_KIND_NON_POLYMER);
        EXPECT_EQ(md_system_part_kind(&sys, md_system_part_count(&sys) - 1), MD_PART_KIND_WATER);
        EXPECT_EQ(md_system_molecule_type_count(&sys), 5);
    }

    {
        // LAMMPS knows its molecules but not what they are: a molecule type per sequence of particle types, named by
        // its formula, and water told by its atoms
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        const char* atom_format = md_lammps_atom_format_strings()[MD_LAMMPS_ATOM_FORMAT_FULL];
        ASSERT_TRUE(md_lammps_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Cubic_Init.data"), atom_format));
        ASSERT_EQ(md_system_molecule_type_count(&sys), 2);
        EXPECT_TRUE(str_eq(md_system_molecule_type_name(&sys, 0), STR_LIT("C2H6")));
        EXPECT_TRUE(str_eq(md_system_molecule_type_name(&sys, 1), STR_LIT("H2O")));
        ASSERT_EQ(md_system_instance_count(&sys), md_system_residue_count(&sys));
        EXPECT_EQ(md_system_instance_molecule_type_idx(&sys, 0), 0);
        EXPECT_EQ(md_system_instance_molecule_type_idx(&sys, md_system_instance_count(&sys) - 1), 1);
        ASSERT_TRUE(md_util_system_infer_parts(&sys));
        EXPECT_EQ(md_system_part_kind(&sys, 0), MD_PART_KIND_NON_POLYMER);
        EXPECT_EQ(md_system_part_kind(&sys, md_system_part_count(&sys) - 1), MD_PART_KIND_WATER);
    }

    md_vm_arena_destroy(alloc);
}

// The molecules of the same make share a type, inferred where no topology gives one: named by their parts, or by their
// formula where they have none. Inferring them again gives the same types, not more of them.
UTEST(util, molecule_type_inferred) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = { .alloc = alloc };
        ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif")));
        ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
        EXPECT_VALID_SYSTEM(&sys);
        // The two chains of the hydrolase, the Mg, the WO4 and the waters
        const size_t num_types = md_system_molecule_type_count(&sys);
        ASSERT_EQ(num_types, 4u);
        for (size_t t = 0; t < num_types; ++t) {
            EXPECT_TRUE(md_system_molecule_type_flags(&sys, t) & MD_MOLECULE_TYPE_FLAG_INFERRED);
        }
        const md_molecule_type_idx_t a = md_system_molecule_type_idx(&sys, (size_t)md_system_molecule_find_by_atom_idx(&sys, md_system_instance_atom_range(&sys, 0).beg));
        const md_molecule_type_idx_t b = md_system_molecule_type_idx(&sys, (size_t)md_system_molecule_find_by_atom_idx(&sys, md_system_instance_atom_range(&sys, 1).beg));
        EXPECT_EQ(a, b);
        EXPECT_TRUE(str_eq(md_system_molecule_type_name(&sys, (size_t)a), STR_LIT("PHOSPHONOACETALDEHYDE HYDROLASE")));
        const md_molecule_type_idx_t w = md_system_molecule_type_idx(&sys, md_system_molecule_count(&sys) - 1);
        EXPECT_TRUE(str_eq(md_system_molecule_type_name(&sys, (size_t)w), STR_LIT("water")));

        ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
        EXPECT_EQ(md_system_molecule_type_count(&sys), num_types);
        EXPECT_VALID_SYSTEM(&sys);
    }

    {
        // No residues, no parts: the type is named by the formula
        md_system_t sys = {.alloc = alloc};
        md_system_state_t state = { .alloc = alloc };
        ASSERT_TRUE(md_xyz_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/C60-Ih.xyz"), MD_XYZ_OPTION_NONE));
        ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
        ASSERT_EQ(md_system_molecule_count(&sys), 1u);
        ASSERT_EQ(md_system_molecule_type_count(&sys), 1u);
        EXPECT_TRUE(str_eq(md_system_molecule_type_name(&sys, 0), STR_LIT("C60")));
    }

    md_vm_arena_destroy(alloc);
}

// A coarse grained strand numbered over a gap, without bonds, in a chain the file names: one part, and so one molecule
// (the hierarchy links follow the part, not the numbering), as every part is within one molecule
UTEST(util, coarse_grained_strand_over_a_gap) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_str_writer_t sb = md_str_writer_create(alloc);
    static const int resids[] = {1, 2, 3, 4, 5, 10, 11, 12, 13, 14};
    for (int i = 0; i < 10; ++i) {
        const float x = 3.8f * i + (i >= 5 ? 20.0f : 0.0f);
        md_str_writer_fmt(&sb, "ATOM  %5d  BB  ALA A%4d    %8.3f%8.3f%8.3f  1.00  0.00            \n", i + 1, resids[i], x, 0.0, 0.0);
    }
    md_str_writer_push_cstr(&sb, "END\n");
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &state, md_str_writer_to_str(&sb), MD_PDB_OPTION_NONE));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    EXPECT_TRUE(md_system_is_coarse_grained(&sys));
    EXPECT_VALID_SYSTEM(&sys);
    ASSERT_EQ(1u, md_system_part_count(&sys));
    EXPECT_EQ(MD_PART_KIND_PEPTIDE, md_system_part_kind(&sys, 0));
    EXPECT_EQ(1u, md_system_molecule_count(&sys));
    md_vm_arena_destroy(alloc);
}

// A topology supplementing a coarse grained GRO splits what the numbering ran together into its molecules: the parts and
// the backbones follow, each backbone range of one part
UTEST(util, topology_supplement_reinfers_backbones) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_str_writer_t sb = md_str_writer_create(alloc);
    md_str_writer_push_cstr(&sb, "cg\n12\n");
    for (int i = 0; i < 12; ++i) {
        const float x = 0.70f * (i % 6) + (i >= 6 ? 3.0f : 0.0f);
        const float y = i >= 6 ? 3.0f : 0.0f;
        md_str_writer_fmt(&sb, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", i + 1, "ALA", "BB", i + 1, x + 1.0f, y + 1.0f, 1.0f);
    }
    md_str_writer_push_cstr(&sb, "   8.00000   8.00000   8.00000\n");
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &state, md_str_writer_to_str(&sb)));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
    EXPECT_EQ(1u, md_system_part_count(&sys));

    md_str_writer_t top = md_str_writer_create(alloc);
    md_str_writer_push_cstr(&top, "[ defaults ]\n1 1\n[ atomtypes ]\nP5 72.0 0.0 A 0.47 5.0\n[ moleculetype ]\nPROT 1\n[ atoms ]\n");
    for (int i = 0; i < 6; ++i) md_str_writer_fmt(&top, "%d P5 %d ALA BB %d 0.0 72.0\n", i + 1, i + 1, i + 1);
    md_str_writer_push_cstr(&top, "[ bonds ]\n");
    for (int i = 0; i < 5; ++i) md_str_writer_fmt(&top, "%d %d 1 0.35 1250\n", i + 1, i + 2);
    md_str_writer_push_cstr(&top, "[ system ]\nx\n[ molecules ]\nPROT 2\n");
    md_itp_data_t data = {0};
    ASSERT_TRUE(md_itp_data_parse_str(&data, md_str_writer_to_str(&top), STR_LIT(""), alloc));
    ASSERT_TRUE(md_itp_system_supplement(&sys, &data));

    EXPECT_VALID_SYSTEM(&sys);
    ASSERT_EQ(2u, md_system_instance_count(&sys));
    EXPECT_EQ(2u, md_system_part_count(&sys));
    EXPECT_EQ(2u, md_system_molecule_count(&sys));
    EXPECT_EQ(1u, md_system_molecule_type_count(&sys));
    const md_protein_backbone_data_t* bb = &sys.protein_backbone;
    ASSERT_EQ(2u, bb->range.count);
    for (size_t r = 0; r < bb->range.count; ++r) {
        const md_part_idx_t part = bb->range.part_idx[r];
        ASSERT_GE(part, 0);
        const md_urange_t residues = md_system_part_residue_range(&sys, (size_t)part);
        for (uint32_t k = bb->range.offset[r]; k < bb->range.offset[r + 1]; ++k) {
            const md_residue_idx_t res = bb->segment.res_idx[k];
            EXPECT_TRUE(residues.beg <= (uint32_t)res && (uint32_t)res < residues.end);
        }
    }
    md_vm_arena_destroy(alloc);
}

// Every change of the topology gives the system a version no consumer has seen
UTEST(util, topology_version) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t st = {.alloc = alloc};
    EXPECT_EQ(sys.topology_version, 0u);

    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro")));
    const uint64_t loaded = sys.topology_version;
    EXPECT_NE(loaded, 0u);

    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    const uint64_t inferred = sys.topology_version;
    EXPECT_NE(inferred, loaded);

    // Inferring again derives the backbones from scratch rather than appending to them
    const size_t num_protein = sys.protein_backbone.segment.count, num_protein_ranges = sys.protein_backbone.range.count;
    const size_t num_nucleic = sys.nucleic_backbone.segment.count, num_nucleic_ranges = sys.nucleic_backbone.range.count;
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    EXPECT_NE(sys.topology_version, inferred);
    EXPECT_EQ(sys.protein_backbone.segment.count, num_protein);
    EXPECT_EQ(sys.protein_backbone.range.count, num_protein_ranges);
    EXPECT_EQ(sys.nucleic_backbone.segment.count, num_nucleic);
    EXPECT_EQ(sys.nucleic_backbone.range.count, num_nucleic_ranges);
    EXPECT_EQ(md_array_size(sys.protein_backbone.range.offset), num_protein_ranges + 1);

    // Bonds added and removed by hand
    const uint64_t before_insert = sys.topology_version;
    md_system_bond_insert(&sys, 0, 100, md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_USER));
    EXPECT_NE(sys.topology_version, before_insert);
    const uint64_t before_remove = sys.topology_version;
    md_system_bond_remove(&sys, (md_bond_idx_t)(sys.bond.count - 1));
    EXPECT_NE(sys.topology_version, before_remove);

    // Read only operations leave it alone
    const uint64_t read_only = sys.topology_version;
    md_bitmask_t mask = md_bitmask_create(0, (uint32_t)sys.atom.count, alloc);
    md_bitmask_set_bit(&mask, 0);
    md_util_mask_grow_by_bonds(&mask, &sys, 2, NULL);
    EXPECT_EQ(sys.topology_version, read_only);

    // A system loaded anew, into the same struct, has a version none of the above had
    const uint64_t last = sys.topology_version;
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro")));
    EXPECT_GT(sys.topology_version, last);

    md_system_free(&sys);
    EXPECT_EQ(sys.topology_version, 0u);
    md_vm_arena_destroy(alloc);
}

// The backbone angles and secondary structure of a frame are the state's, in its attributes
UTEST(util, state_backbone) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t st = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    const size_t n = sys.protein_backbone.segment.count;
    ASSERT_GT(n, 0);

    // Nothing until computed
    EXPECT_TRUE(md_util_state_backbone_angles(&st, &sys) == NULL);
    EXPECT_TRUE(md_util_state_secondary_structure(&st, &sys) == NULL);

    ASSERT_TRUE(md_util_state_backbone_compute(&st, &sys));
    const md_backbone_angles_t*     angle = md_util_state_backbone_angles(&st, &sys);
    const md_secondary_structure_t* ss    = md_util_state_secondary_structure(&st, &sys);
    ASSERT_TRUE(angle != NULL);
    ASSERT_TRUE(ss != NULL);

    // The same as computing them by hand
    md_backbone_angles_t*     ref_angle = md_alloc(alloc, n * sizeof(md_backbone_angles_t));
    md_secondary_structure_t* ref_ss    = md_alloc(alloc, n * sizeof(md_secondary_structure_t));
    md_util_backbone_angles_compute(ref_angle, n, st.xyz, &st.unitcell, &sys.protein_backbone);
    md_dssp_compute(ref_ss, n, st.xyz, &st.unitcell, &sys.protein_backbone);
    size_t structured = 0;
    for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(angle[i].phi, ref_angle[i].phi);
        EXPECT_EQ(angle[i].psi, ref_angle[i].psi);
        EXPECT_EQ(ss[i], ref_ss[i]);
        structured += ss[i] == MD_SECONDARY_STRUCTURE_HELIX_ALPHA || ss[i] == MD_SECONDARY_STRUCTURE_BETA_SHEET;
    }
    EXPECT_GT(structured, n / 4);

    // A producer writing every frame reuses the storage, and every write is a new version
    const md_attribute_t* attr = md_attributes_find(&st.attributes, STR_LIT(MD_BACKBONE_SECONDARY_STRUCTURE_PATH));
    ASSERT_TRUE(attr != NULL);
    const uint64_t version = attr->version;
    md_secondary_structure_t* dst = md_util_state_secondary_structure_write(&st, &sys);
    EXPECT_TRUE(dst == ss);
    EXPECT_GT(md_attributes_version(&st.attributes, attr->id), version);

    // A view has no table to write into
    md_system_state_t view = { .num_atoms = st.num_atoms, .xyz = st.xyz, .unitcell = st.unitcell };
    EXPECT_TRUE(md_util_state_secondary_structure_write(&view, &sys) == NULL);
    EXPECT_FALSE(md_util_state_backbone_compute(&view, &sys));

    // What a state carries for another backbone does not fit this one
    md_system_t other = {.alloc = alloc};
    md_system_state_t other_st = {.alloc = alloc};
    ASSERT_TRUE(md_gro_system_init_from_file(&other, &other_st, STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro")));
    ASSERT_TRUE(md_util_system_infer(&other, &other_st, MD_UTIL_INFER_ALL));
    EXPECT_TRUE(md_util_state_secondary_structure(&st, &other) == NULL);

    md_vm_arena_destroy(alloc);
}

// The binary searches from an atom to its residue and from a residue to its instance
UTEST(util, find_by_index) {
    const uint32_t offset[] = {0, 3, 3, 7, 10};   // [0,3) [3,3) [3,7) [7,10)
    EXPECT_EQ(md_offset_range_find(offset, 4, 0), 0);
    EXPECT_EQ(md_offset_range_find(offset, 4, 2), 0);
    EXPECT_EQ(md_offset_range_find(offset, 4, 3), 2);   // The empty range holds nothing
    EXPECT_EQ(md_offset_range_find(offset, 4, 6), 2);
    EXPECT_EQ(md_offset_range_find(offset, 4, 7), 3);
    EXPECT_EQ(md_offset_range_find(offset, 4, 9), 3);
    EXPECT_EQ(md_offset_range_find(offset, 4, 10), -1);
    EXPECT_EQ(md_offset_range_find(offset, 0, 0), -1);
    EXPECT_EQ(md_offset_range_find(NULL, 4, 0), -1);

    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t sys_state = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro")));
    md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);
    for (size_t c = 0; c < sys.residue.count; ++c) {
        const md_urange_t range = md_system_residue_atom_range(&sys, c);
        EXPECT_EQ(md_system_residue_find_by_atom_idx(&sys, range.beg), (int)c);
        EXPECT_EQ(md_system_residue_find_by_atom_idx(&sys, range.end - 1), (int)c);
    }
    for (size_t i = 0; i < sys.instance.count; ++i) {
        const md_urange_t range = md_system_instance_atom_range(&sys, i);
        EXPECT_EQ(md_system_instance_find_by_atom_idx(&sys, range.beg), (int)i);
        EXPECT_EQ(md_system_instance_find_by_atom_idx(&sys, range.end - 1), (int)i);
    }
    EXPECT_EQ(md_system_residue_find_by_atom_idx(&sys, sys.atom.count), -1);
    md_vm_arena_destroy(alloc);
}

// The walk keeps a per atom depth in temp memory, where zero means 'not reached'. Temp memory is not cleared on
// allocation, so a walk which followed some other use of temp memory used to see stale depths and stop early.
UTEST_F(util, mask_grow_by_bonds_dirty_temp) {
    const md_system_t* sys = &utest_fixture->mol_ala;
    md_bitmask_t mask = md_bitmask_create(0, (uint32_t)sys->atom.count, utest_fixture->alloc);
    md_bitmask_t ref  = md_bitmask_create(0, (uint32_t)sys->atom.count, utest_fixture->alloc);

    // Stale depths of 1 make every neighbour look as if it had already been reached at a lower depth
    const size_t junk_size = MEGABYTES(1);
    {
        md_temp_scope_t temp = md_temp_begin();
        memset(md_temp_alloc(temp, junk_size), 1, junk_size);
        md_temp_end(temp);
    }
    md_bitmask_set_bit(&mask, 0);
    md_util_mask_grow_by_bonds(&mask, sys, 3, NULL);

    {
        md_temp_scope_t temp = md_temp_begin();
        memset(md_temp_alloc(temp, junk_size), 0, junk_size);
        md_temp_end(temp);
    }
    md_bitmask_set_bit(&ref, 0);
    md_util_mask_grow_by_bonds(&ref, sys, 3, NULL);

    EXPECT_GT(md_bitmask_count(&ref), (size_t)1);
    EXPECT_EQ(md_bitmask_count(&ref), md_bitmask_count(&mask));
}

// An atom that crosses the cell boundary between two frames goes the short way, across the
// boundary, and not back through the box. The orthorhombic branch once discarded its
// minimum image and so did exactly that.
UTEST(util, interpolate_linear_across_boundary) {
    const md_unitcell_t cell = md_unitcell_from_extent(10, 10, 10);
    vec3_t a[16] = {0}, b[16] = {0}, out[16] = {0};
    for (int i = 0; i < 16; ++i) {
        a[i] = vec3_set(9.5f, 5.0f, 0.25f);
        b[i] = vec3_set(0.5f, 5.0f, 9.75f);
    }
    const vec3_t* const in[2] = { a, b };
    ASSERT_TRUE(md_util_interpolate_linear(out, in, 16, &cell, 0.5f));
    for (int i = 0; i < 16; ++i) {
        // Halfway is on the boundary: 10 or its image 0, never 5
        EXPECT_NEAR(0.0f, fabsf(fmodf(out[i].x + 5.0f, 10.0f) - 5.0f), 1.0e-4f);
        EXPECT_NEAR(5.0f, out[i].y, 1.0e-4f);
        EXPECT_NEAR(0.0f, fabsf(fmodf(out[i].z + 5.0f, 10.0f) - 5.0f), 1.0e-4f);
    }

    // Without a cell it is a plain blend
    ASSERT_TRUE(md_util_interpolate_linear(out, in, 16, NULL, 0.5f));
    EXPECT_NEAR(5.0f, out[3].x, 1.0e-4f);
}

// Wrapping a subset moves those atoms and no others; the orthorhombic branch once wrote each result
// to the i:th atom instead of the one indexed.
UTEST(util, pbc_indexed) {
    const md_unitcell_t cell = md_unitcell_from_extent(10, 10, 10);
    vec3_t xyz[20];
    for (int i = 0; i < 20; ++i) xyz[i] = vec3_set(12.0f + i, -3.0f, 5.0f);
    const int32_t idx[3] = { 17, 4, 11 };
    ASSERT_TRUE(md_util_pbc(xyz, idx, 3, &cell));
    for (int i = 0; i < 20; ++i) {
        const bool wrapped = (i == 17 || i == 4 || i == 11);
        const float ex = wrapped ? fmodf(12.0f + i, 10.0f) : 12.0f + i;
        EXPECT_NEAR(ex, xyz[i].x, 1.0e-4f);
        EXPECT_NEAR(wrapped ? 7.0f : -3.0f, xyz[i].y, 1.0e-4f);
    }

    // All of them, through the eight wide path and its tail. An atom on the boundary may come out
    // at either end, which is the same place.
    ASSERT_TRUE(md_util_pbc(xyz, NULL, 20, &cell));
    for (int i = 0; i < 20; ++i) {
        EXPECT_GE(xyz[i].x, 0.0f);
        EXPECT_LE(xyz[i].x, 10.0f);
        EXPECT_NEAR(0.0f, remainderf(xyz[i].x - (12.0f + i), 10.0f), 1.0e-4f);
        EXPECT_NEAR(7.0f, xyz[i].y, 1.0e-4f);
        EXPECT_NEAR(5.0f, xyz[i].z, 1.0e-4f);
    }

    // Triclinic: every wrapped atom lands inside the cell
    const md_unitcell_t tri = md_unitcell_from_basis_parameters(10, 10, 10, 2, 1, 3);
    for (int i = 0; i < 20; ++i) xyz[i] = vec3_set(12.0f + i, -3.0f - i, 25.0f);
    ASSERT_TRUE(md_util_pbc(xyz, NULL, 20, &tri));
    mat3_t I = {0};
    md_unitcell_I_extract_float(I.elem, &tri);
    for (int i = 0; i < 20; ++i) {
        const vec3_t f = mat3_mul_vec3(I, xyz[i]);
        for (int k = 0; k < 3; ++k) {
            EXPECT_GE(f.elem[k], -1.0e-4f);
            EXPECT_LE(f.elem[k], 1.0f + 1.0e-4f);
        }
    }
}

// Each of the four paths (plain, radius, index, both) against a scalar reference, on a count that
// leaves a tail after the eight wide part.
UTEST(util, aabb_paths) {
    enum { N = 37 };
    vec3_t xyz[N];
    float r[N];
    int32_t idx[N];
    for (int i = 0; i < N; ++i) {
        const float s = (float)((i * 2654435761u) % 1000) / 100.0f - 5.0f;
        xyz[i] = vec3_set(s, -2.0f * s + i, 0.5f * i - s);
        r[i]   = 0.1f * (i % 7);
        idx[i] = (i * 11) % N;
    }
    for (int variant = 0; variant < 4; ++variant) {
        const float*   rr = (variant & 1) ? r   : NULL;
        const int32_t* ii = (variant & 2) ? idx : NULL;
        const size_t   n  = ii ? 29 : N;
        float ref_min[3] = { FLT_MAX,  FLT_MAX,  FLT_MAX};
        float ref_max[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
        for (size_t k = 0; k < n; ++k) {
            const int32_t j = ii ? ii[k] : (int32_t)k;
            const float rad = rr ? rr[j] : 0.0f;
            for (int c = 0; c < 3; ++c) {
                ref_min[c] = MIN(ref_min[c], xyz[j].elem[c] - rad);
                ref_max[c] = MAX(ref_max[c], xyz[j].elem[c] + rad);
            }
        }
        float got_min[3], got_max[3];
        md_util_aabb_compute(got_min, got_max, xyz, rr, ii, n);
        for (int c = 0; c < 3; ++c) {
            EXPECT_EQ(ref_min[c], got_min[c]);
            EXPECT_EQ(ref_max[c], got_max[c]);
        }
    }
}

// The spline passes through the two middle frames at t = 0 and t = 1, and across a boundary it
// takes the neighbouring frames' images rather than their wrapped positions.
UTEST(util, interpolate_cubic) {
    vec3_t f[4][16], out[16];
    for (int i = 0; i < 16; ++i) {
        for (int k = 0; k < 4; ++k) {
            f[k][i] = vec3_set(1.0f * k + i, 2.0f * k, -0.5f * k * k);
        }
    }
    const vec3_t* const in[4] = { f[0], f[1], f[2], f[3] };
    ASSERT_TRUE(md_util_interpolate_cubic_spline(out, in, 16, NULL, 0.0f, 0.5f));
    for (int i = 0; i < 16; ++i) {
        EXPECT_NEAR(f[1][i].x, out[i].x, 1.0e-4f);
        EXPECT_NEAR(f[1][i].z, out[i].z, 1.0e-4f);
    }
    ASSERT_TRUE(md_util_interpolate_cubic_spline(out, in, 16, NULL, 1.0f, 0.5f));
    for (int i = 0; i < 16; ++i) {
        EXPECT_NEAR(f[2][i].y, out[i].y, 1.0e-4f);
    }

    // Moving +1 per frame along x, wrapped into a box of 10
    const md_unitcell_t cell = md_unitcell_from_extent(10, 10, 10);
    for (int i = 0; i < 16; ++i) {
        for (int k = 0; k < 4; ++k) {
            f[k][i] = vec3_set(fmodf(8.0f + k, 10.0f), 5.0f, 5.0f);   // 8, 9, 0, 1
        }
    }
    ASSERT_TRUE(md_util_interpolate_cubic_spline(out, in, 16, &cell, 0.5f, 0.5f));
    for (int i = 0; i < 16; ++i) {
        EXPECT_NEAR(9.5f, out[i].x, 1.0e-3f);
    }
}

// Both schemes' weights sum to one and reproduce uniform motion. The cardinal spline lands on the middle frames exactly,
// through the weighted sum too; the B-spline is short of them by design, at (p0 + 4 p1 + p2) / 6. Across a periodic
// boundary, orthorhombic or triclinic, the weighted sum takes the frames' images.
UTEST(util, interpolate_cubic_weights) {
    for (int i = 0; i <= 16; ++i) {
        const float t = i / 16.0f;
        float wc[4], wb[4];
        md_util_cubic_spline_weights(wc, t, 0.5f);
        md_util_cubic_bspline_weights(wb, t);
        EXPECT_NEAR(1.0f, wc[0] + wc[1] + wc[2] + wc[3], 1.0e-6f);
        EXPECT_NEAR(1.0f, wb[0] + wb[1] + wb[2] + wb[3], 1.0e-6f);
        // Frames at 0, 1, 2 and 3: the interval between the middle two is at 1 + t
        EXPECT_NEAR(1.0f + t, wc[1] + 2.0f * wc[2] + 3.0f * wc[3], 1.0e-6f);
        EXPECT_NEAR(1.0f + t, wb[1] + 2.0f * wb[2] + 3.0f * wb[3], 1.0e-6f);
    }

    float w[4];
    md_util_cubic_spline_weights(w, 0.0f, 0.5f);
    EXPECT_EQ(0.0f, fabsf(w[0]));
    EXPECT_EQ(1.0f, w[1]);
    EXPECT_EQ(0.0f, fabsf(w[2]));
    EXPECT_EQ(0.0f, fabsf(w[3]));
    md_util_cubic_spline_weights(w, 1.0f, 0.5f);
    EXPECT_EQ(0.0f, fabsf(w[0]));
    EXPECT_EQ(0.0f, fabsf(w[1]));
    EXPECT_EQ(1.0f, w[2]);
    EXPECT_EQ(0.0f, fabsf(w[3]));
    md_util_cubic_bspline_weights(w, 0.0f);
    EXPECT_NEAR(1.0f / 6.0f, w[0], 1.0e-7f);
    EXPECT_NEAR(4.0f / 6.0f, w[1], 1.0e-7f);
    EXPECT_NEAR(1.0f / 6.0f, w[2], 1.0e-7f);
    EXPECT_NEAR(0.0f,        w[3], 1.0e-7f);

    vec3_t f[4][16], out[16];
    for (int i = 0; i < 16; ++i) {
        for (int k = 0; k < 4; ++k) {
            f[k][i] = vec3_set(1.0f * k + i, 2.0f * k * k, -0.5f * k + 0.1f * i);
        }
    }
    const vec3_t* const in[4] = { f[0], f[1], f[2], f[3] };

    md_util_cubic_spline_weights(w, 0.0f, 0.5f);
    ASSERT_TRUE(md_util_interpolate_cubic_weighted(out, in, 16, NULL, w));
    for (int i = 0; i < 16; ++i) {
        EXPECT_EQ(f[1][i].x, out[i].x);
        EXPECT_EQ(f[1][i].y, out[i].y);
        EXPECT_EQ(f[1][i].z, out[i].z);
    }

    md_util_cubic_bspline_weights(w, 0.0f);
    ASSERT_TRUE(md_util_interpolate_cubic_weighted(out, in, 16, NULL, w));
    for (int i = 0; i < 16; ++i) {
        for (int c = 0; c < 3; ++c) {
            const float ref = (f[0][i].elem[c] + 4.0f * f[1][i].elem[c] + f[2][i].elem[c]) / 6.0f;
            EXPECT_NEAR(ref, out[i].elem[c], 1.0e-4f);
        }
    }

    // Moving +1 per frame along x, wrapped into a box of 10: 8, 9, 0, 1
    for (int i = 0; i < 16; ++i) {
        for (int k = 0; k < 4; ++k) {
            f[k][i] = vec3_set(fmodf(8.0f + k, 10.0f), 5.0f, 5.0f);
        }
    }
    const md_unitcell_t cells[2] = {
        md_unitcell_from_extent(10, 10, 10),
        md_unitcell_from_basis_parameters(10, 10, 10, 2, 0, 0),
    };
    for (int c = 0; c < 2; ++c) {
        md_util_cubic_bspline_weights(w, 0.5f);
        ASSERT_TRUE(md_util_interpolate_cubic_weighted(out, in, 16, &cells[c], w));
        for (int i = 0; i < 16; ++i) {
            EXPECT_NEAR(9.5f, out[i].x, 1.0e-3f);
            EXPECT_NEAR(5.0f, out[i].y, 1.0e-3f);
        }
    }
}

// The order is a permutation of the atoms, and atoms that share a position end up together.
UTEST(util, sort_spatial) {
    enum { N = 64 };
    vec3_t xyz[N];
    for (int i = 0; i < N; ++i) {
        const int cluster = i % 4;
        xyz[i] = vec3_set(cluster * 20.0f, cluster * 20.0f, 0.0f);
    }
    uint32_t order[N];
    md_util_sort_spatial(order, xyz, N);
    bool seen[N] = {0};
    for (int i = 0; i < N; ++i) {
        ASSERT_LT(order[i], (uint32_t)N);
        EXPECT_FALSE(seen[order[i]]);
        seen[order[i]] = true;
    }
    int changes = 0;
    for (int i = 1; i < N; ++i) {
        changes += (order[i] % 4) != (order[i - 1] % 4);
    }
    EXPECT_EQ(3, changes);
}

// The component wise kernels step over a few atoms to align their loads, so where an array begins
// decides how they split the work. Every start, against scalar references.
UTEST(util, packed_kernels_any_alignment) {
    enum { N = 37, PAD = 64 };
    const md_unitcell_t cell = md_unitcell_from_extent(10, 12, 14);
    const vec3_t ext = {10, 12, 14};
    vec3_t base_a[N + PAD], base_b[N + PAD], base_o[N + PAD];
    float  base_w[N + PAD];
    for (int off = 0; off < 8; ++off) {
        vec3_t* a = base_a + off;
        vec3_t* b = base_b + off;
        vec3_t* o = base_o + off;
        float*  w = base_w + off;
        for (int i = 0; i < N; ++i) {
            const float s = (float)((i * 2654435761u + off) % 1000) / 1000.0f;
            a[i] = vec3_set(9.0f * s + 0.5f, 11.0f * (1.0f - s) + 0.5f, 13.0f * s * s + 0.5f);
            b[i] = vec3_add(a[i], vec3_set(3.0f - 6.0f * s, 5.0f * s - 2.5f, 1.5f));
            w[i] = 0.5f + s;
        }

        // Plain and weighted COM without a cell
        double ref[3] = {0}, ref_w[3] = {0}, sw = 0;
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < 3; ++k) {
                ref[k]   += a[i].elem[k];
                ref_w[k] += a[i].elem[k] * w[i];
            }
            sw += w[i];
        }
        const vec3_t c  = md_util_com_compute(a, NULL, NULL, N, NULL);
        const vec3_t cw = md_util_com_compute(a, w,    NULL, N, NULL);
        for (int k = 0; k < 3; ++k) {
            EXPECT_NEAR(ref[k] / N, c.elem[k], 1.0e-4);
            EXPECT_NEAR(ref_w[k] / sw, cw.elem[k], 1.0e-4);
        }

        // Periodic COM of a cloud well inside the cell is its plain mean, weighted or not
        vec3_t near_center[N];
        for (int i = 0; i < N; ++i) near_center[i] = vec3_add(vec3_set(5, 6, 7), vec3_mul1(vec3_sub(a[i], vec3_set(5, 6, 7)), 0.1f));
        MEMCPY(o, near_center, sizeof(near_center));
        const vec3_t pc  = md_util_com_compute(o, NULL, NULL, N, &cell);
        const vec3_t pcw = md_util_com_compute(o, w,    NULL, N, &cell);
        const vec3_t mc  = md_util_com_compute(o, NULL, NULL, N, NULL);
        const vec3_t mcw = md_util_com_compute(o, w,    NULL, N, NULL);
        for (int k = 0; k < 3; ++k) {
            EXPECT_NEAR(mc.elem[k],  pc.elem[k],  2.0e-2);
            EXPECT_NEAR(mcw.elem[k], pcw.elem[k], 2.0e-2);
        }

        // AABB
        float mn[3], mx[3];
        md_util_aabb_compute(mn, mx, a, NULL, NULL, N);
        for (int k = 0; k < 3; ++k) {
            float lo = FLT_MAX, hi = -FLT_MAX;
            for (int i = 0; i < N; ++i) { lo = MIN(lo, a[i].elem[k]); hi = MAX(hi, a[i].elem[k]); }
            EXPECT_EQ(lo, mn[k]);
            EXPECT_EQ(hi, mx[k]);
        }

        // Linear interpolation, no cell and orthorhombic, and nothing written past the end
        base_o[off + N] = vec3_set(-1, -1, -1);
        const vec3_t* const in[2] = { a, b };
        ASSERT_TRUE(md_util_interpolate_linear(o, in, N, NULL, 0.25f));
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < 3; ++k) {
                EXPECT_NEAR(a[i].elem[k] + 0.25f * (b[i].elem[k] - a[i].elem[k]), o[i].elem[k], 1.0e-4f);
            }
        }
        ASSERT_TRUE(md_util_interpolate_linear(o, in, N, &cell, 0.25f));
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < 3; ++k) {
                float d = b[i].elem[k] - a[i].elem[k];
                d -= ext.elem[k] * roundf(d / ext.elem[k]);
                EXPECT_NEAR(a[i].elem[k] + 0.25f * d, o[i].elem[k], 1.0e-4f);
            }
        }
        EXPECT_EQ(-1.0f, base_o[off + N].x);

        // Cubic through the middle frames at the ends of the interval
        const vec3_t* const in4[4] = { a, a, b, b };
        ASSERT_TRUE(md_util_interpolate_cubic_spline(o, in4, N, NULL, 1.0f, 0.5f));
        for (int i = 0; i < N; ++i) {
            EXPECT_NEAR(b[i].y, o[i].y, 1.0e-4f);
        }
        EXPECT_EQ(-1.0f, base_o[off + N].x);

        // Wrapping
        MEMCPY(o, b, N * sizeof(vec3_t));
        ASSERT_TRUE(md_util_pbc(o, NULL, N, &cell));
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < 3; ++k) {
                EXPECT_GE(o[i].elem[k], -1.0e-4f);
                EXPECT_LE(o[i].elem[k], ext.elem[k] + 1.0e-4f);
                EXPECT_NEAR(0.0f, remainderf(o[i].elem[k] - b[i].elem[k], ext.elem[k]), 1.0e-4f);
            }
        }
    }
}

// ---- minimum and maximum distance between sets ----------------------------------------------

// The minimum image by exhaustion, in double: reduced into the brick, then every combination of up to three lattice
// vectors either way along the periodic axes. Three is more than any of the cells below needs: a vector A n which
// shortens a reduced d has |n_i| <= 2 |d| |r_i| (r_i a row of the inverse basis), at most 2 for these.
static double dist_ref(vec3_t a, vec3_t b, const md_unitcell_t* cell) {
    double d[3] = { (double)a.x - b.x, (double)a.y - b.y, (double)a.z - b.z };
    double A[3][3] = {{0}};
    md_unitcell_A_extract_double(A, cell);
    int per[3] = { (cell->flags & MD_UNITCELL_PBC_X) != 0, (cell->flags & MD_UNITCELL_PBC_Y) != 0, (cell->flags & MD_UNITCELL_PBC_Z) != 0 };
    for (int i = 0; i < 3; ++i) {
        if (!per[i] || A[i][i] == 0.0) { A[i][0] = A[i][1] = A[i][2] = 0.0; per[i] = 0; }
    }
    for (int i = 2; i >= 0; --i) {
        if (per[i]) {
            const double n = nearbyint(d[i] / A[i][i]);
            for (int j = 0; j < 3; ++j) d[j] -= n * A[i][j];
        }
    }
    const int R = 3;
    double best = DBL_MAX;
    for (int nz = per[2] ? -R : 0; nz <= (per[2] ? R : 0); ++nz) {
        for (int ny = per[1] ? -R : 0; ny <= (per[1] ? R : 0); ++ny) {
            for (int nx = per[0] ? -R : 0; nx <= (per[0] ? R : 0); ++nx) {
                double v[3];
                for (int j = 0; j < 3; ++j) v[j] = d[j] + nx * A[0][j] + ny * A[1][j] + nz * A[2][j];
                best = MIN(best, v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            }
        }
    }
    return sqrt(best);
}

static uint64_t dist_rng_state;
static double dist_rand(void) {
    uint64_t x = dist_rng_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    dist_rng_state = x;
    return (double)(x >> 11) * (1.0 / 9007199254740992.0);
}

// A point at fractional coordinates in [lo, hi) along the periodic axes - several images of the cell, as in an
// unwrapped trajectory - and within span along the others
static vec3_t dist_rand_point(const md_unitcell_t* cell, double span, double lo, double hi) {
    double A[3][3] = {{0}};
    md_unitcell_A_extract_double(A, cell);
    double p[3] = {0};
    for (int i = 0; i < 3; ++i) {
        const double f = lo + (hi - lo) * dist_rand();
        if (A[i][i] != 0.0) {
            for (int j = 0; j < 3; ++j) p[j] += f * A[i][j];
        } else {
            p[i] += f * span;
        }
    }
    return vec3_set((float)p[0], (float)p[1], (float)p[2]);
}

// Groups which are empty, compact, spread over several images, or share points with b, in every kind of cell -
// including cells small enough that most distances are beyond half of them, where reducing a displacement into the
// cell is not the minimum image, and a skewed cell which is not reduced - against the minimum image by exhaustion
UTEST(util, min_max_distance_groups) {
    const double s2 = sqrt(2.0), s6 = sqrt(6.0), d = 30.0, ds = 11.0;
    const md_unitcell_t cells[] = {
        md_unitcell_none(),
        md_unitcell_from_extent(30, 30, 30),
        md_unitcell_from_extent(9, 12, 10),
        md_unitcell_from_basis_parameters(30, 40, 0, 0, 0, 0),                              // periodic in x and y only
        md_unitcell_from_extent_and_angles(30, 30, 30, 70, 70, 70),
        md_unitcell_from_basis_parameters(d, d, d / s2, 0, d / 2, d / 2),                   // rhombic dodecahedron
        md_unitcell_from_basis_parameters(d, 2 * s2 * d / 3, s6 * d / 3, d / 3, -d / 3, s2 * d / 3),  // truncated octahedron
        md_unitcell_from_basis_parameters(ds, ds, ds / s2, 0, ds / 2, ds / 2),              // small dodecahedron
        md_unitcell_from_basis_parameters(30, 20, 25, 27, 20, 15),                          // skewed, not reduced
        md_unitcell_from_basis_parameters(30, 30, 0, 10, 0, 0),                             // triclinic, periodic in x and y only
    };

    dist_rng_state = 0x243F6A8885A308D3ULL;
    vec3_t  a[6 * 12];
    vec3_t  b[60];
    size_t  off[7];
    float   dist[6];
    int64_t ia[6], ib[6];

    for (size_t ci = 0; ci < ARRAY_SIZE(cells); ++ci) {
        const md_unitcell_t* cell = &cells[ci];
        for (int trial = 0; trial < 25; ++trial) {
            const size_t nb = trial == 0 ? 0 : (trial % 5 == 0 ? 1 + (size_t)(dist_rand() * 3) : 1 + (size_t)(dist_rand() * 59));
            const bool blob = trial % 3 == 0;
            const vec3_t bc = dist_rand_point(cell, 40, 0, 1);
            for (size_t j = 0; j < nb; ++j) {
                b[j] = blob ? vec3_add(bc, vec3_set((float)(dist_rand() * 8 - 4), (float)(dist_rand() * 8 - 4), (float)(dist_rand() * 8 - 4)))
                            : dist_rand_point(cell, 40, -1.5, 2.5);
            }

            const size_t ng = 1 + (size_t)(dist_rand() * 6);
            size_t na = 0;
            for (size_t g = 0; g < ng; ++g) {
                off[g] = na;
                const size_t n    = dist_rand() < 0.2 ? 0 : 1 + (size_t)(dist_rand() * 12);
                const int    kind = (int)(dist_rand() * 4);
                const vec3_t cen  = dist_rand_point(cell, 40, -1, 2);
                for (size_t k = 0; k < n; ++k) {
                    if (kind == 0) {
                        a[na++] = dist_rand_point(cell, 40, -1.5, 2.5);
                    } else if (kind == 1 && nb > 0) {
                        a[na++] = b[(size_t)(dist_rand() * nb)];
                    } else {
                        a[na++] = vec3_add(cen, vec3_set((float)(dist_rand() * 5 - 2.5), (float)(dist_rand() * 5 - 2.5), (float)(dist_rand() * 5 - 2.5)));
                    }
                }
            }
            off[ng] = na;

            for (int largest = 0; largest < 2; ++largest) {
                if (largest) {
                    md_util_max_distance_groups(dist, ia, ib, a, off, ng, b, nb, cell);
                } else {
                    md_util_min_distance_groups(dist, ia, ib, a, off, ng, b, nb, cell);
                }
                for (size_t g = 0; g < ng; ++g) {
                    if (off[g] == off[g + 1] || nb == 0) {
                        EXPECT_EQ(0.0f, dist[g]);
                        EXPECT_EQ(-1, ia[g]);
                        EXPECT_EQ(-1, ib[g]);
                        continue;
                    }
                    double ref = largest ? -1.0 : DBL_MAX;
                    for (size_t i = off[g]; i < off[g + 1]; ++i) {
                        for (size_t j = 0; j < nb; ++j) {
                            const double r = dist_ref(a[i], b[j], cell);
                            ref = largest ? MAX(ref, r) : MIN(ref, r);
                        }
                    }
                    EXPECT_NEAR(ref, dist[g], 2.0e-3);
                    // The pair is one which has that distance
                    ASSERT_GE(ia[g], (int64_t)off[g]);
                    ASSERT_LT(ia[g], (int64_t)off[g + 1]);
                    ASSERT_GE(ib[g], 0);
                    ASSERT_LT(ib[g], (int64_t)nb);
                    EXPECT_NEAR(ref, dist_ref(a[ia[g]], b[ib[g]], cell), 2.0e-3);
                }
            }
        }
    }

    // The single set versions: FLT_MAX and 0 for an empty set, the indices untouched then
    const vec3_t p = vec3_set(1, 2, 3);
    const md_unitcell_t cell = md_unitcell_from_extent(10, 10, 10);
    int64_t i0 = 7, i1 = 7;
    EXPECT_EQ(FLT_MAX, md_util_min_distance(&i0, &i1, &p, 0, &p, 1, &cell));
    EXPECT_EQ(0.0f,    md_util_max_distance(&i0, &i1, &p, 1, &p, 0, &cell));
    EXPECT_EQ(7, i0);
    EXPECT_EQ(7, i1);
    const vec3_t q = vec3_set(9.5f, 2, 3);
    EXPECT_NEAR(1.5f, md_util_min_distance(&i0, &i1, &p, 1, &q, 1, &cell), 1.0e-5f);
    EXPECT_NEAR(8.5f, md_util_min_distance(&i0, &i1, &p, 1, &q, 1, NULL), 1.0e-5f);
    EXPECT_EQ(0, i0);
    EXPECT_EQ(0, i1);
}

UTEST(util, ramachandran_classify) {
    // Two chains, both starting with a proline, the second right after a residue which is not one:
    //   chain 0: PRO GLY PRO ALA PRO PRO SER
    //   chain 1: PRO ALA
    const char* names[9] = { "PRO", "GLY", "PRO", "ALA", "PRO", "PRO", "SER", "PRO", "ALA" };
    md_label_t labels[9];
    md_residue_idx_t res_idx[9];
    for (int i = 0; i < 9; ++i) {
        labels[i]   = make_label(str_from_cstr(names[i]));
        res_idx[i] = i;
    }
    uint32_t range_offset[3] = { 0, 7, 9 };

    md_system_t sys = {0};
    sys.residue.count = 9;
    sys.residue.name    = labels;
    sys.protein_backbone.range.count      = 2;
    sys.protein_backbone.range.offset     = range_offset;
    sys.protein_backbone.segment.count    = 9;
    sys.protein_backbone.segment.res_idx = res_idx;

    // MolProbity's types: a glycine or a proline before a proline keeps its own type, and the last residue of a chain
    // is not pre-proline because of the first residue of the next
    enum {
        G = MD_RAMACHANDRAN_TYPE_GENERAL,
        Y = MD_RAMACHANDRAN_TYPE_GLYCINE,
        P = MD_RAMACHANDRAN_TYPE_PROLINE,
        R = MD_RAMACHANDRAN_TYPE_PREPROL,
    };
    const md_ramachandran_type_t expected[9] = { P, Y, P, R, P, P, G, P, G };

    // A guard on each side of the buffer catches a write outside it
    md_ramachandran_type_t buf[1 + 9 + 1];
    MEMSET(buf, 0xAB, sizeof(buf));
    EXPECT_TRUE(md_util_backbone_ramachandran_classify(buf + 1, 9, &sys));
    EXPECT_EQ(0xAB, buf[0]);
    EXPECT_EQ(0xAB, buf[10]);
    for (int i = 0; i < 9; ++i) EXPECT_EQ(expected[i], buf[1 + i]);

    // A capacity short of the segments: the last type written is still decided by the segment after it
    MEMSET(buf, 0xAB, sizeof(buf));
    EXPECT_TRUE(md_util_backbone_ramachandran_classify(buf + 1, 4, &sys));
    EXPECT_EQ(0xAB, buf[0]);
    for (int i = 0; i < 4; ++i) EXPECT_EQ(expected[i], buf[1 + i]);
    EXPECT_EQ(0xAB, buf[5]);
}

// BACKBONES FROM THE BONDS (md_util.c). A file read with its residue and atom names, and the same file as an XYZ file
// would give it, the element and the position of each atom and nothing else, have the same backbones: the same atoms
// in each segment, consecutive along the same chains, and so the same secondary structure and ramachandran types.

// The XYZ text of a system, its atoms in their order, or by element (as VASP and ase.build.sort write them). Atom k of
// the text is atom xyz_to_sys[k] of the system; atoms without an element are left out. Returns the number of atoms.
static size_t xyz_text_of(str_t* out, int32_t* xyz_to_sys, const md_system_t* sys, const md_system_state_t* state, bool by_element, md_allocator_i* alloc) {
    size_t num = 0;
    for (int z = 1; z < 256; ++z) {
        for (size_t i = 0; i < sys->atom.count; ++i) {
            const md_atomic_number_t zi = md_atom_atomic_number(&sys->atom, i);
            if (zi != 0 && (!by_element || zi == z)) xyz_to_sys[num++] = (int32_t)i;
        }
        if (!by_element) break;
    }
    char* buf = md_alloc(alloc, 64 * (num + 1));
    size_t len = (size_t)snprintf(buf, 64, "%zu\nxyz\n", num);
    for (size_t k = 0; k < num; ++k) {
        const str_t sym = md_atomic_number_symbol(md_atom_atomic_number(&sys->atom, (size_t)xyz_to_sys[k]));
        const vec3_t p = md_state_coord(state, (size_t)xyz_to_sys[k]);
        len += (size_t)snprintf(buf + len, 64, "%.*s %.4f %.4f %.4f\n", (int)sym.len, sym.ptr, p.x, p.y, p.z);
    }
    *out = (str_t){ buf, len };
    return num;
}

// The segment of each atom which is the CA (C4') of one, -1 for none
static int32_t* segment_of_atom(size_t num_atoms, md_allocator_i* alloc) {
    int32_t* seg = md_alloc(alloc, sizeof(int32_t) * num_atoms);
    MEMSET(seg, -1, sizeof(int32_t) * num_atoms);
    return seg;
}

static void expect_backbones_from_bonds(int* utest_result, const char* path, bool by_element, bool compare_secondary_structure) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t ref = { .alloc = alloc };
    md_system_state_t ref_state = { .alloc = alloc };
    const str_t file = str_from_cstr(path);
    bool loaded = false;
    if (str_ends_with(file, STR_LIT(".pdb"))) loaded = md_pdb_system_init_from_file(&ref, &ref_state, file, MD_PDB_OPTION_NONE);
    if (str_ends_with(file, STR_LIT(".gro"))) loaded = md_gro_system_init_from_file(&ref, &ref_state, file);
    if (str_ends_with(file, STR_LIT(".cif"))) loaded = md_mmcif_system_init_from_file(&ref, &ref_state, file);
    ASSERT_TRUE(loaded);
    ASSERT_TRUE(md_util_system_infer(&ref, &ref_state, MD_UTIL_INFER_ALL));

    int32_t* xyz_to_ref = md_alloc(alloc, sizeof(int32_t) * ref.atom.count);
    str_t text;
    const size_t num_xyz = xyz_text_of(&text, xyz_to_ref, &ref, &ref_state, by_element, alloc);
    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };
    ASSERT_TRUE(md_xyz_system_init_from_str(&sys, &state, text, MD_XYZ_OPTION_NONE));
    ASSERT_EQ(num_xyz, sys.atom.count);
    ASSERT_EQ(0, sys.residue.count);
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));

    // An atom of the reference as an atom of the XYZ
    int32_t* to_xyz = md_alloc(alloc, sizeof(int32_t) * ref.atom.count);
    MEMSET(to_xyz, -1, sizeof(int32_t) * ref.atom.count);
    for (size_t k = 0; k < num_xyz; ++k) to_xyz[xyz_to_ref[k]] = (int32_t)k;
#define X(a) ((a) < 0 ? -1 : to_xyz[a])
#define SEG(a) (X(a) < 0 ? -1 : seg[X(a)])

    {
        const md_protein_backbone_data_t* rb = &ref.protein_backbone;
        const md_protein_backbone_data_t* bb = &sys.protein_backbone;
        EXPECT_EQ(rb->range.count, bb->range.count);
        EXPECT_EQ(rb->segment.count, bb->segment.count);
        for (size_t s = 0; s < bb->segment.count; ++s) EXPECT_EQ(-1, bb->segment.res_idx[s]);
        for (size_t r = 0; r < bb->range.count; ++r) EXPECT_EQ(-1, bb->range.part_idx[r]);

        int32_t* seg = segment_of_atom(sys.atom.count, alloc);
        for (size_t s = 0; s < bb->segment.count; ++s) seg[bb->segment.atoms[s].ca] = (int32_t)s;
        size_t same = 0;
        for (size_t r = 0; r < rb->range.count; ++r) {
            for (uint32_t j = rb->range.offset[r]; j < rb->range.offset[r + 1]; ++j) {
                const md_amino_acid_atoms_t a = rb->segment.atoms[j];
                const int32_t s = SEG(a.ca);
                if (s < 0) continue;
                const md_amino_acid_atoms_t b = bb->segment.atoms[s];
                // Either oxygen of a C terminal carboxylate
                const bool o = b.o == X(a.o) || (X(a.o) >= 0 && md_bond_find(&sys.bond, b.c, b.o) >= 0 && md_atom_atomic_number(&sys.atom, (size_t)b.o) == MD_Z_O);
                // Consecutive along the chain
                const bool chain = j == rb->range.offset[r] || (s > 0 && SEG(rb->segment.atoms[j - 1].ca) == s - 1);
                same += b.n == X(a.n) && b.c == X(a.c) && o && b.cb == X(a.cb) && chain;
            }
        }
        EXPECT_EQ(rb->segment.count, same);

        if (rb->segment.count && rb->segment.count == bb->segment.count) {
            md_ramachandran_type_t* rama_ref = md_alloc(alloc, sizeof(md_ramachandran_type_t) * rb->segment.count);
            md_ramachandran_type_t* rama = md_alloc(alloc, sizeof(md_ramachandran_type_t) * bb->segment.count);
            EXPECT_TRUE(md_util_backbone_ramachandran_classify(rama_ref, rb->segment.count, &ref));
            EXPECT_TRUE(md_util_backbone_ramachandran_classify(rama, bb->segment.count, &sys));
            md_secondary_structure_t* ss_ref = md_alloc(alloc, sizeof(md_secondary_structure_t) * rb->segment.count);
            md_secondary_structure_t* ss = md_alloc(alloc, sizeof(md_secondary_structure_t) * bb->segment.count);
            if (compare_secondary_structure) {
                EXPECT_TRUE(md_dssp_compute(ss_ref, rb->segment.count, ref_state.xyz, &ref_state.unitcell, rb));
                EXPECT_TRUE(md_dssp_compute(ss, bb->segment.count, state.xyz, &state.unitcell, bb));
            }
            size_t same_rama = 0, same_ss = 0;
            for (size_t j = 0; j < rb->segment.count; ++j) {
                const int32_t s = SEG(rb->segment.atoms[j].ca);
                if (s < 0) continue;
                // A residue whose side chain the crystal did not resolve, not even its CB, is a glycine by its bonds
                // (LYS 156 of 9gfu)
                const bool truncated = rb->segment.atoms[j].cb < 0 && rama_ref[j] != MD_RAMACHANDRAN_TYPE_GLYCINE;
                same_rama += rama[s] == rama_ref[j] || (truncated && rama[s] == MD_RAMACHANDRAN_TYPE_GLYCINE);
                same_ss += !compare_secondary_structure || MEMCMP(&ss[s], &ss_ref[j], sizeof(md_secondary_structure_t)) == 0;
            }
            EXPECT_EQ(rb->segment.count, same_rama);
            EXPECT_EQ(rb->segment.count, same_ss);
        }
    }

    {
        const md_nucleic_backbone_data_t* rb = &ref.nucleic_backbone;
        const md_nucleic_backbone_data_t* bb = &sys.nucleic_backbone;
        EXPECT_EQ(rb->range.count, bb->range.count);
        EXPECT_EQ(rb->segment.count, bb->segment.count);
        for (size_t s = 0; s < bb->segment.count; ++s) EXPECT_EQ(-1, bb->segment.res_idx[s]);

        int32_t* seg = segment_of_atom(sys.atom.count, alloc);
        for (size_t s = 0; s < bb->segment.count; ++s) seg[bb->segment.atoms[s].c4] = (int32_t)s;
        size_t same = 0;
        for (size_t r = 0; r < rb->range.count; ++r) {
            for (uint32_t j = rb->range.offset[r]; j < rb->range.offset[r + 1]; ++j) {
                const md_nucleic_acid_atoms_t a = rb->segment.atoms[j];
                const int32_t s = SEG(a.c4);
                if (s < 0) continue;
                const md_nucleic_acid_atoms_t b = bb->segment.atoms[s];
                const bool chain = j == rb->range.offset[r] || (s > 0 && SEG(rb->segment.atoms[j - 1].c4) == s - 1);
                // The P on the O5' of a 5' triphosphate, which the names do not tell (PA of the GTP which begins a
                // chain of 8g7u)
                const bool p = b.p == X(a.p) || (a.p < 0 && b.p >= 0 && md_bond_find(&sys.bond, b.o5, b.p) >= 0);
                same += p && b.o5 == X(a.o5) && b.c5 == X(a.c5) && b.c3 == X(a.c3) && b.o3 == X(a.o3) &&
                        b.c1 == X(a.c1) && b.c2 == X(a.c2) && b.o4 == X(a.o4) && chain;
            }
        }
        EXPECT_EQ(rb->segment.count, same);
    }
#undef SEG
#undef X

    md_vm_arena_destroy(alloc);
}

// Crystal structures, without hydrogens. In 1k4r the backbone O of HIS 86 of chain C is 1.4 Å from the ND1 of HIS 86 of
// chain A, close enough to be taken for a bond: the residue is still found, and its chain not broken there.
UTEST(util, backbones_from_bonds_protein) {
    const char* files[] = {
        MD_UNITTEST_DATA_DIR "/1k4r.pdb",
        MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb",
        MD_UNITTEST_DATA_DIR "/1fez.cif",
        MD_UNITTEST_DATA_DIR "/2or2.cif",
        MD_UNITTEST_DATA_DIR "/9gfu.cif",
    };
    for (size_t i = 0; i < ARRAY_SIZE(files); ++i) {
        expect_backbones_from_bonds(utest_result, files[i], false, true);
    }
}

// Proteins with DNA (a nucleosome, with hydrogens) and with RNA. With hydrogens the HN of a residue is its hydrogen
// rather than the one the names leave out, which changes the secondary structure DSSP finds now and then.
UTEST(util, backbones_from_bonds_nucleic) {
    expect_backbones_from_bonds(utest_result, MD_UNITTEST_DATA_DIR "/nucl-dna.gro", false, false);
    expect_backbones_from_bonds(utest_result, MD_UNITTEST_DATA_DIR "/8g7u.cif", false, true);
    expect_backbones_from_bonds(utest_result, MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb", false, false);
}

// The order of the atoms in the file plays no part
UTEST(util, backbones_from_bonds_by_element) {
    expect_backbones_from_bonds(utest_result, MD_UNITTEST_DATA_DIR "/nucl-dna.gro", true, false);
    expect_backbones_from_bonds(utest_result, MD_UNITTEST_DATA_DIR "/1k4r.pdb", true, true);
}

// Nothing which is no polymer: fullerenes, and a tryptophan alone (one residue is no chain)
UTEST(util, backbones_from_bonds_none) {
    const char* files[] = {
        MD_UNITTEST_DATA_DIR "/c720.xyz",
        MD_UNITTEST_DATA_DIR "/tryptophan.xyz",
    };
    for (size_t i = 0; i < ARRAY_SIZE(files); ++i) {
        md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
        md_system_t sys = { .alloc = alloc };
        md_system_state_t state = { .alloc = alloc };
        ASSERT_TRUE(md_xyz_system_init_from_file(&sys, &state, str_from_cstr(files[i]), MD_XYZ_OPTION_NONE));
        ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL));
        EXPECT_TRUE(sys.bond.count > 0);
        EXPECT_EQ(0, sys.protein_backbone.range.count);
        EXPECT_EQ(0, sys.protein_backbone.segment.count);
        EXPECT_EQ(0, sys.nucleic_backbone.range.count);
        md_vm_arena_destroy(alloc);
    }
}

UTEST(util, index_data_cleared) {
    // A cleared index data keeps its arrays, emptied: no ranges, as one never pushed to, and none in range
    md_allocator_i* alloc = md_vm_arena_create(MEGABYTES(1));
    md_index_data_t data = { .alloc = alloc };
    EXPECT_EQ(0u, md_index_data_num_ranges(&data));
    const int32_t a[] = { 1, 2, 3 };
    md_index_data_push_arr(&data, a, 3);
    md_index_data_push_arr(&data, a, 2);
    EXPECT_EQ(2u, md_index_data_num_ranges(&data));
    md_index_data_clear(&data);
    EXPECT_EQ(0u, md_index_data_num_ranges(&data));
    EXPECT_TRUE(md_index_range_beg(&data, 0) == NULL);
    EXPECT_EQ(0u, md_index_range_size(&data, 0));
    md_index_data_push_arr(&data, a, 1);
    EXPECT_EQ(1u, md_index_data_num_ranges(&data));
    EXPECT_EQ(1u, md_index_range_size(&data, 0));
    md_vm_arena_destroy(alloc);
}

// The kernels that take packed positions give the same bits wherever the positions start in memory. Each runs on the
// same values copied to every offset of 4 bytes within 64, the alignments a gather or an allocator can hand it. A
// sum whose grouping followed the address (scalar up to an alignment, vector lanes after) fails here.
static double addr_rnd(double lo, double hi) { return lo + (hi - lo) * (double)rand() / (double)RAND_MAX; }

UTEST(util, packed_kernels_independent_of_address) {
    enum { N = 103, OFFSETS = 16 };
    srand(4711);
    vec3_t src[N];
    float w[N];
    for (int i = 0; i < N; ++i) {
        src[i] = (vec3_t){{(float)addr_rnd(-5.0, 45.0), (float)addr_rnd(-5.0, 45.0), (float)addr_rnd(-5.0, 45.0)}};
        w[i] = (float)addr_rnd(1.0, 16.0);
    }
    const md_unitcell_t cell = md_unitcell_from_extent(40.0, 41.0, 42.0);
    const md_unitcell_t none = {0};
    // 64-byte aligned storage with room to slide the positions along it
    float* raw = (float*)md_alloc(md_get_heap_allocator(), (6 * N + 64) * sizeof(float) + 64);
    float* base = (float*)(((uintptr_t)raw + 63) & ~(uintptr_t)63);

    vec3_t ref_com[4] = {0};
    float ref_min[3], ref_max[3];
    vec3_t ref_pbc[N], ref_lin[N];
    for (int o = 0; o < OFFSETS; ++o) {
        vec3_t* p = (vec3_t*)(base + o);
        memcpy(p, src, sizeof(src));
        const vec3_t com[4] = {
            md_util_com_compute(p, NULL, NULL, N, &none),
            md_util_com_compute(p, w, NULL, N, &none),
            md_util_com_compute(p, NULL, NULL, N, &cell),
            md_util_com_compute(p, w, NULL, N, &cell),
        };
        float mn[3], mx[3];
        md_util_aabb_compute(mn, mx, p, NULL, NULL, N);
        // pbc and interpolation write in place or to out: the result at this offset against the first offset's
        memcpy(p, src, sizeof(src));
        md_util_pbc(p, NULL, N, &cell);
        vec3_t pbc[N];
        memcpy(pbc, p, sizeof(pbc));
        memcpy(p, src, sizeof(src));
        vec3_t* q = (vec3_t*)(base + 3 * N + 32 - o);   // a second frame, at another alignment
        for (int i = 0; i < N; ++i) q[i] = vec3_add(src[i], (vec3_t){{0.5f, -0.25f, 0.125f}});
        const vec3_t* frames[2] = {p, q};
        vec3_t lin[N];
        md_util_interpolate_linear(lin, frames, N, &cell, 0.37f);
        if (o == 0) {
            memcpy(ref_com, com, sizeof(com));
            memcpy(ref_min, mn, sizeof(mn));
            memcpy(ref_max, mx, sizeof(mx));
            memcpy(ref_pbc, pbc, sizeof(pbc));
            memcpy(ref_lin, lin, sizeof(lin));
            continue;
        }
        for (int k = 0; k < 4; ++k) {
            if (memcmp(&ref_com[k], &com[k], sizeof(vec3_t))) printf("com %d differs at offset %d bytes: %.9g %.9g %.9g against %.9g %.9g %.9g\n",
                k, 4 * o, com[k].x, com[k].y, com[k].z, ref_com[k].x, ref_com[k].y, ref_com[k].z);
            EXPECT_EQ(0, memcmp(&ref_com[k], &com[k], sizeof(vec3_t)));
        }
        EXPECT_EQ(0, memcmp(ref_min, mn, sizeof(mn)));
        EXPECT_EQ(0, memcmp(ref_max, mx, sizeof(mx)));
        EXPECT_EQ(0, memcmp(ref_pbc, pbc, sizeof(pbc)));
        EXPECT_EQ(0, memcmp(ref_lin, lin, sizeof(lin)));
    }
    md_free(md_get_heap_allocator(), raw, (6 * N + 64) * sizeof(float) + 64);
}
