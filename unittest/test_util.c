#include "utest.h"
#include <string.h>
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

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_str_builder.h>
#include <core/md_bitfield.h>
#include <core/md_hash.h>

#include "rmsd.h"

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

// Appends a component of coarse grained beads. The BB bead of an amino acid component is its backbone.
static void cg_add_component(md_system_t* sys, const char* comp_name, int seq_id, md_component_kind_t comp_kind, const char* atom_names[], size_t num_atoms, md_allocator_i* alloc) {
    if (sys->atom.type.count == 0) {
        // Type 0 is the "not found" sentinel
        md_atom_type_find_or_add(&sys->atom.type, STR_LIT("?"), 0, 0, 0, 0, 0, alloc);
    }
    if (sys->component.count == 0) {
        md_array_push(sys->component.atom_offset, 0, alloc);
    }
    for (size_t i = 0; i < num_atoms; ++i) {
        const md_atom_type_flags_t type_flags = md_atom_type_flags_set_particle_kind(MD_ATOM_TYPE_FLAG_NONE, MD_PARTICLE_BEAD);
        const md_atom_flags_t flags = (comp_kind == MD_COMPONENT_KIND_AMINO_ACID && strcmp(atom_names[i], "BB") == 0) ? MD_ATOM_FLAG_BACKBONE : MD_ATOM_FLAG_NONE;
        md_atom_type_idx_t type = md_atom_type_find_or_add(&sys->atom.type, str_from_cstr(atom_names[i]), 0, 50.0f, 2.35f, 0xFFFFFFFF, type_flags, alloc);
        md_array_push(sys->atom.type_idx, type, alloc);
        md_array_push(sys->atom.flags, flags, alloc);
        sys->atom.count += 1;
    }
    md_array_push(sys->component.name, make_label(str_from_cstr(comp_name)), alloc);
    md_array_push(sys->component.seq_id, seq_id, alloc);
    md_array_push(sys->component.flags, md_component_flags_set_kind(MD_COMPONENT_FLAG_NONE, comp_kind), alloc);
    md_array_push(sys->component.atom_offset, (uint32_t)sys->atom.count, alloc);
    sys->component.count += 1;
}

static void check_structure_invariants(int* utest_result, const md_system_t* sys) {
    size_t total = 0;
    for (size_t s = 0; s < md_structure_count(&sys->structure); ++s) {
        md_structure_t structure = {0};
        md_structure_extract(&structure, &sys->structure, s);
        total += structure.count;
        EXPECT_EQ(structure.parent_idx[0], structure.atom_idx[0]);
        for (size_t k = 0; k < structure.count; ++k) {
            const int32_t atom   = structure.atom_idx[k];
            const int32_t parent = structure.parent_idx[k];
            EXPECT_EQ(md_structure_atom_slot(&sys->structure, atom), (int32_t)(sys->structure.offset[s] + k));
            if (k > 0) {
                EXPECT_NE(parent, atom);
                EXPECT_LT(md_structure_atom_slot(&sys->structure, parent), md_structure_atom_slot(&sys->structure, atom));
            }
        }
    }
    EXPECT_EQ(total, sys->atom.count);
}

// Coarse grained systems frequently carry no bonds. The structure hierarchy must then come from the
// component hierarchy: beads hang off their component's anchor, consecutive polymer components are
// chained, non polymer components stay separate, and nothing is written to sys->bond.
UTEST(util, structure_hierarchy_coarse_grained_without_bonds) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    md_system_t sys = {0};
    sys.alloc = alloc;

    const char* ala[] = {"BB"};
    const char* lys[] = {"BB", "SC1", "SC2"};
    const char* phe[] = {"BB", "SC1", "SC2", "SC3"};
    const char* popc[] = {"NC3", "PO4", "GL1", "GL2", "C1A", "C1B"};
    const char* w[] = {"W"};

    cg_add_component(&sys, "LYS", 1, MD_COMPONENT_KIND_AMINO_ACID, lys, 3, alloc); // atoms 0-2
    cg_add_component(&sys, "ALA", 2, MD_COMPONENT_KIND_AMINO_ACID, ala, 1, alloc); // atom  3
    cg_add_component(&sys, "PHE", 3, MD_COMPONENT_KIND_AMINO_ACID, phe, 4, alloc); // atoms 4-7
    cg_add_component(&sys, "ALA", 7, MD_COMPONENT_KIND_AMINO_ACID, ala, 1, alloc); // atom  8, sequence gap: new chain
    cg_add_component(&sys, "ALA", 8, MD_COMPONENT_KIND_AMINO_ACID, ala, 1, alloc); // atom  9
    cg_add_component(&sys, "POPC", 9, 0, popc, 6, alloc);                // atoms 10-15
    cg_add_component(&sys, "POPC", 10, 0, popc, 6, alloc);               // atoms 16-21, must not join the previous lipid
    cg_add_component(&sys, "W", 11, MD_COMPONENT_KIND_WATER, w, 1, alloc);         // atom  22

    ASSERT_TRUE(md_util_system_infer_structures(&sys));
    EXPECT_EQ(sys.bond.count, 0u);

    ASSERT_EQ(md_structure_count(&sys.structure), 5u);
    check_structure_invariants(utest_result, &sys);

    // The first chain is a path of three anchors, so the middle one is the root
    EXPECT_EQ(sys.structure.offset[1] - sys.structure.offset[0], 8u);
    EXPECT_EQ(sys.structure.atom_idx[0], 3);
    EXPECT_EQ(md_structure_atom_parent(&sys.structure, 0), 3);
    EXPECT_EQ(md_structure_atom_parent(&sys.structure, 1), 0);
    EXPECT_EQ(md_structure_atom_parent(&sys.structure, 2), 0);
    EXPECT_EQ(md_structure_atom_parent(&sys.structure, 4), 3);
    EXPECT_EQ(md_structure_atom_parent(&sys.structure, 7), 4);

    // Second chain, broken off by the sequence gap
    EXPECT_EQ(sys.structure.offset[2] - sys.structure.offset[1], 2u);

    // Each lipid on its own, rooted at its first bead, and the water bead alone
    EXPECT_EQ(sys.structure.offset[3] - sys.structure.offset[2], 6u);
    EXPECT_EQ(sys.structure.offset[4] - sys.structure.offset[3], 6u);
    EXPECT_EQ(sys.structure.offset[5] - sys.structure.offset[4], 1u);
    EXPECT_EQ(md_structure_atom_parent(&sys.structure, 15), 10);
    EXPECT_EQ(md_structure_atom_parent(&sys.structure, 21), 16);

    md_vm_arena_destroy(alloc);
}

// Where bonds exist they are used, and links only fill in what they leave disconnected. A system that is
// not coarse grained gets no links at all.
UTEST(util, structure_hierarchy_coarse_grained_partial_bonds) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    {
        md_system_t sys = {0};
        sys.alloc = alloc;
        const char* lys[] = {"BB", "SC1", "SC2"};
        cg_add_component(&sys, "LYS", 1, MD_COMPONENT_KIND_AMINO_ACID, lys, 3, alloc);
        cg_add_component(&sys, "LYS", 2, MD_COMPONENT_KIND_AMINO_ACID, lys, 3, alloc);

        // SC1-SC2 bonded within each residue, nothing else
        md_atom_pair_t p0 = {{1, 2}};
        md_atom_pair_t p1 = {{4, 5}};
        md_array_push(sys.bond.pairs, p0, alloc);
        md_array_push(sys.bond.pairs, p1, alloc);
        sys.bond.count = 2;
        md_bond_build_connectivity(&sys.bond, sys.atom.count, alloc);

        ASSERT_TRUE(md_util_system_infer_structures(&sys));
        EXPECT_EQ(sys.bond.count, 2u);
        ASSERT_EQ(md_structure_count(&sys.structure), 1u);
        check_structure_invariants(utest_result, &sys);

        // SC2 is reached through its bond to SC1, not through a redundant link to BB
        EXPECT_EQ(md_structure_atom_parent(&sys.structure, 2), 1);
        EXPECT_EQ(md_structure_atom_parent(&sys.structure, 5), 4);
    }

    {
        md_system_t sys = {0};
        sys.alloc = alloc;
        const char* lys[] = {"BB", "SC1", "SC2"};
        cg_add_component(&sys, "LYS", 1, MD_COMPONENT_KIND_AMINO_ACID, lys, 3, alloc);
        cg_add_component(&sys, "LYS", 2, MD_COMPONENT_KIND_AMINO_ACID, lys, 3, alloc);
        for (size_t i = 0; i < sys.atom.type.count; ++i) {
            sys.atom.type.flags[i] = md_atom_type_flags_set_particle_kind(sys.atom.type.flags[i], MD_PARTICLE_ATOM);
        }

        ASSERT_TRUE(md_util_system_infer_structures(&sys));
        EXPECT_EQ(md_structure_count(&sys.structure), 6u);
    }

    md_vm_arena_destroy(alloc);
}

// A bondless coarse grained chain wrapped into a small box must come out whole. The BB beads sit 3.8A
// apart along x in a 10A box and each carries a side chain bead 3A up in y, across the cell boundary.
UTEST(util, unwrap_structure_coarse_grained_without_bonds) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    md_system_t sys = {0};
    sys.alloc = alloc;
    const char* res[] = {"BB", "SC1"};
    enum { N = 7 };
    for (int i = 0; i < N; ++i) {
        cg_add_component(&sys, "LEU", i + 1, MD_COMPONENT_KIND_AMINO_ACID, res, 2, alloc);
    }
    ASSERT_TRUE(md_util_system_infer_structures(&sys));
    ASSERT_EQ(md_structure_count(&sys.structure), 1u);

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

// The coarse grained branch of covalent bond inference used to skip the first component, test every bead
// against itself (always under the cutoff, so a self bond) and emit every other pair twice.
UTEST(util, infer_bonds_coarse_grained) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    md_system_t sys = {0};
    sys.alloc = alloc;
    const char* lys[] = {"BB", "SC1", "SC2"};
    cg_add_component(&sys, "LYS", 1, MD_COMPONENT_KIND_AMINO_ACID, lys, 3, alloc);
    cg_add_component(&sys, "LYS", 2, MD_COMPONENT_KIND_AMINO_ACID, lys, 3, alloc);

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

// Bond inference marks what it produces as of origin MD_BOND_ORIGIN_INFERRED, and re-inference replaces only those.
// A bond of another origin (read from a file, user defined, from a topology) is kept wherever it sits in the arrays,
// is never duplicated by an inferred bond naming the same pair, and survives an inference that cannot run.
UTEST(util, infer_bonds_replaces_only_inferred) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    // Four carbons 1.5A apart along x: inference finds 0-1, 1-2, 2-3
    md_system_t sys = { .alloc = alloc };
    md_atom_type_find_or_add(&sys.atom.type, STR_LIT("?"), 0, 0, 0, 0, 0, alloc);
    md_atom_type_idx_t c = md_atom_type_find_or_add(&sys.atom.type, STR_LIT("C"), 6, 12.011f, 0.76f, 0, 0, alloc);
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
    vec3_t xyz[9];
    md_system_state_t state = { .num_atoms = 9, .xyz = pack_xyz(xyz, x, y, z, 9), .unitcell = md_unitcell_from_extent(10, 10, 10) };

    md_structure_t structure = {0};
    ASSERT_TRUE(md_structure_extract(&structure, &sys.structure, 0));
    md_util_unwrap_structure(&state, &structure);
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
    vec3_t xyz[9];
    md_system_state_t state = { .num_atoms = 9, .xyz = pack_xyz(xyz, x, y, z, 9), .unitcell = cell };

    md_structure_t structure = {0};
    ASSERT_TRUE(md_structure_extract(&structure, &sys.structure, 0));
    md_util_unwrap_structure(&state, &structure);
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

    md_tick_t t0, t1;

    t0 = md_tick_now();
    md_util_sort_radix_uint32(indices, values, N);
    t1 = md_tick_now();

    printf("Time for radix index sort: %.4f ms\n", md_tick_to_milliseconds(t1 - t0));

    t0 = md_tick_now();
    md_util_sort_radix_inplace_uint32(values, N);
    t1 = md_tick_now();

    printf("Time for radix inplace sort: %.4f ms\n", md_tick_to_milliseconds(t1 - t0));
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

static md_entity_flags_t inferred_kind(md_entity_kind_t kind) {
    return md_entity_flags_set_kind(MD_ENTITY_FLAG_INFERRED, kind);
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
        EXPECT_EQ(md_system_entity_flags(&sys, 0), inferred_kind(MD_ENTITY_KIND_PEPTIDE));
        EXPECT_TRUE(str_eq(md_entity_description(&sys.entity, 0), STR_LIT("peptide")));

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
        EXPECT_EQ(md_system_entity_flags(&sys, 0), inferred_kind(MD_ENTITY_KIND_PEPTIDE));

        ASSERT_EQ(md_system_instance_count(&sys), 3);
        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 0), STR_LIT("A")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 0), STR_LIT("A")));

        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 1), STR_LIT("B")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 1), STR_LIT("B")));

        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 2), STR_LIT("C")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 2), STR_LIT("C")));
    }

    {
        // One chain and its 106 waters, each water an instance of its own which shares the id of the others
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1LAF.pdb")));
        EXPECT_GT(md_system_atom_count(&sys),   0);
        ASSERT_EQ(md_system_entity_count(&sys), 2);
        EXPECT_EQ(md_system_entity_flags(&sys, 0), inferred_kind(MD_ENTITY_KIND_PEPTIDE));
        EXPECT_EQ(md_system_entity_flags(&sys, 1), inferred_kind(MD_ENTITY_KIND_WATER));

        ASSERT_EQ(md_system_instance_count(&sys), 1 + 106);
        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 0), STR_LIT("A")));
        EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, 0), STR_LIT("E")));
        EXPECT_EQ(md_system_instance_comp_count(&sys, 0), 239);
        for (size_t i = 1; i < md_system_instance_count(&sys); ++i) {
            EXPECT_TRUE(str_eq(md_system_instance_id(&sys, i), STR_LIT("B")));
            EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, i), STR_LIT("E")));
            EXPECT_EQ(md_system_instance_entity_idx(&sys, i), 1);
            EXPECT_EQ(md_system_instance_comp_count(&sys, i), 1);
        }
    }

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb")));
        EXPECT_GT(md_system_atom_count(&sys),   0);

        // GTP and GDP are named like nucleotides atom by atom, but are not linked into a chain: ligands
        size_t num_nucleotides = 0;
        for (size_t i = 0; i < md_system_component_count(&sys); ++i) {
            num_nucleotides += md_system_component_kind(&sys, i) == MD_COMPONENT_KIND_NUCLEOTIDE;
        }
        EXPECT_EQ(num_nucleotides, 0);

        static const struct { md_entity_kind_t kind; const char* desc; } entities[] = {
            {MD_ENTITY_KIND_PEPTIDE, "peptide"}, {MD_ENTITY_KIND_PEPTIDE, "peptide"}, {MD_ENTITY_KIND_NON_POLYMER, "GTP"}, {MD_ENTITY_KIND_NON_POLYMER, "MG"},
            {MD_ENTITY_KIND_NON_POLYMER, "SO4"}, {MD_ENTITY_KIND_WATER, "water"}, {MD_ENTITY_KIND_NON_POLYMER, "GDP"}, {MD_ENTITY_KIND_NON_POLYMER, "VLB"},
        };
        ASSERT_EQ(md_system_entity_count(&sys), ARRAY_SIZE(entities));
        for (size_t i = 0; i < ARRAY_SIZE(entities); ++i) {
            EXPECT_EQ(md_system_entity_flags(&sys, i), inferred_kind(entities[i].kind));
            EXPECT_TRUE(str_eq(md_entity_description(&sys.entity, i), str_from_cstr(entities[i].desc)));
        }

        // Chains, then the molecules of chain A, of chain B and of chain C. The two SO4 and the five waters of chain A
        // are instances of their own, sharing an id within the chain.
        static const struct { const char* id; const char* auth; int entity; } instances[] = {
            {"A", "A", 0}, {"B", "B", 1},
            {"C", "A", 2}, {"D", "A", 3}, {"E", "A", 4}, {"E", "A", 4}, {"F", "A", 5}, {"F", "A", 5}, {"F", "A", 5}, {"F", "A", 5}, {"F", "A", 5},
            {"G", "B", 6}, {"H", "B", 4}, {"I", "B", 5}, {"I", "B", 5}, {"I", "B", 5},
            {"J", "C", 7},
        };
        ASSERT_EQ(md_system_instance_count(&sys), ARRAY_SIZE(instances));
        for (size_t i = 0; i < ARRAY_SIZE(instances); ++i) {
            EXPECT_TRUE(str_eq(md_system_instance_id(&sys, i),      str_from_cstr(instances[i].id)));
            EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, i), str_from_cstr(instances[i].auth)));
            EXPECT_EQ(md_system_instance_entity_idx(&sys, i),       instances[i].entity);
        }
    }

    {
        // 64 lipids, the first 26 with chain ids of their own, and the waters
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(init_system(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/dppc64.pdb")));
        EXPECT_GT(md_system_atom_count(&sys),   0);

        ASSERT_EQ(md_system_entity_count(&sys), 2);
        EXPECT_EQ(md_system_entity_flags(&sys, 0), inferred_kind(MD_ENTITY_KIND_NON_POLYMER));
        EXPECT_EQ(md_system_entity_flags(&sys, 1), inferred_kind(MD_ENTITY_KIND_WATER));
        EXPECT_TRUE(str_eq(md_entity_description(&sys.entity, 0), STR_LIT("DPP")));

        ASSERT_EQ(md_system_instance_count(&sys), md_system_component_count(&sys));
        for (size_t i = 0; i < 64; ++i) {
            EXPECT_EQ(md_system_instance_entity_idx(&sys, i), 0);
        }
        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 25), STR_LIT("Z")));
        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 26), STR_LIT("AA")));
        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 63), STR_LIT("AA")));
        for (size_t i = 64; i < md_system_instance_count(&sys); ++i) {
            EXPECT_EQ(md_system_instance_entity_idx(&sys, i), 1);
            EXPECT_TRUE(str_eq(md_system_instance_id(&sys, i), STR_LIT("AB")));
        }
    }

    md_temp_end(temp_scope);
}

// Entities the file does not define are inferred on load and flagged as such, file defined ones are not
UTEST(util, entity_derived) {
    md_temp_scope_t temp_scope = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp_scope);

    {
        // Nucleosome: two DNA strands and eight histones, a .gro carries no entities
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_gro_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro")));
        md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);

        ASSERT_EQ(md_system_entity_count(&sys), 6);
        // The DNA strands are DNA (by their residue names), which the nucleic backbone extraction relies upon
        EXPECT_EQ(md_system_entity_flags(&sys, 0), inferred_kind(MD_ENTITY_KIND_DNA));
        EXPECT_EQ(md_system_entity_flags(&sys, 1), inferred_kind(MD_ENTITY_KIND_DNA));
        EXPECT_TRUE(str_eq(md_entity_description(&sys.entity, 0), STR_LIT("DNA")));
        for (size_t i = 2; i < 6; ++i) {
            EXPECT_EQ(md_system_entity_flags(&sys, i), inferred_kind(MD_ENTITY_KIND_PEPTIDE));
        }
        EXPECT_EQ(md_system_instance_count(&sys), 10);
        EXPECT_EQ(sys.nucleic_backbone.range.count, 2);
        EXPECT_EQ(sys.protein_backbone.range.count, 8);

        // Deriving again starts over rather than appending
        ASSERT_TRUE(md_util_system_infer_entity_and_instance(&sys, NULL));
        EXPECT_EQ(md_system_entity_count(&sys), 6);
        EXPECT_EQ(md_system_instance_count(&sys), 10);
    }

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif")));
        md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);

        ASSERT_GT(md_system_entity_count(&sys), 0);
        for (size_t i = 0; i < md_system_entity_count(&sys); ++i) {
            EXPECT_FALSE(md_system_entity_flags(&sys, i) & MD_ENTITY_FLAG_INFERRED);
        }
    }

    md_temp_end(temp_scope);
}

// mmCIF: the entity kinds are the file's, and the molecules of a non-polymer asym (the waters of a chain) are
// instances of their own which share the asym's id
UTEST(util, entity_instance_mmcif) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif")));

        static const struct { md_entity_kind_t kind; const char* desc; } entities[] = {
            {MD_ENTITY_KIND_PEPTIDE, "PHOSPHONOACETALDEHYDE HYDROLASE"}, {MD_ENTITY_KIND_NON_POLYMER, "MAGNESIUM ION"},
            {MD_ENTITY_KIND_NON_POLYMER, "TUNGSTATE(VI)ION"}, {MD_ENTITY_KIND_WATER, "water"},
        };
        ASSERT_EQ(md_system_entity_count(&sys), ARRAY_SIZE(entities));
        for (size_t i = 0; i < ARRAY_SIZE(entities); ++i) {
            EXPECT_EQ(md_system_entity_kind(&sys, i), entities[i].kind);
            EXPECT_TRUE(str_eq(md_entity_description(&sys.entity, i), str_from_cstr(entities[i].desc)));
        }

        static const struct { const char* id; const char* auth; int entity; } instances[] = {
            {"A", "A", 0}, {"B", "B", 0}, {"E", "A", 1}, {"F", "B", 2}, {"G", "B", 1},
            {"K", "A", 3}, {"K", "A", 3}, {"K", "A", 3}, {"K", "A", 3},
            {"L", "B", 3}, {"L", "B", 3}, {"L", "B", 3}, {"L", "B", 3},
        };
        ASSERT_EQ(md_system_instance_count(&sys), ARRAY_SIZE(instances));
        for (size_t i = 0; i < ARRAY_SIZE(instances); ++i) {
            EXPECT_TRUE(str_eq(md_system_instance_id(&sys, i),      str_from_cstr(instances[i].id)));
            EXPECT_TRUE(str_eq(md_system_instance_auth_id(&sys, i), str_from_cstr(instances[i].auth)));
            EXPECT_EQ(md_system_instance_entity_idx(&sys, i),       instances[i].entity);
            if (i >= 2) EXPECT_EQ(md_system_instance_comp_count(&sys, i), 1);
        }
        // MG is an ion as a component, whatever its entity is
        EXPECT_EQ(md_system_component_kind(&sys, md_system_instance_comp_range(&sys, 2).beg), MD_COMPONENT_KIND_ION);
    }

    {
        // RNA strands which begin with a GTP (the 5' triphosphate): linked into the chain, so a nucleotide
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_mmcif_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif")));
        md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL);
        ASSERT_EQ(md_system_entity_count(&sys), 5);
        EXPECT_EQ(md_system_entity_kind(&sys, 2), MD_ENTITY_KIND_RNA);
        EXPECT_EQ(md_system_entity_kind(&sys, 3), MD_ENTITY_KIND_RNA);
        size_t resolved = 0;
        for (size_t i = 0; i < md_system_instance_count(&sys); ++i) {
            if (!md_entity_kind_is_nucleic_acid(md_system_instance_entity_kind(&sys, i))) continue;
            const md_urange_t range = md_system_instance_comp_range(&sys, i);
            EXPECT_TRUE(str_eq(md_component_name(&sys.component, range.beg), STR_LIT("GTP")) || str_eq(md_component_name(&sys.component, range.beg), STR_LIT("UTP")));
            for (uint32_t c = range.beg; c < range.end; ++c) {
                EXPECT_EQ(md_system_component_kind(&sys, c), MD_COMPONENT_KIND_NUCLEOTIDE);
                resolved += (md_system_component_flags(&sys, c) & MD_COMPONENT_FLAG_RESOLVED) != 0;
            }
        }
        EXPECT_EQ(resolved, 48);
        EXPECT_EQ(sys.nucleic_backbone.range.count, 2);
    }

    md_vm_arena_destroy(alloc);
}

// Topologies: each molecule type is an entity, each molecule an instance, and the kinds are classified from the residues
UTEST(util, entity_instance_topology) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/tpr/peptide_tip4p.tpr")));

        static const struct { md_entity_kind_t kind; const char* desc; size_t count; } entities[] = {
            {MD_ENTITY_KIND_PEPTIDE, "Protein_chain_U", 1}, {MD_ENTITY_KIND_WATER, "SOL", 506}, {MD_ENTITY_KIND_NON_POLYMER, "NA", 2}, {MD_ENTITY_KIND_NON_POLYMER, "CL", 2},
        };
        ASSERT_EQ(md_system_entity_count(&sys), ARRAY_SIZE(entities));
        size_t inst = 0;
        for (size_t e = 0; e < ARRAY_SIZE(entities); ++e) {
            EXPECT_EQ(md_system_entity_flags(&sys, e), md_entity_flags_set_kind(MD_ENTITY_FLAG_NONE, entities[e].kind));
            EXPECT_TRUE(str_eq(md_entity_description(&sys.entity, e), str_from_cstr(entities[e].desc)));
            for (size_t k = 0; k < entities[e].count; ++k, ++inst) {
                EXPECT_EQ(md_system_instance_entity_idx(&sys, inst), (int)e);
            }
        }
        ASSERT_EQ(md_system_instance_count(&sys), inst);
        EXPECT_EQ(md_system_instance_comp_count(&sys, 0), 2);
        // Single residue molecules share the id of their block
        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 0),   STR_LIT("A")));
        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 1),   STR_LIT("B")));
        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 506), STR_LIT("B")));
        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 507), STR_LIT("C")));
        EXPECT_TRUE(str_eq(md_system_instance_id(&sys, 510), STR_LIT("D")));

        // Inference leaves the topology's own alone
        ASSERT_TRUE(md_util_system_infer(&sys, &sys_state, MD_UTIL_INFER_ALL & ~MD_UTIL_INFER_BOND_BIT));
        EXPECT_EQ(md_system_entity_count(&sys), ARRAY_SIZE(entities));
        EXPECT_EQ(md_system_instance_count(&sys), inst);
    }

    {
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/tpr/martini3.tpr")));
        // The peptide, two lipids, four ions and 583 waters, as the structures are
        ASSERT_EQ(md_system_entity_count(&sys), 5);
        EXPECT_EQ(md_system_entity_kind(&sys, 0), MD_ENTITY_KIND_PEPTIDE);
        EXPECT_EQ(md_system_entity_kind(&sys, 1), MD_ENTITY_KIND_NON_POLYMER);
        EXPECT_EQ(md_system_entity_kind(&sys, 4), MD_ENTITY_KIND_WATER);
        EXPECT_EQ(md_system_instance_count(&sys), 1 + 2 + 4 + 583);
    }

    {
        // LAMMPS knows its molecules but not what they are: an entity per sequence of atom types, described by its
        // formula, and water told by its atoms
        md_system_t sys = {.alloc = alloc};
        md_system_state_t sys_state = { .alloc = alloc };
        const char* atom_format = md_lammps_atom_format_strings()[MD_LAMMPS_ATOM_FORMAT_FULL];
        ASSERT_TRUE(md_lammps_system_init_from_file(&sys, &sys_state, STR_LIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Cubic_Init.data"), atom_format));
        ASSERT_EQ(md_system_entity_count(&sys), 2);
        EXPECT_TRUE(str_eq(md_entity_description(&sys.entity, 0), STR_LIT("C2H6")));
        EXPECT_TRUE(str_eq(md_entity_description(&sys.entity, 1), STR_LIT("H2O")));
        EXPECT_EQ(md_system_entity_kind(&sys, 0), MD_ENTITY_KIND_NON_POLYMER);
        EXPECT_EQ(md_system_entity_kind(&sys, 1), MD_ENTITY_KIND_WATER);
        ASSERT_EQ(md_system_instance_count(&sys), md_system_component_count(&sys));
        EXPECT_EQ(md_system_instance_entity_idx(&sys, 0), 0);
        EXPECT_EQ(md_system_instance_entity_idx(&sys, md_system_instance_count(&sys) - 1), 1);
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
    md_bitfield_t mask = md_bitfield_create(alloc);
    md_bitfield_set_bit(&mask, 0);
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
    md_util_backbone_secondary_structure_infer(ref_ss, n, st.xyz, &st.unitcell, &sys.protein_backbone);
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

// The binary searches from an atom to its component and from a component to its instance
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
    for (size_t c = 0; c < sys.component.count; ++c) {
        const md_urange_t range = md_system_component_atom_range(&sys, c);
        EXPECT_EQ(md_system_component_find_by_atom_idx(&sys, range.beg), (int)c);
        EXPECT_EQ(md_system_component_find_by_atom_idx(&sys, range.end - 1), (int)c);
    }
    for (size_t i = 0; i < sys.instance.count; ++i) {
        const md_urange_t range = md_system_instance_atom_range(&sys, i);
        EXPECT_EQ(md_system_instance_find_by_atom_idx(&sys, range.beg), (int)i);
        EXPECT_EQ(md_system_instance_find_by_atom_idx(&sys, range.end - 1), (int)i);
    }
    EXPECT_EQ(md_system_component_find_by_atom_idx(&sys, sys.atom.count), -1);
    md_vm_arena_destroy(alloc);
}

// The walk keeps a per atom depth in temp memory, where zero means 'not reached'. Temp memory is not cleared on
// allocation, so a walk which followed some other use of temp memory used to see stale depths and stop early.
UTEST_F(util, mask_grow_by_bonds_dirty_temp) {
    const md_system_t* sys = &utest_fixture->mol_ala;
    md_bitfield_t mask = md_bitfield_create(utest_fixture->alloc);
    md_bitfield_t ref  = md_bitfield_create(utest_fixture->alloc);

    // Stale depths of 1 make every neighbour look as if it had already been reached at a lower depth
    const size_t junk_size = MEGABYTES(1);
    {
        md_temp_scope_t temp = md_temp_begin();
        memset(md_temp_alloc(temp, junk_size), 1, junk_size);
        md_temp_end(temp);
    }
    md_bitfield_set_bit(&mask, 0);
    md_util_mask_grow_by_bonds(&mask, sys, 3, NULL);

    {
        md_temp_scope_t temp = md_temp_begin();
        memset(md_temp_alloc(temp, junk_size), 0, junk_size);
        md_temp_end(temp);
    }
    md_bitfield_set_bit(&ref, 0);
    md_util_mask_grow_by_bonds(&ref, sys, 3, NULL);

    EXPECT_GT(md_bitfield_popcount(&ref), (size_t)1);
    EXPECT_EQ(md_bitfield_popcount(&ref), md_bitfield_popcount(&mask));
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
