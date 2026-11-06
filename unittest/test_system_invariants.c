#include "system_invariants.h"

#include <md_pdb.h>
#include <md_gro.h>
#include <md_mmcif.h>
#include <md_lammps.h>
#include <md_util.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <string.h>

typedef bool (*sys_load_fn)(md_system_t*, md_system_state_t*, str_t);

// Load a file, run full inference, and assert the structural invariants hold.
static void check_file(int* utest_result, sys_load_fn load, str_t path) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };

    if (!load(&sys, &state, path)) {
        printf("failed to load '" STR_FMT "'\n", STR_ARG(path));
        *utest_result = UTEST_TEST_FAILURE;
        md_temp_end(temp);
        return;
    }

    md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL, NULL);
    expect_valid_system(utest_result, __FILE__, __LINE__, &sys);

    md_system_free(&sys);
    md_temp_end(temp);
}

static bool load_pdb(md_system_t* sys, md_system_state_t* state, str_t path) {
    return md_pdb_system_init_from_file(sys, state, path, MD_PDB_OPTION_NONE, NULL);
}

static bool load_gro(md_system_t* sys, md_system_state_t* state, str_t path) {
    return md_gro_system_init_from_file(sys, state, path, NULL);
}

static bool load_lammps(md_system_t* sys, md_system_state_t* state, str_t path) {
    return md_lammps_system_init_from_file(sys, state, path, NULL);
}

UTEST(system_invariants, pdb) {
    check_file(utest_result, load_pdb, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"));
    check_file(utest_result, load_pdb, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"));
    check_file(utest_result, load_pdb, STR_LIT(MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb"));
}

UTEST(system_invariants, mmcif) {
    check_file(utest_result, md_mmcif_system_init_from_file, STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif"));
    check_file(utest_result, md_mmcif_system_init_from_file, STR_LIT(MD_UNITTEST_DATA_DIR "/2or2.cif"));
    check_file(utest_result, md_mmcif_system_init_from_file, STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif"));
}

UTEST(system_invariants, gro) {
    check_file(utest_result, load_gro, STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.gro"));
    check_file(utest_result, load_gro, STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro"));
    check_file(utest_result, load_gro, STR_LIT(MD_UNITTEST_DATA_DIR "/nucleotides.gro"));
}

UTEST(system_invariants, lammps) {
    check_file(utest_result, load_lammps, STR_LIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Cubic_Init.data"));
}

// A copy is the same system in another allocator: valid, and inferred on it gives what inferred on the original gave
UTEST(system_invariants, copy) {
    md_allocator_i* a1 = md_vm_arena_create(GIGABYTES(1));
    md_allocator_i* a2 = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = a1 };
    md_system_state_t state = { .alloc = a1 };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), MD_PDB_OPTION_NONE, NULL));

    // The core, copied and inferred on, as viamd does in the background
    md_system_t core = { .alloc = a2 };
    ASSERT_TRUE(md_system_copy(&core, &sys));
    EXPECT_EQ(sys.topology_version, core.topology_version);
    EXPECT_EQ(sys.atom.count, core.atom.count);
    EXPECT_NE((void*)sys.atom.type_idx, (void*)core.atom.type_idx);
    md_system_state_t core_state = { .alloc = a2 };
    md_system_state_copy(&core_state, &state);
    ASSERT_TRUE(md_util_system_infer(&core, &core_state, MD_UTIL_INFER_ALL, NULL));
    ASSERT_TRUE(md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL, NULL));
    expect_valid_system(utest_result, __FILE__, __LINE__, &core);
    EXPECT_EQ(sys.bond.count, core.bond.count);
    EXPECT_EQ(0, memcmp(sys.bond.pairs, core.bond.pairs, sys.bond.count * sizeof(md_atom_pair_t)));
    EXPECT_EQ(0, memcmp(sys.residue.flags, core.residue.flags, sys.residue.count * sizeof(md_residue_flags_t)));
    EXPECT_EQ(sys.part.count, core.part.count);
    EXPECT_EQ(md_system_molecule_count(&sys), md_system_molecule_count(&core));

    // And the inferred system copied: everything there, nothing shared
    md_allocator_i* a3 = md_vm_arena_create(GIGABYTES(1));
    md_system_t full = { .alloc = a3 };
    ASSERT_TRUE(md_system_copy(&full, &sys));
    expect_valid_system(utest_result, __FILE__, __LINE__, &full);
    EXPECT_EQ(sys.protein_backbone.segment.count, full.protein_backbone.segment.count);
    EXPECT_EQ(md_array_size(sys.ring.rings), md_array_size(full.ring.rings));
    EXPECT_EQ(sys.reference.num_atoms, full.reference.num_atoms);
    EXPECT_NE((void*)sys.reference.xyz, (void*)full.reference.xyz);
    for (size_t t = 0; t < sys.part.type.count; ++t) {
        EXPECT_TRUE(str_eq(md_system_part_type_description(&sys, t), md_system_part_type_description(&full, t)));
    }

    md_vm_arena_destroy(a1);
    md_vm_arena_destroy(a2);
    md_vm_arena_destroy(a3);
}

// ### INFERENCE WHILE THE SYSTEM IS IN USE ###
// md_util_system_infer_prepare shares the CORE of a system with the work and copies the DERIVED (md_system_fields.h).
// That is sound only if inference writes none of the core: these check it of systems of every kind, and that what is
// published is what inferring the system itself gives.

#include <md_system_fields.h>

static bool same_array(const void* a, size_t a_bytes, const void* b, size_t b_bytes) {
    if (a_bytes != b_bytes) return false;
    return a_bytes == 0 || memcmp(a, b, a_bytes) == 0;
}

static bool same_strings(const str_t* a, const str_t* b) {
    if (md_array_size(a) != md_array_size(b)) return false;
    for (size_t i = 0; i < md_array_size(a); ++i) {
        if (!str_eq(a[i], b[i])) return false;
    }
    return true;
}

#define SAME_FIELD(f) \
    if (!same_array(a->f, md_array_size(a->f) * sizeof(*a->f), b->f, md_array_size(b->f) * sizeof(*b->f))) { \
        printf("%s:%i: %s: '%s' differs (%zu and %zu)\n", file, line, what, #f, md_array_size(a->f), md_array_size(b->f)); \
        *utest_result = UTEST_TEST_FAILURE; \
    }
#define SAME_STRINGS(f) \
    if (!same_strings(a->f, b->f)) { \
        printf("%s:%i: %s: '%s' differs\n", file, line, what, #f); \
        *utest_result = UTEST_TEST_FAILURE; \
    }

static void expect_same_core(int* utest_result, const char* file, int line, const char* what, const md_system_t* a, const md_system_t* b) {
    MD_SYSTEM_CORE_ARRAYS(SAME_FIELD)
}

static void expect_same_derived(int* utest_result, const char* file, int line, const char* what, const md_system_t* a, const md_system_t* b) {
    MD_SYSTEM_DERIVED_ARRAYS(SAME_FIELD)
    MD_SYSTEM_DERIVED_STRINGS(SAME_STRINGS)
#define SAME_COUNT(f) \
    if (a->f.count != b->f.count) { \
        printf("%s:%i: %s: '%s.count' differs (%zu and %zu)\n", file, line, what, #f, a->f.count, b->f.count); \
        *utest_result = UTEST_TEST_FAILURE; \
    }
    SAME_COUNT(chain) SAME_COUNT(part) SAME_COUNT(part.type) SAME_COUNT(residue.type) SAME_COUNT(lipid) SAME_COUNT(bond)
    SAME_COUNT(molecule) SAME_COUNT(molecule.type) SAME_COUNT(protein_backbone.range) SAME_COUNT(nucleic_backbone.range)
#undef SAME_COUNT
    if (a->bond.noncovalent_count != b->bond.noncovalent_count) {
        printf("%s:%i: %s: 'bond.noncovalent_count' differs\n", file, line, what);
        *utest_result = UTEST_TEST_FAILURE;
    }
}

#undef SAME_FIELD
#undef SAME_STRINGS

static void check_infer_write_set(int* utest_result, sys_load_fn load, str_t path) {
    md_allocator_i* a_sys  = md_vm_arena_create(GIGABYTES(4));
    md_allocator_i* a_work = md_vm_arena_create(GIGABYTES(4));
    md_allocator_i* a_snap = md_vm_arena_create(GIGABYTES(4));
    md_allocator_i* a_ref  = md_vm_arena_create(GIGABYTES(4));

    md_system_t sys = { .alloc = a_sys };
    md_system_state_t state = { .alloc = a_sys };
    if (!load(&sys, &state, path)) {
        printf("failed to load '" STR_FMT "'\n", STR_ARG(path));
        *utest_result = UTEST_TEST_FAILURE;
    } else {
        const char* name = (const char*)path.ptr;
        // The colors of the types first, as a viewer has them to show the system with: of the types, which the work
        // does not publish (md_util_system_infer_run)
        md_util_system_infer(&sys, &state, MD_UTIL_INFER_COLOR_BIT, NULL);
        md_system_t snap = { .alloc = a_snap };
        md_system_copy(&snap, &sys);

        md_util_infer_work_t work;
        EXPECT_TRUE(md_util_system_infer_prepare(&work, &sys, &state, a_work));
        EXPECT_TRUE(md_util_system_infer_run(&work, MD_UTIL_INFER_ALL, NULL));

        // Nothing of the system was written: not through what the work shares of it
        expect_same_core(utest_result, __FILE__, __LINE__, name, &sys, &snap);
        expect_same_derived(utest_result, __FILE__, __LINE__, name, &sys, &snap);
        // And inference wrote none of the core, of what the work has of its own (the types of the atoms) neither
        expect_same_core(utest_result, __FILE__, __LINE__, name, &work.sys, &sys);

        // What the system inferred itself is: loaded anew, with its attributes (the formal charges of the file)
        md_system_t ref = { .alloc = a_ref };
        md_system_state_t ref_state = { .alloc = a_ref };
        ASSERT_TRUE(load(&ref, &ref_state, path));
        md_util_system_infer(&ref, &ref_state, MD_UTIL_INFER_COLOR_BIT, NULL);
        md_util_system_infer(&ref, &ref_state, MD_UTIL_INFER_ALL & ~MD_UTIL_INFER_COLOR_BIT, NULL);

        EXPECT_TRUE(md_util_system_infer_publish(&sys, &work));
        expect_valid_system(utest_result, __FILE__, __LINE__, &sys);
        expect_same_core(utest_result, __FILE__, __LINE__, name, &sys, &ref);
        expect_same_derived(utest_result, __FILE__, __LINE__, name, &sys, &ref);
        EXPECT_EQ(sys.reference.num_atoms, state.num_atoms);
        if (state.xyz) {
            EXPECT_EQ(0, memcmp(sys.reference.xyz, state.xyz, state.num_atoms * sizeof(vec3_t)));
        }

        // Of a system whose topology changed meanwhile, the work is not published
        const uint64_t version = sys.topology_version;
        EXPECT_TRUE(md_util_system_infer_prepare(&work, &sys, &state, a_work));
        md_system_topology_changed(&sys);
        EXPECT_FALSE(md_util_system_infer_publish(&sys, &work));
        EXPECT_NE(version, sys.topology_version);
    }

    md_vm_arena_destroy(a_sys);
    md_vm_arena_destroy(a_work);
    md_vm_arena_destroy(a_snap);
    md_vm_arena_destroy(a_ref);
}

UTEST(system_invariants, infer_write_set) {
    check_infer_write_set(utest_result, load_pdb, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"));
    check_infer_write_set(utest_result, load_pdb, STR_LIT(MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb"));
    check_infer_write_set(utest_result, load_pdb, STR_LIT(MD_UNITTEST_DATA_DIR "/dppc64.pdb"));
    check_infer_write_set(utest_result, load_pdb, STR_LIT(MD_UNITTEST_DATA_DIR "/ciprofloxacin.pdb"));
    check_infer_write_set(utest_result, md_mmcif_system_init_from_file, STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif"));
    check_infer_write_set(utest_result, md_mmcif_system_init_from_file, STR_LIT(MD_UNITTEST_DATA_DIR "/2or2.cif"));
    check_infer_write_set(utest_result, load_gro, STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro"));
    check_infer_write_set(utest_result, load_gro, STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.gro"));
    check_infer_write_set(utest_result, load_gro, STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini2_membrane_protein.gro"));
    check_infer_write_set(utest_result, load_gro, STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini3_membrane.gro"));
    check_infer_write_set(utest_result, load_lammps, STR_LIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Cubic_Init.data"));
}
