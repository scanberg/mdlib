#include "system_invariants.h"

#include <md_pdb.h>
#include <md_gro.h>
#include <md_mmcif.h>
#include <md_lammps.h>
#include <md_tpr.h>
#include <md_xyz.h>
#include <md_molden.h>
#include <md_bond.h>
#include <md_hierarchy.h>
#include <md_ring.h>
#include <md_fragment.h>
#include <md_system_fields.h>
#include <md_infer.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <string.h>

typedef bool (*sys_load_fn)(md_system_t*, md_system_state_t*, md_str_t);

// Load a file, run full inference, and assert the structural invariants hold.
static void check_file(int* utest_result, sys_load_fn load, md_str_t path) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* alloc = md_temp_allocator(temp);

    md_system_t sys = { .alloc = alloc };
    md_system_state_t state = { .alloc = alloc };

    if (!load(&sys, &state, path)) {
        printf("failed to load '" MD_STR_FMT "'\n", MD_STR_ARG(path));
        *utest_result = UTEST_TEST_FAILURE;
        md_temp_end(temp);
        return;
    }

    md_system_infer(&sys, &state, NULL);
    expect_valid_system(utest_result, __FILE__, __LINE__, &sys);

    md_system_free(&sys);
    md_temp_end(temp);
}

static bool load_pdb(md_system_t* sys, md_system_state_t* state, md_str_t path) {
    return md_pdb_system_init_from_file(sys, state, path, NULL);
}

static bool load_mmcif(md_system_t* sys, md_system_state_t* state, md_str_t path) {
    return md_mmcif_system_init_from_file(sys, state, path, NULL);
}

static bool load_gro(md_system_t* sys, md_system_state_t* state, md_str_t path) {
    return md_gro_system_init_from_file(sys, state, path, NULL);
}

static bool load_lammps(md_system_t* sys, md_system_state_t* state, md_str_t path) {
    return md_lammps_system_init_from_file(sys, state, path, NULL);
}

UTEST(system_invariants, pdb) {
    check_file(utest_result, load_pdb, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"));
    check_file(utest_result, load_pdb, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"));
    check_file(utest_result, load_pdb, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb"));
}

UTEST(system_invariants, mmcif) {
    check_file(utest_result, load_mmcif, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif"));
    check_file(utest_result, load_mmcif, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/2or2.cif"));
    check_file(utest_result, load_mmcif, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif"));
}

UTEST(system_invariants, gro) {
    check_file(utest_result, load_gro, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.gro"));
    check_file(utest_result, load_gro, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro"));
    check_file(utest_result, load_gro, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/nucleotides.gro"));
}

UTEST(system_invariants, lammps) {
    check_file(utest_result, load_lammps, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Cubic_Init.data"));
}

// A copy is the same system in another allocator: valid, and inferred on it gives what inferred on the original gave
UTEST(system_invariants, copy) {
    md_allocator_i* a1 = md_vm_arena_create(MD_GIGABYTES(1));
    md_allocator_i* a2 = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t sys = { .alloc = a1 };
    md_system_state_t state = { .alloc = a1 };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), NULL));

    // The core, copied and inferred on, as viamd does in the background
    md_system_t core = { .alloc = a2 };
    ASSERT_TRUE(md_system_copy(&core, &sys));
    EXPECT_EQ(sys.topology_version, core.topology_version);
    EXPECT_EQ(sys.atom.count, core.atom.count);
    EXPECT_NE((void*)sys.atom.type_idx, (void*)core.atom.type_idx);
    md_system_state_t core_state = { .alloc = a2 };
    md_system_state_copy(&core_state, &state);
    ASSERT_TRUE(md_system_infer(&core, &core_state, NULL));
    ASSERT_TRUE(md_system_infer(&sys, &state, NULL));
    expect_valid_system(utest_result, __FILE__, __LINE__, &core);
    EXPECT_EQ(sys.bond.count, core.bond.count);
    EXPECT_EQ(0, memcmp(sys.bond.pairs, core.bond.pairs, sys.bond.count * sizeof(md_atom_pair_t)));
    EXPECT_EQ(0, memcmp(sys.residue.flags, core.residue.flags, sys.residue.count * sizeof(md_residue_flags_t)));
    EXPECT_EQ(sys.part.count, core.part.count);
    EXPECT_EQ(md_system_molecule_count(&sys), md_system_molecule_count(&core));

    // And the inferred system copied: everything there, nothing shared
    md_allocator_i* a3 = md_vm_arena_create(MD_GIGABYTES(1));
    md_system_t full = { .alloc = a3 };
    ASSERT_TRUE(md_system_copy(&full, &sys));
    expect_valid_system(utest_result, __FILE__, __LINE__, &full);
    EXPECT_EQ(sys.protein_backbone.segment.count, full.protein_backbone.segment.count);
    EXPECT_EQ(md_array_size(sys.ring.rings), md_array_size(full.ring.rings));
    EXPECT_EQ(sys.reference.num_atoms, full.reference.num_atoms);
    EXPECT_NE((void*)sys.reference.xyz, (void*)full.reference.xyz);
    for (size_t t = 0; t < sys.part.type.count; ++t) {
        EXPECT_TRUE(md_str_eq(md_system_part_type_description(&sys, t), md_system_part_type_description(&full, t)));
    }

    md_vm_arena_destroy(a1);
    md_vm_arena_destroy(a2);
    md_vm_arena_destroy(a3);
}

// ### INFERENCE WHILE THE SYSTEM IS IN USE ###
// md_system_infer_prepare shares the CORE of a system with the work and copies the DERIVED (md_system_fields.h).
// That is sound only if inference writes none of the core: these check it of systems of every kind, and that what is
// published is what inferring the system itself gives.

#include <md_system_fields.h>

static bool same_array(const void* a, size_t a_bytes, const void* b, size_t b_bytes) {
    if (a_bytes != b_bytes) return false;
    return a_bytes == 0 || memcmp(a, b, a_bytes) == 0;
}

static bool same_strings(const md_str_t* a, const md_str_t* b) {
    if (md_array_size(a) != md_array_size(b)) return false;
    for (size_t i = 0; i < md_array_size(a); ++i) {
        if (!md_str_eq(a[i], b[i])) return false;
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

static void check_infer_write_set(int* utest_result, sys_load_fn load, md_str_t path) {
    md_allocator_i* a_sys  = md_vm_arena_create(MD_GIGABYTES(4));
    md_allocator_i* a_work = md_vm_arena_create(MD_GIGABYTES(4));
    md_allocator_i* a_snap = md_vm_arena_create(MD_GIGABYTES(4));
    md_allocator_i* a_ref  = md_vm_arena_create(MD_GIGABYTES(4));

    md_system_t sys = { .alloc = a_sys };
    md_system_state_t state = { .alloc = a_sys };
    if (!load(&sys, &state, path)) {
        printf("failed to load '" MD_STR_FMT "'\n", MD_STR_ARG(path));
        *utest_result = UTEST_TEST_FAILURE;
    } else {
        const char* name = (const char*)path.ptr;
        // The colors of the types first, as a viewer has them to show the system with: of the types, which the work
        // does not publish (md_system_infer_run)
        md_system_infer_type_colors(&sys);
        md_system_t snap = { .alloc = a_snap };
        md_system_copy(&snap, &sys);

        md_system_infer_work_t work;
        EXPECT_TRUE(md_system_infer_prepare(&work, &sys, &state, a_work));
        EXPECT_TRUE(md_system_infer_run(&work, NULL));

        // Nothing of the system was written: not through what the work shares of it
        expect_same_core(utest_result, __FILE__, __LINE__, name, &sys, &snap);
        expect_same_derived(utest_result, __FILE__, __LINE__, name, &sys, &snap);
        // And inference wrote none of the core, of what the work has of its own (the types of the atoms) neither
        expect_same_core(utest_result, __FILE__, __LINE__, name, &work.sys, &sys);

        // What the system inferred itself is: loaded anew, with its attributes (the formal charges of the file)
        md_system_t ref = { .alloc = a_ref };
        md_system_state_t ref_state = { .alloc = a_ref };
        ASSERT_TRUE(load(&ref, &ref_state, path));
        md_system_infer_type_colors(&ref);
        md_system_infer(&ref, &ref_state, NULL);

        EXPECT_TRUE(md_system_infer_publish(&sys, &work));
        expect_valid_system(utest_result, __FILE__, __LINE__, &sys);
        expect_same_core(utest_result, __FILE__, __LINE__, name, &sys, &ref);
        expect_same_derived(utest_result, __FILE__, __LINE__, name, &sys, &ref);
        EXPECT_EQ(sys.reference.num_atoms, state.num_atoms);
        if (state.xyz) {
            EXPECT_EQ(0, memcmp(sys.reference.xyz, state.xyz, state.num_atoms * sizeof(md_float3_t)));
        }

        // Of a system whose topology changed meanwhile, the work is not published
        const uint64_t version = sys.topology_version;
        EXPECT_TRUE(md_system_infer_prepare(&work, &sys, &state, a_work));
        md_system_topology_changed(&sys);
        EXPECT_FALSE(md_system_infer_publish(&sys, &work));
        EXPECT_NE(version, sys.topology_version);
    }

    md_vm_arena_destroy(a_sys);
    md_vm_arena_destroy(a_work);
    md_vm_arena_destroy(a_snap);
    md_vm_arena_destroy(a_ref);
}

UTEST(system_invariants, infer_write_set) {
    check_infer_write_set(utest_result, load_pdb, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"));
    check_infer_write_set(utest_result, load_pdb, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb"));
    check_infer_write_set(utest_result, load_pdb, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/dppc64.pdb"));
    check_infer_write_set(utest_result, load_pdb, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/ciprofloxacin.pdb"));
    check_infer_write_set(utest_result, load_mmcif, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif"));
    check_infer_write_set(utest_result, load_mmcif, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/2or2.cif"));
    check_infer_write_set(utest_result, load_gro, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro"));
    check_infer_write_set(utest_result, load_gro, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.gro"));
    check_infer_write_set(utest_result, load_gro, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini2_membrane_protein.gro"));
    check_infer_write_set(utest_result, load_gro, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini3_membrane.gro"));
    check_infer_write_set(utest_result, load_lammps, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Cubic_Init.data"));
}

// ### THE ALLOCATOR OF THE DERIVED ###
//
// What inference writes (MD_SYSTEM_DERIVED_*, and the reference state) is in the system's derived arena, whose pages
// come from where the system's alloc takes its own (md_system_t, ALLOCATORS), and the core is in alloc: an inference
// replaces the derived arena whole. Checked with an arena of the core over a backing which knows its blocks: a loaded
// system has no derived; after inferring, every derived array is within one of its blocks; inferring again, publishing
// a work, adding and removing a bond leave as much of it alive as before; the derived freed, as much as the core alone
// had (inference wrote nothing of the core); freed, the system leaves none of it.

typedef struct tracked_t {
    md_allocator_i alloc;
    uintptr_t* key;         // Open addressing, 0 empty, 1 removed
    size_t*    size;
    size_t     cap;
    size_t     used;        // Keys and removed slots
    size_t     count;
    size_t     live;        // Bytes
} tracked_t;

static size_t tracked_slot(const tracked_t* t, uintptr_t k) {
    uint64_t h = (uint64_t)k * 0x9E3779B97F4A7C15ull;
    return (size_t)(h >> 20) & (t->cap - 1);
}

static void tracked_insert(tracked_t* t, uintptr_t k, size_t size);

static void tracked_grow(tracked_t* t) {
    uintptr_t* key = t->key;
    size_t* size = t->size;
    const size_t cap = t->cap;
    t->cap = cap ? cap * 2 : 1024;
    t->key = (uintptr_t*)calloc(t->cap, sizeof(uintptr_t));
    t->size = (size_t*)calloc(t->cap, sizeof(size_t));
    t->used = t->count = 0;
    for (size_t i = 0; i < cap; ++i) if (key[i] > 1) tracked_insert(t, key[i], size[i]);
    free(key);
    free(size);
}

static void tracked_insert(tracked_t* t, uintptr_t k, size_t size) {
    if ((t->used + 1) * 2 > t->cap) tracked_grow(t);
    size_t i = tracked_slot(t, k);
    while (t->key[i] > 1) i = (i + 1) & (t->cap - 1);
    if (t->key[i] == 0) t->used += 1;
    t->key[i] = k;
    t->size[i] = size;
    t->count += 1;
}

static size_t tracked_find(const tracked_t* t, uintptr_t k) {
    if (!t->cap) return SIZE_MAX;
    size_t i = tracked_slot(t, k);
    while (t->key[i]) {
        if (t->key[i] == k) return i;
        i = (i + 1) & (t->cap - 1);
    }
    return SIZE_MAX;
}

static void* tracked_realloc(md_allocator_o* inst, void* ptr, size_t old_size, size_t new_size, const char* file, size_t line) {
    (void)file; (void)line;
    tracked_t* t = (tracked_t*)inst;
    if (ptr) {
        const size_t i = tracked_find(t, (uintptr_t)ptr);
        MD_ASSERT(i != SIZE_MAX);   // Given back to the allocator it came from
        if (i == SIZE_MAX) return NULL;
        MD_ASSERT(t->size[i] == old_size);
        t->key[i] = 1;
        t->count -= 1;
        t->live -= t->size[i];
    }
    if (new_size == 0) {
        free(ptr);
        return NULL;
    }
    void* out = realloc(ptr, new_size);
    tracked_insert(t, (uintptr_t)out, new_size);
    t->live += new_size;
    return out;
}

static void tracked_init(tracked_t* t) {
    memset(t, 0, sizeof(*t));
    t->alloc.inst = (md_allocator_o*)t;
    t->alloc.realloc = tracked_realloc;
}

static void tracked_free(tracked_t* t) {
    free(t->key);
    free(t->size);
}

// Within one of the blocks (the pages of an arena)
static bool tracked_contains(const tracked_t* t, const void* ptr) {
    const uintptr_t p = (uintptr_t)ptr;
    for (size_t i = 0; i < t->cap; ++i) {
        if (t->key[i] > 1 && p >= t->key[i] && p < t->key[i] + t->size[i]) return true;
    }
    return false;
}

static bool tracked_owns_array(const tracked_t* t, const void* arr) {
    return !arr || tracked_contains(t, md_array_header(arr));
}

// Every derived array of the system a block of t, its strings too, and the reference state's coordinates
static void expect_derived_owned(int* utest_result, const char* what, const md_system_t* sys, const tracked_t* t) {
#define OWNED(f) if (!tracked_owns_array(t, sys->f)) { printf("%s: " #f " is not of the derived allocator\n", what); *utest_result = UTEST_TEST_FAILURE; }
#define OWNED_STRINGS(f) OWNED(f) for (size_t i_ = 0; i_ < md_array_size(sys->f); ++i_) { \
        if (sys->f[i_].ptr && !tracked_contains(t, sys->f[i_].ptr)) { printf("%s: a string of " #f " is not of the derived allocator\n", what); *utest_result = UTEST_TEST_FAILURE; break; } }
    MD_SYSTEM_DERIVED_ARRAYS(OWNED)
    MD_SYSTEM_DERIVED_STRINGS(OWNED_STRINGS)
    OWNED(reference.xyz)
#undef OWNED
#undef OWNED_STRINGS
}

static bool load_tpr(md_system_t* sys, md_system_state_t* state, md_str_t path) {
    return md_tpr_system_init_from_file(sys, state, path, NULL);
}

static bool load_xyz(md_system_t* sys, md_system_state_t* state, md_str_t path) {
    return md_xyz_system_init_from_file(sys, state, path, NULL);
}

static bool load_molden(md_system_t* sys, md_system_state_t* state, md_str_t path) {
    return md_molden_system_init_from_file(sys, state, path, NULL);
}

static void check_derived_arena(int* utest_result, sys_load_fn load, md_str_t path) {
    const char* what = (const char*)path.ptr;
    tracked_t t;
    tracked_init(&t);
    md_allocator_i* core = md_arena_allocator_create(&t.alloc, MD_MEGABYTES(1));
    md_allocator_i* coords = md_vm_arena_create(MD_GIGABYTES(4));
    md_system_t sys = { .alloc = core };
    md_system_state_t state = { .alloc = coords };
    if (!load(&sys, &state, path)) {
        printf("failed to load '%s'\n", what);
        *utest_result = UTEST_TEST_FAILURE;
        md_arena_allocator_destroy(core);
        md_vm_arena_destroy(coords);
        tracked_free(&t);
        return;
    }
    // A loader writes nothing of the derived
    EXPECT_TRUE(sys.derived_arena == NULL);
    const size_t base = t.live;

    md_system_infer(&sys, &state, NULL);
    expect_valid_system(utest_result, __FILE__, __LINE__, &sys);
    expect_derived_owned(utest_result, what, &sys, &t);
    ASSERT_TRUE(sys.derived_arena != NULL);
    EXPECT_TRUE(md_arena_allocator_backing(sys.derived_arena) == &t.alloc);    // Where the core's pages come from
    const size_t live = t.live;
    EXPECT_GT(live, base);

    // Inferred again: the same, in an arena of its own, and the one before given back
    md_system_infer(&sys, &state, NULL);
    expect_derived_owned(utest_result, what, &sys, &t);
    EXPECT_EQ(live, t.live);

    // In the background: the work's arena becomes the system's, the rest of the work's is gone with its own arena
    for (int k = 0; k < 2; ++k) {
        md_allocator_i* a_work = md_vm_arena_create(MD_GIGABYTES(4));
        md_system_infer_work_t work;
        EXPECT_TRUE(md_system_infer_prepare(&work, &sys, &state, a_work));
        EXPECT_TRUE(md_arena_allocator_backing(work.sys.derived_arena) == &t.alloc);    // The system's, not the work's
        EXPECT_TRUE(md_system_infer_run(&work, NULL));
        EXPECT_TRUE(md_system_infer_publish(&sys, &work));
        md_vm_arena_destroy(a_work);
        expect_derived_owned(utest_result, what, &sys, &t);
        expect_valid_system(utest_result, __FILE__, __LINE__, &sys);
        EXPECT_EQ(live, t.live);
    }
    // A work not published is given back whole
    {
        md_allocator_i* a_work = md_vm_arena_create(MD_GIGABYTES(4));
        md_system_infer_work_t work;
        EXPECT_TRUE(md_system_infer_prepare(&work, &sys, &state, a_work));
        EXPECT_TRUE(md_system_infer_run(&work, NULL));
        md_system_infer_discard(&work);
        md_vm_arena_destroy(a_work);
        EXPECT_EQ(live, t.live);
    }

    // The derived freed: as much alive as the core alone had, which inference wrote nothing of
    md_system_derived_free(&sys);
    EXPECT_EQ(base, t.live);
    md_system_infer(&sys, &state, NULL);
    EXPECT_EQ(live, t.live);

    // A bond added by hand between the first and the last atom, and removed: the molecules, rings and fragments with it.
    // The bond is of the core (bond_edit), which may grow by a page; the derived is inferred anew.
    if (sys.atom.count > 1) {
        const md_atom_idx_t a = 0, b = (md_atom_idx_t)(sys.atom.count - 1);
        if (md_system_bond_find(&sys, a, b) == -1) {
            EXPECT_TRUE(md_system_bond_insert(&sys, a, b, md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_USER)));
            md_system_infer(&sys, &state, NULL);
            expect_derived_owned(utest_result, what, &sys, &t);
            EXPECT_NE(-1, md_system_bond_find(&sys, a, b));
            EXPECT_TRUE(md_system_bond_remove(&sys, a, b, md_system_bond_flags(&sys, md_system_bond_find(&sys, a, b))));
            md_system_infer(&sys, &state, NULL);
            expect_derived_owned(utest_result, what, &sys, &t);
            EXPECT_EQ(-1, md_system_bond_find(&sys, a, b));
            EXPECT_LE(t.live, live + MD_MEGABYTES(2));
        }
    }

    // A copy has its derived in an arena of its own, of its own core's pages
    {
        tracked_t t2;
        tracked_init(&t2);
        md_allocator_i* core2 = md_arena_allocator_create(&t2.alloc, MD_MEGABYTES(1));
        md_system_t copy = { .alloc = core2 };
        EXPECT_TRUE(md_system_copy(&copy, &sys));
        EXPECT_TRUE(copy.derived_arena != NULL && copy.derived_arena != sys.derived_arena);
        EXPECT_TRUE(md_arena_allocator_backing(copy.derived_arena) == &t2.alloc);
        expect_derived_owned(utest_result, what, &copy, &t2);
        md_system_free(&copy);
        md_arena_allocator_destroy(core2);
        EXPECT_EQ((size_t)0, t2.live);
        tracked_free(&t2);
    }

    // Freed, and its arena: none of it left
    md_system_free(&sys);
    md_arena_allocator_destroy(core);
    EXPECT_EQ((size_t)0, t.live);
    if (t.live) printf("%s: %zu bytes left in %zu blocks\n", what, t.live, t.count);
    md_vm_arena_destroy(coords);
    tracked_free(&t);
}

UTEST(system_invariants, derived_arena) {
    check_derived_arena(utest_result, load_pdb, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"));
    check_derived_arena(utest_result, load_pdb, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/dppc64.pdb"));
    check_derived_arena(utest_result, load_pdb, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/ciprofloxacin.pdb"));
    check_derived_arena(utest_result, load_mmcif, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif"));
    check_derived_arena(utest_result, load_mmcif, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/2or2.cif"));
    check_derived_arena(utest_result, load_gro, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro"));
    check_derived_arena(utest_result, load_gro, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini3_membrane.gro"));
    check_derived_arena(utest_result, load_lammps, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Cubic_Init.data"));
    check_derived_arena(utest_result, load_tpr, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tpr/peptide_tip3p.tpr"));
    check_derived_arena(utest_result, load_tpr, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini3_membrane.tpr"));
    check_derived_arena(utest_result, load_xyz, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/extended.xyz"));
    check_derived_arena(utest_result, load_molden, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden"));
}

// ### THE DERIVED FROM THE CORE ALONE ###
//
// Everything inference writes (MD_SYSTEM_DERIVED_*) follows from the core: what the data gives (md_system_fields.h) and
// the coordinates. So the derived freed (md_system_derived_free) and inferred again is what it was: nothing the data
// gives is kept only in the derived, where inference had to carry it over from what it replaces. That is what lets the
// derived be replaced whole (an arena of its own, reset) rather than array by array.

static void check_derived_from_core(int* utest_result, sys_load_fn load, md_str_t path) {
    const char* what = (const char*)path.ptr;
    md_allocator_i* a_sys = md_vm_arena_create(MD_GIGABYTES(4));
    md_allocator_i* a_ref = md_vm_arena_create(MD_GIGABYTES(4));
    md_system_t sys = { .alloc = a_sys };
    md_system_state_t state = { .alloc = a_sys };
    md_system_t ref = { .alloc = a_ref };
    md_system_state_t ref_state = { .alloc = a_ref };
    if (!load(&sys, &state, path) || !load(&ref, &ref_state, path)) {
        printf("failed to load '%s'\n", what);
        *utest_result = UTEST_TEST_FAILURE;
    } else {
        md_system_infer(&ref, &ref_state, NULL);
        md_system_derived_free(&sys);
        md_system_infer(&sys, &state, NULL);
        expect_same_core(utest_result, __FILE__, __LINE__, what, &sys, &ref);
        expect_same_derived(utest_result, __FILE__, __LINE__, what, &sys, &ref);
        // And once more, from what was inferred
        md_system_derived_free(&sys);
        md_system_infer(&sys, &state, NULL);
        expect_same_derived(utest_result, __FILE__, __LINE__, what, &sys, &ref);
    }
    md_vm_arena_destroy(a_sys);
    md_vm_arena_destroy(a_ref);
}

UTEST(system_invariants, derived_from_core) {
    check_derived_from_core(utest_result, load_pdb, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"));
    check_derived_from_core(utest_result, load_pdb, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/ciprofloxacin.pdb"));     // CONECT
    check_derived_from_core(utest_result, load_pdb, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/c60.pdb"));
    check_derived_from_core(utest_result, load_mmcif, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/8g7u.cif"));            // Entities
    check_derived_from_core(utest_result, load_mmcif, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1fez.cif"));
    check_derived_from_core(utest_result, load_gro, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro"));
    check_derived_from_core(utest_result, load_gro, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini3_membrane.gro"));
    check_derived_from_core(utest_result, load_lammps, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/Water_Ethane_Cubic_Init.data"));
    check_derived_from_core(utest_result, load_tpr, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tpr/peptide_tip3p.tpr"));     // Topology
    check_derived_from_core(utest_result, load_tpr, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/martini/martini3_membrane.tpr"));
    check_derived_from_core(utest_result, load_xyz, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/extended.xyz"));
    check_derived_from_core(utest_result, load_molden, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden"));  // QM
}
