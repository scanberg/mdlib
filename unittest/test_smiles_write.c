#include "utest.h"

#include <md_smiles.h>
#include <md_match.h>
#include <md_system.h>
#include <md_infer.h>
#include <md_pdb.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_bitset.h>

#include <string.h>
#include <stdlib.h>

// md_smiles_write: canonical whatever the order of the atoms, read back by md_smiles_parse, and found again by the
// matcher, with atom classes for the atoms a pattern is about.

static const md_system_t* tubulin(void) {
    static md_system_t sys;
    static bool loaded, ok;
    if (!loaded) {
        loaded = true;
        md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
        sys = (md_system_t){.alloc = alloc};
        md_system_state_t st = {.alloc = alloc};
        ok = md_pdb_system_init_from_file(&sys, &st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb"), NULL) &&
             md_system_infer(&sys, &st, MD_SYSTEM_INFER_ALL, NULL);
    }
    return ok ? &sys : NULL;
}

// The heavy atoms of the first residue of a name
static size_t residue_atoms(int32_t* out, size_t cap, const md_system_t* sys, const char* name) {
    for (size_t r = 0; r < sys->residue.count; ++r) {
        if (!md_str_eq_cstr(md_system_residue_name(sys, r), name)) continue;
        const md_urange_t ra = md_system_residue_atom_range(sys, r);
        size_t n = 0;
        for (uint32_t i = ra.beg; i < ra.end && n < cap; ++i) if (md_atom_atomic_number(&sys->atom, i) != 1) out[n++] = (int32_t)i;
        return n;
    }
    return 0;
}

static uint32_t next_rand(uint32_t* s) { *s = *s * 1664525u + 1013904223u; return *s >> 8; }

UTEST(smiles_write, canonical_and_found_again) {
    const md_system_t* sys = tubulin();
    ASSERT_TRUE(sys != NULL);
    const char* names[] = {"GTP", "GDP", "VLB", "PHE", "TRP", "HIS", "ARG", "PRO"};
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    uint32_t seed = 7;
    for (size_t k = 0; k < MD_COUNTOF(names); ++k) {
        int32_t atoms[256];
        const size_t n = residue_atoms(atoms, 256, sys, names[k]);
        ASSERT_GT(n, 2u);
        char ref[1024];
        const size_t len = md_smiles_write(ref, sizeof(ref), sys, atoms, n, NULL);
        ASSERT_GT(len, 0u);
        ASSERT_LT(len, sizeof(ref));
        // The same text in any order
        for (int t = 0; t < 20; ++t) {
            for (size_t i = n - 1; i > 0; --i) {
                const size_t j = next_rand(&seed) % (i + 1);
                const int32_t tmp = atoms[i]; atoms[i] = atoms[j]; atoms[j] = tmp;
            }
            char got[1024];
            md_smiles_write(got, sizeof(got), sys, atoms, n, NULL);
            if (strcmp(ref, got)) printf("%s: %s\n   vs %s\n", names[k], ref, got);
            EXPECT_STREQ(ref, got);
        }
        // Read back, with every atom; and the residue is among what it matches
        md_smiles_t parsed = {0};
        md_smiles_error_t err = {0};
        EXPECT_TRUE(md_smiles_parse(&parsed, md_str_from_cstr(ref), alloc, &err));
        EXPECT_EQ(n, parsed.num_atoms);
        md_match_query_t q = {0};
        ASSERT_TRUE(md_match_query_init_smiles(&q, md_str_from_cstr(ref), alloc, NULL));
        md_match_desc_t desc = {.query = &q, .level = MD_MATCH_LEVEL_RESIDUE, .mode = MD_MATCH_MODE_ALL};
        md_match_result_t res = {0};
        ASSERT_TRUE(md_match_find(&res, &desc, sys, alloc));
        bool found = false;
        for (size_t m = 0; m < res.count && !found; ++m) {
            bool all = true;
            for (size_t i = 0; i < n && all; ++i) {
                bool in = false;
                for (size_t c = 0; c < res.width && !in; ++c) in = res.atom_idx[m * res.width + c] == atoms[i];
                all = in;
            }
            found = all;
        }
        if (!found) printf("%s: %s is not found again (%zu matches)\n", names[k], ref, res.count);
        EXPECT_TRUE(found);
    }
    md_arena_allocator_destroy(alloc);
}

UTEST(smiles_write, classes_and_limits) {
    const md_system_t* sys = tubulin();
    ASSERT_TRUE(sys != NULL);
    int32_t atoms[64];
    const size_t n = residue_atoms(atoms, 64, sys, "PHE");
    ASSERT_GT(n, 3u);
    // The first three atoms tagged 1, 2, 3, whatever their order in the text
    uint32_t cls[64] = {0};
    cls[0] = 1; cls[1] = 2; cls[2] = 3;
    char buf[512];
    ASSERT_GT(md_smiles_write(buf, sizeof(buf), sys, atoms, n, cls), 0u);
    EXPECT_TRUE(strstr(buf, ":1]") && strstr(buf, ":2]") && strstr(buf, ":3]"));
    // The length as snprintf gives it, past a small buffer
    char small[4];
    const size_t full = md_smiles_write(small, sizeof(small), sys, atoms, n, NULL);
    EXPECT_GT(full, sizeof(small));
    EXPECT_EQ(strlen(small), sizeof(small) - 1);
    // Nothing to write
    EXPECT_EQ(0u, md_smiles_write(buf, sizeof(buf), sys, atoms, 0, NULL));
    const int32_t twice[2] = {atoms[0], atoms[0]};
    EXPECT_EQ(0u, md_smiles_write(buf, sizeof(buf), sys, twice, 2, NULL));
    const int32_t out[1] = {(int32_t)sys->atom.count};
    EXPECT_EQ(0u, md_smiles_write(buf, sizeof(buf), sys, out, 1, NULL));
    // Two residues not bonded: two components
    int32_t two[128];
    size_t m = residue_atoms(two, 64, sys, "GTP");
    m += residue_atoms(two + m, 64, sys, "SO4");
    ASSERT_GT(md_smiles_write(buf, sizeof(buf), sys, two, m, NULL), 0u);
    EXPECT_TRUE(strchr(buf, '.') != NULL);
}
