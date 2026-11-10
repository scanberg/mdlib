#include "utest.h"

#include <md_mcs.h>
#include <md_system.h>
#include <md_util.h>
#include <md_pdb.h>
#include <md_gro.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

// Maximum common substructure on tubulin (its ligands and residues, and its two chains) and on DNA

typedef struct mcs_fixture_t {
    md_system_t         sys;
    md_system_state_t   state;
    md_allocator_i*     alloc;
    bool                ok;
} mcs_fixture_t;

static mcs_fixture_t* mcs_tubulin(void) {
    static mcs_fixture_t f;
    static bool loaded;
    if (!loaded) {
        loaded = true;
        f.alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
        f.sys = (md_system_t){.alloc = f.alloc};
        f.state = (md_system_state_t){.alloc = f.alloc};
        f.ok = md_pdb_system_init_from_file(&f.sys, &f.state, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tubulin-A-B.pdb"), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE, NULL) &&
               md_util_system_infer(&f.sys, &f.state, MD_UTIL_INFER_ALL, NULL);
    }
    return f.ok ? &f : NULL;
}

// The atoms of the nth residue of that name, into atoms; how many
static size_t mcs_residue(int32_t* atoms, const md_system_t* sys, const char* name, int nth) {
    for (size_t r = 0; r < sys->residue.count; ++r) {
        if (!md_str_eq(md_system_residue_name(sys, r), md_str_from_cstr(name)) || nth-- > 0) continue;
        size_t n = 0;
        for (uint32_t a = sys->residue.atom_offset[r]; a < sys->residue.atom_offset[r + 1]; ++a) atoms[n++] = (int32_t)a;
        return n;
    }
    return 0;
}

static md_mcs_result_t mcs_of(const mcs_fixture_t* f, const char* na, int ia, const char* nb, int ib, bool xyz) {
    int32_t a[256], b[256];
    const size_t ca = mcs_residue(a, &f->sys, na, ia), cb = mcs_residue(b, &f->sys, nb, ib);
    md_mcs_desc_t d = {.a = {.sys = &f->sys, .atoms = a, .count = ca, .xyz = xyz ? f->state.xyz : NULL},
                       .b = {.sys = &f->sys, .atoms = b, .count = cb, .xyz = xyz ? f->state.xyz : NULL}};
    md_mcs_result_t r;
    if (!md_mcs(&r, &d, md_get_heap_allocator())) memset(&r, 0, sizeof(r));
    return r;
}

UTEST(mcs, ligands) {
    const mcs_fixture_t* f = mcs_tubulin();
    ASSERT_TRUE(f != NULL);
    // GDP is GTP less a phosphate: all of it is common
    md_mcs_result_t r = mcs_of(f, "GTP", 0, "GDP", 0, false);
    EXPECT_EQ(28u, (unsigned)r.count);
    EXPECT_EQ(32u, (unsigned)r.a_heavy);
    EXPECT_EQ(28u, (unsigned)r.b_heavy);
    EXPECT_FALSE(r.by_residue);
    EXPECT_FALSE(r.incomplete);
    EXPECT_TRUE(isnan(r.rmsd));
    // Pairs alike by element, a's atoms ascending
    for (size_t i = 0; i < r.count; ++i) {
        EXPECT_EQ(md_atom_atomic_number(&f->sys.atom, r.a[i]), md_atom_atomic_number(&f->sys.atom, r.b[i]));
        if (i) EXPECT_LT(r.a[i - 1], r.a[i]);
    }
    md_mcs_result_free(&r);
    // With positions, of the mappings the phosphates allow, the one of the least RMSD; the first found is not
    md_mcs_result_t g = mcs_of(f, "GTP", 0, "GDP", 0, true);
    ASSERT_EQ(28u, (unsigned)g.count);
    printf("mcs: GTP onto GDP %.2f Å\n", g.rmsd);
    EXPECT_LT(g.rmsd, 0.5);
    md_mcs_result_free(&g);
    // Itself: every heavy atom onto itself
    md_mcs_result_t v = mcs_of(f, "VLB", 0, "VLB", 0, true);
    EXPECT_EQ(59u, (unsigned)v.count);
    for (size_t i = 0; i < v.count; ++i) EXPECT_EQ(v.a[i], v.b[i]);
    EXPECT_LT(v.rmsd, 1e-3);
    md_mcs_result_free(&v);
    // Unrelated molecules still share something; the caller judges how much is enough
    md_mcs_result_t u = mcs_of(f, "GTP", 0, "VLB", 0, false);
    EXPECT_GT(u.count, 3u);
    EXPECT_LT((double)u.count / (double)u.a_heavy, 0.5);
    md_mcs_result_free(&u);
}

UTEST(mcs, residues) {
    const mcs_fixture_t* f = mcs_tubulin();
    ASSERT_TRUE(f != NULL);
    struct { const char* a; const char* b; unsigned common; } cases[] = {
        {"PHE", "TYR", 11}, {"PHE", "TRP", 10}, {"LYS", "ARG", 7}, {"ASP", "GLU", 6}, {"LEU", "ILE", 7}, {"SER", "CYS", 5},
    };
    for (size_t i = 0; i < MD_COUNTOF(cases); ++i) {
        md_mcs_result_t r = mcs_of(f, cases[i].a, 0, cases[i].b, 0, true);
        EXPECT_EQ(cases[i].common, (unsigned)r.count);
        if (r.count != cases[i].common) printf("mcs: %s and %s share %zu\n", cases[i].a, cases[i].b, r.count);
        md_mcs_result_free(&r);
    }
}

UTEST(mcs, chains) {
    const mcs_fixture_t* f = mcs_tubulin();
    ASSERT_TRUE(f != NULL);
    const md_system_t* sys = &f->sys;
    md_allocator_i* heap = md_get_heap_allocator();
    int32_t* a = md_alloc(heap, sys->atom.count * sizeof(int32_t));
    int32_t* b = md_alloc(heap, sys->atom.count * sizeof(int32_t));
    size_t na = 0, nb = 0;
    for (size_t r = 0; r < sys->residue.count; ++r) {
        if (md_system_residue_kind(sys, r) != MD_RESIDUE_KIND_AMINO_ACID) continue;
        const int ch = sys->chain.res_chain[r];
        for (uint32_t x = sys->residue.atom_offset[r]; x < sys->residue.atom_offset[r + 1]; ++x) {
            if (ch == 0) a[na++] = (int32_t)x;
            else if (ch == 1) b[nb++] = (int32_t)x;
        }
    }
    // α against β tubulin: residue by residue, at once
    md_mcs_desc_t d = {.a = {.sys = sys, .atoms = a, .count = na, .xyz = f->state.xyz}, .b = {.sys = sys, .atoms = b, .count = nb, .xyz = f->state.xyz}};
    const md_tick_t t0 = md_tick_now();
    md_mcs_result_t r;
    ASSERT_TRUE(md_mcs(&r, &d, heap));
    const double ms = md_tick_to_milliseconds(md_tick_now() - t0);
    printf("mcs: chain A against B, %zu of %zu and %zu heavy atoms, %.2f Å, %.2f ms\n", r.count, r.a_heavy, r.b_heavy, r.rmsd, ms);
    EXPECT_TRUE(r.by_residue);
    EXPECT_GT(r.count, 2500u);
    EXPECT_LT(r.rmsd, 4.0);
    md_mcs_result_free(&r);
    // As one graph it is too large
    d.flags = MD_MCS_FLAG_GRAPH;
    EXPECT_FALSE(md_mcs(&r, &d, heap));
    // A chain against itself: every heavy atom onto itself
    md_mcs_desc_t s = {.a = {.sys = sys, .atoms = a, .count = na}, .b = {.sys = sys, .atoms = a, .count = na}};
    ASSERT_TRUE(md_mcs(&r, &s, heap));
    EXPECT_EQ(r.a_heavy, r.count);
    bool same = true;
    for (size_t i = 0; i < r.count; ++i) same = same && r.a[i] == r.b[i];
    EXPECT_TRUE(same);
    md_mcs_result_free(&r);
    md_free(heap, a, sys->atom.count * sizeof(int32_t));
    md_free(heap, b, sys->atom.count * sizeof(int32_t));
}

UTEST(mcs, nucleotides) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t st = {.alloc = alloc};
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro"), NULL));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL, NULL));
    // A purine and a pyrimidine share the sugar and what the base's ring holds in common with it: the first DA has no
    // phosphate, 18 heavy atoms
    int32_t a[64], b[64];
    size_t na = 0, nb = 0;
    for (size_t r = 0; r < sys.residue.count && (!na || !nb); ++r) {
        const md_str_t n = md_system_residue_name(&sys, r);
        int32_t* dst = NULL;
        size_t* cnt = NULL;
        if (!na && md_str_begins_with(n, MD_STR_LIT("DA"))) { dst = a; cnt = &na; }
        else if (!nb && md_str_begins_with(n, MD_STR_LIT("DT"))) { dst = b; cnt = &nb; }
        if (!dst) continue;
        for (uint32_t x = sys.residue.atom_offset[r]; x < sys.residue.atom_offset[r + 1] && *cnt < 64; ++x) dst[(*cnt)++] = (int32_t)x;
    }
    ASSERT_GT(na, 0u);
    ASSERT_GT(nb, 0u);
    md_mcs_desc_t d = {.a = {.sys = &sys, .atoms = a, .count = na, .xyz = st.xyz}, .b = {.sys = &sys, .atoms = b, .count = nb, .xyz = st.xyz}};
    md_mcs_result_t r;
    ASSERT_TRUE(md_mcs(&r, &d, alloc));
    printf("mcs: DA and DT share %zu of %zu and %zu heavy atoms, %.2f Å\n", r.count, r.a_heavy, r.b_heavy, r.rmsd);
    EXPECT_FALSE(r.by_residue);
    EXPECT_GE(r.count, 12u);
    EXPECT_LT(r.rmsd, 1.5);
    md_arena_allocator_destroy(alloc);
}
