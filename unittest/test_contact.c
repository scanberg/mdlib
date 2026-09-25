#include "utest.h"

#include <md_contact.h>
#include <md_system.h>
#include <md_pdb.h>
#include <md_gro.h>
#include <md_util.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_bitfield.h>
#include <core/md_vec_math.h>

#include <float.h>
#include <math.h>
#include <string.h>

#include "contact_reference.h"
#include "run_check.h"

static bool check_against_reference(const md_contact_desc_t* desc, const md_system_t* sys, const md_system_state_t* state, md_allocator_i* alloc) {
    md_contact_set_t set;
    if (!md_contact_compute(&set, desc, sys, state, alloc)) {
        printf("md_contact_compute failed\n");
        return false;
    }
    md_array(ref_pair_t) ref = ref_contacts(desc, sys, state, alloc);
    const bool ok = same_contacts(&set, ref, md_array_size(ref));
    if (ok && set.count == 0) {
        printf("Reference found no contacts, the test is vacuous\n");
        return false;
    }
    md_contact_set_free(&set);
    return ok;
}

// One group per component whose flags intersect the mask (or every component if the mask is 0), at most max_count
static md_bitfield_t* make_residue_groups(size_t* out_count, const md_system_t* sys, md_flags_t mask, size_t max_count, md_allocator_i* alloc) {
    md_bitfield_t* groups = md_alloc(alloc, sizeof(md_bitfield_t) * sys->component.count);
    size_t n = 0;
    for (size_t c = 0; c < sys->component.count && n < max_count; ++c) {
        if (mask && !(md_system_component_flags(sys, c) & mask)) continue;
        const md_urange_t range = md_system_component_atom_range(sys, c);
        groups[n] = md_bitfield_create(alloc);
        md_bitfield_set_range(&groups[n], range.beg, range.end);
        ++n;
    }
    *out_count = n;
    return groups;
}

struct contact {
    md_allocator_i* alloc;
    md_system_t ala;
    md_system_state_t ala_state;
    md_system_t npt;
    md_system_state_t npt_state;
};

UTEST_F_SETUP(contact) {
    utest_fixture->alloc = md_vm_arena_create(GIGABYTES(4));

    utest_fixture->ala.alloc = utest_fixture->alloc;
    utest_fixture->ala_state = (md_system_state_t){ .alloc = utest_fixture->alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&utest_fixture->ala, &utest_fixture->ala_state, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    md_util_system_infer(&utest_fixture->ala, &utest_fixture->ala_state, MD_UTIL_INFER_ALL);

    // A triclinic cell, the case which exposes mistakes in the periodic neighbour walk
    utest_fixture->npt.alloc = utest_fixture->alloc;
    utest_fixture->npt_state = (md_system_state_t){ .alloc = utest_fixture->alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&utest_fixture->npt, &utest_fixture->npt_state, STR_LIT(MD_UNITTEST_DATA_DIR "/npt.gro")));
    md_util_system_infer(&utest_fixture->npt, &utest_fixture->npt_state, MD_UTIL_INFER_ALL);
}

UTEST_F_TEARDOWN(contact) {
    md_vm_arena_destroy(utest_fixture->alloc);
}

UTEST_F(contact, self_distance) {
    size_t num = 0;
    md_bitfield_t* res = make_residue_groups(&num, &utest_fixture->ala, 0, SIZE_MAX, utest_fixture->alloc);
    ASSERT_GT(num, (size_t)10);
    md_contact_desc_t desc = { .group_a = res, .num_a = num, .cutoff = 3.5 };
    EXPECT_TRUE(check_against_reference(&desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));
}

UTEST_F(contact, between_sets) {
    size_t num = 0;
    md_bitfield_t* res = make_residue_groups(&num, &utest_fixture->ala, 0, SIZE_MAX, utest_fixture->alloc);
    ASSERT_GT(num, (size_t)10);
    // Overlapping sets: the first two thirds of the residues against the last two thirds
    const size_t third = num / 3;
    md_contact_desc_t desc = { .group_a = res, .num_a = num - third, .group_b = res + third, .num_b = num - third, .cutoff = 4.0 };
    EXPECT_TRUE(check_against_reference(&desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));

    // Protein against everything
    size_t num_prot = 0;
    md_bitfield_t* prot = make_residue_groups(&num_prot, &utest_fixture->ala, MD_FLAG_AMINO_ACID, SIZE_MAX, utest_fixture->alloc);
    ASSERT_GT(num_prot, (size_t)0);
    desc = (md_contact_desc_t){ .group_a = prot, .num_a = num_prot, .group_b = res, .num_b = num, .cutoff = 3.0 };
    EXPECT_TRUE(check_against_reference(&desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));
}

UTEST_F(contact, bond_exclusion_and_separation) {
    size_t num = 0;
    md_bitfield_t* prot = make_residue_groups(&num, &utest_fixture->ala, MD_FLAG_AMINO_ACID, SIZE_MAX, utest_fixture->alloc);
    ASSERT_GT(num, (size_t)3);

    md_contact_desc_t desc = { .group_a = prot, .num_a = num, .cutoff = 4.5, .exclude_bonds = 3 };
    EXPECT_TRUE(check_against_reference(&desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));

    // Neighbouring residues are always in contact, the separation leaves the ones further along the chain
    desc.cutoff = 6.0;
    desc.min_separation = 2;
    EXPECT_TRUE(check_against_reference(&desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));
    {
        md_contact_set_t sep = {0};
        ASSERT_TRUE(md_contact_compute(&sep, &desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));
        for (size_t k = 0; k < sep.count; ++k) {
            EXPECT_GE(sep.j[k] - sep.i[k], 2u);
        }
    }

    // Exclusion removes atom pairs across the peptide bonds, and never adds any
    md_contact_set_t with = {0}, without = {0};
    desc.cutoff = 4.5;
    desc.min_separation = 0;
    ASSERT_TRUE(md_contact_compute(&with, &desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));
    desc.exclude_bonds = 0;
    ASSERT_TRUE(md_contact_compute(&without, &desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));
    size_t sum_with = 0, sum_without = 0;
    for (size_t k = 0; k < with.count; ++k) sum_with += with.atom_pairs[k];
    for (size_t k = 0; k < without.count; ++k) sum_without += without.atom_pairs[k];
    EXPECT_LT(sum_with, sum_without);
    EXPECT_LE(with.count, without.count);
}

UTEST_F(contact, radii) {
    size_t num = 0;
    md_bitfield_t* res = make_residue_groups(&num, &utest_fixture->ala, 0, SIZE_MAX, utest_fixture->alloc);
    md_contact_desc_t desc = { .group_a = res, .num_a = num, .criterion = MD_CONTACT_CRITERION_RADII, .cutoff = 0.3, .exclude_bonds = 3 };
    EXPECT_TRUE(check_against_reference(&desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));

    // A negative tolerance asks for overlap, which bonded atoms have
    desc.cutoff = -0.5;
    desc.exclude_bonds = 0;
    EXPECT_TRUE(check_against_reference(&desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));
}

UTEST_F(contact, triclinic) {
    size_t num = 0;
    md_bitfield_t* all = make_residue_groups(&num, &utest_fixture->npt, 0, SIZE_MAX, utest_fixture->alloc);
    ASSERT_GT(num, (size_t)10);
    // Every third residue: still spread over the whole cell, at a ninth of the cost of the reference
    md_bitfield_t* res = md_alloc(utest_fixture->alloc, sizeof(md_bitfield_t) * num);
    size_t n = 0;
    for (size_t k = 0; k < num; k += 3) res[n++] = all[k];
    num = n;
    md_contact_desc_t desc = { .group_a = res, .num_a = num, .cutoff = 5.0 };
    EXPECT_TRUE(check_against_reference(&desc, &utest_fixture->npt, &utest_fixture->npt_state, utest_fixture->alloc));

    desc = (md_contact_desc_t){ .group_a = res, .num_a = num / 3, .group_b = res + num / 3, .num_b = num - num / 3, .cutoff = 6.0 };
    EXPECT_TRUE(check_against_reference(&desc, &utest_fixture->npt, &utest_fixture->npt_state, utest_fixture->alloc));
}

// Groups sharing atoms, and atoms belonging to no group
UTEST_F(contact, overlapping_groups) {
    size_t num = 0;
    md_bitfield_t* res = make_residue_groups(&num, &utest_fixture->ala, MD_FLAG_AMINO_ACID, SIZE_MAX, utest_fixture->alloc);
    ASSERT_GT(num, (size_t)6);

    // Windows of three consecutive residues, each overlapping the next by two
    const size_t num_win = num - 2;
    md_bitfield_t* win = md_alloc(utest_fixture->alloc, sizeof(md_bitfield_t) * num_win);
    for (size_t w = 0; w < num_win; ++w) {
        win[w] = md_bitfield_create(utest_fixture->alloc);
        md_bitfield_or(&win[w], &res[w], &res[w + 1]);
        md_bitfield_or_inplace(&win[w], &res[w + 2]);
    }
    md_contact_desc_t desc = { .group_a = win, .num_a = num_win, .cutoff = 4.0, .exclude_bonds = 3 };
    EXPECT_TRUE(check_against_reference(&desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));

    desc = (md_contact_desc_t){ .group_a = win, .num_a = num_win, .group_b = res, .num_b = num, .cutoff = 4.0 };
    EXPECT_TRUE(check_against_reference(&desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));
}

// Prepared once, evaluated per frame: the same as preparing for every frame
UTEST_F(contact, query_over_trajectory) {
    md_system_t* sys = &utest_fixture->ala;
    const str_t run = STR_INIT("run/ala");
    ASSERT_TRUE(md_pdb_system_publish_run(sys, STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb"), run, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    const size_t num_frames = run_num_frames(sys, run);
    ASSERT_GT(num_frames, (size_t)1);

    size_t num = 0;
    md_bitfield_t* prot = make_residue_groups(&num, sys, MD_FLAG_AMINO_ACID, SIZE_MAX, utest_fixture->alloc);
    md_contact_desc_t desc = { .group_a = prot, .num_a = num, .cutoff = 4.5, .exclude_bonds = 3, .min_separation = 3 };

    md_contact_query_t query;
    ASSERT_TRUE(md_contact_query_init(&query, &desc, sys, utest_fixture->alloc));

    md_system_state_t state = { .alloc = utest_fixture->alloc };
    ASSERT_TRUE(md_system_state_init(&state, sys->atom.count));

    const str_t paths[] = { STR_INIT("atom/position"), STR_INIT("unitcell") };
    md_system_extract_t* ex = md_system_extract_begin(sys, run, paths, 2, md_get_heap_allocator());
    ASSERT_TRUE(ex != NULL);

    size_t frames_differing = 0;
    size_t prev_count = SIZE_MAX;
    for (size_t f = 0; f < num_frames; ++f) {
        EXPECT_TRUE(md_system_extract_frame(ex, (int64_t)f, &state));
        md_contact_set_t set;
        ASSERT_TRUE(md_contact_query_eval(&set, &query, &state, utest_fixture->alloc));
        EXPECT_EQ(set.flags, (uint32_t)MD_CONTACT_FLAG_SELF);
        EXPECT_EQ(set.num_a, (uint32_t)num);
        for (size_t k = 0; k < set.count; ++k) {
            EXPECT_LT(set.i[k], set.j[k]);
            EXPECT_GE(set.j[k] - set.i[k], 3u);
        }
        if (prev_count != SIZE_MAX && prev_count != set.count) frames_differing += 1;
        prev_count = set.count;

        if (f < 3) {
            md_array(ref_pair_t) ref = ref_contacts(&desc, sys, &state, utest_fixture->alloc);
            EXPECT_TRUE(same_contacts(&set, ref, md_array_size(ref)));
        }
        md_contact_set_free(&set);
    }
    md_system_extract_end(ex);

    // The contacts do change over the trajectory
    EXPECT_GT(frames_differing, (size_t)0);

    // A state of a different system is refused
    md_contact_set_t set;
    EXPECT_FALSE(md_contact_query_eval(&set, &query, &utest_fixture->npt_state, utest_fixture->alloc));

    md_contact_query_free(&query);
}

UTEST_F(contact, empty) {
    md_contact_set_t set;
    // No groups
    md_contact_desc_t desc = { .cutoff = 4.0 };
    EXPECT_TRUE(md_contact_compute(&set, &desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));
    EXPECT_EQ(set.count, (size_t)0);

    // Groups too far apart
    md_bitfield_t g[2] = { md_bitfield_create(utest_fixture->alloc), md_bitfield_create(utest_fixture->alloc) };
    md_bitfield_set_bit(&g[0], 0);
    md_bitfield_set_bit(&g[1], 1);
    desc = (md_contact_desc_t){ .group_a = g, .num_a = 2, .cutoff = 0.01 };
    EXPECT_TRUE(md_contact_compute(&set, &desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));
    EXPECT_EQ(set.count, (size_t)0);

    // An empty B is a set without groups, so no contacts: not the contacts within A, which is what no B means
    static const md_bitfield_t no_groups = {0};
    md_bitfield_t halves[2] = { md_bitfield_create(utest_fixture->alloc), md_bitfield_create(utest_fixture->alloc) };
    md_bitfield_set_range(&halves[0], 0, (uint32_t)utest_fixture->ala.atom.count / 2);
    md_bitfield_set_range(&halves[1], (uint32_t)utest_fixture->ala.atom.count / 2, (uint32_t)utest_fixture->ala.atom.count);
    desc = (md_contact_desc_t){ .group_a = halves, .num_a = 2, .group_b = &no_groups, .num_b = 0, .cutoff = 4.0 };
    EXPECT_TRUE(md_contact_compute(&set, &desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));
    EXPECT_EQ(set.count, (size_t)0);
    EXPECT_EQ(set.num_b, 0u);
    EXPECT_EQ(set.flags, (uint32_t)MD_CONTACT_FLAG_NONE);
    desc.group_b = NULL;
    EXPECT_TRUE(md_contact_compute(&set, &desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));
    EXPECT_GT(set.count, (size_t)0);
    // A count of B groups without the groups is refused
    desc.num_b = 1;
    EXPECT_FALSE(md_contact_compute(&set, &desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));

    // A group referring to atoms beyond the system is refused
    md_bitfield_t bad = md_bitfield_create(utest_fixture->alloc);
    md_bitfield_set_bit(&bad, utest_fixture->ala.atom.count + 5);
    desc = (md_contact_desc_t){ .group_a = &bad, .num_a = 1, .cutoff = 4.0 };
    EXPECT_FALSE(md_contact_compute(&set, &desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));
}


// ### PARTICLE PAIRS ###

typedef struct pair_collect_t {
    md_array(uint64_t) key;     // a << 32 | b
    md_array(float) r;
    md_allocator_i* alloc;
    bool self;
    bool ordered;               // Within one set every pair must come as a < b
} pair_collect_t;

static void collect_pairs(const uint32_t* a, const uint32_t* b, const float* r, size_t count, void* user_param) {
    pair_collect_t* c = (pair_collect_t*)user_param;
    for (size_t k = 0; k < count; ++k) {
        if (c->self && !(a[k] < b[k])) c->ordered = false;
        md_array_push(c->key, ((uint64_t)a[k] << 32) | b[k], c->alloc);
        md_array_push(c->r, r[k], c->alloc);
    }
}

typedef struct key_r_t { uint64_t key; float r; } key_r_t;

static int cmp_key_r(const void* x, const void* y) {
    const uint64_t a = ((const key_r_t*)x)->key, b = ((const key_r_t*)y)->key;
    return (a > b) - (a < b);
}

// Every pair of the stream against the definition: particles of the sets closer than the radius, not excluded
static bool check_pairs(const md_contact_pairs_desc_t* desc, const md_system_t* sys, const md_system_state_t* state, md_allocator_i* alloc) {
    const bool self = desc->set_b == NULL;
    md_contact_pairs_t pairs;
    if (!md_contact_pairs_init(&pairs, desc, sys, alloc)) {
        printf("md_contact_pairs_init failed\n");
        return false;
    }

    // Collected into a temp arena of the caller's, which the stream must leave alone
    md_temp_scope_t temp = md_temp_begin();
    pair_collect_t c = { .alloc = md_temp_allocator(temp), .self = self, .ordered = true };
    bool ok = md_contact_pairs_for_each(&pairs, state, collect_pairs, &c);
    // Anything the stream released beneath the collected pairs gets overwritten here
    memset(md_temp_alloc(temp, MEGABYTES(8)), 0xFF, MEGABYTES(8));
    const size_t n = md_array_size(c.key);

    key_r_t* got = md_alloc(alloc, sizeof(key_r_t) * MAX(n, 1));
    for (size_t k = 0; k < n; ++k) got[k] = (key_r_t){ c.key[k], c.r[k] };
    md_temp_end(temp);
    qsort(got, n, sizeof(key_r_t), cmp_key_r);

    if (!c.ordered) { printf("Pairs within one set not reported as a < b\n"); ok = false; }
    for (size_t k = 1; k < n; ++k) {
        if (got[k].key == got[k-1].key) { printf("Pair reported twice\n"); ok = false; break; }
    }

    md_bitfield_t* excl = md_alloc(alloc, sizeof(md_bitfield_t) * sys->atom.count);
    memset(excl, 0, sizeof(md_bitfield_t) * sys->atom.count);
    const md_bitfield_t* set_b = self ? desc->set_a : desc->set_b;
    const float radius = (float)desc->radius;

    size_t required = 0, found_required = 0, shown = 0;
    md_bitfield_iter_t ia = md_bitfield_iter_create(desc->set_a);
    while (md_bitfield_iter_next(&ia)) {
        const uint32_t a = (uint32_t)md_bitfield_iter_idx(&ia);
        md_bitfield_iter_t jb = md_bitfield_iter_create(set_b);
        while (md_bitfield_iter_next(&jb)) {
            const uint32_t b = (uint32_t)md_bitfield_iter_idx(&jb);
            if (a == b || (self && b < a)) continue;
            if (desc->particle_label && desc->particle_label[a] == desc->particle_label[b]) continue;
            const float d = ref_distance(state, a, b);
            if (!(d < radius + REF_EPS)) continue;
            if (ref_excluded(excl, sys, desc->exclude_bonds, a, b, alloc)) continue;
            const key_r_t probe = { ((uint64_t)a << 32) | b, 0 };
            const key_r_t* hit = bsearch(&probe, got, n, sizeof(key_r_t), cmp_key_r);
            if (d < radius - REF_EPS) {
                required += 1;
                if (hit) found_required += 1;
                else if (shown++ < 10) printf("  missing pair (%u, %u) at %f\n", a, b, d);
            }
            if (hit && fabsf(hit->r - d) > 1.0e-3f) {
                if (shown++ < 10) printf("  pair (%u, %u): distance %f, expected %f\n", a, b, hit->r, d);
                ok = false;
            }
        }
    }
    // Everything reported must be a pair of the definition: check the other way round
    for (size_t k = 0; k < n; ++k) {
        const uint32_t a = (uint32_t)(got[k].key >> 32), b = (uint32_t)(got[k].key & 0xFFFFFFFFu);
        const bool in_sets = md_bitfield_test_bit(desc->set_a, a) && md_bitfield_test_bit(set_b, b) && a != b;
        const float d = ref_distance(state, a, b);
        const bool same_label = desc->particle_label && desc->particle_label[a] == desc->particle_label[b];
        if (!in_sets || same_label || !(d < radius + REF_EPS) || ref_excluded(excl, sys, desc->exclude_bonds, a, b, alloc)) {
            if (shown++ < 10) printf("  unexpected pair (%u, %u) at %f\n", a, b, d);
            ok = false;
        }
    }
    if (found_required != required) ok = false;
    if (required == 0) { printf("No pairs within the radius, the test is vacuous\n"); ok = false; }

    md_contact_pairs_free(&pairs);
    return ok;
}

UTEST_F(contact, pairs) {
    md_system_t* sys = &utest_fixture->ala;
    md_bitfield_t all = md_bitfield_create(utest_fixture->alloc);
    md_bitfield_set_range(&all, 0, sys->atom.count);

    md_contact_pairs_desc_t desc = { .set_a = &all, .radius = 4.0 };
    EXPECT_TRUE(check_pairs(&desc, sys, &utest_fixture->ala_state, utest_fixture->alloc));
    desc.exclude_bonds = 3;
    EXPECT_TRUE(check_pairs(&desc, sys, &utest_fixture->ala_state, utest_fixture->alloc));

    // Two overlapping sets
    md_bitfield_t a = md_bitfield_create(utest_fixture->alloc);
    md_bitfield_t b = md_bitfield_create(utest_fixture->alloc);
    md_bitfield_set_range(&a, 0, sys->atom.count * 2 / 3);
    md_bitfield_set_range(&b, sys->atom.count / 3, sys->atom.count);
    desc = (md_contact_pairs_desc_t){ .set_a = &a, .set_b = &b, .radius = 5.0, .exclude_bonds = 2 };
    EXPECT_TRUE(check_pairs(&desc, sys, &utest_fixture->ala_state, utest_fixture->alloc));
}

UTEST_F(contact, pairs_triclinic) {
    md_system_t* sys = &utest_fixture->npt;
    // Every third particle: spread over the whole cell
    md_bitfield_t sub = md_bitfield_create(utest_fixture->alloc);
    for (size_t k = 0; k < sys->atom.count; k += 3) md_bitfield_set_bit(&sub, k);
    md_contact_pairs_desc_t desc = { .set_a = &sub, .radius = 6.0 };
    EXPECT_TRUE(check_pairs(&desc, sys, &utest_fixture->npt_state, utest_fixture->alloc));
}

// With a parent per group, only groups of the same parent are neighbours
UTEST_F(contact, min_separation_parent) {
    size_t num = 0;
    md_bitfield_t* res = make_residue_groups(&num, &utest_fixture->ala, MD_FLAG_AMINO_ACID, SIZE_MAX, utest_fixture->alloc);
    ASSERT_GT(num, (size_t)8);

    // Two 'chains': residues [0, 7) and [7, num)
    uint32_t* parent = md_alloc(utest_fixture->alloc, sizeof(uint32_t) * num);
    for (size_t k = 0; k < num; ++k) parent[k] = k < 7 ? 0 : 1;

    md_contact_desc_t desc = { .group_a = res, .num_a = num, .cutoff = 4.5, .exclude_bonds = 3, .min_separation = 3, .group_parent = parent };
    EXPECT_TRUE(check_against_reference(&desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));

    md_contact_set_t with_parent = {0}, without_parent = {0};
    ASSERT_TRUE(md_contact_compute(&with_parent, &desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));
    desc.group_parent = NULL;
    ASSERT_TRUE(md_contact_compute(&without_parent, &desc, &utest_fixture->ala, &utest_fixture->ala_state, utest_fixture->alloc));

    // Residues 6 and 7 are consecutive, so in contact, but belong to different parents
    bool found_with = false, found_without = false;
    for (size_t k = 0; k < with_parent.count; ++k)    found_with    |= (with_parent.i[k] == 6 && with_parent.j[k] == 7);
    for (size_t k = 0; k < without_parent.count; ++k) found_without |= (without_parent.i[k] == 6 && without_parent.j[k] == 7);
    EXPECT_TRUE(found_with);
    EXPECT_FALSE(found_without);
}

// ### SAME LABEL EXCLUSION, TYPE PAIR CRITERION ###

// Particles labelled by their residue: the stream only holds pairs between residues
UTEST_F(contact, pairs_label) {
    md_system_t* sys = &utest_fixture->ala;
    md_bitfield_t all = md_bitfield_create(utest_fixture->alloc);
    md_bitfield_set_range(&all, 0, sys->atom.count);

    uint32_t* label = md_alloc(utest_fixture->alloc, sizeof(uint32_t) * sys->atom.count);
    for (size_t c = 0; c < sys->component.count; ++c) {
        const md_urange_t range = md_system_component_atom_range(sys, c);
        for (uint32_t k = range.beg; k < range.end; ++k) label[k] = (uint32_t)c;
    }

    md_contact_pairs_desc_t desc = { .set_a = &all, .radius = 4.5, .particle_label = label };
    EXPECT_TRUE(check_pairs(&desc, sys, &utest_fixture->ala_state, utest_fixture->alloc));

    // Group contacts between those residues are the same with and without the labels: the labels only drop pairs the
    // groups drop anyway, earlier
    size_t num = 0;
    md_bitfield_t* res = make_residue_groups(&num, sys, 0, SIZE_MAX, utest_fixture->alloc);
    md_contact_desc_t gdesc = { .group_a = res, .num_a = num, .cutoff = 4.5, .exclude_bonds = 3 };
    md_contact_set_t plain = {0}, labelled = {0};
    ASSERT_TRUE(md_contact_compute(&plain, &gdesc, sys, &utest_fixture->ala_state, utest_fixture->alloc));
    gdesc.particle_label = label;
    ASSERT_TRUE(md_contact_compute(&labelled, &gdesc, sys, &utest_fixture->ala_state, utest_fixture->alloc));
    ASSERT_EQ(plain.count, labelled.count);
    ASSERT_GT(plain.count, (size_t)0);
    for (size_t k = 0; k < plain.count; ++k) {
        EXPECT_EQ(plain.i[k], labelled.i[k]);
        EXPECT_EQ(plain.j[k], labelled.j[k]);
        EXPECT_EQ(plain.atom_pairs[k], labelled.atom_pairs[k]);
    }
}

// Three made up types with a non additive table: the contact distance of a pair is not a sum of per type sizes
UTEST_F(contact, type_pair) {
    md_system_t* sys = &utest_fixture->ala;
    uint32_t* type = md_alloc(utest_fixture->alloc, sizeof(uint32_t) * sys->atom.count);
    for (size_t k = 0; k < sys->atom.count; ++k) type[k] = (uint32_t)(k % 3);
    const float table[9] = {
        3.0f, 5.5f, 2.5f,
        5.5f, 2.0f, 4.0f,
        2.5f, 4.0f, 4.5f,
    };

    size_t num = 0;
    md_bitfield_t* res = make_residue_groups(&num, sys, 0, SIZE_MAX, utest_fixture->alloc);
    md_contact_desc_t desc = {
        .group_a = res, .num_a = num,
        .criterion = MD_CONTACT_CRITERION_TYPE_PAIR, .particle_type = type, .num_types = 3, .type_cutoff = table,
        .exclude_bonds = 3,
    };
    EXPECT_TRUE(check_against_reference(&desc, sys, &utest_fixture->ala_state, utest_fixture->alloc));

    // Between two sets
    desc.num_a = num / 2;
    desc.group_b = res + num / 3;
    desc.num_b = num - num / 3;
    EXPECT_TRUE(check_against_reference(&desc, sys, &utest_fixture->ala_state, utest_fixture->alloc));

    // A table which is not symmetric, and a type beyond the table, are refused
    md_contact_query_t q;
    float skew[9];
    memcpy(skew, table, sizeof(skew));
    skew[1] = 6.0f;
    desc = (md_contact_desc_t){ .group_a = res, .num_a = num, .criterion = MD_CONTACT_CRITERION_TYPE_PAIR, .particle_type = type, .num_types = 3, .type_cutoff = skew };
    EXPECT_FALSE(md_contact_query_init(&q, &desc, sys, utest_fixture->alloc));
    type[5] = 3;
    desc.type_cutoff = table;
    EXPECT_FALSE(md_contact_query_init(&q, &desc, sys, utest_fixture->alloc));
    type[5] = 2;
    desc.particle_type = NULL;
    EXPECT_FALSE(md_contact_query_init(&q, &desc, sys, utest_fixture->alloc));
}

static int cmp_u64_keys(const void* x, const void* y) {
    const uint64_t a = *(const uint64_t*)x, b = *(const uint64_t*)y;
    return (a > b) - (a < b);
}

// ### CONTACT REGIONS ###

// A contact set made up by hand from unit pairs, as md_contact_query_eval would have produced it
static md_contact_set_t make_contact_set(const uint32_t (*pairs)[2], size_t count, uint32_t num_a, uint32_t num_b, bool self, md_allocator_i* alloc) {
    md_contact_set_t set = { .num_a = num_a, .num_b = num_b, .flags = self ? MD_CONTACT_FLAG_SELF : 0, .count = count, .alloc = alloc };
    set.i = md_alloc(alloc, sizeof(uint32_t) * MAX(count, 1));
    set.j = md_alloc(alloc, sizeof(uint32_t) * MAX(count, 1));
    set.atom_pairs = md_alloc(alloc, sizeof(uint32_t) * MAX(count, 1));
    set.d_min = md_alloc(alloc, sizeof(float) * MAX(count, 1));
    uint64_t* key = md_alloc(alloc, sizeof(uint64_t) * MAX(count, 1));
    for (size_t k = 0; k < count; ++k) {
        uint32_t a = pairs[k][0], b = pairs[k][1];
        if (self && a > b) { uint32_t t = a; a = b; b = t; }
        key[k] = ((uint64_t)a << 32) | b;
    }
    qsort(key, count, sizeof(uint64_t), cmp_u64_keys);
    for (size_t k = 0; k < count; ++k) {
        set.i[k] = (uint32_t)(key[k] >> 32);
        set.j[k] = (uint32_t)(key[k] & 0xFFFFFFFFu);
        set.atom_pairs[k] = 2;
        set.d_min[k] = 1.0f + (float)k;
    }
    return set;
}

// Chains of units in sequence: unit u is in body u / len, and neighbours u - 1 and u + 1 within its body
typedef struct chains_t {
    uint32_t* body;
    uint32_t* adj_off;
    uint32_t* adj;
    md_contact_units_t units;
} chains_t;

static chains_t make_chains(uint32_t num_bodies, uint32_t len, md_allocator_i* alloc) {
    const uint32_t n = num_bodies * len;
    chains_t c = {0};
    c.body = md_alloc(alloc, sizeof(uint32_t) * n);
    c.adj_off = md_alloc(alloc, sizeof(uint32_t) * (n + 1));
    c.adj = md_alloc(alloc, sizeof(uint32_t) * 2 * n);
    uint32_t e = 0;
    for (uint32_t u = 0; u < n; ++u) {
        c.body[u] = u / len;
        c.adj_off[u] = e;
        if (u % len > 0)       c.adj[e++] = u - 1;
        if (u % len < len - 1) c.adj[e++] = u + 1;
    }
    c.adj_off[n] = e;
    c.units = (md_contact_units_t){ .count = n, .body = c.body, .adj_off = c.adj_off, .adj = c.adj };
    return c;
}

// The number of regions, requiring every unit pair to have a region and the per region numbers to add up
static size_t check_regions(const md_contact_region_set_t* r, const md_contact_set_t* set) {
    if (r->num_pairs != set->count) return SIZE_MAX;
    size_t total = 0, atoms = 0;
    for (size_t k = 0; k < r->num_pairs; ++k) {
        if (r->region[k] >= r->count) return SIZE_MAX;
        atoms += set->atom_pairs[k];
    }
    size_t region_atoms = 0;
    for (size_t g = 0; g < r->count; ++g) {
        total += r->unit_pairs[g];
        region_atoms += r->atom_pairs[g];
        if (r->body_i[g] > r->body_j[g] && (set->flags & MD_CONTACT_FLAG_SELF)) return SIZE_MAX;
    }
    if (total != set->count || region_atoms != atoms) return SIZE_MAX;
    return r->count;
}

UTEST(contact_regions, patterns) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const uint32_t L = 20;
    chains_t ch = make_chains(3, L, alloc);
    md_contact_region_set_t r;

    // A parallel streak between chain 0 and chain 1: one region, eight units on either side
    {
        uint32_t p[8][2];
        for (uint32_t k = 0; k < 8; ++k) { p[k][0] = 2 + k; p[k][1] = L + 5 + k; }
        md_contact_set_t set = make_contact_set(p, 8, 3 * L, 3 * L, true, alloc);
        ASSERT_TRUE(md_contact_regions(&r, &set, &ch.units, NULL, 1, alloc));
        EXPECT_EQ((size_t)1, check_regions(&r, &set));
        EXPECT_EQ(0u, r.body_i[0]);
        EXPECT_EQ(1u, r.body_j[0]);
        EXPECT_EQ(8u, r.unit_pairs[0]);
        EXPECT_EQ(8u, r.units_i[0]);
        EXPECT_EQ(8u, r.units_j[0]);
        EXPECT_EQ(16u, r.atom_pairs[0]);
        EXPECT_EQ(1.0f, r.d_min[0]);
    }
    // Antiparallel: an anti-diagonal, still one region
    {
        uint32_t p[8][2];
        for (uint32_t k = 0; k < 8; ++k) { p[k][0] = 2 + k; p[k][1] = L + 15 - k; }
        md_contact_set_t set = make_contact_set(p, 8, 3 * L, 3 * L, true, alloc);
        ASSERT_TRUE(md_contact_regions(&r, &set, &ch.units, NULL, 1, alloc));
        EXPECT_EQ((size_t)1, check_regions(&r, &set));
    }
    // Two crossings of the same two chains: two compact blobs, two regions of two units per side
    {
        const uint32_t p[8][2] = {
            {3, L + 3}, {3, L + 4}, {4, L + 3}, {4, L + 4},
            {12, L + 15}, {12, L + 16}, {13, L + 15}, {13, L + 16},
        };
        md_contact_set_t set = make_contact_set(p, 8, 3 * L, 3 * L, true, alloc);
        ASSERT_TRUE(md_contact_regions(&r, &set, &ch.units, NULL, 1, alloc));
        ASSERT_EQ((size_t)2, check_regions(&r, &set));
        for (size_t g = 0; g < 2; ++g) {
            EXPECT_EQ(4u, r.unit_pairs[g]);
            EXPECT_EQ(2u, r.units_i[g]);
            EXPECT_EQ(2u, r.units_j[g]);
        }
        EXPECT_NE(r.region[0], r.region[7]);
    }
    // A streak with one unit out of contact: split with reach 1, bridged with reach 2
    {
        uint32_t p[7][2];
        const uint32_t along[7] = {2, 3, 4, 5, 7, 8, 9};
        for (uint32_t k = 0; k < 7; ++k) { p[k][0] = along[k]; p[k][1] = L + along[k]; }
        md_contact_set_t set = make_contact_set(p, 7, 3 * L, 3 * L, true, alloc);
        ASSERT_TRUE(md_contact_regions(&r, &set, &ch.units, NULL, 1, alloc));
        EXPECT_EQ((size_t)2, check_regions(&r, &set));
        ASSERT_TRUE(md_contact_regions(&r, &set, &ch.units, NULL, 2, alloc));
        EXPECT_EQ((size_t)1, check_regions(&r, &set));
    }
    // A chain folding onto itself (a hairpin): one region within body 0, six units, which are both sides
    {
        const uint32_t p[3][2] = { {2, 15}, {3, 14}, {4, 13} };
        md_contact_set_t set = make_contact_set(p, 3, 3 * L, 3 * L, true, alloc);
        ASSERT_TRUE(md_contact_regions(&r, &set, &ch.units, NULL, 1, alloc));
        ASSERT_EQ((size_t)1, check_regions(&r, &set));
        EXPECT_EQ(0u, r.body_i[0]);
        EXPECT_EQ(0u, r.body_j[0]);
        EXPECT_EQ(6u, r.units_i[0]);
        EXPECT_EQ(6u, r.units_j[0]);
    }
    // Neighbouring units in different bodies are not neighbours: the last unit of chain 0 and the first of chain 1
    // touching the same unit of chain 2 are contacts of different pairs of bodies
    {
        const uint32_t p[2][2] = { {L - 1, 2 * L + 5}, {L, 2 * L + 5} };
        md_contact_set_t set = make_contact_set(p, 2, 3 * L, 3 * L, true, alloc);
        // An adjacency which does join them across the bodies
        uint32_t* adj_off = md_alloc(alloc, sizeof(uint32_t) * (3 * L + 1));
        uint32_t* adj = md_alloc(alloc, sizeof(uint32_t) * 2);
        for (uint32_t u = 0; u <= 3 * L; ++u) adj_off[u] = (u <= L - 1) ? 0 : (u == L ? 1 : 2);
        adj[0] = L;
        adj[1] = L - 1;
        const md_contact_units_t units = { .count = 3 * L, .body = ch.body, .adj_off = adj_off, .adj = adj };
        ASSERT_TRUE(md_contact_regions(&r, &set, &units, NULL, 3, alloc));
        EXPECT_EQ((size_t)2, check_regions(&r, &set));
    }
    // Without an adjacency, every unit pair is a region of its own
    {
        uint32_t p[4][2];
        for (uint32_t k = 0; k < 4; ++k) { p[k][0] = 2 + k; p[k][1] = L + 2 + k; }
        md_contact_set_t set = make_contact_set(p, 4, 3 * L, 3 * L, true, alloc);
        const md_contact_units_t units = { .count = 3 * L, .body = ch.body };
        ASSERT_TRUE(md_contact_regions(&r, &set, &units, NULL, 1, alloc));
        EXPECT_EQ((size_t)4, check_regions(&r, &set));
    }
    // Between two sets of units, each numbered on its own
    {
        chains_t a = make_chains(1, L, alloc);
        chains_t b = make_chains(2, L, alloc);
        uint32_t p[6][2];
        for (uint32_t k = 0; k < 6; ++k) { p[k][0] = 4 + k; p[k][1] = L + 1 + k; }   // chain 0 of A against chain 1 of B
        md_contact_set_t set = make_contact_set(p, 6, L, 2 * L, false, alloc);
        ASSERT_TRUE(md_contact_regions(&r, &set, &a.units, &b.units, 1, alloc));
        ASSERT_EQ((size_t)1, check_regions(&r, &set));
        EXPECT_EQ(0u, r.body_i[0]);
        EXPECT_EQ(1u, r.body_j[0]);
        EXPECT_EQ(6u, r.units_i[0]);
        EXPECT_EQ(6u, r.units_j[0]);
        // Units which do not match the contact set are refused
        EXPECT_FALSE(md_contact_regions(&r, &set, &b.units, &b.units, 1, alloc));
        EXPECT_FALSE(md_contact_regions(&r, &set, &a.units, &b.units, 0, alloc));
    }

    md_arena_allocator_destroy(alloc);
}

// On a real system: residues as units of their chain, the adjacency from the peptide bonds
UTEST_F(contact, regions_residues) {
    md_system_t* sys = &utest_fixture->ala;
    size_t num = 0;
    md_bitfield_t* res = make_residue_groups(&num, sys, 0, SIZE_MAX, utest_fixture->alloc);
    ASSERT_GT(num, (size_t)3);

    uint32_t* adj_off = 0;
    uint32_t* adj = 0;
    ASSERT_TRUE(md_contact_unit_adjacency_from_bonds(&adj_off, &adj, res, num, sys, utest_fixture->alloc));
    // A chain: the ends have one neighbour, the others their predecessor and successor
    for (uint32_t u = 0; u < num; ++u) {
        const uint32_t expected = (u == 0 || u == num - 1) ? 1 : 2;
        ASSERT_EQ(expected, adj_off[u + 1] - adj_off[u]);
        if (u > 0) EXPECT_EQ(u - 1, adj[adj_off[u]]);
        if (u < num - 1) EXPECT_EQ(u + 1, adj[adj_off[u + 1] - 1]);
    }

    md_contact_desc_t desc = { .group_a = res, .num_a = num, .cutoff = 4.5, .exclude_bonds = 3 };
    md_contact_set_t set;
    ASSERT_TRUE(md_contact_compute(&set, &desc, sys, &utest_fixture->ala_state, utest_fixture->alloc));

    uint32_t* body = md_alloc(utest_fixture->alloc, sizeof(uint32_t) * num);
    memset(body, 0, sizeof(uint32_t) * num);
    const md_contact_units_t units = { .count = num, .body = body, .adj_off = adj_off, .adj = adj };
    md_contact_region_set_t r;
    ASSERT_TRUE(md_contact_regions(&r, &set, &units, NULL, 1, utest_fixture->alloc));
    EXPECT_NE(SIZE_MAX, check_regions(&r, &set));

    // Consecutive residues touch along the whole chain: those contacts form one region
    uint32_t backbone = UINT32_MAX;
    size_t found = 0;
    for (size_t k = 0; k < set.count; ++k) {
        if (set.j[k] == set.i[k] + 1) {
            if (backbone == UINT32_MAX) backbone = r.region[k];
            EXPECT_EQ(backbone, r.region[k]);
            found += 1;
        }
    }
    EXPECT_EQ(num - 1, found);
}
