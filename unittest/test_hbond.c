#include "utest.h"

#include <md_hbond.h>
#include <md_system.h>
#include <md_gro.h>
#include <md_pdb.h>
#include <md_xyz.h>
#include <md_util.h>
#include <md_unitcell.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_bitfield.h>
#include <core/md_str.h>
#include <core/md_vec_math.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

// ### HELPERS ###

static bool load_xyz(md_system_t* sys, md_system_state_t* st, const char* text, const md_unitcell_t* cell, md_allocator_i* alloc) {
    *sys = (md_system_t){ .alloc = alloc };
    *st  = (md_system_state_t){ .alloc = alloc };
    if (!md_xyz_system_init_from_str(sys, st, str_from_cstr(text), MD_XYZ_OPTION_NONE)) return false;
    if (cell) st->unitcell = *cell;
    return md_util_system_infer(sys, st, MD_UTIL_INFER_ALL);
}

static inline uint64_t hb_key(uint32_t h, uint32_t a) { return ((uint64_t)h << 32) | a; }

static int compare_u64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
    return (x > y) - (x < y);
}

static bool contains_u64(const uint64_t* arr, size_t n, uint64_t key) {
    return n && bsearch(&key, arr, n, sizeof(uint64_t), compare_u64) != NULL;
}

static bool set_is_sorted_by_key(const md_hbond_set_t* s) {
    for (size_t k = 1; k < s->count; ++k) {
        if (hb_key(s->hydrogen[k - 1], s->acceptor[k - 1]) >= hb_key(s->hydrogen[k], s->acceptor[k])) return false;
    }
    return true;
}

static md_bitfield_t residue_set(const md_system_t* sys, const char* resname, size_t max_count, md_allocator_i* alloc) {
    md_bitfield_t bf = md_bitfield_create(alloc);
    size_t n = 0;
    for (size_t c = 0; c < sys->component.count && n < max_count; ++c) {
        str_t name = str_trim(md_component_name(&sys->component, c));
        if (resname && !str_eq_cstr(name, resname)) continue;
        md_urange_t r = md_system_component_atom_range(sys, c);
        md_bitfield_set_range(&bf, r.beg, r.end);
        n += 1;
    }
    return bf;
}

static vec3_t mi_diff(const md_system_state_t* st, uint32_t from, uint32_t to) {
    vec3_t d = vec3_sub(st->xyz[to], st->xyz[from]);
    md_util_min_image_vec3(&d, 1, &st->unitcell);
    return d;
}

static double angle_deg(vec3_t a, vec3_t b) {
    const double la = sqrt((double)vec3_dot(a, a));
    const double lb = sqrt((double)vec3_dot(b, b));
    double c = (double)vec3_dot(a, b) / (la * lb);
    c = c < -1 ? -1 : (c > 1 ? 1 : c);
    return acos(c) * 180.0 / 3.14159265358979323846;
}

// Brute force over every 'stride'th donor against every acceptor of the query, gates only (no competition, no
// selection), in double precision with the general minimum image and the bond graph walked anew. Pairs within eps of
// a gate are 'border': either answer is acceptable for them.
static void reference_hbonds(md_array(uint64_t)* out_pass, md_array(uint64_t)* out_border, const md_hbond_query_t* q, const md_system_t* sys, const md_system_state_t* st, size_t stride, md_allocator_i* alloc) {
    const md_hbond_params_t* p = &q->params;
    const double eps_d = 1e-3, eps_a = 1e-2;
    md_array(uint32_t) near = 0;
    md_array(uint32_t) queue = 0;
    uint8_t* depth = md_alloc(alloc, sys->atom.count);
    memset(depth, 0, sys->atom.count);

    for (size_t k = 0; k < q->num_donors; k += stride) {
        const uint32_t D = q->donor_d[k], H = q->donor_h[k];
        // Atoms within exclude_bonds of D
        md_array_shrink(near, 0);
        md_array_shrink(queue, 0);
        if (p->exclude_bonds) {
            depth[D] = 1;
            md_array_push(queue, D, alloc);
            for (size_t head = 0; head < md_array_size(queue); ++head) {
                const uint32_t cur = queue[head];
                if (depth[cur] > p->exclude_bonds) continue;
                md_bond_iter_t it = md_bond_iter(&sys->bond, cur);
                while (md_bond_iter_has_next(&it)) {
                    const uint32_t nx = (uint32_t)md_bond_iter_atom_index(&it);
                    if (!depth[nx]) { depth[nx] = depth[cur] + 1; md_array_push(queue, nx, alloc); md_array_push(near, nx, alloc); }
                    md_bond_iter_next(&it);
                }
            }
            for (size_t v = 0; v < md_array_size(queue); ++v) depth[queue[v]] = 0;
        }
        const vec3_t v_dh = mi_diff(st, D, H);
        for (size_t j = 0; j < q->num_acceptors; ++j) {
            const uint32_t A = q->acceptor[j];
            if (A == D || A == H) continue;
            bool excl = false;
            for (size_t n = 0; n < md_array_size(near); ++n) if (near[n] == A) { excl = true; break; }
            if (excl) continue;
            const vec3_t v_da = mi_diff(st, D, A);
            const double r_da = sqrt((double)vec3_dot(v_da, v_da));
            if (r_da > 6.0) continue;
            const vec3_t v_ha = vec3_sub(v_da, v_dh);
            const double r_ha = sqrt((double)vec3_dot(v_ha, v_ha));
            bool pass = true, border = false;
            if (p->max_ha > 0) { if (r_ha > p->max_ha) pass = false; if (fabs(r_ha - p->max_ha) < eps_d) border = true; }
            if (p->max_da > 0) { if (r_da > p->max_da) pass = false; if (fabs(r_da - p->max_da) < eps_d) border = true; }
            const double dha = angle_deg(vec3_mul1(v_dh, -1.0f), v_ha);
            if (p->min_dha > 0) { if (dha < p->min_dha) pass = false; if (fabs(dha - p->min_dha) < eps_a) border = true; }
            const double hda = angle_deg(v_dh, v_da);
            if (p->max_hda > 0) { if (hda > p->max_hda) pass = false; if (fabs(hda - p->max_hda) < eps_a) border = true; }
            if (p->min_xah > 0) {
                md_bond_iter_t it = md_bond_iter(&sys->bond, A);
                while (md_bond_iter_has_next(&it)) {
                    const uint32_t X = (uint32_t)md_bond_iter_atom_index(&it);
                    md_bond_iter_next(&it);
                    if (X == H || md_atom_atomic_number(&sys->atom, X) == 0) continue;
                    const double xah = angle_deg(mi_diff(st, A, X), vec3_mul1(v_ha, -1.0f));
                    if (xah < p->min_xah) pass = false;
                    if (fabs(xah - p->min_xah) < eps_a) border = true;
                }
            }
            if (border) md_array_push(*out_border, hb_key(H, A), alloc);
            else if (pass) md_array_push(*out_pass, hb_key(H, A), alloc);
        }
    }
    if (md_array_size(*out_pass))   qsort(*out_pass, md_array_size(*out_pass), sizeof(uint64_t), compare_u64);
    if (md_array_size(*out_border)) qsort(*out_border, md_array_size(*out_border), sizeof(uint64_t), compare_u64);
}

// pass ⊆ result ⊆ pass ∪ border, over the bonds of the hydrogens the reference covered
static bool matches_reference(const md_hbond_set_t* s, const uint8_t* covered_h, const uint64_t* pass, size_t num_pass, const uint64_t* border, size_t num_border) {
    size_t found = 0;
    for (size_t k = 0; k < s->count; ++k) {
        if (!covered_h[s->hydrogen[k]]) continue;
        const uint64_t key = hb_key(s->hydrogen[k], s->acceptor[k]);
        if (contains_u64(pass, num_pass, key)) { found += 1; continue; }
        if (contains_u64(border, num_border, key)) continue;
        printf("Bond H %u ... A %u not in the reference (r_ha %.3f, r_da %.3f, dha %.2f)\n", s->hydrogen[k], s->acceptor[k], s->dist_ha[k], s->dist_da[k], s->angle_dha[k]);
        return false;
    }
    if (found != num_pass) {
        printf("%zu of %zu reference bonds missing\n", num_pass - found, num_pass);
        return false;
    }
    return true;
}

// ### FIXTURE ###

struct hbond {
    md_allocator_i* alloc;
    md_system_t npt;            // Triclinic cell, protein, ligand, Cl- in water
    md_system_state_t npt_state;
    md_system_t centered;       // Large orthorhombic water box with a protein
    md_system_state_t centered_state;
};

UTEST_F_SETUP(hbond) {
    utest_fixture->alloc = md_vm_arena_create(GIGABYTES(8));
    md_allocator_i* alloc = utest_fixture->alloc;

    utest_fixture->npt = (md_system_t){ .alloc = alloc };
    utest_fixture->npt_state = (md_system_state_t){ .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&utest_fixture->npt, &utest_fixture->npt_state, STR_LIT(MD_UNITTEST_DATA_DIR "/npt.gro")));
    ASSERT_TRUE(md_util_system_infer(&utest_fixture->npt, &utest_fixture->npt_state, MD_UTIL_INFER_ALL));

    utest_fixture->centered = (md_system_t){ .alloc = alloc };
    utest_fixture->centered_state = (md_system_state_t){ .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&utest_fixture->centered, &utest_fixture->centered_state, STR_LIT(MD_UNITTEST_DATA_DIR "/centered.gro")));
    ASSERT_TRUE(md_util_system_infer(&utest_fixture->centered, &utest_fixture->centered_state, MD_UTIL_INFER_ALL));
}

UTEST_F_TEARDOWN(hbond) {
    md_vm_arena_destroy(utest_fixture->alloc);
}

// ### SMALL SYSTEMS WITH KNOWN ANSWERS ###

// O1-H1 points along +x at O2, 2.9 Å away; the hydrogens of the second water point away from the first.
static const char* water_dimer =
    "6\n"
    "water dimer\n"
    "O  0.000  0.000  0.000\n"
    "H  0.957  0.000  0.000\n"
    "H -0.240  0.927  0.000\n"
    "O  2.900  0.000  0.000\n"
    "H  3.140  0.927  0.000\n"
    "H  3.140 -0.464  0.803\n";

UTEST(hbond, water_dimer) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys;
    md_system_state_t st;
    ASSERT_TRUE(load_xyz(&sys, &st, water_dimer, NULL, alloc));
    ASSERT_EQ(4, (int)sys.bond.count);

    for (int pr = 0; pr < MD_HBOND_PRESET_COUNT; ++pr) {
        md_hbond_params_t params = md_hbond_params_preset((md_hbond_preset_t)pr);
        md_hbond_desc_t desc = { .params = &params };
        md_hbond_set_t set;
        ASSERT_TRUE(md_hbond_compute(&set, &desc, &sys, &st, alloc));
        EXPECT_EQ(1, (int)set.count);
        if (set.count == 1) {
            EXPECT_EQ(0u, set.donor[0]);
            EXPECT_EQ(1u, set.hydrogen[0]);
            EXPECT_EQ(3u, set.acceptor[0]);
            EXPECT_NEAR(2.9f, set.dist_da[0], 1e-3f);
            EXPECT_NEAR(1.943f, set.dist_ha[0], 1e-3f);
            EXPECT_NEAR(180.0f, set.angle_dha[0], 0.1f);
            // Linear and short: strong. Under the 3.0 Å D...A gates (MDAnalysis, VMD) 2.9 Å is near the gate.
            EXPECT_GT(set.strength[0], params.max_da == 3.0f ? 0.4f : 0.9f);
            EXPECT_LE(set.strength[0], 1.0f);
        }
    }
    md_vm_arena_destroy(alloc);
}

// The same dimer straddling the boundary of a periodic cell: the donor's own O-H is split, and so is H...A.
UTEST(hbond, water_dimer_periodic) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    const char* shifted =
        "6\n"
        "water dimer across the boundary at x = 10\n"
        "O  9.600  0.000  5.000\n"
        "H  0.557  0.000  5.000\n"
        "H  9.360  0.927  5.000\n"
        "O  2.500  0.000  5.000\n"
        "H  2.740  0.927  5.000\n"
        "H  2.740 -0.464  5.803\n";
    const md_unitcell_t cell = md_unitcell_from_extent(10, 10, 10);
    md_system_t sys;
    md_system_state_t st;
    ASSERT_TRUE(load_xyz(&sys, &st, shifted, &cell, alloc));
    ASSERT_EQ(4, (int)sys.bond.count);

    md_hbond_set_t set;
    ASSERT_TRUE(md_hbond_compute(&set, NULL, &sys, &st, alloc));
    ASSERT_EQ(1, (int)set.count);
    EXPECT_EQ(0u, set.donor[0]);
    EXPECT_EQ(1u, set.hydrogen[0]);
    EXPECT_EQ(3u, set.acceptor[0]);
    EXPECT_NEAR(2.9f, set.dist_da[0], 1e-3f);
    EXPECT_NEAR(180.0f, set.angle_dha[0], 0.1f);
    md_vm_arena_destroy(alloc);
}

// One hydrogen between two acceptors: a bifurcated bond. One bond per hydrogen keeps the stronger, two keep both, and
// the tolerance drops the second when it is much weaker.
UTEST(hbond, bifurcation) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    const char* text =
        "9\n"
        "water donating to two waters\n"
        "O   0.000  0.000  0.000\n"
        "H   0.957  0.000  0.000\n"
        "H  -0.240  0.927  0.000\n"
        "O   2.700  1.300  0.000\n"
        "H   3.300  1.900  0.000\n"
        "H   3.300  1.300  0.800\n"
        "O   2.900 -1.500  0.000\n"
        "H   3.500 -2.100  0.000\n"
        "H   3.500 -1.500 -0.800\n";
    md_system_t sys;
    md_system_state_t st;
    ASSERT_TRUE(load_xyz(&sys, &st, text, NULL, alloc));

    md_hbond_params_t p = md_hbond_params_preset(MD_HBOND_PRESET_REALISTIC);
    p.max_ha = 2.6f;
    p.min_dha = 110.0f;
    md_hbond_desc_t desc = { .params = &p };
    md_hbond_set_t set;

    ASSERT_TRUE(md_hbond_compute(&set, &desc, &sys, &st, alloc));
    ASSERT_EQ(1, (int)set.count);
    EXPECT_EQ(3u, set.acceptor[0]);     // The closer one

    p.h_capacity = 2;
    ASSERT_TRUE(md_hbond_compute(&set, &desc, &sys, &st, alloc));
    ASSERT_EQ(2, (int)set.count);
    const float s0 = set.strength[0], s1 = set.strength[1];
    const float ratio = MIN(s0, s1) / MAX(s0, s1);

    p.bifurcation_tol = MIN(1.0f, ratio + 0.05f);
    ASSERT_TRUE(md_hbond_compute(&set, &desc, &sys, &st, alloc));
    EXPECT_EQ(1, (int)set.count);

    p.bifurcation_tol = MAX(0.0f, ratio - 0.05f);
    ASSERT_TRUE(md_hbond_compute(&set, &desc, &sys, &st, alloc));
    EXPECT_EQ(2, (int)set.count);
    md_vm_arena_destroy(alloc);
}

// A pyramidal NH3 has a free lone pair, a flattened one does not
UTEST(hbond, roles_planarity) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    const char* text =
        "8\n"
        "pyramidal and planar NH3\n"
        "N   0.000  0.000  0.000\n"
        "H   0.940  0.000 -0.380\n"
        "H  -0.470  0.814 -0.380\n"
        "H  -0.470 -0.814 -0.380\n"
        "N  10.000  0.000  0.000\n"
        "H  11.010  0.000  0.000\n"
        "H   9.495  0.875  0.000\n"
        "H   9.495 -0.875  0.000\n";
    md_system_t sys;
    md_system_state_t st;
    ASSERT_TRUE(load_xyz(&sys, &st, text, NULL, alloc));
    uint8_t role[8], cap[8];
    ASSERT_TRUE(md_hbond_perceive_roles(role, cap, &sys, &st, MD_HBOND_ROLES_DEFAULT));
    EXPECT_EQ(MD_HBOND_ROLE_DONOR | MD_HBOND_ROLE_ACCEPTOR, (int)role[0]);
    EXPECT_EQ(1, (int)cap[0]);
    EXPECT_EQ(MD_HBOND_ROLE_DONOR, (int)role[4]);
    EXPECT_EQ(0, (int)cap[4]);

    // The atom flags of md_util_system_infer carry the same decision to perception without coordinates
    EXPECT_TRUE(sys.atom.flags[0] & MD_FLAG_HBOND_ACCEPTOR);
    EXPECT_FALSE(sys.atom.flags[4] & MD_FLAG_HBOND_ACCEPTOR);
    ASSERT_TRUE(md_hbond_perceive_roles(role, cap, &sys, NULL, MD_HBOND_ROLES_DEFAULT));
    EXPECT_EQ(MD_HBOND_ROLE_DONOR | MD_HBOND_ROLE_ACCEPTOR, (int)role[0]);
    EXPECT_EQ(MD_HBOND_ROLE_DONOR, (int)role[4]);

    // Every N and O
    ASSERT_TRUE(md_hbond_perceive_roles(role, cap, &sys, &st, MD_HBOND_ROLES_ALL_N_O));
    EXPECT_TRUE(role[4] & MD_HBOND_ROLE_ACCEPTOR);
    md_vm_arena_destroy(alloc);
}

UTEST(hbond, no_hydrogens) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/1k4r.pdb"), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    md_hbond_query_t q;
    ASSERT_TRUE(md_hbond_query_init(&q, NULL, &sys, alloc));
    EXPECT_EQ(0, (int)q.num_donors);
    EXPECT_TRUE(q.flags & MD_HBOND_FLAG_NO_HYDROGENS);
    md_hbond_set_t set;
    ASSERT_TRUE(md_hbond_query_eval(&set, &q, &st, alloc));
    EXPECT_EQ(0, (int)set.count);
    EXPECT_TRUE(set.flags & MD_HBOND_FLAG_NO_HYDROGENS);

    // Without hydrogens the standard residues still get their acceptors right: no backbone N
    for (size_t k = 0; k < q.num_acceptors; ++k) {
        const uint32_t a = q.acceptor[k];
        if (md_atom_flags(&sys.atom, a) & MD_FLAG_AMINO_ACID) {
            EXPECT_FALSE(str_eq_cstr(str_trim(md_atom_name(&sys.atom, a)), "N"));
        }
    }
    md_vm_arena_destroy(alloc);
}

UTEST(hbond, invalid_params) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys;
    md_system_state_t st;
    ASSERT_TRUE(load_xyz(&sys, &st, water_dimer, NULL, alloc));
    md_hbond_query_t q;
    md_hbond_params_t p = { 0 };       // No distance gate
    md_hbond_desc_t desc = { .params = &p };
    EXPECT_FALSE(md_hbond_query_init(&q, &desc, &sys, alloc));
    p.max_ha = 2.5f;
    p.min_dha = 200.0f;
    EXPECT_FALSE(md_hbond_query_init(&q, &desc, &sys, alloc));
    p.min_dha = 120.0f;
    EXPECT_TRUE(md_hbond_query_init(&q, &desc, &sys, alloc));
    md_hbond_query_free(&q);
    md_vm_arena_destroy(alloc);
}

// ### ROLES OF REAL SYSTEMS ###

UTEST_F(hbond, roles_protein) {
    const md_system_t* sys = &utest_fixture->centered;
    const size_t N = sys->atom.count;
    md_allocator_i* alloc = utest_fixture->alloc;
    uint8_t* role = md_alloc(alloc, N);
    uint8_t* cap  = md_alloc(alloc, N);
    ASSERT_TRUE(md_hbond_perceive_roles(role, cap, sys, &utest_fixture->centered_state, MD_HBOND_ROLES_DEFAULT));

    size_t his_free_n = 0, his_acc = 0, lys_nz_donor = 0, bb_o = 0, bb_o_acc = 0;
    for (size_t c = 0; c < sys->component.count; ++c) {
        const str_t res = str_trim(md_component_name(&sys->component, c));
        const md_urange_t r = md_system_component_atom_range(sys, c);
        for (uint32_t i = r.beg; i < r.end; ++i) {
            const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, i);
            const str_t name = str_trim(md_atom_name(&sys->atom, i));
            const bool amino = (md_atom_flags(&sys->atom, i) & MD_FLAG_AMINO_ACID) != 0;
            if (amino && z == MD_Z_N) {
                int nh = 0;
                md_bond_iter_t it = md_bond_iter(&sys->bond, i);
                while (md_bond_iter_has_next(&it)) { nh += md_atom_atomic_number(&sys->atom, md_bond_iter_atom_index(&it)) == MD_Z_H; md_bond_iter_next(&it); }
                const bool his_ring = str_eq_cstr(res, "HIS") && (str_eq_cstr(name, "ND1") || str_eq_cstr(name, "NE2"));
                if (role[i] & MD_HBOND_ROLE_ACCEPTOR) {
                    // The only protein N with a free lone pair is the unprotonated ring N of histidine
                    EXPECT_TRUE(his_ring && nh == 0);
                    his_acc += 1;
                }
                if (his_ring && nh == 0) his_free_n += 1;
                if (str_eq_cstr(res, "LYS") && str_eq_cstr(name, "NZ") && (role[i] & MD_HBOND_ROLE_DONOR)) lys_nz_donor += 1;
                EXPECT_EQ(nh > 0, (role[i] & MD_HBOND_ROLE_DONOR) != 0);
            }
            if (amino && z == MD_Z_O && str_eq_cstr(name, "O")) {
                bb_o += 1;
                if ((role[i] & MD_HBOND_ROLE_ACCEPTOR) && cap[i] == 2) bb_o_acc += 1;
            }
        }
    }
    EXPECT_GT(his_free_n, (size_t)0);
    EXPECT_EQ(his_free_n, his_acc);
    EXPECT_GT(lys_nz_donor, (size_t)0);
    EXPECT_GT(bb_o, (size_t)0);
    EXPECT_EQ(bb_o, bb_o_acc);
}

UTEST_F(hbond, roles_water) {
    const md_system_t* sys = &utest_fixture->npt;
    const size_t N = sys->atom.count;
    uint8_t* role = md_alloc(utest_fixture->alloc, N);
    uint8_t* cap  = md_alloc(utest_fixture->alloc, N);
    ASSERT_TRUE(md_hbond_perceive_roles(role, cap, sys, &utest_fixture->npt_state, MD_HBOND_ROLES_DEFAULT));
    size_t water_o = 0;
    for (size_t c = 0; c < sys->component.count; ++c) {
        if (!str_eq_cstr(str_trim(md_component_name(&sys->component, c)), "SOL")) continue;
        const md_urange_t r = md_system_component_atom_range(sys, c);
        for (uint32_t i = r.beg; i < r.end; ++i) {
            if (md_atom_atomic_number(&sys->atom, i) != MD_Z_O) continue;
            EXPECT_EQ(MD_HBOND_ROLE_DONOR | MD_HBOND_ROLE_ACCEPTOR, (int)role[i]);
            EXPECT_EQ(2, (int)cap[i]);
            water_o += 1;
        }
    }
    EXPECT_GT(water_o, (size_t)1000);
}

UTEST(hbond, roles_nucleic) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(2));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/nucl-dna.gro")));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    const size_t N = sys.atom.count;
    uint8_t* role = md_alloc(alloc, N);
    ASSERT_TRUE(md_hbond_perceive_roles(role, NULL, &sys, &st, MD_HBOND_ROLES_DEFAULT));

    // The N acceptors of each base and nothing else: A N1 N3 N7, G N3 N7, C N3, T none
    size_t checked = 0;
    for (size_t c = 0; c < sys.component.count; ++c) {
        const str_t res = str_trim(md_component_name(&sys.component, c));
        if (res.len < 2 || res.ptr[0] != 'D') continue;
        const char base = res.ptr[1];
        const md_urange_t r = md_system_component_atom_range(&sys, c);
        for (uint32_t i = r.beg; i < r.end; ++i) {
            if (md_atom_atomic_number(&sys.atom, i) != MD_Z_N) continue;
            const str_t name = str_trim(md_atom_name(&sys.atom, i));
            bool expect = false;
            switch (base) {
            case 'A': expect = str_eq_cstr(name, "N1") || str_eq_cstr(name, "N3") || str_eq_cstr(name, "N7"); break;
            case 'G': expect = str_eq_cstr(name, "N3") || str_eq_cstr(name, "N7"); break;
            case 'C': expect = str_eq_cstr(name, "N3"); break;
            default:  expect = false; break;
            }
            EXPECT_EQ(expect, (role[i] & MD_HBOND_ROLE_ACCEPTOR) != 0);
            checked += 1;
        }
    }
    EXPECT_GT(checked, (size_t)1000);
    md_vm_arena_destroy(alloc);
}

UTEST_F(hbond, roles_halide) {
    const md_system_t* sys = &utest_fixture->npt;
    md_hbond_query_t q;
    ASSERT_TRUE(md_hbond_query_init(&q, NULL, sys, utest_fixture->alloc));
    size_t num_cl = 0;
    for (size_t k = 0; k < q.num_acceptors; ++k) {
        if (md_atom_atomic_number(&sys->atom, q.acceptor[k]) == MD_Z_Cl) {
            EXPECT_EQ(MD_HBOND_CAPACITY_NO_LIMIT, (int)q.acceptor_cap[k]);
            num_cl += 1;
        }
    }
    EXPECT_GT(num_cl, (size_t)0);

    // A chloride in water takes several bonds
    md_hbond_set_t set;
    ASSERT_TRUE(md_hbond_query_eval(&set, &q, &utest_fixture->npt_state, utest_fixture->alloc));
    int max_per_cl = 0;
    for (size_t k = 0; k < set.count; ) {
        size_t n = 0;
        const uint32_t a = set.acceptor[k];
        for (size_t m = 0; m < set.count; ++m) n += set.acceptor[m] == a;
        if (md_atom_atomic_number(&sys->atom, a) == MD_Z_Cl) max_per_cl = MAX(max_per_cl, (int)n);
        ++k;
    }
    EXPECT_GE(max_per_cl, 3);
}

// ### AGAINST A BRUTE FORCE REFERENCE ###

static void check_preset_against_reference(int* utest_result, const md_system_t* sys, const md_system_state_t* st, md_hbond_params_t params, md_allocator_i* alloc) {
    // The reference has no competition
    params.h_capacity = 0;
    params.acc_capacity_mode = MD_HBOND_CAPACITY_UNLIMITED;
    md_hbond_desc_t desc = { .params = &params };
    md_hbond_query_t q;
    ASSERT_TRUE(md_hbond_query_init(&q, &desc, sys, alloc));
    md_hbond_set_t set;
    ASSERT_TRUE(md_hbond_query_eval(&set, &q, st, alloc));
    EXPECT_TRUE(set_is_sorted_by_key(&set));

    // Every third donor keeps the brute force affordable and still covers thousands of bonds
    const size_t stride = 3;
    uint8_t* covered_h = md_alloc(alloc, sys->atom.count);
    memset(covered_h, 0, sys->atom.count);
    for (size_t k = 0; k < q.num_donors; k += stride) covered_h[q.donor_h[k]] = 1;

    md_array(uint64_t) pass = 0;
    md_array(uint64_t) border = 0;
    reference_hbonds(&pass, &border, &q, sys, st, stride, alloc);
    EXPECT_GT(md_array_size(pass), (size_t)100);
    EXPECT_TRUE(matches_reference(&set, covered_h, pass, md_array_size(pass), border, md_array_size(border)));
}

UTEST_F(hbond, reference_triclinic) {
    for (int pr = 0; pr < MD_HBOND_PRESET_COUNT; ++pr) {
        check_preset_against_reference(utest_result, &utest_fixture->npt, &utest_fixture->npt_state, md_hbond_params_preset((md_hbond_preset_t)pr), utest_fixture->alloc);
    }
}

UTEST_F(hbond, reference_sulfur_fluorine) {
    md_hbond_params_t p = md_hbond_params_preset(MD_HBOND_PRESET_REALISTIC);
    p.roles |= MD_HBOND_ROLES_SULFUR | MD_HBOND_ROLES_FLUORINE;
    p.exclude_bonds = 5;
    p.max_da = 3.3f;
    check_preset_against_reference(utest_result, &utest_fixture->npt, &utest_fixture->npt_state, p, utest_fixture->alloc);
}

// ### COMPETITION ###

UTEST_F(hbond, competition_invariants) {
    const md_system_t* sys = &utest_fixture->npt;
    md_allocator_i* alloc = utest_fixture->alloc;
    md_hbond_params_t p = md_hbond_params_preset(MD_HBOND_PRESET_REALISTIC);
    md_hbond_desc_t desc = { .params = &p };
    md_hbond_query_t q;
    ASSERT_TRUE(md_hbond_query_init(&q, &desc, sys, alloc));
    md_hbond_set_t capped;
    ASSERT_TRUE(md_hbond_query_eval(&capped, &q, &utest_fixture->npt_state, alloc));

    md_hbond_params_t p_all = p;
    p_all.h_capacity = 0;
    p_all.acc_capacity_mode = MD_HBOND_CAPACITY_UNLIMITED;
    md_hbond_desc_t desc_all = { .params = &p_all };
    md_hbond_set_t all;
    ASSERT_TRUE(md_hbond_compute(&all, &desc_all, sys, &utest_fixture->npt_state, alloc));
    EXPECT_GT(capped.count, (size_t)100);
    EXPECT_LE(capped.count, all.count);

    md_array(uint64_t) keys = 0;
    for (size_t k = 0; k < all.count; ++k) md_array_push(keys, hb_key(all.hydrogen[k], all.acceptor[k]), alloc);
    for (size_t k = 0; k < capped.count; ++k) {
        EXPECT_TRUE(contains_u64(keys, md_array_size(keys), hb_key(capped.hydrogen[k], capped.acceptor[k])));
        if (k) EXPECT_NE(capped.hydrogen[k - 1], capped.hydrogen[k]);   // One bond per hydrogen
    }
    // No acceptor beyond its lone pairs
    for (size_t j = 0; j < q.num_acceptors; ++j) {
        if (q.acceptor_cap[j] == MD_HBOND_CAPACITY_NO_LIMIT) continue;
        int n = 0;
        for (size_t k = 0; k < capped.count; ++k) n += capped.acceptor[k] == q.acceptor[j];
        EXPECT_LE(n, (int)q.acceptor_cap[j]);
    }
}

// Every N and O as acceptor through the override equals the ALL_N_O roles
UTEST_F(hbond, role_override) {
    const md_system_t* sys = &utest_fixture->npt;
    md_allocator_i* alloc = utest_fixture->alloc;
    md_bitfield_t n_o = md_bitfield_create(alloc);
    for (size_t i = 0; i < sys->atom.count; ++i) {
        const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, i);
        if (z == MD_Z_N || z == MD_Z_O) md_bitfield_set_bit(&n_o, i);
    }
    md_hbond_params_t p = md_hbond_params_preset(MD_HBOND_PRESET_MDTRAJ);
    md_hbond_desc_t desc = { .params = &p };
    md_hbond_set_t ref;
    ASSERT_TRUE(md_hbond_compute(&ref, &desc, sys, &utest_fixture->npt_state, alloc));

    p.roles = MD_HBOND_ROLES_DEFAULT;
    desc.acceptors = &n_o;
    md_hbond_set_t set;
    ASSERT_TRUE(md_hbond_compute(&set, &desc, sys, &utest_fixture->npt_state, alloc));
    ASSERT_EQ(ref.count, set.count);
    for (size_t k = 0; k < set.count; ++k) {
        EXPECT_EQ(ref.hydrogen[k], set.hydrogen[k]);
        EXPECT_EQ(ref.acceptor[k], set.acceptor[k]);
    }
}

// ### SELECTIONS ###

static void check_selection_independent(int* utest_result, const md_system_t* sys, const md_system_state_t* st, const md_bitfield_t* a, const md_bitfield_t* b, md_allocator_i* alloc) {
    md_hbond_set_t full;
    ASSERT_TRUE(md_hbond_compute(&full, NULL, sys, st, alloc));
    md_hbond_desc_t desc = { .set_a = a, .set_b = b };
    md_hbond_set_t sel;
    ASSERT_TRUE(md_hbond_compute(&sel, &desc, sys, st, alloc));
    EXPECT_TRUE(sel.flags & MD_HBOND_FLAG_SELECTION);

    // The bonds of the full system which the selection reports
    size_t expected = 0, m = 0;
    bool same = true;
    for (size_t k = 0; k < full.count; ++k) {
        const bool da = md_bitfield_test_bit(a, full.donor[k]) || md_bitfield_test_bit(a, full.hydrogen[k]);
        const bool aa = md_bitfield_test_bit(a, full.acceptor[k]);
        bool keep;
        if (b) {
            const bool db = md_bitfield_test_bit(b, full.donor[k]) || md_bitfield_test_bit(b, full.hydrogen[k]);
            const bool ab = md_bitfield_test_bit(b, full.acceptor[k]);
            keep = (da && ab) || (db && aa);
        } else {
            keep = da && aa;
        }
        if (!keep) continue;
        expected += 1;
        if (m < sel.count && sel.hydrogen[m] == full.hydrogen[k] && sel.acceptor[m] == full.acceptor[k] && sel.strength[m] == full.strength[k]) {
            m += 1;
        } else {
            same = false;
        }
    }
    EXPECT_GT(expected, (size_t)0);
    EXPECT_EQ(expected, sel.count);
    EXPECT_TRUE(same);
}

UTEST_F(hbond, selection_between) {
    md_allocator_i* alloc = utest_fixture->alloc;
    const md_system_t* sys = &utest_fixture->centered;
    // A few residues against the rest: a small selection, which exercises the shell of competitors
    md_bitfield_t few  = residue_set(sys, NULL, 60, alloc);
    md_bitfield_t rest = md_bitfield_create(alloc);
    md_bitfield_not(&rest, &few, 0, sys->atom.count);
    check_selection_independent(utest_result, sys, &utest_fixture->centered_state, &few, &rest, alloc);
    // Within a small selection
    check_selection_independent(utest_result, sys, &utest_fixture->centered_state, &few, NULL, alloc);
    // Within most of the system
    check_selection_independent(utest_result, sys, &utest_fixture->centered_state, &rest, NULL, alloc);
}

UTEST_F(hbond, selection_ligand) {
    md_allocator_i* alloc = utest_fixture->alloc;
    const md_system_t* sys = &utest_fixture->npt;
    md_bitfield_t lys = residue_set(sys, "LYS", 1000, alloc);
    md_bitfield_t sol = residue_set(sys, "SOL", 100000, alloc);
    check_selection_independent(utest_result, sys, &utest_fixture->npt_state, &lys, &sol, alloc);
}

// A hydrogen of the selection between two acceptors: the stronger outside the selection, the weaker inside. With one
// bond per hydrogen the full system bonds it to the outside acceptor, so the selection must report nothing, not the
// weaker bond that would win if the competitor were ignored. Distant waters make the selection a small part of the
// system, which is when the shell of competitors is built rather than the whole system taken.
UTEST(hbond, selection_competitor_outside) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    char text[8192];
    int len = 0;
    const int num_far = 20;
    len += snprintf(text + len, sizeof(text) - len, "%d\ncompetitor outside the selection\n", 9 + 3 * num_far);
    len += snprintf(text + len, sizeof(text) - len,
        "O   0.000  0.000  0.000\n"
        "H   0.957  0.000  0.000\n"
        "H  -0.240  0.927  0.000\n"
        "O   2.700  1.300  0.000\n"     // Stronger, outside the selection
        "H   3.300  1.900  0.000\n"
        "H   3.300  1.300  0.800\n"
        "O   2.900 -1.500  0.000\n"     // Weaker, inside
        "H   3.500 -2.100  0.000\n"
        "H   3.500 -1.500 -0.800\n");
    for (int w = 0; w < num_far; ++w) {
        const float x = 30.0f + 6.0f * (float)(w % 5), y = 30.0f + 6.0f * (float)(w / 5);
        len += snprintf(text + len, sizeof(text) - len, "O %.3f %.3f 0.000\nH %.3f %.3f 0.000\nH %.3f %.3f 0.000\n", x, y, x + 0.957f, y, x - 0.240f, y + 0.927f);
    }
    md_system_t sys;
    md_system_state_t st;
    ASSERT_TRUE(load_xyz(&sys, &st, text, NULL, alloc));

    md_hbond_params_t p = md_hbond_params_preset(MD_HBOND_PRESET_REALISTIC);
    p.max_ha = 2.6f;
    p.min_dha = 110.0f;
    md_hbond_desc_t desc = { .params = &p };
    md_hbond_set_t full;
    ASSERT_TRUE(md_hbond_compute(&full, &desc, &sys, &st, alloc));
    ASSERT_EQ(1, (int)full.count);
    EXPECT_EQ(3u, full.acceptor[0]);

    md_bitfield_t sel = md_bitfield_create(alloc);
    md_bitfield_set_range(&sel, 0, 3);
    md_bitfield_set_range(&sel, 6, 9);
    desc.set_a = &sel;
    md_hbond_set_t set;
    ASSERT_TRUE(md_hbond_compute(&set, &desc, &sys, &st, alloc));
    EXPECT_EQ(0, (int)set.count);
    md_vm_arena_destroy(alloc);
}

// A chain of two conflicts. The selection's hydrogen H1 prefers acceptor A2 (outside) over A1 (inside), but A2 takes one
// bond only and a third hydrogen H3, further than one search radius from the selection, bonds to it more strongly.
// So in the full system H1 bonds to A1. Getting that right for the selection needs H3, a competitor's competitor.
UTEST(hbond, selection_competitor_chain) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    char text[8192];
    int len = 0;
    const int num_far = 20;
    len += snprintf(text + len, sizeof(text) - len, "%d\ncompetitor of a competitor\n", 12 + 3 * num_far);
    len += snprintf(text + len, sizeof(text) - len,
        "O   0.000  0.000  0.000\n"     // 0: W1, H1 = 1 in the selection
        "H   0.957  0.000  0.000\n"
        "H  -0.240  0.927  0.000\n"
        "O   2.750  0.900  0.000\n"     // 3: A2, outside
        "H   3.707  0.900  0.000\n"
        "H   2.990  0.900  0.927\n"
        "O   2.600 -1.500  0.000\n"     // 6: A1, in the selection
        "H   2.600 -2.457  0.000\n"
        "H   2.600 -1.740  0.927\n"
        "O   2.750  3.700  0.000\n"     // 9: W4, H3 = 10 points straight at A2
        "H   2.750  2.743  0.000\n"
        "H   3.677  3.940  0.000\n");
    for (int w = 0; w < num_far; ++w) {
        const float x = 30.0f + 6.0f * (float)(w % 5), y = 30.0f + 6.0f * (float)(w / 5);
        len += snprintf(text + len, sizeof(text) - len, "O %.3f %.3f 0.000\nH %.3f %.3f 0.000\nH %.3f %.3f 0.000\n", x, y, x + 0.957f, y, x - 0.240f, y + 0.927f);
    }
    md_system_t sys;
    md_system_state_t st;
    ASSERT_TRUE(load_xyz(&sys, &st, text, NULL, alloc));

    md_hbond_params_t p = md_hbond_params_preset(MD_HBOND_PRESET_REALISTIC);
    p.max_ha = 2.6f;
    p.min_dha = 110.0f;
    p.min_xah = 0.0f;
    p.acc_capacity_mode = MD_HBOND_CAPACITY_FIXED;
    p.acc_capacity_fixed = 1;
    md_hbond_desc_t desc = { .params = &p };
    md_hbond_set_t full;
    ASSERT_TRUE(md_hbond_compute(&full, &desc, &sys, &st, alloc));
    ASSERT_EQ(2, (int)full.count);
    EXPECT_EQ(1u, full.hydrogen[0]);
    EXPECT_EQ(6u, full.acceptor[0]);
    EXPECT_EQ(10u, full.hydrogen[1]);
    EXPECT_EQ(3u, full.acceptor[1]);

    md_bitfield_t sel = md_bitfield_create(alloc);
    md_bitfield_set_range(&sel, 0, 3);
    md_bitfield_set_range(&sel, 6, 9);
    desc.set_a = &sel;
    md_hbond_set_t set;
    ASSERT_TRUE(md_hbond_compute(&set, &desc, &sys, &st, alloc));
    ASSERT_EQ(1, (int)set.count);
    EXPECT_EQ(1u, set.hydrogen[0]);
    EXPECT_EQ(6u, set.acceptor[0]);
    md_vm_arena_destroy(alloc);
}

// Aggressive competition (loose geometry, one bond per hydrogen and per acceptor) creates many conflicts; random
// selections must still report exactly the bonds of the full system.
UTEST_F(hbond, selection_heavy_competition) {
    md_allocator_i* alloc = utest_fixture->alloc;
    const md_system_t* sys = &utest_fixture->npt;
    const md_system_state_t* st = &utest_fixture->npt_state;
    md_hbond_params_t p = md_hbond_params_preset(MD_HBOND_PRESET_REALISTIC);
    p.max_ha = 3.0f;
    p.min_dha = 100.0f;
    p.min_xah = 0.0f;
    p.acc_capacity_mode = MD_HBOND_CAPACITY_FIXED;
    p.acc_capacity_fixed = 1;
    md_hbond_desc_t desc = { .params = &p };
    md_hbond_set_t full;
    ASSERT_TRUE(md_hbond_compute(&full, &desc, sys, st, alloc));

    uint32_t seed = 12345;
    for (int round = 0; round < 8; ++round) {
        // About 5% of the residues
        md_bitfield_t sel = md_bitfield_create(alloc);
        for (size_t c = 0; c < sys->component.count; ++c) {
            seed = seed * 1664525u + 1013904223u;
            if ((seed >> 24) < 13) {
                const md_urange_t r = md_system_component_atom_range(sys, c);
                md_bitfield_set_range(&sel, r.beg, r.end);
            }
        }
        desc.set_a = &sel;
        md_hbond_set_t set;
        ASSERT_TRUE(md_hbond_compute(&set, &desc, sys, st, alloc));
        size_t expected = 0, m = 0;
        bool same = true;
        for (size_t k = 0; k < full.count; ++k) {
            const bool d_in = md_bitfield_test_bit(&sel, full.donor[k]) || md_bitfield_test_bit(&sel, full.hydrogen[k]);
            if (!(d_in && md_bitfield_test_bit(&sel, full.acceptor[k]))) continue;
            expected += 1;
            if (m < set.count && set.hydrogen[m] == full.hydrogen[k] && set.acceptor[m] == full.acceptor[k]) m += 1;
            else same = false;
        }
        EXPECT_EQ(expected, set.count);
        EXPECT_TRUE(same);
    }
}

// ### DETERMINISM ###

UTEST_F(hbond, deterministic) {
    md_allocator_i* alloc = utest_fixture->alloc;
    md_hbond_query_t q;
    ASSERT_TRUE(md_hbond_query_init(&q, NULL, &utest_fixture->npt, alloc));
    md_hbond_set_t a, b;
    ASSERT_TRUE(md_hbond_query_eval(&a, &q, &utest_fixture->npt_state, alloc));
    ASSERT_TRUE(md_hbond_query_eval(&b, &q, &utest_fixture->npt_state, alloc));
    ASSERT_EQ(a.count, b.count);
    EXPECT_TRUE(set_is_sorted_by_key(&a));
    EXPECT_EQ(0, memcmp(a.hydrogen, b.hydrogen, sizeof(uint32_t) * a.count));
    EXPECT_EQ(0, memcmp(a.acceptor, b.acceptor, sizeof(uint32_t) * a.count));
    EXPECT_EQ(0, memcmp(a.strength, b.strength, sizeof(float) * a.count));
    md_hbond_set_free(&a);
    md_hbond_set_free(&b);
    md_hbond_query_free(&q);
}

UTEST_F(hbond, state_mismatch) {
    md_hbond_query_t q;
    ASSERT_TRUE(md_hbond_query_init(&q, NULL, &utest_fixture->npt, utest_fixture->alloc));
    md_hbond_set_t set;
    EXPECT_FALSE(md_hbond_query_eval(&set, &q, &utest_fixture->centered_state, utest_fixture->alloc));
}
