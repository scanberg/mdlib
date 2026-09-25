#include "utest.h"

#include <md_nonbonded.h>
#include <md_tpr.h>
#include <md_system.h>
#include <md_util.h>
#include <core/md_vec_math.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

#include <math.h>
#include <string.h>

#define TPR_DIR MD_UNITTEST_DATA_DIR "/tpr/"

// Relative agreement with a tolerance for values near zero
static bool close_to(double value, double expected, double rel, double abs_tol) {
    return fabs(value - expected) <= rel * fabs(expected) + abs_tol;
}

// Two B particles (sigma 0.4 nm, epsilon 1 kJ/mol) at a distance: LJ (SR) of gmx mdrun -rerun, GROMACS 2023.3,
// vdw cut-off at 1.2 nm, force switch from 0.8 nm and potential shift
UTEST(nonbonded, lj_pair_against_gromacs) {
    const double c6  = 4.0 * pow(0.4, 6);
    const double c12 = 4.0 * pow(0.4, 12);
    const double r[3]         = { 0.5, 0.9, 1.1 };
    const double fswitch[3]   = { -7.60624e-01, -1.77749e-02, -6.05566e-04 };
    const double potshift[3]  = { -7.68219e-01, -2.51123e-02, -3.74752e-03 };

    md_nb_potential_t fsw, psh;
    ASSERT_TRUE(md_nb_potential_init(&fsw, &(md_nb_desc_t){ .lj_modifier = MD_NB_MODIFIER_FORCE_SWITCH, .lj_cutoff = 1.2, .lj_switch = 0.8 }));
    ASSERT_TRUE(md_nb_potential_init(&psh, &(md_nb_desc_t){ .lj_modifier = MD_NB_MODIFIER_POT_SHIFT, .lj_cutoff = 1.2 }));
    for (int i = 0; i < 3; ++i) {
        EXPECT_TRUE(close_to(md_nb_lj_energy(&fsw, c6, c12, r[i] * r[i]), fswitch[i], 1e-5, 1e-9));
        EXPECT_TRUE(close_to(md_nb_lj_energy(&psh, c6, c12, r[i] * r[i]), potshift[i], 1e-5, 1e-9));
    }
    // Before the switch the force switch is the plain potential, shifted
    md_nb_potential_t plain;
    ASSERT_TRUE(md_nb_potential_init(&plain, &(md_nb_desc_t){ .lj_modifier = MD_NB_MODIFIER_NONE, .lj_cutoff = 1.2 }));
    const double d0 = md_nb_lj_energy(&fsw, c6, c12, 0.25) - md_nb_lj_energy(&plain, c6, c12, 0.25);
    const double d1 = md_nb_lj_energy(&fsw, c6, c12, 0.49) - md_nb_lj_energy(&plain, c6, c12, 0.49);
    EXPECT_NEAR(d0, d1, 1e-12);
    // Nothing from the cut-off on
    EXPECT_EQ(0.0, md_nb_lj_energy(&fsw, c6, c12, 1.44));
    EXPECT_EQ(0.0, md_nb_lj_energy(&plain, c6, c12, 2.0));
}

// Shifted and switched potentials reach zero at the cut-off, the switches with a continuous force
UTEST(nonbonded, continuity) {
    const double c6 = 1.0e-2, c12 = 1.0e-5, rc = 1.1, rsw = 0.85;
    const md_nb_modifier_t mods[3] = { MD_NB_MODIFIER_POT_SHIFT, MD_NB_MODIFIER_POT_SWITCH, MD_NB_MODIFIER_FORCE_SWITCH };
    for (int m = 0; m < 3; ++m) {
        md_nb_potential_t pot;
        ASSERT_TRUE(md_nb_potential_init(&pot, &(md_nb_desc_t){ .lj_modifier = mods[m], .lj_cutoff = rc, .lj_switch = rsw }));
        const double e = 1.0e-9;
        EXPECT_NEAR(0.0, md_nb_lj_energy(&pot, c6, c12, (rc - e) * (rc - e)), 1e-9);
        if (mods[m] != MD_NB_MODIFIER_POT_SHIFT) {
            // The force, as one sided differences: zero at the cut-off, and the same on both sides of the switch
            const double h = 1.0e-6;
            #define V(r) md_nb_lj_energy(&pot, c6, c12, (r) * (r))
            EXPECT_NEAR(0.0, -(V(rc - e) - V(rc - e - h)) / h, 1e-5);
            EXPECT_NEAR(-(V(rsw) - V(rsw - h)) / h, -(V(rsw + h) - V(rsw)) / h, 1e-4);
            #undef V
        }
    }

    // Reaction field and shifted Coulomb vanish at the cut-off
    md_nb_potential_t rf, cut, ew;
    ASSERT_TRUE(md_nb_potential_init(&rf, &(md_nb_desc_t){ .lj_cutoff = rc, .coulomb = MD_NB_COULOMB_REACTION_FIELD, .coulomb_cutoff = rc, .epsilon_r = 15, .epsilon_rf = 0 }));
    ASSERT_TRUE(md_nb_potential_init(&cut, &(md_nb_desc_t){ .lj_cutoff = rc, .coulomb = MD_NB_COULOMB_CUTOFF, .coulomb_modifier = MD_NB_MODIFIER_POT_SHIFT, .coulomb_cutoff = rc, .epsilon_r = 1 }));
    ASSERT_TRUE(md_nb_potential_init(&ew, &(md_nb_desc_t){ .lj_cutoff = rc, .coulomb = MD_NB_COULOMB_EWALD, .coulomb_modifier = MD_NB_MODIFIER_POT_SHIFT, .coulomb_cutoff = rc, .epsilon_r = 1, .ewald_rtol = 1e-5 }));
    const double r2 = (rc - 1e-9) * (rc - 1e-9);
    EXPECT_NEAR(0.0, md_nb_coulomb_energy(&rf, 1.0, r2), 1e-6);
    EXPECT_NEAR(0.0, md_nb_coulomb_energy(&cut, 1.0, r2), 1e-6);
    EXPECT_NEAR(0.0, md_nb_coulomb_energy(&ew, 1.0, r2), 1e-6);
    EXPECT_NEAR(1e-5, erfc(ew.ewald_beta * rc), 1e-12);
    // Infinite dielectric: no electrostatics
    md_nb_potential_t none;
    ASSERT_TRUE(md_nb_potential_init(&none, &(md_nb_desc_t){ .lj_cutoff = rc, .coulomb = MD_NB_COULOMB_REACTION_FIELD, .coulomb_cutoff = rc, .epsilon_r = 0 }));
    EXPECT_EQ(0.0, md_nb_coulomb_energy(&none, 1.0, 0.25));
}

// The test systems of test_data/tpr: every non-excluded pair within the cut-off, summed per energy group (the
// three chains, CHN, and the four single particles, SOL), against the group energies of gmx mdrun -rerun.
// Lennard-Jones everywhere. Coulomb between the groups; within a group GROMACS adds the reaction field and
// Ewald corrections of excluded pairs and the self terms, which belong to no pair, except for a plain cut-off.
typedef struct nb_system_case_t {
    const char* file;
    double coul[3];     // CHN-CHN, CHN-SOL, SOL-SOL
    double lj[3];
    bool coul_within_groups;
    double coul_rel;    // Ewald is tabulated in GROMACS, to about 1e-5
} nb_system_case_t;

UTEST(nonbonded, systems_against_gromacs) {
    const nb_system_case_t cases[] = {
        { "nb_fswitch_rf.tpr",   { -0.367957115173, 4.140795707703, -15.993323326111 }, { -0.987358510494, -0.170260623097, -0.001748379553 }, false, 1e-5 },
        { "nb_pswitch_cut.tpr",  { 12.926515579224, 29.370998382568, -13.620864868164 }, { -1.002470016479, -0.167487099767, -0.003118790453 }, true,  1e-5 },
        { "nb_potshift_pme.tpr", { -17.877027511597, 0.242650389671, -80.129722595215 }, { -0.999014973640, -0.182620733976, -0.002159389667 }, false, 1e-4 },
    };
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    for (size_t c = 0; c < sizeof(cases) / sizeof(cases[0]); ++c) {
        char path[512];
        snprintf(path, sizeof(path), "%s%s", TPR_DIR, cases[c].file);
        md_tpr_data_t tpr = {0};
        ASSERT_TRUE(md_tpr_data_parse_file(&tpr, str_from_cstr(path), arena));
        md_nb_potential_t pot;
        ASSERT_TRUE(md_nb_potential_init_from_tpr(&pot, &tpr));
        ASSERT_EQ(19u, tpr.num_atoms);
        ASSERT_TRUE(tpr.box[0][1] == 0 && tpr.box[0][2] == 0 && tpr.box[1][2] == 0);

        // Per particle type and charge, from the molecule blocks
        uint16_t type[19];
        float charge[19];
        size_t n = 0;
        for (size_t b = 0; b < tpr.num_molblocks; ++b) {
            const md_tpr_moltype_t* mt = &tpr.moltypes[tpr.molblocks[b].moltype];
            for (int32_t m = 0; m < tpr.molblocks[b].nmol; ++m) {
                for (size_t a = 0; a < mt->num_atoms; ++a, ++n) {
                    type[n] = mt->atoms[a].type_idx;
                    charge[n] = mt->atoms[a].charge;
                }
            }
        }
        ASSERT_EQ(19u, n);

        double lj[3] = {0}, coul[3] = {0};
        for (size_t i = 0; i < 19; ++i) {
            for (size_t j = i + 1; j < 19; ++j) {
                if (md_tpr_atoms_excluded(&tpr, i, j)) continue;
                double r2 = 0;
                for (int k = 0; k < 3; ++k) {
                    const double L = tpr.box[k][k];
                    double d = (double)tpr.x[i * 3 + k] - (double)tpr.x[j * 3 + k];
                    d -= L * round(d / L);
                    r2 += d * d;
                }
                const int g = (i >= 15) + (j >= 15);
                const md_tpr_lj_t p = md_tpr_lj_pair(&tpr, type[i], type[j]);
                lj[g]   += md_nb_lj_energy(&pot, p.c6, p.c12, r2);
                coul[g] += md_nb_coulomb_energy(&pot, (double)charge[i] * charge[j], r2);
            }
        }
        for (int g = 0; g < 3; ++g) {
            EXPECT_TRUE(close_to(lj[g], cases[c].lj[g], 1e-5, 1e-7));
            if (g == 1 || cases[c].coul_within_groups) {
                EXPECT_TRUE(close_to(coul[g], cases[c].coul[g], cases[c].coul_rel, 1e-7));
            }
            if (!close_to(lj[g], cases[c].lj[g], 1e-5, 1e-7) || ((g == 1 || cases[c].coul_within_groups) && !close_to(coul[g], cases[c].coul[g], cases[c].coul_rel, 1e-7))) {
                printf("  %s group pair %d: LJ %.9g (GROMACS %.9g), Coulomb %.9g (GROMACS %.9g)\n", cases[c].file, g, lj[g], cases[c].lj[g], coul[g], cases[c].coul[g]);
            }
        }
        md_tpr_data_free(&tpr, arena);
    }
    md_arena_allocator_destroy(arena);
}

UTEST(nonbonded, from_tpr) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));

    // LJ-PME has a grid part which is not a sum over pairs
    md_tpr_data_t tpr = {0};
    ASSERT_TRUE(md_tpr_data_parse_file(&tpr, STR_LIT(TPR_DIR "nb_ljpme.tpr"), arena));
    md_nb_potential_t pot;
    EXPECT_FALSE(md_nb_potential_init_from_tpr(&pot, &tpr));
    md_tpr_data_free(&tpr, arena);

    // Martini: reaction field with eps_r 15 and an infinite eps_rf, cut off at 1.1 nm, LJ potential shift at 1.1 nm
    ASSERT_TRUE(md_tpr_data_parse_file(&tpr, STR_LIT(TPR_DIR "martini3.tpr"), arena));
    ASSERT_TRUE(md_nb_potential_init_from_tpr(&pot, &tpr));
    EXPECT_EQ((int)MD_NB_COULOMB_REACTION_FIELD, pot.coulomb);
    EXPECT_EQ((int)MD_NB_MODIFIER_POT_SHIFT, pot.lj_modifier);
    const double rc = 1.1;
    EXPECT_NEAR(1.0 / (2 * rc * rc * rc), pot.k_rf, 1e-6);
    EXPECT_NEAR(1.5 / rc, pot.c_rf, 1e-6);
    EXPECT_NEAR(MD_NB_ONE_4PI_EPS0 / 15.0, pot.epsfac, 1e-6);
    EXPECT_NEAR(1.1, md_nb_potential_cutoff(&pot), 1e-6);
    md_tpr_data_free(&tpr, arena);

    // Without readable simulation parameters there is nothing to go by
    tpr = (md_tpr_data_t){0};
    EXPECT_FALSE(md_nb_potential_init_from_tpr(&pot, &tpr));

    md_arena_allocator_destroy(arena);
}

// The force field a system loaded from a tpr carries: the same energies, from the system's own coordinates (Å)
UTEST(nonbonded, system_forcefield) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_system_t sys = { .alloc = arena };
    md_system_state_t state = { .alloc = arena };
    ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &state, STR_LIT(TPR_DIR "nb_fswitch_rf.tpr")));
    ASSERT_TRUE(sys.nonbonded != NULL);
    const md_nb_forcefield_t* ff = sys.nonbonded;
    ASSERT_EQ((size_t)19, ff->num_atoms);
    EXPECT_TRUE(md_nb_forcefield_excluded(ff, 5, 8));
    EXPECT_FALSE(md_nb_forcefield_excluded(ff, 5, 9));
    EXPECT_FALSE(md_nb_forcefield_excluded(ff, 4, 5));

    const double lj_ref[3]   = { -0.987358510494, -0.170260623097, -0.001748379553 };
    const double coul_ref[3] = { -0.367957115173, 4.140795707703, -15.993323326111 };
    double lj[3] = {0}, coul[3] = {0};
    for (uint32_t i = 0; i < 19; ++i) {
        for (uint32_t j = i + 1; j < 19; ++j) {
            vec3_t d = vec3_sub(state.xyz[i], state.xyz[j]);
            md_util_min_image_vec3(&d, 1, &state.unitcell);
            const double r2 = 0.01 * ((double)d.x * d.x + (double)d.y * d.y + (double)d.z * d.z);
            double elj, ec;
            md_nb_forcefield_pair_energy(ff, i, j, r2, &elj, &ec);
            const int g = (i >= 15) + (j >= 15);
            lj[g] += elj;
            coul[g] += ec;
        }
    }
    for (int g = 0; g < 3; ++g) {
        EXPECT_TRUE(close_to(lj[g], lj_ref[g], 1e-4, 1e-6));
    }
    EXPECT_TRUE(close_to(coul[1], coul_ref[1], 1e-4, 1e-6));
    md_system_free(&sys);

    // LJ-PME cannot be evaluated pair by pair: no force field, so no energies
    md_system_t sys2 = { .alloc = arena };
    md_system_state_t state2 = { .alloc = arena };
    ASSERT_TRUE(md_tpr_system_init_from_file(&sys2, &state2, STR_LIT(TPR_DIR "nb_ljpme.tpr")));
    EXPECT_TRUE(sys2.nonbonded == NULL);

    md_arena_allocator_destroy(arena);
}
