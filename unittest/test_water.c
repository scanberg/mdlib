#include "utest.h"

// Water models with virtual sites: the massless M site of TIP4P (and its variants, and OPC) and the two lone pairs of
// TIP5P. Atomistic models of one molecule each, under the names the common force fields and programs give them.

#include <md_system.h>
#include <md_gro.h>
#include <md_pdb.h>
#include <md_tpr.h>
#include <md_util.h>
#include <md_hbond.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>
#include <core/md_vec_math.h>

#include <math.h>
#include <string.h>

// Two waters each: the first donates to the second, which points its own hydrogens away (and, for TIP5P, a lone pair
// at the incoming hydrogen)
// TIP4P, GROMACS (SOL: OW HW1 HW2 MW)
static const char tip4p_gmx_gro[] =
    "TIP4P gromacs\n"
    "8\n"
    "    1SOL     OW    1   1.000   1.000   1.000\n"
    "    1SOL    HW1    2   1.076   1.059   1.000\n"
    "    1SOL    HW2    3   0.924   1.059   1.000\n"
    "    1SOL     MW    4   1.000   1.015   1.000\n"
    "    2SOL     OW    5   1.225   1.174   1.000\n"
    "    2SOL    HW1    6   1.225   1.233   1.076\n"
    "    2SOL    HW2    7   1.225   1.233   0.924\n"
    "    2SOL     MW    8   1.225   1.189   1.000\n"
    "   3.00000   3.00000   3.00000\n";
// TIP5P, GROMACS (SOL: OW HW1 HW2 LP1 LP2)
static const char tip5p_gmx_gro[] =
    "TIP5P gromacs\n"
    "10\n"
    "    1SOL     OW    1   1.000   1.000   1.000\n"
    "    1SOL    HW1    2   1.076   1.059   1.000\n"
    "    1SOL    HW2    3   0.924   1.059   1.000\n"
    "    1SOL    LP1    4   1.000   0.960   1.057\n"
    "    1SOL    LP2    5   1.000   0.960   0.943\n"
    "    2SOL     OW    6   1.225   1.174   1.000\n"
    "    2SOL    HW1    7   1.225   1.233   1.076\n"
    "    2SOL    HW2    8   1.225   1.233   0.924\n"
    "    2SOL    LP1    9   1.283   1.134   1.000\n"
    "    2SOL    LP2   10   1.168   1.134   1.000\n"
    "   3.00000   3.00000   3.00000\n";
// TIP4P-Ew / OPC, AMBER (WAT: O H1 H2 EPW)
static const char tip4p_amber_pdb[] =
    "HETATM    1  O   WAT A   1      10.000  10.000  10.000  1.00  0.00           O\n"
    "HETATM    2  H1  WAT A   1      10.757  10.586  10.000  1.00  0.00           H\n"
    "HETATM    3  H2  WAT A   1       9.243  10.586  10.000  1.00  0.00           H\n"
    "HETATM    4  EPW WAT A   1      10.000  10.150  10.000  1.00  0.00            \n"
    "HETATM    5  O   WAT A   2      12.254  11.744  10.000  1.00  0.00           O\n"
    "HETATM    6  H1  WAT A   2      12.254  12.330  10.757  1.00  0.00           H\n"
    "HETATM    7  H2  WAT A   2      12.254  12.330   9.243  1.00  0.00           H\n"
    "HETATM    8  EPW WAT A   2      12.254  11.894  10.000  1.00  0.00            \n"
    "END\n";
// TIP4P, CHARMM (TIP4: OH2 OM H1 H2), where OM must not become an oxygen
static const char tip4p_charmm_pdb[] =
    "HETATM    1  OH2 TIP4A   1      10.000  10.000  10.000  1.00  0.00           O\n"
    "HETATM    2  OM  TIP4A   1      10.000  10.150  10.000  1.00  0.00            \n"
    "HETATM    3  H1  TIP4A   1      10.757  10.586  10.000  1.00  0.00           H\n"
    "HETATM    4  H2  TIP4A   1       9.243  10.586  10.000  1.00  0.00           H\n"
    "HETATM    5  OH2 TIP4A   2      12.254  11.744  10.000  1.00  0.00           O\n"
    "HETATM    6  OM  TIP4A   2      12.254  11.894  10.000  1.00  0.00            \n"
    "HETATM    7  H1  TIP4A   2      12.254  12.330  10.757  1.00  0.00           H\n"
    "HETATM    8  H2  TIP4A   2      12.254  12.330   9.243  1.00  0.00           H\n"
    "END\n";
// TIP4P-Ew, OpenMM (HOH: O H1 H2 M)
static const char tip4p_openmm_pdb[] =
    "HETATM    1  O   HOH A   1      10.000  10.000  10.000  1.00  0.00           O\n"
    "HETATM    2  H1  HOH A   1      10.757  10.586  10.000  1.00  0.00           H\n"
    "HETATM    3  H2  HOH A   1       9.243  10.586  10.000  1.00  0.00           H\n"
    "HETATM    4  M   HOH A   1      10.000  10.150  10.000  1.00  0.00            \n"
    "HETATM    5  O   HOH A   2      12.254  11.744  10.000  1.00  0.00           O\n"
    "HETATM    6  H1  HOH A   2      12.254  12.330  10.757  1.00  0.00           H\n"
    "HETATM    7  H2  HOH A   2      12.254  12.330   9.243  1.00  0.00           H\n"
    "HETATM    8  M   HOH A   2      12.254  11.894  10.000  1.00  0.00            \n"
    "END\n";
// TIP5P, AMBER (WAT: O H1 H2 EP1 EP2)
static const char tip5p_amber_pdb[] =
    "HETATM    1  O   WAT A   1      10.000  10.000  10.000  1.00  0.00           O\n"
    "HETATM    2  H1  WAT A   1      10.757  10.586  10.000  1.00  0.00           H\n"
    "HETATM    3  H2  WAT A   1       9.243  10.586  10.000  1.00  0.00           H\n"
    "HETATM    4  EP1 WAT A   1      10.000   9.596  10.572  1.00  0.00            \n"
    "HETATM    5  EP2 WAT A   1      10.000   9.596   9.428  1.00  0.00            \n"
    "HETATM    6  O   WAT A   2      12.254  11.744  10.000  1.00  0.00           O\n"
    "HETATM    7  H1  WAT A   2      12.254  12.330  10.757  1.00  0.00           H\n"
    "HETATM    8  H2  WAT A   2      12.254  12.330   9.243  1.00  0.00           H\n"
    "HETATM    9  EP1 WAT A   2      12.826  11.340  10.000  1.00  0.00            \n"
    "HETATM   10  EP2 WAT A   2      11.682  11.340  10.000  1.00  0.00            \n"
    "END\n";
// TIP5P, CHARMM (TIP5: OH2 H1 H2 LP1 LP2)
static const char tip5p_charmm_pdb[] =
    "HETATM    1  OH2 TIP5A   1      10.000  10.000  10.000  1.00  0.00           O\n"
    "HETATM    2  H1  TIP5A   1      10.757  10.586  10.000  1.00  0.00           H\n"
    "HETATM    3  H2  TIP5A   1       9.243  10.586  10.000  1.00  0.00           H\n"
    "HETATM    4  LP1 TIP5A   1      10.000   9.596  10.572  1.00  0.00            \n"
    "HETATM    5  LP2 TIP5A   1      10.000   9.596   9.428  1.00  0.00            \n"
    "HETATM    6  OH2 TIP5A   2      12.254  11.744  10.000  1.00  0.00           O\n"
    "HETATM    7  H1  TIP5A   2      12.254  12.330  10.757  1.00  0.00           H\n"
    "HETATM    8  H2  TIP5A   2      12.254  12.330   9.243  1.00  0.00           H\n"
    "HETATM    9  LP1 TIP5A   2      12.826  11.340  10.000  1.00  0.00            \n"
    "HETATM   10  LP2 TIP5A   2      11.682  11.340  10.000  1.00  0.00            \n"
    "END\n";
typedef struct water_case_t {
    const char* name;
    const char* text;
    bool gro;
    int num_sites;      // Virtual sites per molecule
} water_case_t;

static const water_case_t water_cases[] = {
    { "tip4p_gmx_gro", tip4p_gmx_gro, true, 1 },
    { "tip5p_gmx_gro", tip5p_gmx_gro, true, 2 },
    { "tip4p_amber_pdb", tip4p_amber_pdb, false, 1 },
    { "tip4p_charmm_pdb", tip4p_charmm_pdb, false, 1 },
    { "tip4p_openmm_pdb", tip4p_openmm_pdb, false, 1 },
    { "tip5p_amber_pdb", tip5p_amber_pdb, false, 2 },
    { "tip5p_charmm_pdb", tip5p_charmm_pdb, false, 2 },
};

static bool load_case(md_system_t* sys, md_system_state_t* st, const water_case_t* c, md_allocator_i* alloc) {
    *sys = (md_system_t){ .alloc = alloc };
    *st  = (md_system_state_t){ .alloc = alloc };
    const str_t text = str_from_cstr(c->text);
    const bool ok = c->gro ? md_gro_system_init_from_str(sys, st, text) : md_pdb_system_init_from_str(sys, st, text, MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE);
    return ok && md_util_system_infer(sys, st, MD_UTIL_INFER_ALL);
}

static inline md_particle_kind_t particle(const md_system_t* sys, size_t i) {
    return md_atom_particle_kind(&sys->atom, i);
}

UTEST(water, virtual_sites) {
    for (size_t k = 0; k < ARRAY_SIZE(water_cases); ++k) {
        const water_case_t* c = &water_cases[k];
        md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
        md_system_t sys;
        md_system_state_t st;
        ASSERT_TRUE(load_case(&sys, &st, c, alloc));

        // Two water molecules, not coarse grained, with only their O-H bonds
        ASSERT_EQ(2, (int)sys.component.count);
        EXPECT_EQ(2 * (3 + c->num_sites), (int)sys.atom.count);
        for (size_t ci = 0; ci < sys.component.count; ++ci) {
            EXPECT_EQ(MD_COMPONENT_KIND_WATER, md_component_kind(&sys.component, ci));
        }
        EXPECT_FALSE(md_system_is_coarse_grained(&sys));
        // Each water molecule is an instance of its own, of one water entity
        ASSERT_EQ(2, (int)sys.instance.count);
        ASSERT_EQ(1, (int)sys.entity.count);
        EXPECT_EQ(MD_ENTITY_KIND_WATER, md_system_instance_entity_kind(&sys, 0));
        EXPECT_EQ(MD_ENTITY_KIND_WATER, md_system_instance_entity_kind(&sys, 1));
        EXPECT_EQ(4, (int)sys.bond.count);

        int num_o = 0, num_h = 0, num_sites = 0;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            const md_atomic_number_t z = md_atom_atomic_number(&sys.atom, i);
            const size_t deg = md_bond_conn_count(&sys.bond, i);
            EXPECT_EQ(MD_COMPONENT_KIND_WATER, md_system_atom_component_kind(&sys, i));
            if (particle(&sys, i) == MD_PARTICLE_VIRTUAL_SITE) {
                EXPECT_EQ(0, (int)z);
                EXPECT_EQ(0.0f, md_atom_mass(&sys.atom, i));
                EXPECT_EQ(0.0f, md_atom_radius(&sys.atom, i));
                EXPECT_EQ(0, (int)deg);
                num_sites += 1;
            } else if (z == MD_Z_O) {
                EXPECT_EQ(2, (int)deg);
                num_o += 1;
            } else {
                EXPECT_EQ(MD_Z_H, (int)z);
                EXPECT_EQ(1, (int)deg);
                num_h += 1;
            }
        }
        EXPECT_EQ(2, num_o);
        EXPECT_EQ(4, num_h);
        EXPECT_EQ(2 * c->num_sites, num_sites);

        // Each molecule is one structure, its sites included
        EXPECT_EQ(2, (int)sys.structure.count);

        // The sites change nothing for hydrogen bonds: one, from the first water to the second
        md_hbond_set_t hb;
        ASSERT_TRUE(md_hbond_compute(&hb, NULL, &sys, &st, alloc));
        ASSERT_EQ(1, (int)hb.count);
        const md_urange_t r0 = md_component_atom_range(&sys.component, 0);
        const md_urange_t r1 = md_component_atom_range(&sys.component, 1);
        EXPECT_TRUE(hb.donor[0] >= r0.beg && hb.donor[0] < r0.end && md_atom_atomic_number(&sys.atom, hb.donor[0]) == MD_Z_O);
        EXPECT_TRUE(hb.acceptor[0] >= r1.beg && hb.acceptor[0] < r1.end && md_atom_atomic_number(&sys.atom, hb.acceptor[0]) == MD_Z_O);

        md_vm_arena_destroy(alloc);
    }
}

// A TIP4P water split across the periodic boundary: unwrapping brings the M site back to its oxygen
UTEST(water, unwrap_virtual_site) {
    const char* text =
        "TIP4P across the boundary at x = 3 nm\n"
        "4\n"
        "    1SOL     OW    1   2.990   1.000   1.000\n"
        "    1SOL    HW1    2   0.049   1.076   1.000\n"
        "    1SOL    HW2    3   0.049   0.924   1.000\n"
        "    1SOL     MW    4   0.005   1.000   1.000\n"
        "   3.00000   3.00000   3.00000\n";
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &st, str_from_cstr(text)));
    ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
    EXPECT_EQ(2, (int)sys.bond.count);
    EXPECT_EQ(1, (int)sys.structure.count);
    md_util_unwrap_system(&st, &sys);
    for (int i = 1; i < 4; ++i) {
        const float d = vec3_length(vec3_sub(st.xyz[i], st.xyz[0]));
        EXPECT_LT(d, 1.0f);
    }
    md_vm_arena_destroy(alloc);
}

// A peptide in TIP4P water, from the GROMACS coordinates and from the run input: the same waters and sites either way
UTEST(water, tip4p_gro_and_tpr) {
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    size_t counts[2][4] = { 0 };
    for (int f = 0; f < 2; ++f) {
        md_system_t sys = { .alloc = alloc };
        md_system_state_t st = { .alloc = alloc };
        const bool ok = f == 0 ? md_gro_system_init_from_file(&sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/tpr/peptide_tip4p.gro"))
                               : md_tpr_system_init_from_file(&sys, &st, STR_LIT(MD_UNITTEST_DATA_DIR "/tpr/peptide_tip4p.tpr"));
        ASSERT_TRUE(ok);
        ASSERT_TRUE(md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL));
        size_t water = 0, sites = 0, bonded_sites = 0, cg = 0;
        for (size_t ci = 0; ci < sys.component.count; ++ci) water += md_component_kind(&sys.component, ci) == MD_COMPONENT_KIND_WATER;
        for (size_t i = 0; i < sys.atom.count; ++i) {
            if (particle(&sys, i) == MD_PARTICLE_VIRTUAL_SITE) {
                sites += 1;
                bonded_sites += md_bond_conn_count(&sys.bond, i) != 0;
            }
        }
        for (size_t t = 0; t < sys.atom.type.count; ++t) cg += md_atom_type_particle_kind(&sys.atom.type, t) == MD_PARTICLE_BEAD;
        counts[f][0] = water;
        counts[f][1] = sites;
        counts[f][2] = bonded_sites;
        counts[f][3] = cg;
    }
    EXPECT_EQ((size_t)506, counts[0][0]);
    EXPECT_EQ((size_t)506, counts[0][1]);
    EXPECT_EQ((size_t)0, counts[0][2]);
    EXPECT_EQ((size_t)0, counts[0][3]);
    for (int k = 0; k < 4; ++k) EXPECT_EQ(counts[0][k], counts[1][k]);
    md_vm_arena_destroy(alloc);
}

// The Martini water bead is still coarse grained: the virtual site entry must not catch it
UTEST(water, martini_bead_unchanged) {
    const char* text =
        "Martini water\n"
        "2\n"
        "    1W        W    1   1.000   1.000   1.000\n"
        "    2W        W    2   1.500   1.000   1.000\n"
        "   3.00000   3.00000   3.00000\n";
    md_allocator_i* alloc = md_vm_arena_create(GIGABYTES(1));
    md_system_t sys = { .alloc = alloc };
    md_system_state_t st = { .alloc = alloc };
    ASSERT_TRUE(md_gro_system_init_from_str(&sys, &st, str_from_cstr(text)));
    EXPECT_EQ(MD_PARTICLE_BEAD, particle(&sys, 0));
    EXPECT_EQ(MD_COMPONENT_KIND_WATER, md_component_kind(&sys.component, 0));
    md_vm_arena_destroy(alloc);
}
