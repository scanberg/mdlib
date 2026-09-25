#include "utest.h"
#include <string.h>
#include <math.h>

#include <md_tpr.h>
#include <md_gro.h>
#include <md_xdr.h>
#include <md_system.h>
#include <md_util.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_os.h>
#include <core/md_str.h>

#define TPR_DIR MD_UNITTEST_DATA_DIR "/tpr/"

// The dipeptide has 46 atoms and 47 bonds, the rest is water and four ions
#define PEPTIDE_ATOMS 46
#define PEPTIDE_BONDS 47

UTEST(xdr, cursor) {
    // 7 (int), -2 (int), 'abcde' as an XDR string, 1.5f, the hyper 0x0000000100000002
    const uint8_t buf[] = {
        0,0,0,7,  0xFF,0xFF,0xFF,0xFE,  0,0,0,5, 'a','b','c','d', 'e',0,0,0,  0x3F,0xC0,0,0,  0,0,0,1, 0,0,0,2,
    };
    md_xdr_t xdr = md_xdr_init(buf, sizeof(buf));
    int32_t i;
    uint32_t u;
    str_t str;
    float f;
    int64_t h;
    EXPECT_TRUE(md_xdr_read_i32(&xdr, &i));
    EXPECT_EQ(7, i);
    EXPECT_TRUE(md_xdr_read_i32(&xdr, &i));
    EXPECT_EQ(-2, i);
    EXPECT_TRUE(md_xdr_read_string(&xdr, &str, 64));
    EXPECT_TRUE(str_eq(str, STR_LIT("abcde")));
    EXPECT_EQ(20u, xdr.pos);   // Padded to a multiple of 4
    EXPECT_TRUE(md_xdr_read_f32(&xdr, &f));
    EXPECT_EQ(1.5f, f);
    EXPECT_TRUE(md_xdr_read_i64(&xdr, &h));
    EXPECT_EQ((int64_t)0x100000002LL, h);
    EXPECT_EQ(0u, md_xdr_remaining(&xdr));

    // Past the end: fails, stays failed, and leaves the value zeroed
    EXPECT_FALSE(md_xdr_read_u32(&xdr, &u));
    EXPECT_EQ(0u, u);
    EXPECT_FALSE(md_xdr_ok(&xdr));

    // A string longer than allowed fails without moving the cursor
    xdr = md_xdr_init(buf + 8, sizeof(buf) - 8);
    EXPECT_FALSE(md_xdr_read_string(&xdr, &str, 4));
    EXPECT_EQ(0u, xdr.pos);
}

// Checks the system built from a tpr against the .gro gmx writes from the same file
static void compare_with_gro(int* utest_result, const md_system_t* sys, const md_system_state_t* state, str_t gro_path) {
    md_allocator_i* alloc = md_get_heap_allocator();
    md_gro_data_t gro = {0};
    ASSERT_TRUE(md_gro_data_parse_file(&gro, gro_path, alloc));
    ASSERT_EQ(gro.num_atoms, sys->atom.count);

    size_t comp = 0;
    size_t res_mismatch = 0, coord_mismatch = 0;
    for (size_t i = 0; i < sys->atom.count; ++i) {
        while (comp + 1 < sys->component.count && sys->component.atom_offset[comp + 1] <= i) comp++;
        const md_gro_atom_t* a = &gro.atom_data[i];
        str_t res_name = str_from_cstrn(a->res_name, sizeof(a->res_name));
        if (!str_eq(md_component_name(&sys->component, comp), res_name) || sys->component.seq_id[comp] != a->res_id) {
            res_mismatch += 1;
        }
        // A .gro has three decimals in nm
        if (fabsf(state->xyz[i].x - a->x * 10.0f) > 0.006f || fabsf(state->xyz[i].y - a->y * 10.0f) > 0.006f || fabsf(state->xyz[i].z - a->z * 10.0f) > 0.006f) {
            coord_mismatch += 1;
        }
    }
    EXPECT_EQ(0u, res_mismatch);
    EXPECT_EQ(0u, coord_mismatch);
    EXPECT_NEAR(gro.box[0][0] * 10.0f, (float)state->unitcell.x, 0.006f);
    EXPECT_NEAR(gro.box[1][1] * 10.0f, (float)state->unitcell.y, 0.006f);
    EXPECT_NEAR(gro.box[2][2] * 10.0f, (float)state->unitcell.z, 0.006f);

    md_gro_data_free(&gro, alloc);
}

UTEST(tpr, parse_tip3p) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));

    md_tpr_data_t data = {0};
    ASSERT_TRUE(md_tpr_data_parse_file(&data, STR_LIT(TPR_DIR "peptide_tip3p.tpr"), arena));
    EXPECT_EQ(129, data.file_version);
    EXPECT_FALSE(data.double_precision);
    EXPECT_TRUE(str_eq(data.name, STR_LIT("Protein in water")));
    EXPECT_EQ(1577u, data.num_atoms);
    EXPECT_EQ(MD_TPR_PBC_XYZ, data.pbc);
    EXPECT_TRUE(data.has_box);
    EXPECT_NEAR(2.5238f, data.box[0][0], 1e-4f);
    EXPECT_TRUE(data.x != NULL);

    ASSERT_EQ(4u, data.num_moltypes);
    ASSERT_EQ(4u, data.num_molblocks);
    const md_tpr_moltype_t* prot = &data.moltypes[0];
    const md_tpr_moltype_t* sol  = &data.moltypes[1];
    EXPECT_TRUE(str_eq(prot->name, STR_LIT("Protein_chain_U")));
    EXPECT_EQ(PEPTIDE_ATOMS, (int)prot->num_atoms);
    EXPECT_EQ(2u, prot->num_residues);
    EXPECT_EQ(PEPTIDE_BONDS, (int)prot->num_bonds);
    EXPECT_TRUE(str_eq(prot->residues[0].name, STR_LIT("TRP")));
    EXPECT_EQ(221, prot->residues[0].nr);
    EXPECT_TRUE(str_eq(prot->atoms[0].name, STR_LIT("N")));
    EXPECT_TRUE(str_eq(prot->atoms[0].type, STR_LIT("N3")));
    EXPECT_EQ(7, prot->atoms[0].atomic_number);
    EXPECT_NEAR(14.01f, prot->atoms[0].mass, 1e-4f);

    // SETTLE water: two O-H bonds and no H-H
    EXPECT_TRUE(str_eq(sol->name, STR_LIT("SOL")));
    ASSERT_EQ(3u, sol->num_atoms);
    ASSERT_EQ(2u, sol->num_bonds);
    EXPECT_EQ(0, sol->bonds[0].idx[0]);
    EXPECT_EQ(1, sol->bonds[0].idx[1]);
    EXPECT_EQ(0, sol->bonds[1].idx[0]);
    EXPECT_EQ(2, sol->bonds[1].idx[1]);
    EXPECT_EQ(509, data.molblocks[1].nmol);

    double total_charge = 0;
    for (size_t b = 0; b < data.num_molblocks; ++b) {
        const md_tpr_moltype_t* mt = &data.moltypes[data.molblocks[b].moltype];
        for (size_t i = 0; i < mt->num_atoms; ++i) total_charge += mt->atoms[i].charge * data.molblocks[b].nmol;
    }
    EXPECT_NEAR(0.0, total_charge, 1e-3);

    md_arena_allocator_destroy(arena);
}

UTEST(tpr, system_tip3p) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_system_t sys = { .alloc = arena };
    md_system_state_t state = { .alloc = arena };

    ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &state, STR_LIT(TPR_DIR "peptide_tip3p.tpr")));
    EXPECT_EQ(1577u, sys.atom.count);
    EXPECT_EQ(1577u, state.num_atoms);
    EXPECT_EQ(2u + 509u + 4u, sys.component.count);
    EXPECT_EQ((size_t)PEPTIDE_BONDS + 509 * 2, sys.bond.count);
    for (size_t i = 0; i < sys.bond.count; ++i) {
        EXPECT_TRUE(sys.bond.flags[i] & MD_BOND_FLAG_TOPOLOGY);
    }

    // Every atom has its element from the topology
    size_t unknown = 0;
    for (size_t i = 0; i < sys.atom.count; ++i) {
        if (md_atom_type_atomic_number(&sys.atom.type, sys.atom.type_idx[i]) == 0) unknown += 1;
    }
    EXPECT_EQ(0u, unknown);

    const md_attribute_t* charge = md_attributes_find(&sys.attributes, STR_LIT("atom/charge"));
    EXPECT_TRUE(charge != NULL);
    // Mass lives on the type, and is exact there: no per atom copy
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("atom/mass")) == NULL);
    EXPECT_NEAR(14.01f, md_atom_mass(&sys.atom, 0), 1e-4f);         // N
    EXPECT_NEAR(16.0f, md_atom_mass(&sys.atom, PEPTIDE_ATOMS), 1e-4f);     // OW
    // The force field type travels with the type; the sentinel type has none
    EXPECT_EQ(sys.atom.type.count, md_array_size(sys.atom.type.ff_type));
    EXPECT_TRUE(str_eq(md_atom_type_ff_type(&sys.atom.type, sys.atom.type_idx[0]), STR_LIT("N3")));
    EXPECT_TRUE(str_empty(md_atom_type_ff_type(&sys.atom.type, 0)));

    // No coarse grained types in an all-atom system
    for (size_t t = 0; t < sys.atom.type.count; ++t) {
        EXPECT_FALSE(sys.atom.type.flags[t] & MD_FLAG_COARSE_GRAINED);
    }

    // Residue names, numbers (the waters and ions renumbered the way gmx does) and coordinates
    compare_with_gro(utest_result, &sys, &state, STR_LIT(TPR_DIR "peptide_tip3p.gro"));

    // What a caller does next: everything but the bonds is derived. One structure per molecule.
    md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL & ~MD_UTIL_INFER_BOND_BIT);
    EXPECT_EQ((size_t)PEPTIDE_BONDS + 509 * 2, sys.bond.count);
    EXPECT_EQ(1u + 509u + 4u, md_structure_count(&sys.structure));

    md_arena_allocator_destroy(arena);
}

UTEST(tpr, double_precision) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_tpr_data_t s = {0};
    md_tpr_data_t d = {0};
    ASSERT_TRUE(md_tpr_data_parse_file(&s, STR_LIT(TPR_DIR "peptide_tip3p.tpr"), arena));
    ASSERT_TRUE(md_tpr_data_parse_file(&d, STR_LIT(TPR_DIR "peptide_tip3p_double.tpr"), arena));
    EXPECT_TRUE(d.double_precision);
    ASSERT_EQ(s.num_atoms, d.num_atoms);
    ASSERT_EQ(s.num_moltypes, d.num_moltypes);
    for (size_t t = 0; t < s.num_moltypes; ++t) {
        EXPECT_EQ(s.moltypes[t].num_atoms, d.moltypes[t].num_atoms);
        EXPECT_EQ(s.moltypes[t].num_bonds, d.moltypes[t].num_bonds);
    }
    float max_diff = 0;
    for (size_t i = 0; i < s.num_atoms * 3; ++i) {
        max_diff = MAX(max_diff, fabsf(s.x[i] - d.x[i]));
    }
    EXPECT_LT(max_diff, 1e-5f);
    md_arena_allocator_destroy(arena);
}

UTEST(tpr, virtual_sites) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_system_t sys = { .alloc = arena };
    md_system_state_t state = { .alloc = arena };

    ASSERT_TRUE(md_tpr_system_init_from_file(&sys, &state, STR_LIT(TPR_DIR "peptide_tip4p.tpr")));
    EXPECT_EQ(2074u, sys.atom.count);
    // TIP4P: the two O-H bonds of the SETTLE, and the virtual site to the oxygen it hangs off
    EXPECT_EQ((size_t)PEPTIDE_BONDS + 506 * 3, sys.bond.count);

    // The first water: OW HW1 HW2 MW, the virtual site with no element
    const size_t ow = PEPTIDE_ATOMS;
    EXPECT_EQ(8, md_atom_type_atomic_number(&sys.atom.type, sys.atom.type_idx[ow]));
    EXPECT_EQ(0, md_atom_type_atomic_number(&sys.atom.type, sys.atom.type_idx[ow + 3]));
    bool found = false;
    for (size_t i = 0; i < sys.bond.count; ++i) {
        if (sys.bond.pairs[i].idx[0] == (md_atom_idx_t)ow && sys.bond.pairs[i].idx[1] == (md_atom_idx_t)(ow + 3)) found = true;
    }
    EXPECT_TRUE(found);

    compare_with_gro(utest_result, &sys, &state, STR_LIT(TPR_DIR "peptide_tip4p.gro"));
    md_arena_allocator_destroy(arena);
}

UTEST(tpr, corrupt) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_allocator_i* heap = md_get_heap_allocator();

    md_file_t file = {0};
    ASSERT_TRUE(md_file_open(&file, STR_LIT(TPR_DIR "peptide_tip4p.tpr"), MD_FILE_READ));
    const size_t size = (size_t)md_file_size(file);
    uint8_t* buf = md_alloc(heap, size);
    ASSERT_EQ(size, md_file_read(file, buf, size));
    md_file_close(&file);

    md_tpr_data_t data = {0};
    EXPECT_TRUE(md_tpr_data_parse_buffer(&data, buf, size, arena));

    // Every truncation fails, and none of them reads out of bounds
    for (size_t len = 0; len < size; len += 97) {
        EXPECT_FALSE(md_tpr_data_parse_buffer(&data, buf, len, arena));
    }

    // Corrupted values anywhere in the body are rejected or survive, but never read out of bounds
    uint8_t* bad = md_alloc(heap, size);
    uint32_t seed = 12345;
    for (int k = 0; k < 400; ++k) {
        MEMCPY(bad, buf, size);
        seed = seed * 1664525u + 1013904223u;
        const size_t pos = 100 + (seed >> 8) % (size - 100);
        bad[pos] ^= (uint8_t)(0x80 | (seed & 0x7F));
        md_tpr_data_parse_buffer(&data, bad, size, arena);
    }
    md_free(heap, bad, size);

    // Not a tpr at all
    EXPECT_FALSE(md_tpr_data_parse_buffer(&data, "This is not a run input file", 28, arena));
    EXPECT_FALSE(md_tpr_data_parse_file(&data, STR_LIT(TPR_DIR "does_not_exist.tpr"), arena));

    md_free(heap, buf, size);
    md_arena_allocator_destroy(arena);
}

// A Martini 3 system (martinized Trp-Ile, two POPC, ions, water): no atomic numbers anywhere, so
// every particle is a bead whose radius comes from its Lennard-Jones parameters
UTEST(tpr, martini) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));

    md_tpr_data_t data = {0};
    ASSERT_TRUE(md_tpr_data_parse_file(&data, STR_LIT(TPR_DIR "martini3.tpr"), arena));
    EXPECT_EQ(12u, data.num_nb_types);
    // W: sigma 0.47 nm, epsilon 4.65 kJ/mol
    const md_tpr_moltype_t* w = &data.moltypes[4];
    ASSERT_TRUE(str_eq(w->name, STR_LIT("W")));
    const md_tpr_lj_t lj = md_tpr_lj_pair(&data, w->atoms[0].type_idx, w->atoms[0].type_idx);
    EXPECT_NEAR(0.47f, powf(lj.c12 / lj.c6, 1.0f / 6.0f), 1e-4f);
    EXPECT_NEAR(0.5f * powf(2.0f, 1.0f / 6.0f) * 4.7f, md_tpr_lj_vdw_radius(lj), 1e-3f);

    // The whole table of type pairs is kept: symmetric, and Martini 3 sets its pairs explicitly rather than
    // mixing the diagonal, so some pair differs from the combination of its two types
    EXPECT_TRUE(data.nb_is_lj);
    const size_t T = data.num_nb_types;
    size_t not_mixed = 0;
    for (size_t a = 0; a < T; ++a) {
        for (size_t b = 0; b < T; ++b) {
            const md_tpr_lj_t ab = md_tpr_lj_pair(&data, a, b);
            const md_tpr_lj_t ba = md_tpr_lj_pair(&data, b, a);
            EXPECT_EQ(ab.c6, ba.c6);
            EXPECT_EQ(ab.c12, ba.c12);
            const md_tpr_lj_t aa = md_tpr_lj_pair(&data, a, a);
            const md_tpr_lj_t bb = md_tpr_lj_pair(&data, b, b);
            // Geometric mixing of c6 and c12 (combination rule 1 and 3 for sigma and epsilon alike)
            const float c6  = sqrtf(aa.c6 * bb.c6);
            const float c12 = sqrtf(aa.c12 * bb.c12);
            if (a != b && (fabsf(ab.c6 - c6) > 1e-3f * c6 || fabsf(ab.c12 - c12) > 1e-3f * c12)) not_mixed += 1;
        }
    }
    EXPECT_GT(not_mixed, (size_t)0);
    // Beyond the table
    EXPECT_EQ(0.0f, md_tpr_lj_pair(&data, T, 0).c6);

    md_system_t sys = { .alloc = arena };
    md_system_state_t state = { .alloc = arena };
    ASSERT_TRUE(md_tpr_system_init_from_data(&sys, &state, &data));
    ASSERT_EQ(619u, sys.atom.count);

    // Every type is a bead with a radius from the force field, and no element
    for (size_t t = 1; t < sys.atom.type.count; ++t) {
        EXPECT_EQ(0, sys.atom.type.z[t]);
        EXPECT_TRUE(sys.atom.type.flags[t] & MD_FLAG_COARSE_GRAINED);
        EXPECT_GT(sys.atom.type.radius[t], 1.8f);    // Tiny beads: sigma 0.34 nm
        EXPECT_LT(sys.atom.type.radius[t], 3.0f);
    }

    // BB is the same particle (Q5, 72) in both residues and shares a type. SC1 is a TC4 bead of 36 in
    // Trp and an SC2 bead of 54 in Ile: same name, different particles, different types.
    const md_atom_type_idx_t trp_bb = sys.atom.type_idx[0], trp_sc1 = sys.atom.type_idx[1];
    const md_atom_type_idx_t ile_bb = sys.atom.type_idx[6], ile_sc1 = sys.atom.type_idx[7];
    EXPECT_EQ(trp_bb, ile_bb);
    EXPECT_NE(trp_sc1, ile_sc1);
    EXPECT_TRUE(str_eq(md_atom_type_ff_type(&sys.atom.type, trp_sc1), STR_LIT("TC4")));
    EXPECT_TRUE(str_eq(md_atom_type_ff_type(&sys.atom.type, ile_sc1), STR_LIT("SC2")));
    EXPECT_TRUE(str_eq(md_atom_type_ff_type(&sys.atom.type, trp_bb),  STR_LIT("Q5")));
    EXPECT_NEAR(36.0f, sys.atom.type.mass[trp_sc1], 1e-4f);
    EXPECT_NEAR(54.0f, sys.atom.type.mass[ile_sc1], 1e-4f);
    EXPECT_NEAR(72.0f, md_atom_mass(&sys.atom, 618), 1e-4f);    // The last water

    // Trp's SC3 is a massless virtual site that interacts: still a bead with a radius
    EXPECT_EQ(0.0f, md_atom_mass(&sys.atom, 3));
    EXPECT_TRUE(sys.atom.type.flags[sys.atom.type_idx[3]] & MD_FLAG_COARSE_GRAINED);

    // The bead tables know the protein backbone
    EXPECT_TRUE(sys.atom.type.flags[trp_bb] & MD_FLAG_BACKBONE);
    EXPECT_TRUE(sys.component.flags[0] & MD_FLAG_AMINO_ACID);

    compare_with_gro(utest_result, &sys, &state, STR_LIT(TPR_DIR "martini3.gro"));

    // One structure per molecule: the peptide, two lipids, four ions, 583 waters
    md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL & ~MD_UTIL_INFER_BOND_BIT);
    EXPECT_EQ(1u + 2u + 4u + 583u, md_structure_count(&sys.structure));

    md_arena_allocator_destroy(arena);
}


// ### NON-BONDED INTERACTIONS ###
// The expected values are what 'gmx dump -s' (GROMACS 2023.3) prints for each file. epsilon-rf 'inf' is
// stored as 0.

typedef struct nb_expect_t {
    const char* file;
    int32_t scheme;
    float rlist;
    int32_t coulomb_type, coulomb_modifier;
    float rcoulomb_switch, rcoulomb, epsilon_r, epsilon_rf;
    int32_t vdw_type, vdw_modifier;
    float rvdw_switch, rvdw;
    int32_t disp_corr;
    float ewald_rtol, ewald_rtol_lj;
    int32_t ljpme_comb_rule;
    float fudge_qq;
} nb_expect_t;

UTEST(tpr, nonbonded_settings) {
    const nb_expect_t cases[] = {
        { "martini3.tpr", MD_TPR_CUTOFF_SCHEME_VERLET, 1.21f, MD_TPR_COULOMB_RF, MD_TPR_MODIFIER_POT_SHIFT, 0.0f, 1.1f, 15.0f, 0.0f,
          MD_TPR_VDW_CUT, MD_TPR_MODIFIER_POT_SHIFT, 0.0f, 1.1f, MD_TPR_DISP_CORR_NO, 1e-5f, 1e-3f, 0, 1.0f },
        { "peptide_tip3p.tpr", MD_TPR_CUTOFF_SCHEME_VERLET, 0.99f, MD_TPR_COULOMB_PME, MD_TPR_MODIFIER_POT_SHIFT, 0.0f, 0.9f, 1.0f, 0.0f,
          MD_TPR_VDW_CUT, MD_TPR_MODIFIER_POT_SHIFT, 0.0f, 0.9f, MD_TPR_DISP_CORR_NO, 1e-5f, 1e-3f, 0, 0.8333f },
        { "peptide_tip3p_double.tpr", MD_TPR_CUTOFF_SCHEME_VERLET, 0.99f, MD_TPR_COULOMB_PME, MD_TPR_MODIFIER_POT_SHIFT, 0.0f, 0.9f, 1.0f, 0.0f,
          MD_TPR_VDW_CUT, MD_TPR_MODIFIER_POT_SHIFT, 0.0f, 0.9f, MD_TPR_DISP_CORR_NO, 1e-5f, 1e-3f, 0, 0.8333f },
        { "peptide_tip4p.tpr", MD_TPR_CUTOFF_SCHEME_VERLET, 0.99f, MD_TPR_COULOMB_PME, MD_TPR_MODIFIER_POT_SHIFT, 0.0f, 0.9f, 1.0f, 0.0f,
          MD_TPR_VDW_CUT, MD_TPR_MODIFIER_POT_SHIFT, 0.0f, 0.9f, MD_TPR_DISP_CORR_NO, 1e-5f, 1e-3f, 0, 0.8333f },
        { "nb_fswitch_rf.tpr", MD_TPR_CUTOFF_SCHEME_VERLET, 1.32f, MD_TPR_COULOMB_RF, MD_TPR_MODIFIER_POT_SHIFT, 0.0f, 1.2f, 2.0f, 78.0f,
          MD_TPR_VDW_CUT, MD_TPR_MODIFIER_FORCE_SWITCH, 0.8f, 1.2f, MD_TPR_DISP_CORR_ENER_PRES, 1e-5f, 1e-3f, 0, 0.5f },
        { "nb_pswitch_cut.tpr", MD_TPR_CUTOFF_SCHEME_VERLET, 1.1f, MD_TPR_COULOMB_CUT, MD_TPR_MODIFIER_NONE, 0.0f, 1.0f, 1.0f, 0.0f,
          MD_TPR_VDW_CUT, MD_TPR_MODIFIER_POT_SWITCH, 0.9f, 1.0f, MD_TPR_DISP_CORR_NO, 1e-5f, 1e-3f, 0, 0.5f },
        { "nb_ljpme.tpr", MD_TPR_CUTOFF_SCHEME_VERLET, 1.21f, MD_TPR_COULOMB_PME, MD_TPR_MODIFIER_POT_SHIFT, 0.0f, 1.1f, 1.0f, 0.0f,
          MD_TPR_VDW_PME, MD_TPR_MODIFIER_POT_SHIFT, 0.0f, 1.1f, MD_TPR_DISP_CORR_NO, 1e-6f, 1e-4f, 1, 0.5f },
        { "nb_potshift_pme.tpr", MD_TPR_CUTOFF_SCHEME_VERLET, 1.21f, MD_TPR_COULOMB_PME, MD_TPR_MODIFIER_POT_SHIFT, 0.0f, 1.1f, 1.0f, 0.0f,
          MD_TPR_VDW_CUT, MD_TPR_MODIFIER_POT_SHIFT, 0.0f, 1.1f, MD_TPR_DISP_CORR_NO, 1e-5f, 1e-3f, 0, 0.5f },
    };

    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    for (size_t i = 0; i < ARRAY_SIZE(cases); ++i) {
        const nb_expect_t* e = &cases[i];
        char path[512];
        snprintf(path, sizeof(path), "%s%s", TPR_DIR, e->file);
        md_tpr_data_t data = {0};
        ASSERT_TRUE(md_tpr_data_parse_file(&data, str_from_cstr(path), arena));
        const md_tpr_nonbonded_t* nb = &data.nonbonded;
        EXPECT_TRUE(nb->valid);
        EXPECT_EQ(e->scheme, nb->cutoff_scheme);
        EXPECT_NEAR(e->rlist, nb->rlist, 1e-5f);
        EXPECT_EQ(e->coulomb_type, nb->coulomb_type);
        EXPECT_EQ(e->coulomb_modifier, nb->coulomb_modifier);
        EXPECT_NEAR(e->rcoulomb_switch, nb->rcoulomb_switch, 1e-6f);
        EXPECT_NEAR(e->rcoulomb, nb->rcoulomb, 1e-6f);
        EXPECT_NEAR(e->epsilon_r, nb->epsilon_r, 1e-6f);
        EXPECT_NEAR(e->epsilon_rf, nb->epsilon_rf, 1e-6f);
        EXPECT_EQ(e->vdw_type, nb->vdw_type);
        EXPECT_EQ(e->vdw_modifier, nb->vdw_modifier);
        EXPECT_NEAR(e->rvdw_switch, nb->rvdw_switch, 1e-6f);
        EXPECT_NEAR(e->rvdw, nb->rvdw, 1e-6f);
        EXPECT_EQ(e->disp_corr, nb->disp_corr);
        EXPECT_NEAR(e->ewald_rtol, nb->ewald_rtol, 1e-9f);
        EXPECT_NEAR(e->ewald_rtol_lj, nb->ewald_rtol_lj, 1e-8f);
        EXPECT_EQ(e->ljpme_comb_rule, nb->ljpme_comb_rule);
        EXPECT_NEAR(e->fudge_qq, data.fudge_qq, 1e-6f);
        EXPECT_EQ(12.0, data.repulsion_power);
        // Reading the simulation parameters leaves the rest as it was
        EXPECT_EQ(MD_TPR_PBC_XYZ, data.pbc);
        md_tpr_data_free(&data, arena);
    }
    md_arena_allocator_destroy(arena);
}

// nb_topol.top: types A (sigma 0.30, epsilon 0.5) and B (0.40, 1.0) mixed by combination rule 2, except for
// A-B which [ nonbond_params ] sets to sigma 0.25, epsilon 2.0. CHAIN is A B A B A bonded in a row with
// nrexcl 2 and an explicit [ exclusions ] 1 4; three of them, then four single atom molecules: two SOLP, two SOLN.
UTEST(tpr, nonbonded_parameters_and_exclusions) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_tpr_data_t data = {0};
    ASSERT_TRUE(md_tpr_data_parse_file(&data, STR_LIT(TPR_DIR "nb_fswitch_rf.tpr"), arena));

    // c6 = 4 epsilon sigma^6, c12 = 4 epsilon sigma^12
    ASSERT_EQ(2u, data.num_nb_types);
    ASSERT_TRUE(data.nb_is_lj);
    const md_tpr_lj_t aa = md_tpr_lj_pair(&data, 0, 0);
    const md_tpr_lj_t bb = md_tpr_lj_pair(&data, 1, 1);
    const md_tpr_lj_t ab = md_tpr_lj_pair(&data, 0, 1);
    const md_tpr_lj_t ba = md_tpr_lj_pair(&data, 1, 0);
    EXPECT_NEAR(4.0 * 0.5 * pow(0.30, 6),  aa.c6,  1e-9);
    EXPECT_NEAR(4.0 * 0.5 * pow(0.30, 12), aa.c12, 1e-12);
    EXPECT_NEAR(4.0 * 1.0 * pow(0.40, 6),  bb.c6,  1e-8);
    // The override, not the combination (sigma 0.35, epsilon 0.707)
    EXPECT_NEAR(4.0 * 2.0 * pow(0.25, 6),  ab.c6,  1e-9);
    EXPECT_NEAR(4.0 * 2.0 * pow(0.25, 12), ab.c12, 1e-12);
    EXPECT_EQ(ab.c6, ba.c6);
    EXPECT_EQ(ab.c12, ba.c12);

    ASSERT_EQ(3u, data.num_moltypes);
    const md_tpr_moltype_t* chain = &data.moltypes[0];
    const md_tpr_moltype_t* solo  = &data.moltypes[1];
    ASSERT_TRUE(str_eq(solo->name, STR_LIT("SOLP")));
    ASSERT_TRUE(str_eq(chain->name, STR_LIT("CHAIN")));
    ASSERT_EQ(5u, chain->num_atoms);
    ASSERT_TRUE(chain->excl_offset && chain->excl);
    // Within two bonds, and 1-4 given explicitly. The atom itself is not listed.
    const uint32_t expected[5][4] = { {1, 2, 3}, {0, 2, 3}, {0, 1, 3, 4}, {0, 1, 2, 4}, {2, 3} };
    const uint32_t expected_len[5] = { 3, 3, 4, 4, 2 };
    for (uint32_t i = 0; i < 5; ++i) {
        ASSERT_EQ(expected_len[i], chain->excl_offset[i + 1] - chain->excl_offset[i]);
        for (uint32_t k = 0; k < expected_len[i]; ++k) {
            EXPECT_EQ(expected[i][k], chain->excl[chain->excl_offset[i] + k]);
        }
    }
    // A single atom is only excluded from itself, which is not listed
    EXPECT_TRUE(solo->excl_offset == NULL);

    // Global: the second chain is atoms 5..9, the solos 15..18
    EXPECT_TRUE(md_tpr_atoms_excluded(&data, 5, 8));
    EXPECT_TRUE(md_tpr_atoms_excluded(&data, 8, 5));
    EXPECT_FALSE(md_tpr_atoms_excluded(&data, 5, 9));
    EXPECT_TRUE(md_tpr_atoms_excluded(&data, 7, 9));
    EXPECT_FALSE(md_tpr_atoms_excluded(&data, 4, 5));     // Last of one chain, first of the next
    EXPECT_FALSE(md_tpr_atoms_excluded(&data, 0, 15));
    EXPECT_FALSE(md_tpr_atoms_excluded(&data, 15, 16));
    EXPECT_TRUE(md_tpr_atoms_excluded(&data, 16, 16));
    EXPECT_FALSE(md_tpr_atoms_excluded(&data, 0, 19));    // Beyond the system

    // SETTLE water excludes its three atoms from each other
    md_tpr_data_t pep = {0};
    ASSERT_TRUE(md_tpr_data_parse_file(&pep, STR_LIT(TPR_DIR "peptide_tip3p.tpr"), arena));
    const md_tpr_moltype_t* sol = &pep.moltypes[1];
    ASSERT_TRUE(str_eq(sol->name, STR_LIT("SOL")));
    ASSERT_TRUE(sol->excl_offset != NULL);
    for (uint32_t i = 0; i < 3; ++i) {
        EXPECT_EQ(2u, sol->excl_offset[i + 1] - sol->excl_offset[i]);
    }
    // The peptide excludes up to three bonds (amber: nrexcl 3): the atoms bonded to atom 0 and theirs
    const md_tpr_moltype_t* prot = &pep.moltypes[0];
    ASSERT_TRUE(prot->excl_offset != NULL);
    for (size_t b = 0; b < prot->num_bonds; ++b) {
        const int32_t i = prot->bonds[b].idx[0], j = prot->bonds[b].idx[1];
        bool found = false;
        for (uint32_t k = prot->excl_offset[i]; k < prot->excl_offset[i + 1]; ++k) found |= prot->excl[k] == (uint32_t)j;
        EXPECT_TRUE(found);
    }

    md_tpr_data_free(&pep, arena);
    md_tpr_data_free(&data, arena);
    md_arena_allocator_destroy(arena);
}
