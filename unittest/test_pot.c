#include "utest.h"

#include <md_pot.h>
#include <core/md_allocator.h>
#include <core/md_str.h>

UTEST(pot, water_pe_npe) {
    md_pot_t pot = {0};
    md_allocator_i* alloc = md_get_heap_allocator();

    bool result = md_pot_parse_file(&pot, STR_LIT(MD_UNITTEST_DATA_DIR "/pot/water_pe_npe.pot"), alloc);
    ASSERT_TRUE(result);

    EXPECT_EQ((int)MD_POT_UNIT_ANGSTROM, (int)pot.unit);
    EXPECT_NEAR(1.0, md_pot_unit_to_angstrom(pot.unit), 1e-12);

    // 50 water molecules, 39 polarizably embedded (HOH_pe) and 11 not (HOH_npe)
    ASSERT_EQ(150u, (unsigned)pot.num_sites);

    // First site
    EXPECT_STREQ("O", pot.sites[0].element);
    EXPECT_NEAR(-5.672, pot.sites[0].coord[0], 1e-6);
    EXPECT_NEAR( 2.390, pot.sites[0].coord[1], 1e-6);
    EXPECT_NEAR(-4.911, pot.sites[0].coord[2], 1e-6);
    EXPECT_EQ(1, pot.sites[0].fragment_id);
    EXPECT_TRUE(str_eq(pot.sites[0].fragment_name, STR_LIT("HOH_pe")));
    EXPECT_TRUE(str_eq(pot.sites[0].atom_name,     STR_LIT("OW")));

    EXPECT_STREQ("H", pot.sites[2].element);
    EXPECT_TRUE(str_eq(pot.sites[2].atom_name, STR_LIT("H2")));

    // First non-polarizable site (line 121 of the file -> site index 117)
    EXPECT_TRUE(str_eq(pot.sites[117].fragment_name, STR_LIT("HOH_npe")));
    EXPECT_EQ(8, pot.sites[117].fragment_id);
    EXPECT_NEAR(-11.932, pot.sites[117].coord[0], 1e-6);

    // Last site
    EXPECT_STREQ("H", pot.sites[149].element);
    EXPECT_NEAR(-5.908, pot.sites[149].coord[0], 1e-6);
    EXPECT_NEAR( 6.663, pot.sites[149].coord[1], 1e-6);
    EXPECT_NEAR( 2.605, pot.sites[149].coord[2], 1e-6);
    EXPECT_EQ(49, pot.sites[149].fragment_id);
    EXPECT_TRUE(str_eq(pot.sites[149].fragment_name, STR_LIT("HOH_npe")));

    // Fragment names are interned, so equal names share storage
    EXPECT_EQ(pot.sites[0].fragment_name.ptr, pot.sites[3].fragment_name.ptr);
    EXPECT_NE(pot.sites[0].fragment_name.ptr, pot.sites[117].fragment_name.ptr);
    // HOH_pe, HOH_npe, OW, H1, H2
    EXPECT_EQ(5u, (unsigned)pot.num_strings);

    // Charges: two fragment types, three atoms each
    ASSERT_EQ(6u, (unsigned)pot.num_charges);
    EXPECT_STREQ("O", pot.charges[0].element);
    EXPECT_NEAR(-0.67444, pot.charges[0].charge, 1e-8);
    EXPECT_TRUE(str_eq(pot.charges[0].fragment_name, STR_LIT("HOH_pe")));
    EXPECT_NEAR( 0.33722, pot.charges[1].charge, 1e-8);
    EXPECT_NEAR( 0.33722, pot.charges[2].charge, 1e-8);
    EXPECT_STREQ("O", pot.charges[3].element);
    EXPECT_NEAR(-0.83400, pot.charges[3].charge, 1e-8);
    EXPECT_TRUE(str_eq(pot.charges[3].fragment_name, STR_LIT("HOH_npe")));
    EXPECT_NEAR( 0.41700, pot.charges[5].charge, 1e-8);

    // Polarizabilities: only the HOH_pe fragment type is present
    ASSERT_EQ(3u, (unsigned)pot.num_polarizabilities);
    EXPECT_STREQ("O", pot.polarizabilities[0].element);
    EXPECT_TRUE(str_eq(pot.polarizabilities[0].fragment_name, STR_LIT("HOH_pe")));
    EXPECT_NEAR(5.73935, pot.polarizabilities[0].alpha[MD_POT_ALPHA_XX], 1e-8);
    EXPECT_NEAR(0.0,     pot.polarizabilities[0].alpha[MD_POT_ALPHA_XY], 1e-8);
    EXPECT_NEAR(0.0,     pot.polarizabilities[0].alpha[MD_POT_ALPHA_XZ], 1e-8);
    EXPECT_NEAR(5.73935, pot.polarizabilities[0].alpha[MD_POT_ALPHA_YY], 1e-8);
    EXPECT_NEAR(0.0,     pot.polarizabilities[0].alpha[MD_POT_ALPHA_YZ], 1e-8);
    EXPECT_NEAR(5.73935, pot.polarizabilities[0].alpha[MD_POT_ALPHA_ZZ], 1e-8);
    EXPECT_NEAR(2.30839, pot.polarizabilities[1].alpha[MD_POT_ALPHA_XX], 1e-8);
    EXPECT_NEAR(2.30839, pot.polarizabilities[2].alpha[MD_POT_ALPHA_ZZ], 1e-8);

    double m[3][3];
    md_pot_alpha_to_mat3(m, pot.polarizabilities[0].alpha);
    EXPECT_NEAR(5.73935, m[0][0], 1e-8);
    EXPECT_NEAR(5.73935, m[1][1], 1e-8);
    EXPECT_NEAR(5.73935, m[2][2], 1e-8);
    EXPECT_NEAR(0.0,     m[0][1], 1e-8);
    EXPECT_NEAR(0.0,     m[2][0], 1e-8);

    md_pot_free(&pot, alloc);
}

UTEST(pot, minimal_and_units) {
    md_pot_t pot = {0};
    md_allocator_i* alloc = md_get_heap_allocator();

    // No atom-name column, bohr units, a comment and an unsupported section
    str_t src = STR_LIT(
        "# a comment\n"
        "@environment\n"
        "units: au\n"
        "xyz:\n"
        "O   0.000000  0.000000  0.000000  water  1\n"
        "H   1.000000  0.000000  0.000000  water  1\n"
        "@end\n"
        "\n"
        "@exclusion_lists\n"
        "1 2\n"
        "@end\n"
        "\n"
        "@charges\n"
        "O  -0.60000000  water\n"
        "H   0.30000000  water\n"
        "@end\n");

    ASSERT_TRUE(md_pot_parse_str(&pot, src, alloc));
    EXPECT_EQ((int)MD_POT_UNIT_BOHR, (int)pot.unit);
    EXPECT_NEAR(0.52917721, md_pot_unit_to_angstrom(pot.unit), 1e-8);

    ASSERT_EQ(2u, (unsigned)pot.num_sites);
    EXPECT_EQ(1, pot.sites[0].fragment_id);
    EXPECT_TRUE(str_empty(pot.sites[0].atom_name));
    EXPECT_TRUE(str_eq(pot.sites[1].fragment_name, STR_LIT("water")));
    EXPECT_NEAR(1.0, pot.sites[1].coord[0], 1e-8);

    ASSERT_EQ(2u, (unsigned)pot.num_charges);
    EXPECT_NEAR(-0.6, pot.charges[0].charge, 1e-8);
    EXPECT_EQ(0u, (unsigned)pot.num_polarizabilities);

    md_pot_free(&pot, alloc);
}

UTEST(pot, malformed) {
    md_allocator_i* alloc = md_get_heap_allocator();

    {   // Too few fields in @environment
        md_pot_t pot = {0};
        str_t src = STR_LIT("@environment\nunits: angstrom\nxyz:\nO 0.0 0.0\n@end\n");
        EXPECT_FALSE(md_pot_parse_str(&pot, src, alloc));
        EXPECT_EQ(0u, (unsigned)pot.num_sites);
    }
    {   // Coordinate that is not a number
        md_pot_t pot = {0};
        str_t src = STR_LIT("@environment\nunits: angstrom\nxyz:\nO 0.0 0.0 abc\n@end\n");
        EXPECT_FALSE(md_pot_parse_str(&pot, src, alloc));
    }
    {   // Missing 'xyz:' directive
        md_pot_t pot = {0};
        str_t src = STR_LIT("@environment\nunits: angstrom\nO 0.0 0.0 0.0\n@end\n");
        EXPECT_FALSE(md_pot_parse_str(&pot, src, alloc));
    }
    {   // No @environment at all
        md_pot_t pot = {0};
        str_t src = STR_LIT("@charges\nO -0.6 water\n@end\n");
        EXPECT_FALSE(md_pot_parse_str(&pot, src, alloc));
    }
}
