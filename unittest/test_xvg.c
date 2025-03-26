#include "utest.h"
#include <string.h>

#include <md_xvg.h>
#include <core/md_allocator.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include <core/md_array.h>

UTEST(xvg, rdf) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/rdf.xvg");
    md_xvg_t xvg = {0};
    bool result = md_xvg_parse_file(&xvg, path, md_get_heap_allocator());
    ASSERT_TRUE(result);
    
    EXPECT_EQ(2,    xvg.num_fields);
    EXPECT_EQ(1024, xvg.num_values);

    EXPECT_NEAR(0.010, xvg.fields[0][10], 1.0e-6f);
    EXPECT_NEAR(0.000, xvg.fields[1][10], 1.0e-6f);
    
    EXPECT_NEAR(0.049, xvg.fields[0][50], 1.0e-6f);
    EXPECT_NEAR(0.000, xvg.fields[1][50], 1.0e-6f);

    ASSERT_EQ(1, xvg.header_info.num_legends);
    EXPECT_TRUE(str_eq(xvg.header_info.legends[0], STR_LIT("OW")));

    EXPECT_TRUE(str_eq(xvg.header_info.title, STR_LIT("Radial distribution")));
    EXPECT_TRUE(str_eq(xvg.header_info.xaxis_label, STR_LIT("r (nm)")));
    EXPECT_TRUE(str_eq(xvg.header_info.yaxis_label, STR_LIT("g(r)")));

    md_xvg_free(&xvg, md_get_heap_allocator());
}

UTEST(xvg, energy) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/energy.xvg");
    md_xvg_t xvg = {0};
    bool result = md_xvg_parse_file(&xvg, path, md_get_heap_allocator());
    ASSERT_TRUE(result);

    EXPECT_EQ(5,   xvg.num_fields);
    EXPECT_EQ(356, xvg.num_values);

    //    18.000000  -427.499512  -170.860535  -3487.963135  -2470.888916
    EXPECT_NEAR(   18.000000, xvg.fields[0][9], 1.0e-6f);
    EXPECT_NEAR( -427.499512, xvg.fields[1][9], 1.0e-6f);
    EXPECT_NEAR( -170.860535, xvg.fields[2][9], 1.0e-6f);
    EXPECT_NEAR(-3487.963135, xvg.fields[3][9], 1.0e-6f);
    EXPECT_NEAR(-2470.888916, xvg.fields[4][9], 1.0e-6f);

    //   710.000000  -414.190063  -272.062927  -3277.823975  -2500.806641
    EXPECT_NEAR(  710.000000, xvg.fields[0][355], 1.0e-6f);
    EXPECT_NEAR( -414.190063, xvg.fields[1][355], 1.0e-6f);
    EXPECT_NEAR( -272.062927, xvg.fields[2][355], 1.0e-6f);
    EXPECT_NEAR(-3277.823975, xvg.fields[3][355], 1.0e-6f);
    EXPECT_NEAR(-2500.806641, xvg.fields[4][355], 1.0e-6f);

    ASSERT_EQ(4, xvg.header_info.num_legends);
    EXPECT_TRUE(str_eq(xvg.header_info.legends[0], STR_LIT("Coul-SR:2S29-2S29")));
    EXPECT_TRUE(str_eq(xvg.header_info.legends[1], STR_LIT("LJ-SR:2S29-2S29")));
    EXPECT_TRUE(str_eq(xvg.header_info.legends[2], STR_LIT("Coul-SR:2S29-SOL")));
    EXPECT_TRUE(str_eq(xvg.header_info.legends[3], STR_LIT("LJ-SR:2S29-SOL")));

    EXPECT_TRUE(str_eq(xvg.header_info.title, STR_LIT("GROMACS Energies")));
    EXPECT_TRUE(str_eq(xvg.header_info.xaxis_label, STR_LIT("Time (ps)")));
    EXPECT_TRUE(str_eq(xvg.header_info.yaxis_label, STR_LIT("(kJ/mol)")));

    md_xvg_free(&xvg, md_get_heap_allocator());
}
