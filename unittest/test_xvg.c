#include "utest.h"
#include <string.h>

#include <md_xvg.h>
#include <core/md_allocator.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include <core/md_array.h>

UTEST(xvg, rdf) {
    str_t path = STR(MD_UNITTEST_DATA_DIR "/rdf.xvg");
    md_xvg_t xvg = {0};
    bool result = md_xvg_init_from_file(&xvg, path, default_allocator);
    ASSERT_TRUE(result);
    
    EXPECT_EQ(2,    xvg.num_cols);
    EXPECT_EQ(1024, xvg.num_rows);

    const float* row10 = md_xvg_row(&xvg, 10);
    EXPECT_NEAR(0.010, row10[0], 1.0e-6f);
    EXPECT_NEAR(0.000, row10[1], 1.0e-6f);
    
    const float* row50 = md_xvg_row(&xvg, 50);
    EXPECT_NEAR(0.049, row50[0], 1.0e-6f);
    EXPECT_NEAR(0.000, row50[1], 1.0e-6f);

    md_xvg_free(&xvg, default_allocator);
}

UTEST(xvg, energy) {
    str_t path = STR(MD_UNITTEST_DATA_DIR "/energy.xvg");
    md_xvg_t xvg = {0};
    bool result = md_xvg_init_from_file(&xvg, path, default_allocator);
    ASSERT_TRUE(result);

    EXPECT_EQ(5,   xvg.num_cols);
    EXPECT_EQ(356, xvg.num_rows);

    const float* row9 = md_xvg_row(&xvg, 9);
    //    18.000000  -427.499512  -170.860535  -3487.963135  -2470.888916
    EXPECT_NEAR(   18.000000, row9[0], 1.0e-6f);
    EXPECT_NEAR( -427.499512, row9[1], 1.0e-6f);
    EXPECT_NEAR( -170.860535, row9[2], 1.0e-6f);
    EXPECT_NEAR(-3487.963135, row9[3], 1.0e-6f);
    EXPECT_NEAR(-2470.888916, row9[4], 1.0e-6f);

    //   710.000000  -414.190063  -272.062927  -3277.823975  -2500.806641
    const float* row355 = md_xvg_row(&xvg, 355);
    EXPECT_NEAR(  710.000000, row355[0], 1.0e-6f);
    EXPECT_NEAR( -414.190063, row355[1], 1.0e-6f);
    EXPECT_NEAR( -272.062927, row355[2], 1.0e-6f);
    EXPECT_NEAR(-3277.823975, row355[3], 1.0e-6f);
    EXPECT_NEAR(-2500.806641, row355[4], 1.0e-6f);

    md_xvg_free(&xvg, default_allocator);
}
