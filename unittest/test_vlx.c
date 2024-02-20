#include "utest.h"

#include <md_vlx.h>
#include <core/md_allocator.h>
#include <core/md_str.h>

static const int ref_iter[] = {
    1,2,3,4,5,6,7,8,9,10,11,12
};
static const double ref_ener_tot[] = {
    -444.502770501781,
    -444.502826985951,
    -444.518269328072,
    -444.518408360201,
    -444.518497179843,
    -444.518500090662,
    -444.518500772981,
    -444.518500782306,
    -444.518500783079,
    -444.518500783141,
    -444.518500783156,
    -444.518500783158,
};
static const double ref_ener_change[] = {
    0.0000000000,
    -0.0000564842,
    -0.0154423421,
    -0.0001390321,
    -0.0000888196,
    -0.0000029108,
    -0.0000006823,
    -0.0000000093,
    -0.0000000008,
    -0.0000000001,
    -0.0000000000,
    -0.0000000000,
};
static const double ref_grad_norm[] = {
    0.33661898,
    0.33937059,
    0.04132164,
    0.02546891,
    0.00508750,
    0.00230342,
    0.00028359,
    0.00008491,
    0.00002482,
    0.00000991,
    0.00000406,
    0.00000094,
};
static const double ref_max_grad[] = {
    0.01168705,
    0.01156775,
    0.00163688,
    0.00092076,
    0.00023894,
    0.00009228,
    0.00000913,
    0.00000327,
    0.00000086,
    0.00000045,
    0.00000018,
    0.00000003,
};
static const double ref_density_change[] = {
    0.00000000,
    0.31704185,
    0.17166350,
    0.02943676,
    0.01204202,
    0.00304038,
    0.00107142,
    0.00015964,
    0.00003900,
    0.00001305,
    0.00000550,
    0.00000235,
};

UTEST(vlx, vlx_parse) {
    str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/vlx/mol.out");
    md_vlx_data_t vlx = {0};
    bool result = md_vlx_data_parse_file(&vlx, path, md_heap_allocator);
    ASSERT_TRUE(result);

    EXPECT_EQ(0,  vlx.geom.molecular_charge);
    EXPECT_EQ(26, vlx.geom.num_atoms);
    EXPECT_EQ(1,  vlx.geom.spin_multiplicity);
    EXPECT_EQ(41, vlx.geom.num_alpha_electrons);
    EXPECT_EQ(41, vlx.geom.num_beta_electrons);

    EXPECT_NEAR(-3.259400000000, vlx.geom.coord_x[0], 1.0e-5);
    EXPECT_NEAR( 0.145200000000, vlx.geom.coord_y[0], 1.0e-5);
    EXPECT_NEAR(-0.048400000000, vlx.geom.coord_z[0], 1.0e-5);

    EXPECT_TRUE(str_eq(vlx.basis.ident, STR_LIT("DEF2-SVP")));
    EXPECT_EQ(229, vlx.basis.num_contracted_basis_functions);
    EXPECT_EQ(369, vlx.basis.num_primitive_basis_functions);

    EXPECT_EQ(-444.5185007832,  vlx.scf.total_energy);
    EXPECT_EQ(-1041.3439889091, vlx.scf.electronic_energy);
    EXPECT_EQ(596.8254881260,   vlx.scf.nuclear_repulsion_energy);
    EXPECT_EQ(0.0000009436,     vlx.scf.gradient_norm);

    ASSERT_EQ(12, vlx.scf.iter.count);
    for (size_t i = 0; i < vlx.scf.iter.count; ++i) {
        EXPECT_EQ(ref_iter[i], vlx.scf.iter.iteration[i]);
        EXPECT_NEAR(ref_ener_tot[i], vlx.scf.iter.energy_total[i], 1.0e-5);
        EXPECT_NEAR(ref_ener_change[i], vlx.scf.iter.energy_change[i], 1.0e-5);
        EXPECT_NEAR(ref_grad_norm[i], vlx.scf.iter.gradient_norm[i], 1.0e-5);
        EXPECT_NEAR(ref_max_grad[i], vlx.scf.iter.max_gradient[i], 1.0e-5);
        EXPECT_NEAR(ref_density_change[i], vlx.scf.iter.density_change[i], 1.0e-5);
    }

    // @TODO: Test RSP

    md_vlx_data_free(&vlx);
}
