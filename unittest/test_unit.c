#include "utest.h"

#include <core/md_unit.h>
#include <core/md_str.h>

#define UNIT_KILOGRAM       {.base = {.dim = {.mass = 1,}},     .mult = 1.0}
#define UNIT_NANOMETER      {.base = {.dim = {.length = 1,}},   .mult = 1e-9}

#define UNIT_M_S2			{.base = {.dim = {.length = 1, .time = -2}},    .mult = 1.0}
#define UNIT_KJ_MOL			{.base = {.dim = {.mass = 1, .length =  2, .time = -2, .mole = -1}}, .mult = 1000.0}

#define UNIT_PASCAL         {.base = {.dim = {.mass = 1, .length = -1, .time = -2,}}, .mult = 1.0}
#define UNIT_BAR            {.base = {.dim = {.mass = 1, .length = -1, .time = -2,}}, .mult = 1.0e5}
#define UNIT_JOULE          {.base = {.dim = {.mass = 1, .length =  2, .time = -2,}}, .mult = 1.0}
#define UNIT_NEWTON         {.base = {.dim = {.mass = 1, .length =  1, .time = -2,}}, .mult = 1.0}

#define UNIT_VOLT           {.base = {.dim = {.mass = 1, .length =  2, .time = -3, .current = -1}}, .mult = 1.0}
#define UNIT_COULOMB        {.base = {.dim = {.time = 1, .current = 1}}, .mult = 1.0}
#define UNIT_AMPEREHOUR     {.base = {.dim = {.time = 1, .current = 1}}, .mult = 3600.0}

static md_unit_t parse(str_t str) {
	md_unit_t unit;
	md_unit_parse(&unit, str);
	return unit;
}

UTEST(unit, print) {
	char buf[512] = {0};

	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_KILOGRAM);		EXPECT_STREQ("kg", buf);
	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_PASCAL);		EXPECT_STREQ("Pa", buf);
	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_NANOMETER);		EXPECT_STREQ("nm", buf);
	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_BAR);			EXPECT_STREQ("bar", buf);
	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_M_S2);			EXPECT_STREQ("m/s^2", buf);
	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_KJ_MOL);		EXPECT_STREQ("kJ/mol", buf);
	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_NEWTON);		EXPECT_STREQ("N", buf);
	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_VOLT);			EXPECT_STREQ("V", buf);
	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_COULOMB);		EXPECT_STREQ("C", buf);
	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_AMPEREHOUR);	EXPECT_STREQ("Ah", buf);

	// A unit without dimensions has no representation
	md_unit_print(buf, sizeof(buf), (md_unit_t){0});					EXPECT_STREQ("", buf);
	md_unit_print(buf, sizeof(buf), md_unit_none());					EXPECT_STREQ("", buf);
	md_unit_print(buf, sizeof(buf), md_unit_scl(md_unit_none(), 1e3));	EXPECT_STREQ("", buf);

	// The prefix binds before the exponent, nm^2 is (1e-9 m)^2
	md_unit_print(buf, sizeof(buf), md_unit_inv(md_unit_nanometer()));					EXPECT_STREQ("1/nm", buf);
	md_unit_print(buf, sizeof(buf), md_unit_pow(md_unit_nanometer(), 2));				EXPECT_STREQ("nm^2", buf);
	md_unit_print(buf, sizeof(buf), md_unit_inv(md_unit_pow(md_unit_nanosecond(), 2)));	EXPECT_STREQ("1/ns^2", buf);
	md_unit_print(buf, sizeof(buf), md_unit_pow(md_unit_angstrom(), 3));				EXPECT_STREQ("\xc3\x85^3", buf);

	// Factors are separated and a denominator holding several factors is parenthesized
	md_unit_print(buf, sizeof(buf), md_unit_div(md_unit_joule(), md_unit_mul(md_unit_mole(), md_unit_kelvin())));
	EXPECT_STREQ("J/(mol*K)", buf);
	md_unit_print(buf, sizeof(buf), md_unit_div(md_unit_count(), md_unit_pow(md_unit_angstrom(), 3)));
	EXPECT_STREQ("count/\xc3\x85^3", buf);
}

UTEST(unit, parse) {
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_KILOGRAM,		parse(STR_LIT("kg"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_PASCAL,		parse(STR_LIT("Pa"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_NANOMETER,	parse(STR_LIT("nm"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_BAR,			parse(STR_LIT("bar"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_M_S2,			parse(STR_LIT("m/s^2"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_M_S2,			parse(STR_LIT("m*s^-2"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_KJ_MOL,		parse(STR_LIT("kJ/mol"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_NEWTON,		parse(STR_LIT("kg*m/s^2"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_NEWTON,		parse(STR_LIT("m kg/s^2"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_NEWTON,		parse(STR_LIT("N"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_VOLT,			parse(STR_LIT("V"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_COULOMB,		parse(STR_LIT("C"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_AMPEREHOUR,	parse(STR_LIT("Ah"))));
	EXPECT_TRUE(md_unit_equal(md_unit_inv(md_unit_second()),parse(STR_LIT("1/s"))));
	EXPECT_TRUE(md_unit_equal(md_unit_scl(md_unit_meter(), 1e-6), parse(STR_LIT("um"))));

	// Parenthesis, explicit numeric scale and left associative division
	const md_unit_t J_mol_K = md_unit_div(md_unit_joule(), md_unit_mul(md_unit_mole(), md_unit_kelvin()));
	EXPECT_TRUE(md_unit_equal(J_mol_K,						parse(STR_LIT("J/(mol*K)"))));
	EXPECT_TRUE(md_unit_equal(J_mol_K,						parse(STR_LIT("J/mol/K"))));
	EXPECT_TRUE(md_unit_equal(md_unit_angstrom(),			parse(STR_LIT("1e-10*m"))));
	EXPECT_TRUE(md_unit_equal(md_unit_scl(md_unit_joule(), 4.184), parse(STR_LIT("4.184*J"))));
}

UTEST(unit, parse_reports_failure) {
	md_unit_t unit;

	EXPECT_TRUE(md_unit_parse(&unit, STR_LIT("nm")));
	EXPECT_TRUE(md_unit_equal(unit, md_unit_nanometer()));

	// An empty string is a successfully parsed 'no unit'
	EXPECT_TRUE(md_unit_parse(&unit, STR_LIT("")));
	EXPECT_TRUE(md_unit_is_none(unit));

	EXPECT_FALSE(md_unit_parse(&unit, STR_LIT("banana")));
	EXPECT_TRUE(md_unit_is_none(unit));

	EXPECT_FALSE(md_unit_parse(&unit, STR_LIT("q m")));
	EXPECT_TRUE(md_unit_is_none(unit));
}

UTEST(unit, print_parse_roundtrip) {
	const md_unit_t base[] = {
		md_unit_meter(), md_unit_kilogram(), md_unit_second(), md_unit_ampere(), md_unit_mole(),
		md_unit_kelvin(), md_unit_radian(), md_unit_count(), md_unit_angstrom(), md_unit_nanometer(),
		md_unit_picosecond(), md_unit_femtosecond(), md_unit_joule(), md_unit_electronvolt(),
		md_unit_hertz(), md_unit_pascal(), md_unit_bar(), md_unit_degree(),
	};

	// Anything we print has to read back as the same unit
	for (size_t i = 0; i < ARRAY_SIZE(base); ++i) {
		for (size_t j = 0; j < ARRAY_SIZE(base); ++j) {
			for (int exp = -2; exp <= 2; ++exp) {
				const md_unit_t unit = md_unit_mul(base[i], md_unit_pow(base[j], exp));
				if (md_unit_is_none(unit)) continue;

				char buf[256];
				const size_t len = md_unit_print(buf, sizeof(buf), unit);
				EXPECT_GT(len, 0U);

				md_unit_t parsed;
				EXPECT_TRUE(md_unit_parse(&parsed, (str_t){buf, len}));
				EXPECT_TRUE(md_unit_equal(parsed, unit));
			}
		}
	}
}

UTEST(unit, conversion_factor) {
	double factor = 0.0;
	EXPECT_TRUE(md_unit_conversion_factor(&factor, md_unit_degree(), md_unit_radian()));
	EXPECT_NEAR(factor, DEG_TO_RAD(1.0), 1e-15);

	EXPECT_TRUE(md_unit_conversion_factor(&factor, md_unit_angstrom(), md_unit_meter()));
	EXPECT_NEAR(factor, 1e-10, 1e-25);

	EXPECT_TRUE(md_unit_conversion_factor(&factor, md_unit_nanometer(), md_unit_angstrom()));
	EXPECT_NEAR(factor, 10.0, 1e-9);

	EXPECT_TRUE(md_unit_conversion_factor(&factor, md_unit_none(), md_unit_none()));
	EXPECT_EQ(factor, 1.0);

	// Mismatching dimensions are not convertible
	EXPECT_FALSE(md_unit_conversion_factor(&factor, md_unit_meter(), md_unit_second()));
	EXPECT_FALSE(md_unit_conversion_factor(&factor, md_unit_none(), md_unit_radian()));
}

UTEST(unit, zero_init_sentinel) {
	const md_unit_t zero = {0};

	// The zero initialized unit is a valid unit without dimensions and with a scale of one
	EXPECT_EQ(md_unit_scale(zero), 1.0);
	EXPECT_TRUE(md_unit_is_none(zero));
	EXPECT_TRUE(md_unit_equal(zero, md_unit_none()));

	// It is the identity element under multiplication and division
	EXPECT_TRUE(md_unit_equal(md_unit_mul(zero, md_unit_meter()), md_unit_meter()));
	EXPECT_TRUE(md_unit_equal(md_unit_mul(md_unit_meter(), zero), md_unit_meter()));
	EXPECT_TRUE(md_unit_equal(md_unit_div(md_unit_meter(), zero), md_unit_meter()));
	EXPECT_TRUE(md_unit_equal(md_unit_div(zero, md_unit_second()), md_unit_inv(md_unit_second())));
	EXPECT_TRUE(md_unit_equal(md_unit_pow(zero, 3), zero));
	EXPECT_TRUE(md_unit_equal(md_unit_inv(zero), zero));

	// A unit without dimensions may still carry a scale
	const md_unit_t kilo = md_unit_scl(zero, 1e3);
	EXPECT_EQ(md_unit_scale(kilo), 1e3);
	EXPECT_FALSE(md_unit_equal(kilo, zero));
	EXPECT_TRUE(md_unit_equal(md_unit_div(md_unit_meter(), md_unit_scl(md_unit_meter(), 1e-3)), kilo));
}

UTEST(unit, inverse_and_power_scale) {
	md_unit_t inv_ns = md_unit_inv(md_unit_nanosecond());
	EXPECT_TRUE(md_unit_base_equal(inv_ns, md_unit_pow(md_unit_second(), -1)));
	EXPECT_NEAR(md_unit_scale(inv_ns), 1e9, 1e-3);

	md_unit_t nm3 = md_unit_pow(md_unit_nanometer(), 3);
	EXPECT_TRUE(md_unit_base_equal(nm3, md_unit_pow(md_unit_meter(), 3)));
	EXPECT_NEAR(md_unit_scale(nm3), 1e-27, 1e-36);
}

UTEST(unit, add_sub_require_exact_units) {
	md_unit_t m = md_unit_meter();
	md_unit_t mm = md_unit_scl(md_unit_meter(), 1e-3);

	EXPECT_TRUE(md_unit_is_none(md_unit_add(m, mm)));
	EXPECT_TRUE(md_unit_is_none(md_unit_sub(m, mm)));
	EXPECT_TRUE(md_unit_equal(md_unit_add(m, m), m));
	EXPECT_TRUE(md_unit_equal(md_unit_sub(m, m), m));
}
