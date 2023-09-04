#include "utest.h"

#include <core/md_unit.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_str.h>

typedef struct {
	md_unit_t unit;
	str_t     str;
} item_t;

#define ADD(unit, cstr) md_array_push(table, ((item_t){unit, STR(cstr)}), alloc)

#define UNIT_KILOGRAM       {.base = {.dim = {.mass = 1,}},     .mult = 1.0}
#define UNIT_NANOMETER      {.base = {.dim = {.length = 1,}},   .mult = 1e-9}

#define UNIT_M_S2			{.base = {.dim = {.length = 1, .time = -2}},    .mult = 1.0}
#define UNIT_KJ_MOL			{.base = {.dim = {.mass = 1, .length =  2, .time = -2, .mole = -1}}, .mult = 1000.0}

#define UNIT_PASCAL         {.base = {.dim = {.mass = 1, .length = -1, .time = -2,}}, .mult = 1.0}
#define UNIT_BAR            {.base = {.dim = {.mass = 1, .length = -1, .time = -2,}}, .mult = 1.0e5}
#define UNIT_JOULE          {.base = {.dim = {.mass = 1, .length =  2, .time = -2,}}, .mult = 1.0}
#define UNIT_NEWTON         {.base = {.dim = {.mass = 1, .length =  1, .time = -2,}}, .mult = 1.0}

#define UNIT_VOLT           {.base = {.dim = {.mass = 1, .length =  2, .time = -3, .current = -1}}, .mult = 1.0}
#define UNIT_WATT           {.base = {.dim = {.mass = 1, .length =  2, .time = -3}}, .mult = 1.0}
#define UNIT_COULOMB        {.base = {.dim = {.time = 1, .current = 1}}, .mult = 1.0}
#define UNIT_AMPEREHOUR     {.base = {.dim = {.time = 1, .current = 1}}, .mult = 3600.0}

item_t* build_test_table(md_allocator_i* alloc) {
	item_t* table = 0;

	ADD(UNIT_KILOGRAM,		"kg");
	ADD(UNIT_PASCAL,		"Pa");
	ADD(UNIT_NANOMETER,		"nm");
	ADD(UNIT_BAR,			"bar");
	ADD(UNIT_M_S2,			"m/s^2");
	ADD(UNIT_KJ_MOL,		"kJ/mol");
	ADD(UNIT_VOLT,			"V");
	ADD(UNIT_COULOMB,		"C");
	ADD(UNIT_AMPEREHOUR,	"Ah");	

	return table;
}

UTEST(unit, print) {
	char buf[512] = {0};
	
	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_KILOGRAM);		EXPECT_STREQ("kg", buf);
	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_PASCAL);		EXPECT_STREQ("Pa", buf);
	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_NANOMETER);		EXPECT_STREQ("nm", buf);
	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_BAR);			EXPECT_STREQ("bar", buf);
//	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_M_S2);			EXPECT_STREQ("m/s^2", buf);
	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_KJ_MOL);		EXPECT_STREQ("kJ/mol", buf);
	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_NEWTON);		EXPECT_STREQ("N", buf);
	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_VOLT);			EXPECT_STREQ("V", buf);
	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_COULOMB);		EXPECT_STREQ("C", buf);
	md_unit_print(buf, sizeof(buf), (md_unit_t)UNIT_AMPEREHOUR);	EXPECT_STREQ("Ah", buf);
}

UTEST(unit, from_string) {
	
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_KILOGRAM,		md_unit_from_string(STR("kg"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_PASCAL,		md_unit_from_string(STR("Pa"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_NANOMETER,	md_unit_from_string(STR("nm"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_BAR,			md_unit_from_string(STR("bar"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_M_S2,			md_unit_from_string(STR("m/s^2"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_M_S2,			md_unit_from_string(STR("m*s^-2"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_KJ_MOL,		md_unit_from_string(STR("kJ/mol"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_NEWTON,		md_unit_from_string(STR("kg*m/s^2"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_NEWTON,		md_unit_from_string(STR("m kg/s^2"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_NEWTON,		md_unit_from_string(STR("N"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_VOLT,			md_unit_from_string(STR("V"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_COULOMB,		md_unit_from_string(STR("C"))));
	EXPECT_TRUE(md_unit_equal((md_unit_t)UNIT_AMPEREHOUR,	md_unit_from_string(STR("Ah"))));
	
}

UTEST(unit, convert) {
	// This is a bit clunky, but the interface is designed for converting larget batches of values
	double angle_deg = 1.0;
	md_unit_convert_inplace_d(&angle_deg, 1, md_unit_degree(), md_unit_radian());
	EXPECT_EQ(angle_deg, DEG_TO_RAD(1.0));
}
