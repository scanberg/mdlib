#include "utest.h"

#include <core/md_unit.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_str.h>

typedef struct {
	md_unit_t unit;
	str_t     str;
} item_t;

#define ADD(unit, cstr) md_array_push(table, ((item_t){unit, MAKE_STR(cstr)}), alloc)

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
	
	unit_print(buf, sizeof(buf), (md_unit_t)UNIT_KILOGRAM);		EXPECT_STREQ("kg", buf);
	unit_print(buf, sizeof(buf), (md_unit_t)UNIT_PASCAL);		EXPECT_STREQ("Pa", buf);
	unit_print(buf, sizeof(buf), (md_unit_t)UNIT_NANOMETER);	EXPECT_STREQ("nm", buf);
	unit_print(buf, sizeof(buf), (md_unit_t)UNIT_BAR);			EXPECT_STREQ("bar", buf);
	unit_print(buf, sizeof(buf), (md_unit_t)UNIT_M_S2);			EXPECT_STREQ("m/s^2", buf);
	//unit_print(buf, sizeof(buf), (md_unit_t)UNIT_KJ_MOL);		EXPECT_STREQ("kJ/mol", buf);
	unit_print(buf, sizeof(buf), (md_unit_t)UNIT_VOLT);			EXPECT_STREQ("V", buf);
	unit_print(buf, sizeof(buf), (md_unit_t)UNIT_COULOMB);		EXPECT_STREQ("C", buf);
	unit_print(buf, sizeof(buf), (md_unit_t)UNIT_AMPEREHOUR);	EXPECT_STREQ("Ah", buf);
}

UTEST(unit, from_string) {
	/*
	EXPECT_TRUE(unit_equal((md_unit_t)UNIT_KILOGRAM,	unit_from_string(MAKE_STR("Kg"))));
	EXPECT_TRUE(unit_equal((md_unit_t)UNIT_PASCAL,		unit_from_string(MAKE_STR("Pa"))));
	EXPECT_TRUE(unit_equal((md_unit_t)UNIT_NANOMETER,	unit_from_string(MAKE_STR("nm"))));
	EXPECT_TRUE(unit_equal((md_unit_t)UNIT_BAR,			unit_from_string(MAKE_STR("bar"))));
	EXPECT_TRUE(unit_equal((md_unit_t)UNIT_M_S2,		unit_from_string(MAKE_STR("m/s^2"))));
	EXPECT_TRUE(unit_equal((md_unit_t)UNIT_KJ_MOL,		unit_from_string(MAKE_STR("kJ/mol"))));
	EXPECT_TRUE(unit_equal((md_unit_t)UNIT_VOLT,		unit_from_string(MAKE_STR("V"))));
	EXPECT_TRUE(unit_equal((md_unit_t)UNIT_COULOMB,		unit_from_string(MAKE_STR("C"))));
	EXPECT_TRUE(unit_equal((md_unit_t)UNIT_AMPEREHOUR,	unit_from_string(MAKE_STR("Ah"))));
	*/
}

UTEST(unit, convert) {
	// This is a bit clunky, but the interface is designed for converting larget batches of values
	double angle_deg = 1.0;
	unit_convert_d_inplace(&angle_deg, 1, unit_degree(), unit_radian());
	EXPECT_EQ(angle_deg, DEG_TO_RAD(1.0));
}
