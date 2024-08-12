#pragma once

#include <core/md_str.h>

#include <stdint.h>
#include <stdbool.h>

struct md_allocator_i;
struct md_strb_t;

/*

This is a very simplistic unit system to cover our needs.
There is a single unit structure which represents a fixed set of units with the corresponding
dimensions. Some of which are dimensionless, but are needed to be propageted through the system.
*/

// Base represents the data of the SI base units, excluding luminous intensity (candela) and with the addition of dimensionless units angle and count
// dim represents the exponent in each base unit
// Examples:
// Newton (N) (kg * m / s^2) would have base { mass=1, length=1 and time=-2 } and mult = 1.0
// Ångström (Å) (1.0e-9 m) would have base { length=1 } and mult = 1.0e-10
typedef union md_unit_base_t {
    struct {
        int8_t length;      // Meters
        int8_t mass;        // Kilogram
        int8_t time;        // Seconds
        int8_t current;     // Ampere
        int8_t mole;        // Mole            
        int8_t temp;        // Kelvin
        int8_t angle;       // Radians
        int8_t count;
    } dim;
    int8_t arr[8];
    uint64_t raw_bits;
} md_unit_base_t;

typedef struct md_unit_t {
    md_unit_base_t base;
    double         mult;  // Multiplier (e.g. kilo=1e3, nano=1e-9, Mega=1e6)
} md_unit_t;

#ifdef __cplusplus
extern "C" {
#endif

static inline bool md_unit_empty(md_unit_t unit) {
    return unit.base.raw_bits == 0 && unit.mult == 0;
}

static inline bool md_unit_unitless(md_unit_t unit) {
    return unit.base.raw_bits == 0 && unit.mult != 0;
}

static inline bool md_unit_base_equal(md_unit_t a, md_unit_t b) {
    return a.base.raw_bits == b.base.raw_bits;
}

bool md_unit_equal(md_unit_t a, md_unit_t b);

md_unit_t md_unit_mul(md_unit_t a, md_unit_t b);
md_unit_t md_unit_div(md_unit_t a, md_unit_t b);
// These are kind of pointless, but if the units differ, the result will be unitless
md_unit_t md_unit_add(md_unit_t a, md_unit_t b);
md_unit_t md_unit_sub(md_unit_t a, md_unit_t b);

md_unit_t md_unit_inv(md_unit_t unit);
md_unit_t md_unit_pow(md_unit_t unit, int pow);

// Scales the multiplier with supplied value
md_unit_t md_unit_scl(md_unit_t unit, double scl);

// Will convert values in its current unit into values in the new unit
bool md_unit_convert_inplace_d(double* values, int64_t num_values, md_unit_t cur_unit, md_unit_t new_unit);
bool md_unit_convert_inplace_f(float*  values, int64_t num_values, md_unit_t cur_unit, md_unit_t new_unit);

//void unit_convert_d(double* dst_values, const double* src_values, int64_t num_values, md_unit_t* dst_unit, md_unit_t src_unit, md_unit_base_t new_base);
//void unit_convert_f(float* dst_values, const float* src_values, int64_t num_values, md_unit_t* dst_unit, md_unit_t src_unit, md_unit_base_t new_base);

// Print a unit into a supplied buffer
size_t md_unit_print(char* buf, size_t cap, md_unit_t unit);

// Write unit to a string allocated by supplied allocator
str_t md_unit_to_string(md_unit_t unit, struct md_allocator_i* alloc);

// Parse unit from string
md_unit_t md_unit_from_string(str_t str);

// Some defined helper
md_unit_t md_unit_none(void);  // Unitless quantity
 
// Base units
md_unit_t md_unit_meter(void);
md_unit_t md_unit_nanometer(void);
md_unit_t md_unit_angstrom(void);

md_unit_t md_unit_kilogram(void);

md_unit_t md_unit_second(void);
md_unit_t md_unit_nanosecond(void);
md_unit_t md_unit_pikosecond(void);
md_unit_t md_unit_femtosecond(void);

md_unit_t md_unit_ampere(void);
md_unit_t md_unit_mole(void);
md_unit_t md_unit_kelvin(void);

// Not real base units, but required to represent some quantities
md_unit_t md_unit_radian(void);
md_unit_t md_unit_degree(void);
md_unit_t md_unit_count(void);

// Common units
md_unit_t md_unit_joule(void);
md_unit_t md_unit_pascal(void);
md_unit_t md_unit_bar(void);

#ifdef __cplusplus
}
#endif
