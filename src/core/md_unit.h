#pragma once

#include <core/md_str.h>

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

/*

A very simplistic unit system to cover our needs.

A unit is represented by a set of integer exponents over a fixed set of base dimensions
together with a scale factor which expresses the unit in terms of the SI base units:

    value_in_SI = value * md_unit_scale(unit)

Examples:
    Newton (N)      dim { mass = 1, length = 1, time = -2 },  scale = 1.0
    Angstrom (A)    dim { length = 1 },                       scale = 1.0e-10
    Degree (deg)    dim { angle = 1 },                        scale = PI / 180

A zero initialized md_unit_t is a valid unit which represents 'no unit' (dimensionless):
all dimensions are zero and a stored scale of zero is interpreted as 1.0.
This means that (md_unit_t){0} can be used as a sentinel, that it is the identity element
under multiplication and division, and that there is no invalid state to guard against.
Operations which may fail (parsing, conversion) report this through their return value.

*/

// Base represents the data of the SI base units, excluding luminous intensity (candela) and with the addition of dimensionless units angle and count
// dim represents the exponent in each base unit
typedef union md_unit_base_t {
    struct {
        int8_t length;      // Meters
        int8_t mass;        // Kilogram
        int8_t time;        // Seconds
        int8_t current;     // Ampere
        int8_t mole;        // Mole
        int8_t temp;        // Kelvin
        int8_t angle;       // Radians
        int8_t count;       // Dimensionless count
    } dim;
    int8_t   arr[8];
    uint64_t raw_bits;
} md_unit_base_t;

typedef struct md_unit_t {
    md_unit_base_t base;
    double         mult;  // Scale relative to the SI base units. Zero is interpreted as 1.0, use md_unit_scale() to read it.
} md_unit_t;

#ifdef __cplusplus
extern "C" {
#endif

// Scale of the unit expressed in SI base units. A stored multiplier of zero is interpreted as 1.0
// such that a zero initialized unit is dimensionless with a scale of one.
static inline double md_unit_scale(md_unit_t unit) {
    return unit.mult == 0.0 ? 1.0 : unit.mult;
}

// True if the unit carries no dimensions and therefore has no representation.
// Any scale it may carry belongs to the value itself (e.g. the 1e3 of km/m).
static inline bool md_unit_is_none(md_unit_t unit) {
    return unit.base.raw_bits == 0;
}

// True if the units share the same dimensions and are therefore convertible.
static inline bool md_unit_base_equal(md_unit_t a, md_unit_t b) {
    return a.base.raw_bits == b.base.raw_bits;
}

// Equal dimensions and equal scale (within a small relative epsilon).
bool md_unit_equal(md_unit_t a, md_unit_t b);

md_unit_t md_unit_mul(md_unit_t a, md_unit_t b);
md_unit_t md_unit_div(md_unit_t a, md_unit_t b);

// Addition and subtraction only preserve the unit if both operands agree, otherwise the result is 'none'
md_unit_t md_unit_add(md_unit_t a, md_unit_t b);
md_unit_t md_unit_sub(md_unit_t a, md_unit_t b);

md_unit_t md_unit_inv(md_unit_t unit);
md_unit_t md_unit_pow(md_unit_t unit, int exponent);

// Scales the unit with the supplied factor (e.g. md_unit_scl(md_unit_meter(), 1e-3) -> millimeter)
md_unit_t md_unit_scl(md_unit_t unit, double scl);

// Computes the factor which converts values expressed in from_unit into to_unit.
// Returns false if the units do not share the same dimensions, in which case no conversion is possible.
bool md_unit_conversion_factor(double* factor, md_unit_t from_unit, md_unit_t to_unit);

// Print a unit into the supplied buffer, returns the length written (excluding null terminator).
// The result is a valid input to md_unit_parse(). A unit which is 'none' prints nothing.
size_t md_unit_print(char* buf, size_t cap, md_unit_t unit);

// Parse a unit from a string, the unit is set to 'none' if it could not be parsed.
// An empty string is a successfully parsed 'none'.
bool md_unit_parse(md_unit_t* unit, str_t str);

// Some defined helpers

md_unit_t md_unit_none(void);   // No unit (dimensionless), equivalent to (md_unit_t){0}

// Base units

// Length units
md_unit_t md_unit_meter(void);
md_unit_t md_unit_nanometer(void);
md_unit_t md_unit_angstrom(void);
md_unit_t md_unit_bohr_radius(void);

// Mass units
md_unit_t md_unit_kilogram(void);

// Time units
md_unit_t md_unit_second(void);
md_unit_t md_unit_nanosecond(void);
md_unit_t md_unit_picosecond(void);
md_unit_t md_unit_femtosecond(void);

// Electric current units
md_unit_t md_unit_ampere(void);

md_unit_t md_unit_mole(void);

// Temperature units
md_unit_t md_unit_kelvin(void);

// Not real base units, but required to represent some quantities
md_unit_t md_unit_radian(void);
md_unit_t md_unit_degree(void);
md_unit_t md_unit_count(void);

// Common units
md_unit_t md_unit_joule(void);
md_unit_t md_unit_electronvolt(void);
md_unit_t md_unit_hertz(void);
md_unit_t md_unit_pascal(void);
md_unit_t md_unit_bar(void);

#ifdef __cplusplus
}
#endif
