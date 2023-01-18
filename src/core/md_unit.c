#include "md_unit.h"

#include "md_log.h"
#include "md_common.h"
#include "md_str.h"

#include <stdio.h>
#include <math.h>
#include <string.h>

#define S(str) {str"", sizeof(str)-1}

typedef struct {
    float value;
    str_t prefix;
} si_prefix_t;


static const si_prefix_t si_prefixes[] = {
    {1e-3F,  S(u8"m")},  {1.0F / 1e3F,   S(u8"m")},
    {1e3F,   S(u8"k")},  {1.0F / 1e-3F,  S(u8"k")},
    {1e-6F,  S(u8"μ")},  {1.0F / 1e6F,   S(u8"μ")},
    {1e-6F,  S(u8"u")},  {1.0F / 1e6F,   S(u8"u")},
    {1e-2F,  S(u8"c")},  {1.0F / 1e2F,   S(u8"c")},
    {1e6F,   S(u8"M")},  {1.0F / 1e-6F,  S(u8"M")},
    {1e9F,   S(u8"G")},  {1.0F / 1e-9F,  S(u8"G")},
    {1e-9F,  S(u8"n")},  {1.0F / 1e9F,   S(u8"n")},
    {1e-12F, S(u8"p")},  {1.0F / 1e12F,  S(u8"p")},
    {1e-15F, S(u8"f")},  {1.0F / 1e15F,  S(u8"f")},
    {1e-18F, S(u8"a")},  {1.0F / 1e18F,  S(u8"a")},
    {1e-21F, S(u8"z")},  {1.0F / 1e21F,  S(u8"z")},
    {1e-24F, S(u8"y")},  {1.0F / 1e24F,  S(u8"y")},
    {1e12F,  S(u8"T")},  {1.0F / 1e-12F, S(u8"T")},
    {1e15F,  S(u8"P")},  {1.0F / 1e-15F, S(u8"P")},
    {1e18F,  S(u8"E")},  {1.0F / 1e-18F, S(u8"E")},
    {1e21F,  S(u8"Z")},  {1.0F / 1e-21F, S(u8"Z")},
    {1e24F,  S(u8"Y")},  {1.0F / 1e-24F, S(u8"Y")},
};

static str_t find_prefix_str_from_value(float value) {
    for (int i = 0; i < (int)ARRAY_SIZE(si_prefixes); ++i) {
        if (si_prefixes[i].value == value) return si_prefixes[i].prefix;
    }
    return (str_t){0};
}

static float find_prefix_value_from_str(str_t str) {
    for (int i = 0; i < (int)ARRAY_SIZE(si_prefixes); ++i) {
        if (str_equal(str, si_prefixes[i].prefix)) return si_prefixes[i].value;
    }
    return 0;
}

typedef struct {
    md_unit_t   unit;
    str_t       str;
} unit_name_t;

#define UNIT_NONE           {.mult = 1.0}

// Base
#define UNIT_METER          {.base = {.dim = {.length = 1,}},   .mult = 1.0}
#define UNIT_KILOGRAM       {.base = {.dim = {.mass = 1,}},     .mult = 1.0}
#define UNIT_SECOND         {.base = {.dim = {.time = 1,}},     .mult = 1.0}
#define UNIT_AMPERE         {.base = {.dim = {.current = 1,}},  .mult = 1.0}
#define UNIT_MOLE           {.base = {.dim = {.mole = 1,}},     .mult = 1.0}
#define UNIT_KELVIN         {.base = {.dim = {.temp = 1,}},     .mult = 1.0}
#define UNIT_RADIAN         {.base = {.dim = {.angle = 1,}},    .mult = 1.0}
#define UNIT_COUNT          {.base = {.dim = {.count = 1,}},    .mult = 1.0}

// Extra simple
#define UNIT_MILLIMETER     {.base = {.dim = {.length = 1,}},   .mult = 1e-3}
#define UNIT_MICROMETER     {.base = {.dim = {.length = 1,}},   .mult = 1e-6}
#define UNIT_NANOMETER      {.base = {.dim = {.length = 1,}},   .mult = 1e-9}
#define UNIT_ANGSTROM       {.base = {.dim = {.length = 1,}},   .mult = 1e-10}
#define UNIT_MILLISECOND    {.base = {.dim = {.time = 1,}},     .mult = 1e-3}
#define UNIT_NANOSECOND     {.base = {.dim = {.time = 1,}},     .mult = 1e-9}
#define UNIT_PIKOSECOND     {.base = {.dim = {.time = 1,}},     .mult = 1e-12}

#define UNIT_DEGREE         {.base = {.dim = {.angle = 1,}},    .mult = RAD_TO_DEG(1),}

#define UNIT_METER_SQ       {.base = {.dim = {.length = 2,}},   .mult = 1.0}

// Composite
#define UNIT_PASCAL         {.base = {.dim = {.mass = 1, .length = -1, .time = -2,}}, .mult = 1.0}
#define UNIT_BAR            {.base = {.dim = {.mass = 1, .length = -1, .time = -2,}}, .mult = 1.0e5}
#define UNIT_JOULE          {.base = {.dim = {.mass = 1, .length =  2, .time = -2,}}, .mult = 1.0}
#define UNIT_NEWTON         {.base = {.dim = {.mass = 1, .length =  1, .time = -2,}}, .mult = 1.0}
//#define UNIT_NEWTONMETER    {.base = {.dim = {.mass = 1, .length =  2, .time = -2,}}, .mult = 1.0}

#define UNIT_VOLT           {.base = {.dim = {.mass = 1, .length =  2, .time = -3, .current = -1}}, .mult = 1.0}
#define UNIT_WATT           {.base = {.dim = {.mass = 1, .length =  2, .time = -3}}, .mult = 1.0}
#define UNIT_COULOMB        {.base = {.dim = {.time = 1, .current = 1}}, .mult = 1.0}
#define UNIT_AMPEREHOUR     {.base = {.dim = {.time = 1, .current = 1}}, .mult = 3600.0}

static const unit_name_t predefined_units[] = {
    // Base
    {UNIT_METER,        S(u8"m")},
    {UNIT_KILOGRAM,     S(u8"kg")},
    {UNIT_SECOND,       S(u8"s")},
    {UNIT_AMPERE,       S(u8"A")},
    {UNIT_MOLE,         S(u8"mol")},
    {UNIT_KELVIN,       S(u8"K")},
    {UNIT_RADIAN,       S(u8"rad")},
    {UNIT_COUNT,        S(u8"count")},

    // Extra extra!
    {UNIT_MILLIMETER,   S(u8"mm")},
    {UNIT_MICROMETER,   S(u8"μm")},
    {UNIT_METER_SQ,     S(u8"m²")},
    {UNIT_METER_SQ,     S(u8"m^2")},
    {UNIT_NANOSECOND,   S(u8"ns")},
    {UNIT_PIKOSECOND,   S(u8"ps")},
    {UNIT_NANOMETER,    S(u8"nm")},
    {UNIT_ANGSTROM,     S(u8"Å")},
    {UNIT_DEGREE,       S(u8"deg")},
    {UNIT_DEGREE,       S(u8"\u00B0")},
    {UNIT_PASCAL,       S(u8"Pa")},
    {UNIT_BAR,          S(u8"bar")},
    {UNIT_JOULE,        S(u8"J")},
    {UNIT_NEWTON,       S(u8"N")},
    {UNIT_VOLT,         S(u8"V")},
    {UNIT_WATT,         S(u8"W")},
    {UNIT_COULOMB,      S(u8"C")},
    {UNIT_AMPEREHOUR,   S(u8"Ah")},
};

static md_unit_t find_unit_from_predefined(str_t name) {
    for (int i = 0; i < (int)ARRAY_SIZE(predefined_units); ++i) {
        if (str_equal(name, predefined_units[i].str)) {
            return predefined_units[i].unit;
        }
    }
    return (md_unit_t){0};
}

static str_t find_name_from_predefined(md_unit_t unit) {
    for (int i = 0; i < (int)ARRAY_SIZE(predefined_units); ++i) {
        if (unit_equal(predefined_units[i].unit, unit)) {
            return predefined_units[i].str;
        }
    }
    return (str_t){0};
}

// b is allowed to match or be a subset of a (with the same sign), but not vice-versa
static inline int compute_match_score(md_unit_t a, md_unit_t b) {
    int score = 0;
    for (int i = 0; i < (int)ARRAY_SIZE(b.base.arr); ++i) {
        
        int sign_a = (a.base.arr[i] > 0) - (a.base.arr[i] < 0);
        int sign_b = (b.base.arr[i] > 0) - (b.base.arr[i] < 0);

        if (sign_b) {
            score += (sign_b == sign_a);
        }
    }
    return score;
}

static md_unit_t find_best_matching_predefined(md_unit_t unit) {
    int best_match_score = 0;
    int best_match_idx   = -1;
    for (int i = 0; i < (int)ARRAY_SIZE(predefined_units); ++i) {
        int match_score = compute_match_score(unit, predefined_units[i].unit);
        if (match_score > best_match_score) {
            best_match_score = match_score;
            best_match_idx = i;
        }
    }
    if (best_match_idx != -1) {
        return predefined_units[best_match_idx].unit;
    }
    return (md_unit_t){0};
}

static int print_exponent(char* buf, int cap, int exp) {
    switch (exp) {
    case 0: return 0;
    case 1: return 0;
    //case 2: return snprintf(buf, cap, u8"²");
    //case 3: return snprintf(buf, cap, u8"³");
    //case -1: return snprintf(buf, cap, u8"⁻¹");
    //case -2: return snprintf(buf, cap, u8"⁻²");
    //case -3: return snprintf(buf, cap, u8"⁻³");
    default: return snprintf(buf, cap, "^%i",exp);
    }
}


bool unit_equal(md_unit_t a, md_unit_t b) {
    if (!unit_base_equal(a, b)) return false;

    // Compare multiplier
    if (a.mult == b.mult) return true;

    // Compare magnitude within reasonable threshold
    double d = a.mult - b.mult;
    if (d == 0 || fpclassify(d) == FP_SUBNORMAL) {
        return true;
    }

    // @TODO: This is not really correct
    // Improve this to a more robust comparison
    return fabs(d) < 5e-13;
}

md_unit_t unit_mul(md_unit_t a, md_unit_t b) {
    md_unit_t result = {
        .base = {
            .dim.length     = a.base.dim.length  + b.base.dim.length,
            .dim.mass       = a.base.dim.mass    + b.base.dim.mass,
            .dim.time       = a.base.dim.time    + b.base.dim.time,
            .dim.current    = a.base.dim.current + b.base.dim.current,
            .dim.mole       = a.base.dim.mole    + b.base.dim.mole,
            .dim.temp       = a.base.dim.temp    + b.base.dim.temp,
            .dim.angle      = a.base.dim.angle   + b.base.dim.angle,
            .dim.count      = a.base.dim.count   + b.base.dim.count,
    },
    .mult = a.mult * b.mult,
    };

    return result;
}

md_unit_t unit_div(md_unit_t a, md_unit_t b) {
    md_unit_t result = {
        .base = {
            .dim.length     = a.base.dim.length  - b.base.dim.length,
            .dim.mass       = a.base.dim.mass    - b.base.dim.mass,
            .dim.time       = a.base.dim.time    - b.base.dim.time,
            .dim.current    = a.base.dim.current - b.base.dim.current,
            .dim.mole       = a.base.dim.mole    - b.base.dim.mole,
            .dim.temp       = a.base.dim.temp    - b.base.dim.temp,
            .dim.angle      = a.base.dim.angle   - b.base.dim.angle,
            .dim.count      = a.base.dim.count   - b.base.dim.count,
    },
    .mult = a.mult / b.mult,
    };

    return result;
}

md_unit_t unit_add(md_unit_t a, md_unit_t b) {
    md_unit_t result = {0};
    if (a.base.raw_bits == b.base.raw_bits) {
        result = a;
    }
    return result;
}

md_unit_t unit_sub(md_unit_t a, md_unit_t b) {
    md_unit_t result = {0};
    if (a.base.raw_bits == b.base.raw_bits) {
        result = a;
    }
    return result;
}

md_unit_t unit_inv(md_unit_t unit) {
    md_unit_t inv = {
        .base = {
            .dim = {
                .length     = -unit.base.dim.length,
                .mass       = -unit.base.dim.mass,
                .time       = -unit.base.dim.time,
                .current    = -unit.base.dim.current,
                .mole       = -unit.base.dim.mole,
                .temp       = -unit.base.dim.temp,
                .angle      = -unit.base.dim.angle,
                .count      = -unit.base.dim.count,
            },
        },
    };
    return inv;
}

md_unit_t unit_pow(md_unit_t unit, int pow) {
    md_unit_t unit_pow = {
        .base = {
            .dim = {
                .length     = (int8_t)(unit.base.dim.length  * pow),
                .mass       = (int8_t)(unit.base.dim.mass    * pow),
                .time       = (int8_t)(unit.base.dim.time    * pow),
                .current    = (int8_t)(unit.base.dim.current * pow),
                .mole       = (int8_t)(unit.base.dim.mole    * pow),
                .temp       = (int8_t)(unit.base.dim.temp    * pow),
                .angle      = (int8_t)(unit.base.dim.angle   * pow),
                .count      = (int8_t)(unit.base.dim.count   * pow),
            },
        },
        .mult = unit.mult,
    };
    return unit_pow;
}

md_unit_t unit_scl(md_unit_t unit, double scl) {
    md_unit_t unit_scl = unit;
    unit_scl.mult *= scl;
    return unit_scl;
}

bool unit_convert_inplace_d(double* values, int64_t num_values, md_unit_t cur_unit, md_unit_t new_unit) {
    if (!unit_base_equal(cur_unit, new_unit)) {
        md_log(MD_LOG_TYPE_ERROR, "Failed to perform unit conversion, the current unit cannot be converted into new unit");
        return false;
    }

    if (unit_equal(cur_unit, new_unit)) {
        // Easy peasy, do nothing
        return true;
    }

    if (cur_unit.mult == 0) {
        md_log(MD_LOG_TYPE_INFO, "Failed to perform unit conversion, multiplier in current unit was zero");
        return false;
    }
    if (new_unit.mult == 0) {
        md_log(MD_LOG_TYPE_INFO, "Failed to perform unit conversion, multiplier in new unit was zero");
        return false;
    }

    // unit base is equal, compute conversion factor from multipliers
    double mult_factor = new_unit.mult / cur_unit.mult;
    for (int64_t i = 0; i < num_values; ++i) {
        values[i] *= mult_factor;
    }

    return true;
}

bool unit_convert_inplace_f(float* values, int64_t num_values, md_unit_t cur_unit, md_unit_t new_unit) {
    if (!unit_base_equal(cur_unit, new_unit)) {
        md_log(MD_LOG_TYPE_ERROR, "Failed to perform unit conversion, the current unit cannot be converted into new unit");
        return false;
    }

    if (unit_equal(cur_unit, new_unit)) {
        // Easy peasy, do nothing
        return true;
    }

    if (cur_unit.mult == 0) {
        md_log(MD_LOG_TYPE_INFO, "Failed to perform unit conversion, multiplier in current unit was zero");
        return false;
    }
    if (new_unit.mult == 0) {
        md_log(MD_LOG_TYPE_INFO, "Failed to perform unit conversion, multiplier in new unit was zero");
        return false;
    }

    // unit base is equal, compute conversion factor from multipliers
    double mult_factor = new_unit.mult / cur_unit.mult;
    for (int64_t i = 0; i < num_values; ++i) {
        values[i] = (float)(mult_factor * values[i]);
    }

    return true;
}


#define PRINT(...) len += snprintf(buf + len, MAX(0, cap - len), ##__VA_ARGS__)

md_unit_t mask_pos(md_unit_t unit) {
    md_unit_t unit_pos = {
        .base = {
            .dim = {
                .length  = unit.base.dim.length  > 0 ? unit.base.dim.length  : 0,
                .mass    = unit.base.dim.mass    > 0 ? unit.base.dim.mass    : 0,
                .time    = unit.base.dim.time    > 0 ? unit.base.dim.time    : 0,
                .current = unit.base.dim.current > 0 ? unit.base.dim.current : 0,
                .mole    = unit.base.dim.mole    > 0 ? unit.base.dim.mole    : 0,
                .temp    = unit.base.dim.temp    > 0 ? unit.base.dim.temp    : 0,
                .angle   = unit.base.dim.angle   > 0 ? unit.base.dim.angle   : 0,
                .count   = unit.base.dim.count   > 0 ? unit.base.dim.count   : 0,
            },
        },
        .mult = unit.mult,
    };
    return unit_pos;
}

md_unit_t mask_neg(md_unit_t unit) {
    md_unit_t unit_neg = {
        .base = {
            .dim = {
                .length  = unit.base.dim.length  < 0 ? unit.base.dim.length  : 0,
                .mass    = unit.base.dim.mass    < 0 ? unit.base.dim.mass    : 0,
                .time    = unit.base.dim.time    < 0 ? unit.base.dim.time    : 0,
                .current = unit.base.dim.current < 0 ? unit.base.dim.current : 0,
                .mole    = unit.base.dim.mole    < 0 ? unit.base.dim.mole    : 0,
                .temp    = unit.base.dim.temp    < 0 ? unit.base.dim.temp    : 0,
                .angle   = unit.base.dim.angle   < 0 ? unit.base.dim.angle   : 0,
                .count   = unit.base.dim.count   < 0 ? unit.base.dim.count   : 0,
            },
        },
        .mult = unit.mult,
    };
    return unit_neg;
}

// A bit wierd name, counts the number of base unit dimensions which is not zero
uint8_t base_count(md_unit_t unit) {
    uint8_t count = 0;
    count += unit.base.dim.length  != 0;
    count += unit.base.dim.mass    != 0;
    count += unit.base.dim.time    != 0;
    count += unit.base.dim.current != 0;
    count += unit.base.dim.mole    != 0;
    count += unit.base.dim.temp    != 0;
    count += unit.base.dim.angle   != 0;
    count += unit.base.dim.count   != 0;
    return count;
}

int internal_print_dims_SI(char* buf, int cap, md_unit_t unit) {
    int len = 0;
    if (unit.base.dim.length != 0) {
        PRINT("m");
        len += print_exponent(buf + len, cap - len, unit.base.dim.length);
    }
    if (unit.base.dim.mass != 0) {
        PRINT("kg");
        len += print_exponent(buf + len, cap - len, unit.base.dim.mass);
    }
    if (unit.base.dim.time != 0) {
        PRINT("s");
        len += print_exponent(buf + len, cap - len, unit.base.dim.time);
    }
    if (unit.base.dim.current != 0) {
        PRINT("A");
        len += print_exponent(buf + len, cap - len, unit.base.dim.current);
    }
    if (unit.base.dim.mole != 0) {
        PRINT("mol");
        len += print_exponent(buf + len, cap - len, unit.base.dim.mole);
    }
    if (unit.base.dim.temp != 0) {
        PRINT("K");
        len += print_exponent(buf + len, cap - len, unit.base.dim.temp);
    }
    if (unit.base.dim.angle != 0) {
        PRINT("rad");
        len += print_exponent(buf + len, cap - len, unit.base.dim.angle);
    }
    if (unit.base.dim.count != 0) {
        PRINT("count");
        len += print_exponent(buf + len, cap - len, unit.base.dim.count);
    }
    return len;
}

int internal_print(char* buf, int cap, md_unit_t unit, int depth) {
    int len = 0;

    if (unit_empty(unit) || unit_unitless(unit)) {
        return 0;
    }

    if (unit.mult == 0) {
        PRINT("0*");
        unit.mult = 1;
    }

    // Try to match against predefined units (exact)
    str_t exact = find_name_from_predefined(unit);
    if (!str_empty(exact)) {
        PRINT("%.*s", (int)exact.len, exact.ptr);
        return len;
    }

    // Try to match inverted against predefined units in table
    {
        str_t str = find_name_from_predefined(unit_inv(unit));
        if (!str_empty(str)) {
            PRINT("1/%.*s", (int)str.len, str.ptr);
            return len;
        }
    }

    // Try to match find an exact matching base and see if we can apply a prefix
    for (int i = 0; i < (int)ARRAY_SIZE(predefined_units); ++i) {
        if (unit_base_equal(unit, predefined_units[i].unit)) {
            // We found a matching base, we just need to work out a matching prefix scaling
            str_t str  = predefined_units[i].str;
            double mult  = unit.mult / predefined_units[i].unit.mult;
            str_t prefix = STR("");

            if (mult != 1.0) {
                prefix = find_prefix_str_from_value((float)mult);
                if (str_empty(prefix)) {
                    md_log(MD_LOG_TYPE_DEBUG, "Unable to find suitable prefix?");
                    return 0;
                }
            }

            PRINT("%.*s%.*s", (int)prefix.len, prefix.ptr, (int)str.len, str.ptr);
            return len;
        }
        // Test for inverted
        if (unit_base_equal(unit_inv(unit), predefined_units[i].unit)) {
            // We found a matching unit, we just need to work out a matching prefix scaling
            str_t str  = predefined_units[i].str;
            double mult  = unit.mult / predefined_units[i].unit.mult;
            str_t prefix = STR("");

            if (mult != 1.0) {
                prefix = find_prefix_str_from_value((float)mult);
                if (str_empty(prefix)) {
                    md_log(MD_LOG_TYPE_DEBUG, "Unable to find suitable prefix?");
                    return 0;
                }
            }

            PRINT("1/%.*s%.*s", (int)prefix.len, prefix.ptr, (int)str.len, str.ptr);
            return len;
        }
    }

    // Try to find best matching base and see where we can go from there
    {
        md_unit_t best_match = find_best_matching_predefined(unit);
        if (!unit_empty(best_match)) {
            double mult  = unit.mult / best_match.mult;
            str_t prefix = STR("");

            if (mult != 1.0) {
                prefix = find_prefix_str_from_value((float)mult);
                if (str_empty(prefix)) {
                    // Error, could not determine prefix
                    md_log(MD_LOG_TYPE_DEBUG, "Could not determine prefix");
                    return 0;
                }
            }

            if (unit.base.raw_bits == best_match.base.raw_bits) {
                // Easy peasy, exact match with optional prefix
                
                PRINT("%.*s", (int)prefix.len, prefix.ptr);
            }
        }
    }

    // Last resort, just print out whatever there is in metric units
    {
        str_t prefix = find_prefix_str_from_value((float)unit.mult);
        if (!str_empty(prefix)) {
            PRINT("%.*s", (int)prefix.len, prefix.ptr);
        }
        else if (unit.mult != 1.0) {
            PRINT("%f*", unit.mult);
        }
        unit.mult = 1.0;
        len += internal_print_dims_SI(buf + len, cap - len, mask_pos(unit));
        if (base_count(mask_neg(unit))) {
            PRINT("/");
            len += internal_print_dims_SI(buf + len, cap - len, unit_inv(mask_neg(unit)));
        }   
    }

    return len;
}

int unit_print(char* buf, int cap, md_unit_t unit) {
    return internal_print(buf, cap, unit, 0);
}

str_t unit_to_string(md_unit_t unit, struct md_allocator_i* alloc) {
    str_t str = {0};

    char buf[128];
    int len = unit_print(buf, sizeof(buf), unit);
    if (len > 0) {
        str = str_copy((str_t){buf, len}, alloc);
    }

    return str;
}

md_unit_t unit_from_string(str_t str) {
    str = str_trim(str);
    if (str_empty(str)) {
        return (md_unit_t) {0};
    }

    int64_t delim;
    
    delim = str_find_char(str, '/');
    if (delim != -1) {
        str_t num = str_substr(str, 0, delim);
        str_t den = str_substr(str, delim + 1, str.len);
        return unit_div(unit_from_string(num), unit_from_string(den));
    }

    delim = str_find_char(str, '*');
    if (delim != -1) {
        str_t p1 = str_substr(str, 0, delim);
        str_t p2 = str_substr(str, delim + 1, str.len);
        return unit_mul(unit_from_string(p1), unit_from_string(p2));
    }
    
    delim = str_find_char(str, ' ');
    if (delim != -1) {
        str_t p1 = str_substr(str, 0, delim);
        str_t p2 = str_substr(str, delim + 1, str.len);
        return unit_mul(unit_from_string(p1), unit_from_string(p2));
    }
    
    // Try to match against defined units (exact)
    md_unit_t exact = find_unit_from_predefined(str);
    if (!unit_empty(exact)) {
        return exact;
    }
        
    // Try to match against longest matching predefined unit (from right to left)
    delim = str_find_char(str, '^');
    str_t base = str;
    str_t exp = {0};
    str_t rest = {0};
    if (delim != -1) {
        base = str_substr(str, 0, delim);
        int64_t i = delim + 1;
        while (i < str.len && (str.ptr[i] == '-' || str.ptr[i] == '+' || is_alpha(str.ptr[i]))) ++i;
        exp = str_substr(str, delim + 1, i);
        rest = str_substr(str, i + 1, str.len);
    }
        
    int best_match_idx = -1;
    for (int i = 0; i < (int)ARRAY_SIZE(predefined_units); ++i) {
        if (str_ends_with(base, predefined_units[i].str)) {
            if (best_match_idx == -1 || predefined_units[i].str.len > predefined_units[best_match_idx].str.len) {
                best_match_idx = i;
            }
        }
    }
    
    if (best_match_idx != -1) {
        str_t prefix = str_substr(base, 0, base.len - predefined_units[best_match_idx].str.len);
        md_unit_t unit = predefined_units[best_match_idx].unit;
        if (!str_empty(prefix)) {
            unit.mult = unit.mult * find_prefix_value_from_str(prefix);
        }
            
        if (!str_empty(exp)) {
            int pow;
            if (str_extract_i32(&pow, &exp)) {
                unit = unit_pow(unit, pow);
            }
        }

        if (!str_empty(rest)) {
            unit = unit_mul(unit, unit_from_string(rest));
        }

        return unit;
    }

    md_logf(MD_LOG_TYPE_ERROR, "Could not convert string to unit: '%.*s'", (int)str.len, str.ptr);
    return (md_unit_t) {0};
}

md_unit_t unit_none() {
    return (md_unit_t)UNIT_NONE;
}

md_unit_t unit_meter() {
    return (md_unit_t)UNIT_METER;
}

md_unit_t unit_nanometer() {
    return (md_unit_t)UNIT_NANOMETER;
}

md_unit_t unit_angstrom() {
    return (md_unit_t)UNIT_ANGSTROM;
}

md_unit_t unit_kilogram() {
    return (md_unit_t)UNIT_KILOGRAM;
}

md_unit_t unit_second() {
    return (md_unit_t)UNIT_SECOND;
}

md_unit_t unit_nanosecond() {
    return (md_unit_t)UNIT_NANOSECOND;
}

md_unit_t unit_pikosecond() {
    return (md_unit_t)UNIT_PIKOSECOND;
}

md_unit_t unit_ampere() {
    return (md_unit_t)UNIT_AMPERE;
}

md_unit_t unit_mole() {
    return (md_unit_t)UNIT_MOLE;
}

md_unit_t unit_kelvin() {
    return (md_unit_t)UNIT_KELVIN;
}

md_unit_t unit_radian() {
    return (md_unit_t)UNIT_RADIAN;
}

md_unit_t unit_degree() {
    return (md_unit_t)UNIT_DEGREE;
}

md_unit_t unit_count() {
    return (md_unit_t)UNIT_COUNT;
}

md_unit_t unit_joule() {
    return (md_unit_t)UNIT_JOULE;
}

md_unit_t unit_pascal() {
    return (md_unit_t)UNIT_PASCAL;
}

md_unit_t unit_bar() {
    return (md_unit_t)UNIT_BAR;
}
