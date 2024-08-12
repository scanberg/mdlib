#include <core/md_unit.h>

#include <core/md_log.h>
#include <core/md_common.h>
#include <core/md_str.h>

#include <math.h>  // fpclassify
#include <stdio.h> // snprintf

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
        if (str_eq(str, si_prefixes[i].prefix)) return si_prefixes[i].value;
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
#define UNIT_FEMTOSECOND	{.base = {.dim = {.time = 1,}},     .mult = 1e-15}

#define UNIT_DEGREE         {.base = {.dim = {.angle = 1,}},    .mult = RAD_TO_DEG(1),}

#define UNIT_METER_SQ       {.base = {.dim = {.length = 2,}},   .mult = 1.0}

// Composite
#define UNIT_PASCAL         {.base = {.dim = {.mass = 1, .length = -1, .time = -2,}}, .mult = 1.0}
#define UNIT_BAR            {.base = {.dim = {.mass = 1, .length = -1, .time = -2,}}, .mult = 1.0e5}
#define UNIT_JOULE          {.base = {.dim = {.mass = 1, .length =  2, .time = -2,}}, .mult = 1.0}
#define UNIT_NEWTON         {.base = {.dim = {.mass = 1, .length =  1, .time = -2,}}, .mult = 1.0}
#define UNIT_NEWTONMETER    {.base = {.dim = {.mass = 1, .length =  2, .time = -2,}}, .mult = 1.0}

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
    {UNIT_NEWTONMETER,  S(u8"Nm")},
    {UNIT_VOLT,         S(u8"V")},
    {UNIT_WATT,         S(u8"W")},
    {UNIT_COULOMB,      S(u8"C")},
    {UNIT_AMPEREHOUR,   S(u8"Ah")},
};

static md_unit_t find_unit_from_predefined(str_t name) {
    for (int i = 0; i < (int)ARRAY_SIZE(predefined_units); ++i) {
        if (str_eq(name, predefined_units[i].str)) {
            return predefined_units[i].unit;
        }
    }
    return (md_unit_t){0};
}

static str_t find_name_from_predefined(md_unit_t unit) {
    for (int i = 0; i < (int)ARRAY_SIZE(predefined_units); ++i) {
        if (md_unit_equal(predefined_units[i].unit, unit)) {
            return predefined_units[i].str;
        }
    }
    return (str_t){0};
}

static inline int compute_match_score(md_unit_t unit, md_unit_t ref) {
    int score = 0;
    for (int i = 0; i < (int)ARRAY_SIZE(ref.base.arr); ++i) {
        if (ref.base.arr[i] != 0 && unit.base.arr[i] == ref.base.arr[i]) {
            score += 10;
        }
    }
    return score;
}

static int find_best_matching_predefined(md_unit_t unit) {
    int best_match_score = 0;
    int best_match_idx   = -1;
    for (int i = 0; i < (int)ARRAY_SIZE(predefined_units); ++i) {
        int match_score = compute_match_score(unit, predefined_units[i].unit);
        if (match_score > best_match_score) {
            best_match_score = match_score;
            best_match_idx = i;
        }
    }
    return best_match_idx;
}

static size_t print_exponent(char* buf, size_t cap, int exp) {
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


bool md_unit_equal(md_unit_t a, md_unit_t b) {
    if (!md_unit_base_equal(a, b)) return false;

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

md_unit_t md_unit_mul(md_unit_t a, md_unit_t b) {
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

md_unit_t md_unit_div(md_unit_t a, md_unit_t b) {
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

md_unit_t md_unit_add(md_unit_t a, md_unit_t b) {
    if (a.base.raw_bits == b.base.raw_bits) {
        return a;
    }
    return md_unit_none();
}

md_unit_t md_unit_sub(md_unit_t a, md_unit_t b) {
    if (a.base.raw_bits == b.base.raw_bits) {
        return a;
    }
    return md_unit_none();
}

md_unit_t md_unit_inv(md_unit_t unit) {
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

md_unit_t md_unit_pow(md_unit_t unit, int pow) {
    md_unit_t md_unit_pow = {
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
    return md_unit_pow;
}

md_unit_t md_unit_scl(md_unit_t unit, double scl) {
    md_unit_t md_unit_scl = unit;
    md_unit_scl.mult *= scl;
    return md_unit_scl;
}

bool md_unit_convert_inplace_d(double* values, int64_t num_values, md_unit_t cur_unit, md_unit_t new_unit) {
    if (!md_unit_base_equal(cur_unit, new_unit)) {
        MD_LOG_ERROR("Failed to perform unit conversion, the current unit cannot be converted into new unit");
        return false;
    }

    if (md_unit_equal(cur_unit, new_unit)) {
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

bool md_unit_convert_inplace_f(float* values, int64_t num_values, md_unit_t cur_unit, md_unit_t new_unit) {
    if (!md_unit_base_equal(cur_unit, new_unit)) {
        MD_LOG_ERROR("Failed to perform unit conversion, the current unit cannot be converted into new unit");
        return false;
    }

    if (md_unit_equal(cur_unit, new_unit)) {
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


#define PRINT(...) len += snprintf(buf + len, cap - MIN(len, cap), ##__VA_ARGS__)

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

// A bit wierd name, counts the number of base unit dimensions that are not zero
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

size_t internal_print_dims_SI(char* buf, size_t cap, md_unit_t unit) {
    size_t len = 0;
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

size_t internal_print(char* buf, size_t cap, md_unit_t unit) {
    size_t len = 0;

    if (md_unit_empty(unit) || md_unit_unitless(unit)) {
        return 0;
    }

    if (unit.mult == 0) {
        PRINT("0*");
        unit.mult = 1;
    }

    // Try to match against predefined units (exact)
    str_t exact = find_name_from_predefined(unit);
    if (!str_empty(exact)) {
        PRINT(STR_FMT, STR_ARG(exact));
        return len;
    }

    // Try to match inverted against predefined units (exact)
    {
        str_t str = find_name_from_predefined(md_unit_inv(unit));
        if (!str_empty(str)) {
            PRINT("1/"STR_FMT, STR_ARG(str));
            return len;
        }
    }

    // Try to match find an exact matching base and see if we can apply a prefix
    for (size_t i = 0; i < ARRAY_SIZE(predefined_units); ++i) {
        if (md_unit_base_equal(unit, predefined_units[i].unit)) {
            // We found a matching base, we just need to work out a matching prefix scaling
            str_t str  = predefined_units[i].str;
            double mult  = unit.mult / predefined_units[i].unit.mult;
            str_t prefix = STR_LIT("");

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
        if (md_unit_base_equal(md_unit_inv(unit), predefined_units[i].unit)) {
            // We found a matching unit, we just need to work out a matching prefix scaling
            str_t str  = predefined_units[i].str;
            double mult  = unit.mult / predefined_units[i].unit.mult;
            str_t prefix = STR_LIT("");

            if (mult != 1.0) {
                prefix = find_prefix_str_from_value((float)mult);
                if (str_empty(prefix)) {
                    MD_LOG_DEBUG("Unable to find suitable prefix?");
                    return 0;
                }
            }

            PRINT("1/"STR_FMT""STR_FMT, STR_ARG(prefix), STR_ARG(str));
            return len;
        }
    }

    // Try to find best matching base and see where we can go from there
    {
        md_unit_t original_unit = unit;

        size_t num_len = 0;
        size_t den_len = 0;
        int num_arr[32];
        int den_arr[32];
        for (size_t i = 0; i < ARRAY_SIZE(num_arr); ++i) {
            int idx = find_best_matching_predefined(unit);
            if (idx == -1) {
                break;
            }
            num_arr[num_len++] = idx;
            unit = md_unit_div(unit, predefined_units[idx].unit);
            if (unit.base.raw_bits == 0) {
                break;
            }
        }
        if (unit.base.raw_bits != 0) {
            for (size_t i = 0; i < ARRAY_SIZE(den_arr); ++i) {
                int idx = find_best_matching_predefined(md_unit_inv(unit));
                if (idx == -1) {
                    break;
                }
                den_arr[den_len++] = idx;
                unit = md_unit_mul(unit, predefined_units[idx].unit);
                if (unit.base.raw_bits == 0) {
                    break;
                }
            }
        }

        if (unit.base.raw_bits == 0) {
            str_t prefix = STR_LIT("");
            if (unit.mult != 1.0) {
                prefix = find_prefix_str_from_value((float)unit.mult);
                if (str_empty(prefix)) {
                    MD_LOG_DEBUG("Unable to find suitable prefix?");
                    return 0;
                }
            }

            PRINT(STR_FMT, STR_ARG(prefix));
            for (size_t i = 0; i < num_len; ++i) {
                int idx = num_arr[i];
                str_t str = predefined_units[idx].str;
                PRINT(STR_FMT, STR_ARG(str));
            }
            if (den_len > 0) {
                PRINT("/");
                for (size_t i = 0; i < den_len; ++i) {
                    int idx = den_arr[i];
                    str_t str = predefined_units[idx].str;
                    PRINT(STR_FMT, STR_ARG(str));
                }
            }
            return len;
        }
        
        unit = original_unit;
    }

    // Last resort, just print out whatever there is in SI units
    {
        str_t prefix = find_prefix_str_from_value((float)unit.mult);
        if (!str_empty(prefix)) {
            PRINT(STR_FMT, STR_ARG(prefix));
        }
        else if (unit.mult != 1.0) {
            PRINT("%f*", unit.mult);
        }
        unit.mult = 1.0;
        len += internal_print_dims_SI(buf + len, cap - len, mask_pos(unit));
        if (base_count(mask_neg(unit))) {
            PRINT("/");
            len += internal_print_dims_SI(buf + len, cap - len, md_unit_inv(mask_neg(unit)));
        }   
    }

    return len;
}

size_t md_unit_print(char* buf, size_t cap, md_unit_t unit) {
    return internal_print(buf, cap, unit);
}

str_t md_unit_to_string(md_unit_t unit, struct md_allocator_i* alloc) {
    str_t str = {0};

    char buf[128];
    size_t len = md_unit_print(buf, sizeof(buf), unit);
    if (len > 0) {
        str = str_copy((str_t){buf, len}, alloc);
    }

    return str;
}

md_unit_t md_unit_from_string(str_t str) {
    str = str_trim(str);
    if (str_empty(str)) {
        return (md_unit_t) {0};
    }

    size_t loc;
    
    if (str_find_char(&loc, str, '/')) {
        str_t num = str_substr(str, 0, loc);
        str_t den = str_substr(str, loc + 1, str.len);
        return md_unit_div(md_unit_from_string(num), md_unit_from_string(den));
    }

    if (str_find_char(&loc, str, '*')) {
        str_t p1 = str_substr(str, 0, loc);
        str_t p2 = str_substr(str, loc + 1, str.len);
        return md_unit_mul(md_unit_from_string(p1), md_unit_from_string(p2));
    }
    
    if (str_find_char(&loc, str, ' ')) {
        str_t p1 = str_substr(str, 0, loc);
        str_t p2 = str_substr(str, loc + 1, str.len);
        return md_unit_mul(md_unit_from_string(p1), md_unit_from_string(p2));
    }
    
    // Try to match against defined units (exact)
    md_unit_t exact = find_unit_from_predefined(str);
    if (!md_unit_empty(exact)) {
        return exact;
    }
        
    // Try to match against longest matching predefined unit (from right to left)
    str_t base = str;
    str_t exp = {0};
    str_t rest = {0};
    if (str_find_char(&loc, str, '^')) {
        base = str_substr(str, 0, loc);
        size_t i = loc + 1;
        while (i < str.len && (str.ptr[i] == '-' || str.ptr[i] == '+' || is_alpha(str.ptr[i]))) ++i;
        exp = str_substr(str, loc + 1, i);
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
                unit = md_unit_pow(unit, pow);
            }
        }

        if (!str_empty(rest)) {
            unit = md_unit_mul(unit, md_unit_from_string(rest));
        }

        return unit;
    }

    MD_LOG_ERROR("Could not convert string to unit: '%.*s'", (int)str.len, str.ptr);
    return (md_unit_t) {0};
}

md_unit_t md_unit_none(void) {
    return (md_unit_t)UNIT_NONE;
}

md_unit_t md_unit_meter(void) {
    return (md_unit_t)UNIT_METER;
}

md_unit_t md_unit_nanometer(void) {
    return (md_unit_t)UNIT_NANOMETER;
}

md_unit_t md_unit_angstrom(void) {
    return (md_unit_t)UNIT_ANGSTROM;
}

md_unit_t md_unit_kilogram(void) {
    return (md_unit_t)UNIT_KILOGRAM;
}

md_unit_t md_unit_second(void) {
    return (md_unit_t)UNIT_SECOND;
}

md_unit_t md_unit_nanosecond(void) {
    return (md_unit_t)UNIT_NANOSECOND;
}

md_unit_t md_unit_pikosecond(void) {
    return (md_unit_t)UNIT_PIKOSECOND;
}

md_unit_t md_unit_femtosecond(void) {
	return (md_unit_t)UNIT_FEMTOSECOND;
}

md_unit_t md_unit_ampere(void) {
    return (md_unit_t)UNIT_AMPERE;
}

md_unit_t md_unit_mole(void) {
    return (md_unit_t)UNIT_MOLE;
}

md_unit_t md_unit_kelvin(void) {
    return (md_unit_t)UNIT_KELVIN;
}

md_unit_t md_unit_radian(void) {
    return (md_unit_t)UNIT_RADIAN;
}

md_unit_t md_unit_degree(void) {
    return (md_unit_t)UNIT_DEGREE;
}

md_unit_t md_unit_count(void) {
    return (md_unit_t)UNIT_COUNT;
}

md_unit_t md_unit_joule(void) {
    return (md_unit_t)UNIT_JOULE;
}

md_unit_t md_unit_pascal(void) {
    return (md_unit_t)UNIT_PASCAL;
}

md_unit_t md_unit_bar(void) {
    return (md_unit_t)UNIT_BAR;
}
