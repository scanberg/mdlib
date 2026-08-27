#include <core/md_unit.h>

#include <core/md_log.h>
#include <core/md_common.h>
#include <core/md_str.h>
#include <core/md_parse.h>

#include <float.h>
#include <math.h>
#include <stdio.h> // snprintf

#define S(str) {(const char*)str"", sizeof(str)-1}

typedef struct {
    double value;
    str_t  str;
    bool   alias;   // Alternative spelling, only considered when parsing
} si_prefix_t;

static const si_prefix_t si_prefixes[] = {
    {1e-24, S(u8"y"), false},
    {1e-21, S(u8"z"), false},
    {1e-18, S(u8"a"), false},
    {1e-15, S(u8"f"), false},
    {1e-12, S(u8"p"), false},
    {1e-9,  S(u8"n"), false},
    {1e-6,  S(u8"μ"), false},
    {1e-6,  S(u8"u"), true },
    {1e-3,  S(u8"m"), false},
    {1e-2,  S(u8"c"), false},
    {1e3,   S(u8"k"), false},
    {1e6,   S(u8"M"), false},
    {1e9,   S(u8"G"), false},
    {1e12,  S(u8"T"), false},
    {1e15,  S(u8"P"), false},
    {1e18,  S(u8"E"), false},
    {1e21,  S(u8"Z"), false},
    {1e24,  S(u8"Y"), false},
};

// Relative comparison, the scale factors are the result of floating point arithmetic
// and will not always land exactly on the tabulated values.
// The comparison has to remain relative over the full range, the scale factors span
// many orders of magnitude and 1e-30 is not the same thing as 1e-24.
static inline bool value_equal(double a, double b) {
    if (a == b) return true;
    const double d = fabs(a - b);
    const double scale = MAX(fabs(a), fabs(b));
    return d <= 16.0 * DBL_EPSILON * scale;
}

static str_t find_prefix_str_from_value(double value) {
    for (size_t i = 0; i < ARRAY_SIZE(si_prefixes); ++i) {
        if (si_prefixes[i].alias) continue;
        if (value_equal(value, si_prefixes[i].value)) return si_prefixes[i].str;
    }
    return (str_t){0};
}

static double find_prefix_value_from_str(str_t str) {
    for (size_t i = 0; i < ARRAY_SIZE(si_prefixes); ++i) {
        if (str_eq(str, si_prefixes[i].str)) return si_prefixes[i].value;
    }
    return 0;
}

typedef struct {
    md_unit_t   unit;
    str_t       str;
} unit_name_t;

// Physical constants. The elementary charge is exact by the 2019 SI redefinition; the Bohr radius
// is CODATA 2022; the Debye is exactly 1e-21 / c coulomb metre.
#define ELEMENTARY_CHARGE_COULOMB 1.602176634e-19
#define BOHR_RADIUS_METER         5.29177210544e-11
#define DEBYE_COULOMB_METER       3.33564095198152e-30

#define UNIT_NONE           {0}

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
#define UNIT_BOHR           {.base = {.dim = {.length = 1,}},   .mult = BOHR_RADIUS_METER}
#define UNIT_MILLISECOND    {.base = {.dim = {.time = 1,}},     .mult = 1e-3}
#define UNIT_NANOSECOND     {.base = {.dim = {.time = 1,}},     .mult = 1e-9}
#define UNIT_PICOSECOND     {.base = {.dim = {.time = 1,}},     .mult = 1e-12}
#define UNIT_FEMTOSECOND    {.base = {.dim = {.time = 1,}},     .mult = 1e-15}

// The scale is expressed in SI base units, one degree is PI/180 radians
#define UNIT_DEGREE         {.base = {.dim = {.angle = 1,}},    .mult = DEG_TO_RAD(1.0)}

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
#define UNIT_ELEMENTARY_CHARGE {.base = {.dim = {.time = 1, .current = 1}}, .mult = ELEMENTARY_CHARGE_COULOMB}

// Dipole moment: charge times length. The Debye and the atomic unit (e a0) differ only in scale,
// so md_unit_conversion_factor moves between them and between either and C m.
#define UNIT_DEBYE          {.base = {.dim = {.length = 1, .time = 1, .current = 1}}, .mult = DEBYE_COULOMB_METER}
#define UNIT_ELEMENTARY_CHARGE_BOHR {.base = {.dim = {.length = 1, .time = 1, .current = 1}}, .mult = ELEMENTARY_CHARGE_COULOMB * BOHR_RADIUS_METER}
#define UNIT_ELECTRONVOLT   {.base = {.dim = {.mass = 1, .length =  2, .time = -2,}}, .mult = ELEMENTARY_CHARGE_COULOMB}
#define UNIT_HERTZ          {.base = {.dim = {.time = -1,}}, .mult = 1.0}

// The order matters, the first entry which matches is the one used when printing
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
    {UNIT_PICOSECOND,   S(u8"ps")},
    {UNIT_NANOMETER,    S(u8"nm")},
    {UNIT_ANGSTROM,     S(u8"Å")},
    {UNIT_BOHR,         S(u8"bohr")},
    {UNIT_DEGREE,       S(u8"deg")},
    {UNIT_DEGREE,       S(u8"°")},
    {UNIT_PASCAL,       S(u8"Pa")},
    {UNIT_BAR,          S(u8"bar")},
    {UNIT_JOULE,        S(u8"J")},
    {UNIT_NEWTON,       S(u8"N")},
    {UNIT_NEWTONMETER,  S(u8"Nm")},
    {UNIT_VOLT,         S(u8"V")},
    {UNIT_WATT,         S(u8"W")},
    {UNIT_COULOMB,      S(u8"C")},
    {UNIT_AMPEREHOUR,   S(u8"Ah")},
    {UNIT_ELEMENTARY_CHARGE, S(u8"e")},
    {UNIT_DEBYE,        S(u8"D")},
    {UNIT_ELECTRONVOLT, S(u8"eV")},
    {UNIT_HERTZ,        S(u8"Hz")},
};

static bool find_unit_from_predefined(md_unit_t* unit, str_t name) {
    for (size_t i = 0; i < ARRAY_SIZE(predefined_units); ++i) {
        if (str_eq(name, predefined_units[i].str)) {
            *unit = predefined_units[i].unit;
            return true;
        }
    }
    return false;
}

static str_t find_name_from_predefined(md_unit_t unit) {
    for (size_t i = 0; i < ARRAY_SIZE(predefined_units); ++i) {
        if (md_unit_equal(predefined_units[i].unit, unit)) {
            return predefined_units[i].str;
        }
    }
    return (str_t){0};
}

// Exponents considered when matching a unit against a power of a predefined unit
#define MAX_UNIT_EXPONENT 12

// True if the name can carry an exponent, a name which already contains one cannot ('m²^6')
static bool name_can_be_raised(str_t str) {
    for (size_t i = 0; i < str.len; ++i) {
        const unsigned char c = (unsigned char)str.ptr[i];
        // '^' or the trailing byte of the superscript two and three in UTF-8
        if (c == '^' || c == 0xB2 || c == 0xB3) return false;
    }
    return true;
}

// True if none of the dimensions are inverted, we prefer to express leftover dimensions
// using such units, i.e. 'm/s^2' rather than 'mHz^2'
static bool base_is_positive(md_unit_t unit) {
    for (size_t i = 0; i < ARRAY_SIZE(unit.base.arr); ++i) {
        if (unit.base.arr[i] < 0) return false;
    }
    return true;
}

// Every dimension of ref, raised to exponent, has to line up with unit, otherwise there is no match.
// The score is proportional to the number of dimensions covered, plus a bonus for candidates which
// leave a scale we can express. This is what decides between otherwise equal candidates,
// e.g. 'count/Å^3' rather than 'Zcount/mm^3'.
static int compute_match_score(md_unit_t unit, md_unit_t ref, int exponent) {
    int score = 0;
    for (size_t i = 0; i < ARRAY_SIZE(ref.base.arr); ++i) {
        if (ref.base.arr[i] == 0) continue;
        if (unit.base.arr[i] != ref.base.arr[i] * exponent) return 0;
        score += 10;
    }
    if (score == 0) {
        return 0;
    }

    const double leftover = md_unit_scale(unit) / md_unit_scale(md_unit_pow(ref, exponent));
    if (value_equal(leftover, 1.0)) {
        score += 8;
    } else if (!str_empty(find_prefix_str_from_value(leftover))) {
        score += 5;
    }

    if (base_is_positive(ref)) {
        score += 2;
    }

    return score;
}

// Finds the predefined unit (and exponent) which covers the most dimensions of unit
static int find_best_matching_predefined(int* exponent, md_unit_t unit) {
    int best_score = 0;
    int best_idx   = -1;
    int best_exp   = 0;

    // Prefer the lowest magnitude exponents
    for (int mag = 1; mag <= MAX_UNIT_EXPONENT; ++mag) {
        for (int sign = 1; sign >= -1; sign -= 2) {
            const int exp = mag * sign;
            for (size_t i = 0; i < ARRAY_SIZE(predefined_units); ++i) {
                if (mag > 1 && !name_can_be_raised(predefined_units[i].str)) continue;
                int score = compute_match_score(unit, predefined_units[i].unit, exp);
                if (score > best_score) {
                    best_score = score;
                    best_idx = (int)i;
                    best_exp = exp;
                }
            }
        }
    }

    if (exponent) *exponent = best_exp;
    return best_idx;
}

bool md_unit_equal(md_unit_t a, md_unit_t b) {
    return md_unit_base_equal(a, b) && value_equal(md_unit_scale(a), md_unit_scale(b));
}

md_unit_t md_unit_mul(md_unit_t a, md_unit_t b) {
    md_unit_t result = {
        .base = {
            .dim = {
                .length     = (int8_t)(a.base.dim.length  + b.base.dim.length),
                .mass       = (int8_t)(a.base.dim.mass    + b.base.dim.mass),
                .time       = (int8_t)(a.base.dim.time    + b.base.dim.time),
                .current    = (int8_t)(a.base.dim.current + b.base.dim.current),
                .mole       = (int8_t)(a.base.dim.mole    + b.base.dim.mole),
                .temp       = (int8_t)(a.base.dim.temp    + b.base.dim.temp),
                .angle      = (int8_t)(a.base.dim.angle   + b.base.dim.angle),
                .count      = (int8_t)(a.base.dim.count   + b.base.dim.count),
            },
        },
        .mult = md_unit_scale(a) * md_unit_scale(b),
    };
    return result;
}

md_unit_t md_unit_div(md_unit_t a, md_unit_t b) {
    md_unit_t result = {
        .base = {
            .dim = {
                .length     = (int8_t)(a.base.dim.length  - b.base.dim.length),
                .mass       = (int8_t)(a.base.dim.mass    - b.base.dim.mass),
                .time       = (int8_t)(a.base.dim.time    - b.base.dim.time),
                .current    = (int8_t)(a.base.dim.current - b.base.dim.current),
                .mole       = (int8_t)(a.base.dim.mole    - b.base.dim.mole),
                .temp       = (int8_t)(a.base.dim.temp    - b.base.dim.temp),
                .angle      = (int8_t)(a.base.dim.angle   - b.base.dim.angle),
                .count      = (int8_t)(a.base.dim.count   - b.base.dim.count),
            },
        },
        .mult = md_unit_scale(a) / md_unit_scale(b),
    };
    return result;
}

md_unit_t md_unit_add(md_unit_t a, md_unit_t b) {
    if (md_unit_equal(a, b)) {
        return a;
    }
    return md_unit_none();
}

md_unit_t md_unit_sub(md_unit_t a, md_unit_t b) {
    if (md_unit_equal(a, b)) {
        return a;
    }
    return md_unit_none();
}

md_unit_t md_unit_inv(md_unit_t unit) {
    md_unit_t inv = {
        .base = {
            .dim = {
                .length     = (int8_t)(-unit.base.dim.length),
                .mass       = (int8_t)(-unit.base.dim.mass),
                .time       = (int8_t)(-unit.base.dim.time),
                .current    = (int8_t)(-unit.base.dim.current),
                .mole       = (int8_t)(-unit.base.dim.mole),
                .temp       = (int8_t)(-unit.base.dim.temp),
                .angle      = (int8_t)(-unit.base.dim.angle),
                .count      = (int8_t)(-unit.base.dim.count),
            },
        },
        .mult = 1.0 / md_unit_scale(unit),
    };
    return inv;
}

md_unit_t md_unit_pow(md_unit_t unit, int exponent) {
    md_unit_t result = {
        .base = {
            .dim = {
                .length     = (int8_t)(unit.base.dim.length  * exponent),
                .mass       = (int8_t)(unit.base.dim.mass    * exponent),
                .time       = (int8_t)(unit.base.dim.time    * exponent),
                .current    = (int8_t)(unit.base.dim.current * exponent),
                .mole       = (int8_t)(unit.base.dim.mole    * exponent),
                .temp       = (int8_t)(unit.base.dim.temp    * exponent),
                .angle      = (int8_t)(unit.base.dim.angle   * exponent),
                .count      = (int8_t)(unit.base.dim.count   * exponent),
            },
        },
        .mult = exponent == 0 ? 1.0 : pow(md_unit_scale(unit), exponent),
    };
    return result;
}

md_unit_t md_unit_scl(md_unit_t unit, double scl) {
    md_unit_t result = unit;
    result.mult = md_unit_scale(unit) * scl;
    return result;
}

bool md_unit_conversion_factor(double* factor, md_unit_t from_unit, md_unit_t to_unit) {
    ASSERT(factor);

    if (!md_unit_base_equal(from_unit, to_unit)) {
        return false;
    }

    // value_in_SI = value * scale, thus value_in_to = value_in_from * (from_scale / to_scale)
    *factor = md_unit_scale(from_unit) / md_unit_scale(to_unit);
    return true;
}

#define PRINT(...) len += (size_t)snprintf(buf + MIN(len, cap), cap - MIN(len, cap), ##__VA_ARGS__)

// Each print_* procedure writes from the beginning of buf and returns the length written.
// A return value of zero means that the representation was not applicable and the next one should be tried.
// The printed form is always parseable, factors are separated by '*' and a denominator holding more
// than a single factor is parenthesized.

// A scale which cannot be folded into a prefix is printed verbatim, using the shortest
// representation which still reads back as the same value
static size_t print_scale(char* buf, size_t cap, double scale) {
    size_t len = 0;
    char tmp[64];
    for (int precision = 6; precision <= 17; ++precision) {
        snprintf(tmp, sizeof(tmp), "%.*g", precision, scale);
        if (value_equal(parse_float(str_from_cstr(tmp)), scale)) {
            break;
        }
    }
    PRINT("%s*", tmp);
    return len;
}

// Exact match against a predefined unit, e.g. 'Pa'
static size_t print_predefined(char* buf, size_t cap, md_unit_t unit) {
    size_t len = 0;

    str_t str = find_name_from_predefined(unit);
    if (!str_empty(str)) {
        PRINT(STR_FMT, STR_ARG(str));
        return len;
    }

    str = find_name_from_predefined(md_unit_inv(unit));
    if (!str_empty(str)) {
        PRINT("1/"STR_FMT, STR_ARG(str));
        return len;
    }

    return 0;
}

// Match a predefined unit raised to a power and express the remaining scale as an SI prefix,
// e.g. 'nm', '1/ns' or 'Å^3'. The prefix is applied before the exponent, nm^2 is (1e-9 m)^2.
static size_t print_predefined_pow(char* buf, size_t cap, md_unit_t unit) {
    size_t len = 0;

    // The first pass only accepts units which match without a prefix, we would rather
    // print 'nm^2' than 'am²' even though the two are equivalent
    for (int pass = 0; pass < 2; ++pass) {
        for (int mag = 1; mag <= MAX_UNIT_EXPONENT; ++mag) {
            for (int sign = 1; sign >= -1; sign -= 2) {
                const int exp = mag * sign;
                for (size_t i = 0; i < ARRAY_SIZE(predefined_units); ++i) {
                    if (mag > 1 && !name_can_be_raised(predefined_units[i].str)) continue;

                    const md_unit_t ref = md_unit_pow(predefined_units[i].unit, exp);
                    if (!md_unit_base_equal(unit, ref)) {
                        continue;
                    }

                    str_t prefix = STR_LIT("");
                    const double leftover = md_unit_scale(unit) / md_unit_scale(ref);
                    if (!value_equal(leftover, 1.0)) {
                        if (pass == 0) {
                            continue;
                        }
                        // The prefix binds before the exponent, i.e. nm^2 is (1e-9 m)^2
                        prefix = find_prefix_str_from_value(pow(leftover, 1.0 / exp));
                        if (str_empty(prefix)) {
                            continue;
                        }
                    }

                    if (exp < 0) {
                        PRINT("1/");
                    }
                    PRINT(STR_FMT""STR_FMT, STR_ARG(prefix), STR_ARG(predefined_units[i].str));
                    if (mag > 1) {
                        PRINT("^%i", mag);
                    }
                    return len;
                }
            }
        }
    }

    return 0;
}

// Decompose into a product/quotient of predefined units, e.g. 'kJ/mol', 'm/s^2' or 'count/Å^3'
static size_t print_decomposed(char* buf, size_t cap, md_unit_t unit) {
    size_t len = 0;

    struct {
        int idx;
        int exp;
    } term[8];
    size_t num_terms = 0;

    md_unit_t rem = unit;
    while (rem.base.raw_bits != 0 && num_terms < ARRAY_SIZE(term)) {
        int exp = 0;
        int idx = find_best_matching_predefined(&exp, rem);
        if (idx == -1) break;
        term[num_terms].idx = idx;
        term[num_terms].exp = exp;
        num_terms += 1;
        rem = md_unit_div(rem, md_unit_pow(predefined_units[idx].unit, exp));
    }

    if (rem.base.raw_bits != 0) {
        return 0;
    }

    // The leftover scale is folded into the prefix of a single term, taking its exponent into
    // account since the prefix binds tighter, i.e. 'ns^2' is (1e-9 s)^2 and not 1e-9 * s^2.
    // If no term can absorb it, the scale is printed as an explicit factor.
    const double leftover = md_unit_scale(rem);
    str_t prefix = STR_LIT("");
    size_t prefix_term = num_terms;
    if (!value_equal(leftover, 1.0)) {
        for (size_t i = 0; i < num_terms; ++i) {
            str_t str = find_prefix_str_from_value(pow(leftover, 1.0 / term[i].exp));
            if (!str_empty(str)) {
                prefix = str;
                prefix_term = i;
                break;
            }
        }
        if (prefix_term == num_terms) {
            len += print_scale(buf, cap, leftover);
        }
    }

    size_t num_count = 0;
    size_t den_count = 0;
    for (size_t i = 0; i < num_terms; ++i) {
        if (term[i].exp > 0) num_count += 1;
        else den_count += 1;
    }

    if (num_count == 0) {
        PRINT("1");
    }
    for (size_t i = 0, printed = 0; i < num_terms; ++i) {
        if (term[i].exp < 0) continue;
        if (printed++ > 0) PRINT("*");
        if (i == prefix_term) PRINT(STR_FMT, STR_ARG(prefix));
        PRINT(STR_FMT, STR_ARG(predefined_units[term[i].idx].str));
        if (term[i].exp > 1) PRINT("^%i", term[i].exp);
    }

    if (den_count > 0) {
        PRINT("/");
        if (den_count > 1) PRINT("(");
        for (size_t i = 0, printed = 0; i < num_terms; ++i) {
            if (term[i].exp > 0) continue;
            if (printed++ > 0) PRINT("*");
            if (i == prefix_term) PRINT(STR_FMT, STR_ARG(prefix));
            PRINT(STR_FMT, STR_ARG(predefined_units[term[i].idx].str));
            if (term[i].exp < -1) PRINT("^%i", -term[i].exp);
        }
        if (den_count > 1) PRINT(")");
    }

    return len;
}

// Last resort for units which cannot be decomposed, spell out the SI base dimensions,
// e.g. '1e-30*m^12*s^-3'
static size_t print_SI(char* buf, size_t cap, md_unit_t unit) {
    static const str_t dim_str[8] = {
        S(u8"m"), S(u8"kg"), S(u8"s"), S(u8"A"), S(u8"mol"), S(u8"K"), S(u8"rad"), S(u8"count"),
    };

    size_t len = 0;

    const double scale = md_unit_scale(unit);
    if (!value_equal(scale, 1.0)) {
        len += print_scale(buf, cap, scale);
    }

    for (size_t i = 0, printed = 0; i < ARRAY_SIZE(unit.base.arr); ++i) {
        const int exp = unit.base.arr[i];
        if (exp == 0) continue;
        if (printed++ > 0) PRINT("*");
        PRINT(STR_FMT, STR_ARG(dim_str[i]));
        if (exp != 1) {
            PRINT("^%i", exp);
        }
    }

    return len;
}

size_t md_unit_print(char* buf, size_t cap, md_unit_t unit) {
    if (!buf || cap == 0) {
        return 0;
    }

    buf[0] = '\0';

    // A dimensionless unit has no representation, any scale it carries belongs to the value itself
    if (md_unit_is_none(unit)) {
        return 0;
    }

    size_t len = 0;
    if ((len = print_predefined(buf, cap, unit)) > 0) return len;
    if ((len = print_predefined_pow(buf, cap, unit)) > 0) return len;
    if ((len = print_decomposed(buf, cap, unit)) > 0) return len;
    return print_SI(buf, cap, unit);
}

#undef PRINT

static bool parse_unit(md_unit_t* unit, str_t str);

static bool parse_term(md_unit_t* unit, str_t str) {
    // A numeric factor, either as the numerator of an inverse unit ('1/s') or as an explicit scale ('1e-10*m')
    if (is_digit(str.ptr[0]) || str.ptr[0] == '.') {
        if (!is_float(str)) {
            return false;
        }
        *unit = md_unit_scl(md_unit_none(), parse_float(str));
        return true;
    }

    // Try to match against the predefined units (exact)
    if (find_unit_from_predefined(unit, str)) {
        return true;
    }

    // Split off a potential exponent, e.g. 's^-2' or 'm^3kg'
    str_t base = str;
    str_t exp  = {0};
    str_t rest = {0};
    size_t loc;
    if (str_find_char(&loc, str, '^')) {
        base = str_substr(str, 0, loc);
        size_t i = loc + 1;
        if (i < str.len && (str.ptr[i] == '-' || str.ptr[i] == '+')) ++i;
        while (i < str.len && is_digit(str.ptr[i])) ++i;
        exp  = str_substr(str, loc + 1, i - (loc + 1));
        rest = str_substr(str, i, str.len);
    }

    // Match against the longest predefined unit which the base ends with, the remainder is an SI prefix
    int best_match_idx = -1;
    for (size_t i = 0; i < ARRAY_SIZE(predefined_units); ++i) {
        if (str_ends_with(base, predefined_units[i].str)) {
            if (best_match_idx == -1 || predefined_units[i].str.len > predefined_units[best_match_idx].str.len) {
                best_match_idx = (int)i;
            }
        }
    }

    if (best_match_idx == -1) {
        return false;
    }

    md_unit_t result = predefined_units[best_match_idx].unit;

    str_t prefix = str_substr(base, 0, base.len - predefined_units[best_match_idx].str.len);
    if (!str_empty(prefix)) {
        const double prefix_scale = find_prefix_value_from_str(prefix);
        if (prefix_scale == 0) {
            return false;
        }
        result = md_unit_scl(result, prefix_scale);
    }

    if (!str_empty(exp)) {
        int32_t exponent;
        if (!str_extract_i32(&exponent, &exp)) {
            return false;
        }
        result = md_unit_pow(result, exponent);
    }

    if (!str_empty(rest)) {
        md_unit_t rest_unit;
        if (!parse_unit(&rest_unit, rest)) {
            return false;
        }
        result = md_unit_mul(result, rest_unit);
    }

    *unit = result;
    return true;
}

// Finds the last operator which is not enclosed in parenthesis.
// Searching from the right yields left associativity, 'a/b/c' is (a/b)/c.
static bool find_last_operator(size_t* loc, str_t str) {
    int depth = 0;
    for (size_t i = str.len; i > 0; --i) {
        const char c = str.ptr[i-1];
        if (c == ')') {
            depth += 1;
        } else if (c == '(') {
            depth -= 1;
        } else if (depth == 0 && (c == '*' || c == '/' || c == ' ')) {
            *loc = i - 1;
            return true;
        }
    }
    return false;
}

static bool parse_unit(md_unit_t* unit, str_t str) {
    str = str_trim(str);
    if (str_empty(str)) {
        *unit = md_unit_none();
        return true;
    }

    size_t loc;
    if (find_last_operator(&loc, str)) {
        md_unit_t lhs, rhs;
        if (!parse_unit(&lhs, str_substr(str, 0, loc)) || !parse_unit(&rhs, str_substr(str, loc + 1, str.len))) {
            return false;
        }
        *unit = (str.ptr[loc] == '/') ? md_unit_div(lhs, rhs) : md_unit_mul(lhs, rhs);
        return true;
    }

    if (str.ptr[0] == '(') {
        if (str.ptr[str.len - 1] != ')') {
            return false;
        }
        return parse_unit(unit, str_substr(str, 1, str.len - 2));
    }

    return parse_term(unit, str);
}

bool md_unit_parse(md_unit_t* unit, str_t str) {
    ASSERT(unit);

    if (!parse_unit(unit, str)) {
        *unit = md_unit_none();
        return false;
    }

    return true;
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

md_unit_t md_unit_bohr_radius(void) {
    return (md_unit_t)UNIT_BOHR;
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

md_unit_t md_unit_picosecond(void) {
    return (md_unit_t)UNIT_PICOSECOND;
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

md_unit_t md_unit_electronvolt(void) {
    return (md_unit_t)UNIT_ELECTRONVOLT;
}

md_unit_t md_unit_hertz(void) {
    return (md_unit_t)UNIT_HERTZ;
}

md_unit_t md_unit_pascal(void) {
    return (md_unit_t)UNIT_PASCAL;
}

md_unit_t md_unit_bar(void) {
    return (md_unit_t)UNIT_BAR;
}

md_unit_t md_unit_coulomb(void) {
    return (md_unit_t)UNIT_COULOMB;
}

md_unit_t md_unit_elementary_charge(void) {
    return (md_unit_t)UNIT_ELEMENTARY_CHARGE;
}

md_unit_t md_unit_debye(void) {
    return (md_unit_t)UNIT_DEBYE;
}

md_unit_t md_unit_elementary_charge_bohr(void) {
    return (md_unit_t)UNIT_ELEMENTARY_CHARGE_BOHR;
}
