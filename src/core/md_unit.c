#include "md_unit.h"

#include "md_log.h"
#include "md_common.h"

#include <stdio.h>
#include <math.h>

md_unit_t unit_mul(md_unit_t a, md_unit_t b) {
    md_unit_t result = {
        .dim.length     = a.dim.length  + b.dim.length,
        .dim.mass       = a.dim.mass    + b.dim.mass,
        .dim.time       = a.dim.time    + b.dim.time,
        .dim.current    = a.dim.current + b.dim.current,
        .dim.mole       = a.dim.mole    + b.dim.mole,
        .dim.temp       = a.dim.temp    + b.dim.temp,
        .dim.count      = a.dim.count   + b.dim.count,
        .dim.angle      = a.dim.angle   + b.dim.angle,
    };

    return result;
}

md_unit_t unit_div(md_unit_t a, md_unit_t b) {
    md_unit_t result = {
        .dim.length     = a.dim.length  - b.dim.length,
        .dim.mass       = a.dim.mass    - b.dim.mass,
        .dim.time       = a.dim.time    - b.dim.time,
        .dim.current    = a.dim.current - b.dim.current,
        .dim.mole       = a.dim.mole    - b.dim.mole,
        .dim.temp       = a.dim.temp    - b.dim.temp,
        .dim.count      = a.dim.count   - b.dim.count,
        .dim.angle      = a.dim.angle   - b.dim.angle,
    };

    return result;
}

md_unit_t unit_add(md_unit_t a, md_unit_t b) {
    md_unit_t result = {0};
    if (a.dim.raw_bits == b.dim.raw_bits) {
        result = a;
    }
    return result;
}

md_unit_t unit_sub(md_unit_t a, md_unit_t b) {
    md_unit_t result = {0};
    if (a.dim.raw_bits == b.dim.raw_bits) {
        result = a;
    }
    return result;
}

static bool validate_unit_base(md_unit_base_t unit_base) {
    return
        (unit_base.length == UNIT_LENGTH_ANGSTROM || unit_base.length == UNIT_LENGTH_NANOMETER) &&
        (unit_base.mass == 0) &&
        (unit_base.time == UNIT_TIME_PIKOSECONDS || unit_base.time == UNIT_TIME_NANOSECONDS) &&
        (unit_base.current == 0) &&
        (unit_base.mole == 0) &&
        (unit_base.temp == 0) &&
        (unit_base.count == 0) &&
        (unit_base.angle == UNIT_ANGLE_RADIANS || unit_base.angle == UNIT_ANGLE_DEGREES);
}

static double conversion_table_length[2][2] = {
    {1.0, 10.0},
    {10.0, 1.0},
};

static double conversion_table_time[2][2] = {
    {1.0,  0.001},
    {1000.0, 1.0},
};

static double conversion_table_angle[2][2] = {
    {1.0, 180.0 / PI},
    {PI / 180.0, 1.0},
};

static double get_conversion_factor(md_unit_dim_t dim, md_unit_base_t from, md_unit_base_t to) {
    double factor = 1.0;

    if (dim.length != 0) {
        const double x = conversion_table_length[from.length][to.length];
        factor *= pow(x, (double)dim.length);
    }

    if (dim.time != 0) {
        const double x = conversion_table_time[from.time][to.time];
        factor *= pow(x, (double)dim.time);
    }

    if (dim.angle != 0) {
        const double x = conversion_table_angle[from.time][to.time];
        factor *= pow(x, (double)dim.angle);
    }

    return factor;
}

void unit_convert_d(double* values, int64_t num_values, md_unit_t* unit, md_unit_base_t new_base) {
    if (new_base.raw_bits == 0) return;

    if (!validate_unit_base(new_base)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to perform unit conversion, the new base was invalid.");
        return;
    }

    const double factor = get_conversion_factor(unit->dim, unit->base, new_base);
    if (factor != 1.0) {
        for (int64_t i = 0; i < num_values; ++i) {
            values[i] *= factor;
        }
    }
}

void unit_convert_f(float* values, int64_t num_values, md_unit_t* unit, md_unit_base_t new_base) {
    if (new_base.raw_bits == 0) return;

    if (!validate_unit_base(new_base)) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to perform unit conversion, the new base was invalid.");
        return;
    }

    const double factor = get_conversion_factor(unit->dim, unit->base, new_base);
    if (factor != 1.0) {
        for (int64_t i = 0; i < num_values; ++i) {
            values[i] = (float)(values[i] * factor);
        }
    }
}

static int print_exponent(char* buf, int cap, int exp) {
    switch (exp) {
    case 0: return 0;
    case 1: return 0;
    case 2: return snprintf(buf, cap, u8"²");
    case 3: return snprintf(buf, cap, u8"³");
    case -1: return snprintf(buf, cap, u8"⁻¹");
    case -2: return snprintf(buf, cap, u8"⁻²");
    case -3: return snprintf(buf, cap, u8"⁻³");
    default: return snprintf(buf, cap, "^%i",exp);
    }
}

#define PRINT(str) len += snprintf(buf + len, MAX(0, cap - len), str)

int unit_print_long(char* buf, int cap, md_unit_t unit) {
    const bool have_neg_dim =
        unit.dim.length     < 0 ||
        unit.dim.mass       < 0 ||
        unit.dim.time       < 0 ||
        unit.dim.current    < 0 ||
        unit.dim.temp       < 0 ||
        unit.dim.mole       < 0 ||
        unit.dim.count      < 0 ||
        unit.dim.angle      < 0;

    int len = 0;
    // Write positive indices first
    if (unit.dim.length > 0) {
        switch (unit.base.length) {
        case UNIT_LENGTH_ANGSTROM:  PRINT(u8"Ångström");  break;
        case UNIT_LENGTH_NANOMETER: PRINT(u8"Nanometer"); break;
        default:
            ASSERT(false);
        }
        len += print_exponent(buf + len, MAX(0, cap - len), unit.dim.length);
    }
    if (unit.dim.mass > 0) {
        switch (unit.base.length) {
        case UNIT_LENGTH_ANGSTROM:  PRINT(u8"Kilogram");  break;
        default:
            ASSERT(false);
        }
        len += print_exponent(buf + len, MAX(0, cap - len), unit.dim.mass);
    }
    if (unit.dim.time > 0) {
        switch (unit.base.time) {
        case UNIT_TIME_PIKOSECONDS: PRINT(u8"Picoseconds"); break;
        case UNIT_TIME_NANOSECONDS: PRINT(u8"Nanoseconds"); break;
        default:
            ASSERT(false);
        }
        len += print_exponent(buf + len, MAX(0, cap - len), unit.dim.time);
    }
    if (unit.dim.current > 0) {
        PRINT(u8"Ampere");
        len += print_exponent(buf + len, MAX(0, cap - len), unit.dim.current);
    }
    if (unit.dim.temp > 0) {
        PRINT(u8"Kelvin");
        len += print_exponent(buf + len, MAX(0, cap - len), unit.dim.temp);
    }
    if (unit.dim.count > 0) {
        PRINT(u8"Count");
        len += print_exponent(buf + len, MAX(0, cap - len), unit.dim.count); // Should we print exponent here???
    }
    if (unit.dim.angle > 0) {
        switch (unit.base.angle) {
        case UNIT_ANGLE_RADIANS: PRINT(u8"Radians"); break;
        case UNIT_ANGLE_DEGREES: PRINT(u8"Degrees"); break;
        default:
            ASSERT(false);
        }
        len += print_exponent(buf + len, MAX(0, cap - len), unit.dim.angle); // Should we print exponent here???
    }

    if (have_neg_dim) {
        len += snprintf(buf + len, MAX(0, cap - len), "/");

        // Negative indices
        if (unit.dim.length < 0) {
            switch (unit.base.length) {
            case UNIT_LENGTH_ANGSTROM:  PRINT(u8"Ångström");  break;
            case UNIT_LENGTH_NANOMETER: PRINT(u8"Nanometer"); break;
            default:
                ASSERT(false);
            }
            len += print_exponent(buf + len, MAX(0, cap - len), -unit.dim.length);
        }
        if (unit.dim.mass < 0) {
            switch (unit.base.mass) {
            case UNIT_LENGTH_ANGSTROM:  PRINT(u8"Kilogram");  break;
            default:
                ASSERT(false);
            }
            len += print_exponent(buf + len, MAX(0, cap - len), -unit.dim.mass);
        }
        if (unit.dim.time < 0) {
            switch (unit.base.time) {
            case UNIT_TIME_PIKOSECONDS: PRINT(u8"Picoseconds"); break;
            case UNIT_TIME_NANOSECONDS: PRINT(u8"Nanoseconds"); break;
            default:
                ASSERT(false);
            }
            len += print_exponent(buf + len, MAX(0, cap - len), -unit.dim.time);
        }
        if (unit.dim.current < 0) {
            PRINT(u8"Ampere");
            len += print_exponent(buf + len, MAX(0, cap - len), -unit.dim.current);
        }
        if (unit.dim.temp < 0) {
            PRINT(u8"Kelvin");
            len += print_exponent(buf + len, MAX(0, cap - len), -unit.dim.temp);
        }
        if (unit.dim.count < 0) {
            PRINT(u8"Count");
            len += print_exponent(buf + len, MAX(0, cap - len), -unit.dim.count); // Should we print exponent here???
        }
        if (unit.dim.angle < 0) {
            switch (unit.base.angle) {
            case UNIT_ANGLE_RADIANS: PRINT(u8"Radians"); break;
            case UNIT_ANGLE_DEGREES: PRINT(u8"Degrees"); break;
            default:
                ASSERT(false);
            }
            len += print_exponent(buf + len, MAX(0, cap - len), -unit.dim.angle); // Should we print exponent here???
        }
    }
    return len;
}

int unit_print(char* buf, int cap, md_unit_t unit) {
    const bool have_neg_dim =
        unit.dim.length     < 0 ||
        unit.dim.mass       < 0 ||
        unit.dim.time       < 0 ||
        unit.dim.current    < 0 ||
        unit.dim.temp       < 0 ||
        unit.dim.mole       < 0 ||
        unit.dim.count      < 0 ||
        unit.dim.angle      < 0;

    int len = 0;
    // Write positive indices first
    if (unit.dim.length > 0) {
        switch (unit.base.length) {
        case UNIT_LENGTH_ANGSTROM:  PRINT(u8"Å");  break;
        case UNIT_LENGTH_NANOMETER: PRINT(u8"nm"); break;
        default:
            ASSERT(false);
        }
        len += print_exponent(buf + len, MAX(0, cap - len), unit.dim.length);
    }
    if (unit.dim.mass > 0) {
        PRINT(u8"Kg");
        len += print_exponent(buf + len, MAX(0, cap - len), unit.dim.mass);
    }
    if (unit.dim.time > 0) {
        switch (unit.base.length) {
        case UNIT_TIME_PIKOSECONDS: PRINT(u8"ps"); break;
        case UNIT_TIME_NANOSECONDS: PRINT(u8"ns"); break;
        default:
            ASSERT(false);
        }
        len += print_exponent(buf + len, MAX(0, cap - len), unit.dim.time);
    }
    if (unit.dim.current > 0) {
        PRINT(u8"A");
        len += print_exponent(buf + len, MAX(0, cap - len), unit.dim.current);
    }
    if (unit.dim.temp > 0) {
        PRINT(u8"K");
        len += print_exponent(buf + len, MAX(0, cap - len), unit.dim.temp);
    }
    if (unit.dim.angle > 0) {
        PRINT(u8"rad");
        len += print_exponent(buf + len, MAX(0, cap - len), unit.dim.angle); // Should we print exponent here???
    }

    if (have_neg_dim) {
        PRINT("/");

        // Negative indices
        if (unit.dim.length < 0) {
            PRINT(u8"Å");
            len += print_exponent(buf + len, MAX(0, cap - len), -unit.dim.length);
        }
        if (unit.dim.mass < 0) {
            PRINT(u8"Kg");
            len += print_exponent(buf + len, MAX(0, cap - len), -unit.dim.mass);
        }
        if (unit.dim.time < 0) {
            switch (unit.base.time) {
            case UNIT_TIME_PIKOSECONDS: PRINT(u8"ps"); break;
            case UNIT_TIME_NANOSECONDS: PRINT(u8"ns"); break;
            default:
                ASSERT(false);
            }
            len += print_exponent(buf + len, MAX(0, cap - len), -unit.dim.time);
        }
        if (unit.dim.current < 0) {
            PRINT(u8"A");
            len += print_exponent(buf + len, MAX(0, cap - len), -unit.dim.current);
        }
        if (unit.dim.temp < 0) {
            PRINT(u8"K");
            len += print_exponent(buf + len, MAX(0, cap - len), -unit.dim.temp);
        }
        if (unit.dim.angle < 0) {
            switch (unit.base.angle) {
            case UNIT_ANGLE_RADIANS: PRINT(u8"rad"); break;
            case UNIT_ANGLE_DEGREES: PRINT(u8"°"); break;
            default:
                ASSERT(false);
            }
            len += print_exponent(buf + len, MAX(0, cap - len), -unit.dim.angle); // Should we print exponent here???
        }
    }
    return len;
}
