#include <md_pot.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_parse.h>
#include <core/md_str.h>

#define BOHR_TO_ANGSTROM 0.5291772109029999

typedef enum {
    SECTION_NONE = 0,
    SECTION_ENVIRONMENT,
    SECTION_CHARGES,
    SECTION_POLARIZABILITIES,
    SECTION_SKIP,
} section_t;

// Unique strings only. The count of distinct fragment/atom names in a .pot is small
// (a handful of fragment types times their atom count), so a linear scan is fine and
// keeps the returned str_t:s pointer-comparable.
static str_t pot_intern(md_pot_t* pot, str_t str, md_allocator_i* alloc) {
    str_t empty = {0};
    if (str_empty(str)) return empty;
    for (size_t i = 0; i < md_array_size(pot->strings); ++i) {
        if (str_eq(pot->strings[i], str)) {
            return pot->strings[i];
        }
    }
    str_t cpy = str_copy(str, alloc);
    md_array_push(pot->strings, cpy, alloc);
    return cpy;
}

static void pot_copy_symbol(char dst[4], str_t tok, int line) {
    size_t len = MIN(tok.len, (size_t)3);
    if (tok.len > 3) {
        MD_LOG_INFO("POT: line %i: element symbol '%.*s' is longer than 3 characters and was truncated", line, (int)tok.len, tok.ptr);
    }
    MEMCPY(dst, tok.ptr, len);
    dst[len] = '\0';
}

static md_pot_unit_t pot_parse_unit(str_t str) {
    if (str_eq_cstr_ignore_case(str, "angstrom") ||
        str_eq_cstr_ignore_case(str, "angstroms") ||
        str_eq_cstr_ignore_case(str, "aa") ||
        str_eq_cstr_ignore_case(str, "a")) {
        return MD_POT_UNIT_ANGSTROM;
    }
    if (str_eq_cstr_ignore_case(str, "au") ||
        str_eq_cstr_ignore_case(str, "bohr") ||
        str_eq_cstr_ignore_case(str, "bohrs") ||
        str_eq_cstr_ignore_case(str, "atomic")) {
        return MD_POT_UNIT_BOHR;
    }
    return MD_POT_UNIT_UNKNOWN;
}

double md_pot_unit_to_angstrom(md_pot_unit_t unit) {
    // The format's implicit default is angstrom, so UNKNOWN is treated as such.
    return (unit == MD_POT_UNIT_BOHR) ? BOHR_TO_ANGSTROM : 1.0;
}

bool md_pot_parse_str(md_pot_t* pot, str_t str, md_allocator_i* alloc) {
    if (!pot) {
        MD_LOG_ERROR("POT: pot object was NULL");
        return false;
    }
    if (!alloc) {
        MD_LOG_ERROR("POT: allocator was NULL");
        return false;
    }

    MEMSET(pot, 0, sizeof(md_pot_t));

    section_t section = SECTION_NONE;
    bool in_xyz = false;    // Set once the 'xyz:' directive has been seen within @environment
    int  line_num = 0;
    str_t line;

    while (str_extract_line(&line, &str)) {
        line_num += 1;

        str_t l = str_trim(line);
        if (str_empty(l))  continue;
        if (l.ptr[0] == '#' || l.ptr[0] == '!') continue;

        if (l.ptr[0] == '@') {
            str_t kw  = {0};
            str_t rem = str_substr(l, 1, SIZE_MAX);
            extract_token(&kw, &rem);

            if (str_eq_cstr_ignore_case(kw, "end")) {
                if (section == SECTION_NONE) {
                    MD_LOG_INFO("POT: line %i: '@end' without an open section", line_num);
                }
                section = SECTION_NONE;
                in_xyz  = false;
            } else if (str_eq_cstr_ignore_case(kw, "environment")) {
                section = SECTION_ENVIRONMENT;
                in_xyz  = false;
            } else if (str_eq_cstr_ignore_case(kw, "charges")) {
                section = SECTION_CHARGES;
            } else if (str_eq_cstr_ignore_case(kw, "polarizabilities")) {
                section = SECTION_POLARIZABILITIES;
            } else {
                MD_LOG_INFO("POT: line %i: skipping unsupported section '@%.*s'", line_num, (int)kw.len, kw.ptr);
                section = SECTION_SKIP;
            }
            continue;
        }

        if (section == SECTION_SKIP) {
            continue;
        }

        if (section == SECTION_NONE) {
            MD_LOG_ERROR("POT: line %i: data outside of any section: '%.*s'", line_num, (int)l.len, l.ptr);
            goto fail;
        }

        str_t  tok[10];
        str_t  rem = l;
        size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &rem);

        if (section == SECTION_ENVIRONMENT && !in_xyz) {
            // Directives of the form 'key: value'. Coordinate rows never contain a colon.
            size_t loc;
            if (!str_find_char(&loc, l, ':')) {
                MD_LOG_ERROR("POT: line %i: expected a directive ('units:' or 'xyz:') within @environment, got '%.*s'", line_num, (int)l.len, l.ptr);
                goto fail;
            }
            str_t key = str_trim(str_substr(l, 0, loc));
            str_t val = str_trim(str_substr(l, loc + 1, SIZE_MAX));

            if (str_eq_cstr_ignore_case(key, "xyz")) {
                in_xyz = true;
                if (!str_empty(val)) {
                    MD_LOG_ERROR("POT: line %i: unexpected trailing content after 'xyz:': '%.*s'", line_num, (int)val.len, val.ptr);
                    goto fail;
                }
            } else if (str_eq_cstr_ignore_case(key, "units")) {
                md_pot_unit_t unit = pot_parse_unit(val);
                if (unit == MD_POT_UNIT_UNKNOWN) {
                    MD_LOG_ERROR("POT: line %i: unrecognized unit '%.*s'", line_num, (int)val.len, val.ptr);
                    goto fail;
                }
                pot->unit = unit;
            } else {
                MD_LOG_INFO("POT: line %i: ignoring unknown directive '%.*s'", line_num, (int)key.len, key.ptr);
            }
            continue;
        }

        switch (section) {
        case SECTION_ENVIRONMENT: {
            // element x y z [fragment_name [fragment_id [atom_name]]]
            if (num_tok < 4) {
                MD_LOG_ERROR("POT: line %i: expected at least 4 fields in @environment, got %i", line_num, (int)num_tok);
                goto fail;
            }
            for (int i = 1; i < 4; ++i) {
                if (!is_float(tok[i])) {
                    MD_LOG_ERROR("POT: line %i: coordinate field %i is not a number: '%.*s'", line_num, i, (int)tok[i].len, tok[i].ptr);
                    goto fail;
                }
            }

            md_pot_site_t site = {0};
            pot_copy_symbol(site.element, tok[0], line_num);
            site.coord[0]   = parse_float(tok[1]);
            site.coord[1]   = parse_float(tok[2]);
            site.coord[2]   = parse_float(tok[3]);
            site.fragment_id = -1;

            if (num_tok > 4) {
                site.fragment_name = pot_intern(pot, tok[4], alloc);
            }
            if (num_tok > 5) {
                if (!is_int(tok[5])) {
                    MD_LOG_ERROR("POT: line %i: fragment number is not an integer: '%.*s'", line_num, (int)tok[5].len, tok[5].ptr);
                    goto fail;
                }
                site.fragment_id = (int)parse_int(tok[5]);
            }
            if (num_tok > 6) {
                site.atom_name = pot_intern(pot, tok[6], alloc);
            }
            if (num_tok > 7) {
                MD_LOG_INFO("POT: line %i: ignoring %i trailing field(s) in @environment", line_num, (int)(num_tok - 7));
            }

            md_array_push(pot->sites, site, alloc);
            break;
        }
        case SECTION_CHARGES: {
            // element charge [fragment_name]
            if (num_tok < 2) {
                MD_LOG_ERROR("POT: line %i: expected at least 2 fields in @charges, got %i", line_num, (int)num_tok);
                goto fail;
            }
            if (!is_float(tok[1])) {
                MD_LOG_ERROR("POT: line %i: charge is not a number: '%.*s'", line_num, (int)tok[1].len, tok[1].ptr);
                goto fail;
            }

            md_pot_charge_t charge = {0};
            pot_copy_symbol(charge.element, tok[0], line_num);
            charge.charge = parse_float(tok[1]);
            if (num_tok > 2) {
                charge.fragment_name = pot_intern(pot, tok[2], alloc);
            }
            if (num_tok > 3) {
                MD_LOG_INFO("POT: line %i: ignoring %i trailing field(s) in @charges", line_num, (int)(num_tok - 3));
            }

            md_array_push(pot->charges, charge, alloc);
            break;
        }
        case SECTION_POLARIZABILITIES: {
            // element a_xx a_xy a_xz a_yy a_yz a_zz [fragment_name]
            if (num_tok < 7) {
                MD_LOG_ERROR("POT: line %i: expected at least 7 fields in @polarizabilities, got %i", line_num, (int)num_tok);
                goto fail;
            }
            for (int i = 1; i < 7; ++i) {
                if (!is_float(tok[i])) {
                    MD_LOG_ERROR("POT: line %i: polarizability field %i is not a number: '%.*s'", line_num, i, (int)tok[i].len, tok[i].ptr);
                    goto fail;
                }
            }

            md_pot_polarizability_t pol = {0};
            pot_copy_symbol(pol.element, tok[0], line_num);
            for (int i = 0; i < 6; ++i) {
                pol.alpha[i] = parse_float(tok[i + 1]);
            }
            if (num_tok > 7) {
                pol.fragment_name = pot_intern(pot, tok[7], alloc);
            }
            if (num_tok > 8) {
                MD_LOG_INFO("POT: line %i: ignoring %i trailing field(s) in @polarizabilities", line_num, (int)(num_tok - 8));
            }

            md_array_push(pot->polarizabilities, pol, alloc);
            break;
        }
        default:
            ASSERT(false);
            break;
        }
    }

    if (section != SECTION_NONE) {
        MD_LOG_INFO("POT: file ended with an unterminated section (missing '@end')");
    }

    if (md_array_size(pot->sites) == 0) {
        MD_LOG_ERROR("POT: no sites were parsed, missing or empty @environment section");
        goto fail;
    }

    pot->num_sites            = md_array_size(pot->sites);
    pot->num_charges          = md_array_size(pot->charges);
    pot->num_polarizabilities = md_array_size(pot->polarizabilities);
    pot->num_strings          = md_array_size(pot->strings);

    return true;

fail:
    md_pot_free(pot, alloc);
    return false;
}

bool md_pot_parse_file(md_pot_t* pot, str_t filename, md_allocator_i* alloc) {
    if (!pot) {
        MD_LOG_ERROR("POT: pot object was NULL");
        return false;
    }
    if (!alloc) {
        MD_LOG_ERROR("POT: allocator was NULL");
        return false;
    }

    bool success = false;
    md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);

    str_t str = load_textfile(filename, temp_alloc);
    if (str_empty(str)) {
        MD_LOG_ERROR("POT: failed to load file '%.*s'", (int)filename.len, filename.ptr);
    } else {
        success = md_pot_parse_str(pot, str, alloc);
    }

    md_temp_end(temp_scope);
    return success;
}

void md_pot_free(md_pot_t* pot, md_allocator_i* alloc) {
    ASSERT(pot);
    ASSERT(alloc);

    md_array_free(pot->sites, alloc);
    md_array_free(pot->charges, alloc);
    md_array_free(pot->polarizabilities, alloc);

    for (size_t i = 0; i < md_array_size(pot->strings); ++i) {
        str_free(pot->strings[i], alloc);
    }
    md_array_free(pot->strings, alloc);

    MEMSET(pot, 0, sizeof(md_pot_t));
}
