#include <md_amber_energy.h>

#include <md_system.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_parse.h>
#include <core/md_str.h>
#include <core/md_unit.h>

#include <math.h>
#include <string.h>
#include <stdio.h>

// A table of records: one row per record, one column per term, NaN where a record lacks the term
typedef struct amber_term_t {
    md_str_t   name;       // as mdout writes it
    char    slug[64];
    md_array(double) values;
} amber_term_t;

typedef struct amber_table_t {
    md_array(double)       time;
    md_array(amber_term_t) terms;
    md_allocator_i*        alloc;
} amber_table_t;

static md_str_t amber_next_line(md_str_t* text) {
    size_t nl = 0;
    md_str_t line;
    if (!md_str_find_char(&nl, *text, '\n')) {
        line = *text;
        *text = (md_str_t){0};
    } else {
        line = md_str_substr(*text, 0, nl);
        *text = md_str_substr(*text, nl + 1, SIZE_MAX);
    }
    if (line.len && line.ptr[line.len - 1] == '\r') line.len -= 1;
    return line;
}

// Lowercase, every run of other characters one '_', none at the ends
static size_t amber_slug(char* out, size_t cap, md_str_t name) {
    size_t n = 0;
    bool sep = false;
    for (size_t i = 0; i < name.len && n + 2 < cap; ++i) {
        char c = name.ptr[i];
        if (c >= 'A' && c <= 'Z') c = (char)(c - 'A' + 'a');
        const bool alnum = (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9');
        if (alnum) {
            if (sep && n) out[n++] = '_';
            out[n++] = c;
            sep = false;
        } else {
            sep = true;
        }
    }
    out[n] = '\0';
    return n;
}

static size_t amber_term(amber_table_t* t, md_str_t name) {
    char slug[64];
    amber_slug(slug, sizeof(slug), name);
    for (size_t i = 0; i < md_array_size(t->terms); ++i) {
        if (strcmp(t->terms[i].slug, slug) == 0) return i;
    }
    amber_term_t term = { .name = md_str_copy(name, t->alloc) };
    md_memcpy(term.slug, slug, sizeof(slug));
    // Rows recorded before the term appeared lack it
    for (size_t r = 0; r < md_array_size(t->time); ++r) md_array_push(term.values, NAN, t->alloc);
    md_array_push(t->terms, term, t->alloc);
    return md_array_size(t->terms) - 1;
}

static void amber_begin_row(amber_table_t* t, double time) {
    md_array_push(t->time, time, t->alloc);
    for (size_t i = 0; i < md_array_size(t->terms); ++i) md_array_push(t->terms[i].values, NAN, t->alloc);
}

static void amber_set(amber_table_t* t, md_str_t name, double value) {
    if (md_array_size(t->time) == 0) return;
    const size_t i = amber_term(t, name);
    *md_array_last(t->terms[i].values) = value;
}

// ## mdout
//
//  NSTEP =      500   TIME(PS) =       0.500  TEMP(K) =   182.95  PRESS =     0.0
//  Etot   =       -16.7421  EKtot   =         5.4532  EPtot      =       -22.1954
//  BOND   =         2.4475  ANGLE   =         2.8681  DIHED      =         3.3755
//  1-4 NB =         0.5154  1-4 EEL =       -19.8735  VDWAALS    =        -0.2550
//  EELEC  =       -11.2733  EHBOND  =         0.0000  RESTRAINT  =         0.0000
//  ------------------------------------------------------------------------------
//
// A record runs from its NSTEP line to the dashes. Every "<name> = <value>" in it is a term; a name is what lies
// between the value before and the '='.

// 'region' > 1: the terms are those of that TI region, named "<name> (TI region <region>)"
static void amber_parse_pairs(amber_table_t* t, md_str_t line, bool first_line, double* out_time, int region) {
    size_t pos = 0;
    while (pos < line.len) {
        size_t eq = 0;
        if (!md_str_find_char(&eq, md_str_substr(line, pos, SIZE_MAX), '=')) break;
        eq += pos;
        const md_str_t name = md_str_trim(md_str_substr(line, pos, eq - pos));
        size_t v = eq + 1;
        while (v < line.len && line.ptr[v] == ' ') ++v;
        size_t e = v;
        while (e < line.len && line.ptr[e] != ' ') ++e;
        const md_str_t token = md_str_substr(line, v, e - v);
        double value = NAN;
        if (md_parse_f64(&value, token) != token.len || token.len == 0) value = NAN;   // ****** overflows
        pos = e;
        if (md_str_empty(name)) continue;
        if (first_line && md_str_eq_cstr_ignore_case(name, "NSTEP")) continue;
        if (first_line && md_str_eq_cstr_ignore_case(name, "TIME(PS)")) {
            *out_time = value;
            continue;
        }
        if (region > 1) {
            char buf[128];
            const int len = snprintf(buf, sizeof(buf), "%.*s (TI region %d)", (int)MD_MIN(name.len, 64), name.ptr, region);
            amber_set(t, (md_str_t){ buf, (size_t)MD_MAX(len, 0) }, value);
        } else {
            amber_set(t, name, value);
        }
    }
}

// A thermodynamic integration run states every record once per TI region, each opened by "| TI region <n>", with the
// same step and time and energies of their own (each end state's). They are one row: region 1's terms under their own
// names, the others' as "<name> (TI region <n>)". Read as rows of their own they would put every time on the axis twice.
static bool amber_parse_mdout(amber_table_t* t, md_str_t text, md_str_t path) {
    bool in_results = false;
    bool any_section = false;
    int ti_region = 0;
    while (text.len) {
        md_str_t line = amber_next_line(&text);
        const md_str_t trimmed = md_str_trim(line);
        if (md_str_begins_with(trimmed, MD_STR_LIT("| TI region"))) {
            int64_t region = 0;
            if (md_parse_i64(&region, md_str_trim(md_str_substr(trimmed, 11, SIZE_MAX))) > 0) ti_region = (int)region;
            continue;
        }
        // "   4.  RESULTS", "   5.  TIMINGS": the numbered sections of the log
        if (trimmed.len > 3 && trimmed.ptr[0] >= '1' && trimmed.ptr[0] <= '9' && trimmed.ptr[1] == '.' && trimmed.ptr[2] == ' ') {
            any_section = true;
            in_results = md_str_eq_cstr(md_str_trim(md_str_substr(trimmed, 2, SIZE_MAX)), "RESULTS");
            continue;
        }
        // The averages and fluctuations at the end of the results, and between the stages of a run of several; and in a
        // TI run the averages of DV/DL ("DV/DL, AVERAGES OVER 50000 STEPS"), which are a record of the last step to look at
        if (md_str_find_str(NULL, trimmed, MD_STR_LIT("A V E R A G E S")) || md_str_find_str(NULL, trimmed, MD_STR_LIT("F L U C T U A T I O N S")) ||
            md_str_find_str(NULL, trimmed, MD_STR_LIT("AVERAGES OVER"))) {
            ti_region = 0;
            // Skip their record, up to its dashes
            while (text.len) {
                const md_str_t l = md_str_trim(amber_next_line(&text));
                if (md_str_begins_with(l, MD_STR_LIT("-----"))) break;
            }
            continue;
        }
        if (any_section && !in_results) continue;
        if (!md_str_begins_with(trimmed, MD_STR_LIT("NSTEP")) || !md_str_find_char(NULL, trimmed, '=')) continue;

        // A record
        const int region = ti_region;
        ti_region = 0;
        if (region > 1 && md_array_size(t->time) > 0) {
            // Another TI region's record of the step just read, into its row, when it is that step
            double time = NAN;
            md_str_t probe = line;
            amber_parse_pairs(&(amber_table_t){ .alloc = t->alloc }, probe, true, &time, 0);
            if (time == *md_array_last(t->time)) {
                amber_parse_pairs(t, line, true, &time, region);
                while (text.len) {
                    line = amber_next_line(&text);
                    const md_str_t l = md_str_trim(line);
                    if (md_str_begins_with(l, MD_STR_LIT("-----")) || md_str_empty(l)) break;
                    amber_parse_pairs(t, line, false, &time, region);
                }
                continue;
            }
        }
        double time = NAN;
        const size_t row = md_array_size(t->time);
        amber_begin_row(t, NAN);
        amber_parse_pairs(t, line, true, &time, 0);
        while (text.len) {
            line = amber_next_line(&text);
            const md_str_t l = md_str_trim(line);
            if (md_str_begins_with(l, MD_STR_LIT("-----")) || md_str_empty(l)) break;
            amber_parse_pairs(t, line, false, &time, 0);
        }
        if (!isfinite(time)) {
            // Not a record of molecular dynamics (minimization writes no time): dropped
            md_array_pop(t->time);
            for (size_t i = 0; i < md_array_size(t->terms); ++i) md_array_pop(t->terms[i].values);
            continue;
        }
        t->time[row] = time;
    }
    if (md_array_size(t->time) == 0) {
        MD_LOG_ERROR("Amber energy: '" MD_STR_FMT "' holds no records of molecular dynamics (NSTEP = ... TIME(PS) = ...)", MD_STR_ARG(path));
        return false;
    }
    return true;
}

// ## mden
//
//     L0  Nsteps           time(ps)         Etot             EKinetic
//     L1  Temp             T_solute         T_solv           Pres_scal_solu
//     ...
//     L9  AV_permMoment    AV_indMoment     AV_totMoment     Density          dV/dlambda
//     L0        0  0.0000000000E+00  -.1000000000E+02  0.5000000000E+01
//     L1 0.3000000000E+03  ...
//
// The first ten lines name the columns, each record is ten lines of the same labels. The names are mapped onto those
// of mdout where it has the term.

static md_str_t amber_mden_name(md_str_t name) {
    static const char* const map[][2] = {
        { "Etot", "Etot" }, { "EKinetic", "EKtot" }, { "Temp", "TEMP(K)" }, { "volume", "VOLUME" }, { "Pressure", "PRESS" },
        { "EKComTot", "EKCMT" }, { "VIRIAL_tot", "VIRIAL" }, { "E_pot", "EPtot" }, { "E_vdw", "VDWAALS" }, { "E_el", "EELEC" },
        { "E_hbon", "EHBOND" }, { "E_bon", "BOND" }, { "E_angle", "ANGLE" }, { "E_dih", "DIHED" }, { "E_14vdw", "1-4 NB" },
        { "E_14el", "1-4 EEL" }, { "E_const", "RESTRAINT" }, { "E_pol", "EPOLZ" }, { "Density", "Density" }, { "dV/dlambda", "DV/DL" },
    };
    for (size_t i = 0; i < MD_COUNTOF(map); ++i) {
        if (md_str_eq_cstr(name, map[i][0])) return md_str_from_cstr(map[i][1]);
    }
    return name;
}

static size_t amber_tokens(md_str_t* out, size_t cap, md_str_t line) {
    size_t n = 0, i = 0;
    while (i < line.len && n < cap) {
        while (i < line.len && (line.ptr[i] == ' ' || line.ptr[i] == '\t')) ++i;
        size_t j = i;
        while (j < line.len && line.ptr[j] != ' ' && line.ptr[j] != '\t') ++j;
        if (j > i) out[n++] = (md_str_t){line.ptr + i, j - i};
        i = j;
    }
    return n;
}

static bool amber_parse_mden(amber_table_t* t, md_str_t text, md_str_t path) {
    enum { MAX_COLS = 16, MAX_LINES = 32 };
    md_str_t names[MAX_LINES][MAX_COLS];
    size_t num_names[MAX_LINES] = {0};
    size_t num_header = 0;

    // The header: labelled lines of names, until the first labelled line of numbers
    md_str_t rest = text;
    while (rest.len) {
        md_str_t peek = rest;
        const md_str_t line = amber_next_line(&peek);
        md_str_t tok[MAX_COLS + 1];
        const size_t n = amber_tokens(tok, MAX_COLS + 1, line);
        if (n < 2 || tok[0].len < 2 || tok[0].ptr[0] != 'L') { rest = peek; continue; }
        double v;
        if (md_parse_f64(&v, tok[1]) == tok[1].len) break;  // numbers: the first record
        if (num_header == MAX_LINES) return false;
        num_names[num_header] = n - 1;
        for (size_t k = 1; k < n && k <= MAX_COLS; ++k) names[num_header][k - 1] = tok[k];
        num_header += 1;
        rest = peek;
    }
    if (num_header == 0 || num_names[0] < 2 || !md_str_eq_cstr(names[0][0], "Nsteps")) {
        MD_LOG_ERROR("Amber energy: '" MD_STR_FMT "' does not start with the column names of an mden file", MD_STR_ARG(path));
        return false;
    }

    // The records: num_header lines each, L0 first
    size_t line_in_record = 0;
    while (rest.len) {
        const md_str_t line = amber_next_line(&rest);
        md_str_t tok[MAX_COLS + 1];
        const size_t n = amber_tokens(tok, MAX_COLS + 1, line);
        if (n == 0) continue;
        if (tok[0].len < 2 || tok[0].ptr[0] != 'L') break;
        int64_t label = -1;
        md_parse_i64(&label, md_str_substr(tok[0], 1, SIZE_MAX));
        if (label == 0) {
            line_in_record = 0;
        }
        if (label != (int64_t)line_in_record || (size_t)label >= num_header || n - 1 != num_names[label]) {
            MD_LOG_ERROR("Amber energy: '" MD_STR_FMT "' has a record line out of order, or of another number of values", MD_STR_ARG(path));
            return false;
        }
        for (size_t k = 1; k < n; ++k) {
            double v = NAN;
            if (md_parse_f64(&v, tok[k]) != tok[k].len) v = NAN;
            const md_str_t name = names[label][k - 1];
            if (label == 0 && k == 1) {
                continue;   // Nsteps
            }
            if (label == 0 && k == 2) {
                amber_begin_row(t, v);   // time(ps)
                continue;
            }
            amber_set(t, amber_mden_name(name), v);
        }
        line_in_record += 1;
    }
    if (md_array_size(t->time) == 0) {
        MD_LOG_ERROR("Amber energy: '" MD_STR_FMT "' holds no records", MD_STR_ARG(path));
        return false;
    }
    return true;
}

// ## Units by name

static md_unit_t amber_term_unit(const char* slug) {
    const md_unit_t kcal_per_mol = md_unit_div(md_unit_scl(md_unit_calorie(), 1.0e3), md_unit_mole());
    static const char* const kelvin[] = { "temp_k", "t_solute", "t_solv", "temp0" };
    static const char* const bar[] = { "press", "pres_x", "pres_y", "pres_z" };
    static const char* const angstrom[] = { "boxx", "boxy", "boxz" };
    static const char* const debye[] = { "av_permmoment", "av_indmoment", "av_totmoment" };
    static const char* const none[] = { "pres_scal_solu", "pres_scal_solv", "ewald_error_estimate", "surften", "solvent_accessible_surface" };
    for (size_t i = 0; i < MD_COUNTOF(kelvin); ++i)   if (!strcmp(slug, kelvin[i]))   return md_unit_kelvin();
    for (size_t i = 0; i < MD_COUNTOF(bar); ++i)      if (!strcmp(slug, bar[i]))      return md_unit_bar();
    for (size_t i = 0; i < MD_COUNTOF(angstrom); ++i) if (!strcmp(slug, angstrom[i])) return md_unit_angstrom();
    for (size_t i = 0; i < MD_COUNTOF(debye); ++i)    if (!strcmp(slug, debye[i]))    return md_unit_debye();
    for (size_t i = 0; i < MD_COUNTOF(none); ++i)     if (!strcmp(slug, none[i]))     return md_unit_none();
    if (!strcmp(slug, "volume")) return md_unit_pow(md_unit_angstrom(), 3);
    if (!strcmp(slug, "density")) {
        // g/cm^3
        return md_unit_div(md_unit_scl(md_unit_kilogram(), 1.0e-3), md_unit_pow(md_unit_scl(md_unit_meter(), 1.0e-2), 3));
    }
    return kcal_per_mol;
}

// ## Publishing

bool md_amber_energy_file_is_mdout(md_str_t filename) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) return false;
    char buf[8192];
    const size_t len = md_file_read(file, buf, sizeof(buf));
    md_file_close(&file);
    const md_str_t head = {buf, len};
    // The file assignments at the top of every mdout ("|  MDIN: mdin"), and the program that wrote it
    return (md_str_find_str(NULL, head, MD_STR_LIT("MDIN:")) || md_str_find_str(NULL, head, MD_STR_LIT("MDOUT:"))) &&
           (md_str_find_str(NULL, head, MD_STR_LIT("PMEMD")) || md_str_find_str(NULL, head, MD_STR_LIT("SANDER")) || md_str_find_str(NULL, head, MD_STR_LIT("Amber")));
}

bool md_amber_energy_file_is_mden(md_str_t filename) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) return false;
    char buf[256];
    const size_t len = md_file_read(file, buf, sizeof(buf));
    md_file_close(&file);
    md_str_t tok[3];
    const size_t n = amber_tokens(tok, 3, amber_next_line(&(md_str_t){buf, len}));
    return n >= 2 && md_str_eq_cstr(tok[0], "L0") && md_str_eq_cstr(tok[1], "Nsteps");
}

static md_str_t amber_join(char* buf, size_t cap, md_str_t a, const char* b) {
    const size_t b_len = strlen(b);
    if (a.len + 1 + b_len + 1 > cap) return (md_str_t){0};
    md_memcpy(buf, a.ptr, a.len);
    buf[a.len] = '/';
    md_memcpy(buf + a.len + 1, b, b_len);
    buf[a.len + 1 + b_len] = '\0';
    return (md_str_t){buf, a.len + 1 + b_len};
}

bool md_amber_energy_system_supplement_from_file(md_system_t* sys, md_str_t filename, md_str_t run) {
    MD_ASSERT(sys);
    if (md_str_empty(run)) {
        MD_LOG_ERROR("Amber energy: no run to publish the energies under; load the trajectory first");
        return false;
    }
    md_attributes_t* attributes = &sys->attributes;
    if (!attributes->alloc) attributes->alloc = sys->alloc;

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    bool result = false;
    md_attributes_t probe = { .alloc = temp_alloc };

    amber_table_t t = { .alloc = temp_alloc };
    const md_str_t text = md_load_textfile(filename, temp_alloc);
    if (md_str_empty(text)) {
        MD_LOG_ERROR("Amber energy: cannot read '" MD_STR_FMT "'", MD_STR_ARG(filename));
        goto done;
    }
    const bool mden = md_amber_energy_file_is_mden(filename);
    if (!(mden ? amber_parse_mden(&t, text, filename) : amber_parse_mdout(&t, text, filename))) goto done;

    const size_t R = md_array_size(t.time);
    for (size_t r = 1; r < R; ++r) {
        if (!(t.time[r] >= t.time[r - 1])) {
            MD_LOG_ERROR("Amber energy: the times of '" MD_STR_FMT "' decrease at record %zu (%g ps after %g ps); a restarted run is to be loaded file by file",
                MD_STR_ARG(filename), r, t.time[r], t.time[r - 1]);
            goto done;
        }
    }

    // The file's own axis, against a scratch table first: whether it covers the run is a question about the values
    const md_attribute_desc_t time_desc = {
        .path   = MD_STR_INIT("time"),
        .format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)R } },
        .flags  = MD_ATTRIBUTE_FLAG_TEMPORAL,
        .unit   = md_unit_picosecond(),
        .label  = MD_STR_INIT("Time"),
        .data   = t.time,
        .byte_size = R * sizeof(double),
    };
    const md_attribute_id_t probe_id = md_attributes_create(&probe, &time_desc);
    if (!probe_id) goto done;

    char path_buf[512], group_buf[512];
    const md_attribute_t* run_axis = md_attributes_find(attributes, amber_join(path_buf, sizeof(path_buf), run, "time"));
    if (run_axis && md_attributes_axis(attributes, run_axis) == run_axis) {
        const md_attribute_t* axis = md_attributes_get(&probe, probe_id);
        for (size_t f = 0; f < run_axis->format.shape[0]; ++f) {
            size_t row;
            if (!md_attribute_axis_map(&row, run_axis, f, axis)) {
                double time = 0;
                md_attribute_extract_f64(&time, 1, run_axis, md_attribute_slice_1((uint32_t)f), md_unit_picosecond());
                MD_LOG_ERROR("Amber energy: frame %zu of '" MD_STR_FMT "' (%g ps) has no record in '" MD_STR_FMT "' (%g - %g ps); the file belongs to another run, or was written at another interval",
                    f, MD_STR_ARG(run), time, MD_STR_ARG(filename), t.time[0], t.time[R - 1]);
                goto done;
            }
        }
    }

    const md_str_t group = amber_join(group_buf, sizeof(group_buf), run, "amber_energy");
    if (md_str_empty(group)) goto done;
    md_attributes_remove_prefix(attributes, group);

    {
        md_attribute_desc_t desc = time_desc;
        desc.path = amber_join(path_buf, sizeof(path_buf), group, "time");
        if (!md_attributes_create(attributes, &desc)) goto fail;
    }
    size_t published = 0;
    for (size_t i = 0; i < md_array_size(t.terms); ++i) {
        const amber_term_t* term = &t.terms[i];
        bool complete = term->slug[0] != '\0' && strcmp(term->slug, "time") && strcmp(term->slug, "source");
        for (size_t r = 0; complete && r < R; ++r) complete = isfinite(term->values[r]);
        if (!complete) {
            MD_LOG_INFO("Amber energy: '" MD_STR_FMT "' is left out, it is not in every record of '" MD_STR_FMT "'", MD_STR_ARG(term->name), MD_STR_ARG(filename));
            continue;
        }
        const md_attribute_desc_t desc = {
            .path   = amber_join(path_buf, sizeof(path_buf), group, term->slug),
            .format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)R } },
            .flags  = MD_ATTRIBUTE_FLAG_TEMPORAL,
            .unit   = amber_term_unit(term->slug),
            .label  = term->name,
            .data   = term->values,
            .byte_size = R * sizeof(double),
        };
        if (md_str_empty(desc.path) || !md_attributes_create(attributes, &desc)) goto fail;
        published += 1;
    }
    {
        char canonical[4096];
        const size_t len = md_path_write_canonical(canonical, sizeof(canonical), filename);
        const md_str_t source = len ? md_str_copy((md_str_t){canonical, len}, attributes->alloc) : md_str_copy(filename, attributes->alloc);
        const md_attribute_desc_t desc = {
            .path   = amber_join(path_buf, sizeof(path_buf), group, "source"),
            .format = { .type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 0 },
            .unit   = md_unit_none(),
            .label  = MD_STR_INIT("Source"),
            .data   = &source,
            .byte_size = sizeof(md_str_t),
        };
        if (!md_attributes_create(attributes, &desc)) goto fail;
    }
    MD_LOG_INFO("Amber energy: %zu terms of %zu records from '" MD_STR_FMT "'", published, R, MD_STR_ARG(filename));
    result = true;
    goto done;

fail:
    MD_LOG_ERROR("Amber energy: failed to publish the energies under '" MD_STR_FMT "'", MD_STR_ARG(group));
    md_attributes_remove_prefix(attributes, group);

done:
    md_attributes_free(&probe);
    md_temp_end(temp);
    return result;
}
