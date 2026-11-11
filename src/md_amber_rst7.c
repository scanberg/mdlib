#include <md_amber_rst7.h>

#include <md_system.h>
#include <md_unitcell.h>

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_file.h>
#include <core/md_parse.h>
#include <core/md_str.h>

#include <math.h>
#include <string.h>

#define RST7_FIELD 12

static md_str_t rst7_next_line(md_str_t* text) {
    md_str_t line;
    md_str_cut(&line, text, *text, md_str_find_char(*text, '\n'));
    if (line.len && line.ptr[line.len - 1] == '\r') line.len -= 1;
    return line;
}

// NATOM and, optionally, the time on the second line: "   2304  0.1000000E+04"
static bool rst7_header(md_str_t line, size_t* out_atoms) {
    line = md_str_trim(line);
    int64_t n = 0;
    const size_t len = md_parse_i64(&n, line);
    if (!len || n <= 0 || n > INT32_MAX) return false;
    const md_str_t rest = md_str_trim(md_str_substr(line, len, SIZE_MAX));
    if (!md_str_empty(rest)) {
        double t;
        if (md_parse_f64(&t, rest) == 0) return false;
    }
    *out_atoms = (size_t)n;
    return true;
}

// The values of a line of 12 wide fields. Returns the count, or SIZE_MAX when a field is not a number (an overflowed
// field is written as asterisks).
static size_t rst7_fields(double* out, size_t cap, md_str_t line) {
    size_t n = 0;
    for (size_t beg = 0; beg < line.len; beg += RST7_FIELD) {
        const md_str_t field = md_str_trim(md_str_substr(line, beg, RST7_FIELD));
        if (md_str_empty(field)) continue;
        if (n == cap) return SIZE_MAX;
        if (md_parse_f64(&out[n], field) != field.len || !isfinite(out[n])) return SIZE_MAX;
        n += 1;
    }
    return n;
}

static bool rst7_peek(size_t* out_atoms, md_str_t filename) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) return false;
    char buf[1024];
    const size_t len = md_file_read(file, buf, sizeof(buf));
    md_file_close(&file);
    md_str_t text = {buf, len};
    rst7_next_line(&text);                                  // title
    if (!rst7_header(rst7_next_line(&text), out_atoms)) return false;
    if (!text.len) return false;
    // The first line of coordinates: whole fields of 12, at most 6 of them, and numbers
    const md_str_t line = rst7_next_line(&text);
    if (!text.ptr) return false;                            // the line may have been cut by the buffer
    double v[6];
    const size_t n = rst7_fields(v, 6, line);
    return n != SIZE_MAX && n >= 3 && line.len <= 6 * RST7_FIELD + 2;
}

bool md_amber_rst7_file_is_restart(md_str_t filename) {
    size_t n;
    return rst7_peek(&n, filename);
}

size_t md_amber_rst7_file_atom_count(md_str_t filename) {
    size_t n = 0;
    return rst7_peek(&n, filename) ? n : 0;
}

bool md_amber_rst7_state_init_from_file(md_system_state_t* state, md_str_t filename, size_t num_atoms) {
    MD_ASSERT(state);
    if (!state->alloc) {
        MD_LOG_ERROR("Amber restart: the state has no allocator");
        return false;
    }
    md_temp_scope_t temp = md_temp_begin_avoid(state->alloc);
    md_allocator_i* alloc = md_temp_allocator(temp);
    bool result = false;

    md_str_t text = md_load_textfile(filename, alloc);
    size_t N = 0;
    if (md_str_empty(text)) {
        MD_LOG_ERROR("Amber restart: cannot read '" MD_STR_FMT "'", MD_STR_ARG(filename));
        goto done;
    }
    rst7_next_line(&text);
    if (!rst7_header(rst7_next_line(&text), &N)) {
        MD_LOG_ERROR("Amber restart: '" MD_STR_FMT "' has no atom count on its second line", MD_STR_ARG(filename));
        goto done;
    }
    if (num_atoms && N != num_atoms) {
        MD_LOG_ERROR("Amber restart: '" MD_STR_FMT "' holds %zu atoms, the system %zu", MD_STR_ARG(filename), N, num_atoms);
        goto done;
    }

    // Two atoms a line: the coordinates take ceil(3N / 6) lines, velocities as many again, the box one
    const size_t lines_per_block = (3 * N + 5) / 6;
    float* xyz = md_temp_alloc_array(temp, float, 3 * N);
    size_t read = 0;
    double v[6];
    for (size_t l = 0; l < lines_per_block; ++l) {
        const md_str_t line = rst7_next_line(&text);
        const size_t want = MD_MIN((size_t)6, 3 * N - read);
        if (rst7_fields(v, 6, line) != want) {
            MD_LOG_ERROR("Amber restart: line %zu of '" MD_STR_FMT "' does not hold %zu coordinates", l + 3, MD_STR_ARG(filename), want);
            goto done;
        }
        for (size_t k = 0; k < want; ++k) xyz[read++] = (float)v[k];
    }

    // What follows: nothing, a box line, velocities, or velocities and a box line
    size_t rest = 0;
    md_str_t lines[2] = {0};   // the last two non blank lines
    for (md_str_t t = text; t.len;) {
        const md_str_t line = rst7_next_line(&t);
        if (md_str_empty(md_str_trim(line))) continue;
        lines[0] = lines[1];
        lines[1] = line;
        rest += 1;
    }
    md_unitcell_t cell = md_unitcell_none();
    const bool has_box = rest == 1 || rest == lines_per_block + 1;
    if (!(rest == 0 || rest == 1 || rest == lines_per_block || rest == lines_per_block + 1)) {
        MD_LOG_ERROR("Amber restart: '" MD_STR_FMT "' has %zu lines after its coordinates, which are neither velocities nor a box", MD_STR_ARG(filename), rest);
        goto done;
    }
    if (has_box) {
        const size_t n = rst7_fields(v, 6, lines[1]);
        if (n == 3 || n == 6) {
            if (n == 3) { v[3] = v[4] = v[5] = 90.0; }
            if (v[0] > 0 && v[1] > 0 && v[2] > 0) {
                cell = md_unitcell_from_extent_and_angles(v[0], v[1], v[2], v[3], v[4], v[5]);
            }
        } else {
            MD_LOG_ERROR("Amber restart: the box of '" MD_STR_FMT "' is not 3 or 6 numbers", MD_STR_ARG(filename));
            goto done;
        }
    }

    if (!md_system_state_init(state, N)) goto done;
    md_memcpy(state->xyz, xyz, 3 * N * sizeof(float));
    state->num_atoms = N;
    state->unitcell = cell;
    result = true;

done:
    md_temp_end(temp);
    return result;
}
