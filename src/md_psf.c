#include <md_psf.h>

#include <md_system.h>
#include <md_topology_build.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_parse.h>
#include <core/md_str.h>

#include <math.h>
#include <string.h>

static str_t psf_next_line(str_t* text) {
    size_t nl = 0;
    str_t line;
    if (!str_find_char(&nl, *text, '\n')) {
        line = *text;
        *text = (str_t){0};
    } else {
        line = str_substr(*text, 0, nl);
        *text = str_substr(*text, nl + 1, SIZE_MAX);
    }
    if (line.len && line.ptr[line.len - 1] == '\r') line.len -= 1;
    return line;
}

static size_t psf_tokens(str_t* out, size_t cap, str_t line) {
    size_t n = 0, i = 0;
    while (i < line.len && n < cap) {
        while (i < line.len && (line.ptr[i] == ' ' || line.ptr[i] == '\t')) ++i;
        size_t j = i;
        while (j < line.len && line.ptr[j] != ' ' && line.ptr[j] != '\t') ++j;
        if (j > i) out[n++] = (str_t){line.ptr + i, j - i};
        i = j;
    }
    return n;
}

// A section header: "    3341 !NATOM", "    3365 !NBOND: bonds", "       2       0 !NGRP NST2". The first count, and
// the name after the '!' up to a ':' or a space.
static bool psf_section(str_t line, int64_t* out_count, str_t* out_name) {
    size_t bang = 0;
    if (!str_find_char(&bang, line, '!')) return false;
    const str_t head = str_trim(str_substr(line, 0, bang));
    int64_t n = 0;
    if (head.len == 0 || md_parse_i64(&n, head) == 0 || n < 0) return false;
    str_t name = str_substr(line, bang + 1, SIZE_MAX);
    size_t end = 0;
    while (end < name.len && name.ptr[end] != ':' && name.ptr[end] != ' ') ++end;
    *out_count = n;
    *out_name = str_substr(name, 0, end);
    return true;
}

static bool psf_is_header(str_t text) {
    // "PSF", possibly with flags; some writers put blank lines first
    while (text.len) {
        const str_t line = str_trim(psf_next_line(&text));
        if (str_empty(line)) continue;
        return str_begins_with(line, STR_LIT("PSF")) && (line.len == 3 || line.ptr[3] == ' ' || line.ptr[3] == '\t');
    }
    return false;
}

bool md_psf_file_is_psf(str_t filename) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) return false;
    char buf[512];
    const size_t len = md_file_read(file, buf, sizeof(buf));
    md_file_close(&file);
    return psf_is_header((str_t){buf, len});
}

size_t md_psf_file_atom_count(str_t filename) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) return 0;
    // The title comes first and is short; !NATOM follows within the first block
    const size_t cap = 64 * 1024;
    char* buf = md_alloc(md_get_heap_allocator(), cap);
    const size_t len = md_file_read(file, buf, cap);
    md_file_close(&file);
    size_t result = 0;
    str_t text = {buf, len};
    if (psf_is_header(text)) {
        while (text.len) {
            int64_t n;
            str_t name;
            if (psf_section(psf_next_line(&text), &n, &name) && str_eq_cstr(name, "NATOM")) {
                result = (size_t)n;
                break;
            }
        }
    }
    md_free(md_get_heap_allocator(), buf, cap);
    return result;
}

// The header's flags: "PSF EXT CMAP XPLOR", "PSF NAMD"
static bool psf_has_flag(str_t text, const char* flag) {
    while (text.len) {
        const str_t line = str_trim(psf_next_line(&text));
        if (str_empty(line)) continue;
        str_t tok[16];
        const size_t n = psf_tokens(tok, ARRAY_SIZE(tok), line);
        for (size_t i = 1; i < n; ++i) {
            if (str_eq_cstr(tok[i], flag)) return true;
        }
        return false;
    }
    return false;
}

// The fields of an atom line: id segment residue resname name, then type charge mass and whatever a flavour adds
typedef struct psf_atom_fields_t {
    str_t id, seg, resid, resname, name, type, charge, mass;
} psf_atom_fields_t;

static bool psf_blank(str_t line, size_t col) {
    return col >= line.len || line.ptr[col] == ' ' || line.ptr[col] == '\t';
}

// The atom line as CHARMM writes it, in columns: (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,...) or with EXT (I10,1X,A8,1X,A8,1X,
// A8,1X,A8,1X,...), which VMD and ParmEd keep to as well. The type, charge and mass after the name are taken as
// words: the type is 4 or 6 wide and the numbers are in different columns from writer to writer, but none of them
// is ever blank. False when the line is not laid out so (a blank between the columns, the id a number).
static bool psf_atom_columns(psf_atom_fields_t* f, str_t line, bool ext) {
    const size_t id_w = ext ? 10 : 8;
    const size_t w    = ext ? 8 : 4;
    const size_t name_end = id_w + 4 * (w + 1);
    if (line.len <= name_end) return false;
    for (size_t k = 0; k < 4; ++k) {
        if (!psf_blank(line, id_w + k * (w + 1))) return false;
    }
    if (!psf_blank(line, name_end)) return false;
    f->id      = str_trim(str_substr(line, 0, id_w));
    f->seg     = str_trim(str_substr(line, id_w + 1, w));
    f->resid   = str_trim(str_substr(line, id_w + 1 + (w + 1), w));
    f->resname = str_trim(str_substr(line, id_w + 1 + 2 * (w + 1), w));
    f->name    = str_trim(str_substr(line, id_w + 1 + 3 * (w + 1), w));
    str_t tok[3];
    if (psf_tokens(tok, 3, str_substr(line, name_end, SIZE_MAX)) != 3) return false;
    f->type = tok[0]; f->charge = tok[1]; f->mass = tok[2];
    int64_t id;
    return !str_empty(f->id) && md_parse_i64(&id, f->id) == f->id.len && !str_empty(f->resid) && !str_empty(f->resname) && !str_empty(f->name);
}

// The atom line as words: NAMD's format, and any writer that does not keep to the columns. A blank segment cannot
// be told apart here, so the columns are tried first.
static bool psf_atom_words(psf_atom_fields_t* f, str_t line) {
    str_t tok[8];
    if (psf_tokens(tok, ARRAY_SIZE(tok), line) < 8) return false;
    f->id = tok[0]; f->seg = tok[1]; f->resid = tok[2]; f->resname = tok[3];
    f->name = tok[4]; f->type = tok[5]; f->charge = tok[6]; f->mass = tok[7];
    return true;
}

typedef struct psf_data_t {
    size_t num_atoms;
    str_t*  atom_name;
    str_t*  atom_type;
    double* charge;
    double* mass;
    md_array(str_t)    res_name;
    md_array(uint32_t) res_offset;
    md_array(int64_t)  res_number;
    md_array(str_t)    res_chain;
    md_array(md_atom_pair_t) bonds;
    str_t title;
} psf_data_t;

static bool psf_read(psf_data_t* d, str_t path, md_allocator_i* alloc) {
    MEMSET(d, 0, sizeof(*d));
    str_t text = load_textfile(path, alloc);
    if (!psf_is_header(text)) {
        MD_LOG_ERROR("PSF: '" STR_FMT "' does not start with 'PSF'", STR_ARG(path));
        return false;
    }
    // EXT widens the columns of the atoms; NAMD's are words of any width
    const bool ext  = psf_has_flag(text, "EXT");
    const bool namd = psf_has_flag(text, "NAMD");
    bool have_atoms = false;
    while (text.len) {
        const str_t line = psf_next_line(&text);
        int64_t count;
        str_t name;
        if (!psf_section(line, &count, &name)) continue;

        if (str_eq_cstr(name, "NTITLE")) {
            const char* beg = text.ptr;
            for (int64_t i = 0; i < count && text.len; ++i) psf_next_line(&text);
            if (beg && text.ptr) d->title = str_trim((str_t){beg, (size_t)(text.ptr - beg)});
        } else if (str_eq_cstr(name, "NATOM")) {
            if (count <= 0 || count > INT32_MAX) {
                MD_LOG_ERROR("PSF: '" STR_FMT "' has %lld atoms", STR_ARG(path), (long long)count);
                return false;
            }
            const size_t N = (size_t)count;
            d->num_atoms = N;
            d->atom_name = md_alloc(alloc, N * sizeof(str_t));
            d->atom_type = md_alloc(alloc, N * sizeof(str_t));
            d->charge    = md_alloc(alloc, N * sizeof(double));
            d->mass      = md_alloc(alloc, N * sizeof(double));
            str_t prev_seg = {0}, prev_resid = {0}, prev_resname = {0};
            for (size_t i = 0; i < N; ++i) {
                if (!text.len) {
                    MD_LOG_ERROR("PSF: '" STR_FMT "' ends within its atoms", STR_ARG(path));
                    return false;
                }
                const str_t line = psf_next_line(&text);
                psf_atom_fields_t f;
                int64_t id = 0;
                const bool fields = (!namd && psf_atom_columns(&f, line, ext)) || psf_atom_words(&f, line);
                if (!fields || md_parse_i64(&id, f.id) != f.id.len ||
                    md_parse_f64(&d->charge[i], f.charge) != f.charge.len || md_parse_f64(&d->mass[i], f.mass) != f.mass.len) {
                    MD_LOG_ERROR("PSF: atom %zu of '" STR_FMT "' is not 'id segment residue resname name type charge mass'", i + 1, STR_ARG(path));
                    return false;
                }
                if (id != (int64_t)i + 1) {
                    MD_LOG_ERROR("PSF: the atoms of '" STR_FMT "' are not numbered 1, 2, 3, ... (atom %zu is %lld)", STR_ARG(path), i + 1, (long long)id);
                    return false;
                }
                const str_t seg = f.seg, resid = f.resid, resname = f.resname;
                if (i == 0 || !str_eq(seg, prev_seg) || !str_eq(resid, prev_resid) || !str_eq(resname, prev_resname)) {
                    int64_t number = 0;
                    md_parse_i64(&number, resid);   // an insertion code may follow: 52A
                    md_array_push(d->res_offset, (uint32_t)i, alloc);
                    md_array_push(d->res_name, resname, alloc);
                    md_array_push(d->res_number, number, alloc);
                    md_array_push(d->res_chain, seg, alloc);
                    prev_seg = seg; prev_resid = resid; prev_resname = resname;
                }
                d->atom_name[i] = f.name;
                d->atom_type[i] = f.type;
            }
            md_array_push(d->res_offset, (uint32_t)N, alloc);
            have_atoms = true;
        } else if (str_eq_cstr(name, "NBOND")) {
            if (!have_atoms) break;
            md_array_ensure(d->bonds, (size_t)count, alloc);
            int64_t pending = -1;
            int64_t read = 0;
            while (read < 2 * count && text.len) {
                str_t tok[16];
                const size_t n = psf_tokens(tok, ARRAY_SIZE(tok), psf_next_line(&text));
                for (size_t k = 0; k < n; ++k) {
                    int64_t v = 0;
                    if (md_parse_i64(&v, tok[k]) != tok[k].len || v < 1 || v > (int64_t)d->num_atoms) {
                        MD_LOG_ERROR("PSF: '" STR_FMT "' has a bond to an atom that is not there", STR_ARG(path));
                        return false;
                    }
                    read += 1;
                    if (pending < 0) {
                        pending = v - 1;
                    } else {
                        if (pending != v - 1) {
                            const md_atom_pair_t pair = {{ (md_atom_idx_t)pending, (md_atom_idx_t)(v - 1) }};
                            md_array_push(d->bonds, pair, alloc);
                        }
                        pending = -1;
                    }
                }
            }
            break;  // nothing after the bonds is read
        }
    }
    if (!have_atoms) {
        MD_LOG_ERROR("PSF: '" STR_FMT "' has no !NATOM section", STR_ARG(path));
        return false;
    }
    return true;
}

static md_topology_desc_t psf_desc(const psf_data_t* d) {
    return (md_topology_desc_t){
        .format       = "PSF",
        .num_atoms    = d->num_atoms,
        .atom_name    = d->atom_name,
        .atom_type    = d->atom_type,
        .charge       = d->charge,
        .mass         = d->mass,
        .num_residues = md_array_size(d->res_name),
        .res_name     = d->res_name,
        .res_offset   = d->res_offset,
        .res_number   = d->res_number,
        .res_chain    = d->res_chain,
        .num_bonds    = md_array_size(d->bonds),
        .bonds        = d->bonds,
        .title        = d->title,
    };
}

bool md_psf_system_init_from_file(md_system_t* sys, md_system_state_t* state, str_t filename) {
    ASSERT(sys);
    ASSERT(state);
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    psf_data_t d;
    bool result = psf_read(&d, filename, md_temp_allocator(temp));
    if (result) {
        const md_topology_desc_t desc = psf_desc(&d);
        result = md_topology_system_init(sys, state, &desc);
    }
    md_temp_end(temp);
    return result;
}

bool md_psf_system_supplement_from_file(md_system_t* sys, str_t filename) {
    ASSERT(sys);
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    psf_data_t d;
    bool result = psf_read(&d, filename, md_temp_allocator(temp));
    if (result) {
        const md_topology_desc_t desc = psf_desc(&d);
        result = md_topology_system_supplement(sys, &desc, filename);
    }
    md_temp_end(temp);
    return result;
}
