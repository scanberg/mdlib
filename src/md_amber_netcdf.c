#include <md_amber_netcdf.h>

#include <md_system.h>
#include <md_unitcell.h>
#include <md_xdr.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include <core/md_unit.h>

#include <math.h>
#include <string.h>

// ## THE CLASSIC NETCDF FORMAT
//
// https://docs.unidata.ucar.edu/netcdf-c/current/file_format_specifications.html
//
//     file     = header data
//     header   = magic numrecs dim_list gatt_list var_list
//     magic    = 'C' 'D' 'F' version            1: CDF-1, 2: CDF-2 (64 bit offsets), 5: CDF-5 (64 bit data)
//     dim_list = ABSENT | NC_DIMENSION count dim*           dim = name length       (length 0: the record dimension)
//     att_list = ABSENT | NC_ATTRIBUTE count attr*          attr = name type count values (padded to 4)
//     var_list = ABSENT | NC_VARIABLE count var*            var = name rank dimid* att_list type vsize begin
//
// Everything is big endian (XDR). Counts, lengths and dimension ids are 4 bytes, 8 in CDF-5; begin is an
// offset of 4 bytes in CDF-1, 8 otherwise. A variable without the record dimension is one array at begin. One
// with it (its first dimension) has one slab per record, and the records follow each other after the fixed
// variables, every record holding a slab of each record variable in turn: slab r of a variable is at
// begin + r * recsize. A frame of an Amber trajectory is therefore one read at a known offset.

enum {
    CDF_BYTE = 1, CDF_CHAR, CDF_SHORT, CDF_INT, CDF_FLOAT, CDF_DOUBLE,
    CDF_UBYTE, CDF_USHORT, CDF_UINT, CDF_INT64, CDF_UINT64,
};

enum {
    CDF_TAG_DIMENSION = 0x0A,
    CDF_TAG_VARIABLE  = 0x0B,
    CDF_TAG_ATTRIBUTE = 0x0C,
};

#define CDF_MAX_RANK 8
#define CDF_MAX_HEADER_BYTES (256ull << 20)

typedef struct cdf_dim_t {
    str_t    name;
    uint64_t len;           // 0 for the record dimension
} cdf_dim_t;

typedef struct cdf_att_t {
    str_t          name;
    uint32_t       type;
    uint64_t       count;
    const uint8_t* data;    // into the header buffer
} cdf_att_t;

typedef struct cdf_var_t {
    str_t      name;
    uint32_t   type;
    uint32_t   rank;
    uint32_t   dim[CDF_MAX_RANK];
    cdf_att_t* att;
    size_t     num_att;
    uint64_t   begin;
    uint64_t   slab;        // bytes of one record (a record variable) or of all of it, unpadded
    bool       record;
} cdf_var_t;

typedef struct cdf_t {
    uint32_t   version;
    uint64_t   numrecs;
    bool       streaming;   // numrecs unknown: the writer did not say, count the records from the size
    cdf_dim_t* dim;
    size_t     num_dim;
    cdf_att_t* gatt;
    size_t     num_gatt;
    cdf_var_t* var;
    size_t     num_var;
    uint64_t   recsize;
    uint64_t   file_size;
} cdf_t;

static size_t cdf_type_size(uint32_t type) {
    switch (type) {
    case CDF_BYTE: case CDF_CHAR: case CDF_UBYTE:   return 1;
    case CDF_SHORT: case CDF_USHORT:                return 2;
    case CDF_INT: case CDF_UINT: case CDF_FLOAT:    return 4;
    case CDF_DOUBLE: case CDF_INT64: case CDF_UINT64: return 8;
    default: return 0;
    }
}

typedef struct cdf_cursor_t {
    md_xdr_t xdr;
    uint32_t version;
} cdf_cursor_t;

// A count, a length or a dimension id: 8 bytes in CDF-5, 4 otherwise
static uint64_t cdf_read_count(cdf_cursor_t* c) {
    if (c->version == 5) {
        uint64_t v = 0;
        md_xdr_read_u64(&c->xdr, &v);
        return v;
    }
    uint32_t v = 0;
    md_xdr_read_u32(&c->xdr, &v);
    return v;
}

static uint64_t cdf_read_offset(cdf_cursor_t* c) {
    if (c->version == 1) {
        uint32_t v = 0;
        md_xdr_read_u32(&c->xdr, &v);
        return v;
    }
    uint64_t v = 0;
    md_xdr_read_u64(&c->xdr, &v);
    return v;
}

static str_t cdf_read_name(cdf_cursor_t* c) {
    const uint64_t len = cdf_read_count(c);
    const uint8_t* ptr = NULL;
    if (len > 4096 || !md_xdr_read_opaque(&c->xdr, &ptr, (size_t)len)) {
        c->xdr.error = true;
        return (str_t){0};
    }
    return (str_t){(const char*)ptr, (size_t)len};
}

// The list header: its tag, which must be the one given or zero (absent), and the count
static uint64_t cdf_read_list(cdf_cursor_t* c, uint32_t tag) {
    uint32_t t = 0;
    md_xdr_read_u32(&c->xdr, &t);
    const uint64_t n = cdf_read_count(c);
    if (t != tag && !(t == 0 && n == 0)) {
        c->xdr.error = true;
        return 0;
    }
    return n;
}

static bool cdf_read_atts(cdf_cursor_t* c, cdf_att_t** out, size_t* out_count, md_allocator_i* alloc) {
    const uint64_t n = cdf_read_list(c, CDF_TAG_ATTRIBUTE);
    if (!md_xdr_ok(&c->xdr) || n > md_xdr_remaining(&c->xdr)) return false;
    cdf_att_t* att = n ? md_alloc(alloc, n * sizeof(cdf_att_t)) : NULL;
    for (uint64_t i = 0; i < n; ++i) {
        att[i].name = cdf_read_name(c);
        md_xdr_read_u32(&c->xdr, &att[i].type);
        att[i].count = cdf_read_count(c);
        const size_t elem = cdf_type_size(att[i].type);
        if (!md_xdr_ok(&c->xdr) || !elem || att[i].count > md_xdr_remaining(&c->xdr) / elem) return false;
        if (!md_xdr_read_opaque(&c->xdr, &att[i].data, (size_t)(att[i].count * elem))) return false;
    }
    *out = att;
    *out_count = (size_t)n;
    return md_xdr_ok(&c->xdr);
}

typedef enum cdf_parse_result_t {
    CDF_PARSE_OK,
    CDF_PARSE_NOT_CDF,
    CDF_PARSE_HDF5,
    CDF_PARSE_TRUNCATED,    // the buffer ended inside the header: read more of the file and try again
    CDF_PARSE_INVALID,
} cdf_parse_result_t;

static cdf_parse_result_t cdf_parse_buffer(cdf_t* cdf, const uint8_t* data, size_t size, md_allocator_i* alloc) {
    MEMSET(cdf, 0, sizeof(*cdf));
    if (size >= 8 && memcmp(data, "\211HDF\r\n\032\n", 8) == 0) return CDF_PARSE_HDF5;
    if (size < 4 || data[0] != 'C' || data[1] != 'D' || data[2] != 'F') return CDF_PARSE_NOT_CDF;
    const uint32_t version = data[3];
    if (version != 1 && version != 2 && version != 5) return CDF_PARSE_NOT_CDF;

    cdf_cursor_t c = { .xdr = md_xdr_init(data + 4, size - 4), .version = version };
    cdf->version = version;

    if (version == 5) {
        uint64_t v = 0;
        md_xdr_read_u64(&c.xdr, &v);
        cdf->streaming = (v == UINT64_MAX);
        cdf->numrecs = v;
    } else {
        uint32_t v = 0;
        md_xdr_read_u32(&c.xdr, &v);
        cdf->streaming = (v == UINT32_MAX);
        cdf->numrecs = v;
    }

    const uint64_t num_dim = cdf_read_list(&c, CDF_TAG_DIMENSION);
    if (!md_xdr_ok(&c.xdr)) goto fail;
    if (num_dim > md_xdr_remaining(&c.xdr)) goto fail;
    cdf->dim = num_dim ? md_alloc(alloc, num_dim * sizeof(cdf_dim_t)) : NULL;
    cdf->num_dim = (size_t)num_dim;
    int record_dim = -1;
    for (uint64_t i = 0; i < num_dim; ++i) {
        cdf->dim[i].name = cdf_read_name(&c);
        cdf->dim[i].len  = cdf_read_count(&c);
        if (cdf->dim[i].len == 0) {
            if (record_dim >= 0) goto invalid;  // one record dimension at most
            record_dim = (int)i;
        }
    }
    if (!md_xdr_ok(&c.xdr)) goto fail;

    if (!cdf_read_atts(&c, &cdf->gatt, &cdf->num_gatt, alloc)) goto fail;

    const uint64_t num_var = cdf_read_list(&c, CDF_TAG_VARIABLE);
    if (!md_xdr_ok(&c.xdr)) goto fail;
    if (num_var > md_xdr_remaining(&c.xdr)) goto fail;
    cdf->var = num_var ? md_alloc(alloc, num_var * sizeof(cdf_var_t)) : NULL;
    cdf->num_var = (size_t)num_var;
    size_t num_record_vars = 0;
    const cdf_var_t* first_record_var = NULL;
    for (uint64_t i = 0; i < num_var; ++i) {
        cdf_var_t* v = &cdf->var[i];
        MEMSET(v, 0, sizeof(*v));
        v->name = cdf_read_name(&c);
        const uint64_t rank = cdf_read_count(&c);
        if (!md_xdr_ok(&c.xdr)) goto fail;
        if (rank > CDF_MAX_RANK) goto invalid;
        v->rank = (uint32_t)rank;
        for (uint32_t d = 0; d < v->rank; ++d) {
            const uint64_t id = cdf_read_count(&c);
            if (md_xdr_ok(&c.xdr) && id >= num_dim) goto invalid;
            v->dim[d] = (uint32_t)id;
        }
        if (!cdf_read_atts(&c, &v->att, &v->num_att, alloc)) goto fail;
        md_xdr_read_u32(&c.xdr, &v->type);
        (void)cdf_read_count(&c);   // vsize: computed below, it saturates for variables past 4 GB
        v->begin = cdf_read_offset(&c);
        if (!md_xdr_ok(&c.xdr)) goto fail;

        const size_t elem = cdf_type_size(v->type);
        if (!elem) goto invalid;
        uint64_t bytes = elem;
        for (uint32_t d = 0; d < v->rank; ++d) {
            if ((int)v->dim[d] == record_dim) {
                if (d != 0) goto invalid;   // the record dimension can only be the first
                v->record = true;
                continue;
            }
            const uint64_t len = cdf->dim[v->dim[d]].len;
            if (len && bytes > UINT64_MAX / len) goto invalid;
            bytes *= len;
        }
        v->slab = bytes;
        if (v->record) {
            num_record_vars += 1;
            if (!first_record_var) first_record_var = v;
            cdf->recsize += (bytes + 3) & ~(uint64_t)3;
        }
    }
    // One record variable alone is not padded (NC_computeshapes in NetCDF-C)
    if (num_record_vars == 1) {
        cdf->recsize = first_record_var->slab;
    }
    return CDF_PARSE_OK;

fail:
    // A read past the end of the buffer: the header may go on in the file (cdf_open reads more and tries again)
    return CDF_PARSE_TRUNCATED;
invalid:
    return CDF_PARSE_INVALID;
}

// Reads and parses the header of the file. Its size is not known in advance: the first block holds it for any
// file of the Amber convention, and a larger header is read again with a larger block.
static cdf_parse_result_t cdf_open(cdf_t* cdf, md_file_t* out_file, str_t path, md_allocator_i* alloc, bool quiet) {
    md_file_t file = {0};
    if (!md_file_open(&file, path, MD_FILE_READ)) {
        if (!quiet) MD_LOG_ERROR("Amber NetCDF: cannot open '" STR_FMT "'", STR_ARG(path));
        return CDF_PARSE_NOT_CDF;
    }
    const uint64_t file_size = (uint64_t)md_file_size(file);
    size_t block = 64 * 1024;
    cdf_parse_result_t result = CDF_PARSE_INVALID;
    for (;;) {
        const size_t want = (size_t)MIN((uint64_t)block, file_size);
        uint8_t* buf = md_alloc(alloc, want + 1);
        const size_t got = md_file_read_at(file, 0, buf, want);
        result = cdf_parse_buffer(cdf, buf, got, alloc);
        if (result != CDF_PARSE_TRUNCATED || want >= file_size || block >= CDF_MAX_HEADER_BYTES) break;
        block *= 4;
    }
    if (result == CDF_PARSE_TRUNCATED) result = CDF_PARSE_INVALID;
    cdf->file_size = file_size;
    if (result == CDF_PARSE_OK && out_file) {
        *out_file = file;
    } else {
        md_file_close(&file);
    }
    return result;
}

static const cdf_att_t* cdf_find_att(const cdf_att_t* att, size_t count, const char* name) {
    for (size_t i = 0; i < count; ++i) {
        if (str_eq_cstr(att[i].name, name)) return &att[i];
    }
    return NULL;
}

static str_t cdf_att_text(const cdf_att_t* att) {
    if (!att || att->type != CDF_CHAR) return (str_t){0};
    str_t s = {(const char*)att->data, (size_t)att->count};
    // C writers include the terminator in the count
    while (s.len && s.ptr[s.len - 1] == '\0') s.len -= 1;
    return str_trim(s);
}

static bool cdf_att_number(const cdf_att_t* att, double* out) {
    if (!att || att->count < 1) return false;
    const uint8_t* p = att->data;
    switch (att->type) {
    case CDF_BYTE:   *out = (double)(int8_t)p[0]; return true;
    case CDF_UBYTE:  *out = (double)p[0]; return true;
    case CDF_SHORT:  *out = (double)(int16_t)((p[0] << 8) | p[1]); return true;
    case CDF_USHORT: *out = (double)(uint16_t)((p[0] << 8) | p[1]); return true;
    case CDF_INT:    *out = (double)md_xdr_load_i32(p); return true;
    case CDF_UINT:   *out = (double)md_xdr_load_u32(p); return true;
    case CDF_FLOAT:  *out = (double)md_xdr_load_f32(p); return true;
    case CDF_DOUBLE: *out = md_xdr_load_f64(p); return true;
    case CDF_INT64:  *out = (double)md_xdr_load_i64(p); return true;
    case CDF_UINT64: *out = (double)md_xdr_load_u64(p); return true;
    default: return false;
    }
}

static const cdf_var_t* cdf_find_var(const cdf_t* cdf, const char* name) {
    for (size_t i = 0; i < cdf->num_var; ++i) {
        if (str_eq_cstr(cdf->var[i].name, name)) return &cdf->var[i];
    }
    return NULL;
}

static int cdf_find_dim(const cdf_t* cdf, const char* name) {
    for (size_t i = 0; i < cdf->num_dim; ++i) {
        if (str_eq_cstr(cdf->dim[i].name, name)) return (int)i;
    }
    return -1;
}

// ## THE AMBER CONVENTION

// A token of Conventions, which may list several separated by commas or spaces
static bool amber_has_convention(str_t conventions, const char* name) {
    size_t i = 0;
    while (i < conventions.len) {
        while (i < conventions.len && (conventions.ptr[i] == ',' || conventions.ptr[i] == ' ')) ++i;
        size_t j = i;
        while (j < conventions.len && conventions.ptr[j] != ',' && conventions.ptr[j] != ' ') ++j;
        if (j > i && str_eq_cstr((str_t){conventions.ptr + i, j - i}, name)) return true;
        i = j;
    }
    return false;
}

static md_amber_netcdf_kind_t amber_kind(const cdf_t* cdf) {
    const str_t conventions = cdf_att_text(cdf_find_att(cdf->gatt, cdf->num_gatt, "Conventions"));
    if (amber_has_convention(conventions, "AMBERRESTART")) return MD_AMBER_NETCDF_RESTART;
    if (amber_has_convention(conventions, "AMBER"))        return MD_AMBER_NETCDF_TRAJECTORY;
    return MD_AMBER_NETCDF_NONE;
}

// How one variable is laid out frame by frame: frame f is at begin + f * stride, bytes long
typedef struct amber_section_t {
    bool     present;
    uint32_t type;          // CDF_FLOAT, CDF_DOUBLE, ...
    uint64_t begin;
    uint64_t stride;        // 0 in a restart, which has one frame
    uint64_t bytes;         // of one frame
    double   scale;         // the variable's scale_factor, 1 without one
    str_t    units;
} amber_section_t;

// What the file is to the reader, from its header
typedef struct amber_layout_t {
    md_amber_netcdf_kind_t kind;
    size_t          num_atoms;
    size_t          num_frames;     // the frames complete in the file
    amber_section_t coordinates;
    amber_section_t velocities;
    amber_section_t forces;
    amber_section_t time;
    amber_section_t cell_lengths;
    amber_section_t cell_angles;
    amber_section_t temp0;
    amber_section_t repidx;
    amber_section_t crdidx;
} amber_layout_t;

// The section of variable name, whose dimensions after the frame (in a trajectory) must be inner[]. False with
// the section absent when the variable is not there; false and logged when it is there but not as the
// convention has it.
static bool amber_section(amber_section_t* out, const cdf_t* cdf, const char* name, bool restart, int frame_dim, const int inner[], size_t num_inner, str_t path) {
    MEMSET(out, 0, sizeof(*out));
    const cdf_var_t* v = cdf_find_var(cdf, name);
    if (!v) return true;

    const uint32_t lead = restart ? 0 : 1;
    bool ok = v->rank == lead + num_inner;
    if (ok && !restart) ok = (int)v->dim[0] == frame_dim;
    for (size_t d = 0; ok && d < num_inner; ++d) {
        ok = (int)v->dim[lead + d] == inner[d];
    }
    switch (v->type) {
    case CDF_FLOAT: case CDF_DOUBLE: case CDF_SHORT: case CDF_INT: break;
    default: ok = false;
    }
    if (!ok) {
        MD_LOG_ERROR("Amber NetCDF: the variable '%s' in '" STR_FMT "' is not laid out as the Amber convention has it", name, STR_ARG(path));
        return false;
    }

    out->present = true;
    out->type  = v->type;
    out->begin = v->begin;
    out->bytes = v->slab;
    if (!restart) {
        if (v->record) {
            out->stride = cdf->recsize;
        } else {
            // A fixed frame dimension: the frames are consecutive slabs of the one array
            const uint64_t frames = cdf->dim[frame_dim].len;
            out->bytes  = frames ? v->slab / frames : 0;
            out->stride = out->bytes;
        }
    }
    out->scale = 1.0;
    double scale;
    if (cdf_att_number(cdf_find_att(v->att, v->num_att, "scale_factor"), &scale) && isfinite(scale) && scale != 0.0) {
        out->scale = scale;
    }
    out->units = cdf_att_text(cdf_find_att(v->att, v->num_att, "units"));
    return true;
}

// Frames of the section complete in a file of the given size
static uint64_t amber_section_frames(const amber_section_t* s, uint64_t file_size) {
    if (!s->present) return UINT64_MAX;
    if (s->begin + s->bytes > file_size) return 0;
    if (s->stride == 0) return 1;
    return (file_size - s->begin - s->bytes) / s->stride + 1;
}

static bool amber_layout(amber_layout_t* out, const cdf_t* cdf, str_t path) {
    MEMSET(out, 0, sizeof(*out));
    out->kind = amber_kind(cdf);
    if (out->kind == MD_AMBER_NETCDF_NONE) {
        MD_LOG_ERROR("Amber NetCDF: '" STR_FMT "' is NetCDF but does not follow the Amber convention (no Conventions \"AMBER\" or \"AMBERRESTART\")", STR_ARG(path));
        return false;
    }
    const bool restart = out->kind == MD_AMBER_NETCDF_RESTART;

    const str_t version = cdf_att_text(cdf_find_att(cdf->gatt, cdf->num_gatt, "ConventionVersion"));
    if (!str_empty(version) && !str_eq_cstr(version, "1.0")) {
        MD_LOG_INFO("Amber NetCDF: '" STR_FMT "' declares convention version '" STR_FMT "', read as 1.0", STR_ARG(path), STR_ARG(version));
    }

    const int atom    = cdf_find_dim(cdf, "atom");
    const int spatial = cdf_find_dim(cdf, "spatial");
    const int frame   = restart ? -1 : cdf_find_dim(cdf, "frame");
    if (atom < 0 || spatial < 0 || (!restart && frame < 0) || cdf->dim[spatial].len != 3 || cdf->dim[atom].len == 0) {
        MD_LOG_ERROR("Amber NetCDF: '" STR_FMT "' lacks the atom, spatial (3) or frame dimensions of the convention", STR_ARG(path));
        return false;
    }
    out->num_atoms = (size_t)cdf->dim[atom].len;

    const int xyz[2] = { atom, spatial };
    const int cell_spatial = cdf_find_dim(cdf, "cell_spatial");
    const int cell_angular = cdf_find_dim(cdf, "cell_angular");

    bool ok = amber_section(&out->coordinates, cdf, "coordinates", restart, frame, xyz, 2, path)
           && amber_section(&out->velocities,  cdf, "velocities",  restart, frame, xyz, 2, path)
           && amber_section(&out->forces,      cdf, "forces",      restart, frame, xyz, 2, path)
           && amber_section(&out->time,        cdf, "time",        restart, frame, NULL, 0, path);
    if (ok && cell_spatial >= 0) ok = amber_section(&out->cell_lengths, cdf, "cell_lengths", restart, frame, &cell_spatial, 1, path);
    if (ok && cell_angular >= 0) ok = amber_section(&out->cell_angles,  cdf, "cell_angles",  restart, frame, &cell_angular, 1, path);
    if (ok && !restart) {
        // Replica exchange, per frame. Anything else about it is left out.
        ok = amber_section(&out->temp0,  cdf, "temp0",       restart, frame, NULL, 0, path)
          && amber_section(&out->repidx, cdf, "remd_repidx", restart, frame, NULL, 0, path)
          && amber_section(&out->crdidx, cdf, "remd_crdidx", restart, frame, NULL, 0, path);
    }
    if (!ok) return false;

    if (!out->coordinates.present) {
        MD_LOG_ERROR("Amber NetCDF: '" STR_FMT "' has no coordinates", STR_ARG(path));
        return false;
    }
    if (out->cell_lengths.present && (cell_spatial < 0 || cdf->dim[cell_spatial].len != 3)) out->cell_lengths.present = false;
    if (out->cell_angles.present  && (cell_angular < 0 || cdf->dim[cell_angular].len != 3)) out->cell_angles.present  = false;

    if (!str_empty(out->coordinates.units) && !str_eq_cstr_ignore_case(out->coordinates.units, "angstrom") && !str_eq_cstr_ignore_case(out->coordinates.units, "angstroms")) {
        MD_LOG_ERROR("Amber NetCDF: the coordinates of '" STR_FMT "' are in '" STR_FMT "', the convention has angstrom", STR_ARG(path), STR_ARG(out->coordinates.units));
        return false;
    }

    // ## Frames: as many as the header says, of those as many as the file holds completely. A trajectory still
    // being written has its last record coming, and one cut short lacks some.
    uint64_t frames = 1;
    if (!restart) {
        frames = cdf->dim[frame].len;
        if (frames == 0) {
            frames = cdf->streaming ? UINT64_MAX : cdf->numrecs;
        }
    }
    const amber_section_t* sections[] = { &out->coordinates, &out->velocities, &out->forces, &out->time, &out->cell_lengths, &out->cell_angles, &out->temp0, &out->repidx, &out->crdidx };
    uint64_t complete = UINT64_MAX;
    for (size_t i = 0; i < ARRAY_SIZE(sections); ++i) {
        complete = MIN(complete, amber_section_frames(sections[i], cdf->file_size));
    }
    if (complete < frames) {
        if (frames != UINT64_MAX) {
            MD_LOG_INFO("Amber NetCDF: '" STR_FMT "' holds %llu complete frames of the %llu its header says", STR_ARG(path),
                (unsigned long long)complete, (unsigned long long)frames);
        }
        frames = complete;
    }
    if (frames == 0 || frames > UINT32_MAX) {
        MD_LOG_ERROR("Amber NetCDF: '" STR_FMT "' holds no complete frame", STR_ARG(path));
        return false;
    }
    out->num_frames = (size_t)frames;
    return true;
}

static md_amber_netcdf_kind_t amber_open(amber_layout_t* layout, md_file_t* out_file, str_t path, md_allocator_i* alloc, bool quiet) {
    cdf_t cdf;
    const cdf_parse_result_t r = cdf_open(&cdf, out_file, path, alloc, quiet);
    if (r == CDF_PARSE_HDF5) {
        if (!quiet) {
            MD_LOG_ERROR("Amber NetCDF: '" STR_FMT "' is NetCDF-4 (HDF5), which the Amber convention does not allow and which is not supported. "
                "Convert it with cpptraj or MDAnalysis to the classic format", STR_ARG(path));
        }
        return MD_AMBER_NETCDF_HDF5;
    }
    if (r != CDF_PARSE_OK) {
        if (!quiet && r != CDF_PARSE_NOT_CDF) MD_LOG_ERROR("Amber NetCDF: '" STR_FMT "' is not a NetCDF file, or its header is damaged", STR_ARG(path));
        else if (!quiet && md_path_is_valid(path)) MD_LOG_ERROR("Amber NetCDF: '" STR_FMT "' is not a NetCDF file", STR_ARG(path));
        return MD_AMBER_NETCDF_NONE;
    }
    if (quiet) {
        // The kind alone, without validating the rest
        MEMSET(layout, 0, sizeof(*layout));
        layout->kind = amber_kind(&cdf);
        const int atom = cdf_find_dim(&cdf, "atom");
        layout->num_atoms = atom >= 0 ? (size_t)cdf.dim[atom].len : 0;
        if (out_file) md_file_close(out_file);
        return layout->kind;
    }
    if (!amber_layout(layout, &cdf, path)) {
        if (out_file) md_file_close(out_file);
        return MD_AMBER_NETCDF_NONE;
    }
    return layout->kind;
}

md_amber_netcdf_kind_t md_amber_netcdf_file_kind(str_t filename) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(256));
    amber_layout_t layout;
    const md_amber_netcdf_kind_t kind = amber_open(&layout, NULL, filename, arena, true);
    md_arena_allocator_destroy(arena);
    return kind;
}

size_t md_amber_netcdf_file_atom_count(str_t filename) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(256));
    amber_layout_t layout;
    const md_amber_netcdf_kind_t kind = amber_open(&layout, NULL, filename, arena, true);
    md_arena_allocator_destroy(arena);
    return (kind == MD_AMBER_NETCDF_TRAJECTORY || kind == MD_AMBER_NETCDF_RESTART) ? layout.num_atoms : 0;
}

// ## READING VALUES

// count values of the section's type, converted to double and scaled
static void amber_convert_f64(double* out, const uint8_t* raw, uint32_t type, size_t count, double scale) {
    switch (type) {
    case CDF_FLOAT:  for (size_t i = 0; i < count; ++i) out[i] = (double)md_xdr_load_f32(raw + i * 4); break;
    case CDF_DOUBLE: for (size_t i = 0; i < count; ++i) out[i] = md_xdr_load_f64(raw + i * 8); break;
    case CDF_INT:    for (size_t i = 0; i < count; ++i) out[i] = (double)md_xdr_load_i32(raw + i * 4); break;
    case CDF_SHORT:  for (size_t i = 0; i < count; ++i) out[i] = (double)(int16_t)((raw[i * 2] << 8) | raw[i * 2 + 1]); break;
    default:         for (size_t i = 0; i < count; ++i) out[i] = 0.0; break;
    }
    if (scale != 1.0) {
        for (size_t i = 0; i < count; ++i) out[i] *= scale;
    }
}

static void amber_convert_f32(float* out, const uint8_t* raw, uint32_t type, size_t count, double scale) {
    const float s = (float)scale;
    switch (type) {
    case CDF_FLOAT:
        if (scale == 1.0) {
            md_xdr_load_f32_array(out, raw, count);
        } else {
            for (size_t i = 0; i < count; ++i) out[i] = md_xdr_load_f32(raw + i * 4) * s;
        }
        break;
    case CDF_DOUBLE: for (size_t i = 0; i < count; ++i) out[i] = (float)(md_xdr_load_f64(raw + i * 8) * scale); break;
    case CDF_INT:    for (size_t i = 0; i < count; ++i) out[i] = (float)((double)md_xdr_load_i32(raw + i * 4) * scale); break;
    case CDF_SHORT:  for (size_t i = 0; i < count; ++i) out[i] = (float)((double)(int16_t)((raw[i * 2] << 8) | raw[i * 2 + 1]) * scale); break;
    default:         MEMSET(out, 0, count * sizeof(float)); break;
    }
}

// The values of a section of one value per frame (or a few), frames [0, F), as doubles. One read for a
// variable along a fixed frame dimension, one per frame for a record variable.
static bool amber_read_series(double* out, size_t values_per_frame, md_file_t file, const amber_section_t* s, size_t F, md_allocator_i* alloc) {
    const size_t elem = cdf_type_size(s->type);
    if (s->bytes != values_per_frame * elem) return false;
    if (s->stride == s->bytes) {
        const size_t bytes = F * (size_t)s->bytes;
        uint8_t* raw = md_alloc(alloc, bytes);
        if (md_file_read_at(file, (int64_t)s->begin, raw, bytes) != bytes) return false;
        amber_convert_f64(out, raw, s->type, F * values_per_frame, s->scale);
        md_free(alloc, raw, bytes);
        return true;
    }
    uint8_t raw[64];
    if (s->bytes > sizeof(raw)) return false;
    for (size_t f = 0; f < F; ++f) {
        if (md_file_read_at(file, (int64_t)(s->begin + f * s->stride), raw, (size_t)s->bytes) != s->bytes) return false;
        amber_convert_f64(out + f * values_per_frame, raw, s->type, values_per_frame, s->scale);
    }
    return true;
}

// The cell from lengths and angles. Zero lengths are no box (some writers keep the variables for a run
// without one); angles default to 90 when the file has lengths alone.
static md_unitcell_t amber_cell(const double lengths[3], const double* angles) {
    if (!(lengths[0] > 0.0 && lengths[1] > 0.0 && lengths[2] > 0.0) || !isfinite(lengths[0] + lengths[1] + lengths[2])) {
        return md_unitcell_none();
    }
    double a[3] = { 90.0, 90.0, 90.0 };
    if (angles && angles[0] > 0.0 && angles[0] < 180.0 && angles[1] > 0.0 && angles[1] < 180.0 && angles[2] > 0.0 && angles[2] < 180.0) {
        a[0] = angles[0]; a[1] = angles[1]; a[2] = angles[2];
    }
    return md_unitcell_from_extent_and_angles(lengths[0], lengths[1], lengths[2], a[0], a[1], a[2]);
}

// ## INITIAL STATE

bool md_amber_netcdf_state_init_from_file(md_system_state_t* state, str_t filename, size_t num_atoms) {
    ASSERT(state);
    if (!state->alloc) {
        MD_LOG_ERROR("Amber NetCDF: the state has no allocator");
        return false;
    }
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    bool result = false;
    md_file_t file = {0};
    amber_layout_t layout;
    const md_amber_netcdf_kind_t kind = amber_open(&layout, &file, filename, arena, false);
    if (kind != MD_AMBER_NETCDF_TRAJECTORY && kind != MD_AMBER_NETCDF_RESTART) goto done;

    const size_t N = layout.num_atoms;
    if (num_atoms && N != num_atoms) {
        MD_LOG_ERROR("Amber NetCDF: '" STR_FMT "' holds %zu atoms, the system %zu", STR_ARG(filename), N, num_atoms);
        goto done;
    }

    const amber_section_t* xyz = &layout.coordinates;
    const size_t bytes = (size_t)xyz->bytes;
    uint8_t* raw = md_alloc(arena, bytes);
    if (bytes != N * 3 * cdf_type_size(xyz->type) || md_file_read_at(file, (int64_t)xyz->begin, raw, bytes) != bytes) {
        MD_LOG_ERROR("Amber NetCDF: failed to read the coordinates of '" STR_FMT "'", STR_ARG(filename));
        goto done;
    }

    md_unitcell_t cell = md_unitcell_none();
    if (layout.cell_lengths.present) {
        double lengths[3], angles[3];
        const bool has_angles = layout.cell_angles.present && amber_read_series(angles, 3, file, &layout.cell_angles, 1, arena);
        if (amber_read_series(lengths, 3, file, &layout.cell_lengths, 1, arena)) {
            cell = amber_cell(lengths, has_angles ? angles : NULL);
        }
    }

    if (!md_system_state_init(state, N)) goto done;
    amber_convert_f32((float*)state->xyz, raw, xyz->type, N * 3, xyz->scale);
    state->num_atoms = N;
    state->unitcell = cell;
    result = true;

done:
    if (md_file_valid(file)) md_file_close(&file);
    md_arena_allocator_destroy(arena);
    return result;
}

// ## RUN

// The per atom sections a run publishes
static const struct {
    const char* leaf;
} amber_run_sections[] = {
    { "atom/position" },
    { "atom/velocity" },
    { "atom/force" },
};

// <run>/source/layout: per section begin, stride, type and 0; <run>/source/scale: per section its scale_factor
enum { AMBER_LAYOUT_STRIDE = 4 };

// <run>/atom/{position,velocity,force}: frame f of the section is at begin + f * stride
static size_t amber_section_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    const md_system_t* sys = (const md_system_t*)user_data;
    ASSERT(sys);
    const md_attributes_t* attributes = &sys->attributes;
    if (!slice || slice->num_idx == 0 || slice->num_idx > 2) return 0;

    int kind = -1;
    md_run_source_t src = {0};
    for (int i = 0; i < (int)ARRAY_SIZE(amber_run_sections); ++i) {
        if (md_run_source(&src, attributes, attr, str_from_cstr(amber_run_sections[i].leaf))) {
            kind = i;
            break;
        }
    }
    if (kind < 0) return 0;

    const md_attribute_t* layout_attr = md_attributes_find_in(attributes, src.run, STR_LIT("source/layout"));
    const md_attribute_t* scale_attr  = md_attributes_find_in(attributes, src.run, STR_LIT("source/scale"));
    if (!layout_attr || !scale_attr || layout_attr->format.type != MD_ATTRIBUTE_TYPE_I64 || scale_attr->format.type != MD_ATTRIBUTE_TYPE_F64 ||
        md_attribute_value_count(&layout_attr->format) != ARRAY_SIZE(amber_run_sections) * AMBER_LAYOUT_STRIDE ||
        md_attribute_value_count(&scale_attr->format) != ARRAY_SIZE(amber_run_sections) || !layout_attr->data || !scale_attr->data) {
        MD_LOG_ERROR("Amber NetCDF: the run '" STR_FMT "' has lost its source attributes", STR_ARG(src.run));
        return 0;
    }
    const int64_t* layout = (const int64_t*)layout_attr->data + kind * AMBER_LAYOUT_STRIDE;
    const double scale = ((const double*)scale_attr->data)[kind];
    const int64_t begin  = layout[0];
    const int64_t stride = layout[1];
    const uint32_t type  = (uint32_t)layout[2];
    const size_t elem = cdf_type_size(type);
    if (!elem) return 0;

    const uint32_t frame = slice->idx[0];
    const size_t N = attr->format.shape[1];
    if (frame >= attr->format.shape[0]) return 0;
    size_t first = 0, count = N;
    if (slice->num_idx == 2) {
        if (slice->idx[1] >= N) return 0;
        first = slice->idx[1];
        count = 1;
    }
    if (cap != count * 3) return 0;

    const int64_t offset = begin + (int64_t)frame * stride + (int64_t)(first * 3 * elem);
    const size_t bytes = count * 3 * elem;
    md_temp_scope_t temp = md_temp_begin();
    size_t written = 0;
    uint8_t* raw = md_temp_alloc(temp, bytes);
    if (raw && md_attribute_io_read_at(io, src.path, offset, raw, bytes) == bytes) {
        amber_convert_f32((float*)dst, raw, type, count * 3, scale);
        written = cap;
    } else {
        MD_LOG_ERROR("Amber NetCDF: failed to read frame %u from '" STR_FMT "'", frame, STR_ARG(src.path));
    }
    md_temp_end(temp);
    return written;
}

static md_unit_t amber_time_unit(str_t units, bool* known) {
    *known = true;
    if (str_empty(units) || str_eq_cstr_ignore_case(units, "picosecond") || str_eq_cstr_ignore_case(units, "picoseconds") || str_eq_cstr_ignore_case(units, "ps")) {
        return md_unit_picosecond();
    }
    if (str_eq_cstr_ignore_case(units, "nanosecond") || str_eq_cstr_ignore_case(units, "nanoseconds") || str_eq_cstr_ignore_case(units, "ns")) {
        return md_unit_nanosecond();
    }
    if (str_eq_cstr_ignore_case(units, "femtosecond") || str_eq_cstr_ignore_case(units, "femtoseconds") || str_eq_cstr_ignore_case(units, "fs")) {
        return md_unit_femtosecond();
    }
    *known = false;
    return md_unit_none();
}

bool md_amber_netcdf_system_publish_run(md_system_t* sys, str_t filename, str_t run, uint32_t flags) {
    ASSERT(sys);
    (void)flags;    // Nothing to index: the header says where every frame is
    char path_buf[4096];
    const size_t path_len = md_path_write_canonical(path_buf, sizeof(path_buf), filename);
    const str_t path = {path_buf, path_len};

    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    bool result = false;
    md_file_t file = {0};
    amber_layout_t layout;
    if (path_len == 0) {
        MD_LOG_ERROR("Amber NetCDF: cannot resolve '" STR_FMT "'", STR_ARG(filename));
        goto done;
    }
    const md_amber_netcdf_kind_t kind = amber_open(&layout, &file, path, arena, false);
    if (kind == MD_AMBER_NETCDF_RESTART) {
        MD_LOG_ERROR("Amber NetCDF: '" STR_FMT "' is a restart, a single set of coordinates rather than a trajectory", STR_ARG(path));
        goto done;
    }
    if (kind != MD_AMBER_NETCDF_TRAJECTORY) goto done;

    const size_t F = layout.num_frames;
    const size_t N = layout.num_atoms;
    if (layout.coordinates.bytes != N * 3 * cdf_type_size(layout.coordinates.type)) {
        MD_LOG_ERROR("Amber NetCDF: the coordinates of '" STR_FMT "' do not hold 3 values per atom", STR_ARG(path));
        goto done;
    }

    // ## Time. Without it, or with one that does not increase (concatenated runs which each start over), the
    // frames are ordinals: the frame axis has to be non decreasing to be matched against.
    double* times = md_alloc(arena, F * sizeof(double));
    md_unit_t time_unit = md_unit_none();
    bool use_time = false;
    if (layout.time.present) {
        bool known = false;
        time_unit = amber_time_unit(layout.time.units, &known);
        if (!known) {
            MD_LOG_INFO("Amber NetCDF: the time of '" STR_FMT "' is in '" STR_FMT "', which is not known here; the frames are numbered instead", STR_ARG(path), STR_ARG(layout.time.units));
        } else if (amber_read_series(times, 1, file, &layout.time, F, arena)) {
            use_time = true;
            for (size_t f = 0; f < F; ++f) {
                if (!isfinite(times[f]) || (f > 0 && times[f] <= times[f - 1])) {
                    MD_LOG_INFO("Amber NetCDF: the time of '" STR_FMT "' does not increase from frame to frame (at frame %zu, a file written without time has zeros); the frames are numbered instead", STR_ARG(path), f);
                    use_time = false;
                    break;
                }
            }
        } else {
            MD_LOG_ERROR("Amber NetCDF: failed to read the time of '" STR_FMT "'", STR_ARG(path));
            goto done;
        }
    }
    if (!use_time) {
        for (size_t f = 0; f < F; ++f) times[f] = (double)f;
        time_unit = md_unit_none();
    }

    // ## Cells
    float* boxes = md_alloc(arena, F * 9 * sizeof(float));
    MEMSET(boxes, 0, F * 9 * sizeof(float));
    if (layout.cell_lengths.present) {
        double* lengths = md_alloc(arena, F * 3 * sizeof(double));
        double* angles  = layout.cell_angles.present ? md_alloc(arena, F * 3 * sizeof(double)) : NULL;
        if (!amber_read_series(lengths, 3, file, &layout.cell_lengths, F, arena) ||
            (angles && !amber_read_series(angles, 3, file, &layout.cell_angles, F, arena))) {
            MD_LOG_ERROR("Amber NetCDF: failed to read the cells of '" STR_FMT "'", STR_ARG(path));
            goto done;
        }
        for (size_t f = 0; f < F; ++f) {
            const md_unitcell_t cell = amber_cell(lengths + f * 3, angles ? angles + f * 3 : NULL);
            md_unitcell_A_extract_float((float(*)[3])(boxes + f * 9), &cell);
        }
    }

    // ## Frame table: where each frame's coordinates are
    int64_t* offsets = md_alloc(arena, F * sizeof(int64_t));
    int64_t* sizes   = md_alloc(arena, F * sizeof(int64_t));
    for (size_t f = 0; f < F; ++f) {
        offsets[f] = (int64_t)(layout.coordinates.begin + f * layout.coordinates.stride);
        sizes[f]   = (int64_t)layout.coordinates.bytes;
    }

    const md_attribute_virtual_t virt = { .provider = amber_section_provider, .user_data = sys };
    const md_run_desc_t desc = {
        .num_frames    = F,
        .num_atoms     = N,
        .time          = times,
        .time_unit     = time_unit,
        .unitcell      = boxes,     // zero where the file has no box, as every run has it
        .source_path   = path,
        .source_offset = offsets,
        .source_size   = sizes,
        .position_virt = &virt,
    };
    if (!md_run_publish(sys, run, &desc)) goto done;

    md_attributes_t* attributes = &sys->attributes;
    char buf[512];
    const amber_section_t* sections[ARRAY_SIZE(amber_run_sections)] = { &layout.coordinates, &layout.velocities, &layout.forces };
    int64_t layout_values[ARRAY_SIZE(amber_run_sections) * AMBER_LAYOUT_STRIDE] = {0};
    double scales[ARRAY_SIZE(amber_run_sections)] = {0};
    for (size_t k = 0; k < ARRAY_SIZE(amber_run_sections); ++k) {
        const amber_section_t* s = sections[k];
        if (!s->present) continue;
        layout_values[k * AMBER_LAYOUT_STRIDE + 0] = (int64_t)s->begin;
        layout_values[k * AMBER_LAYOUT_STRIDE + 1] = (int64_t)s->stride;
        layout_values[k * AMBER_LAYOUT_STRIDE + 2] = (int64_t)s->type;
        scales[k] = s->scale;
    }

    bool ok = md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), run, STR_LIT("source/layout")),
        .format = { .type = MD_ATTRIBUTE_TYPE_I64, .components = 1, .rank = 1, .shape = { (uint32_t)ARRAY_SIZE(layout_values) } },
        .unit = md_unit_none(),
        .description = STR_INIT("Per section (position, velocity, force): byte offset of frame 0, bytes from frame to frame, NetCDF type, 0"),
        .data = layout_values, .byte_size = sizeof(layout_values)});
    ok = ok && md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), run, STR_LIT("source/scale")),
        .format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)ARRAY_SIZE(scales) } },
        .unit = md_unit_none(),
        .description = STR_INIT("Per section, the scale_factor its values are multiplied by"),
        .data = scales, .byte_size = sizeof(scales)});

    // ## Velocities, forces
    const md_unit_t kcal_per_mol = md_unit_div(md_unit_scl(md_unit_calorie(), 1.0e3), md_unit_mole());
    const struct { int k; md_unit_t unit; str_t label; } per_atom[] = {
        { 1, md_unit_div(md_unit_angstrom(), md_unit_picosecond()), STR_INIT("Velocity") },
        { 2, md_unit_div(kcal_per_mol, md_unit_angstrom()),         STR_INIT("Force") },
    };
    for (size_t i = 0; ok && i < ARRAY_SIZE(per_atom); ++i) {
        const int k = per_atom[i].k;
        if (!sections[k]->present) continue;
        if (sections[k]->bytes != N * 3 * cdf_type_size(sections[k]->type)) continue;
        ok = md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), run, str_from_cstr(amber_run_sections[k].leaf)),
            .format = { .type = MD_ATTRIBUTE_TYPE_F32, .components = 3, .rank = 2, .shape = { (uint32_t)F, (uint32_t)N } },
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = per_atom[i].unit, .label = per_atom[i].label,
            .virt = &virt});
    }

    // ## Replica exchange
    if (ok && layout.temp0.present) {
        double* temp0 = md_alloc(arena, F * sizeof(double));
        if (amber_read_series(temp0, 1, file, &layout.temp0, F, arena)) {
            ok = md_attributes_replace(attributes, &(md_attribute_desc_t){
                .path = md_run_path(buf, sizeof(buf), run, STR_LIT("remd/temp0")),
                .format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)F } },
                .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = md_unit_kelvin(), .label = STR_INIT("Target temperature"),
                .data = temp0, .byte_size = F * sizeof(double)});
        }
    }
    const struct { const amber_section_t* s; const char* leaf; str_t label; } indices[] = {
        { &layout.repidx, "remd/repidx", STR_INIT("Replica index") },
        { &layout.crdidx, "remd/crdidx", STR_INIT("Coordinate index") },
    };
    for (size_t i = 0; ok && i < ARRAY_SIZE(indices); ++i) {
        if (!indices[i].s->present) continue;
        double*  values = md_alloc(arena, F * sizeof(double));
        int32_t* ints   = md_alloc(arena, F * sizeof(int32_t));
        if (!amber_read_series(values, 1, file, indices[i].s, F, arena)) continue;
        for (size_t f = 0; f < F; ++f) ints[f] = (int32_t)values[f];
        ok = md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), run, str_from_cstr(indices[i].leaf)),
            .format = { .type = MD_ATTRIBUTE_TYPE_I32, .components = 1, .rank = 1, .shape = { (uint32_t)F } },
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = md_unit_none(), .label = indices[i].label,
            .data = ints, .byte_size = F * sizeof(int32_t)});
    }

    if (!ok) {
        MD_LOG_ERROR("Amber NetCDF: failed to publish '" STR_FMT "' as '" STR_FMT "'", STR_ARG(path), STR_ARG(run));
        md_attributes_remove_prefix(attributes, run);
        goto done;
    }
    result = true;

done:
    if (md_file_valid(file)) md_file_close(&file);
    md_arena_allocator_destroy(arena);
    return result;
}
