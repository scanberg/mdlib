#include <md_trr.h>
#include <md_xdr.h>
#include <md_system.h>

#include <md_util.h>

#include <core/md_common.h>
#include <core/md_array.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_vec_math.h>
#include <core/md_str_builder.h>

#include <string.h>
#include <stdio.h>


#define MD_TRR_CACHE_MAGIC 0x67b7cbab123452
#define MD_TRR_CACHE_VERSION 5

#define TRR_MAGIC 1993

#define TRR_MIN_HEADER_SIZE 72

#define DIM 3

// This file cherry picks bits and pieces from the provided xdrfile_trr.c implementation.
// See xdrfile_trr.c for the copyright specific to that file.
// The implementation is modified from its original to use mdlibs allocator and error logging for clarity.

// Taken from xdrfile_trr.c
typedef struct trr_header_t {
    bool use_double; /* Double precision?                    */
    int ir_size;     /* Backward compatibility               */
    int e_size;      /* Backward compatibility               */
    int box_size;    /* Non zero if a box is present         */
    int vir_size;    /* Backward compatibility               */
    int pres_size;   /* Backward compatibility               */
    int top_size;    /* Backward compatibility               */
    int sym_size;    /* Backward compatibility               */
    int x_size;      /* Non zero if coordinates are present  */
    int v_size;      /* Non zero if velocities are present   */
    int f_size;      /* Non zero if forces are present       */

    int natoms;     /* The total number of atoms            */
    int step;       /* Current step number                  */
    int nre;        /* Backward compatibility               */
    //float tf;       /* Current time                         */
    //float lambdaf;  /* Current value of lambda              */
    //double td;      /* Current time                         */
    //double lambdad; /* Current value of lambda              */
    // We only store the double version and if a float version is required, we cast it
    double t;
    double lambda;
} trr_header_t;

static bool n_float_size(const trr_header_t* sh, int* nflsz) {
    int nflsize = 0;

    if (sh->box_size) {
        nflsize = sh->box_size / (DIM * DIM);
    } else if (sh->x_size) {
        nflsize = sh->x_size / (sh->natoms * DIM);
    } else if (sh->v_size) {
        nflsize = sh->v_size / (sh->natoms * DIM);
    } else if (sh->f_size) {
        nflsize = sh->f_size / (sh->natoms * DIM);
    } else {
        MD_LOG_ERROR("TRR: Unexpected size in header");
        return false;
    }

    if (((nflsize != sizeof(float)) && (nflsize != sizeof(double)))) {
        MD_LOG_ERROR("TRR: Unexpected float size in header");
        return false;
    }

    *nflsz = nflsize;

    return true;
}

static int calc_framebytes(const trr_header_t* sh) {
    return sh->ir_size + sh->e_size + sh->box_size + sh->vir_size + sh->pres_size + sh->top_size +
        sh->sym_size + sh->x_size + sh->v_size + sh->f_size;
}

// Largest possible frame header: magic, the version string (length, then XDR string of 12
// characters), 13 ints, and time and lambda as doubles
#define TRR_MAX_HEADER_SIZE (4 + 4 + 4 + 12 + 13 * 4 + 2 * 8)

static bool trr_read_frame_header_buf(md_xdr_t* xdr, trr_header_t* sh) {
    ASSERT(xdr);
    ASSERT(sh);

    const char version[] = "GMX_trn_file";
    int32_t magic, slen, nflsz;

    if (!md_xdr_read_i32(xdr, &magic)) {
        MD_LOG_ERROR("TRR: Failed to read header magic number");
        return false;
    }
    if (magic != TRR_MAGIC) {
        MD_LOG_ERROR("TRR: Magic number did not match");
        return false;
    }

    // The length including the terminating zero, then the string as XDR (without it)
    str_t str;
    if (!md_xdr_read_i32(xdr, &slen) || !md_xdr_read_string(xdr, &str, 128)) {
        MD_LOG_ERROR("TRR: Failed to read header version string");
        return false;
    }
    if (slen != (int32_t)sizeof(version)) {
        MD_LOG_ERROR("TRR: Incorrect version string length");
        return false;
    }

    int32_t fields[13];
    if (!md_xdr_read_i32_array(xdr, fields, 13)) {
        MD_LOG_ERROR("TRR: Failed to read header fields");
        return false;
    }
    sh->ir_size   = fields[0];
    sh->e_size    = fields[1];
    sh->box_size  = fields[2];
    sh->vir_size  = fields[3];
    sh->pres_size = fields[4];
    sh->top_size  = fields[5];
    sh->sym_size  = fields[6];
    sh->x_size    = fields[7];
    sh->v_size    = fields[8];
    sh->f_size    = fields[9];
    sh->natoms    = fields[10];
    sh->step      = fields[11];
    sh->nre       = fields[12];

    if (!n_float_size(sh, &nflsz)) {
        return false;
    }
    sh->use_double = (nflsz == sizeof(double));

    if (sh->use_double) {
        md_xdr_read_f64(xdr, &sh->t);
        md_xdr_read_f64(xdr, &sh->lambda);
    } else {
        float tf, lf;
        md_xdr_read_f32(xdr, &tf);
        md_xdr_read_f32(xdr, &lf);
        sh->t      = tf;
        sh->lambda = lf;
    }
    if (!md_xdr_ok(xdr)) {
        MD_LOG_ERROR("TRR: Failed to read header time and lambda");
        return false;
    }
    return true;
}

// Which per frame sections a frame carries, as the cache stores it
enum {
    TRR_SECTION_X   = 1,
    TRR_SECTION_V   = 2,
    TRR_SECTION_F   = 4,
    TRR_SECTION_BOX = 8,
};

static uint32_t trr_sections(const trr_header_t* sh) {
    return (sh->x_size   ? TRR_SECTION_X   : 0) |
           (sh->v_size   ? TRR_SECTION_V   : 0) |
           (sh->f_size   ? TRR_SECTION_F   : 0) |
           (sh->box_size ? TRR_SECTION_BOX : 0);
}

// Walks the frame headers: where each frame is, and the time, step, sections and box it states.
// One positional read per frame, of the header and the box behind it. A trailing frame cut short -
// a run still being written - is left out rather than failing the file.
static size_t trr_scan(md_file_t file, int* out_natoms, md_array(int64_t)* offsets, md_array(double)* times, md_array(int64_t)* steps,
    md_array(uint32_t)* sections, md_array(float)* boxes, md_allocator_i* alloc) {
    const int64_t filesize = (int64_t)md_file_size(file);
    int64_t off = 0;
    int natoms = 0;
    size_t num_frames = 0;

    while (off < filesize) {
        uint8_t buf[TRR_MAX_HEADER_SIZE + DIM * DIM * 8];
        const size_t got = md_file_read_at(file, off, buf, sizeof(buf));
        md_xdr_t xdr = md_xdr_init(buf, got);
        trr_header_t sh;
        if (!trr_read_frame_header_buf(&xdr, &sh)) {
            if (num_frames == 0) return 0;
            MD_LOG_INFO("TRR: frame %zu has no readable header, the file is taken to end there", num_frames);
            break;
        }
        if (natoms == 0) {
            natoms = sh.natoms;
        } else if (sh.natoms != natoms) {
            MD_LOG_ERROR("TRR: frame %zu holds %i atoms, the first frame %i", num_frames, sh.natoms, natoms);
            return 0;
        }
        const int64_t end = off + (int64_t)xdr.pos + calc_framebytes(&sh);
        if (end > filesize) {
            MD_LOG_INFO("TRR: frame %zu is incomplete and left out", num_frames);
            break;
        }

        float box[DIM * DIM] = {0};
        if (sh.box_size) {
            if (sh.use_double) {
                double pv[DIM * DIM];
                md_xdr_read_f64_array(&xdr, pv, DIM * DIM);
                for (int i = 0; i < DIM * DIM; ++i) box[i] = (float)pv[i];
            } else {
                md_xdr_read_f32_array(&xdr, box, DIM * DIM);
            }
            if (!md_xdr_ok(&xdr)) {
                MD_LOG_ERROR("TRR: failed to read the box of frame %zu", num_frames);
                return 0;
            }
        }

        md_array_push(*offsets, off, alloc);
        md_array_push(*times, sh.t, alloc);
        md_array_push(*steps, (int64_t)sh.step, alloc);
        md_array_push(*sections, trr_sections(&sh), alloc);
        md_array_push_array(*boxes, box, DIM * DIM, alloc);
        num_frames += 1;
        off = end;
    }

    if (num_frames) {
        md_array_push(*offsets, off, alloc);
    }
    *out_natoms = natoms;
    return num_frames;
}

// Everything the scan learns about the file, short of the per atom data. Kept beside the
// trajectory as '<file>.cache' so the scan is paid once per file rather than once per load.
typedef struct trr_cache_t {
    md_run_cache_header_t header;
    int64_t*  frame_offsets;   // num_frames + 1: the last one is where the last frame ends
    double*   frame_times;     // ps
    int64_t*  frame_steps;
    uint32_t* frame_sections;  // TRR_SECTION_* present in each frame
    float*    frame_boxes;     // 9 per frame, nm, row major; zero without a box
} trr_cache_t;

static bool cache_read_block(md_file_t file, void** dst, size_t bytes, md_allocator_i* alloc) {
    *dst = md_alloc(alloc, bytes);
    return *dst && md_file_read(file, *dst, bytes) == bytes;
}

// On failure nothing is freed here: alloc is an arena owned by the caller, which goes as a whole.
static bool try_read_cache(trr_cache_t* cache, str_t path, md_allocator_i* alloc) {
    md_file_t file = {0};
    if (!md_run_cache_open(&file, &cache->header, path, MD_TRR_CACHE_MAGIC, MD_TRR_CACHE_VERSION)) {
        return false;
    }
    const size_t n = cache->header.num_frames;
    const bool ok =
        cache_read_block(file, (void**)&cache->frame_offsets,  (n + 1) * sizeof(int64_t), alloc) &&
        cache_read_block(file, (void**)&cache->frame_times,    n * sizeof(double),        alloc) &&
        cache_read_block(file, (void**)&cache->frame_steps,    n * sizeof(int64_t),       alloc) &&
        cache_read_block(file, (void**)&cache->frame_sections, n * sizeof(uint32_t),      alloc) &&
        cache_read_block(file, (void**)&cache->frame_boxes,    n * 9 * sizeof(float),     alloc) &&
        md_file_tell(file) == (int64_t)md_file_size(file);
    if (!ok) {
        MD_LOG_ERROR("TRR: the cache of '" STR_FMT "' is incomplete", STR_ARG(path));
    }
    md_file_close(&file);
    return ok;
}

static bool write_cache(const trr_cache_t* cache, str_t path, const md_file_info_t* scanned) {
    md_file_t file = {0};
    const size_t n = cache->header.num_frames;
    if (!md_run_cache_create(&file, path, scanned, MD_TRR_CACHE_MAGIC, MD_TRR_CACHE_VERSION, cache->header.num_atoms, n)) {
        return false;
    }
    const bool ok =
        md_file_write(file, cache->frame_offsets,  (n + 1) * sizeof(int64_t)) == (n + 1) * sizeof(int64_t) &&
        md_file_write(file, cache->frame_times,    n * sizeof(double))        == n * sizeof(double) &&
        md_file_write(file, cache->frame_steps,    n * sizeof(int64_t))       == n * sizeof(int64_t) &&
        md_file_write(file, cache->frame_sections, n * sizeof(uint32_t))      == n * sizeof(uint32_t) &&
        md_file_write(file, cache->frame_boxes,    n * 9 * sizeof(float))     == n * 9 * sizeof(float);
    if (!ok) {
        MD_LOG_ERROR("TRR: failed to write the cache of '" STR_FMT "'", STR_ARG(path));
    }
    md_file_close(&file);
    return ok;
}

// What the file is, from its cache when that is current and from a scan otherwise (writing the
// cache unless told not to). path must be canonical; every array is allocated from alloc.
static bool trr_index_load(trr_cache_t* cache, str_t path, uint32_t flags, md_allocator_i* alloc) {
    MEMSET(cache, 0, sizeof(*cache));

    md_file_t file = {0};
    if (!md_file_open(&file, path, MD_FILE_READ)) {
        MD_LOG_ERROR("TRR: Failed to open file '" STR_FMT "'", STR_ARG(path));
        return false;
    }

    bool result = false;
    md_file_info_t scanned = {0};
    md_file_info_extract(file, &scanned);

    if (try_read_cache(cache, path, alloc)) {
        result = true;
        goto done;
    }

    MEMSET(cache, 0, sizeof(*cache));
    md_array(int64_t)  offsets  = 0;
    md_array(double)   times    = 0;
    md_array(int64_t)  steps    = 0;
    md_array(uint32_t) sections = 0;
    md_array(float)    boxes    = 0;
    int natoms = 0;
    const size_t num_frames = trr_scan(file, &natoms, &offsets, &times, &steps, &sections, &boxes, alloc);
    if (num_frames == 0 || natoms <= 0) {
        MD_LOG_ERROR("TRR: '" STR_FMT "' holds no readable frames", STR_ARG(path));
        goto done;
    }

    cache->header.num_atoms  = (uint64_t)natoms;
    cache->header.num_frames = num_frames;
    cache->frame_offsets  = offsets;
    cache->frame_times    = times;
    cache->frame_steps    = steps;
    cache->frame_sections = sections;
    cache->frame_boxes    = boxes;

    if (!(flags & MD_RUN_FLAG_DISABLE_CACHE_WRITE)) {
        if (write_cache(cache, path, &scanned)) {
            MD_LOG_INFO("TRR: Successfully created cache file for '" STR_FMT "'", STR_ARG(path));
        }
    }
    result = true;

done:
    md_file_close(&file);
    return result;
}

// ### RUN ###

// The per atom sections a run publishes, in the order they sit in a frame
static const struct {
    const char* leaf;       // below the run
    uint32_t    section;
    float       scale;      // coordinates and velocities nm to Angstrom; forces as stored
} trr_run_sections[] = {
    { "atom/position", TRR_SECTION_X, 10.0f },
    { "atom/velocity", TRR_SECTION_V, 10.0f },
    { "atom/force",    TRR_SECTION_F, 1.0f  },
};

// <run>/atom/{position,velocity,force}. Everything it needs beyond the slice is in the table: the
// run is the path minus the leaf, the file, offset and size of the frame are its source attributes.
// Two reads per frame: the header, which says where the section is, then the section alone - a
// frame that carries all three reads a third of it for the positions.
static size_t trr_section_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    const md_system_t* sys = (const md_system_t*)user_data;
    ASSERT(sys);
    const md_attributes_t* attributes = &sys->attributes;

    // A whole temporal virtual attribute is refused before it gets here, so the frame is fixed.
    if (!slice || slice->num_idx == 0 || slice->num_idx > 2) {
        return 0;
    }

    int kind = -1;
    md_run_source_t src = {0};
    for (int i = 0; i < (int)ARRAY_SIZE(trr_run_sections); ++i) {
        if (md_run_source(&src, attributes, attr, str_from_cstr(trr_run_sections[i].leaf))) {
            kind = i;
            break;
        }
    }
    if (kind < 0) {
        return 0;
    }
    const int64_t* offsets = src.offset;
    const int64_t* sizes   = src.size;

    const uint32_t frame = slice->idx[0];
    const size_t num_atoms = attr->format.shape[1];
    const str_t file_path = src.path;

    uint8_t head[TRR_MAX_HEADER_SIZE];
    const size_t head_bytes = MIN(sizeof(head), (size_t)sizes[frame]);
    if (md_attribute_io_read_at(io, file_path, offsets[frame], head, head_bytes) != head_bytes) {
        MD_LOG_ERROR("TRR: Failed to read the header of frame %u from '" STR_FMT "'", frame, STR_ARG(file_path));
        return 0;
    }
    md_xdr_t xdr = md_xdr_init(head, head_bytes);
    trr_header_t sh;
    if (!trr_read_frame_header_buf(&xdr, &sh) || (size_t)sh.natoms != num_atoms) {
        return 0;
    }

    const size_t elem = sh.use_double ? 8 : 4;
    int64_t at = offsets[frame] + (int64_t)xdr.pos + sh.box_size + sh.vir_size + sh.pres_size;
    int section_size = sh.x_size;
    if (kind >= 1) { at += sh.x_size; section_size = sh.v_size; }
    if (kind >= 2) { at += sh.v_size; section_size = sh.f_size; }
    if ((size_t)section_size != num_atoms * 3 * elem) {
        // Not in this frame. The run only publishes what every frame carries, so this is a file
        // that changed underneath it.
        MD_LOG_ERROR("TRR: frame %u of '" STR_FMT "' has no %s", frame, STR_ARG(file_path), trr_run_sections[kind].leaf);
        return 0;
    }

    size_t first = 0;
    size_t count = num_atoms;
    if (slice->num_idx == 2) {
        if (slice->idx[1] >= num_atoms) return 0;
        first = slice->idx[1];
        count = 1;
    }
    if (cap != count * 3) {
        return 0;
    }

    md_temp_scope_t temp = md_temp_begin();
    size_t written = 0;
    const size_t bytes = count * 3 * elem;
    uint8_t* raw = md_temp_alloc(temp, bytes);
    if (raw && md_attribute_io_read_at(io, file_path, at + (int64_t)(first * 3 * elem), raw, bytes) == bytes) {
        float* out = (float*)dst;
        const float scale = trr_run_sections[kind].scale;
        if (sh.use_double) {
            for (size_t i = 0; i < count * 3; ++i) out[i] = (float)(md_xdr_load_f64(raw + i * 8) * scale);
        } else {
            for (size_t i = 0; i < count * 3; ++i) out[i] = md_xdr_load_f32(raw + i * 4) * scale;
        }
        written = cap;
    } else {
        MD_LOG_ERROR("TRR: Failed to read frame %u from '" STR_FMT "'", frame, STR_ARG(file_path));
    }
    md_temp_end(temp);
    return written;
}

// A section that is not in exactly the frames that have positions: its own group below the run,
// "<run>/trr/<name>", with its own frame axis and frame table, so each of its values sits at the
// time it was written.
static bool trr_publish_sparse_section(md_system_t* sys, str_t run, str_t path, const trr_cache_t* index, int kind, const md_unit_t unit, const str_t label,
    const md_attribute_virtual_t* virt, md_allocator_i* arena) {
    const size_t F = index->header.num_frames;
    const size_t N = index->header.num_atoms;
    const uint32_t section = trr_run_sections[kind].section;

    size_t n = 0;
    double*  times   = md_alloc(arena, F * sizeof(double));
    int64_t* steps   = md_alloc(arena, F * sizeof(int64_t));
    int64_t* offsets = md_alloc(arena, F * sizeof(int64_t));
    int64_t* sizes   = md_alloc(arena, F * sizeof(int64_t));
    for (size_t i = 0; i < F; ++i) {
        if (index->frame_sections[i] & section) {
            times[n]   = index->frame_times[i];
            steps[n]   = index->frame_steps[i];
            offsets[n] = index->frame_offsets[i];
            sizes[n]   = index->frame_offsets[i + 1] - index->frame_offsets[i];
            n += 1;
        }
    }
    if (n == 0) {
        return true;
    }

    // "<run>/trr/velocity" from "atom/velocity"
    const str_t leaf = str_from_cstr(trr_run_sections[kind].leaf);
    const str_t name = str_substr(leaf, 5, SIZE_MAX);
    char group_buf[512];
    const int group_len = snprintf(group_buf, sizeof(group_buf), STR_FMT "/trr/" STR_FMT, STR_ARG(run), STR_ARG(name));
    if (group_len <= 0 || group_len >= (int)sizeof(group_buf)) {
        return false;
    }
    const str_t group = { group_buf, (size_t)group_len };

    const md_attribute_format_t series_f64 = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)n } };
    const md_attribute_format_t series_i64 = { .type = MD_ATTRIBUTE_TYPE_I64, .components = 1, .rank = 1, .shape = { (uint32_t)n } };
    md_attributes_t* attributes = &sys->attributes;
    char buf[512];

    bool ok = md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), group, STR_LIT("time")), .format = series_f64, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
        .unit = md_unit_picosecond(), .label = STR_LIT("Time"), .data = times, .byte_size = n * sizeof(double)});
    ok = ok && md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), group, STR_LIT("step")), .format = series_i64, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
        .unit = md_unit_none(), .label = STR_LIT("Step"), .data = steps, .byte_size = n * sizeof(int64_t)});
    ok = ok && md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), group, STR_LIT("source/path")),
        .format = { .type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 0 },
        .unit = md_unit_none(), .data = &path, .byte_size = sizeof(str_t)});
    ok = ok && md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), group, STR_LIT("source/offset")), .format = series_i64, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
        .unit = md_unit_none(), .data = offsets, .byte_size = n * sizeof(int64_t)});
    ok = ok && md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), group, STR_LIT("source/size")), .format = series_i64, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
        .unit = md_unit_none(), .data = sizes, .byte_size = n * sizeof(int64_t)});
    ok = ok && md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), group, leaf),
        .format = { .type = MD_ATTRIBUTE_TYPE_F32, .components = 3, .rank = 2, .shape = { (uint32_t)n, (uint32_t)N } },
        .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = unit, .label = label, .virt = virt});
    return ok;
}

bool md_trr_system_publish_run(md_system_t* sys, str_t filename, str_t run, uint32_t flags) {
    ASSERT(sys);
    char path_buf[4096];
    const size_t path_len = md_path_write_canonical(path_buf, sizeof(path_buf), filename);
    const str_t path = {path_buf, path_len};

    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    bool result = false;

    trr_cache_t index;
    if (path_len == 0 || !trr_index_load(&index, path, flags, arena)) {
        goto done;
    }

    const size_t F = index.header.num_frames;
    const size_t N = index.header.num_atoms;

    // The run's frames are the frames with positions.
    size_t R = 0;
    double*  times   = md_alloc(arena, F * sizeof(double));
    int64_t* steps   = md_alloc(arena, F * sizeof(int64_t));
    int64_t* offsets = md_alloc(arena, F * sizeof(int64_t));
    int64_t* sizes   = md_alloc(arena, F * sizeof(int64_t));
    float*   boxes   = md_alloc(arena, F * 9 * sizeof(float));
    for (size_t i = 0; i < F; ++i) {
        if (!(index.frame_sections[i] & TRR_SECTION_X)) continue;
        times[R]   = index.frame_times[i];
        steps[R]   = index.frame_steps[i];
        offsets[R] = index.frame_offsets[i];
        sizes[R]   = index.frame_offsets[i + 1] - index.frame_offsets[i];
        for (int k = 0; k < 9; ++k) {
            boxes[R * 9 + k] = index.frame_boxes[i * 9 + k] * 10.0f;   // nm to Angstrom, as the coordinates
        }
        R += 1;
    }
    if (R == 0) {
        MD_LOG_ERROR("TRR: '" STR_FMT "' has no frame with coordinates", STR_ARG(path));
        goto done;
    }

    const md_attribute_virtual_t virt = { .provider = trr_section_provider, .user_data = sys };
    const md_run_desc_t desc = {
        .num_frames    = R,
        .num_atoms     = N,
        .time          = times,
        .time_unit     = md_unit_picosecond(),
        .step          = steps,
        .unitcell      = boxes,
        .source_path   = path,
        .source_offset = offsets,
        .source_size   = sizes,
        .position_virt = &virt,
    };
    if (!md_run_publish(sys, run, &desc)) {
        goto done;
    }

    const md_unit_t units[] = {
        md_unit_angstrom(),
        md_unit_div(md_unit_angstrom(), md_unit_picosecond()),
        md_unit_div(md_unit_div(md_unit_scl(md_unit_joule(), 1.0e3), md_unit_mole()), md_unit_nanometer()),
    };
    const str_t labels[] = { STR_LIT("Position"), STR_LIT("Velocity"), STR_LIT("Force") };
    char buf[512];

    // Velocities and forces: beside the positions when written in exactly the frames the positions
    // are, in a group of their own otherwise (GROMACS' nstvout and nstfout against nstxout).
    for (int k = 1; k < (int)ARRAY_SIZE(trr_run_sections); ++k) {
        const uint32_t section = trr_run_sections[k].section;
        bool any = false;
        bool aligned = true;
        for (size_t i = 0; i < F; ++i) {
            const bool has_x = (index.frame_sections[i] & TRR_SECTION_X) != 0;
            const bool has_s = (index.frame_sections[i] & section) != 0;
            any |= has_s;
            aligned &= (has_x == has_s);
        }
        if (!any) continue;

        bool ok;
        if (aligned) {
            ok = md_attributes_replace(&sys->attributes, &(md_attribute_desc_t){
                .path = md_run_path(buf, sizeof(buf), run, str_from_cstr(trr_run_sections[k].leaf)),
                .format = { .type = MD_ATTRIBUTE_TYPE_F32, .components = 3, .rank = 2, .shape = { (uint32_t)R, (uint32_t)N } },
                .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = units[k], .label = labels[k],
                .virt = &virt});
        } else {
            ok = trr_publish_sparse_section(sys, run, path, &index, k, units[k], labels[k], &virt, arena);
        }
        if (!ok) {
            MD_LOG_ERROR("TRR: failed to publish '" STR_FMT "' as '" STR_FMT "'", STR_ARG(path), STR_ARG(run));
            md_attributes_remove_prefix(&sys->attributes, run);
            goto done;
        }
    }
    result = true;

done:
    md_arena_allocator_destroy(arena);
    return result;
}
