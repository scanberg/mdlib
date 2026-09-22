#include <md_trr.h>
#include <md_xdr.h>
#include <md_system.h>

#include <md_util.h>
#include <md_trajectory.h>

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

#define MD_TRR_TRAJ_MAGIC 0x75CF81728AB71723
#define MD_TRR_TRAJ_READER_MAGIC 0x75CF81728AB71724

#define MD_TRR_CACHE_MAGIC 0x67b7cbab123452
#define MD_TRR_CACHE_VERSION 3

#define TRR_MAGIC 1993

#define TRR_MIN_HEADER_SIZE 72

#define DIM 3
typedef float matrix[DIM][DIM];

// This file cherry picks bits and pieces from the provided xdrfile_trr.c implementation.
// See xdrfile_trr.c for the copyright specific to that file.
// The implementation is modified from its original to use mdlibs allocator and error logging for clarity.

typedef struct trr_t {
    uint64_t magic;
    str_t filepath;
    int64_t* frame_offsets;
    md_trajectory_header_t header;
    md_allocator_i* allocator;
} trr_t;

typedef struct trr_reader_t {
    uint64_t magic;
    md_file_t file;
    const trr_t* traj;
    md_array(uint8_t) frame_data;
    md_allocator_i* arena;
} trr_reader_t;

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

// Reads the frame header at the current file position and leaves the file positioned right after it
static bool trr_read_frame_header(md_file_t file, trr_header_t* sh) {
    const int64_t beg = md_file_tell(file);
    uint8_t buf[TRR_MAX_HEADER_SIZE];
    const size_t bytes = md_file_read(file, buf, sizeof(buf));

    md_xdr_t xdr = md_xdr_init(buf, bytes);
    if (!trr_read_frame_header_buf(&xdr, sh)) {
        return false;
    }
    return md_file_seek(file, beg + (int64_t)xdr.pos, MD_FILE_BEG);
}

// Reads one per atom section (natoms 3-vectors) into planar arrays, scaled; skips it if dst is NULL
static bool trr_read_vec_section(md_xdr_t* xdr, const trr_header_t* sh, int section_size, float* dst[3], float scale, const char* what) {
    if (section_size == 0) {
        return true;
    }
    if (!dst) {
        if (!md_xdr_skip(xdr, (size_t)section_size)) {
            MD_LOG_ERROR("TRR: Failed to skip %s section in frame", what);
            return false;
        }
        return true;
    }
    const size_t elem = sh->use_double ? 8 : 4;
    const uint8_t* p = md_xdr_take(xdr, (size_t)sh->natoms * 3 * elem);
    if (!p) {
        MD_LOG_ERROR("TRR: Failed to read %s section in frame", what);
        return false;
    }
    for (int i = 0; i < sh->natoms; ++i, p += 3 * elem) {
        for (int k = 0; k < 3; ++k) {
            const double c = sh->use_double ? md_xdr_load_f64(p + k * 8) : md_xdr_load_f32(p + k * 4);
            dst[k][i] = (float)(c * scale);
        }
    }
    return true;
}

static bool trr_read_frame_data(md_xdr_t* xdr, const trr_header_t* sh, matrix box, float* x[3], float* v[3], float* f[3]) {
    if (sh->box_size != 0) {
        double pv[DIM * DIM];
        if (sh->use_double) {
            md_xdr_read_f64_array(xdr, pv, DIM * DIM);
        } else {
            float pvf[DIM * DIM];
            md_xdr_read_f32_array(xdr, pvf, DIM * DIM);
            for (int i = 0; i < DIM * DIM; ++i) pv[i] = pvf[i];
        }
        if (!md_xdr_ok(xdr)) {
            MD_LOG_ERROR("TRR: Failed to read frame box");
            return false;
        }
        if (box) {
            for (int i = 0; i < DIM; i++)
                for (int j = 0; j < DIM; j++)
                    box[i][j] = (float)pv[i * DIM + j];
        }
    }

    if (!md_xdr_skip(xdr, (size_t)sh->vir_size) || !md_xdr_skip(xdr, (size_t)sh->pres_size)) {
        return false;
    }

    // Coordinates and velocities nm -> Ångström, forces as stored
    return trr_read_vec_section(xdr, sh, sh->x_size, x, 10.0f, "coordinate")
        && trr_read_vec_section(xdr, sh, sh->v_size, v, 10.0f, "velocity")
        && trr_read_vec_section(xdr, sh, sh->f_size, f, 1.0f,  "force");
}

static int64_t trr_read_frame_offsets_and_times(md_file_t xd, md_array(int64_t)* offsets, md_array(double)* times, md_allocator_i* alloc) {
    trr_header_t sh = {0};

    int64_t filesize = md_file_size(xd);
    if (filesize <= 0) return 0;

    /* Go to file beg */
    if (!md_file_seek(xd, 0L, MD_FILE_BEG)) {
        return 0;
    }

    if (!trr_read_frame_header(xd, &sh)) {
        return 0;
    }

    int framebytes = calc_framebytes(&sh);
    int est_nframes = (int)(filesize / ((int64_t)(framebytes + TRR_MIN_HEADER_SIZE)) + 1); /* must be at least 1 for successful growth */

    /* Allocate memory for the frame index array */
    md_array_ensure(*offsets, (size_t)est_nframes, alloc);
    md_array_ensure(*times,   (size_t)est_nframes, alloc);
    md_array_push(*offsets, 0, alloc);

    int64_t num_frames = 1;

    while (1) {
        /* Skip `framebytes` */
        if (!md_file_seek(xd, (int64_t)(framebytes), MD_FILE_CUR)) {
            goto fail;
        }

        /* Store position in `offsets` */
        int64_t curr_offset = md_file_tell(xd);
        md_array_push(*offsets, curr_offset, alloc);
        md_array_push(*times,   sh.t, alloc);

        if (curr_offset == filesize) {
            return num_frames;
        }

        /* Read header and calculate how much to skip next time */
        if (!trr_read_frame_header(xd, &sh)) {
            goto fail;
        }
        framebytes = calc_framebytes(&sh);
        num_frames += 1;
    }

    return num_frames;

fail:
    md_array_free(*offsets, alloc);
    md_array_free(*times, alloc);
    return 0;
}

bool trr_get_header(struct md_trajectory_o* inst, md_trajectory_header_t* header) {
    trr_t* trr = (trr_t*)inst;
    ASSERT(trr);
    ASSERT(trr->magic == MD_TRR_TRAJ_MAGIC);
    ASSERT(header);

    *header = trr->header;
    return true;
}

// This is lowlevel cruft for enabling parallel loading and decoding of frames
static size_t trr_fetch_frame_data(const trr_t* trr, md_file_t file, int64_t frame_idx, void* frame_data_ptr) {
    ASSERT(trr);
    ASSERT(trr->magic == MD_TRR_TRAJ_MAGIC);

    if (!md_file_valid(file)) {
        MD_LOG_ERROR("TRR: File handle is NULL");
        return 0;
    }

    if (!trr->frame_offsets) {
        MD_LOG_ERROR("TRR: Frame offsets is empty");
        return 0;
    }

    if (frame_idx < 0 || (int64_t)trr->header.num_frames <= frame_idx) {
        MD_LOG_ERROR("TRR: Frame index is out of range");
        return 0;
    }

    const int64_t beg = trr->frame_offsets[frame_idx];
    const int64_t end = trr->frame_offsets[frame_idx + 1];
    const size_t frame_size = (size_t)(end - beg);

    if (frame_data_ptr) {
        if (!md_file_seek(file, beg, MD_FILE_BEG)) {
            MD_LOG_ERROR("TRR: Failed to seek to frame %lld", frame_idx);
            return 0;
        }
        const size_t bytes_read = md_file_read(file, frame_data_ptr, frame_size);
        if (bytes_read != frame_size) {
            MD_LOG_ERROR("TRR: Failed to read frame %lld", frame_idx);
            return 0;
        }
    }
    return frame_size;
}

// The optional per atom sections a TRR frame may carry beside its coordinates. Requested by handing
// over planar destination arrays; each 'present' flag reports whether this particular frame actually
// had that section, which is not something the caller can infer from the values it gets back.
typedef struct trr_frame_extras_t {
    float* v[3];        // NULL to skip the section
    float* f[3];
    bool   v_present;
    bool   f_present;
} trr_frame_extras_t;

static bool trr_decode_frame_data(const trr_t* trr, const void* frame_data_ptr, size_t frame_data_size, size_t* num_atoms, md_unitcell_t* cell, float* x, float* y, float* z, trr_frame_extras_t* extras) {
    ASSERT(frame_data_ptr);
    ASSERT(frame_data_size);

    bool result = true;

    if (trr->magic != MD_TRR_TRAJ_MAGIC) {
        MD_LOG_ERROR("TRR: Error when decoding frame coord, trr magic did not match");
        return false;
    }

    if ((x || y || z) && !(x && y && z)) {
        MD_LOG_ERROR("TRR: User supplied coordinates (x,y,z) cannot be partially supplied");
        return false;
    }

    md_xdr_t buf = md_xdr_init(frame_data_ptr, frame_data_size);

    // Get header
    trr_header_t sh;
    float box[3][3];
    float* coords[3] = { x, y, z };
    result = trr_read_frame_header_buf(&buf, &sh)
          && trr_read_frame_data(&buf, &sh, box, coords, extras ? extras->v : 0, extras ? extras->f : 0);
    if (result) {
        if (extras) {
            // The section sizes come off the frame header, so this is what the file says rather
            // than a guess from whether the values look written.
            extras->v_present = (sh.v_size != 0);
            extras->f_present = (sh.f_size != 0);
        }
        if (num_atoms) {
            *num_atoms = sh.natoms;
        }

        if (cell) {
            // @TODO: This scaling should be moved out of the core parts and into the loader which ties it with viamd.
            // nm -> Ångström
            for (int i = 0; i < 3; ++i) {
                box[i][0] *= 10.0f;
                box[i][1] *= 10.0f;
                box[i][2] *= 10.0f;
            }
            *cell = md_unitcell_from_basis_parameters(box[0][0], box[1][1], box[2][2], box[0][1], box[0][2], box[1][2]);
        }
    }

    return result;
}

static bool trr_reader_load_frame_raw(struct md_trajectory_reader_o* inst, int64_t frame_idx, size_t* num_atoms, md_unitcell_t* cell, float* x, float* y, float* z, trr_frame_extras_t* extras) {
    ASSERT(inst);

    trr_reader_t* reader = (trr_reader_t*)inst;
    ASSERT(reader->magic == MD_TRR_TRAJ_READER_MAGIC);

    const trr_t* trr = reader->traj;
    if (trr->magic != MD_TRR_TRAJ_MAGIC) {
        MD_LOG_ERROR("TRR: Error when decoding frame coord, trr magic did not match");
        return false;
    }

    bool result = false;
    const size_t frame_size = trr_fetch_frame_data(trr, reader->file, frame_idx, NULL);
    if (frame_size > 0) {
        md_array_ensure(reader->frame_data, frame_size, reader->arena);

        const size_t read_size = trr_fetch_frame_data(trr, reader->file, frame_idx, reader->frame_data);
        if (read_size != frame_size) {
            MD_LOG_ERROR("TRR: Failed to read the expected size");
            return false;
        }

        result = trr_decode_frame_data(trr, reader->frame_data, frame_size, num_atoms, cell, x, y, z, extras);
    }

    return result;
}

static void trr_trajectory_reader_free(struct md_trajectory_reader_i* reader) {
    if (!reader) {
        return;
    }

    trr_reader_t* inst = (trr_reader_t*)reader->inst;
    if (inst) {
        ASSERT(inst->magic == MD_TRR_TRAJ_READER_MAGIC);
        if (md_file_valid(inst->file)) {
            md_file_close(&inst->file);
        }
        md_arena_allocator_destroy(inst->arena);
    }

    MEMSET(reader, 0, sizeof(*reader));
}

// Adapts the raw reader to the state based interface. Everything the frame yields lands on the one
// state, which is what makes a metadata/coordinate mismatch unrepresentable here.
// @NOTE: state->frame is stamped by md_trajectory_reader_load_frame, not here.
static bool trr_reader_load_frame(struct md_trajectory_reader_o* inst, int64_t idx, md_system_state_t* state) {
    size_t num_atoms = 0;
    md_unitcell_t cell = {0};
    float* x = state ? state->x : NULL;
    float* y = state ? state->y : NULL;
    float* z = state ? state->z : NULL;

    // A TRR frame may carry velocities and forces beside its coordinates. They are per FRAME, so
    // they land on the state's own attribute table rather than on the system - which is the whole
    // reason that table exists. Only asked for when the state can hold them: a view state owns no
    // allocator and a raw coordinate fetch wants no extra decoding work.
    md_temp_scope_t temp = md_temp_begin();
    trr_frame_extras_t extras = {0};
    const size_t n = state ? state->num_atoms : 0;
    const bool want_extras = (state && state->attributes.alloc && n > 0);
    float* planar = NULL;

    if (want_extras) {
        planar = (float*)md_temp_alloc(temp, sizeof(float) * n * 6);
        if (planar) {
            for (int i = 0; i < 3; ++i) {
                extras.v[i] = planar + n * i;
                extras.f[i] = planar + n * (3 + i);
            }
        }
    }

    bool result = trr_reader_load_frame_raw(inst, idx, &num_atoms, &cell, x, y, z, planar ? &extras : NULL);
    if (result && state) {
        state->unitcell = cell;
        if (state->num_atoms == 0) {
            state->num_atoms = num_atoms;
        }

        if (planar && num_atoms == n) {
            // Planar out of the decoder, interleaved into the table: an attribute of 3 component
            // values stores xyz together, which is also what every consumer of a vector wants.
            float* xyz = (float*)md_temp_alloc(temp, sizeof(float) * n * 3);
            if (xyz) {
                if (extras.v_present) {
                    for (size_t i = 0; i < n; ++i) {
                        xyz[i * 3 + 0] = extras.v[0][i];
                        xyz[i * 3 + 1] = extras.v[1][i];
                        xyz[i * 3 + 2] = extras.v[2][i];
                    }
                    // The decoder already scaled nm/ps to Angstrom/ps, matching the coordinates.
                    md_attributes_publish_atom_column(&state->attributes, STR_LIT("atom/velocity"),
                        md_unit_div(md_unit_angstrom(), md_unit_picosecond()), 3, xyz, n);
                }
                if (extras.f_present) {
                    for (size_t i = 0; i < n; ++i) {
                        xyz[i * 3 + 0] = extras.f[0][i];
                        xyz[i * 3 + 1] = extras.f[1][i];
                        xyz[i * 3 + 2] = extras.f[2][i];
                    }
                    // Forces are NOT rescaled by the decoder, so they stay in GROMACS' kJ/mol/nm.
                    const md_unit_t kj_per_mol_nm = md_unit_div(
                        md_unit_div(md_unit_scl(md_unit_joule(), 1.0e3), md_unit_mole()), md_unit_nanometer());
                    md_attributes_publish_atom_column(&state->attributes, STR_LIT("atom/force"), kj_per_mol_nm, 3, xyz, n);
                }
            }
        }
    }

    md_temp_end(temp);
    return result;
}

static bool trr_trajectory_reader_init(md_trajectory_reader_i* reader, struct md_trajectory_o* traj_inst) {
    ASSERT(reader);
    ASSERT(traj_inst);

    trr_t* trr = (trr_t*)traj_inst;
    ASSERT(trr->magic == MD_TRR_TRAJ_MAGIC);

    md_file_t file = {0};
    if (!md_file_open(&file, trr->filepath, MD_FILE_READ)) {
        MD_LOG_ERROR("TRR: Failed to open '" STR_FMT "'", STR_ARG(trr->filepath));
        return false;
    }

    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    trr_reader_t* inst = md_alloc(arena, sizeof(trr_reader_t));
    MEMSET(inst, 0, sizeof(trr_reader_t));
    inst->magic = MD_TRR_TRAJ_READER_MAGIC;
    inst->file = file;
    inst->traj = trr;
    inst->arena = arena;

    MEMSET(reader, 0, sizeof(*reader));
    reader->inst = (struct md_trajectory_reader_o*)inst;
    reader->free = trr_trajectory_reader_free;
    reader->load_frame = trr_reader_load_frame;

    return true;
}

typedef struct trr_cache_t {
    md_trajectory_cache_header_t header;
    int64_t* frame_offsets;
    double*  frame_times;
} trr_cache_t;

static bool try_read_cache(trr_cache_t* cache, str_t cache_file, size_t traj_num_bytes, md_file_time_t traj_last_modified, md_allocator_i* alloc) {
    ASSERT(cache);
    ASSERT(alloc);

    bool result = false;
    md_file_t file = {0};
    if (md_file_open(&file, cache_file, MD_FILE_READ)) {
        if (md_file_read(file, &cache->header, sizeof(cache->header)) != sizeof(cache->header)) {
            MD_LOG_ERROR("TRR trajectory cache: failed to read header");
            goto done;
        }

        if (cache->header.magic != MD_TRR_CACHE_MAGIC) {
            MD_LOG_ERROR("TRR trajectory cache: magic was incorrect or corrupt");
            goto done;
        }
        if (cache->header.version != MD_TRR_CACHE_VERSION) {
            MD_LOG_INFO("TRR trajectory cache: version mismatch, expected %i, got %i", MD_TRR_CACHE_VERSION, (int)cache->header.version);
            goto done;
        }
        if (cache->header.num_bytes != traj_num_bytes) {
            MD_LOG_INFO("TRR trajectory cache: trajectory size mismatch, expected %zu, got %zu", traj_num_bytes, cache->header.num_bytes);
        }
        if (traj_last_modified != 0 && cache->header.last_modified != traj_last_modified) {
            MD_LOG_INFO("TRR trajectory cache: source file has been modified, cache is stale");
            goto done;
        }
        if (cache->header.num_atoms == 0) {
            MD_LOG_ERROR("TRR trajectory cache: num atoms was zero");
            goto done;
        }
        if (cache->header.num_frames == 0) {
            MD_LOG_ERROR("TRR trajectory cache: num frames was zero");
            goto done;
        }

        const size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
        cache->frame_offsets = md_alloc(alloc, offset_bytes);
        if (md_file_read(file, cache->frame_offsets, offset_bytes) != offset_bytes) {
            MD_LOG_ERROR("TRR trajectory cache: Failed to read offset data");
            md_free(alloc, cache->frame_offsets, offset_bytes);
            goto done;
        }

        const size_t time_bytes = cache->header.num_frames * sizeof(double);
        cache->frame_times = md_alloc(alloc, time_bytes);
        if (md_file_read(file, cache->frame_times, time_bytes) != time_bytes) {
            MD_LOG_ERROR("TRR trajectory cache: times are incomplete");
            md_free(alloc, cache->frame_offsets, offset_bytes);
            md_free(alloc, cache->frame_times, time_bytes);
            goto done;
        }

        // Test position in file, we expect to be at the end of the file
        if (md_file_tell(file) != (int64_t)md_file_size(file)) {
            MD_LOG_ERROR("TRR trajectory cache: file position was not at the end of the file");
            md_free(alloc, cache->frame_offsets, offset_bytes);
            md_free(alloc, cache->frame_times, time_bytes);
            goto done;
        }

        result = true;
    done:
        md_file_close(&file);
    }
    return result;
}

static bool write_cache(const trr_cache_t* cache, str_t cache_file) {
    bool result = false;

    md_file_t file = {0};
    if (!md_file_open(&file, cache_file, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
        MD_LOG_INFO("TRR trajectory cache: could not open file '"STR_FMT"'", STR_ARG(cache_file));
        return false;
    }

    if (md_file_write(file, &cache->header, sizeof(cache->header)) != sizeof(cache->header)) {
        MD_LOG_ERROR("TRR trajectory cache: failed to write header");
        goto done;
    }

    const size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
    if (md_file_write(file, cache->frame_offsets, offset_bytes) != offset_bytes) {
        MD_LOG_ERROR("TRR trajectory cache: failed to write offsets");
        goto done;
    }

    const size_t time_bytes = cache->header.num_frames * sizeof(double);
    if (md_file_write(file, cache->frame_times, time_bytes) != time_bytes) {
        MD_LOG_ERROR("TRR trajectory cache: failed to write times");
        goto done;
    }

    result = true;

done:
    md_file_close(&file);
    return result;
}

void md_trr_trajectory_free(md_trajectory_i* traj) {
    ASSERT(traj);
    ASSERT(traj->inst);
    trr_t* trr = (trr_t*)traj->inst;
    if (trr->magic != MD_TRR_TRAJ_MAGIC) {
        MD_LOG_ERROR("TRR: Cannot free trajectory, is not a valid trr trajectory.");
        ASSERT(false);
        return;
    }
    md_arena_allocator_destroy(trr->allocator);
}

md_trajectory_i* md_trr_trajectory_create(str_t filename, md_allocator_i* ext_alloc, uint32_t flags) {
    ASSERT(ext_alloc);
    md_allocator_i* alloc = md_arena_allocator_create(ext_alloc, MEGABYTES(1));

    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        md_arena_allocator_destroy(alloc);
        return NULL;
    }

    const int64_t filesize = md_file_size(file);

    md_file_info_t file_info = {0};
    md_file_info_extract_from_path(filename, &file_info);

    trr_header_t sh;
    if (!trr_read_frame_header(file, &sh)) {
        goto fail;
    }

    if (sh.natoms == 0) {
        MD_LOG_ERROR("TRR: Number of atoms is zero");
        goto fail;
    }

    {
        md_temp_scope_t temp = md_temp_begin_avoid(ext_alloc);
        md_allocator_i* temp_arena = md_temp_allocator(temp);
        md_strb_t sb = md_strb_create(temp_arena);
        md_strb_push_str(&sb, filename);
        md_strb_push_cstr(&sb, ".cache");
        str_t cache_file = md_strb_to_str(sb);

        trr_cache_t cache = {0};
        if (!try_read_cache(&cache, cache_file, filesize, file_info.modified_time, alloc)) {
            cache.header.magic = MD_TRR_CACHE_MAGIC;
            cache.header.version = MD_TRR_CACHE_VERSION;
            cache.header.num_bytes = filesize;
            cache.header.num_atoms = sh.natoms;
            cache.header.last_modified = file_info.modified_time;
            cache.header.num_frames = trr_read_frame_offsets_and_times(file, &cache.frame_offsets, &cache.frame_times, alloc);
            if (!cache.header.num_frames) {
                goto fail;
            }

            if (!(flags & MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE)) {
                // If we fail to write the cache, that's ok, we can inform about it, but do not halt
                if (write_cache(&cache, cache_file)) {
                    MD_LOG_INFO("TRR: Successfully created cache file for '" STR_FMT "'", STR_ARG(cache_file));
                }
            }
        }

        if (!cache.frame_offsets || !cache.frame_times) {
            MD_LOG_DEBUG("TRR: frame offsets or frame times was empty");
            goto fail;
        }

        void* mem = md_alloc(alloc, sizeof(md_trajectory_i) + sizeof(trr_t));
        ASSERT(mem);
        MEMSET(mem, 0, sizeof(md_trajectory_i) + sizeof(trr_t));

        md_trajectory_i* traj = mem;
        trr_t* trr = (trr_t*)(traj + 1);

        trr->magic = MD_TRR_TRAJ_MAGIC;
        trr->allocator = alloc;
        trr->filepath = str_copy(filename, alloc);
        trr->frame_offsets = cache.frame_offsets;

        trr->header = (md_trajectory_header_t) {
            .num_frames = cache.header.num_frames,
            .num_atoms = sh.natoms,
            .time_unit = md_unit_picosecond(),
            .frame_times = cache.frame_times,
        };

        traj->inst = (struct md_trajectory_o*)trr;
        traj->free = md_trr_trajectory_free;
        traj->get_header = trr_get_header;
        traj->init_reader = trr_trajectory_reader_init;

        md_file_close(&file);
        md_temp_end(temp);
        return traj;
    }

fail:
    md_file_close(&file);
    md_arena_allocator_destroy(alloc);
    return NULL;
}

// Attach convenience wrapper: create trajectory and attach to system
bool md_trr_attach_from_file(struct md_system_t* sys, str_t filename, uint32_t flags) {
    ASSERT(sys);
    if (!sys->alloc) {
        MD_LOG_ERROR("System allocator not set");
        return false;
    }

    md_trajectory_i* traj = md_trr_trajectory_create(filename, sys->alloc, flags);
    if (!traj) return false;
    md_system_attach_trajectory(sys, traj);
    return true;
}
