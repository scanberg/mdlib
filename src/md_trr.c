#include <md_trr.h>
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

// =====================================================================
// Internal helpers: big-endian primitive extraction
// TRR files use the XDR encoding which is big-endian.
// =====================================================================

static inline bool trr_read_int32(int32_t* out, md_file_t f) {
    uint32_t raw;
    if (md_file_read(f, &raw, sizeof(raw)) != sizeof(raw)) return false;
#if __LITTLE_ENDIAN__
    raw = BSWAP32(raw);
#endif
    MEMCPY(out, &raw, sizeof(int32_t));
    return true;
}

static inline bool trr_read_float(float* out, md_file_t f) {
    uint32_t raw;
    if (md_file_read(f, &raw, sizeof(raw)) != sizeof(raw)) return false;
#if __LITTLE_ENDIAN__
    raw = BSWAP32(raw);
#endif
    MEMCPY(out, &raw, sizeof(float));
    return true;
}

static inline bool trr_read_double(double* out, md_file_t f) {
    uint64_t raw;
    if (md_file_read(f, &raw, sizeof(raw)) != sizeof(raw)) return false;
#if __LITTLE_ENDIAN__
    raw = BSWAP64(raw);
#endif
    MEMCPY(out, &raw, sizeof(double));
    return true;
}

// XDR string on disk: uint32 char-count (big-endian), then `count` bytes,
// padded up to the next 4-byte boundary. Returns false on error.
static inline bool trr_read_string(char* buf, int maxlen, md_file_t f) {
    int32_t len;
    if (!trr_read_int32(&len, f)) return false;
    if (len <= 0 || len > maxlen) return false;
    if (md_file_read(f, buf, (size_t)len) != (size_t)len) return false;
    // consume padding bytes
    int pad = (4 - (len & 3)) & 3;
    if (pad > 0) {
        char tmp[3];
        if (md_file_read(f, tmp, (size_t)pad) != (size_t)pad) return false;
    }
    return true;
}

// =====================================================================
// Buffer-cursor helpers for in-memory decoding (replaces xdrfile_mem).
// =====================================================================

typedef struct trr_buf_t {
    const uint8_t* data;
    size_t         size;
    size_t         pos;
} trr_buf_t;

static inline bool buf_read_int32(int32_t* out, trr_buf_t* b) {
    if (b->pos + 4 > b->size) return false;
    uint32_t raw;
    MEMCPY(&raw, b->data + b->pos, 4);
    b->pos += 4;
#if __LITTLE_ENDIAN__
    raw = BSWAP32(raw);
#endif
    MEMCPY(out, &raw, sizeof(int32_t));
    return true;
}

static inline bool buf_read_float(float* out, trr_buf_t* b) {
    if (b->pos + 4 > b->size) return false;
    uint32_t raw;
    MEMCPY(&raw, b->data + b->pos, 4);
    b->pos += 4;
#if __LITTLE_ENDIAN__
    raw = BSWAP32(raw);
#endif
    MEMCPY(out, &raw, sizeof(float));
    return true;
}

static inline bool buf_read_double(double* out, trr_buf_t* b) {
    if (b->pos + 8 > b->size) return false;
    uint64_t raw;
    MEMCPY(&raw, b->data + b->pos, 8);
    b->pos += 8;
#if __LITTLE_ENDIAN__
    raw = BSWAP64(raw);
#endif
    MEMCPY(out, &raw, sizeof(double));
    return true;
}

static inline bool buf_read_string(char* buf, int maxlen, trr_buf_t* b) {
    int32_t len;
    if (!buf_read_int32(&len, b)) return false;
    if (len <= 0 || len > maxlen) return false;
    if (b->pos + (size_t)len > b->size) return false;
    MEMCPY(buf, b->data + b->pos, (size_t)len);
    b->pos += (size_t)len;
    int pad = (4 - (len & 3)) & 3;
    if (b->pos + (size_t)pad > b->size) return false;
    b->pos += (size_t)pad;
    return true;
}

static inline bool buf_skip(trr_buf_t* b, size_t count) {
    if (b->pos + count > b->size) return false;
    b->pos += count;
    return true;
}

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

static bool trr_read_frame_header(md_file_t xd, trr_header_t* sh) {
    ASSERT(sh);

    const char version[] = "GMX_trn_file";
    char buf[128];
    int32_t magic, slen, nflsz;

    if (!trr_read_int32(&magic, xd)) {
        MD_LOG_ERROR("TRR: Failed to read header magic number");
        return false;
    }
    if (magic != TRR_MAGIC) {
        MD_LOG_ERROR("TRR: Magic number did not match");
        return false;
    }

    if (!trr_read_int32(&slen, xd)) {
        MD_LOG_ERROR("TRR: Failed to read header version string length");
        return false;
    }
    if (slen != (int32_t)sizeof(version)) {
        MD_LOG_ERROR("TRR: Incorrect version string length");
        return false;
    }
    if (!trr_read_string(buf, (int)sizeof(buf), xd)) {
        MD_LOG_ERROR("TRR: Failed to read header version string");
        return false;
    }

    int32_t fields[11];
    for (int i = 0; i < 11; ++i) {
        if (!trr_read_int32(&fields[i], xd)) {
            MD_LOG_ERROR("TRR: Failed to read header fields");
            return false;
        }
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

    if (!n_float_size(sh, &nflsz)) {
        return false;
    }

    sh->use_double = (nflsz == sizeof(double));

    int32_t step, nre;
    if (!trr_read_int32(&step, xd)) {
        MD_LOG_ERROR("TRR: Failed to read header step");
        return false;
    }
    sh->step = step;

    if (!trr_read_int32(&nre, xd)) {
        MD_LOG_ERROR("TRR: Failed to read header \"nre\"");
        return false;
    }
    sh->nre = nre;

    if (sh->use_double) {
        if (!trr_read_double(&sh->t, xd)) {
            MD_LOG_ERROR("TRR: Failed to read header t");
            return false;
        }
        if (!trr_read_double(&sh->lambda, xd)) {
            MD_LOG_ERROR("TRR: Failed to read header lambda");
            return false;
        }
    } else {
        float tf, lf;
        if (!trr_read_float(&tf, xd)) {
            MD_LOG_ERROR("TRR: Failed to read header t");
            return false;
        }
        if (!trr_read_float(&lf, xd)) {
            MD_LOG_ERROR("TRR: Failed to read header lambda");
            return false;
        }
        sh->t = tf;
        sh->lambda = lf;
    }

    return true;
}

static bool trr_read_frame_header_buf(trr_buf_t* buf, trr_header_t* sh) {
    ASSERT(buf);
    ASSERT(sh);

    const char version[] = "GMX_trn_file";
    char strbuf[128];
    int32_t magic, slen, nflsz;

    if (!buf_read_int32(&magic, buf)) {
        MD_LOG_ERROR("TRR: Failed to read header magic number");
        return false;
    }
    if (magic != TRR_MAGIC) {
        MD_LOG_ERROR("TRR: Magic number did not match");
        return false;
    }

    if (!buf_read_int32(&slen, buf)) {
        MD_LOG_ERROR("TRR: Failed to read header version string length");
        return false;
    }
    if (slen != (int32_t)sizeof(version)) {
        MD_LOG_ERROR("TRR: Incorrect version string length");
        return false;
    }
    if (!buf_read_string(strbuf, (int)sizeof(strbuf), buf)) {
        MD_LOG_ERROR("TRR: Failed to read header version string");
        return false;
    }

    int32_t fields[11];
    for (int i = 0; i < 11; ++i) {
        if (!buf_read_int32(&fields[i], buf)) {
            MD_LOG_ERROR("TRR: Failed to read header fields");
            return false;
        }
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

    if (!n_float_size(sh, &nflsz)) {
        return false;
    }
    sh->use_double = (nflsz == sizeof(double));

    int32_t step, nre;
    if (!buf_read_int32(&step, buf)) { MD_LOG_ERROR("TRR: Failed to read header step"); return false; }
    if (!buf_read_int32(&nre,  buf)) { MD_LOG_ERROR("TRR: Failed to read header \"nre\""); return false; }
    sh->step = step;
    sh->nre  = nre;

    if (sh->use_double) {
        if (!buf_read_double(&sh->t,      buf)) { MD_LOG_ERROR("TRR: Failed to read header t");      return false; }
        if (!buf_read_double(&sh->lambda, buf)) { MD_LOG_ERROR("TRR: Failed to read header lambda"); return false; }
    } else {
        float tf, lf;
        if (!buf_read_float(&tf, buf)) { MD_LOG_ERROR("TRR: Failed to read header t");      return false; }
        if (!buf_read_float(&lf, buf)) { MD_LOG_ERROR("TRR: Failed to read header lambda"); return false; }
        sh->t      = tf;
        sh->lambda = lf;
    }
    return true;
}

static bool trr_read_frame_data(trr_buf_t* buf, const trr_header_t* sh, matrix box, float* x[3], float* v[3], float* f[3]) {
    /* Double */
    if (sh->use_double) {
        if (sh->box_size != 0) {
            double pvd[DIM * DIM];
            for (int i = 0; i < DIM * DIM; ++i) {
                if (!buf_read_double(&pvd[i], buf)) {
                    MD_LOG_ERROR("TRR: Failed to read frame box");
                    return false;
                }
            }
            if (box) {
                for (int i = 0; i < DIM; i++)
                    for (int j = 0; j < DIM; j++)
                        box[i][j] = (float)pvd[i * DIM + j];
            }
        }

        if (sh->vir_size != 0) {
            if (!buf_skip(buf, (size_t)sh->vir_size)) return false;
        }
        if (sh->pres_size != 0) {
            if (!buf_skip(buf, (size_t)sh->pres_size)) return false;
        }

        if (sh->x_size != 0) {
            if (x) {
                for (int i = 0; i < sh->natoms; ++i) {
                    double c[3];
                    if (!buf_read_double(&c[0], buf) || !buf_read_double(&c[1], buf) || !buf_read_double(&c[2], buf)) {
                        MD_LOG_ERROR("TRR: Failed to read coordinate entry in frame");
                        return false;
                    }
                    x[0][i] = (float)(c[0] * 10);
                    x[1][i] = (float)(c[1] * 10);
                    x[2][i] = (float)(c[2] * 10);
                }
            } else {
                if (!buf_skip(buf, (size_t)sh->x_size)) {
                    MD_LOG_ERROR("TRR: Failed to skip coordinate section in frame");
                    return false;
                }
            }
        }
        if (sh->v_size != 0) {
            if (v) {
                for (int i = 0; i < sh->natoms; ++i) {
                    double c[3];
                    if (!buf_read_double(&c[0], buf) || !buf_read_double(&c[1], buf) || !buf_read_double(&c[2], buf)) {
                        MD_LOG_ERROR("TRR: Failed to read velocity entry in frame");
                        return false;
                    }
                    v[0][i] = (float)(c[0] * 10);
                    v[1][i] = (float)(c[1] * 10);
                    v[2][i] = (float)(c[2] * 10);
                }
            } else {
                if (!buf_skip(buf, (size_t)sh->v_size)) {
                    MD_LOG_ERROR("TRR: Failed to skip velocity section in frame");
                    return false;
                }
            }
        }
        if (sh->f_size != 0) {
            if (f) {
                for (int i = 0; i < sh->natoms; ++i) {
                    double c[3];
                    if (!buf_read_double(&c[0], buf) || !buf_read_double(&c[1], buf) || !buf_read_double(&c[2], buf)) {
                        MD_LOG_ERROR("TRR: Failed to read force entry in frame");
                        return false;
                    }
                    f[0][i] = (float)c[0];
                    f[1][i] = (float)c[1];
                    f[2][i] = (float)c[2];
                }
            } else {
                if (!buf_skip(buf, (size_t)sh->f_size)) {
                    MD_LOG_ERROR("TRR: Failed to skip force section in frame");
                    return false;
                }
            }
        }
    }
    else
    /* Float */
    {
        if (sh->box_size != 0) {
            float pvf[DIM * DIM];
            for (int i = 0; i < DIM * DIM; ++i) {
                if (!buf_read_float(&pvf[i], buf)) {
                    MD_LOG_ERROR("TRR: Failed to read frame box");
                    return false;
                }
            }
            if (box) {
                for (int i = 0; i < DIM; i++)
                    for (int j = 0; j < DIM; j++)
                        box[i][j] = pvf[i * DIM + j];
            }
        }

        if (sh->vir_size != 0) {
            if (!buf_skip(buf, (size_t)sh->vir_size)) return false;
        }
        if (sh->pres_size != 0) {
            if (!buf_skip(buf, (size_t)sh->pres_size)) return false;
        }

        if (sh->x_size != 0) {
            if (x) {
                for (int i = 0; i < sh->natoms; ++i) {
                    float c[3];
                    if (!buf_read_float(&c[0], buf) || !buf_read_float(&c[1], buf) || !buf_read_float(&c[2], buf)) {
                        MD_LOG_ERROR("TRR: Failed to read coordinate entry in frame");
                        return false;
                    }
                    x[0][i] = c[0] * 10;
                    x[1][i] = c[1] * 10;
                    x[2][i] = c[2] * 10;
                }
            } else {
                if (!buf_skip(buf, (size_t)sh->x_size)) {
                    MD_LOG_ERROR("TRR: Failed to skip coordinate section in frame");
                    return false;
                }
            }
        }
        if (sh->v_size != 0) {
            if (v) {
                for (int i = 0; i < sh->natoms; ++i) {
                    float c[3];
                    if (!buf_read_float(&c[0], buf) || !buf_read_float(&c[1], buf) || !buf_read_float(&c[2], buf)) {
                        MD_LOG_ERROR("TRR: Failed to read velocity entry in frame");
                        return false;
                    }
                    v[0][i] = c[0] * 10;
                    v[1][i] = c[1] * 10;
                    v[2][i] = c[2] * 10;
                }
            } else {
                if (!buf_skip(buf, (size_t)sh->v_size)) {
                    MD_LOG_ERROR("TRR: Failed to skip velocity section in frame");
                    return false;
                }
            }
        }
        if (sh->f_size != 0) {
            if (f) {
                for (int i = 0; i < sh->natoms; ++i) {
                    float c[3];
                    if (!buf_read_float(&c[0], buf) || !buf_read_float(&c[1], buf) || !buf_read_float(&c[2], buf)) {
                        MD_LOG_ERROR("TRR: Failed to read force entry in frame");
                        return false;
                    }
                    f[0][i] = c[0];
                    f[1][i] = c[1];
                    f[2][i] = c[2];
                }
            } else {
                if (!buf_skip(buf, (size_t)sh->f_size)) {
                    MD_LOG_ERROR("TRR: Failed to skip force section in frame");
                    return false;
                }
            }
        }
    }
    return true;
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

static bool trr_decode_frame_data(const trr_t* trr, const void* frame_data_ptr, size_t frame_data_size, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
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

    trr_buf_t buf = { (const uint8_t*)frame_data_ptr, frame_data_size, 0 };

    // Get header
    trr_header_t sh;
    float box[3][3];
    float* coords[3] = { x, y, z };
    result = trr_read_frame_header_buf(&buf, &sh) && trr_read_frame_data(&buf, &sh, box, coords, 0, 0);
    if (result && header) {
        // @TODO: This scaling should be moved out of the core parts and into the loader which ties it with viamd.
        // nm -> Ångström
        for (int i = 0; i < 3; ++i) {
            box[i][0] *= 10.0f;
            box[i][1] *= 10.0f;
            box[i][2] *= 10.0f;
        }
        header->num_atoms = sh.natoms;
        header->index = sh.step;
        header->timestamp = sh.t;
        header->unitcell = md_unitcell_from_basis_parameters(box[0][0], box[1][1], box[2][2], box[0][1], box[0][2], box[1][2]);
    }

    return result;
}

static bool trr_reader_load_frame(struct md_trajectory_reader_o* inst, int64_t frame_idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
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

        result = trr_decode_frame_data(trr, reader->frame_data, frame_size, header, x, y, z);
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
        md_strb_t sb = md_strb_create(md_get_temp_arena());
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
