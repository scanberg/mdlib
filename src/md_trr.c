#include <md_trr.h>

#include <md_util.h>
#include <md_trajectory.h>

#include <core/md_common.h>
#include <core/md_array.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_vec_math.h>

#include <xdrfile.h>

#include <string.h>
#include <stdio.h>

#define MD_TRR_TRAJ_MAGIC 0x75CF81728AB71723

#define MD_TRR_CACHE_MAGIC 0x67b7cbab
#define MD_TRR_CACHE_VERSION 2

#define TRR_MAGIC 1993

#define TRR_MIN_HEADER_SIZE 72

// This file cherry picks bits and pieces from the provided xdrfile_trr.c implementation.
// See xdrfile_trr.c for the copyright specific to that file.
// The implementation is modified from its original to use mdlibs allocator and error logging for clarity.

typedef struct trr_t {
    uint64_t magic;
    XDRFILE* file;
    int64_t* frame_offsets;
    md_trajectory_header_t header;
    md_allocator_i* allocator;
    md_mutex_t mutex;
} trr_t;

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

static bool trr_read_frame_header(XDRFILE* xd, trr_header_t* sh) {
    ASSERT(xd);
    ASSERT(sh);

    const char version[] = "GMX_trn_file";
    char buf[128];
    int magic, slen, nflsz;

    if (xdrfile_read_int(&magic, 1, xd) != 1) {
        MD_LOG_ERROR("TRR: Failed to read header magic number");
        return false;
    }
    if (magic != TRR_MAGIC) {
        MD_LOG_ERROR("TRR: Magic number did not match");
        return false;
    }

    if (xdrfile_read_int(&slen, 1, xd) != 1) {
        MD_LOG_ERROR("TRR: Failed to read header version string length");
        return false;
    }
    if (slen != sizeof(version)) {
        MD_LOG_ERROR("TRR: Incorrect version string length");
        return false;
    }
    if (xdrfile_read_string(buf, sizeof(buf), xd) <= 0) {
        MD_LOG_ERROR("TRR: Failed to read header version string");
        return false;
    }

    if (xdrfile_read_int(&sh->ir_size,   1, xd) != 1 ||
        xdrfile_read_int(&sh->e_size,    1, xd) != 1 ||
        xdrfile_read_int(&sh->box_size,  1, xd) != 1 ||
        xdrfile_read_int(&sh->vir_size,  1, xd) != 1 ||
        xdrfile_read_int(&sh->pres_size, 1, xd) != 1 ||
        xdrfile_read_int(&sh->top_size,  1, xd) != 1 ||
        xdrfile_read_int(&sh->sym_size,  1, xd) != 1 ||
        xdrfile_read_int(&sh->x_size,    1, xd) != 1 ||
        xdrfile_read_int(&sh->v_size,    1, xd) != 1 ||
        xdrfile_read_int(&sh->f_size,    1, xd) != 1 ||
        xdrfile_read_int(&sh->natoms,    1, xd) != 1)
    {
        MD_LOG_ERROR("TRR: Failed to read header fields");
        return false;
    }

    if (!n_float_size(sh, &nflsz)) {
        return false;
    }

    sh->use_double = (nflsz == sizeof(double));

    if (xdrfile_read_int(&sh->step, 1, xd) != 1) {
        MD_LOG_ERROR("TRR: Failed to read header step");
        return false;
    }
    if (xdrfile_read_int(&sh->nre, 1, xd) != 1) {
        MD_LOG_ERROR("TRR: Failed to read header \"nre\"");
        return false;
    }

    if (sh->use_double) {
        if (xdrfile_read_double(&sh->t, 1, xd) != 1) {
            MD_LOG_ERROR("TRR: Failed to read header t");
            return false;
        }
        if (xdrfile_read_double(&sh->lambda, 1, xd) != 1) {
            MD_LOG_ERROR("TRR: Failed to read header lambda");
            return false;
        }
    } else {
        float tf, lf;
        if (xdrfile_read_float(&tf, 1, xd) != 1) {
            MD_LOG_ERROR("TRR: Failed to read header t");
            return false;
        }
        if (xdrfile_read_float(&lf, 1, xd) != 1) {
            MD_LOG_ERROR("TRR: Failed to read header lambda");
            return false;
        }
        sh->t = tf;
        sh->lambda = lf;
    }

    return true;
}

static bool trr_read_frame_data(XDRFILE* xd, const trr_header_t* sh, matrix box, md_vec3_soa_t* x, md_vec3_soa_t* v, md_vec3_soa_t* f) {
    bool allocate_temp = (sh->x_size != 0) || (sh->v_size != 0) || (sh->f_size != 0);

    /* Double */
    if (sh->use_double) {
        double pvd[DIM * DIM];
        double* dx = NULL;

        if (sh->box_size != 0) {
            if (xdrfile_read_double(pvd, DIM * DIM, xd) == DIM * DIM) {
                if (box) {
                    for (int i = 0; (i < DIM); i++) {
                        for (int j = 0; (j < DIM); j++) {
                            box[i][j] = (float)pvd[i * DIM + j];
                        }
                    }
                }
            } else {
                MD_LOG_ERROR("TRR: Failed to read frame box");
                return false;
            }
        }

        if (sh->vir_size != 0) {
            if (xdrfile_read_double(pvd, DIM * DIM, xd) != DIM * DIM) {
                return false;
            }
        }

        if (sh->pres_size != 0) {
            if (xdrfile_read_double(pvd, DIM * DIM, xd) != DIM * DIM) {
                return false;
            }
        }

        if (allocate_temp) {
            dx = (double*)calloc(sh->natoms * DIM, sizeof(dx[0]));
            if (NULL == dx) {
                return false;
            }
        }
        if (sh->x_size != 0) {
            if (x) {
                double c[3];
                for (int i = 0; i < sh->natoms; ++i) {
                    if (xdrfile_read_double(c, DIM, xd) == DIM) {
                        // @TODO: This scaling should be moved out of the core parts and into the loader which ties it with viamd.
                        x->x[i] = (float)(c[0] * 10);
                        x->y[i] = (float)(c[1] * 10);
                        x->z[i] = (float)(c[2] * 10);
                    } else {
                        MD_LOG_ERROR("TRR: Failed to read coordinate entry in frame");
                        return false;
                    }
                }
            } else {
                if (xdr_seek(xd, sh->x_size, SEEK_CUR) != exdrOK) {
                    MD_LOG_ERROR("TRR: Failed to skip coordinate section in frame");
                    return false;
                }
            }
        }
        if (sh->v_size != 0) {
            if (v) {
                double c[3];
                for (int i = 0; i < sh->natoms; ++i) {
                    if (xdrfile_read_double(c, DIM, xd) == DIM) {
                        // @TODO: This scaling should be moved out of the core parts and into the loader which ties it with viamd.
                        v->x[i] = (float)(c[0] * 10);
                        v->y[i] = (float)(c[1] * 10);
                        v->z[i] = (float)(c[2] * 10);
                    } else {
                        MD_LOG_ERROR("TRR: Failed to read velocity entry in frame");
                        return false;
                    }
                }
            } else {
                if (xdr_seek(xd, sh->v_size, SEEK_CUR) != exdrOK) {
                    MD_LOG_ERROR("TRR: Failed to skip velocity section in frame");
                    return false;
                }
            }
        }
        if (sh->f_size != 0) {
            if (f) {
                double c[3];
                for (int i = 0; i < sh->natoms; ++i) {
                    if (xdrfile_read_double(c, DIM, xd) == DIM) {
                        f->x[i] = (float)c[0];
                        f->y[i] = (float)c[1];
                        f->z[i] = (float)c[2];
                    } else {
                        MD_LOG_ERROR("TRR: Failed to read force entry in frame");
                        return false;
                    }
                }
            } else {
                if (xdr_seek(xd, sh->f_size, SEEK_CUR) != exdrOK) {
                    MD_LOG_ERROR("TRR: Failed to skip force section in frame");
                    return false;
                }
            }
        }
        if (dx) {
            free(dx);
        }
    }
    else
    /* Float */
    {
        float pvf[DIM * DIM];
        float* fx = NULL;

        if (sh->box_size != 0) {
            if (xdrfile_read_float(pvf, DIM * DIM, xd) == DIM * DIM) {
                if (NULL != box) {
                    for (int i = 0; (i < DIM); i++) {
                        for (int j = 0; (j < DIM); j++) {
                            box[i][j] = pvf[i * DIM + j];
                        }
                    }
                }
            } else {
                MD_LOG_ERROR("TRR: Failed to read frame box");
                return false;
            }
        }

        if (sh->vir_size != 0) {
            if (xdrfile_read_float(pvf, DIM * DIM, xd) != DIM * DIM) {
                return false;
            }
        }

        if (sh->pres_size != 0) {
            if (xdrfile_read_float(pvf, DIM * DIM, xd) != DIM * DIM) {
                return false;
            }
        }

        if (allocate_temp) {
            fx = (float*)calloc(sh->natoms * DIM, sizeof(fx[0]));
            if (NULL == fx) {
                return false;
            }
        }
        if (sh->x_size != 0) {
            if (x) {
                float c[3];
                for (int i = 0; i < sh->natoms; ++i) {
                    if (xdrfile_read_float(c, DIM, xd) == DIM) {
                        // @TODO: This scaling should be moved out of the core parts and into the loader which ties it with viamd.
                        x->x[i] = c[0] * 10;
                        x->y[i] = c[1] * 10;
                        x->z[i] = c[2] * 10;
                    } else {
                        MD_LOG_ERROR("TRR: Failed to read coordinate entry in frame");
                        return false;
                    }
                }
            } else {
                if (xdr_seek(xd, sh->x_size, SEEK_CUR) != exdrOK) {
                    MD_LOG_ERROR("TRR: Failed to skip coordinate section in frame");
                    return false;
                }
            }
        }
        if (sh->v_size != 0) {
            if (v) {
                float c[3];
                for (int i = 0; i < sh->natoms; ++i) {
                    if (xdrfile_read_float(c, DIM, xd) == DIM) {
                        // @TODO: This scaling should be moved out of the core parts and into the loader which ties it with viamd.
                        v->x[i] = c[0] * 10;
                        v->y[i] = c[1] * 10;
                        v->z[i] = c[2] * 10;
                    } else {
                        MD_LOG_ERROR("TRR: Failed to read velocity entry in frame");
                        return false;
                    }
                }
            } else {
                if (xdr_seek(xd, sh->v_size, SEEK_CUR) != exdrOK) {
                    MD_LOG_ERROR("TRR: Failed to skip velocity section in frame");
                    return false;
                }
            }
        }
        if (sh->f_size != 0) {
            if (f) {
                float c[3];
                for (int i = 0; i < sh->natoms; ++i) {
                    if (xdrfile_read_float(c, DIM, xd) == DIM) {
                        f->x[i] = c[0];
                        f->y[i] = c[1];
                        f->z[i] = c[2];
                    } else {
                        MD_LOG_ERROR("TRR: Failed to read force entry in frame");
                        return false;
                    }
                }
            } else {
                if (xdr_seek(xd, sh->f_size, SEEK_CUR) != exdrOK) {
                    MD_LOG_ERROR("TRR: Failed to skip force section in frame");
                    return false;
                }
            }
        }
        if (fx) {
            free(fx);
        }
    }
    return true;
}

static bool trr_read_frame_offsets_and_times(XDRFILE* xd, md_array(int64_t)* offsets, md_array(double)* times, md_allocator_i* alloc) {
    int result = 0;
    trr_header_t sh = {0};

    /* Go to file end */
    if (xdr_seek(xd, 0L, SEEK_END) != exdrOK) {
        return false;
    }

    int64_t filesize = xdr_tell(xd);

    /* Go to file beg */
    if (xdr_seek(xd, 0L, SEEK_SET) != exdrOK) {
        return false;
    }

    if (!trr_read_frame_header(xd, &sh)) {
        return false;
    }

    int framebytes = calc_framebytes(&sh);
    int est_nframes = (int)(filesize / ((int64_t)(framebytes + TRR_MIN_HEADER_SIZE)) + 1); /* must be at least 1 for successful growth */

    /* Allocate memory for the frame index array */
    md_array_ensure(*offsets, est_nframes, alloc);
    md_array_ensure(*times,   est_nframes, alloc);
    md_array_push(*offsets,    0, alloc);

    while (1) {
        /* Skip `framebytes` */
        result = xdr_seek(xd, (int64_t)(framebytes), SEEK_CUR);
        if (result != exdrOK) {
            break;
        }

        /* Store position in `offsets` */
        int64_t curr_offset = xdr_tell(xd);
        md_array_push(*offsets, curr_offset, alloc);
        md_array_push(*times,   sh.t, alloc);

        if (curr_offset == filesize) {
            return true;
        }

        /* Read header and calculate how much to skip next time */
        if (!trr_read_frame_header(xd, &sh)) {
            return false;
        }
        framebytes = calc_framebytes(&sh);
    }

    return false;
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
static int64_t trr_fetch_frame_data(struct md_trajectory_o* inst, int64_t frame_idx, void* frame_data_ptr) {
    trr_t* trr = (trr_t*)inst;
    ASSERT(trr);
    ASSERT(trr->magic == MD_TRR_TRAJ_MAGIC);

    if (!trr->file) {
        MD_LOG_ERROR("TRR: File handle is NULL");
        return 0;
    }

    if (!trr->frame_offsets) {
        MD_LOG_ERROR("TRR: Frame offsets is empty");
        return 0;
    }

    if (frame_idx < 0 || md_array_size(trr->frame_offsets) <= frame_idx + 1) {
        MD_LOG_ERROR("TRR: Frame index is out of range");
        return 0;
    }

    const int64_t beg = trr->frame_offsets[frame_idx];
    const int64_t end = trr->frame_offsets[frame_idx + 1];
    const int64_t frame_size = end - beg;

    if (frame_data_ptr) {
        ASSERT(trr->file);
        md_mutex_lock(&trr->mutex);
        // Seek and read must be an atomic operation to avoid race conditions
        // Since we use a shared file handle internally
        xdr_seek(trr->file, beg, SEEK_SET);
        const int64_t bytes_read = xdr_read(trr->file, frame_data_ptr, frame_size);
        md_mutex_unlock(&trr->mutex);
        ASSERT(frame_size == bytes_read);
    }
    return frame_size;
}

static bool trr_decode_frame_data(struct md_trajectory_o* inst, const void* frame_data_ptr, int64_t frame_data_size, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(inst);
    ASSERT(frame_data_ptr);
    ASSERT(frame_data_size);

    bool result = true;

    trr_t* trr = (trr_t*)inst;
    if (trr->magic != MD_TRR_TRAJ_MAGIC) {
        MD_LOG_ERROR("TRR: Error when decoding frame coord, trr magic did not match");
        return false;
    }

    if ((x || y || z) && !(x && y && z)) {
        MD_LOG_ERROR("TRR: User supplied coordinates (x,y,z) cannot be partially supplied");
        return false;
    }

    // There is a warning for ignoring const qualifier for frame_data_ptr, but it is ok since we only perform read operations "r" with the data.
    XDRFILE* file = xdrfile_mem((void*)frame_data_ptr, frame_data_size, "r");
    ASSERT(file);

    // Get header
    trr_header_t sh;
    mat3_t box = {0};
    md_vec3_soa_t coords = { x, y, z };
    result = trr_read_frame_header(file, &sh) && trr_read_frame_data(file, &sh, box.elem, &coords, 0, 0);
    if (result && header) {
        // @TODO: This scaling should be moved out of the core parts and into the loader which ties it with viamd.
        box = mat3_mul_f(box, 10.0f); // nm -> Ångström
        header->num_atoms = sh.natoms;
        header->index = sh.step;
        header->timestamp = sh.t;
        header->unit_cell = md_util_unit_cell_mat3(box);
    }

    xdrfile_close(file);

    return result;
}

bool trr_load_frame(struct md_trajectory_o* inst, int64_t frame_idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(inst);

    trr_t* trr = (trr_t*)inst;
    if (trr->magic != MD_TRR_TRAJ_MAGIC) {
        MD_LOG_ERROR("TRR: Error when decoding frame coord, trr magic did not match");
        return false;
    }

    // Should this be exposed?
    // This is a borderline case if one should use the md_temp_allocator as the raw frame size could potentially be several megabytes...
    md_allocator_i* alloc = md_heap_allocator;

    bool result = true;
    const int64_t frame_size = trr_fetch_frame_data(inst, frame_idx, NULL);
    if (frame_size > 0) {
        void* frame_data = md_alloc(alloc, frame_size);
        const int64_t read_size = trr_fetch_frame_data(inst, frame_idx, frame_data);
        ASSERT(read_size == frame_size);

        result = trr_decode_frame_data(inst, frame_data, frame_size, header, x, y, z);

        md_free(alloc, frame_data, frame_size);
    }

    return result;
}


static bool try_read_cache(md_array(int64_t)* frame_offsets, md_array(double)* frame_times, str_t cache_file, md_allocator_i* alloc) {
    bool result = false;
    md_file_o* file = md_file_open(cache_file, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        uint32_t magic = 0;
        uint32_t version = 0;
        int64_t  num_offsets = 0;
        int64_t  num_times = 0;

        if (md_file_read(file, &magic, sizeof(magic)) != sizeof(magic) || magic != MD_TRR_CACHE_MAGIC) {
            MD_LOG_DEBUG("TRR: Failed to read offset cache, magic was incorrect");
            goto done;
        }

        if (md_file_read(file, &version, sizeof(version)) != sizeof(version) || version != MD_TRR_CACHE_VERSION) {
            MD_LOG_DEBUG("TRR: Failed to read offset cache, version was incorrect");
            goto done;
        }

        if (md_file_read(file, &num_offsets, sizeof(num_offsets)) != sizeof(num_offsets) || num_offsets == 0) {
            MD_LOG_DEBUG("TRR: Failed to read offset cache, number of offsets was zero or corrupted");
            goto done;
        }

        md_array_resize(*frame_offsets, num_offsets, alloc);
        const int64_t offset_bytes = md_array_bytes(*frame_offsets);
        if (md_file_read(file, *frame_offsets, offset_bytes) != offset_bytes) {
            MD_LOG_DEBUG("TRR: Failed to read offsets from cache, incomplete");
            md_array_free(*frame_offsets, alloc);
            *frame_offsets = 0;
            goto done;
        }

        if (md_file_read(file, &num_times, sizeof(num_times)) != sizeof(num_times) || num_times == 0) {
            MD_LOG_DEBUG("TRR: Failed to read offset cache, number of times was zero or corrupted");
            goto done;
        }

        md_array_resize(*frame_times, num_times, alloc);
        const int64_t time_bytes = md_array_bytes(*frame_times);
        if (md_file_read(file, *frame_times, time_bytes) != time_bytes) {
            MD_LOG_DEBUG("TRR: Failed to read times from cache, incomplete");
            md_array_free(*frame_times, alloc);
            *frame_times = 0;
            goto done;
        }

        result = true;
    done:
        md_file_close(file);
    }

    return result;
}

static bool write_cache(str_t cache_file, const md_array(int64_t) frame_offsets, const md_array(double) frame_times) {
    bool result = false;

    if (md_array_size(frame_offsets) == 0 || md_array_size(frame_times) == 0 ||
        md_array_size(frame_offsets) != md_array_size(frame_times) + 1) {
        MD_LOG_DEBUG("XTC: cache offsets and times are not correct in length");
        return false;
    }
    
    md_file_o* file = md_file_open(cache_file, MD_FILE_WRITE | MD_FILE_BINARY);
    if (file) {
        const uint32_t magic = MD_TRR_CACHE_MAGIC;
        const uint32_t version = MD_TRR_CACHE_VERSION;
        const int64_t offset_bytes = md_array_bytes(frame_offsets);
        const int64_t num_offsets = md_array_size(frame_offsets);
        const int64_t time_bytes = md_array_bytes(frame_times);
        const int64_t num_times = md_array_size(frame_times);

        if (md_file_write(file, &magic,         sizeof(magic))          != sizeof(magic) ||
            md_file_write(file, &version,       sizeof(version))        != sizeof(version) ||
            md_file_write(file, &num_offsets,   sizeof(num_offsets))    != sizeof(num_offsets) ||
            md_file_write(file, frame_offsets,  offset_bytes)           != offset_bytes ||
            md_file_write(file, &num_times,     sizeof(num_times))      != sizeof(num_times) ||
            md_file_write(file, frame_times,    time_bytes)             != time_bytes)
        {
            MD_LOG_ERROR("TRR: Failed to write cache data");
            goto done;
        }
        result = true;

        done:
        md_file_close(file);
    } else {
        MD_LOG_ERROR("TRR: Failed to write cache, could not open file '%.*s", (int)cache_file.len, cache_file.ptr);
    }

    return result;
}

void trr_free(struct md_trajectory_o* inst) {
    trr_t* trr = (trr_t*)inst;
    if (trr->file) xdrfile_close(trr->file);
    if (trr->frame_offsets) md_array_free(trr->frame_offsets, trr->allocator);
    if (trr->header.frame_times) md_array_free(trr->header.frame_times, trr->allocator);
    md_mutex_destroy(&trr->mutex);
}

md_trajectory_i* md_trr_trajectory_create(str_t filename, md_allocator_i* alloc) {
    ASSERT(alloc);

    // Ensure that filename is zero terminated
    filename = str_copy(filename, md_temp_allocator);
    XDRFILE* file = xdrfile_open(filename.ptr, "r");

    if (file) {
        trr_header_t sh;
        if (!trr_read_frame_header(file, &sh)) {
            xdrfile_close(file);
            return false;
        }

        if (sh.natoms == 0) {
            MD_LOG_ERROR("TRR: Number of atoms is zero");
            xdrfile_close(file);
            return false;
        }

        char buf[2048];
        int len = snprintf(buf, sizeof(buf), "%.*s.cache", (int)filename.len, filename.ptr);
        str_t cache_file = {buf, len};

        md_array(int64_t) frame_offsets = 0;
        md_array(double)  frame_times   = 0;
        if (!try_read_cache(&frame_offsets, &frame_times, cache_file, alloc)) {
            md_array_shrink(frame_offsets, 0);
            md_array_shrink(frame_times, 0);

            if (!trr_read_frame_offsets_and_times(file, &frame_offsets, &frame_times, alloc)) {
                xdrfile_close(file);
                return false;
            }
            write_cache(cache_file, frame_offsets, frame_times);
        }

        if (!frame_offsets || !frame_times) {
            MD_LOG_ERROR("TRR: Failed to read frame offsets and times");
            xdrfile_close(file);
            return false;
        }

        const int64_t num_frames = md_array_size(frame_offsets) - 1;

        int64_t max_frame_size = 0;
        for (int64_t i = 0; i < num_frames; ++i) {
            const int64_t beg = frame_offsets[i + 0];
            const int64_t end = frame_offsets[i + 1];
            const int64_t frame_size = end - beg;
            max_frame_size = MAX(max_frame_size, frame_size);
        }

        void* mem = md_alloc(alloc, sizeof(md_trajectory_i) + sizeof(trr_t));
        ASSERT(mem);
        MEMSET(mem, 0, sizeof(md_trajectory_i) + sizeof(trr_t));

        md_trajectory_i* traj = mem;
        trr_t* trr = (trr_t*)(traj + 1);

        trr->magic = MD_TRR_TRAJ_MAGIC;
        trr->allocator = alloc;
        trr->file = file;
        trr->frame_offsets = frame_offsets;
        trr->mutex = md_mutex_create();

        trr->header = (md_trajectory_header_t) {
            .num_frames = num_frames,
            .num_atoms = sh.natoms,
            .max_frame_data_size = max_frame_size,
            .time_unit = md_unit_pikosecond(),
            .frame_times = frame_times,
        };

        traj->inst = (struct md_trajectory_o*)trr;
        traj->get_header = trr_get_header;
        traj->load_frame = trr_load_frame;
        traj->fetch_frame_data = trr_fetch_frame_data;
        traj->decode_frame_data = trr_decode_frame_data;

        return traj;
    }
    
    return NULL;
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

    md_allocator_i* alloc = trr->allocator;
    trr_free(traj->inst);
    md_free(alloc, traj, sizeof(md_trajectory_i) + sizeof(trr_t));
}

static md_trajectory_loader_i trr_loader = {
    md_trr_trajectory_create,
    md_trr_trajectory_free,
};

md_trajectory_loader_i* md_trr_trajectory_loader(void) {
    return &trr_loader;
}
