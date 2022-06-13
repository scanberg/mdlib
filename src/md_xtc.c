#include "md_xtc.h"

#include <core/md_common.h>
#include <core/md_array.inl>
#include <core/md_file.h>
#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_sync.h>
#include <md_trajectory.h>

#include <ext/xtc/xdrfile.h>

#include <string.h>
#include <stdio.h>

#define MD_XTC_TRAJ_MAGIC 0x162365dac721995

#define XTC_MAGIC 1995

/* XTC small header size (natoms<=9).
*  > int(4) magic
*  > int(4) natoms
*  > int(4) step
*  > float(4) time
*  > 9xfloat(4) box
*  > int(4) natoms (again)
*/
#define XTC_SMALL_HEADER_SIZE 56

/* Size of uncompressed coordinates for one atom.
* 3xfloat(4) x
*/
#define XTC_SMALL_COORDS_SIZE 12

/* XTC header size (natoms>=10).
* Compressed trajectories contain some additional values:
*  > float(4) precision
*  > 3xint(4) minint
*  > 3xint(4) maxint
*  > int(4) smallidx
* See `xdrfile_compress_coord_double()`.
*/
#define XTC_HEADER_SIZE (XTC_SMALL_HEADER_SIZE + 32)

typedef struct xtc_t {
    uint64_t magic;
    XDRFILE* file;
    int64_t* frame_offsets;
    md_trajectory_header_t header;
    md_allocator_i* allocator;
    md_mutex_t mutex;
} xtc_t;

static bool xtc_header(XDRFILE* xd, int* natoms, int* step, float* time, matrix box) {
    int result, magic;

    if ((result = xdrfile_read_int(&magic, 1, xd)) != 1) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to read magic number");
        return false;
    }
    if (magic != XTC_MAGIC) {
        md_print(MD_LOG_TYPE_ERROR, "Magic number did not match");
        return false;
    }
    if ((result = xdrfile_read_int(natoms, 1, xd)) != 1) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to read number of atoms");
        return false;
    }
    if ((result = xdrfile_read_int(step, 1, xd)) != 1) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to read step");
        return false;
    }
    if ((result = xdrfile_read_float(time, 1, xd)) != 1) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to read timestamp");
        return false;
    }
    if ((result = xdrfile_read_float(box[0], DIM * DIM, xd)) != DIM * DIM) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to read box dimensions");
        return false;
    }
    /*
    if ((result = xdrfile_read_int(&natoms2, 1, xd)) != 1) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to read number of atoms (second)");
        return false;
    }
    if (*natoms != natoms2) {
        md_print(MD_LOG_TYPE_ERROR, "Number of atoms did not match second number of atoms, header might be corrupt");
        return false;
    }
    */

    return true;
}

static bool xtc_coord(XDRFILE* xd, int natoms, rvec* pos) {
    int result;
    float prec;

    result = xdrfile_decompress_coord_float(pos[0], &natoms, &prec, xd);
    if (result != natoms) {
        md_print(MD_LOG_TYPE_ERROR, "Failed to read coordinates");
        return false;
    }

    return true;
}

static bool xtc_offsets(XDRFILE* xd, int64_t** offsets, md_allocator_i* alloc) {
    int i, result, est_nframes, framebytes, step, natoms, nframes;
    float time;
    float box[3][3];

    int64_t* arr = 0;

    /* Go to file beg */
    if (xdr_seek(xd, 0L, SEEK_SET) != exdrOK) {
        return false;
    }

    if (!xtc_header(xd, &natoms, &step, &time, box)) {
        return false;
    }

    /* Go to file end */
    if (xdr_seek(xd, 0L, SEEK_END) != exdrOK) {
        return false;
    }
    /* Cursor position is equivalent to file size */
    int64_t filesize = xdr_tell(xd);

    /* Dont bother with compression for nine atoms or less */
    if (natoms <= 9) {
        xdrfile_close(xd);
        framebytes = XTC_SMALL_HEADER_SIZE + XTC_SMALL_COORDS_SIZE * natoms;
        nframes = (int)(filesize / framebytes); /* Should we complain if framesize doesn't divide filesize? */
        for (i = 0; i < nframes; i++) {
            md_array_push(arr, (int64_t)i * (int64_t)framebytes, alloc);
        }
    } else {
        /* Go back to the beginning of the file */
        if (xdr_seek(xd, (int64_t)XTC_HEADER_SIZE, SEEK_SET) != exdrOK) {
            return false;
        }

        if (xdrfile_read_int(&framebytes, 1, xd) == 0) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to read framebytes");
            return false;
        }
        framebytes = (framebytes + 3) & ~0x03; /* Rounding to the next 32-bit boundary */
        est_nframes = (int)(filesize / ((int64_t)(framebytes + XTC_HEADER_SIZE)) +
            1); /* must be at least 1 for successful growth */
                /* First `framebytes` might be larger than average, so we would underestimate `est_nframes`
                */
        est_nframes += est_nframes / 5;
        md_array_ensure(arr, est_nframes, alloc);
        md_array_push(arr, 0, alloc);

        while (1) {
            /* Skip `framebytes` and next header */
            result = xdr_seek(xd, (int64_t)(framebytes + XTC_HEADER_SIZE), SEEK_CUR);
            if (result != exdrOK) {
                md_print(MD_LOG_TYPE_ERROR, "Failed to skip framebytes");
                md_array_free(arr, alloc);
                return false;
            }

            /* Store position in `offsets`, adjust for header */
            md_array_push(arr, xdr_tell(xd) - (int64_t)(XTC_HEADER_SIZE), alloc);

            /* Read how much to skip next time */
            if (xdrfile_read_int(&framebytes, 1, xd) == 0) {
                *offsets = arr;
                return true;
            }
            framebytes = (framebytes + 3) & ~0x03; /* Rounding to the next 32-bit boundary */
        }
    }

    return false;
}

bool xtc_get_header(struct md_trajectory_o* inst, md_trajectory_header_t* header) {
    xtc_t* xtc = (xtc_t*)inst;
    ASSERT(xtc);
    ASSERT(xtc->magic == MD_XTC_TRAJ_MAGIC);
    ASSERT(header);

    *header = xtc->header;
    return true;
}

// This is lowlevel cruft for enabling parallel loading and decoding of frames
static int64_t xtc_fetch_frame_data(struct md_trajectory_o* inst, int64_t frame_idx, void* frame_data_ptr) {
    xtc_t* xtc = (xtc_t*)inst;
    ASSERT(xtc);
    ASSERT(xtc->magic == MD_XTC_TRAJ_MAGIC);

    if (!xtc->file) {
        md_print(MD_LOG_TYPE_ERROR, "File handle is NULL");
        return 0;
    }

    if (!xtc->frame_offsets) {
        md_print(MD_LOG_TYPE_ERROR, "Frame offsets is empty");
        return 0;
    }

    if (frame_idx < 0 || md_array_size(xtc->frame_offsets) <= frame_idx + 1) {
        md_print(MD_LOG_TYPE_ERROR, "Frame index is out of range");
        return 0;
    }

    const int64_t beg = xtc->frame_offsets[frame_idx];
    const int64_t end = xtc->frame_offsets[frame_idx + 1];
    const int64_t frame_size = end - beg;

    if (frame_data_ptr) {
        ASSERT(xtc->file);
        md_mutex_lock(&xtc->mutex);
        // Seek and read must be an atomic operation to avoid race conditions
        xdr_seek(xtc->file, beg, SEEK_SET);
        const int64_t bytes_read = xdr_read(xtc->file, frame_data_ptr, frame_size);
        md_mutex_unlock(&xtc->mutex);
        ASSERT(frame_size == bytes_read);
    }
    return frame_size;
}

static bool xtc_decode_frame_data(struct md_trajectory_o* inst, const void* frame_data_ptr, int64_t frame_data_size, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(inst);
    ASSERT(frame_data_ptr);
    ASSERT(frame_data_size);

    bool result = true;

    xtc_t* xtc = (xtc_t*)inst;
    if (xtc->magic != MD_XTC_TRAJ_MAGIC) {
        md_print(MD_LOG_TYPE_ERROR, "Error when decoding frame coord, xtc magic did not match");
        return false;
    }

    // There is a warning for ignoring const qualifier for frame_data_ptr, but it is ok since we only perform read operations "r" with the data.
    XDRFILE* file = xdrfile_mem((void*)frame_data_ptr, frame_data_size, "r");
    ASSERT(file);

    // Get header
    int natoms = 0, step = 0;
    float time = 0;
    float box[3][3] = {0};
    result = xtc_header(file, &natoms, &step, &time, box);
    if (result) {
        if (header) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    box[i][j] *= 10.0f;
                }
            }
            header->num_atoms = natoms;
            header->index = step;
            header->timestamp = time;
            memcpy(header->box, box, sizeof(header->box));
        }

        if (x || y || z) {
            int64_t byte_size = natoms * sizeof(rvec);
            rvec* pos = md_alloc(default_allocator, byte_size);
            result = xtc_coord(file, natoms, pos);
            if (result) {            
                // nm -> Ångström
                for (int64_t i = 0; i < natoms; ++i) {
                    if (x) x[i] = pos[i][0] * 10.0f;
                    if (y) y[i] = pos[i][1] * 10.0f;
                    if (z) z[i] = pos[i][2] * 10.0f;
                }
            }
            md_free(default_allocator, pos, byte_size);
        }
    }

    xdrfile_close(file);

    return result;
}

bool xtc_load_frame(struct md_trajectory_o* inst, int64_t frame_idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(inst);

    xtc_t* xtc = (xtc_t*)inst;
    if (xtc->magic != MD_XTC_TRAJ_MAGIC) {
        md_print(MD_LOG_TYPE_ERROR, "Error when decoding frame coord, xtc magic did not match");
        return false;
    }

    // Should this be exposed?
    md_allocator_i* alloc = default_temp_allocator;

    bool result = true;
    const int64_t frame_size = xtc_fetch_frame_data(inst, frame_idx, NULL);
    if (frame_size > 0) {
        // This is a borderline case if one should use the default_temp_allocator as the raw frame size could potentially be several megabytes...
        void* frame_data = md_alloc(alloc, frame_size);
        const int64_t read_size = xtc_fetch_frame_data(inst, frame_idx, frame_data);
        ASSERT(read_size == frame_size);

        result = xtc_decode_frame_data(inst, frame_data, frame_size, header, x, y, z);

        md_free(alloc, frame_data, frame_size);
    }

    return result;
}

#define OFFSET_MAGIC 0x828123727bcbabc
static int64_t* try_read_offset_cache(str_t cache_file, md_allocator_i* alloc) {
    md_file_o* file = md_file_open(cache_file, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        int64_t* offsets = 0;
        int64_t read_bytes = 0;
        int64_t num_frames = 0;
        uint64_t magic = 0;

        read_bytes = md_file_read(file, &magic, sizeof(magic));
        if (read_bytes != sizeof(magic) || magic != OFFSET_MAGIC) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to read offset cache, magic was incorrect");
            goto done;
        }

        read_bytes = md_file_read(file, &num_frames, sizeof(num_frames));
        if (read_bytes != sizeof(num_frames) || num_frames == 0) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to read offset cache, number of frames was zero or corrupted");
            goto done;
        }

        md_array_resize(offsets, num_frames, alloc);

        read_bytes = md_file_read(file, offsets, num_frames * sizeof(int64_t));
        if (read_bytes != (int64_t)(num_frames * sizeof(int64_t))) {
            md_print(MD_LOG_TYPE_ERROR, "Failed to read offset cache, offsets are incomplete");
            md_array_free(offsets, alloc);
            offsets = 0;
            goto done;
        }
    done:
        md_file_close(file);
        return offsets;
    }

    return 0;
}

static bool write_offset_cache(str_t cache_file, int64_t* offsets) {
    md_file_o* file = md_file_open(cache_file, MD_FILE_WRITE | MD_FILE_BINARY);
    if (file) {
        int64_t num_offsets = md_array_size(offsets);
        int64_t written_bytes = 0;
        const uint64_t magic = OFFSET_MAGIC;
        written_bytes = md_file_write(file, &magic, sizeof(magic));
        ASSERT(written_bytes == sizeof(OFFSET_MAGIC));

        written_bytes = md_file_write(file, &num_offsets, sizeof(num_offsets));
        ASSERT(written_bytes == sizeof(num_offsets));

        written_bytes = md_file_write(file, offsets, num_offsets * sizeof(int64_t));
        ASSERT(written_bytes == (int64_t)(num_offsets * sizeof(int64_t)));

        md_file_close(file);
        return true;
    }

    md_printf(MD_LOG_TYPE_ERROR, "Failed to write offset cache, could not open file '%.*s", (int)cache_file.len, cache_file.ptr);
    return false;
}

void xtc_free(struct md_trajectory_o* inst) {
    xtc_t* xtc = (xtc_t*)inst;
    if (xtc->file) xdrfile_close(xtc->file);
    if (xtc->frame_offsets) md_array_free(xtc->frame_offsets, xtc->allocator);
    md_mutex_destroy(&xtc->mutex);
}

md_trajectory_i* md_xtc_trajectory_create(str_t filename, md_allocator_i* alloc) {
    ASSERT(alloc);

    // Ensure that filename is zero terminated
    filename = copy_str(filename, default_temp_allocator);
    XDRFILE* file = xdrfile_open(filename.ptr, "r");

    if (file) {
        int num_atoms, step;
        float time;
        float box[3][3];
        if (!xtc_header(file, &num_atoms, &step, &time, box)) {
            xdrfile_close(file);
            return false;
        }

        char buf[1024];
        int len = snprintf(buf, sizeof(buf), "%.*s.cache", (int)filename.len, filename.ptr);
        str_t cache_file = {buf, len};

        int64_t* offsets = try_read_offset_cache(cache_file, alloc);
        if (!offsets) {
            if (!xtc_offsets(file, &offsets, alloc)) {
                xdrfile_close(file);
                return false;
            }
            write_offset_cache(cache_file, offsets);
        }

        int64_t max_frame_size = 0;
        for (int64_t i = 0; i < md_array_size(offsets) - 1; ++i) {
            const int64_t beg = offsets[i + 0];
            const int64_t end = offsets[i + 1];
            const int64_t frame_size = end - beg;
            max_frame_size = MAX(max_frame_size, frame_size);
        }
        
        ASSERT(offsets);
        ASSERT(num_atoms > 0);

        void* mem = md_alloc(alloc, sizeof(md_trajectory_i) + sizeof(xtc_t));
        ASSERT(mem);
        memset(mem, 0, sizeof(md_trajectory_i) + sizeof(xtc_t));

        md_trajectory_i* traj = mem;
        xtc_t* xtc = (xtc_t*)(traj + 1);

        xtc->magic = MD_XTC_TRAJ_MAGIC;
        xtc->allocator = alloc;
        xtc->file = file;
        xtc->frame_offsets = offsets;
        xtc->mutex = md_mutex_create();

        xtc->header = (md_trajectory_header_t) {
            .num_frames = md_array_size(offsets) - 1,
            .num_atoms = num_atoms,
            .max_frame_data_size = max_frame_size,
            .time_unit = {.base = {.time = UNIT_TIME_PIKOSECONDS}, .dim = {.time = 1}},
        };

        traj->inst = (struct md_trajectory_o*)xtc;
        traj->get_header = xtc_get_header;
        traj->load_frame = xtc_load_frame;
        traj->fetch_frame_data = xtc_fetch_frame_data;
        traj->decode_frame_data = xtc_decode_frame_data;

        return traj;
    }
    
    return NULL;
}

void md_xtc_trajectory_free(md_trajectory_i* traj) {
    ASSERT(traj);
    ASSERT(traj->inst);
    xtc_t* xtc = (xtc_t*)traj->inst;
    if (xtc->magic != MD_XTC_TRAJ_MAGIC) {
        md_printf(MD_LOG_TYPE_ERROR, "Trajectory is not a valid XTC trajectory.");
        ASSERT(false);
        return;
    }

    md_allocator_i* alloc = xtc->allocator;
    xtc_free(traj->inst);
    md_free(alloc, traj, sizeof(md_trajectory_i) + sizeof(xtc_t));
}

static md_trajectory_api xtc_api = {
    md_xtc_trajectory_create,
    md_xtc_trajectory_free,
};

md_trajectory_api* md_xtc_trajectory_api(void) {
    return &xtc_api;
}