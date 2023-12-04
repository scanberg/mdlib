#include <md_xtc.h>

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

#define MD_XTC_CACHE_MAGIC   0x8281237
#define MD_XTC_CACHE_VERSION 2

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

static bool xtc_frame_header(XDRFILE* xd, int* natoms, int* step, float* time, matrix box) {
    int magic;

    if ((xdrfile_read_int(&magic, 1, xd)) != 1) {
        MD_LOG_ERROR("XTC: Failed to read magic number in header");
        return false;
    }
    if (magic != XTC_MAGIC) {
        MD_LOG_ERROR("XTC: Magic number did not match");
        return false;
    }
    if ((xdrfile_read_int(natoms, 1, xd)) != 1) {
        MD_LOG_ERROR("XTC: Failed to read number of atoms");
        return false;
    }
    if ((xdrfile_read_int(step, 1, xd)) != 1) {
        MD_LOG_ERROR("XTC: Failed to read step");
        return false;
    }
    if ((xdrfile_read_float(time, 1, xd)) != 1) {
        MD_LOG_ERROR("XTC: Failed to read timestamp");
        return false;
    }
    if ((xdrfile_read_float(box[0], DIM * DIM, xd)) != DIM * DIM) {
        MD_LOG_ERROR("XTC: Failed to read box dimensions");
        return false;
    }
    /*
    if ((result = xdrfile_read_int(&natoms2, 1, xd)) != 1) {
        MD_LOG_ERROR("Failed to read number of atoms (second)");
        return false;
    }
    if (*natoms != natoms2) {
        MD_LOG_ERROR("Number of atoms did not match second number of atoms, header might be corrupt");
        return false;
    }
    */

    return true;
}

static bool xtc_frame_coords(XDRFILE* xd, int natoms, rvec* x) {
    int result;
    float prec;

    result = xdrfile_decompress_coord_float(x[0], &natoms, &prec, xd);
    if (result != natoms) {
        MD_LOG_ERROR("XTC: Failed to read coordinates");
        return false;
    }

    return true;
}

static bool xtc_frame_offsets_and_times(XDRFILE* xd, md_array(int64_t)* frame_offsets, md_array(double)* frame_times, md_allocator_i* alloc) {
    int step, natoms;
    float time;
    float box[3][3];

    /* Go to file beg */
    if (xdr_seek(xd, 0L, SEEK_SET) != exdrOK) {
        MD_LOG_ERROR("XTC: Failed to seek to beginning of file");
        return false;
    }

    if (!xtc_frame_header(xd, &natoms, &step, &time, box)) {
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
        const int64_t framebytes = XTC_SMALL_HEADER_SIZE + XTC_SMALL_COORDS_SIZE * natoms;
        const int64_t num_frames = (filesize / framebytes); /* Should we complain if framesize doesn't divide filesize? */

        md_array_push(*frame_offsets, 0, alloc);
        md_array_push(*frame_times, time, alloc);

        for (int64_t i = 1; i < num_frames; i++) {
            const int64_t offset = i * framebytes;

            if (xdr_seek(xd, offset, SEEK_SET) != exdrOK || !xtc_frame_header(xd, &natoms, &step, &time, box)) {
                md_array_free(*frame_offsets, alloc);
                md_array_free(*frame_times, alloc);
                *frame_offsets = 0;
                *frame_times = 0;
                return false;
            }

            md_array_push(*frame_offsets, i * framebytes, alloc);
            md_array_push(*frame_times, time, alloc);
        }
        md_array_push(*frame_offsets, num_frames * framebytes, alloc);

        return true;
    } else {
        int framebytes, est_nframes;

        /* Move pos back to end of first header */
        if (xdr_seek(xd, (int64_t)XTC_HEADER_SIZE, SEEK_SET) != exdrOK) {
            return false;
        }
        
        if (xdrfile_read_int(&framebytes, 1, xd) == 0) {
            MD_LOG_ERROR("XTC: Failed to read framebytes");
            return false;
        }
        framebytes = (framebytes + 3) & ~0x03; /* Rounding to the next 32-bit boundary */
        est_nframes = (int)(filesize / ((int64_t)(framebytes + XTC_HEADER_SIZE)) +
            1); /* must be at least 1 for successful growth */
                /* First `framebytes` might be larger than average, so we would underestimate `est_nframes`
                */
        est_nframes += est_nframes / 5;

        /* Skip `framebytes` */
        if (xdr_seek(xd, (int64_t)(framebytes), SEEK_CUR) != exdrOK) {
            goto fail;
        }
        
        md_array_ensure(*frame_offsets, est_nframes, alloc);
        md_array_ensure(*frame_times,   est_nframes, alloc);

        md_array_push(*frame_offsets, 0, alloc);
        md_array_push(*frame_times, time, alloc);

        while (1) {
            const int64_t offset = xdr_tell(xd);
            if (offset == filesize) {
                // Done, add final filesize as offset
                md_array_push(*frame_offsets, offset, alloc);
                return true;
            }
            if (!xtc_frame_header(xd, &natoms, &step, &time, box)) {
                goto fail;
            }

            /* Store position in `offsets`, adjust for header */
            md_array_push(*frame_offsets, offset, alloc);
            md_array_push(*frame_times, time, alloc);

            if (xdr_seek(xd, offset + XTC_HEADER_SIZE, SEEK_SET) != exdrOK) {
                goto fail;
            }
            /* Read how much to skip */
            if (xdrfile_read_int(&framebytes, 1, xd) == 0) {
                goto fail;
            }
            framebytes = (framebytes + 3) & ~0x03; /* Rounding to the next 32-bit boundary */

            /* Skip `framebytes` to next header */
            if (xdr_seek(xd, framebytes, SEEK_CUR) != exdrOK) {
                goto fail;
            }
        }

    fail:
        md_array_free(*frame_offsets, alloc);
        md_array_free(*frame_times, alloc);
        *frame_offsets = 0;
        *frame_times = 0;
        return false;
    }
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
        MD_LOG_ERROR("XTC: File handle is NULL");
        return 0;
    }

    if (!xtc->frame_offsets) {
        MD_LOG_ERROR("XTC: Frame offsets is empty");
        return 0;
    }

    if (frame_idx < 0 || md_array_size(xtc->frame_offsets) <= frame_idx + 1) {
        MD_LOG_ERROR("XTC: Frame index is out of range");
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
        (void)bytes_read;
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
        MD_LOG_ERROR("XTC: Error when decoding frame coord, xtc magic did not match");
        return false;
    }

    if ((x || y || z) && !(x && y && z)) {
        MD_LOG_ERROR("XTC: User supplied coordinates (x,y,z) cannot be partially supplied");
        return false;
    }

    // There is a warning for ignoring const qualifier for frame_data_ptr, but it is ok since we only perform read operations "r" with the data.
    XDRFILE* file = xdrfile_mem((void*)frame_data_ptr, frame_data_size, "r");
    ASSERT(file);

    // Get header
    int natoms = 0, step = 0;
    float time = 0;
    float box[3][3];
    result = xtc_frame_header(file, &natoms, &step, &time, box);
    if (result) {
        if (header) {
            // nm -> Ångström
            for (int i = 0; i < 3; ++i) {
                box[i][0] *= 10.0f;
                box[i][1] *= 10.0f;
                box[i][2] *= 10.0f;
            }
            header->num_atoms = natoms;
            header->index = step;
            header->timestamp = time;
            header->unit_cell = md_util_unit_cell_from_matrix(box);
        }

        if (x || y || z) {
            int64_t byte_size = natoms * sizeof(rvec);
            rvec* pos = md_alloc(md_heap_allocator, byte_size);
            result = xtc_frame_coords(file, natoms, pos);
            if (result) {            
                // nm -> Ångström
                for (int64_t i = 0; i < natoms; ++i) {
                    if (x) x[i] = pos[i][0] * 10.0f;
                    if (y) y[i] = pos[i][1] * 10.0f;
                    if (z) z[i] = pos[i][2] * 10.0f;
                }
            }
            md_free(md_heap_allocator, pos, byte_size);
        }
    }

    xdrfile_close(file);

    return result;
}

bool xtc_load_frame(struct md_trajectory_o* inst, int64_t frame_idx, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
    ASSERT(inst);

    xtc_t* xtc = (xtc_t*)inst;
    if (xtc->magic != MD_XTC_TRAJ_MAGIC) {
        MD_LOG_ERROR("XTC: Error when decoding frame coord, xtc magic did not match");
        return false;
    }

    // Should this be exposed?
    md_allocator_i* alloc = md_temp_allocator;

    bool result = true;
    const int64_t frame_size = xtc_fetch_frame_data(inst, frame_idx, NULL);
    if (frame_size > 0) {
        // This is a borderline case if one should use the md_temp_allocator as the raw frame size could potentially be several megabytes...
        void* frame_data = md_alloc(alloc, frame_size);
        const int64_t read_size = xtc_fetch_frame_data(inst, frame_idx, frame_data);
        ASSERT(read_size == frame_size);

        result = xtc_decode_frame_data(inst, frame_data, frame_size, header, x, y, z);

        md_free(alloc, frame_data, frame_size);
    }

    return result;
}

static bool try_read_cache(md_array(int64_t)* frame_offsets, md_array(double)* frame_times, str_t cache_file, md_allocator_i* alloc) {
    bool result = false;
    md_file_o* file = md_file_open(cache_file, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        int magic = 0;
        int version = 0;
        int64_t num_offsets = 0;
        int64_t num_times = 0;

        if (md_file_read(file, &magic, sizeof(magic)) != sizeof(magic)) {
            MD_LOG_ERROR("XTC: Failed to read offset cache magic");
            goto done;
        }

        if (magic != MD_XTC_CACHE_MAGIC) {
            MD_LOG_ERROR("XTC: Incorrect cache magic");
            goto done;
        }

        if (md_file_read(file, &version, sizeof(version)) != sizeof(version)) {
            MD_LOG_ERROR("XTC: Failed to read cache version");
            goto done;
        }

        if (version != MD_XTC_CACHE_VERSION) {
            MD_LOG_ERROR("XTC: Cache version is incorrect");
            goto done;
        }

        if (md_file_read(file, &num_offsets, sizeof(num_offsets)) != sizeof(num_offsets) || num_offsets == 0) {
            MD_LOG_ERROR("XTC: Failed to read offset cache number of frames");
            goto done;
        }

        md_array_resize(*frame_offsets, num_offsets, alloc);
        const int64_t offset_bytes = md_array_bytes(*frame_offsets);
        if (md_file_read(file, *frame_offsets, offset_bytes) != offset_bytes) {
            MD_LOG_ERROR("XTC: Failed to read offset cache, offsets are incomplete");
            goto done;
        }

        if (md_file_read(file, &num_times, sizeof(num_times)) != sizeof(num_times) || num_times == 0) {
            MD_LOG_ERROR("XTC: Failed to read offset cache number of times");
            goto done;
        }

        md_array_resize(*frame_times, num_times, alloc);
        const int64_t time_bytes = md_array_bytes(*frame_times);
        if (md_file_read(file, *frame_times, time_bytes) != time_bytes) {
            MD_LOG_ERROR("XTC: Failed to read cache, times are incomplete");
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

    const int64_t num_offsets = md_array_size(frame_offsets);
    const int64_t num_times = md_array_size(frame_times);

    if (num_offsets == 0 || num_times == 0 || num_offsets != num_times + 1) {
        MD_LOG_DEBUG("XTC: cache offsets and times are not correct in length");
        return false;
    }

    md_file_o* file = md_file_open(cache_file, MD_FILE_WRITE | MD_FILE_BINARY);
    if (file) {
        const int magic   = MD_XTC_CACHE_MAGIC;
        const int version = MD_XTC_CACHE_VERSION;
        const int64_t offset_bytes = md_array_bytes(frame_offsets);
        const int64_t time_bytes = md_array_bytes(frame_times);
        
        if (md_file_write(file, &magic,         sizeof(magic))          != sizeof(magic) ||
            md_file_write(file, &version,       sizeof(version))        != sizeof(version) ||
            md_file_write(file, &num_offsets,   sizeof(num_offsets))    != sizeof(num_offsets) ||
            md_file_write(file, frame_offsets,  offset_bytes)           != offset_bytes ||
            md_file_write(file, &num_times,     sizeof(num_times))      != sizeof(num_times) ||
            md_file_write(file, frame_times,    time_bytes)             != time_bytes)
        {
            MD_LOG_ERROR("XTC: Failed to write cache data");
            goto done;
        }
        result = true;

    done:
        md_file_close(file);
    } else {
        MD_LOG_ERROR("XTC: Failed to write cache, could not open file '%.*s", (int)cache_file.len, cache_file.ptr);
    }

    return result;
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
    filename = str_copy(filename, md_temp_allocator);
    XDRFILE* file = xdrfile_open(filename.ptr, "r");

    if (file) {
        int num_atoms, step;
        float time;
        float box[3][3];
        if (!xtc_frame_header(file, &num_atoms, &step, &time, box)) {
            xdrfile_close(file);
            return false;
        }

        if (num_atoms == 0) {
            MD_LOG_ERROR("XTC: Number of atoms in trajectory was zero");
            xdrfile_close(file);
            return false;
        }

        char buf[1024];
        int len = snprintf(buf, sizeof(buf), "%.*s.cache", (int)filename.len, filename.ptr);
        str_t cache_file = {buf, len};

        md_array(int64_t) frame_offsets = 0;
        md_array(double)  frame_times   = 0;
        if (!try_read_cache(&frame_offsets, &frame_times, cache_file, alloc)) {
            md_array_shrink(frame_offsets, 0);
            md_array_shrink(frame_times, 0);

            MD_LOG_INFO("XTC: Attempting to create new cache file...");
            if (!xtc_frame_offsets_and_times(file, &frame_offsets, &frame_times, alloc)) {
                xdrfile_close(file);
                return false;
            }
            if (write_cache(cache_file, frame_offsets, frame_times)) {
                MD_LOG_INFO("XTC: Successfully created cache file for trajectory");
            }
        }

        if (!frame_offsets || !frame_times) {
            MD_LOG_DEBUG("XTC: frame offsets or frame times was empty");
            xdrfile_close(file);
            return false;
        }

        int64_t max_frame_size = 0;
        for (int64_t i = 0; i < md_array_size(frame_offsets) - 1; ++i) {
            const int64_t beg = frame_offsets[i + 0];
            const int64_t end = frame_offsets[i + 1];
            const int64_t frame_size = end - beg;
            max_frame_size = MAX(max_frame_size, frame_size);
        }
        
        const int64_t num_frames = md_array_size(frame_offsets) - 1;

        void* mem = md_alloc(alloc, sizeof(md_trajectory_i) + sizeof(xtc_t));
        ASSERT(mem);
        MEMSET(mem, 0, sizeof(md_trajectory_i) + sizeof(xtc_t));

        md_trajectory_i* traj = mem;
        xtc_t* xtc = (xtc_t*)(traj + 1);

        xtc->magic = MD_XTC_TRAJ_MAGIC;
        xtc->allocator = alloc;
        xtc->file = file;
        xtc->frame_offsets = frame_offsets;
        xtc->mutex = md_mutex_create();

        xtc->header = (md_trajectory_header_t) {
            .num_frames = num_frames,
            .num_atoms = num_atoms,
            .max_frame_data_size = max_frame_size,
            .time_unit = md_unit_pikosecond(),
            .frame_times = frame_times,
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
        MD_LOG_ERROR("XTC: Cannot free trajectory, is not a valid XTC trajectory.");
        ASSERT(false);
        return;
    }

    md_allocator_i* alloc = xtc->allocator;
    xtc_free(traj->inst);
    md_free(alloc, traj, sizeof(md_trajectory_i) + sizeof(xtc_t));
}

static md_trajectory_loader_i xtc_loader = {
    md_xtc_trajectory_create,
    md_xtc_trajectory_free,
};

md_trajectory_loader_i* md_xtc_trajectory_loader(void) {
    return &xtc_loader;
}
