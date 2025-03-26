#include <md_xtc.h>

#include <md_util.h>
#include <md_trajectory.h>

#include <core/md_common.h>
#include <core/md_array.h>
#include <core/md_str_builder.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_vec_math.h>

#include <xdrfile.h>

#include <string.h>
#include <stdio.h>

#define MD_XTC_CACHE_MAGIC   0x8281237612371
#define MD_XTC_CACHE_VERSION 3

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
    md_array(int64_t) frame_offsets;
    md_trajectory_header_t header;
    md_allocator_i* allocator;
    md_mutex_t mutex;    
} xtc_t;


// XDR Specific stuff
static inline bool xdr_read_bytes (uint8_t* ptr, size_t count, md_file_o* xdr_file) {
    ASSERT(ptr);
    ASSERT(xdr_file);

    size_t bytes = count;
    if (md_file_read(xdr_file, ptr, bytes) == bytes) {
        return true;
    }

    return false;
}

static inline bool xdr_read_int16(int16_t* ptr, size_t count, md_file_o* xdr_file) {
    ASSERT(ptr);
    ASSERT(xdr_file);

    size_t bytes = count * sizeof(int16_t);
    if (md_file_read(xdr_file, ptr, bytes) == bytes) {
#if __LITTLE_ENDIAN__
        for (size_t i = 0; i < count; ++i) {
            ptr[i] = BSWAP16(ptr[i]);
        }
#endif
        return true;
    }

    return false;
}

static inline bool xdr_read_int32(int32_t* ptr, size_t count, md_file_o* xdr_file) {
    ASSERT(ptr);
    ASSERT(xdr_file);

    size_t bytes = count * sizeof(int32_t);
    if (md_file_read(xdr_file, ptr, bytes) == bytes) {
#if __LITTLE_ENDIAN__
        for (size_t i = 0; i < count; ++i) {
            ptr[i] = BSWAP32(ptr[i]);
        }
#endif
        return true;
    }

    return false;
}

static inline bool xdr_read_float (float* ptr, size_t count, md_file_o* xdr_file) {
    ASSERT(ptr);
    ASSERT(xdr_file);
    return xdr_read_int32((int32_t*)ptr, count, xdr_file);
}



static inline int sizeofint(int size) {
    unsigned int num = 1;
    int num_of_bits = 0;

    while (size >= num && num_of_bits < 32) {
        num_of_bits++;
        num <<= 1;
    }
    return num_of_bits;
}

static inline int sizeofints(int num_of_ints, unsigned int sizes[]) {
    int i, num;
    unsigned int num_of_bytes, num_of_bits, bytes[32], bytecnt, tmp;
    num_of_bytes = 1;
    bytes[0] = 1;
    num_of_bits = 0;
    for (i = 0; i < num_of_ints; i++) {
        tmp = 0;
        for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++) {
            tmp = bytes[bytecnt] * sizes[i] + tmp;
            bytes[bytecnt] = tmp & 0xff;
            tmp >>= 8;
        }
        while (tmp != 0) {
            bytes[bytecnt++] = tmp & 0xff;
            tmp >>= 8;
        }
        num_of_bytes = bytecnt;
    }
    num = 1;
    num_of_bytes--;
    while (bytes[num_of_bytes] >= num) {
        num_of_bits++;
        num *= 2;
    }
    return num_of_bits + num_of_bytes * 8;
}

static inline int decodebits(int buf[3], int num_of_bits) {
    int cnt, num;
    unsigned int lastbits, lastbyte;
    unsigned char* cbuf;
    int mask = (1 << num_of_bits) - 1;

    cbuf = ((unsigned char*)buf) + 3 * sizeof(*buf);
    cnt = buf[0];
    lastbits = (unsigned int)buf[1];
    lastbyte = (unsigned int)buf[2];

    num = 0;
    while (num_of_bits >= 8) {
        lastbyte = (lastbyte << 8) | cbuf[cnt++];
        num |= (lastbyte >> lastbits) << (num_of_bits - 8);
        num_of_bits -= 8;
    }
    if (num_of_bits > 0) {
        if (lastbits < num_of_bits) {
            lastbits += 8;
            lastbyte = (lastbyte << 8) | cbuf[cnt++];
        }
        lastbits -= num_of_bits;
        num |= (lastbyte >> lastbits) & ((1 << num_of_bits) - 1);
    }
    num &= mask;
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    return num;
}

/*
* decodeints - decode 'small' integers from the buf array
*
* this routine is the inverse from encodeints() and decodes the small integers
* written to buf by calculating the remainder and doing divisions with
* the given sizes[]. You need to specify the total number of bits to be
* used from buf in num_of_bits.
*
*/

static inline void decodeints(int buf[], int num_of_ints, int num_of_bits, unsigned int sizes[], int nums[3]) {

    int bytes[16];
    int i, j, num_of_bytes, p, num, size;

    bytes[1] = bytes[2] = bytes[3] = 0;
    num_of_bytes = 0;
    while (num_of_bits > 8) {
        bytes[num_of_bytes++] = decodebits(buf, 8);
        num_of_bits -= 8;
    }
    if (num_of_bits > 0) {
        bytes[num_of_bytes++] = decodebits(buf, num_of_bits);
    }
    //printf("num_bytes %i\n", num_of_bytes);
    for (i = num_of_ints - 1; i > 0; i--) {
        num = 0;
        size = sizes[i];
        for (j = num_of_bytes - 1; j >= 0; j--) {
            num = (num << 8) | bytes[j];
            p = num / size;
            num = num - p * size;
            bytes[j] = p;
        }
        nums[i] = num;
    }
    nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
}

static const int magicints[] = {
    0,        0,        0,       0,       0,       0,       0,       0,       0,       8,
    10,       12,       16,      20,      25,      32,      40,      50,      64,      80,
    101,      128,      161,     203,     256,     322,     406,     512,     645,     812,
    1024,     1290,     1625,    2048,    2580,    3250,    4096,    5060,    6501,    8192,
    10321,    13003,    16384,   20642,   26007,   32768,   41285,   52015,   65536,   82570,
    104031,   131072,   165140,  208063,  262144,  330280,  416127,  524287,  660561,  832255,
    1048576,  1321122,  1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042, 8388607,
    10568983, 13316085, 16777216};

#define FIRSTIDX 9
/* note that magicints[FIRSTIDX-1] == 0 */
#define LASTIDX (sizeof(magicints) / sizeof(*magicints))

/* Compressed coordinate routines - modified from the original
* implementation by Frans v. Hoesel to make them threadsafe.
*/
bool xdr_decompress_coord_float(float* ptr, int* size, float* precision, md_file_o* xfp) {
    int minint[3], maxint[3], *lip;
    int smallidx;
    unsigned sizeint[3], sizesmall[3], bitsizeint[3], size3;
    int k, *buf1, *buf2, lsize, flag;
    int smallnum, smaller, i, is_smaller, run;
    float *lfp, inv_precision;
    int tmp, *thiscoord, prevcoord[3];
    unsigned int bitsize;
    const float* ptrstart = ptr;

    bitsizeint[0] = 0;
    bitsizeint[1] = 0;
    bitsizeint[2] = 0;

    if (xfp == NULL || ptr == NULL) {
        return false;
    }
    if (!xdr_read_int32(&lsize, 1, xfp)) {
        return false; /* return if we could not read size */
    }
    *size = lsize;
    size3 = *size * 3;

    /* Dont bother with compression for three atoms or less */
    if (*size <= 9) {
        return xdr_read_float(ptr, size3, xfp);
    }


    /* Compression-time if we got here. Read precision first */
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));

    buf1 = md_vm_arena_push(arena, sizeof(int) * size3);
    buf2 = md_vm_arena_push(arena, sizeof(int) * (int)(size3 * 1.2));

    bool result = false;

    if (!xdr_read_float(precision, 1, xfp)) {
        goto done;
    }
    if (!xdr_read_int32(minint, 3, xfp) || !xdr_read_int32(maxint, 3, xfp)) {
        goto done;
    }

    /* buf2[0-2] are special and do not contain actual data */
    buf2[0] = buf2[1] = buf2[2] = 0;

    sizeint[0] = maxint[0] - minint[0] + 1;
    sizeint[1] = maxint[1] - minint[1] + 1;
    sizeint[2] = maxint[2] - minint[2] + 1;

    /* check if one of the sizes is to big to be multiplied */
    if ((sizeint[0] | sizeint[1] | sizeint[2]) > 0xffffff) {
        bitsizeint[0] = sizeofint(sizeint[0]);
        bitsizeint[1] = sizeofint(sizeint[1]);
        bitsizeint[2] = sizeofint(sizeint[2]);
        bitsize = 0; /* flag the use of large sizes */
    } else {
        bitsize = sizeofints(3, sizeint);
    }

    if (!xdr_read_int32(&smallidx, 1, xfp)) {
        return false; /* not sure what has happened here or why we return... */
    }
    tmp = smallidx + 8;
    tmp = smallidx - 1;
    tmp = (FIRSTIDX > tmp) ? FIRSTIDX : tmp;
    smaller = magicints[tmp] / 2;
    smallnum = magicints[smallidx] / 2;
    sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];

    /* buf2[0] holds the length in bytes */
    uint32_t num_bytes;
    if (!xdr_read_int32((int32_t*)&num_bytes, 1, xfp)) {
        goto done;
    }
    if (!xdr_read_bytes((uint8_t*)&(buf2[3]), num_bytes, xfp)) {
        goto done;
    }
    buf2[0] = buf2[1] = buf2[2] = 0;

    lfp = ptr;
    inv_precision = 1.0f / *precision;
    run = 0;
    i = 0;
    lip = buf1;
    while (i < lsize) {
        thiscoord = (int*)(lip) + i * 3;

        if (bitsize == 0) {
            thiscoord[0] = decodebits(buf2, bitsizeint[0]);
            thiscoord[1] = decodebits(buf2, bitsizeint[1]);
            thiscoord[2] = decodebits(buf2, bitsizeint[2]);
        } else {
            decodeints(buf2, 3, bitsize, sizeint, thiscoord);
        }

        i++;
        thiscoord[0] += minint[0];
        thiscoord[1] += minint[1];
        thiscoord[2] += minint[2];

        prevcoord[0] = thiscoord[0];
        prevcoord[1] = thiscoord[1];
        prevcoord[2] = thiscoord[2];

        flag = decodebits(buf2, 1);
        is_smaller = 0;
        if (flag == 1) {
            run = decodebits(buf2, 5);
            is_smaller = run % 3;
            run -= is_smaller;
            is_smaller--;
        }
        if ((lfp - ptrstart) + run > size3) {
            MD_LOG_ERROR("XTC: Buffer overrun during decompression.");
            goto done;
        }
        if (run > 0) {
            thiscoord += 3;
            for (k = 0; k < run; k += 3) {
                decodeints(buf2, 3, smallidx, sizesmall, thiscoord);
                i++;
                thiscoord[0] += prevcoord[0] - smallnum;
                thiscoord[1] += prevcoord[1] - smallnum;
                thiscoord[2] += prevcoord[2] - smallnum;
                if (k == 0) {
                    /* interchange first with second atom for better
                    * compression of water molecules
                    */
                    tmp = thiscoord[0];
                    thiscoord[0] = prevcoord[0];
                    prevcoord[0] = tmp;
                    tmp = thiscoord[1];
                    thiscoord[1] = prevcoord[1];
                    prevcoord[1] = tmp;
                    tmp = thiscoord[2];
                    thiscoord[2] = prevcoord[2];
                    prevcoord[2] = tmp;
                    *lfp++ = prevcoord[0] * inv_precision;
                    *lfp++ = prevcoord[1] * inv_precision;
                    *lfp++ = prevcoord[2] * inv_precision;
                } else {
                    prevcoord[0] = thiscoord[0];
                    prevcoord[1] = thiscoord[1];
                    prevcoord[2] = thiscoord[2];
                }
                *lfp++ = thiscoord[0] * inv_precision;
                *lfp++ = thiscoord[1] * inv_precision;
                *lfp++ = thiscoord[2] * inv_precision;
            }
        } else {
            *lfp++ = thiscoord[0] * inv_precision;
            *lfp++ = thiscoord[1] * inv_precision;
            *lfp++ = thiscoord[2] * inv_precision;
        }
        smallidx += is_smaller;
        if (is_smaller < 0) {
            smallnum = smaller;

            if (smallidx > FIRSTIDX) {
                smaller = magicints[smallidx - 1] / 2;
            } else {
                smaller = 0;
            }
        } else if (is_smaller > 0) {
            smaller = smallnum;
            smallnum = magicints[smallidx] / 2;
        }
        sizesmall[0] = sizesmall[1] = sizesmall[2] = magicints[smallidx];
        if (sizesmall[0] == 0 || sizesmall[1] == 0 || sizesmall[2] == 0) {
            MD_LOG_ERROR("XTC: Invalid size found in 'xdrfile_decompress_coord_float'.");
            goto done;
        }
    }
    result = true;
done:
    md_vm_arena_destroy(arena);
    return result;
}




static bool xtc_frame_header(XDRFILE* xd, int* natoms, int* step, float* time, float box[3][3]) {
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
    if ((xdrfile_read_float((float*)box, DIM * DIM, xd)) != DIM * DIM) {
        MD_LOG_ERROR("XTC: Failed to read box dimensions");
        return false;
    }

    return true;
}

bool md_xtc_read_frame_header(md_file_o* xdr, int* natoms, int* step, float* time, float box[3][3]) {
    int magic;

    if (!xdr_read_int32(&magic, 1, xdr)) {
        MD_LOG_ERROR("XTC: Failed to read magic number in header");
        return false;
    }
    if (magic != XTC_MAGIC) {
        MD_LOG_ERROR("XTC: Magic number did not match");
        return false;
    }
    if (!xdr_read_int32(natoms, 1, xdr)) {
        MD_LOG_ERROR("XTC: Failed to read number of atoms");
        return false;
    }
    if (!xdr_read_int32(step, 1, xdr)) {
        MD_LOG_ERROR("XTC: Failed to read step");
        return false;
    }
    if (!xdr_read_float(time, 1, xdr)) {
        MD_LOG_ERROR("XTC: Failed to read timestamp");
        return false;
    }
    if (!xdr_read_float((float*)box, 9, xdr)) {
        MD_LOG_ERROR("XTC: Failed to read box dimensions");
        return false;
    }

    return true;
}

bool md_xtc_read_frame_offsets_and_times(md_file_o* xdr, md_array(int64_t)* frame_offsets, md_array(double)* frame_times, md_allocator_i* alloc) {
    int step, natoms;
    float time;
    float box[3][3];

    size_t filesize = (size_t)md_file_size(xdr);

    if (filesize == 0) {
        MD_LOG_ERROR("XTC: Failed extract filesize");
        return false;
    }

    /* Go to file beg */
    if (!md_file_seek(xdr, 0, MD_FILE_BEG)) {
        MD_LOG_ERROR("XTC: Failed to seek to beginning of file");
        return false;
    }

    if (!md_xtc_read_frame_header(xdr, &natoms, &step, &time, box)) {
        return false;
    }

    /* Dont bother with compression for nine atoms or less */
    if (natoms <= 9) {
        const size_t framebytes = XTC_SMALL_HEADER_SIZE + XTC_SMALL_COORDS_SIZE * natoms;
        const size_t num_frames = (filesize / framebytes); /* Should we complain if framesize doesn't divide filesize? */

        md_array_ensure(*frame_offsets, num_frames, alloc);
        md_array_ensure(*frame_times,   num_frames, alloc);

        md_array_push(*frame_offsets, 0, alloc);
        md_array_push(*frame_times, time, alloc);

        for (size_t i = 1; i < num_frames; i++) {
            const size_t offset = i * framebytes;

            if (!md_file_seek(xdr, offset, MD_FILE_BEG) || !md_xtc_read_frame_header(xdr, &natoms, &step, &time, box)) {
                goto fail;
            }

            md_array_push(*frame_offsets, i * framebytes, alloc);
            md_array_push(*frame_times, time, alloc);
        }
        md_array_push(*frame_offsets, num_frames * framebytes, alloc);

        return true;
    } else {
        int framebytes, est_nframes;

        /* Move pos back to end of first header */
        if (!md_file_seek(xdr, (int64_t)XTC_HEADER_SIZE, MD_FILE_BEG)) {
            return false;
        }

        if (!xdr_read_int32(&framebytes, 1, xdr)) {
            MD_LOG_ERROR("XTC: Failed to read framebytes");
            return false;
        }
        framebytes = (framebytes + 3) & ~0x03; /* Rounding to the next 32-bit boundary */
        est_nframes = (int)(filesize / ((int64_t)(framebytes + XTC_HEADER_SIZE)) + 1); /* must be at least 1 for successful growth */
        /* First `framebytes` might be larger than average, so we would underestimate `est_nframes`
        */
        est_nframes += est_nframes / 5;

        /* Skip `framebytes` */
        if (!md_file_seek(xdr, (int64_t)(framebytes), MD_FILE_CUR)) {
            goto fail;
        }

        md_array_ensure(*frame_offsets, (size_t)est_nframes, alloc);
        md_array_ensure(*frame_times,   (size_t)est_nframes, alloc);

        md_array_push(*frame_offsets, 0, alloc);
        md_array_push(*frame_times, time, alloc);

        while (true) {
            const int64_t offset = md_file_tell(xdr);
            if (offset == (int64_t)filesize) {
                // Add last offset
                md_array_push(*frame_offsets, offset, alloc);
                return true;
            }
            if (!md_xtc_read_frame_header(xdr, &natoms, &step, &time, box)) {
                goto fail;
            }

            /* Store position in `offsets`, adjust for header */
            md_array_push(*frame_offsets, offset, alloc);
            md_array_push(*frame_times, time, alloc);

            if (!md_file_seek(xdr, offset + XTC_HEADER_SIZE, MD_FILE_BEG)) {
                goto fail;
            }
            /* Read how much to skip */
            if (!xdr_read_int32(&framebytes, 1, xdr)) {
                goto fail;
            }
            
            framebytes = (framebytes + 3) & ~0x03; /* Rounding to the next 32-bit boundary */

            /* Skip `framebytes` to next header */
            if (!md_file_seek(xdr, framebytes, MD_FILE_CUR)) {
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

size_t md_xtc_read_frame_coords(md_file_o* xdr_file, float* coords, size_t capacity) {
    int natoms;
    float prec;

    bool result = xdr_decompress_coord_float(coords, &natoms, &prec, xdr_file);
    if (!result) {
        MD_LOG_ERROR("XTC: Failed to read coordinates");
        return false;
    }

    return true;
}

static size_t xtc_frame_offsets_and_times(XDRFILE* xd, md_array(int64_t)* frame_offsets, md_array(double)* frame_times, md_allocator_i* alloc) {
    int step, natoms;
    float time;
    float box[3][3];

    /* Go to file beg */
    if (xdr_seek(xd, 0L, SEEK_SET) != exdrOK) {
        MD_LOG_ERROR("XTC: Failed to seek to beginning of file");
        return 0;
    }

    if (!xtc_frame_header(xd, &natoms, &step, &time, box)) {
        return 0;
    }

    /* Go to file end */
    if (xdr_seek(xd, 0L, SEEK_END) != exdrOK) {
        return 0;
    }
    /* Cursor position is equivalent to file size */
    size_t filesize = (size_t)xdr_tell(xd);

    /* Dont bother with compression for nine atoms or less */
    if (natoms <= 9) {
        const size_t framebytes = XTC_SMALL_HEADER_SIZE + XTC_SMALL_COORDS_SIZE * natoms;
        const size_t num_frames = (filesize / framebytes); /* Should we complain if framesize doesn't divide filesize? */

        md_array_push(*frame_offsets, 0, alloc);
        md_array_push(*frame_times, time, alloc);

        for (size_t i = 1; i < num_frames; i++) {
            const size_t offset = i * framebytes;

            if (xdr_seek(xd, offset, SEEK_SET) != exdrOK || !xtc_frame_header(xd, &natoms, &step, &time, box)) {
                md_array_free(*frame_offsets, alloc);
                md_array_free(*frame_times, alloc);
                *frame_offsets = 0;
                *frame_times = 0;
                return 0;
            }

            md_array_push(*frame_offsets, i * framebytes, alloc);
            md_array_push(*frame_times, time, alloc);
        }
        md_array_push(*frame_offsets, num_frames * framebytes, alloc);

        return num_frames;
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
        
        md_array_ensure(*frame_offsets, (size_t)est_nframes, alloc);
        md_array_ensure(*frame_times,   (size_t)est_nframes, alloc);

        md_array_push(*frame_offsets, 0, alloc);
        md_array_push(*frame_times, time, alloc);
        size_t num_frames = 1;

        while (1) {
            const int64_t offset = xdr_tell(xd);
            if (offset == (int64_t)filesize) {
                // Add last offset
                md_array_push(*frame_offsets, offset, alloc);
                return num_frames;
            }
            if (!xtc_frame_header(xd, &natoms, &step, &time, box)) {
                goto fail;
            }

            /* Store position in `offsets`, adjust for header */
            md_array_push(*frame_offsets, offset, alloc);
            md_array_push(*frame_times, time, alloc);
            num_frames += 1;

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
        return 0;
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
static size_t xtc_fetch_frame_data(struct md_trajectory_o* inst, int64_t frame_idx, void* frame_data_ptr) {
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

    if (frame_idx < 0 || (int64_t)xtc->header.num_frames <= frame_idx) {
        MD_LOG_ERROR("XTC: Frame index is out of range");
        return 0;
    }

    const int64_t beg = xtc->frame_offsets[frame_idx];
    const int64_t end = xtc->frame_offsets[frame_idx + 1];
    const size_t frame_size = (size_t)MAX(0, end - beg);

    if (frame_data_ptr) {
        ASSERT(xtc->file);
        md_mutex_lock(&xtc->mutex);
        // Seek and read must be an atomic operation to avoid race conditions
        xdr_seek(xtc->file, beg, SEEK_SET);
        const size_t bytes_read = xdr_read(xtc->file, frame_data_ptr, frame_size);
        md_mutex_unlock(&xtc->mutex);
        (void)bytes_read;
        ASSERT(frame_size == bytes_read);
    }
    return frame_size;
}

static bool xtc_decode_frame_data(struct md_trajectory_o* inst, const void* frame_data_ptr, size_t frame_data_size, md_trajectory_frame_header_t* header, float* x, float* y, float* z) {
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

        if (x && y && z) {
            size_t byte_size = natoms * sizeof(rvec);
            rvec* pos = md_alloc(md_get_heap_allocator(), byte_size);
            result = xtc_frame_coords(file, natoms, pos);
            if (result) {            
                // nm -> Ångström
                for (int i = 0; i < natoms; ++i) {
                    x[i] = pos[i][0] * 10.0f;
                    y[i] = pos[i][1] * 10.0f;
                    z[i] = pos[i][2] * 10.0f;
                }
            }
            md_free(md_get_heap_allocator(), pos, byte_size);
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

    bool result = false;
    const size_t frame_size = xtc_fetch_frame_data(inst, frame_idx, NULL);
    if (frame_size > 0) {
        md_allocator_i* alloc = md_get_heap_allocator();

        void* frame_data = md_alloc(alloc, frame_size);
        ASSERT(frame_data);

        const size_t read_size = xtc_fetch_frame_data(inst, frame_idx, frame_data);
        if (read_size != frame_size) {
            MD_LOG_ERROR("Failed to read the expected size");
            goto done;
        }
        result = xtc_decode_frame_data(inst, frame_data, frame_size, header, x, y, z);
    done:
        md_free(alloc, frame_data, frame_size);
    }

    return result;
}

typedef struct xtc_cache_t {
    md_trajectory_cache_header_t header;
    int64_t* frame_offsets;
    double*  frame_times;
} xtc_cache_t;

static bool try_read_cache(xtc_cache_t* cache, str_t cache_file, size_t traj_num_bytes, md_allocator_i* alloc) {
    ASSERT(cache);
    ASSERT(alloc);

    bool result = false;
    md_file_o* file = md_file_open(cache_file, MD_FILE_READ | MD_FILE_BINARY);
    if (file) {
        if (md_file_read(file, &cache->header, sizeof(cache->header)) != sizeof(cache->header)) {
            MD_LOG_ERROR("XTC trajectory cache: failed to read header");
            goto done;
        }

        if (cache->header.magic != MD_XTC_CACHE_MAGIC) {
            MD_LOG_ERROR("XTC trajectory cache: magic was incorrect or corrupt");
            goto done;
        }
        if (cache->header.version != MD_XTC_CACHE_VERSION) {
            MD_LOG_INFO("XTC trajectory cache: version mismatch, expected %i, got %i", MD_XTC_CACHE_VERSION, (int)cache->header.version);
        }
        if (cache->header.num_bytes != traj_num_bytes) {
            MD_LOG_INFO("XTC trajectory cache: trajectory size mismatch, expected %zu, got %zu", traj_num_bytes, cache->header.num_bytes);
        }
        if (cache->header.num_atoms == 0) {
            MD_LOG_ERROR("XTC trajectory cache: num atoms was zero");
            goto done;
        }
        if (cache->header.num_frames == 0) {
            MD_LOG_ERROR("XTC trajectory cache: num frames was zero");
            goto done;
        }

        const size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
        cache->frame_offsets = md_alloc(alloc, offset_bytes);
        if (md_file_read(file, cache->frame_offsets, offset_bytes) != offset_bytes) {
            MD_LOG_ERROR("XTC trajectory cache: Failed to read offset data");
            md_free(alloc, cache->frame_offsets, offset_bytes);
            goto done;
        }

        const size_t time_bytes = cache->header.num_frames * sizeof(double);
        cache->frame_times = md_alloc(alloc, time_bytes);
        if (md_file_read(file, cache->frame_times, time_bytes) != time_bytes) {
        	MD_LOG_ERROR("XTC trajectory cache: times are incomplete");
        	md_free(alloc, cache->frame_offsets, offset_bytes);
        	md_free(alloc, cache->frame_times, time_bytes);
        	goto done;
        }

        // Test position in file, we expect to be at the end of the file
        if (md_file_tell(file) != (int64_t)md_file_size(file)) {
        	MD_LOG_ERROR("XTC trajectory cache: file position was not at the end of the file");
        	md_free(alloc, cache->frame_offsets, offset_bytes);
        	md_free(alloc, cache->frame_times, time_bytes);
        	goto done;
        }

        result = true;
    done:
        md_file_close(file);
    }
    return result;
}

static bool write_cache(const xtc_cache_t* cache, str_t cache_file) {
    bool result = false;

    md_file_o* file = md_file_open(cache_file, MD_FILE_WRITE | MD_FILE_BINARY);
    if (!file) {
        MD_LOG_INFO("XTC trajectory cache: could not open file '"STR_FMT"'", STR_ARG(cache_file));
        return false;
    }

    if (md_file_write(file, &cache->header, sizeof(cache->header)) != sizeof(cache->header)) {
        MD_LOG_ERROR("XTC trajectory cache: failed to write header");
        goto done;
    }

    const size_t offset_bytes = (cache->header.num_frames + 1) * sizeof(int64_t);
    if (md_file_write(file, cache->frame_offsets, offset_bytes) != offset_bytes) {
        MD_LOG_ERROR("Failed to write offset cache, offsets");
        goto done;
    }

    const size_t time_bytes = cache->header.num_frames * sizeof(double);
    if (md_file_write(file, cache->frame_times, time_bytes) != time_bytes) {
    	MD_LOG_ERROR("Failed to write offset cache, times");
    	goto done;
    }

    result = true;

done:
    md_file_close(file);
    return result;
}

md_trajectory_i* md_xtc_trajectory_create(str_t filename, md_allocator_i* ext_alloc, uint32_t flags) {
    ASSERT(ext_alloc);
    md_allocator_i* alloc = md_arena_allocator_create(ext_alloc, MEGABYTES(1));

    // Ensure that the path is zero terminated (not guaranteed by str_t)
    md_strb_t sb = md_strb_create(md_get_temp_allocator());
    md_strb_push_str(&sb, filename);
    XDRFILE* file = xdrfile_open(md_strb_to_cstr(sb), "r");

    if (file) {
        xdr_seek(file, 0L, SEEK_END);
        const size_t filesize = (size_t)xdr_tell(file);
        xdr_seek(file, 0L, SEEK_SET);

        int num_atoms, step;
        float time;
        float box[3][3];
        if (!xtc_frame_header(file, &num_atoms, &step, &time, box)) {
            goto fail;
        }

        if (num_atoms == 0) {
            MD_LOG_ERROR("XTC: Number of atoms in trajectory was zero");
            goto fail;
        }

        md_strb_push_str(&sb, STR_LIT(".cache"));
        str_t path = md_strb_to_str(sb);

        xtc_cache_t cache = {0};
        if (!try_read_cache(&cache, path, filesize, alloc)) {
            cache.header.magic     = MD_XTC_CACHE_MAGIC;
            cache.header.version   = MD_XTC_CACHE_VERSION;
            cache.header.num_bytes = filesize;
            cache.header.num_atoms = num_atoms;
            cache.header.num_frames = xtc_frame_offsets_and_times(file, &cache.frame_offsets, &cache.frame_times, alloc);
            if (!cache.header.num_frames) {
                goto fail;
            }
            if (!cache.frame_offsets || !cache.frame_times) {
                MD_LOG_DEBUG("XTC: frame offsets or frame times was empty");
                goto fail;
            }

            if (!(flags & MD_TRAJECTORY_FLAG_DISABLE_CACHE_WRITE)) {
                // If we fail to write the cache, that's ok, we can inform about it, but do not halt
                if (write_cache(&cache, path)) {
                    MD_LOG_INFO("XTC: Successfully created cache file for '" STR_FMT "'", STR_ARG(path));
                }
            }
        }

        size_t max_frame_size = 0;
        for (size_t i = 0; i < cache.header.num_frames; ++i) {
            const size_t frame_size = (size_t)MAX(0, cache.frame_offsets[i + 1] - cache.frame_offsets[i]);
            max_frame_size = MAX(max_frame_size, frame_size);
        }

        void* mem = md_alloc(alloc, sizeof(md_trajectory_i) + sizeof(xtc_t));
        ASSERT(mem);
        MEMSET(mem, 0, sizeof(md_trajectory_i) + sizeof(xtc_t));

        md_trajectory_i* traj = mem;
        xtc_t* xtc = (xtc_t*)(traj + 1);

        xtc->magic = MD_XTC_TRAJ_MAGIC;
        xtc->allocator = alloc;
        xtc->file = file;
        xtc->frame_offsets = cache.frame_offsets;
        xtc->mutex = md_mutex_create();

        xtc->header = (md_trajectory_header_t) {
            .num_frames = cache.header.num_frames,
            .num_atoms = num_atoms,
            .max_frame_data_size = max_frame_size,
            .time_unit = md_unit_pikosecond(),
            .frame_times = cache.frame_times,
        };

        traj->inst = (struct md_trajectory_o*)xtc;
        traj->get_header = xtc_get_header;
        traj->load_frame = xtc_load_frame;
        //traj->fetch_frame_data = xtc_fetch_frame_data;
        //traj->decode_frame_data = xtc_decode_frame_data;

        return traj;
    }
fail:
    if (file) xdrfile_close(file);
    md_arena_allocator_destroy(alloc);
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

    if (xtc->file) xdrfile_close(xtc->file);
    md_mutex_destroy(&xtc->mutex);
    md_arena_allocator_destroy(xtc->allocator);
}

static md_trajectory_loader_i xtc_loader = {
    md_xtc_trajectory_create,
    md_xtc_trajectory_free,
};

md_trajectory_loader_i* md_xtc_trajectory_loader(void) {
    return &xtc_loader;
}
