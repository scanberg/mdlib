#include <md_dcd.h>
#include <md_trajectory.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_unit.h>

#include <math.h>

#define MD_DCD_TRAJ_MAGIC 0x3F5C8A1B7E902D46ULL

// 1 AKMA time unit = 48.88821 fs = 0.04888821 ps
// NAMD's TIMEFACTOR for converting DCD timestamps to femtoseconds
#define DCD_AKMA_TO_PS 0.04888821

// CHARMm / X-PLOR DCD flags
#define DCD_IS_CHARMM       0x01
#define DCD_HAS_4DIMS       0x02
#define DCD_HAS_EXTRA_BLOCK 0x04

// ==================== Internal DCD file header ====================

// Holds metadata parsed directly from the DCD binary header.
typedef struct dcd_file_header_t {
    int32_t  nset;           // Number of frames declared in the header (may be 0 if file was not closed properly)
    int32_t  istart;         // First timestep number saved
    int32_t  nsavc;          // Integration steps between successive saves
    int32_t  nfixed;         // Number of fixed atoms (0 = all atoms move)
    int32_t  natoms;         // Total number of atoms
    int32_t  charmm;         // CHARMm flags (DCD_IS_CHARMM, DCD_HAS_EXTRA_BLOCK, DCD_HAS_4DIMS)
    double   delta;          // Integration timestep (AKMA units)
    bool     reverse_endian; // File byte order differs from host
    int64_t  header_size;    // Byte offset to the first frame data
    int32_t* free_indices;   // 1-based indices of free (non-fixed) atoms; NULL when nfixed == 0
} dcd_file_header_t;

// ==================== Trajectory instance ====================

typedef struct dcd_t {
    uint64_t magic;
    md_file_t file;
    md_allocator_i* allocator;

	dcd_file_header_t header;
	size_t num_frames;           // Actual frame count (may differ from header.nset when file was not closed cleanly)

	// Parameters for analytical frame offset computation.
	// offset[0] = file_header_size
	// offset[i] = file_header_size + first_frame_size + (i-1) * frame_size   (i >= 1)
	size_t file_header_size;
	size_t first_frame_size;
	size_t frame_size;

	double* frame_times;          // Precomputed frame timestamps (picoseconds), length = num_frames.
	md_unit_t time_unit;          // Time unit for frame_times (we try to convert to picoseconds if time data is available, otherwise we fallback to indices)

    // Full coordinate snapshot of frame 0, used for fixed-atom reconstruction in subsequent frames.
    // Only allocated when nfixed > 0.
    float* first_frame_x;
    float* first_frame_y;
    float* first_frame_z;
} dcd_t;

// ==================== Byte-swap helpers ====================

static inline int32_t dcd_swap32(int32_t v, bool swap) {
    return swap ? (int32_t)BSWAP32((uint32_t)v) : v;
}

static inline void dcd_swap32_array(int32_t* data, int count, bool swap) {
    if (!swap) return;
    for (int i = 0; i < count; ++i) {
        data[i] = BSWAP32(data[i]);
    }
}

// ==================== Analytical frame offset ====================

// DCD frames have fixed sizes so the file offset for any frame is a closed-form
// expression.  No offset table is needed.
static inline int64_t dcd_frame_offset(const dcd_t* dcd, int64_t frame_idx) {
    if (frame_idx == 0) return (int64_t)dcd->file_header_size;
    return (int64_t)dcd->file_header_size + (int64_t)dcd->first_frame_size + (frame_idx - 1) * (int64_t)dcd->frame_size;
}

// Compute the byte sizes of the first and subsequent frames from header fields.
// Each coordinate dimension is a Fortran record: [int32: N*4][float[N]][int32: N*4]
// = (N + 2) * 4 bytes per dimension.
static void dcd_compute_frame_sizes(const dcd_file_header_t* h,
                                    int64_t* out_first,
                                    int64_t* out_subsequent)
{
    const int64_t extrablocksize = (h->charmm & DCD_HAS_EXTRA_BLOCK) ? (48 + 8) : 0;
    const int64_t ndims          = (h->charmm & DCD_HAS_4DIMS)       ? 4 : 3;
    *out_first      = (int64_t)(h->natoms              + 2) * ndims * 4 + extrablocksize;
    *out_subsequent = (int64_t)(h->natoms - h->nfixed  + 2) * ndims * 4 + extrablocksize;
}

// ==================== DCD file header parsing ====================

// Parses the full DCD file header starting from the beginning of the file.
// On success, populates *out and leaves the file positioned at the first frame.
// When nfixed > 0, allocates out->free_indices from alloc.
static bool dcd_parse_file_header(md_file_t file, dcd_file_header_t* out, md_allocator_i* alloc) {
    ASSERT(md_file_valid(file) && out && alloc);
    MEMSET(out, 0, sizeof(*out));

    // The very first 4 bytes are a Fortran record marker that must equal 84.
    // If the bytes compare as 84 when swapped, the file uses reverse endianness.
    int32_t marker = 0;
    if (md_file_read(file, &marker, 4) != 4) {
        MD_LOG_ERROR("DCD: Failed to read initial record marker");
        return false;
    }
    if (marker != 84) {
        if ((int32_t)BSWAP32((uint32_t)marker) == 84) {
            out->reverse_endian = true;
        } else {
            MD_LOG_ERROR("DCD: Invalid initial record marker (got %d, expected 84)", marker);
            return false;
        }
    }

    const bool rev = out->reverse_endian;

    // Read the 84-byte main header block
    char hdrbuf[84];
    if (md_file_read(file, hdrbuf, 84) != 84) {
        MD_LOG_ERROR("DCD: Failed to read 84-byte header block");
        return false;
    }

    // Trailing record marker must also equal 84
    int32_t trailer = 0;
    if (md_file_read(file, &trailer, 4) != 4) {
        MD_LOG_ERROR("DCD: Failed to read trailing record marker");
        return false;
    }
    if (dcd_swap32(trailer, rev) != 84) {
        MD_LOG_ERROR("DCD: Trailing record marker mismatch");
        return false;
    }

    // Verify the 4-byte file type identifier "CORD"
    if (hdrbuf[0] != 'C' || hdrbuf[1] != 'O' || hdrbuf[2] != 'R' || hdrbuf[3] != 'D') {
        MD_LOG_ERROR("DCD: File identifier is not 'CORD'");
        return false;
    }

    // Extract the integer fields packed into the 84-byte block.
    // Offsets are relative to the start of hdrbuf (i.e. 4 bytes after the initial marker).
    int32_t nset, istart, nsavc, namnf, charmm_version;
    MEMCPY(&nset,           hdrbuf +  4, 4);
    MEMCPY(&istart,         hdrbuf +  8, 4);
    MEMCPY(&nsavc,          hdrbuf + 12, 4);
    MEMCPY(&namnf,          hdrbuf + 36, 4);
    MEMCPY(&charmm_version, hdrbuf + 80, 4);
    nset           = dcd_swap32(nset,           rev);
    istart         = dcd_swap32(istart,         rev);
    nsavc          = dcd_swap32(nsavc,          rev);
    namnf          = dcd_swap32(namnf,          rev);
    charmm_version = dcd_swap32(charmm_version, rev);

    // A non-zero value in the last int slot of the block marks this as a CHARMm file.
    int charmm = 0;
    if (charmm_version != 0) {
        charmm = DCD_IS_CHARMM;

        int32_t has_extra, has_4dims;
        MEMCPY(&has_extra, hdrbuf + 44, 4);
        MEMCPY(&has_4dims, hdrbuf + 48, 4);
        has_extra = dcd_swap32(has_extra, rev);
        has_4dims = dcd_swap32(has_4dims, rev);

        if (has_extra)    charmm |= DCD_HAS_EXTRA_BLOCK;
        if (has_4dims == 1) charmm |= DCD_HAS_4DIMS;
    }

    // DELTA is stored as a 32-bit float for CHARMm and as a 64-bit double for X-PLOR.
    // In the CHARMm case the float sits at hdrbuf[40]; the two adjacent CHARMm flags
    // at [44] and [48] would overlap a 64-bit double at [40], hence the split.
    double delta = 0.0;
    if (charmm & DCD_IS_CHARMM) {
        uint32_t raw;
        MEMCPY(&raw, hdrbuf + 40, 4);
        if (rev) raw = BSWAP32(raw);
        float fdelta;
        MEMCPY(&fdelta, &raw, 4);
        delta = (double)fdelta;
    } else {
        uint64_t raw;
        MEMCPY(&raw, hdrbuf + 40, 8);
        if (rev) raw = BSWAP64(raw);
        MEMCPY(&delta, &raw, 8);
    }

    // Title block: [int32: size][size bytes: titles][int32: size]
    // We don't need the title content, so we skip over it entirely.
    int32_t beg_title_size;
    if (md_file_read(file, &beg_title_size, 4) != 4) {
        MD_LOG_ERROR("DCD: Failed to read title block size");
        return false;
    }
    beg_title_size = dcd_swap32(beg_title_size, rev);

    if ((beg_title_size - 4) % 80 == 0) {
        int32_t num_title;
		md_file_read(file, &num_title, 4);
        num_title = dcd_swap32(num_title, rev);
        if (num_title <= 0) {
            MD_LOG_ERROR("DCD: Title block has non-positive title count");
            return false;
		}
		// Skip the title lines, each 80 bytes long.
		md_file_seek(file, (int64_t)num_title * 80, MD_FILE_CUR);
	}
	else {
        // Malformed title block
        return false;
    }

	int32_t end_title_size;
    if (md_file_read(file, &end_title_size, 4) != 4) {
        MD_LOG_ERROR("DCD: Failed to read title block trailing size");
        return false;
	}
	end_title_size = dcd_swap32(end_title_size, rev);

	if (end_title_size != beg_title_size) {
        MD_LOG_ERROR("DCD: Title block size mismatch");
        return false;
    }


    // Natoms block: [int32: 4][int32: natoms][int32: 4]
    int32_t natoms_sz;
    if (md_file_read(file, &natoms_sz, 4) != 4) {
        MD_LOG_ERROR("DCD: Failed to read natoms block size");
        return false;
    }
    if (dcd_swap32(natoms_sz, rev) != 4) {
        MD_LOG_ERROR("DCD: Natoms record has unexpected block size");
        return false;
    }

    int32_t natoms = 0;
    if (md_file_read(file, &natoms, 4) != 4) {
        MD_LOG_ERROR("DCD: Failed to read natoms");
        return false;
    }
    natoms = dcd_swap32(natoms, rev);

    int32_t natoms_sz2 = 0;
    if (md_file_read(file, &natoms_sz2, 4) != 4) {
        MD_LOG_ERROR("DCD: Failed to read natoms block trailing size");
        return false;
    }

    // Free-atom index block (only present when fixed atoms exist):
    // [int32: nfree*4][int32[nfree]: 1-based indices][int32: nfree*4]
    out->free_indices = NULL;
    if (namnf > 0) {
        const int32_t nfree = natoms - namnf;
        if (nfree <= 0) {
            MD_LOG_ERROR("DCD: Fixed atom count (%d) is >= total atom count (%d)", namnf, natoms);
            return false;
        }

        int32_t idx_sz;
        if (md_file_read(file, &idx_sz, 4) != 4) {
            MD_LOG_ERROR("DCD: Failed to read free-atom index block size");
            return false;
        }
        idx_sz = dcd_swap32(idx_sz, rev);
        if (idx_sz != nfree * 4) {
            MD_LOG_ERROR("DCD: Free-atom index block has unexpected size (%d, expected %d)", idx_sz, nfree * 4);
            return false;
        }

        out->free_indices = (int32_t*)md_alloc(alloc, (size_t)nfree * sizeof(int32_t));
        if (md_file_read(file, out->free_indices, (size_t)nfree * 4) != (size_t)(nfree * 4)) {
            MD_LOG_ERROR("DCD: Failed to read free atom indices");
            md_free(alloc, out->free_indices, (size_t)nfree * sizeof(int32_t));
            out->free_indices = NULL;
            return false;
        }
        if (rev) {
            for (int32_t i = 0; i < nfree; ++i)
                out->free_indices[i] = dcd_swap32(out->free_indices[i], true);
        }

        int32_t idx_sz2 = 0;
        if (md_file_read(file, &idx_sz2, 4) != 4) {
            MD_LOG_ERROR("DCD: Failed to read free-atom index block trailing size");
            return false;
        }
    }

    out->nset        = nset;
    out->istart      = istart;
    out->nsavc       = MAX(nsavc, 1);
    out->nfixed      = namnf;
    out->natoms      = natoms;
    out->charmm      = charmm;
    out->delta       = delta;
    out->header_size = md_file_tell(file);

    return true;
}

// ==================== Unit cell conversion ====================

// The DCD periodic-cell extra block stores 6 doubles in this order:
//   [0] A        length of a-vector
//   [1] cos(γ)   angle between a and b  (or γ in degrees for old NAMD ≤ 2.5)
//   [2] B        length of b-vector
//   [3] cos(β)   angle between a and c  (or β in degrees)
//   [4] cos(α)   angle between b and c  (or α in degrees)
//   [5] C        length of c-vector
//
// We distinguish cosines from degree values by checking whether the angle
// slot at index 1 lies within the valid cosine range [-1, 1].
static md_unitcell_t dcd_unitcell_from_params(const double uc[6]) {
    double a = uc[0], b = uc[2], c = uc[5];
    double alpha_deg, beta_deg, gamma_deg;

    if (fabs(uc[1]) <= 1.0) {
        gamma_deg = RAD_TO_DEG(acos(uc[1]));
        beta_deg  = RAD_TO_DEG(acos(uc[3]));
        alpha_deg = RAD_TO_DEG(acos(uc[4]));
    } else {
        gamma_deg = uc[1];
        beta_deg  = uc[3];
        alpha_deg = uc[4];
    }

    return md_unitcell_from_extent_and_angles(a, b, c, alpha_deg, beta_deg, gamma_deg);
}

// ==================== Low-level frame reading ====================

static bool dcd_read_bytes_at(md_file_t file, md_file_offset_t* offset, void* ptr, size_t num_bytes) {
    ASSERT(md_file_valid(file));
    ASSERT(offset);

    if (md_file_read_at(file, *offset, ptr, num_bytes) != num_bytes) {
        return false;
    }

    *offset += (md_file_offset_t)num_bytes;
    return true;
}

// Read one Fortran coordinate record: [int32: N*4][float[N]][int32: N*4].
// When out is NULL the data is skipped.  Returns the atom count on success, -1 on error.
// expected_count is used to validate the record size.
static int dcd_read_coord_record(md_file_t file, md_file_offset_t* offset, int expected_count, bool rev, float* out) {
    int32_t sz = 0;
    if (!dcd_read_bytes_at(file, offset, &sz, 4)) return -1;
    sz = dcd_swap32(sz, rev);

    const int n = sz / 4;
    if (n != expected_count) {
        MD_LOG_ERROR("DCD: Coordinate record size mismatch (got %d atoms, expected %d)", n, expected_count);
        return -1;
    }

    if (out) {
        if (!dcd_read_bytes_at(file, offset, out, (size_t)sz)) return -1;
        dcd_swap32_array((int32_t*)out, n, rev);
    } else {
        *offset += sz;
    }

    int32_t sz2 = 0;
    if (!dcd_read_bytes_at(file, offset, &sz2, 4)) return -1;

    return n;
}

// Read a complete frame from the current file position.
//
// For the first frame (or when nfixed == 0) all natoms coordinates are read directly
// into out_x/y/z.  For subsequent frames with fixed atoms, only nfree = natoms - nfixed
// values are read and scattered into out_x/y/z at the positions given by free_indices.
// The caller is responsible for pre-filling out_x/y/z with the first-frame coordinates
// before calling this function in that case.
//
// Any output pointer may be NULL (the corresponding data will be skipped).
static bool dcd_read_frame_at(md_file_t file,
                               md_file_offset_t* offset,
                               int natoms, int nfixed, int charmm, bool rev,
                               const int32_t* free_indices,
                               bool is_first_frame,
                               md_unitcell_t* out_unitcell,
                               float* out_x, float* out_y, float* out_z)
{
    ASSERT(md_file_valid(file));
    ASSERT(offset);
    const bool coords_requested = (out_x && out_y && out_z);
    const int  nfree            = natoms - nfixed;
    const bool all_atoms        = (nfixed == 0 || is_first_frame);
    const int  atoms_in_frame   = all_atoms ? natoms : nfree;

    // --- Optional periodic cell block (CHARMm only) ---
    if ((charmm & DCD_IS_CHARMM) && (charmm & DCD_HAS_EXTRA_BLOCK)) {
        int32_t block_sz = 0;
        if (!dcd_read_bytes_at(file, offset, &block_sz, 4)) {
            MD_LOG_ERROR("DCD: Failed to read extra block size");
            return false;
        }
        block_sz = dcd_swap32(block_sz, rev);

        if (block_sz == 48) {
            double uc[6];
            if (!dcd_read_bytes_at(file, offset, uc, 48)) {
                MD_LOG_ERROR("DCD: Failed to read unit cell data");
                return false;
            }
            if (rev) {
                for (int i = 0; i < 6; ++i) {
                    uint64_t tmp;
                    MEMCPY(&tmp, &uc[i], 8);
                    tmp = BSWAP64(tmp);
                    MEMCPY(&uc[i], &tmp, 8);
                }
            }
            if (out_unitcell)
                *out_unitcell = dcd_unitcell_from_params(uc);
        } else {
            // Unknown block size: skip the data
            *offset += block_sz;
        }

        int32_t block_sz2 = 0;
        if (!dcd_read_bytes_at(file, offset, &block_sz2, 4)) {
            MD_LOG_ERROR("DCD: Failed to read extra block trailing size");
            return false;
        }
    }

    // --- Coordinate records ---
    if (all_atoms) {
        // First frame or no fixed atoms: read all atom coordinates directly.
        if (dcd_read_coord_record(file, offset, natoms, rev, coords_requested ? out_x : NULL) < 0) {
            MD_LOG_ERROR("DCD: Failed to read X coordinate record");
            return false;
        }
        if (dcd_read_coord_record(file, offset, natoms, rev, coords_requested ? out_y : NULL) < 0) {
            MD_LOG_ERROR("DCD: Failed to read Y coordinate record");
            return false;
        }
        if (dcd_read_coord_record(file, offset, natoms, rev, coords_requested ? out_z : NULL) < 0) {
            MD_LOG_ERROR("DCD: Failed to read Z coordinate record");
            return false;
        }
    } else {
        // Subsequent frame with fixed atoms: read only the nfree mobile coordinates
        // and scatter them into the correct positions in the output buffers.
        // The caller has already filled out_x/y/z with the first-frame values.
        if (coords_requested) {
            md_allocator_i* heap = md_get_heap_allocator();
            float* tmp = (float*)md_alloc(heap, (size_t)nfree * sizeof(float));

            for (int dim = 0; dim < 3; ++dim) {
                float* dst = (dim == 0) ? out_x : (dim == 1) ? out_y : out_z;
                if (dcd_read_coord_record(file, offset, nfree, rev, tmp) < 0) {
                    MD_LOG_ERROR("DCD: Failed to read free-atom coordinate record (dim %d)", dim);
                    md_free(heap, tmp, (size_t)nfree * sizeof(float));
                    return false;
                }
                for (int i = 0; i < nfree; ++i)
                    dst[free_indices[i] - 1] = tmp[i];
            }

            md_free(heap, tmp, (size_t)nfree * sizeof(float));
        } else {
            for (int dim = 0; dim < 3; ++dim) {
                if (dcd_read_coord_record(file, offset, nfree, rev, NULL) < 0) {
                    MD_LOG_ERROR("DCD: Failed to skip coordinate record (dim %d)", dim);
                    return false;
                }
            }
        }
    }

    // --- Optional 4th-dimension record (CHARMm only) ---
    if (charmm & DCD_HAS_4DIMS) {
        if (dcd_read_coord_record(file, offset, atoms_in_frame, rev, NULL) < 0) {
            MD_LOG_ERROR("DCD: Failed to skip 4th-dimension record");
            return false;
        }
    }

    return true;
}

static bool dcd_get_header(struct md_trajectory_o* inst, md_trajectory_header_t* out) {
    ASSERT(inst && out);
    dcd_t* dcd = (dcd_t*)inst;
    if (dcd->magic != MD_DCD_TRAJ_MAGIC) {
        MD_LOG_ERROR("DCD: Invalid trajectory instance in get_header");
        return false;
    }
	*out = (md_trajectory_header_t){
		.num_frames  = dcd->num_frames,
		.num_atoms   = (size_t)dcd->header.natoms,
		.time_unit   = dcd->time_unit,
		.frame_times = dcd->frame_times,
	};
    return true;
}

static bool dcd_load_frame(struct md_trajectory_o* inst, int64_t idx, md_trajectory_frame_header_t* out_hdr, float* x, float* y, float* z) {
    ASSERT(inst);
    dcd_t* dcd = (dcd_t*)inst;
    if (dcd->magic != MD_DCD_TRAJ_MAGIC) {
        MD_LOG_ERROR("DCD: Invalid trajectory instance in load_frame");
        return false;
    }
	if (idx < 0 || idx >= (int64_t)dcd->num_frames) {
		MD_LOG_ERROR("DCD: Frame index out of range (got %lld, valid range is [0, %zu))", idx, dcd->num_frames);
		return false;
	}
	if ((x || y || z) && !(x && y && z)) {
		MD_LOG_ERROR("DCD: Coordinate arrays must all be provided or all be NULL");
		return false;
	}
    if (!md_file_valid(dcd->file)) {
        MD_LOG_ERROR("DCD: Invalid file handle in trajectory instance");
		return false;
	}
    md_file_offset_t offset = dcd_frame_offset(dcd, idx);
	const bool is_first_frame = (idx == 0);
	const int natoms = (int)dcd->header.natoms;
	// Pre-fill output with first-frame coords so fixed-atom slots are correct
	// before dcd_read_frame_at scatters only the free-atom coordinates.
	if (x && dcd->header.nfixed > 0 && !is_first_frame) {
		MEMCPY(x, dcd->first_frame_x, (size_t)natoms * sizeof(float));
		MEMCPY(y, dcd->first_frame_y, (size_t)natoms * sizeof(float));
		MEMCPY(z, dcd->first_frame_z, (size_t)natoms * sizeof(float));
	}
	md_unitcell_t unitcell = md_unitcell_none();
    bool success = dcd_read_frame_at(dcd->file, &offset, natoms, dcd->header.nfixed, dcd->header.charmm, dcd->header.reverse_endian,
									 dcd->header.free_indices,
									 is_first_frame,
									 &unitcell,
									 x, y, z);
    if (!success) {
        MD_LOG_ERROR("DCD: Failed to read frame %lld", idx);
        return false;
    }
	if (out_hdr) {
		*out_hdr = (md_trajectory_frame_header_t){
			.num_atoms = (size_t)natoms,
			.index     = idx,
			.timestamp = dcd->frame_times[idx],
			.unitcell  = unitcell,
		};
	}
	return true;
}

// ==================== Creation and destruction ====================
md_trajectory_i* md_dcd_trajectory_create(str_t filename, md_allocator_i* ext_alloc, md_trajectory_flags_t flags) {
    ASSERT(ext_alloc);
    (void)flags;

    md_allocator_i* alloc = md_arena_allocator_create(ext_alloc, MEGABYTES(1));

    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("DCD: Failed to open '" STR_FMT "'", STR_ARG(filename));
        goto fail;
    }

    const int64_t filesize = (int64_t)md_file_size(file);

    md_file_seek(file, 0, MD_FILE_BEG);

    // Always parse the raw DCD header so we have ground-truth metadata.
    dcd_file_header_t fhdr;
    if (!dcd_parse_file_header(file, &fhdr, alloc)) {
        MD_LOG_ERROR("DCD: Failed to parse header of '" STR_FMT "'", STR_ARG(filename));
        goto fail;
    }

    if (fhdr.natoms <= 0) {
        MD_LOG_ERROR("DCD: Atom count is zero or negative");
        goto fail;
    }

    // Compute fixed per-frame byte sizes analytically from the header fields.
    int64_t first_frame_size, frame_size;
    dcd_compute_frame_sizes(&fhdr, &first_frame_size, &frame_size);

    // Determine how many complete frames fit in the file payload.
    const int64_t payload   = filesize - fhdr.header_size;
    size_t num_frames = 0;
    if (payload >= first_frame_size) {
        num_frames = 1 + (size_t)MAX(0LL, (payload - first_frame_size) / frame_size);
    }

    // If the header declares fewer frames, trust that (the file may be truncated or
    // still being written).
    if (fhdr.nset > 0 && (size_t)(fhdr.nset + 1) < num_frames) {
        num_frames = (size_t)fhdr.nset + 1;
    }
    if (num_frames == 0) {
        MD_LOG_ERROR("DCD: No complete frames found in '" STR_FMT "'", STR_ARG(filename));
        goto fail;
    }

    // Frame timestamps: convert from AKMA to picoseconds.
    double* frame_times = (double*)md_alloc(alloc, num_frames * sizeof(double));
	md_unit_t time_unit = md_unit_none();
    if (fhdr.delta > 0) {
        for (size_t i = 0; i < num_frames; ++i) {
            frame_times[i] = (fhdr.istart + (double)i * fhdr.nsavc) * fhdr.delta * DCD_AKMA_TO_PS;
        }
		time_unit = md_unit_picosecond();
    }
    else {
        // No frame time available, default to frame index
        for (size_t i = 0; i < num_frames; ++i) {
            frame_times[i] = (double)i;
		}
    }

    // Allocate trajectory and internal structs from the arena.
    void* mem = md_alloc(alloc, sizeof(md_trajectory_i) + sizeof(dcd_t));
    ASSERT(mem);
    MEMSET(mem, 0, sizeof(md_trajectory_i) + sizeof(dcd_t));

    md_trajectory_i* traj = (md_trajectory_i*)mem;
    dcd_t*           dcd  = (dcd_t*)(traj + 1);

    dcd->magic            = MD_DCD_TRAJ_MAGIC;
    dcd->file             = file;
    dcd->allocator        = alloc;
    dcd->header           = fhdr;
    dcd->num_frames       = num_frames;
    dcd->file_header_size = (size_t)fhdr.header_size;
    dcd->first_frame_size = (size_t)first_frame_size;
    dcd->frame_size       = (size_t)frame_size;
    dcd->frame_times      = frame_times;
	dcd->time_unit        = time_unit;

    // When fixed atoms are present, store the full first-frame coordinates.
    // These are blended with subsequent frames that only contain free-atom deltas.
    dcd->first_frame_x = NULL;
    dcd->first_frame_y = NULL;
    dcd->first_frame_z = NULL;
    if (fhdr.nfixed > 0) {
        const int natoms   = fhdr.natoms;
		dcd->first_frame_x = (float*)md_alloc(alloc, (size_t)natoms * sizeof(float));
		dcd->first_frame_y = (float*)md_alloc(alloc, (size_t)natoms * sizeof(float));
		dcd->first_frame_z = (float*)md_alloc(alloc, (size_t)natoms * sizeof(float));
		md_file_offset_t first_frame_offset = dcd_frame_offset(dcd, 0);
        // Read the first frame treating nfixed as 0 so that all atoms are read.
		if (!dcd_read_frame_at(file, &first_frame_offset, natoms, 0, fhdr.charmm, fhdr.reverse_endian,
                               NULL, true, NULL,
                               dcd->first_frame_x, dcd->first_frame_y, dcd->first_frame_z))
        {
            MD_LOG_ERROR("DCD: Failed to read first frame for fixed-atom initialisation");
            goto fail;
        }
    }

    traj->inst       = (struct md_trajectory_o*)dcd;
    traj->get_header = dcd_get_header;
    traj->load_frame = dcd_load_frame;

    return traj;

fail:
    if (md_file_valid(file)) md_file_close(&file);
    md_arena_allocator_destroy(alloc);
    return NULL;
}

void md_dcd_trajectory_free(md_trajectory_i* traj) {
    ASSERT(traj && traj->inst);
    dcd_t* dcd = (dcd_t*)traj->inst;
    if (dcd->magic != MD_DCD_TRAJ_MAGIC) {
        MD_LOG_ERROR("DCD: Cannot free trajectory, invalid magic");
        ASSERT(false);
        return;
    }
    if (md_file_valid(dcd->file)) md_file_close(&dcd->file);
    md_arena_allocator_destroy(dcd->allocator);
}

static md_trajectory_loader_i dcd_loader = {
    md_dcd_trajectory_create,
    md_dcd_trajectory_free,
};

md_trajectory_loader_i* md_dcd_trajectory_loader(void) {
    return &dcd_loader;
}
