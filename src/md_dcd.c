#include <md_dcd.h>
#include <md_system.h>
#include <md_util.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include <core/md_unit.h>

#include <math.h>
#include <stdio.h>


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
    int32_t  version;        // The CHARMM version the header states: 24 for NAMD, LAMMPS, VMD and MDAnalysis, which
                             // pretend to be CHARMM 24; CHARMM itself writes its own (35, 36, ...)
    int32_t  marker;         // The size of the Fortran record markers: 4, or 8 (gfortran -frecord-marker=8, g77 on
                             // some 64 bit systems)
    double   delta;          // Integration timestep (AKMA units)
    bool     reverse_endian; // File byte order differs from host
    int64_t  header_size;    // Byte offset to the first frame data
    int32_t* free_indices;   // 1-based indices of free (non-fixed) atoms; NULL when nfixed == 0
} dcd_file_header_t;

// ==================== What the file is ====================

typedef struct dcd_t {

	dcd_file_header_t header;
	size_t num_frames;           // Actual frame count (may differ from header.nset when file was not closed cleanly)

	// Parameters for analytical frame offset computation.
	// offset[0] = file_header_size
	// offset[i] = file_header_size + first_frame_size + (i-1) * frame_size   (i >= 1)
	size_t file_header_size;
	size_t first_frame_size;
	size_t frame_size;

	double* frame_times;          // Precomputed frame timestamps (picoseconds), length = num_frames.
	int64_t* frame_steps;         // Simulation step per frame (istart + i*nsavc), length = num_frames.
	md_unit_t time_unit;          // Time unit for frame_times. Empty when the file carried no timestep and the times are frame ordinals standing in for real time.

    // Full coordinate snapshot of frame 0, used for fixed-atom reconstruction in subsequent frames.
    // Only allocated when nfixed > 0.
    float* first_frame_x;
    float* first_frame_y;
    float* first_frame_z;
} dcd_t;

// ==================== Byte-swap helpers ====================

static inline int32_t dcd_swap32(int32_t v, bool swap) {
    return swap ? (int32_t)MD_BSWAP32((uint32_t)v) : v;
}

static inline void dcd_swap32_array(int32_t* data, int count, bool swap) {
    if (!swap) return;
    for (int i = 0; i < count; ++i) {
        data[i] = MD_BSWAP32(data[i]);
    }
}

// A Fortran record marker of 'marker' bytes (4 or 8)
static inline int64_t dcd_load_marker(const uint8_t* p, int marker, bool rev) {
    if (marker == 8) {
        uint64_t u;
        MD_MEMCPY(&u, p, 8);
        if (rev) u = MD_BSWAP64(u);
        return (int64_t)u;
    }
    uint32_t u;
    MD_MEMCPY(&u, p, 4);
    if (rev) u = MD_BSWAP32(u);
    return (int32_t)u;
}

static bool dcd_read_marker(md_file_t file, int marker, bool rev, int64_t* out) {
    uint8_t buf[8];
    if (md_file_read(file, buf, (size_t)marker) != (size_t)marker) return false;
    *out = dcd_load_marker(buf, marker, rev);
    return true;
}

// ==================== Analytical frame offset ====================

// DCD frames have fixed sizes so the file offset for any frame is a closed-form
// expression.  No offset table is needed.
static inline int64_t dcd_frame_offset(const dcd_t* dcd, int64_t frame_idx) {
    if (frame_idx == 0) return (int64_t)dcd->file_header_size;
    return (int64_t)dcd->file_header_size + (int64_t)dcd->first_frame_size + (frame_idx - 1) * (int64_t)dcd->frame_size;
}

// Compute the byte sizes of the first and subsequent frames from header fields.
// Each coordinate dimension is a Fortran record: [marker: N*4][float[N]][marker: N*4]
// = N * 4 + 2 * marker bytes per dimension.
static void dcd_compute_frame_sizes(const dcd_file_header_t* h,
                                    int64_t* out_first,
                                    int64_t* out_subsequent)
{
    const int64_t m              = h->marker;
    const int64_t extrablocksize = (h->charmm & DCD_HAS_EXTRA_BLOCK) ? (48 + 2 * m) : 0;
    const int64_t ndims          = (h->charmm & DCD_HAS_4DIMS)       ? 4 : 3;
    *out_first      = ((int64_t)h->natoms              * 4 + 2 * m) * ndims + extrablocksize;
    *out_subsequent = ((int64_t)(h->natoms - h->nfixed) * 4 + 2 * m) * ndims + extrablocksize;
}

// ==================== DCD file header parsing ====================

// Parses the full DCD file header starting from the beginning of the file.
// On success, populates *out and leaves the file positioned at the first frame.
// When nfixed > 0, allocates out->free_indices from alloc.
static bool dcd_parse_file_header(md_file_t file, dcd_file_header_t* out, md_allocator_i* alloc) {
    MD_ASSERT(md_file_valid(file) && out && alloc);
    MD_MEMSET(out, 0, sizeof(*out));

    // The very first bytes are a Fortran record marker that must equal 84: 4 bytes of it, or 8 in a file written
    // with 64 bit markers. If they compare as 84 when swapped, the file uses reverse endianness.
    uint8_t first[8];
    if (md_file_read(file, first, 8) != 8) {
        MD_LOG_ERROR("DCD: Failed to read initial record marker");
        return false;
    }
    // A little endian 64 bit marker begins as a 32 bit one does: what follows it (the "CORD" of the header, or the
    // zeros of the marker's upper half) tells them apart
    const bool cord = first[4] == 'C' && first[5] == 'O' && first[6] == 'R' && first[7] == 'D';
    int m = 0;
    if (cord && dcd_load_marker(first, 4, false) == 84)      { m = 4; }
    else if (cord && dcd_load_marker(first, 4, true) == 84)  { m = 4; out->reverse_endian = true; }
    else if (dcd_load_marker(first, 8, false) == 84) { m = 8; }
    else if (dcd_load_marker(first, 8, true) == 84)  { m = 8; out->reverse_endian = true; }
    else {
        MD_LOG_ERROR("DCD: Invalid initial record marker (expected 84)");
        return false;
    }
    out->marker = m;
    md_file_seek(file, m, MD_FILE_BEG);

    const bool rev = out->reverse_endian;

    // Read the 84-byte main header block
    char hdrbuf[84];
    if (md_file_read(file, hdrbuf, 84) != 84) {
        MD_LOG_ERROR("DCD: Failed to read 84-byte header block");
        return false;
    }

    // Trailing record marker must also equal 84
    int64_t trailer = 0;
    if (!dcd_read_marker(file, m, rev, &trailer)) {
        MD_LOG_ERROR("DCD: Failed to read trailing record marker");
        return false;
    }
    if (trailer != 84) {
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
    MD_MEMCPY(&nset,           hdrbuf +  4, 4);
    MD_MEMCPY(&istart,         hdrbuf +  8, 4);
    MD_MEMCPY(&nsavc,          hdrbuf + 12, 4);
    MD_MEMCPY(&namnf,          hdrbuf + 36, 4);
    MD_MEMCPY(&charmm_version, hdrbuf + 80, 4);
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
        MD_MEMCPY(&has_extra, hdrbuf + 44, 4);
        MD_MEMCPY(&has_4dims, hdrbuf + 48, 4);
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
        MD_MEMCPY(&raw, hdrbuf + 40, 4);
        if (rev) raw = MD_BSWAP32(raw);
        float fdelta;
        MD_MEMCPY(&fdelta, &raw, 4);
        delta = (double)fdelta;
    } else {
        uint64_t raw;
        MD_MEMCPY(&raw, hdrbuf + 40, 8);
        if (rev) raw = MD_BSWAP64(raw);
        MD_MEMCPY(&delta, &raw, 8);
    }

    // Title block: [int32: size][size bytes: titles][int32: size]
    // We don't need the title content, so we skip over it entirely.
    int64_t beg_title_size;
    if (!dcd_read_marker(file, m, rev, &beg_title_size)) {
        MD_LOG_ERROR("DCD: Failed to read title block size");
        return false;
    }

    if ((beg_title_size - 4) % 80 == 0) {
        int32_t num_title;
		if (md_file_read(file, &num_title, 4) != 4) {
            MD_LOG_ERROR("DCD: Failed to read the title count");
            return false;
        }
        num_title = dcd_swap32(num_title, rev);
        // The titles are skipped by their count, which is not always what the size of the block says: some
        // writers declare a block of two titles and write three (the trailing size then follows the third)
        if (num_title < 0) {
            MD_LOG_ERROR("DCD: Title block has a negative title count");
            return false;
		}
		// Skip the title lines, each 80 bytes long.
		md_file_seek(file, (int64_t)num_title * 80, MD_FILE_CUR);
	}
	else {
        // Malformed title block
        return false;
    }

	int64_t end_title_size;
    if (!dcd_read_marker(file, m, rev, &end_title_size)) {
        MD_LOG_ERROR("DCD: Failed to read title block trailing size");
        return false;
	}

	if (end_title_size != beg_title_size) {
        MD_LOG_ERROR("DCD: Title block size mismatch");
        return false;
    }


    // Natoms block: [int32: 4][int32: natoms][int32: 4]
    int64_t natoms_sz;
    if (!dcd_read_marker(file, m, rev, &natoms_sz)) {
        MD_LOG_ERROR("DCD: Failed to read natoms block size");
        return false;
    }
    if (natoms_sz != 4) {
        MD_LOG_ERROR("DCD: Natoms record has unexpected block size");
        return false;
    }

    int32_t natoms = 0;
    if (md_file_read(file, &natoms, 4) != 4) {
        MD_LOG_ERROR("DCD: Failed to read natoms");
        return false;
    }
    natoms = dcd_swap32(natoms, rev);

    int64_t natoms_sz2 = 0;
    if (!dcd_read_marker(file, m, rev, &natoms_sz2)) {
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

        int64_t idx_sz;
        if (!dcd_read_marker(file, m, rev, &idx_sz)) {
            MD_LOG_ERROR("DCD: Failed to read free-atom index block size");
            return false;
        }
        if (idx_sz != (int64_t)nfree * 4) {
            MD_LOG_ERROR("DCD: Free-atom index block has unexpected size (%lld, expected %d)", (long long)idx_sz, nfree * 4);
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

        int64_t idx_sz2 = 0;
        if (!dcd_read_marker(file, m, rev, &idx_sz2)) {
            MD_LOG_ERROR("DCD: Failed to read free-atom index block trailing size");
            return false;
        }
    }

    out->nset        = nset;
    out->istart      = istart;
    out->nsavc       = MD_MAX(nsavc, 1);
    out->nfixed      = namnf;
    out->natoms      = natoms;
    out->charmm      = charmm;
    out->version     = charmm_version;
    out->delta       = delta;
    out->header_size = md_file_tell(file);

    return true;
}

// ==================== Unit cell conversion ====================

// The DCD periodic-cell extra block stores 6 doubles, written three ways:
//
//   NAMD > 2.5 (and VMD, MDAnalysis): [A, cos(gamma), B, cos(beta), cos(alpha), C]
//   NAMD <= 2.5:                      [A, gamma, B, beta, alpha, C], angles in degrees
//   CHARMM:                           the symmetric box (shape) matrix, lower triangle [h11, h21, h22, h31, h32, h33],
//                                     whose rows are the cell vectors, in the frame of the coordinates
//
// The header's CHARMM version tells CHARMM's own files (35, 36, ...) from the others, which all state 24; within
// those, the first two are told by whether the angle slots hold cosines, either taken only when it makes a cell
// (angles that span a volume), and the matrix otherwise - what MDAnalysis does, with the difference that a reading
// which makes no cell is not taken.

// Angles in degrees that a cell can have: each strictly between 0 and 180, and the three spanning a volume
static bool dcd_cell_angles_valid(double alpha, double beta, double gamma) {
    if (!(alpha > 0.0 && alpha < 180.0 && beta > 0.0 && beta < 180.0 && gamma > 0.0 && gamma < 180.0)) return false;
    const double ca = cos(MD_DEG_TO_RAD(alpha)), cb = cos(MD_DEG_TO_RAD(beta)), cg = cos(MD_DEG_TO_RAD(gamma));
    return 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg > 1.0e-6;
}

static double dcd_vec_angle(const double u[3], const double v[3]) {
    const double uu = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    const double vv = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    const double c = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (uu * vv);
    return MD_RAD_TO_DEG(acos(MD_CLAMP(c, -1.0, 1.0)));
}

// out_R (may be NULL): the rotation of the frame's positions into the orientation of the cell, p' = R p, and the
// result true when there is one. Only the box matrix has an orientation of its own: CHARMM's coordinates are in the
// frame of its symmetric box, whose vectors are not lower triangular, and they turn with it into an md_unitcell_t's.
// Taken without the turn (as MDAnalysis does), a triclinic CHARMM box held its atoms up to 0.74 of a box length from
// where they are in it.
// The box matrix of CHARMM: its rows (= columns, it is symmetric) are the cell vectors. None when they make none.
static bool dcd_cell_from_box_matrix(md_unitcell_t* out, const double uc[6], double out_R[3][3], bool* out_rotated) {
    const double e1[3] = { uc[0], uc[1], uc[3] };
    const double e2[3] = { uc[1], uc[2], uc[4] };
    const double e3[3] = { uc[3], uc[4], uc[5] };
    const double la = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
    const double lb = sqrt(e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2]);
    const double lc = sqrt(e3[0] * e3[0] + e3[1] * e3[1] + e3[2] * e3[2]);
    if (!(la > 0.0 && lb > 0.0 && lc > 0.0 && isfinite(la + lb + lc))) return false;
    if (!dcd_cell_angles_valid(dcd_vec_angle(e2, e3), dcd_vec_angle(e1, e3), dcd_vec_angle(e1, e2))) return false;
    const double basis[3][3] = { { e1[0], e1[1], e1[2] }, { e2[0], e2[1], e2[2] }, { e3[0], e3[1], e3[2] } };
    md_unitcell_t cell;
    double R[3][3];
    const bool rotated = md_unitcell_from_basis_vectors(&cell, R, basis);
    if (md_unitcell_flags(&cell) == 0) return false;
    *out = cell;
    if (out_R) memcpy(out_R, R, sizeof(R));
    if (out_rotated) *out_rotated = rotated;
    return true;
}

// version is the CHARMM version of the header. A file CHARMM wrote states its own (35, 36, ...) and holds the box
// matrix, which is then read as one first: the readings by values alone cannot tell a matrix of a box sheared by less
// than an Angstrom from cosines. NAMD, LAMMPS, VMD and MDAnalysis write 24 (pretending to be CHARMM 24), and their
// files are read by the values.
static md_unitcell_t dcd_unitcell_from_params(const double uc[6], int32_t version, double out_R[3][3], bool* out_rotated) {
    if (out_R) {
        for (int i = 0; i < 9; ++i) out_R[i / 3][i % 3] = (i % 4 == 0) ? 1.0 : 0.0;
    }
    if (out_rotated) *out_rotated = false;
    md_unitcell_t matrix_cell;
    if (version > 24 && dcd_cell_from_box_matrix(&matrix_cell, uc, out_R, out_rotated)) {
        return matrix_cell;
    }
    const double a = uc[0], b = uc[2], c = uc[5];
    const bool lengths = a > 0.0 && b > 0.0 && c > 0.0 && isfinite(a + b + c);

    if (lengths && fabs(uc[1]) <= 1.0 && fabs(uc[3]) <= 1.0 && fabs(uc[4]) <= 1.0) {
        const double gamma = MD_RAD_TO_DEG(acos(uc[1]));
        const double beta  = MD_RAD_TO_DEG(acos(uc[3]));
        const double alpha = MD_RAD_TO_DEG(acos(uc[4]));
        if (dcd_cell_angles_valid(alpha, beta, gamma)) {
            return md_unitcell_from_extent_and_angles(a, b, c, alpha, beta, gamma);
        }
    }
    if (lengths && dcd_cell_angles_valid(uc[4], uc[3], uc[1])) {
        return md_unitcell_from_extent_and_angles(a, b, c, uc[4], uc[3], uc[1]);
    }

    if (version <= 24 && dcd_cell_from_box_matrix(&matrix_cell, uc, out_R, out_rotated)) {
        return matrix_cell;
    }
    return md_unitcell_none();
}

// ==================== Low-level frame reading ====================

static bool dcd_read_bytes_at(md_file_t file, md_file_offset_t* offset, void* ptr, size_t num_bytes) {
    MD_ASSERT(md_file_valid(file));
    MD_ASSERT(offset);

    if (md_file_read_at(file, *offset, ptr, num_bytes) != num_bytes) {
        return false;
    }

    *offset += (md_file_offset_t)num_bytes;
    return true;
}

// Read one Fortran coordinate record: [int32: N*4][float[N]][int32: N*4].
// When out is NULL the data is skipped.  Returns the atom count on success, -1 on error.
// expected_count is used to validate the record size.
static int dcd_read_coord_record(md_file_t file, md_file_offset_t* offset, int expected_count, bool rev, int marker, float* out) {
    uint8_t mbuf[8];
    if (!dcd_read_bytes_at(file, offset, mbuf, (size_t)marker)) return -1;
    const int64_t sz = dcd_load_marker(mbuf, marker, rev);

    const int n = (int)(sz / 4);
    if (sz != (int64_t)expected_count * 4) {
        MD_LOG_ERROR("DCD: Coordinate record size mismatch (got %d atoms, expected %d)", n, expected_count);
        return -1;
    }

    if (out) {
        if (!dcd_read_bytes_at(file, offset, out, (size_t)sz)) return -1;
        dcd_swap32_array((int32_t*)out, n, rev);
    } else {
        *offset += sz;
    }

    if (!dcd_read_bytes_at(file, offset, mbuf, (size_t)marker)) return -1;

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
                              int natoms, int nfixed, int charmm, int version, bool rev, int marker,
                              const int32_t* free_indices,
                              bool is_first_frame,
                              md_unitcell_t* out_unitcell,
                              float* out_x, float* out_y, float* out_z)
{
    MD_ASSERT(md_file_valid(file));
    MD_ASSERT(offset);
    const bool coords_requested = (out_x && out_y && out_z);
    const int  nfree            = natoms - nfixed;
    const bool all_atoms        = (nfixed == 0 || is_first_frame);
    const int  atoms_in_frame   = all_atoms ? natoms : nfree;

    // --- Optional periodic cell block (CHARMm only) ---
    if ((charmm & DCD_IS_CHARMM) && (charmm & DCD_HAS_EXTRA_BLOCK)) {
        uint8_t mbuf[8];
        if (!dcd_read_bytes_at(file, offset, mbuf, (size_t)marker)) {
            MD_LOG_ERROR("DCD: Failed to read extra block size");
            return false;
        }
        const int64_t block_sz = dcd_load_marker(mbuf, marker, rev);

        if (block_sz == 48) {
            double uc[6];
            if (!dcd_read_bytes_at(file, offset, uc, 48)) {
                MD_LOG_ERROR("DCD: Failed to read unit cell data");
                return false;
            }
            if (rev) {
                for (int i = 0; i < 6; ++i) {
                    uint64_t tmp;
                    MD_MEMCPY(&tmp, &uc[i], 8);
                    tmp = MD_BSWAP64(tmp);
                    MD_MEMCPY(&uc[i], &tmp, 8);
                }
            }
            if (out_unitcell)
                *out_unitcell = dcd_unitcell_from_params(uc, version, NULL, NULL);
        } else {
            // Unknown block size: skip the data
            *offset += block_sz;
        }

        if (!dcd_read_bytes_at(file, offset, mbuf, (size_t)marker)) {
            MD_LOG_ERROR("DCD: Failed to read extra block trailing size");
            return false;
        }
    }

    // --- Coordinate records ---
    if (all_atoms) {
        // First frame or no fixed atoms: read all atom coordinates directly.
        if (dcd_read_coord_record(file, offset, natoms, rev, marker, coords_requested ? out_x : NULL) < 0) {
            MD_LOG_ERROR("DCD: Failed to read X coordinate record");
            return false;
        }
        if (dcd_read_coord_record(file, offset, natoms, rev, marker, coords_requested ? out_y : NULL) < 0) {
            MD_LOG_ERROR("DCD: Failed to read Y coordinate record");
            return false;
        }
        if (dcd_read_coord_record(file, offset, natoms, rev, marker, coords_requested ? out_z : NULL) < 0) {
            MD_LOG_ERROR("DCD: Failed to read Z coordinate record");
            return false;
        }
    } else {
        // Subsequent frame with fixed atoms: read only the nfree mobile coordinates
        // and scatter them into the correct positions in the output buffers.
        // The caller has already filled out_x/y/z with the first-frame values.
        if (coords_requested) {
            bool success = false;
            md_temp_scope_t temp_scope = md_temp_begin();
            float* tmp = (float*)md_temp_alloc(temp_scope, (size_t)nfree * sizeof(float));

            for (int dim = 0; dim < 3; ++dim) {
                float* dst = (dim == 0) ? out_x : (dim == 1) ? out_y : out_z;
                if (dcd_read_coord_record(file, offset, nfree, rev, marker, tmp) < 0) {
                    MD_LOG_ERROR("DCD: Failed to read free-atom coordinate record (dim %d)", dim);
                    goto done;
                }
                for (int i = 0; i < nfree; ++i)
                    dst[free_indices[i] - 1] = tmp[i];
            }

            success = true;

done:
            md_temp_end(temp_scope);
            if (!success) return false;
        } else {
            for (int dim = 0; dim < 3; ++dim) {
                if (dcd_read_coord_record(file, offset, nfree, rev, marker, NULL) < 0) {
                    MD_LOG_ERROR("DCD: Failed to skip coordinate record (dim %d)", dim);
                    return false;
                }
            }
        }
    }

    // --- Optional 4th-dimension record (CHARMm only) ---
    if (charmm & DCD_HAS_4DIMS) {
        if (dcd_read_coord_record(file, offset, atoms_in_frame, rev, marker, NULL) < 0) {
            MD_LOG_ERROR("DCD: Failed to skip 4th-dimension record");
            return false;
        }
    }

    return true;
}

// ==================== Index ====================

// Everything known about the file short of its coordinates: the header, where the frames are and when
// they were taken, the first frame in full (which the fixed atoms of every later frame take their
// positions from). Allocated from alloc.
static bool dcd_index_load(dcd_t* dcd, md_str_t path, md_allocator_i* alloc) {
    MD_MEMSET(dcd, 0, sizeof(*dcd));

    md_file_t file = {0};
    if (!md_file_open(&file, path, MD_FILE_READ)) {
        MD_LOG_ERROR("DCD: Failed to open '" MD_STR_FMT "'", MD_STR_ARG(path));
        return false;
    }
    bool result = false;

    const int64_t filesize = (int64_t)md_file_size(file);
    md_file_seek(file, 0, MD_FILE_BEG);

    // Always parse the raw DCD header so we have ground-truth metadata.
    dcd_file_header_t fhdr;
    if (!dcd_parse_file_header(file, &fhdr, alloc)) {
        MD_LOG_ERROR("DCD: Failed to parse header of '" MD_STR_FMT "'", MD_STR_ARG(path));
        goto done;
    }
    if (fhdr.natoms <= 0) {
        MD_LOG_ERROR("DCD: Atom count is zero or negative");
        goto done;
    }

    // Compute fixed per-frame byte sizes analytically from the header fields.
    int64_t first_frame_size, frame_size;
    dcd_compute_frame_sizes(&fhdr, &first_frame_size, &frame_size);

    // Determine how many complete frames fit in the file payload.
    const int64_t payload = filesize - fhdr.header_size;
    size_t num_frames = 0;
    if (payload >= first_frame_size) {
        num_frames = 1 + (size_t)MD_MAX(0LL, (payload - first_frame_size) / frame_size);
    }

    // If the header declares fewer frames, trust that (the file may be truncated or
    // still being written).
    if (fhdr.nset > 0 && (size_t)(fhdr.nset + 1) < num_frames) {
        num_frames = (size_t)fhdr.nset + 1;
    }
    if (num_frames == 0) {
        MD_LOG_ERROR("DCD: No complete frames found in '" MD_STR_FMT "'", MD_STR_ARG(path));
        goto done;
    }

    // Simulation step per frame. DCD records the first step and the dump interval, so the step of
    // every frame follows without reading any of them.
    int64_t* frame_steps = (int64_t*)md_alloc(alloc, num_frames * sizeof(int64_t));
    for (size_t i = 0; i < num_frames; ++i) {
        frame_steps[i] = (int64_t)fhdr.istart + (int64_t)i * (int64_t)fhdr.nsavc;
    }

    // Frame timestamps: convert from AKMA to picoseconds.
    double* frame_times = (double*)md_alloc(alloc, num_frames * sizeof(double));
    md_unit_t time_unit = md_unit_none();
    if (fhdr.delta > 0) {
        for (size_t i = 0; i < num_frames; ++i) {
            frame_times[i] = (double)frame_steps[i] * fhdr.delta * DCD_AKMA_TO_PS;
        }
        time_unit = md_unit_picosecond();
    } else {
        // No integration timestep in the file, so real time is unknowable. Fall back to the frame
        // ordinal and leave time_unit empty rather than labelling a count as picoseconds.
        for (size_t i = 0; i < num_frames; ++i) {
            frame_times[i] = (double)i;
        }
    }

    dcd->header           = fhdr;
    dcd->num_frames       = num_frames;
    dcd->file_header_size = (size_t)fhdr.header_size;
    dcd->first_frame_size = (size_t)first_frame_size;
    dcd->frame_size       = (size_t)frame_size;
    dcd->frame_times      = frame_times;
    dcd->frame_steps      = frame_steps;
    dcd->time_unit        = time_unit;

    // The full first frame. Later frames with fixed atoms carry only the free ones; the rest keep
    // these positions.
    {
        const int natoms   = fhdr.natoms;
        dcd->first_frame_x = (float*)md_alloc(alloc, (size_t)natoms * sizeof(float));
        dcd->first_frame_y = (float*)md_alloc(alloc, (size_t)natoms * sizeof(float));
        dcd->first_frame_z = (float*)md_alloc(alloc, (size_t)natoms * sizeof(float));
        md_file_offset_t first_frame_offset = dcd_frame_offset(dcd, 0);
        md_unitcell_t unitcell = { 0 };
        // Read the first frame treating nfixed as 0 so that all atoms are read.
        if (!dcd_read_frame_at(file, &first_frame_offset, natoms, 0, fhdr.charmm, fhdr.version, fhdr.reverse_endian, fhdr.marker,
                               NULL, true, &unitcell,
                               dcd->first_frame_x, dcd->first_frame_y, dcd->first_frame_z))
        {
            MD_LOG_ERROR("DCD: Failed to read first frame for fixed-atom initialisation");
            goto done;
        }

        // The coordinates are taken as the file holds them, as every other trajectory's are: a structure
        // of the same simulation (a PDB of NAMD or CHARMM-GUI, centred on the origin like the DCD) lines up
        // with its frames. (A shift by the centre of the cell was made here when the centre of mass of the
        // first frame was nearer the centre of the cell than the origin: which moved atoms that sat in the
        // cell to its corner, and, the centre of mass taken across the periodic boundaries of a box full
        // of atoms being any point in it, was a coin toss for a solvated system.)
    }
    result = true;

done:
    md_file_close(&file);
    return result;
}

// ==================== Run ====================

static const void* dcd_resident(const md_attributes_t* attributes, md_str_t run, const char* leaf, md_attribute_type_t type, size_t* out_count) {
    const md_attribute_t* a = md_attributes_find_in(attributes, run, md_str_from_cstr(leaf));
    if (!a || a->format.type != type || !a->data) return NULL;
    if (out_count) *out_count = md_attribute_value_count(&a->format);
    return a->data;
}

static inline float dcd_load_f32(const uint8_t* p, bool rev) {
    uint32_t u;
    MD_MEMCPY(&u, p, 4);
    if (rev) u = MD_BSWAP32(u);
    float f;
    MD_MEMCPY(&f, &u, 4);
    return f;
}

static inline int32_t dcd_load_i32(const uint8_t* p, bool rev) {
    int32_t v;
    MD_MEMCPY(&v, p, 4);
    return dcd_swap32(v, rev);
}

// source/layout: how a frame of this file is laid out
enum { DCD_LAYOUT_CHARMM, DCD_LAYOUT_REVERSE_ENDIAN, DCD_LAYOUT_NFIXED, DCD_LAYOUT_MARKER, DCD_LAYOUT_VERSION, DCD_LAYOUT_COUNT };

static bool dcd_has_cell_block(int32_t charmm) {
    return (charmm & DCD_IS_CHARMM) && (charmm & DCD_HAS_EXTRA_BLOCK);
}

// <run>/atom/position. One read of the frame; the file stores x, y and z as three planes, which are
// interleaved here into the packed layout the attribute has. Fixed atoms keep their first frame
// position, which the run holds as source/first_frame.
static size_t dcd_position_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    const md_system_t* sys = (const md_system_t*)user_data;
    MD_ASSERT(sys);
    const md_attributes_t* attributes = &sys->attributes;
    if (!slice || slice->num_idx == 0 || slice->num_idx > 2) return 0;

    md_run_source_t src;
    if (!md_run_source(&src, attributes, attr, MD_STR_LIT("atom/position"))) return 0;
    const md_str_t run = src.run;
    const int64_t* offsets = src.offset;
    const int64_t* sizes   = src.size;
    size_t layout_count = 0;
    const int32_t* layout      = dcd_resident(attributes, run, "source/layout",      MD_ATTRIBUTE_TYPE_I32, &layout_count);
    if (!layout || layout_count != DCD_LAYOUT_COUNT) {
        MD_LOG_ERROR("DCD: the run '" MD_STR_FMT "' has lost its source attributes", MD_STR_ARG(run));
        return 0;
    }

    const uint32_t frame = slice->idx[0];
    const size_t N = attr->format.shape[1];
    const bool rev = layout[DCD_LAYOUT_REVERSE_ENDIAN] != 0;
    const size_t nfixed = (size_t)layout[DCD_LAYOUT_NFIXED];
    const bool all_atoms = (nfixed == 0 || frame == 0);
    const size_t n = all_atoms ? N : N - nfixed;

    const int32_t* free_atoms  = NULL;
    const float*   first_frame = NULL;
    if (!all_atoms) {
        size_t free_count = 0, first_count = 0;
        free_atoms  = dcd_resident(attributes, run, "source/free_atoms",  MD_ATTRIBUTE_TYPE_I32, &free_count);
        first_frame = dcd_resident(attributes, run, "source/first_frame", MD_ATTRIBUTE_TYPE_F32, &first_count);
        if (!free_atoms || !first_frame || free_count != n || first_count != N) {
            MD_LOG_ERROR("DCD: the run '" MD_STR_FMT "' has lost its fixed atoms", MD_STR_ARG(run));
            return 0;
        }
    }

    size_t first = 0, count = N;
    if (slice->num_idx == 2) {
        if (slice->idx[1] >= N) return 0;
        first = slice->idx[1];
        count = 1;
    }
    if (cap != count * 3) return 0;

    const size_t frame_size = (size_t)sizes[frame];
    const int    marker = layout[DCD_LAYOUT_MARKER];
    const size_t skip = dcd_has_cell_block(layout[DCD_LAYOUT_CHARMM]) ? 48 + 2 * (size_t)marker : 0;
    const size_t record = n * 4 + 2 * (size_t)marker;
    if (skip + 3 * record > frame_size) return 0;

    md_temp_scope_t temp = md_temp_begin();
    size_t written = 0;
    const md_str_t file_path = src.path;
    // The cell block, if any, is read along: one read of the frame rather than one per plane.
    uint8_t* raw = md_temp_alloc(temp, skip + 3 * record);
    float* xyz = (count == N) ? (float*)dst : md_temp_alloc(temp, N * 3 * sizeof(float));

    if (raw && xyz && md_attribute_io_read_at(io, file_path, offsets[frame], raw, skip + 3 * record) == skip + 3 * record) {
        bool ok = true;
        if (!all_atoms) {
            MD_MEMCPY(xyz, first_frame, N * 3 * sizeof(float));
        }
        for (size_t d = 0; d < 3 && ok; ++d) {
            const uint8_t* rec = raw + skip + d * record;
            if (dcd_load_marker(rec, marker, rev) != (int64_t)(n * 4)) {
                MD_LOG_ERROR("DCD: frame %u of '" MD_STR_FMT "' is not laid out as its header says", frame, MD_STR_ARG(file_path));
                ok = false;
                break;
            }
            const uint8_t* data = rec + marker;
            if (all_atoms) {
                for (size_t i = 0; i < N; ++i) xyz[i * 3 + d] = dcd_load_f32(data + i * 4, rev);
            } else {
                for (size_t i = 0; i < n; ++i) {
                    const int32_t idx = free_atoms[i];
                    if (idx < 0 || (size_t)idx >= N) { ok = false; break; }
                    xyz[(size_t)idx * 3 + d] = dcd_load_f32(data + i * 4, rev);
                }
            }
        }
        if (ok && skip && dcd_load_marker(raw, marker, rev) == 48) {
            // A box matrix not in the orientation of an md_unitcell_t: the positions turn with it
            double uc[6], R[3][3];
            bool rotated = false;
            for (int i = 0; i < 6; ++i) {
                uint64_t u;
                MD_MEMCPY(&u, raw + marker + i * 8, 8);
                if (rev) u = MD_BSWAP64(u);
                MD_MEMCPY(&uc[i], &u, 8);
            }
            dcd_unitcell_from_params(uc, layout[DCD_LAYOUT_VERSION], R, &rotated);
            if (rotated) {
                for (size_t i = 0; i < N; ++i) {
                    float* q = xyz + i * 3;
                    const double x = q[0], y = q[1], z = q[2];
                    q[0] = (float)(R[0][0] * x + R[0][1] * y + R[0][2] * z);
                    q[1] = (float)(R[1][0] * x + R[1][1] * y + R[1][2] * z);
                    q[2] = (float)(R[2][0] * x + R[2][1] * y + R[2][2] * z);
                }
            }
        }
        if (ok) {
            if (count != N) {
                MD_MEMCPY(dst, xyz + first * 3, 3 * sizeof(float));
            }
            written = cap;
        }
    } else {
        MD_LOG_ERROR("DCD: Failed to read frame %u from '" MD_STR_FMT "'", frame, MD_STR_ARG(file_path));
    }
    md_temp_end(temp);
    return written;
}

// <run>/unitcell, for a file whose frames carry a cell block: 48 bytes at the start of the frame
static size_t dcd_cell_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io) {
    const md_system_t* sys = (const md_system_t*)user_data;
    MD_ASSERT(sys);
    const md_attributes_t* attributes = &sys->attributes;
    if (!slice || slice->num_idx == 0 || slice->num_idx > 2) return 0;

    md_run_source_t src;
    if (!md_run_source(&src, attributes, attr, MD_STR_LIT("unitcell"))) return 0;
    const md_str_t run = src.run;
    const int64_t* offsets = src.offset;
    size_t layout_count = 0;
    const int32_t* layout  = dcd_resident(attributes, run, "source/layout", MD_ATTRIBUTE_TYPE_I32, &layout_count);
    if (!layout || layout_count != DCD_LAYOUT_COUNT) {
        MD_LOG_ERROR("DCD: the run '" MD_STR_FMT "' has lost its source attributes", MD_STR_ARG(run));
        return 0;
    }
    if (slice->num_idx == 1 ? cap != 9 : (cap != 3 || slice->idx[1] >= 3)) return 0;

    const bool rev = layout[DCD_LAYOUT_REVERSE_ENDIAN] != 0;
    const int  marker = layout[DCD_LAYOUT_MARKER];
    const md_str_t file_path = src.path;
    uint8_t raw[8 + 48];
    const size_t raw_size = (size_t)marker + 48;
    if (md_attribute_io_read_at(io, file_path, offsets[slice->idx[0]], raw, raw_size) != raw_size) {
        MD_LOG_ERROR("DCD: Failed to read the cell of frame %u from '" MD_STR_FMT "'", slice->idx[0], MD_STR_ARG(file_path));
        return 0;
    }

    float A[3][3] = {0};
    if (dcd_load_marker(raw, marker, rev) == 48) {
        double uc[6];
        for (int i = 0; i < 6; ++i) {
            uint64_t u;
            MD_MEMCPY(&u, raw + marker + i * 8, 8);
            if (rev) u = MD_BSWAP64(u);
            MD_MEMCPY(&uc[i], &u, 8);
        }
        const md_unitcell_t cell = dcd_unitcell_from_params(uc, layout[DCD_LAYOUT_VERSION], NULL, NULL);
        md_unitcell_A_extract_float(A, &cell);
    }
    // A block of another size carries no cell known here: the frame has none.
    if (slice->num_idx == 1) {
        MD_MEMCPY(dst, A, sizeof(A));
    } else {
        MD_MEMCPY(dst, A[slice->idx[1]], sizeof(A[0]));
    }
    return cap;
}

bool md_dcd_system_publish_run(md_system_t* sys, md_str_t filename, md_str_t run, uint32_t flags) {
    MD_ASSERT(sys);
    (void)flags;
    char path_buf[4096];
    const size_t path_len = md_path_write_canonical(path_buf, sizeof(path_buf), filename);
    const md_str_t path = {path_buf, path_len};

    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    bool result = false;

    dcd_t dcd;
    if (path_len == 0 || !dcd_index_load(&dcd, path, arena)) {
        goto done;
    }

    const size_t F = dcd.num_frames;
    const size_t N = (size_t)dcd.header.natoms;
    const size_t nfixed = (size_t)dcd.header.nfixed;

    int64_t* offsets = md_alloc(arena, F * sizeof(int64_t));
    int64_t* sizes   = md_alloc(arena, F * sizeof(int64_t));
    for (size_t i = 0; i < F; ++i) {
        offsets[i] = dcd_frame_offset(&dcd, (int64_t)i);
        sizes[i]   = (int64_t)(i == 0 ? dcd.first_frame_size : dcd.frame_size);
    }

    // The cell: read from each frame when the file keeps one per frame; otherwise the system's own at
    // every frame.
    const md_attribute_virtual_t cell_virt = { .provider = dcd_cell_provider, .user_data = sys };
    float* boxes = NULL;
    if (!dcd_has_cell_block(dcd.header.charmm)) {
        float A[3][3] = {0};
        md_unitcell_A_extract_float(A, &sys->reference.unitcell);
        boxes = md_alloc(arena, F * 9 * sizeof(float));
        for (size_t i = 0; i < F; ++i) {
            MD_MEMCPY(boxes + i * 9, A, sizeof(A));
        }
    }

    const md_attribute_virtual_t pos_virt = { .provider = dcd_position_provider, .user_data = sys };
    const md_run_desc_t desc = {
        .num_frames    = F,
        .num_atoms     = N,
        .time          = dcd.frame_times,
        .time_unit     = dcd.time_unit,   // none without a timestep in the file: time is then ordinals
        .step          = dcd.frame_steps,
        .unitcell      = boxes,
        .unitcell_virt = boxes ? NULL : &cell_virt,
        .source_path   = path,
        .source_offset = offsets,
        .source_size   = sizes,
        .position_virt = &pos_virt,
    };
    if (!md_run_publish(sys, run, &desc)) {
        goto done;
    }

    // How to read a frame of this file: what the providers need beyond the frame table
    md_attributes_t* attributes = &sys->attributes;
    const int32_t layout[DCD_LAYOUT_COUNT] = {
        [DCD_LAYOUT_CHARMM]         = dcd.header.charmm,
        [DCD_LAYOUT_REVERSE_ENDIAN] = dcd.header.reverse_endian ? 1 : 0,
        [DCD_LAYOUT_NFIXED]         = dcd.header.nfixed,
        [DCD_LAYOUT_MARKER]         = dcd.header.marker,
        [DCD_LAYOUT_VERSION]        = dcd.header.version,
    };
    char buf[512];

    bool ok = md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), run, MD_STR_LIT("source/layout")),
        .format = { .type = MD_ATTRIBUTE_TYPE_I32, .components = 1, .rank = 1, .shape = { DCD_LAYOUT_COUNT } },
        .unit = md_unit_none(),
        .description = MD_STR_INIT("CHARMM flags, byte order reversed, number of fixed atoms"),
        .data = layout, .byte_size = sizeof(layout)});


    if (ok && nfixed > 0) {
        // Every frame after the first carries only the free atoms; the fixed ones stay where the
        // first frame put them. Both are small next to the trajectory, so they are held.
        const size_t nfree = N - nfixed;
        int32_t* free_atoms = md_alloc(arena, nfree * sizeof(int32_t));
        for (size_t i = 0; i < nfree; ++i) {
            free_atoms[i] = dcd.header.free_indices[i] - 1;
        }
        float* first_frame = md_alloc(arena, N * 3 * sizeof(float));
        for (size_t i = 0; i < N; ++i) {
            first_frame[i * 3 + 0] = dcd.first_frame_x[i];
            first_frame[i * 3 + 1] = dcd.first_frame_y[i];
            first_frame[i * 3 + 2] = dcd.first_frame_z[i];
        }
        ok = md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), run, MD_STR_LIT("source/free_atoms")),
            .format = { .type = MD_ATTRIBUTE_TYPE_I32, .components = 1, .rank = 1, .shape = { (uint32_t)nfree } },
            .unit = md_unit_none(), .description = MD_STR_INIT("The atoms each frame after the first carries, zero based"),
            .data = free_atoms, .byte_size = nfree * sizeof(int32_t)});
        ok = ok && md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), run, MD_STR_LIT("source/first_frame")),
            .format = { .type = MD_ATTRIBUTE_TYPE_F32, .components = 3, .rank = 1, .shape = { (uint32_t)N } },
            .unit = md_unit_angstrom(), .description = MD_STR_INIT("The first frame as stored"),
            .data = first_frame, .byte_size = N * 3 * sizeof(float)});
    }

    if (!ok) {
        MD_LOG_ERROR("DCD: failed to publish '" MD_STR_FMT "' as '" MD_STR_FMT "'", MD_STR_ARG(path), MD_STR_ARG(run));
        md_attributes_remove_prefix(attributes, run);
        goto done;
    }
    result = true;

done:
    md_arena_allocator_destroy(arena);
    return result;
}

size_t md_dcd_file_atom_count(md_str_t filename) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) return 0;
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_KILOBYTES(64));
    dcd_file_header_t header;
    const size_t n = dcd_parse_file_header(file, &header, arena) && header.natoms > 0 ? (size_t)header.natoms : 0;
    md_arena_allocator_destroy(arena);
    md_file_close(&file);
    return n;
}
