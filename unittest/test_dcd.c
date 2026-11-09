#include "utest.h"
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <md_unitcell.h>

#include <md_dcd.h>
#include <md_system.h>
#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>

// DCD files written by hand, so that what the format allows - a cell per frame, fixed atoms, the
// other byte order - is each exercised without a sample file for every combination.

typedef struct dcd_writer_t {
    md_file_t file;
    bool      rev;
    int       marker;   // the size of the Fortran record markers, 4 or 8
} dcd_writer_t;

static void dcd_put_u32(dcd_writer_t* w, uint32_t u) {
    if (w->rev) u = MD_BSWAP32(u);
    md_file_write(w->file, &u, 4);
}

static void dcd_put_i32(dcd_writer_t* w, int32_t v) { dcd_put_u32(w, (uint32_t)v); }

static void dcd_put_f32(dcd_writer_t* w, float f) {
    uint32_t u;
    MD_MEMCPY(&u, &f, 4);
    dcd_put_u32(w, u);
}

static void dcd_put_f64(dcd_writer_t* w, double d) {
    uint64_t u;
    MD_MEMCPY(&u, &d, 8);
    if (w->rev) u = MD_BSWAP64(u);
    md_file_write(w->file, &u, 8);
}

// A Fortran record marker
static void dcd_put_marker(dcd_writer_t* w, int64_t n) {
    if (w->marker == 8) {
        uint64_t u = (uint64_t)n;
        if (w->rev) u = MD_BSWAP64(u);
        md_file_write(w->file, &u, 8);
    } else {
        dcd_put_i32(w, (int32_t)n);
    }
}

typedef struct dcd_spec_t {
    bool    charmm;       // CHARMM (float delta, flags) rather than X-PLOR (double delta)
    bool    cell_block;   // a unit cell block in every frame; CHARMM only
    int     cell_style;   // how the block holds it: 0 lengths and cosines, 1 lengths and degrees, 2 CHARMM's box matrix,
                          // 3 a box matrix sheared by less than an Angstrom (its off diagonal values pass for cosines)
    int     version;      // the CHARMM version of the header; 0 is 24, what everything but CHARMM writes
    int     marker;       // the size of the record markers; 0 is 4
    bool    rev;          // the other byte order
    int     natoms;
    int     nfixed;       // the last nfixed atoms are fixed
    int     nframes;
    int     istart;
    int     nsavc;
    double  delta;        // AKMA
} dcd_spec_t;

// Atom i of frame f. Near (15, 12, 10), the middle of the cell the CHARMM file declares: where they stay.
static float dcd_coord(int f, int i, int d) {
    const float base[3] = { 15.0f, 12.0f, 10.0f };
    return base[d] + 0.25f * (float)i - 0.5f * (float)d + 0.125f * (float)f;
}

static bool dcd_write(str_t path, const dcd_spec_t* s) {
    dcd_writer_t w = { .rev = s->rev, .marker = s->marker ? s->marker : 4 };
    if (!md_file_open(&w.file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) return false;

    // The 84 byte header record
    uint8_t hdr[84] = {0};
    int32_t ints[21] = {0};
    ints[1] = s->nframes;
    ints[2] = s->istart;
    ints[3] = s->nsavc;
    ints[9] = s->nfixed;
    if (s->charmm) {
        ints[11] = s->cell_block ? 1 : 0;
        ints[20] = s->version ? s->version : 24;
    }
    for (int i = 1; i < 21; ++i) {
        uint32_t u = (uint32_t)ints[i];
        if (s->rev) u = MD_BSWAP32(u);
        MD_MEMCPY(hdr + i * 4, &u, 4);
    }
    MD_MEMCPY(hdr, "CORD", 4);
    if (s->charmm) {
        const float d = (float)s->delta;
        uint32_t u;
        MD_MEMCPY(&u, &d, 4);
        if (s->rev) u = MD_BSWAP32(u);
        MD_MEMCPY(hdr + 40, &u, 4);
    } else {
        uint64_t u;
        MD_MEMCPY(&u, &s->delta, 8);
        if (s->rev) u = MD_BSWAP64(u);
        MD_MEMCPY(hdr + 40, &u, 8);
    }
    dcd_put_marker(&w, 84);
    md_file_write(w.file, hdr, 84);
    dcd_put_marker(&w, 84);

    // One title line
    char title[80];
    MD_MEMSET(title, ' ', sizeof(title));
    MD_MEMCPY(title, "md_unittest", 11);
    dcd_put_marker(&w, 84);
    dcd_put_i32(&w, 1);
    md_file_write(w.file, title, 80);
    dcd_put_marker(&w, 84);

    dcd_put_marker(&w, 4);
    dcd_put_i32(&w, s->natoms);
    dcd_put_marker(&w, 4);

    const int nfree = s->natoms - s->nfixed;
    if (s->nfixed > 0) {
        dcd_put_marker(&w, nfree * 4);
        for (int i = 0; i < nfree; ++i) dcd_put_i32(&w, i + 1);
        dcd_put_marker(&w, nfree * 4);
    }

    for (int f = 0; f < s->nframes; ++f) {
        if (s->charmm && s->cell_block) {
            // A, cos(gamma), B, cos(beta), cos(alpha), C: a cell that grows a little per frame. In degrees as NAMD
            // 2.5 wrote it, or as the box matrix of CHARMM c36: h11, h21, h22, h31, h32, h33
            const double styles[4][6] = {
                { 30.0 + f, 0.0, 24.0, 0.0, 0.0, 20.0 },
                { 30.0 + f, 90.0, 24.0, 90.0, 90.0, 20.0 },
                { 30.0 + f, 5.0, 24.0, 3.0, -2.0, 20.0 },
                { 30.0 + f, 0.5, 24.0, 0.3, -0.2, 20.0 },
            };
            const double* uc = styles[s->cell_style];
            dcd_put_marker(&w, 48);
            for (int i = 0; i < 6; ++i) dcd_put_f64(&w, uc[i]);
            dcd_put_marker(&w, 48);
        }
        const int n = (f == 0 || s->nfixed == 0) ? s->natoms : nfree;
        for (int d = 0; d < 3; ++d) {
            dcd_put_marker(&w, n * 4);
            for (int i = 0; i < n; ++i) {
                // Fixed atoms stay where the first frame put them; a later frame never has them.
                dcd_put_f32(&w, dcd_coord(f, i, d));
            }
            dcd_put_marker(&w, n * 4);
        }
    }
    md_file_close(&w.file);
    return true;
}

static const str_t dcd_paths[] = { STR_INIT("atom/position"), STR_INIT("unitcell") };

// Every frame through one context against what the writer put there: free atoms at their frame's
// position, fixed atoms at their first frame's, and the cell
// the frame's block describes or, without blocks, the system's own.
static void dcd_check_frames(int* utest_result, const md_system_t* sys, str_t run, const dcd_spec_t* spec, const md_unitcell_t* fixed_cell, md_allocator_i* arena) {
    const size_t N = (size_t)spec->natoms;
    const int nfree = spec->natoms - spec->nfixed;
    md_system_state_t got = {.alloc = arena};
    md_system_state_init(&got, N);

    md_system_extract_t* ex = md_system_extract_begin(sys, run, dcd_paths, MD_COUNTOF(dcd_paths), md_get_heap_allocator());
    ASSERT_TRUE(ex != NULL);
    for (int f = 0; f < spec->nframes; ++f) {
        ASSERT_TRUE(md_system_extract_frame(ex, f, &got));
        for (int i = 0; i < spec->natoms; ++i) {
            const int src = (f > 0 && spec->nfixed > 0 && i >= nfree) ? 0 : f;
            EXPECT_EQ(dcd_coord(src, i, 0), got.xyz[i].x);
            EXPECT_EQ(dcd_coord(src, i, 1), got.xyz[i].y);
            EXPECT_EQ(dcd_coord(src, i, 2), got.xyz[i].z);
        }
        if (fixed_cell) {
            EXPECT_EQ(fixed_cell->flags, got.unitcell.flags);
            EXPECT_NEAR(fixed_cell->x, got.unitcell.x, 1.0e-5);
            EXPECT_NEAR(fixed_cell->y, got.unitcell.y, 1.0e-5);
            EXPECT_NEAR(fixed_cell->z, got.unitcell.z, 1.0e-5);
        } else {
            EXPECT_NEAR(30.0 + f, got.unitcell.x, 1.0e-4);
            EXPECT_NEAR(24.0,     got.unitcell.y, 1.0e-4);
            EXPECT_NEAR(20.0,     got.unitcell.z, 1.0e-4);
            EXPECT_NEAR(0.0,      got.unitcell.xy, 1.0e-4);
        }
    }
    md_system_extract_end(ex);
}

#define DCD_RUN STR_LIT("run/dcd")

// A cell in every frame, fixed atoms and the other byte order: everything the reader has to piece
// together, which the run must piece together the same way.
UTEST(dcd, run_matches_the_trajectory) {
    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    const str_t path = STR_INIT("md_unittest_run.dcd");
    const dcd_spec_t spec = { .charmm = true, .cell_block = true, .rev = true, .natoms = 7, .nfixed = 2, .nframes = 5, .istart = 1000, .nsavc = 50, .delta = 0.5 };
    ASSERT_TRUE(dcd_write(path, &spec));

    md_system_t sys = {.alloc = arena};
    ASSERT_TRUE(md_dcd_system_publish_run(&sys, path, DCD_RUN, 0));

    const md_attributes_t* t = &sys.attributes;
    const md_attribute_t* time = md_attributes_find(t, STR_LIT("run/dcd/time"));
    const md_attribute_t* step = md_attributes_find(t, STR_LIT("run/dcd/step"));
    const md_attribute_t* cell = md_attributes_find(t, STR_LIT("run/dcd/unitcell"));
    const md_attribute_t* pos  = md_attributes_find(t, STR_LIT("run/dcd/atom/position"));
    ASSERT_TRUE(time && step && cell && pos);
    EXPECT_EQ((uint32_t)spec.nframes, time->format.shape[0]);
    EXPECT_TRUE(md_unit_equal(time->unit, md_unit_picosecond()));
    // istart + frame * nsavc steps of delta AKMA units
    EXPECT_NEAR(1100 * 0.5 * 0.04888821, ((const double*)time->data)[2], 1.0e-9);
    EXPECT_EQ(1100, ((const int64_t*)step->data)[2]);
    EXPECT_TRUE(md_attributes_find(t, STR_LIT("run/dcd/source/free_atoms")) != NULL);
    // The coordinates as the file holds them: the atoms sit around the middle of the cell, and stay there
    EXPECT_TRUE(md_attributes_find(t, STR_LIT("run/dcd/source/translation")) == NULL);
    dcd_check_frames(utest_result, &sys, DCD_RUN, &spec, NULL, arena);

    // A fixed atom at a later frame, alone: where the first frame put it.
    float xyz[3];
    md_attribute_slice_t one = md_attribute_slice_2(4, 6);
    ASSERT_EQ(3u, md_attribute_extract_f32(xyz, 3, pos, one, md_unit_none()));
    EXPECT_EQ(dcd_coord(0, 6, 0), xyz[0]);

    // The cell of one frame on its own
    float box[9];
    md_attribute_slice_t fr = md_attribute_slice_1(3);
    ASSERT_EQ(9u, md_attribute_extract_f32(box, 9, cell, fr, md_unit_none()));
    EXPECT_NEAR(33.0f, box[0], 1.0e-4f);
    EXPECT_NEAR(20.0f, box[8], 1.0e-4f);

    md_system_free(&sys);
    md_vm_arena_destroy(arena);
    remove(path.ptr);
}

// X-PLOR, no cell in the file and no timestep: the frames are ordinals and every frame has the
// system's own cell, both as the reader has them.
UTEST(dcd, run_without_cell_or_timestep) {
    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    const str_t path = STR_INIT("md_unittest_run_xplor.dcd");
    const dcd_spec_t spec = { .charmm = false, .natoms = 9, .nframes = 3, .istart = 0, .nsavc = 10, .delta = 0.0 };
    ASSERT_TRUE(dcd_write(path, &spec));

    md_system_t sys = {.alloc = arena};
    sys.reference.unitcell = md_unitcell_from_extent(40.0, 30.0, 20.0);
    ASSERT_TRUE(md_dcd_system_publish_run(&sys, path, DCD_RUN, 0));

    const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/dcd/time"));
    ASSERT_TRUE(time != NULL);
    EXPECT_TRUE(md_unit_is_none(time->unit));
    EXPECT_EQ(2.0, ((const double*)time->data)[2]);
    EXPECT_TRUE(md_attributes_find(&sys.attributes, STR_LIT("run/dcd/source/free_atoms")) == NULL);

    dcd_check_frames(utest_result, &sys, DCD_RUN, &spec, &sys.reference.unitcell, arena);

    md_system_state_t st = {.alloc = arena};
    md_system_state_init(&st, spec.natoms);
    const str_t cell_only[] = { STR_INIT("unitcell") };
    md_system_extract_t* ex = md_system_extract_begin(&sys, DCD_RUN, cell_only, 1, md_get_heap_allocator());
    ASSERT_TRUE(ex != NULL);
    ASSERT_TRUE(md_system_extract_frame(ex, 1, &st));
    EXPECT_NEAR(40.0, st.unitcell.x, 1.0e-5);
    EXPECT_NEAR(20.0, st.unitcell.z, 1.0e-5);
    md_system_extract_end(ex);

    md_system_free(&sys);
    md_vm_arena_destroy(arena);
    remove(path.ptr);
}

UTEST(dcd, run_refuses_a_different_system) {
    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    const str_t path = STR_INIT("md_unittest_run_other.dcd");
    const dcd_spec_t spec = { .charmm = true, .natoms = 4, .nframes = 2, .nsavc = 1, .delta = 1.0 };
    ASSERT_TRUE(dcd_write(path, &spec));

    md_system_t sys = {.alloc = arena};
    sys.atom.count = 5;
    EXPECT_FALSE(md_dcd_system_publish_run(&sys, path, DCD_RUN, 0));

    md_vm_arena_destroy(arena);
    remove(path.ptr);
}

UTEST(dcd, atom_count_from_the_header) {
    const str_t path = STR_INIT("md_unittest_count.dcd");
    const dcd_spec_t spec = { .charmm = true, .natoms = 11, .nframes = 2, .istart = 0, .nsavc = 1, .delta = 1.0 };
    ASSERT_TRUE(dcd_write(path, &spec));
    EXPECT_EQ(11u, md_dcd_file_atom_count(path));
    remove(path.ptr);
}

// The three ways a cell block holds the cell, the box matrix of CHARMM c36 among them: its rows are the cell vectors
UTEST(dcd, cell_in_every_convention) {
    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    const str_t path = STR_INIT("md_unittest_cell.dcd");
    for (int style = 0; style < 3; ++style) {
        const dcd_spec_t spec = { .charmm = true, .cell_block = true, .cell_style = style, .natoms = 4, .nframes = 2, .istart = 0, .nsavc = 1, .delta = 1.0 };
        ASSERT_TRUE(dcd_write(path, &spec));
        md_system_t sys = {.alloc = arena};
        ASSERT_TRUE(md_dcd_system_publish_run(&sys, path, DCD_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
        const md_attribute_t* cell = md_attributes_find(&sys.attributes, STR_LIT("run/dcd/unitcell"));
        ASSERT_TRUE(cell != NULL);
        float A[9];
        ASSERT_EQ(9u, md_attribute_extract_f32(A, 9, cell, md_attribute_slice_1(1), md_unit_none()));
        const md_unitcell_t uc = md_unitcell_from_matrix_float((const float(*)[3])A);
        double a, b, c, alpha, beta, gamma;
        md_unitcell_extract_extent_angles(&a, &b, &c, &alpha, &beta, &gamma, &uc);
        if (style < 2) {
            EXPECT_NEAR(31.0, a, 1e-4);
            EXPECT_NEAR(24.0, b, 1e-4);
            EXPECT_NEAR(20.0, c, 1e-4);
            EXPECT_NEAR(90.0, alpha, 1e-3);
            EXPECT_NEAR(90.0, gamma, 1e-3);
        } else {
            // e1 = (31, 5, 3), e2 = (5, 24, -2), e3 = (3, -2, 20)
            const double e1[3] = {31, 5, 3}, e2[3] = {5, 24, -2}, e3[3] = {3, -2, 20};
            const double l1 = sqrt(31.0*31 + 25 + 9), l2 = sqrt(25.0 + 576 + 4), l3 = sqrt(9.0 + 4 + 400);
            EXPECT_NEAR(l1, a, 1e-3);
            EXPECT_NEAR(l2, b, 1e-3);
            EXPECT_NEAR(l3, c, 1e-3);
            EXPECT_NEAR(acos((e2[0]*e3[0] + e2[1]*e3[1] + e2[2]*e3[2]) / (l2 * l3)) * 180.0 / 3.14159265358979, alpha, 1e-3);
            EXPECT_NEAR(acos((e1[0]*e3[0] + e1[1]*e3[1] + e1[2]*e3[2]) / (l1 * l3)) * 180.0 / 3.14159265358979, beta, 1e-3);
            EXPECT_NEAR(acos((e1[0]*e2[0] + e1[1]*e2[1] + e1[2]*e2[2]) / (l1 * l2)) * 180.0 / 3.14159265358979, gamma, 1e-3);
        }
        md_system_free(&sys);
    }
    md_vm_arena_destroy(arena);
    remove(path.ptr);
}

// Fortran record markers of 8 bytes (gfortran -frecord-marker=8, g77 on some 64 bit systems), in both byte orders,
// with a cell in every frame and fixed atoms: read as the same trajectory. They were refused at the first marker.
UTEST(dcd, record_markers_of_8_bytes) {
    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    const str_t path = STR_INIT("md_unittest_marker8.dcd");
    for (int rev = 0; rev < 2; ++rev) {
        const dcd_spec_t spec = { .charmm = true, .cell_block = true, .rev = rev != 0, .marker = 8, .natoms = 7, .nfixed = 2, .nframes = 4, .istart = 0, .nsavc = 10, .delta = 0.5 };
        ASSERT_TRUE(dcd_write(path, &spec));
        EXPECT_EQ(7u, md_dcd_file_atom_count(path));
        md_system_t sys = {.alloc = arena};
        ASSERT_TRUE(md_dcd_system_publish_run(&sys, path, DCD_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
        const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/dcd/time"));
        ASSERT_TRUE(time != NULL);
        EXPECT_EQ(4u, time->format.shape[0]);
        dcd_check_frames(utest_result, &sys, DCD_RUN, &spec, NULL, arena);
        md_system_free(&sys);
    }
    md_vm_arena_destroy(arena);
    remove(path.ptr);
}

// CHARMM's box matrix is the cell in the frame of the coordinates, and its vectors are not in the orientation of an
// md_unitcell_t (a along x, b in the xy plane): the positions turn with it, so that each atom is where it was in the
// cell. They were taken unturned (as MDAnalysis takes them), and an atom of a triclinic CHARMM box was up to 0.74 of a
// box length from its place. A file of CHARMM itself (version 36) is read as a box matrix even when its shear is so
// small that its values pass for cosines, which made a cell of 60 degrees of this box.
UTEST(dcd, charmm_box_matrix_turns_the_positions) {
    md_allocator_i* arena = md_vm_arena_create(MD_GIGABYTES(1));
    const str_t path = STR_INIT("md_unittest_charmm_box.dcd");
    for (int style = 2; style <= 3; ++style) {
        const dcd_spec_t spec = { .charmm = true, .cell_block = true, .cell_style = style, .version = 36, .natoms = 5, .nframes = 2, .istart = 0, .nsavc = 1, .delta = 1.0 };
        ASSERT_TRUE(dcd_write(path, &spec));
        md_system_t sys = {.alloc = arena};
        ASSERT_TRUE(md_dcd_system_publish_run(&sys, path, DCD_RUN, MD_RUN_FLAG_DISABLE_CACHE_WRITE));
        md_system_state_t got = {.alloc = arena};
        md_system_state_init(&got, spec.natoms);
        md_system_extract_t* ex = md_system_extract_begin(&sys, DCD_RUN, dcd_paths, MD_COUNTOF(dcd_paths), md_get_heap_allocator());
        ASSERT_TRUE(ex != NULL);
        for (int f = 0; f < spec.nframes; ++f) {
            ASSERT_TRUE(md_system_extract_frame(ex, f, &got));
            // The box matrix as written, rows the cell vectors, and the fractional coordinates of each atom in it
            const double off[2][3] = { { 5.0, 3.0, -2.0 }, { 0.5, 0.3, -0.2 } };
            const double* o = off[style - 2];
            const double H[3][3] = { { 30.0 + f, o[0], o[1] }, { o[0], 24.0, o[2] }, { o[1], o[2], 20.0 } };
            float Lf[3][3];
            md_unitcell_A_extract_float(Lf, &got.unitcell);
            double L[3][3];
            for (int i = 0; i < 9; ++i) L[i / 3][i % 3] = Lf[i / 3][i % 3];
            // Same lengths: |row i| of L and of H
            for (int r = 0; r < 3; ++r) {
                EXPECT_NEAR(sqrt(H[r][0] * H[r][0] + H[r][1] * H[r][1] + H[r][2] * H[r][2]), sqrt(L[r][0] * L[r][0] + L[r][1] * L[r][1] + L[r][2] * L[r][2]), 1.0e-3);
            }
            for (int i = 0; i < spec.natoms; ++i) {
                // s with x = s H (file) and with x' = s L (read): the same s
                const double x[3]  = { dcd_coord(f, i, 0), dcd_coord(f, i, 1), dcd_coord(f, i, 2) };
                const double xr[3] = { got.xyz[i].x, got.xyz[i].y, got.xyz[i].z };
                // L is lower triangular: solve s L = x' by substitution from the last column
                const double s2 = xr[2] / L[2][2];
                const double s1 = (xr[1] - s2 * L[2][1]) / L[1][1];
                const double s0 = (xr[0] - s1 * L[1][0] - s2 * L[2][0]) / L[0][0];
                // and s H against x
                for (int k = 0; k < 3; ++k) {
                    EXPECT_NEAR(x[k], s0 * H[0][k] + s1 * H[1][k] + s2 * H[2][k], 1.0e-3);
                }
            }
        }
        md_system_extract_end(ex);
        md_system_free(&sys);
    }
    md_vm_arena_destroy(arena);
    remove(path.ptr);
}
