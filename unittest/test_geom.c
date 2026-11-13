#include "utest.h"

#include <md_system.h>
#include <md_run.h>
#include <core/md_geom.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "rmsd.h"
#include "system_util.h"

// md_geom_com over planar inputs: packed first. With indices the source is as long as the
// largest index reaches.
static md_float3_t com_planar(const float* x, const float* y, const float* z, const float* w, const int32_t* idx, size_t n, const md_unitcell_t* cell) {
    size_t len = n;
    if (idx) {
        len = 0;
        for (size_t i = 0; i < n; ++i) len = MD_MAX(len, (size_t)idx[i] + 1);
    }
    md_temp_scope_t temp = md_temp_begin();
    md_float3_t* xyz = md_temp_alloc_array(temp, md_float3_t, MD_ALIGN_TO(len, 16));
    pack_xyz(xyz, x, y, z, len);
    const md_float3_t com = md_geom_com((const float*)xyz, 0, w, 0, idx, n, cell);
    md_temp_end(temp);
    return com;
}

static void unpack_xyz(float* x, float* y, float* z, const md_float3_t* src, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        x[i] = src[i].x;
        y[i] = src[i].y;
        z[i] = src[i].z;
    }
}

UTEST(geom, rmsd) {
    reference_systems_t refsys = reference_systems_load();
    md_allocator_i* alloc = refsys.alloc;
    md_system_t* mol = &refsys.mol_ala;
    ASSERT_TRUE(mol);

    const size_t N = mol->atom.count;
    const size_t cap = MD_ALIGN_TO(N, 16);
    md_float3_t* xyz[2] = {
        md_alloc(alloc, cap * sizeof(md_float3_t)),
        md_alloc(alloc, cap * sizeof(md_float3_t)),
    };
    double* xyz0 = md_alloc(alloc, N * 3 * sizeof(double));
    double* xyz1 = md_alloc(alloc, N * 3 * sizeof(double));

    float* w = md_alloc(alloc, sizeof(float) * N);
	md_atom_extract_masses(w, 0, mol->atom.count, &mol->atom);

    {
        const md_str_t paths[] = { MD_STR_INIT("atom/position") };
        md_system_extract_t* ex = md_system_extract_begin(mol, MD_STR_LIT("run/ala"), paths, 1, md_get_heap_allocator());
        ASSERT_TRUE(ex != NULL);
        ASSERT_TRUE(md_system_extract_frame(ex, 0, &(md_system_state_t){ .xyz = xyz[0] }));
        ASSERT_TRUE(md_system_extract_frame(ex, 1, &(md_system_state_t){ .xyz = xyz[1] }));
        md_system_extract_end(ex);
    }

    for (size_t i = 0; i < N; ++i) {
        xyz0[i * 3 + 0] = xyz[0][i].x;
        xyz0[i * 3 + 1] = xyz[0][i].y;
        xyz0[i * 3 + 2] = xyz[0][i].z;

        xyz1[i * 3 + 0] = xyz[1][i].x;
        xyz1[i * 3 + 1] = xyz[1][i].y;
        xyz1[i * 3 + 2] = xyz[1][i].z;
    }

    // Reference
    double ref_rmsd;
    fast_rmsd((double(*)[3])xyz0, (double(*)[3])xyz1, (int)N, &ref_rmsd);

    // Our implementation
    md_float4_t* xyzw[2] = {
        md_alloc(alloc, N * sizeof(md_float4_t)),
        md_alloc(alloc, N * sizeof(md_float4_t)),
    };
    for (size_t i = 0; i < N; ++i) {
        xyzw[0][i] = md_float4_from_float3(xyz[0][i], w[i]);
        xyzw[1][i] = md_float4_from_float3(xyz[1][i], w[i]);
    }
    const md_float4_t* const cxyzw[2] = { xyzw[0], xyzw[1] };

    md_float3_t com[2] = {
        md_geom_com((const float*)xyz[0], 0, w, 0, 0, N, 0),
        md_geom_com((const float*)xyz[1], 0, w, 0, 0, N, 0),
    };
    double rmsd = md_geom_rmsd(cxyzw, 0, N, com);

    EXPECT_LE(fabs(ref_rmsd - rmsd), 0.1);

    md_free(alloc, xyzw[0], N * sizeof(md_float4_t));
    md_free(alloc, xyzw[1], N * sizeof(md_float4_t));

    md_free(alloc, xyz[0], cap * sizeof(md_float3_t));
    md_free(alloc, xyz[1], cap * sizeof(md_float3_t));
    md_free(alloc, xyz0, N * 3 * sizeof(double));
    md_free(alloc, xyz1, N * 3 * sizeof(double));
    md_free(alloc, w, N * sizeof(float));
    reference_systems_free(&refsys);
}

UTEST(geom, com) {
    /* DISCLAIMER
        Computing the center of mass is a bit tricky, because of the periodic boundary conditions.
        The problem occurs when we have structures which have an extent which covers more than half the box.
        There are many different variations of how to handle this and it seems none of them are perfect.
        Unless you handpick your algorithm based on some external knowledge of the structure you are dealing with.

        In this case, we have settled on the trigonometric approach presented in:
        Bai, Linge, and David Breen. "Calculating center of mass in an unbounded 2D environment." Journal of Graphics Tools 13.4 (2008): 53-60.
    */
    md_float3_t pbc_ext = {5,0,0};
    md_unitcell_t cell = md_unitcell_from_extent(5,0,0);
    {
        const md_float4_t xyzw[] = {
            {1,0,0,1},
            {2,0,0,1},
            {3,0,0,1},
            {4,0,0,1},
        };

        md_float3_t com = md_geom_com((const float*)xyzw, sizeof(md_float4_t), (const float*)xyzw + 3, sizeof(md_float4_t), 0, MD_COUNTOF(xyzw), &cell);
        EXPECT_NEAR(2.5f, com.x, 1.0E-5F);
        EXPECT_EQ(0, com.y);
        EXPECT_EQ(0, com.z);
    }
    
    {
        const md_float4_t xyzw[] = {
            {0,0,0,1},
            {5,0,0,1},
        };

        md_float3_t com = md_geom_com((const float*)xyzw, sizeof(md_float4_t), (const float*)xyzw + 3, sizeof(md_float4_t), 0, MD_COUNTOF(xyzw), &cell);
		com = md_float3_deperiodize_ortho(com, (md_float3_t){ 0,0,0 }, pbc_ext);
        EXPECT_NEAR(0, com.x, 1.0E-5F);
        EXPECT_EQ(0, com.y);
        EXPECT_EQ(0, com.z);
    }

    {
        const md_float4_t xyzw[] = {
            {0,0,0,1},
            {4,0,0,1},
        };

        /*
        Here we expect the 4 to wrap around to -1,
        then added to the 0, producing a center of mass of -0.5.
        which is then placed within the period to 4.5.
        */

        md_float3_t com = md_geom_com((const float*)xyzw, sizeof(md_float4_t), (const float*)xyzw + 3, sizeof(md_float4_t), 0, MD_COUNTOF(xyzw), &cell);
        com = md_float3_deperiodize_ortho(com, md_float3_mul1(pbc_ext, 0.5f), pbc_ext);
        EXPECT_NEAR(4.5f, com.x, 1.0E-5F);
        EXPECT_EQ(0, com.y);
        EXPECT_EQ(0, com.z);
    }

    {
        const md_float3_t pbc_ext = { 5,0,0 };

        const md_float4_t pos0[] = {
            {4,0,0, 1},
            {5,0,0, 1},
            {6,0,0, 1},
            {7,0,0, 1},
        };

        const md_float4_t pos1[] = {
            {4,0,0, 1},
            {0,0,0, 1},
            {1,0,0, 1},
            {2,0,0, 1},
        };

        const md_float4_t pos2[] = {
            {-1,0, 0, 1},
            {0 ,0, 0, 1},
            {1 ,0, 0, 1},
            {2 ,0, 0, 1},
        };

        md_float3_t com0 = md_geom_com((const float*)pos0, sizeof(md_float4_t), (const float*)pos0 + 3, sizeof(md_float4_t), 0, MD_COUNTOF(pos0), &cell);
        md_float3_t com1 = md_geom_com((const float*)pos1, sizeof(md_float4_t), (const float*)pos1 + 3, sizeof(md_float4_t), 0, MD_COUNTOF(pos1), &cell);
        md_float3_t com2 = md_geom_com((const float*)pos2, sizeof(md_float4_t), (const float*)pos2 + 3, sizeof(md_float4_t), 0, MD_COUNTOF(pos2), &cell);

        com0 = md_float3_deperiodize_ortho(com0, (md_float3_t){ 0,0,0 }, pbc_ext);
        com1 = md_float3_deperiodize_ortho(com1, (md_float3_t){ 0,0,0 }, pbc_ext);
        com2 = md_float3_deperiodize_ortho(com2, (md_float3_t){ 0,0,0 }, pbc_ext);
        
        EXPECT_NEAR(0.5f, com0.x, 1.0E-5F);
        EXPECT_NEAR(0.5f, com1.x, 1.0E-5F);
        EXPECT_NEAR(0.5f, com2.x, 1.0E-5F);
    }
}

// A bondless coarse grained chain wrapped into a small box must come out whole. The BB beads sit 3.8A
// apart along x in a 10A box and each carries a side chain bead 3A up in y, across the cell boundary.
UTEST(geom, unwrap_molecule_coarse_grained_without_bonds) {
    md_allocator_i* alloc = md_vm_arena_create(MD_GIGABYTES(1));

    md_system_t sys = {0};
    sys.alloc = alloc;
    const char* res[] = {"BB", "SC1"};
    enum { N = 7 };
    for (int i = 0; i < N; ++i) {
        cg_add_residue(&sys, "LEU", i + 1, MD_RESIDUE_KIND_AMINO_ACID, res, 2, alloc);
    }
    ASSERT_TRUE(test_infer_molecules(&sys));
    ASSERT_EQ(md_system_molecule_count(&sys), 1u);

    float x[2 * N], y[2 * N], z[2 * N];
    for (int i = 0; i < N; ++i) {
        x[2 * i]     = fmodf(1.0f + 3.8f * i, 10.0f);
        y[2 * i]     = 8.5f;
        z[2 * i]     = 5.0f;
        x[2 * i + 1] = x[2 * i];
        y[2 * i + 1] = fmodf(8.5f + 3.0f, 10.0f); // 3A up, wrapped to the bottom of the cell
        z[2 * i + 1] = 5.0f;
    }
    md_float3_t xyz[2 * N];
    md_system_state_t state = { .num_atoms = 2 * N, .xyz = pack_xyz(xyz, x, y, z, 2 * N), .unitcell = md_unitcell_from_extent(10, 10, 10) };

    md_system_unwrap(&state, &sys);
    unpack_xyz(x, y, z, xyz, state.num_atoms);

    for (int i = 0; i + 1 < N; ++i) {
        EXPECT_NEAR(x[2 * (i + 1)] - x[2 * i], 3.8f, 1.0e-4f);
    }
    for (int i = 0; i < N; ++i) {
        EXPECT_NEAR(x[2 * i + 1] - x[2 * i], 0.0f, 1.0e-4f);
        EXPECT_NEAR(y[2 * i + 1] - y[2 * i], 3.0f, 1.0e-4f);
    }

    md_vm_arena_destroy(alloc);
}

// A chain of nine atoms spaced 2A along x inside a 10A box, so the wrapped input folds back on
// itself twice. md_geom_unwrap must recover a straight line.
//
// @NOTE: this is the regression test for two bugs that were live together here. min_image_ortho was
// being fed a half extent where it wants a reciprocal one, and the result was written as x + dx
// instead of parent + dx, which is wrong for every non root atom even with no image correction at
// all. Either one alone destroys the spacing checked below.
UTEST(geom, unwrap_molecule_ortho) {
    md_allocator_i* alloc = md_get_heap_allocator();

    md_system_t sys = {0};
    make_chain_system(&sys, alloc, 1, 9);
    ASSERT_TRUE(test_infer_molecules(&sys));

    float x[9], y[9] = {0}, z[9] = {0};
    for (int i = 0; i < 9; ++i) {
        x[i] = fmodf(1.0f + 2.0f * i, 10.0f);
    }
    md_float3_t xyz[9];
    md_system_state_t state = { .num_atoms = 9, .xyz = pack_xyz(xyz, x, y, z, 9), .unitcell = md_unitcell_from_extent(10, 10, 10) };

    const md_molecule_t mol = md_system_molecule(&sys, 0);
    ASSERT_LT(0u, (unsigned)mol.count);
    md_geom_unwrap((float*)state.xyz, sizeof(md_float3_t), mol.atom_idx, mol.parent_idx, mol.count, &state.unitcell);
    unpack_xyz(x, y, z, xyz, state.num_atoms);

    for (int i = 0; i < 8; ++i) {
        EXPECT_NEAR(x[i + 1] - x[i], 2.0f, 1.0e-4f);
        EXPECT_NEAR(y[i], 0.0f, 1.0e-4f);
        EXPECT_NEAR(z[i], 0.0f, 1.0e-4f);
    }

    md_system_free(&sys);
}

// The same chain folded back by whole lattice vectors of a sheared cell. Unwrapping must recover the
// original straight line up to a single rigid lattice translation of the whole structure.
UTEST(geom, unwrap_molecule_triclinic) {
    md_allocator_i* alloc = md_get_heap_allocator();

    md_system_t sys = {0};
    make_chain_system(&sys, alloc, 1, 9);
    ASSERT_TRUE(test_infer_molecules(&sys));

    md_unitcell_t cell = md_unitcell_from_basis_parameters(10, 10, 10, 2, 1, 3);
    float A[3][3] = {0};
    md_unitcell_A_extract_float(A, &cell);

    float tx[9], ty[9], tz[9];
    float x[9], y[9], z[9];
    for (int i = 0; i < 9; ++i) {
        tx[i] = 1.0f + 2.00f * i;
        ty[i] = 1.0f + 0.50f * i;
        tz[i] = 1.0f + 0.25f * i;

        const int na = (i >= 5) ? -1 : 0;   // displace the tail by one lattice vector along a
        const int nb = (i >= 7) ? -1 : 0;   // and the last two also along b

        x[i] = tx[i] + na * A[0][0] + nb * A[1][0];
        y[i] = ty[i] + na * A[0][1] + nb * A[1][1];
        z[i] = tz[i] + na * A[0][2] + nb * A[1][2];
    }
    md_float3_t xyz[9];
    md_system_state_t state = { .num_atoms = 9, .xyz = pack_xyz(xyz, x, y, z, 9), .unitcell = cell };

    const md_molecule_t mol = md_system_molecule(&sys, 0);
    ASSERT_LT(0u, (unsigned)mol.count);
    md_geom_unwrap((float*)state.xyz, sizeof(md_float3_t), mol.atom_idx, mol.parent_idx, mol.count, &state.unitcell);
    unpack_xyz(x, y, z, xyz, state.num_atoms);

    const float ox = x[0] - tx[0];
    const float oy = y[0] - ty[0];
    const float oz = z[0] - tz[0];

    for (int i = 0; i < 9; ++i) {
        EXPECT_NEAR(x[i], tx[i] + ox, 1.0e-3f);
        EXPECT_NEAR(y[i], ty[i] + oy, 1.0e-3f);
        EXPECT_NEAR(z[i], tz[i] + oz, 1.0e-3f);
    }

    md_system_free(&sys);
}

// md_geom_min_image, packed and of md_float4_t, reduces a separation vector to its shortest periodic image, so no
// residue may exceed half the box and the result must differ from the input by whole box lengths.
UTEST(geom, min_image) {
    const float e = 10.0f;
    md_unitcell_t cell = md_unitcell_from_extent(e, e, e);

    md_float3_t dx3[4] = {
        { 1.0f, -2.0f,  3.0f},      // already minimal, must be left alone
        {12.0f,  0.0f,  0.0f},
        {-27.0f, 41.0f, -8.0f},
        { 4.9f, -4.9f,  4.9f},      // just inside the half box
    };
    const md_float3_t in3[4] = { dx3[0], dx3[1], dx3[2], dx3[3] };

    md_geom_min_image((float*)dx3, 0, 4, &cell);

    for (int i = 0; i < 4; ++i) {
        for (int c = 0; c < 3; ++c) {
            EXPECT_LE(fabsf(dx3[i].elem[c]), 0.5f * e + 1.0e-4f);
            const float shift = in3[i].elem[c] - dx3[i].elem[c];
            EXPECT_NEAR(shift, e * nearbyintf(shift / e), 1.0e-3f);
        }
    }

    EXPECT_NEAR(dx3[0].x,  1.0f, 1.0e-4f);
    EXPECT_NEAR(dx3[0].y, -2.0f, 1.0e-4f);
    EXPECT_NEAR(dx3[0].z,  3.0f, 1.0e-4f);
    EXPECT_NEAR(dx3[1].x,  2.0f, 1.0e-4f);
    EXPECT_NEAR(dx3[3].x,  4.9f, 1.0e-4f);

    md_float4_t dx4[2] = { {12.0f, 0.0f, 0.0f, 7.0f}, {-27.0f, 41.0f, -8.0f, 8.0f} };
    md_geom_min_image((float*)dx4, sizeof(md_float4_t), 2, &cell);
    EXPECT_NEAR(dx4[0].x, 2.0f, 1.0e-4f);
    EXPECT_NEAR(dx4[1].x, 3.0f, 1.0e-4f);
    EXPECT_NEAR(dx4[1].y, 1.0f, 1.0e-4f);
    EXPECT_NEAR(dx4[1].z, 2.0f, 1.0e-4f);
    EXPECT_NEAR(dx4[0].w, 7.0f, 1.0e-4f);   // w must be left untouched
    EXPECT_NEAR(dx4[1].w, 8.0f, 1.0e-4f);
}

// ---- periodic image selection -------------------------------------------------------------

#define PBC_NP 12

// A lumpy, chiral-ish shape so no two points are interchangeable and the fit is well determined.
static void pbc_make_shape(md_float4_t out[PBC_NP], float radius) {
    for (int i = 0; i < PBC_NP; ++i) {
        const float t = (float)i / PBC_NP * 6.2831853f;
        const float u = (float)(i % 5) / 5.0f;
        out[i] = md_float4_set(radius * cosf(t) * (0.4f + u),
                          radius * sinf(t) * (0.5f + 0.5f * u),
                          radius * (u - 0.5f),
                          1.0f + 0.05f * i);
    }
}

static md_float3x3_t pbc_rot(float angle) {
    const float c = cosf(angle), s = sinf(angle);
    const md_float3_t k = md_float3_normalize(md_float3_set(0.3f, 0.6f, 0.74f));
    md_float3x3_t m;
    m.elem[0][0]=c+k.x*k.x*(1-c);     m.elem[1][0]=k.x*k.y*(1-c)-k.z*s; m.elem[2][0]=k.x*k.z*(1-c)+k.y*s;
    m.elem[0][1]=k.y*k.x*(1-c)+k.z*s; m.elem[1][1]=c+k.y*k.y*(1-c);     m.elem[2][1]=k.y*k.z*(1-c)-k.x*s;
    m.elem[0][2]=k.z*k.x*(1-c)-k.y*s; m.elem[1][2]=k.z*k.y*(1-c)+k.x*s; m.elem[2][2]=c+k.z*k.z*(1-c);
    return m;
}

static float pbc_wrap(float v, float L) { v = fmodf(v, L); return v < 0.0f ? v + L : v; }

// Angle between R and the inverse of R_true, in degrees. Zero when R undoes R_true.
static float pbc_angle_err(md_float3x3_t R, md_float3x3_t R_true) {
    const md_float3x3_t E = md_float3x3_mul(R, R_true);
    const float tr = E.elem[0][0] + E.elem[1][1] + E.elem[2][2];
    return acosf(MD_CLAMP((tr - 1.0f) * 0.5f, -1.0f, 1.0f)) * 180.0f / 3.14159265f;
}

// A set broken across a boundary must come back whole, with a centre that does not depend on which
// images the input happened to arrive in.
UTEST(geom, pbc_deperiodize_self) {
    const float L = 20.0f;
    md_unitcell_t cell = md_unitcell_from_extent(L, L, L);

    md_float4_t shape[PBC_NP];
    pbc_make_shape(shape, 3.0f);

    // Place it straddling the origin corner, then wrap every point independently
    md_float4_t pts[PBC_NP];
    for (int i = 0; i < PBC_NP; ++i) {
        const md_float3_t v = md_float3_from_float4(shape[i]);
        pts[i] = md_float4_set(pbc_wrap(v.x, L), pbc_wrap(v.y, L), pbc_wrap(v.z, L), shape[i].w);
    }

    md_float3_t com = {0};
    md_geom_deperiodize_self(pts, PBC_NP, &cell, &com);

    // Every pairwise separation must match the original shape: the set is whole again
    for (int i = 0; i < PBC_NP; ++i) {
        for (int j = i + 1; j < PBC_NP; ++j) {
            const float got = md_float3_distance(md_float3_from_float4(pts[i]), md_float3_from_float4(pts[j]));
            const float ref = md_float3_distance(md_float3_from_float4(shape[i]), md_float3_from_float4(shape[j]));
            EXPECT_NEAR(got, ref, 1.0e-3f);
        }
    }

    // The centre must sit inside the reassembled set, not at some average of scattered images
    const md_float3_t shape_com = md_geom_com((const float*)shape, sizeof(md_float4_t), (const float*)shape + 3, sizeof(md_float4_t), 0, PBC_NP, 0);
    for (int i = 0; i < PBC_NP; ++i) {
        const md_float3_t d = md_float3_sub(md_float3_from_float4(pts[i]), com);
        const md_float3_t r = md_float3_sub(md_float3_from_float4(shape[i]), shape_com);
        EXPECT_NEAR(d.x, r.x, 1.0e-3f);
        EXPECT_NEAR(d.y, r.y, 1.0e-3f);
        EXPECT_NEAR(d.z, r.z, 1.0e-3f);
    }

    // Already consistent input must be left alone
    md_float4_t again[PBC_NP];
    md_memcpy(again, pts, sizeof(pts));
    md_float3_t com2 = {0};
    md_geom_deperiodize_self(again, PBC_NP, &cell, &com2);
    for (int i = 0; i < PBC_NP; ++i) {
        EXPECT_NEAR(again[i].x, pts[i].x, 1.0e-4f);
        EXPECT_NEAR(again[i].w, pts[i].w, 1.0e-6f);  // weights carried through untouched
    }
}

// A set that arrives whole, in an image other than the reference one, must be left in THAT image.
//
// The alternation seeds from the circular mean, which always lands in the reference cell, so
// without an explicit correction the whole set is quietly translated into image (0,0,0). That is
// invisible to a caller reading relative geometry and a whole cell vector wrong for one that
// combines the reported centre with coordinates it never passed in - which is how viamd's recenter
// came to place a structure in the middle of the WRONG image on nojump trajectories.
UTEST(geom, pbc_deperiodize_self_keeps_image) {
    const float L = 20.0f;
    md_unitcell_t cell = md_unitcell_from_extent(L, L, L);

    md_float4_t shape[PBC_NP];
    pbc_make_shape(shape, 3.0f);
    const md_float3_t shape_com = md_geom_com((const float*)shape, sizeof(md_float4_t), (const float*)shape + 3, sizeof(md_float4_t), 0, PBC_NP, 0);

    for (int nx = -2; nx <= 2; ++nx) {
        for (int nz = -1; nz <= 1; ++nz) {
            const md_float3_t off = md_float3_set(nx * L, 10.0f, 10.0f + nz * L);

            md_float4_t pts[PBC_NP];
            for (int i = 0; i < PBC_NP; ++i) {
                pts[i] = md_float4_from_float3(md_float3_add(md_float3_from_float4(shape[i]), off), shape[i].w);
            }

            md_float3_t com = {0};
            md_geom_deperiodize_self(pts, PBC_NP, &cell, &com);

            // Nothing was broken to begin with, so nothing may move
            for (int i = 0; i < PBC_NP; ++i) {
                const md_float3_t want = md_float3_add(md_float3_from_float4(shape[i]), off);
                EXPECT_NEAR(pts[i].x, want.x, 1.0e-3f);
                EXPECT_NEAR(pts[i].y, want.y, 1.0e-3f);
                EXPECT_NEAR(pts[i].z, want.z, 1.0e-3f);
            }

            // and the centre must be reported in the image the points actually occupy
            const md_float3_t want_com = md_float3_add(shape_com, off);
            EXPECT_NEAR(com.x, want_com.x, 1.0e-3f);
            EXPECT_NEAR(com.y, want_com.y, 1.0e-3f);
            EXPECT_NEAR(com.z, want_com.z, 1.0e-3f);
        }
    }
}

// Same contract for the joint fit: out_com and the placed points come back in the target's image,
// so a caller can build a transform against coordinates it did not hand over. With a rotation in
// play the old behaviour was not merely one cell vector off - it was off by R times a cell vector,
// which is not a periodic image at all.
UTEST(geom, pbc_optimal_rotation_keeps_image) {
    const float L = 20.0f;
    md_unitcell_t cell = md_unitcell_from_extent(L, L, L);

    md_float4_t ref[PBC_NP];
    pbc_make_shape(ref, 3.0f);
    const md_float3_t ref_com = md_geom_com((const float*)ref, sizeof(md_float4_t), (const float*)ref + 3, sizeof(md_float4_t), 0, PBC_NP, 0);

    const md_float3x3_t R_true = pbc_rot(40.0f * 3.14159265f / 180.0f);

    for (int nx = -1; nx <= 2; ++nx) {
        const md_float3_t c_true = md_float3_set(9.0f + nx * L, 11.0f, 8.0f - nx * L);

        md_float4_t trg[PBC_NP];
        for (int i = 0; i < PBC_NP; ++i) {
            const md_float3_t v = md_float3_add(md_float3x3_mul_float3(R_true, md_float3_sub(md_float3_from_float4(ref[i]), ref_com)), c_true);
            trg[i] = md_float4_from_float3(v, ref[i].w);
        }

        md_float3x3_t R = {0};
        md_float3_t com = {0};
        md_float4_t placed[PBC_NP];
        const float residual = md_geom_optimal_rotation_pbc(&R, &com, placed, ref, ref_com, trg, PBC_NP, &cell, 8, 1.0e-6f);

        EXPECT_LT(residual, 1.0e-2f);
        EXPECT_NEAR(com.x, c_true.x, 1.0e-2f);
        EXPECT_NEAR(com.y, c_true.y, 1.0e-2f);
        EXPECT_NEAR(com.z, c_true.z, 1.0e-2f);

        // Untouched input, so the placed points must be exactly where they came in
        for (int i = 0; i < PBC_NP; ++i) {
            EXPECT_NEAR(placed[i].x, trg[i].x, 1.0e-2f);
            EXPECT_NEAR(placed[i].y, trg[i].y, 1.0e-2f);
            EXPECT_NEAR(placed[i].z, trg[i].z, 1.0e-2f);
        }

        // The recentering transform viamd builds must land the set on the cell centre, in image 0.
        const md_float3_t centre = md_float3_set(0.5f * L, 0.5f * L, 0.5f * L);
        md_float3_t sum = md_float3_zero();
        float  wsum = 0.0f;
        for (int i = 0; i < PBC_NP; ++i) {
            const md_float3_t q = md_float3_add(centre, md_float3x3_mul_float3(R, md_float3_sub(md_float3_from_float4(trg[i]), com)));
            sum = md_float3_add(sum, md_float3_mul1(q, trg[i].w));
            wsum += trg[i].w;
        }
        sum = md_float3_div1(sum, wsum);
        EXPECT_NEAR(sum.x, centre.x, 1.0e-2f);
        EXPECT_NEAR(sum.y, centre.y, 1.0e-2f);
        EXPECT_NEAR(sum.z, centre.z, 1.0e-2f);
    }
}

// The image the set is left in has to survive a basis with off diagonal terms, where the lattice
// vector being undone is not axis aligned.
UTEST(geom, pbc_deperiodize_self_keeps_image_triclinic) {
    const float a = 47.8497f;
    const float c = 33.8349f;
    md_unitcell_t cell = md_unitcell_from_basis_parameters(a, a, c, 0.0, a * 0.5f, a * 0.5f);
    ASSERT_TRUE(md_unitcell_is_triclinic(&cell));

    md_float3x3_t A = {0};
    md_unitcell_A_extract_float(A.elem, &cell);

    md_float4_t shape[PBC_NP];
    pbc_make_shape(shape, 5.0f);
    const md_float3_t shape_com = md_geom_com((const float*)shape, sizeof(md_float4_t), (const float*)shape + 3, sizeof(md_float4_t), 0, PBC_NP, 0);
    const md_float3_t base = md_float3x3_mul_float3(A, md_float3_set(0.5f, 0.5f, 0.5f));

    for (int nx = -1; nx <= 1; ++nx) {
        for (int ny = -1; ny <= 1; ++ny) {
            const md_float3_t off = md_float3_add(base, md_float3x3_mul_float3(A, md_float3_set((float)nx, (float)ny, 0.0f)));

            md_float4_t pts[PBC_NP];
            for (int i = 0; i < PBC_NP; ++i) {
                pts[i] = md_float4_from_float3(md_float3_add(md_float3_from_float4(shape[i]), off), shape[i].w);
            }

            md_float3_t com = {0};
            md_geom_deperiodize_self(pts, PBC_NP, &cell, &com);

            const md_float3_t want_com = md_float3_add(shape_com, off);
            EXPECT_NEAR(com.x, want_com.x, 1.0e-2f);
            EXPECT_NEAR(com.y, want_com.y, 1.0e-2f);
            EXPECT_NEAR(com.z, want_com.z, 1.0e-2f);
        }
    }
}

// Every other test in this block uses a cube, which hides anything that only goes wrong when the
// basis has off diagonal terms. This one uses a rhombic dodecahedron - the shape GROMACS writes for
// -bt dodecahedron, and the shape that exposed the bug this test now guards.
//
// The circular mean has to be taken in FRACTIONAL space and carried back through the full basis.
// Handling each Cartesian axis on its own is only valid when the basis is diagonal; do it on a
// triclinic cell and the centre comes back displaced by a lattice vector, which then drags the
// whole deperiodize / align chain into the wrong periodic image.
UTEST(geom, pbc_com_triclinic) {
    // a = b, c = a/sqrt(2), third vector leaning by half a cell in x and y
    const float a = 47.8497f;
    const float c = 33.8349f;
    md_unitcell_t cell = md_unitcell_from_basis_parameters(a, a, c, 0.0, a * 0.5f, a * 0.5f);
    ASSERT_TRUE(md_unitcell_is_triclinic(&cell));

    md_float3x3_t A = {0};
    md_unitcell_A_extract_float(A.elem, &cell);
    md_float3x3_t I = {0};
    md_unitcell_I_extract_float(I.elem, &cell);

    // A compact blob sitting well inside the cell, nowhere near a boundary.
    md_float4_t pts[PBC_NP];
    const md_float3_t centre = md_float3x3_mul_float3(A, md_float3_set(0.41f, 0.63f, 0.48f));
    for (int i = 0; i < PBC_NP; ++i) {
        const float t = (float)i / PBC_NP * 6.2831853f;
        pts[i] = md_float4_set(centre.x + 6.0f * cosf(t), centre.y + 6.0f * sinf(t), centre.z + 3.0f * cosf(t * 2.0f), 1.0f + (float)(i % 3));
    }

    // Nothing is wrapped, so the periodic mean must agree with the plain mean. It is the same set
    // of points either way.
    const md_float3_t com_plain = md_geom_com((const float*)pts, sizeof(md_float4_t), (const float*)pts + 3, sizeof(md_float4_t), 0, PBC_NP, 0);
    const md_float3_t com_pbc   = md_geom_com((const float*)pts, sizeof(md_float4_t), (const float*)pts + 3, sizeof(md_float4_t), 0, PBC_NP, &cell);
    EXPECT_NEAR(com_pbc.x, com_plain.x, 0.2f);
    EXPECT_NEAR(com_pbc.y, com_plain.y, 0.2f);
    EXPECT_NEAR(com_pbc.z, com_plain.z, 0.2f);

    // Same blob, every point folded into the primary cell the way a trajectory stores it. The
    // recovered centre must be the same one, in the same image - not a lattice vector away.
    md_float4_t wrapped[PBC_NP];
    for (int i = 0; i < PBC_NP; ++i) {
        md_float3_t f = md_float3x3_mul_float3(I, md_float3_from_float4(pts[i]));
        f.x -= floorf(f.x);
        f.y -= floorf(f.y);
        f.z -= floorf(f.z);
        wrapped[i] = md_float4_from_float3(md_float3x3_mul_float3(A, f), pts[i].w);
    }

    md_float3_t com_wrapped = {0};
    md_geom_deperiodize_self(wrapped, PBC_NP, &cell, &com_wrapped);
    EXPECT_NEAR(com_wrapped.x, com_plain.x, 0.05f);
    EXPECT_NEAR(com_wrapped.y, com_plain.y, 0.05f);
    EXPECT_NEAR(com_wrapped.z, com_plain.z, 0.05f);

    // and the blob itself must be back in one piece, in that same image
    for (int i = 0; i < PBC_NP; ++i) {
        EXPECT_NEAR(wrapped[i].x, pts[i].x, 0.05f);
        EXPECT_NEAR(wrapped[i].y, pts[i].y, 0.05f);
        EXPECT_NEAR(wrapped[i].z, pts[i].z, 0.05f);
    }
}

// ------------------------------------------------------------------------------------------------
// Periodic centre of mass.
//
// The circular mean has to be taken in FRACTIONAL space and carried back out through the full
// basis. Resolving each Cartesian axis against itself alone is only valid when the basis is
// diagonal, so an ORTHORHOMBIC cell hides that mistake completely and a triclinic one does not.
// Everything below therefore runs against both a cube and a rhombic dodecahedron - the shape
// GROMACS writes for -bt dodecahedron, and the shape that surfaced this.
//
// md_geom_com (packed) dispatches to one of four internal variants depending on whether weights and
// indices are present, and each of those has an AVX512, an AVX2, an SSE2 and a scalar remainder
// path. The counts below straddle the block sizes (4 / 8 / 16) so that on any given build both the
// vector body and the remainder loop are exercised, and every variant is called for each count.

#define COM_MAX_PTS 253

static md_unitcell_t com_cell_ortho(void) {
    return md_unitcell_from_extent(47.8497, 47.8497, 33.8349);
}

static md_unitcell_t com_cell_triclinic(void) {
    const double a = 47.8497;   // a = b, c = a/sqrt(2), third vector leaning half a cell in x and y
    const double c = 33.8349;
    return md_unitcell_from_basis_parameters(a, a, c, 0.0, a * 0.5, a * 0.5);
}

// A compact blob sitting well inside the cell, nowhere near a boundary, so the periodic mean has an
// unambiguous answer: the ordinary arithmetic mean of the very same points.
static void com_make_blob(float* x, float* y, float* z, float* w, size_t count, const md_unitcell_t* cell) {
    md_float3x3_t A = {0};
    md_unitcell_A_extract_float(A.elem, cell);
    const md_float3_t centre = md_float3x3_mul_float3(A, md_float3_set(0.41f, 0.63f, 0.48f));
    for (size_t i = 0; i < count; ++i) {
        const float t = (float)i * 0.7913f;
        x[i] = centre.x + 3.0f * cosf(t);
        y[i] = centre.y + 3.0f * sinf(t * 1.3f);
        z[i] = centre.z + 2.0f * cosf(t * 0.7f);
        w[i] = 1.0f + (float)(i % 7);
    }
}

static const size_t COM_COUNTS[] = { 1, 2, 3, 5, 7, 8, 15, 16, 17, 33, 61, 64, 253 };

// With nothing wrapped, the periodic mean is just the mean. All four variants have to say so, for
// both cell shapes, at every count - and strided md_float4_t points have to agree with packed ones.
UTEST(geom, com_pbc_matches_plain_mean) {
    md_unitcell_t cells[2];
    cells[0] = com_cell_ortho();
    cells[1] = com_cell_triclinic();
    ASSERT_TRUE(md_unitcell_is_orthorhombic(&cells[0]));
    ASSERT_TRUE(md_unitcell_is_triclinic(&cells[1]));

    for (int ci = 0; ci < 2; ++ci) {
        for (size_t k = 0; k < MD_COUNTOF(COM_COUNTS); ++k) {
            const size_t n = COM_COUNTS[k];
            float x[COM_MAX_PTS], y[COM_MAX_PTS], z[COM_MAX_PTS], w[COM_MAX_PTS];
            int32_t idx[COM_MAX_PTS];
            md_float4_t xyzw[COM_MAX_PTS];
            com_make_blob(x, y, z, w, n, &cells[ci]);
            for (size_t i = 0; i < n; ++i) {
                idx[i]  = (int32_t)i;
                xyzw[i] = md_float4_set(x[i], y[i], z[i], w[i]);
            }

            const md_float3_t plain_u = com_planar(x, y, z, NULL, NULL, n, NULL);
            const md_float3_t plain_w = com_planar(x, y, z, w,    NULL, n, NULL);

            // with weights and an index, and without: the packed kernel and the general one
            const md_float3_t pbc    = com_planar(x, y, z, NULL, NULL, n, &cells[ci]);
            const md_float3_t pbc_w  = com_planar(x, y, z, w,    NULL, n, &cells[ci]);
            const md_float3_t pbc_i  = com_planar(x, y, z, NULL, idx,  n, &cells[ci]);
            const md_float3_t pbc_iw = com_planar(x, y, z, w,    idx,  n, &cells[ci]);

            // and strided, points of md_float4_t carrying their weight in w
            const md_float3_t pbc_v4 = md_geom_com((const float*)xyzw, sizeof(md_float4_t), (const float*)xyzw + 3, sizeof(md_float4_t), NULL, n, &cells[ci]);
            const md_float3_t pbc_v4i = md_geom_com((const float*)xyzw, sizeof(md_float4_t), (const float*)xyzw + 3, sizeof(md_float4_t), idx, n, &cells[ci]);

            for (int e = 0; e < 3; ++e) {
                EXPECT_NEAR(pbc.elem[e],     plain_u.elem[e], 0.1f);
                EXPECT_NEAR(pbc_i.elem[e],   plain_u.elem[e], 0.1f);
                EXPECT_NEAR(pbc_w.elem[e],   plain_w.elem[e], 0.1f);
                EXPECT_NEAR(pbc_iw.elem[e],  plain_w.elem[e], 0.1f);
                EXPECT_NEAR(pbc_v4.elem[e],  plain_w.elem[e], 0.1f);
                EXPECT_NEAR(pbc_v4i.elem[e], plain_w.elem[e], 0.1f);
            }

            // indexed and contiguous walk the same points, so they must agree exactly, not merely
            // to within the tolerance above
            for (int e = 0; e < 3; ++e) {
                EXPECT_NEAR(pbc_i.elem[e],  pbc.elem[e],   1.0e-4f);
                EXPECT_NEAR(pbc_iw.elem[e], pbc_w.elem[e], 1.0e-4f);
            }
        }
    }
}

// The defining property, and the one that needs no tolerance argument: moving individual points by
// whole lattice vectors does not change the configuration, so it must not move the centre. A
// per axis circular mean is periodic in the cell EXTENT rather than in the lattice, which is the
// same thing for a cube and is not the same thing at all for a triclinic cell.
UTEST(geom, com_pbc_invariant_under_lattice_shift) {
    md_unitcell_t cells[2];
    cells[0] = com_cell_ortho();
    cells[1] = com_cell_triclinic();

    for (int ci = 0; ci < 2; ++ci) {
        md_float3x3_t A = {0};
        md_unitcell_A_extract_float(A.elem, &cells[ci]);

        for (size_t k = 0; k < MD_COUNTOF(COM_COUNTS); ++k) {
            const size_t n = COM_COUNTS[k];
            float x[COM_MAX_PTS], y[COM_MAX_PTS], z[COM_MAX_PTS], w[COM_MAX_PTS];
            int32_t idx[COM_MAX_PTS];
            com_make_blob(x, y, z, w, n, &cells[ci]);
            for (size_t i = 0; i < n; ++i) idx[i] = (int32_t)i;

            const md_float3_t before    = com_planar(x, y, z, NULL, NULL, n, &cells[ci]);
            const md_float3_t before_w  = com_planar(x, y, z, w,    NULL, n, &cells[ci]);
            const md_float3_t before_i  = com_planar(x, y, z, NULL, idx,  n, &cells[ci]);
            const md_float3_t before_iw = com_planar(x, y, z, w,    idx,  n, &cells[ci]);

            // scatter the points across images - a different lattice vector for every third one
            for (size_t i = 0; i < n; ++i) {
                if (i % 3 != 0) continue;
                // NOTE: cast before subtracting - i is size_t and the subtraction would wrap
                const md_float3_t nvec = md_float3_set((float)((int)(i % 5) - 2), (float)((int)(i % 3) - 1), (float)((int)(i % 7) - 3));
                const md_float3_t shift = md_float3x3_mul_float3(A, nvec);
                x[i] += shift.x;
                y[i] += shift.y;
                z[i] += shift.z;
            }

            const md_float3_t after    = com_planar(x, y, z, NULL, NULL, n, &cells[ci]);
            const md_float3_t after_w  = com_planar(x, y, z, w,    NULL, n, &cells[ci]);
            const md_float3_t after_i  = com_planar(x, y, z, NULL, idx,  n, &cells[ci]);
            const md_float3_t after_iw = com_planar(x, y, z, w,    idx,  n, &cells[ci]);

            for (int e = 0; e < 3; ++e) {
                EXPECT_NEAR(after.elem[e],    before.elem[e],    0.02f);
                EXPECT_NEAR(after_w.elem[e],  before_w.elem[e],  0.02f);
                EXPECT_NEAR(after_i.elem[e],  before_i.elem[e],  0.02f);
                EXPECT_NEAR(after_iw.elem[e], before_iw.elem[e], 0.02f);
            }
        }
    }
}

// The indexed variants must read through the index list and nothing else. Bury the real points in a
// larger array whose unselected slots hold coordinates far away, and the answer must not move.
UTEST(geom, com_pbc_indices_select) {
    md_unitcell_t cells[2];
    cells[0] = com_cell_ortho();
    cells[1] = com_cell_triclinic();

    for (int ci = 0; ci < 2; ++ci) {
        const size_t n = 61;
        float xs[COM_MAX_PTS], ys[COM_MAX_PTS], zs[COM_MAX_PTS], ws[COM_MAX_PTS];
        com_make_blob(xs, ys, zs, ws, n, &cells[ci]);

        const md_float3_t want   = com_planar(xs, ys, zs, NULL, NULL, n, &cells[ci]);
        const md_float3_t want_w = com_planar(xs, ys, zs, ws,   NULL, n, &cells[ci]);

        // stride the real points through a 3x larger array, junk in between
        float bx[COM_MAX_PTS * 3], by[COM_MAX_PTS * 3], bz[COM_MAX_PTS * 3], bw[COM_MAX_PTS * 3];
        int32_t idx[COM_MAX_PTS];
        for (size_t i = 0; i < n * 3; ++i) {
            bx[i] = -931.0f + (float)i;
            by[i] =  757.0f - (float)i;
            bz[i] =  613.0f + (float)(i * 2);
            bw[i] =  99.0f;
        }
        for (size_t i = 0; i < n; ++i) {
            const size_t slot = i * 3 + 2;
            bx[slot] = xs[i]; by[slot] = ys[i]; bz[slot] = zs[i]; bw[slot] = ws[i];
            idx[i] = (int32_t)slot;
        }

        const md_float3_t got   = com_planar(bx, by, bz, NULL, idx, n, &cells[ci]);
        const md_float3_t got_w = com_planar(bx, by, bz, bw,   idx, n, &cells[ci]);

        for (int e = 0; e < 3; ++e) {
            EXPECT_NEAR(got.elem[e],   want.elem[e],   1.0e-4f);
            EXPECT_NEAR(got_w.elem[e], want_w.elem[e], 1.0e-4f);
        }
    }
}

// Weights have to actually weight. Two points, one ten times heavier, in a triclinic cell: the
// centre belongs near the heavy one, and nowhere near the midpoint.
UTEST(geom, com_pbc_weights_apply) {
    md_unitcell_t cell = com_cell_triclinic();
    md_float3x3_t A = {0};
    md_unitcell_A_extract_float(A.elem, &cell);

    const md_float3_t p0 = md_float3x3_mul_float3(A, md_float3_set(0.30f, 0.30f, 0.30f));
    const md_float3_t p1 = md_float3x3_mul_float3(A, md_float3_set(0.34f, 0.36f, 0.38f));

    float x[2] = { p0.x, p1.x };
    float y[2] = { p0.y, p1.y };
    float z[2] = { p0.z, p1.z };
    float w[2] = { 10.0f, 1.0f };

    const md_float3_t expect = md_float3_div1(md_float3_add(md_float3_mul1(p0, 10.0f), p1), 11.0f);
    const md_float3_t got    = com_planar(x, y, z, w, NULL, 2, &cell);

    EXPECT_NEAR(got.x, expect.x, 0.05f);
    EXPECT_NEAR(got.y, expect.y, 0.05f);
    EXPECT_NEAR(got.z, expect.z, 0.05f);
}

// The rotation must be recovered from a target whose points arrive scattered across images, for any
// rotation, without ever consulting topology.
UTEST(geom, pbc_optimal_rotation) {
    const float L = 20.0f;
    md_unitcell_t cell = md_unitcell_from_extent(L, L, L);

    for (int deg = 0; deg <= 180; deg += 20) {
        md_float4_t ref[PBC_NP];
        pbc_make_shape(ref, 3.0f);
        const md_float3_t ref_com = md_geom_com((const float*)ref, sizeof(md_float4_t), (const float*)ref + 3, sizeof(md_float4_t), 0, PBC_NP, 0);

        const md_float3x3_t R_true = pbc_rot(deg * 3.14159265f / 180.0f);
        const md_float3_t c_true = { 3.0f, 17.0f, 9.0f };   // deliberately near two boundaries

        md_float4_t trg[PBC_NP];
        for (int i = 0; i < PBC_NP; ++i) {
            const md_float3_t v = md_float3_add(md_float3x3_mul_float3(R_true, md_float3_sub(md_float3_from_float4(ref[i]), ref_com)), c_true);
            trg[i] = md_float4_set(pbc_wrap(v.x, L), pbc_wrap(v.y, L), pbc_wrap(v.z, L), ref[i].w);
        }

        md_float3x3_t R = {0};
        md_float3_t com = {0};
        md_float4_t placed[PBC_NP];
        const float residual = md_geom_optimal_rotation_pbc(&R, &com, placed, ref, ref_com, trg, PBC_NP, &cell, 8, 1.0e-6f);

        EXPECT_LT(residual, 1.0e-2f);
        EXPECT_LT(pbc_angle_err(R, R_true), 0.5f);

        // The placed points must satisfy the transform the call reported
        for (int i = 0; i < PBC_NP; ++i) {
            const md_float3_t p = md_float3_sub(md_float3_from_float4(ref[i]), ref_com);
            const md_float3_t q = md_float3x3_mul_float3(R, md_float3_sub(md_float3_from_float4(placed[i]), com));
            EXPECT_NEAR(q.x, p.x, 1.0e-2f);
            EXPECT_NEAR(q.y, p.y, 1.0e-2f);
            EXPECT_NEAR(q.z, p.z, 1.0e-2f);
            EXPECT_NEAR(placed[i].w, ref[i].w, 1.0e-6f);
        }
    }
}

// Beyond the regime where any per point image choice can work, the fit must SAY so rather than
// return a confident wrong answer. A set of radius 8.5 in a 20A cell flipped end for end is past it.
UTEST(geom, pbc_optimal_rotation_reports_ambiguity) {
    const float L = 20.0f;
    md_unitcell_t cell = md_unitcell_from_extent(L, L, L);

    md_float4_t ref[PBC_NP];
    pbc_make_shape(ref, 8.5f);
    const md_float3_t ref_com = md_geom_com((const float*)ref, sizeof(md_float4_t), (const float*)ref + 3, sizeof(md_float4_t), 0, PBC_NP, 0);

    const md_float3x3_t R_true = pbc_rot(3.14159265f);
    const md_float3_t c_true = { 3.0f, 17.0f, 9.0f };

    md_float4_t trg[PBC_NP];
    for (int i = 0; i < PBC_NP; ++i) {
        const md_float3_t v = md_float3_add(md_float3x3_mul_float3(R_true, md_float3_sub(md_float3_from_float4(ref[i]), ref_com)), c_true);
        trg[i] = md_float4_set(pbc_wrap(v.x, L), pbc_wrap(v.y, L), pbc_wrap(v.z, L), ref[i].w);
    }

    md_float3x3_t R = {0};
    md_float3_t com = {0};
    const float residual = md_geom_optimal_rotation_pbc(&R, &com, NULL, ref, ref_com, trg, PBC_NP, &cell, 8, 1.0e-6f);

    // The certificate is the whole point: a residual on the order of half a cell means do not trust it
    EXPECT_GT(residual, 0.25f * L);
}

// Degenerate inputs must not misbehave
UTEST(geom, pbc_optimal_rotation_degenerate) {
    md_unitcell_t cell = md_unitcell_from_extent(20, 20, 20);
    md_float4_t ref[PBC_NP];
    pbc_make_shape(ref, 3.0f);
    const md_float3_t ref_com = md_geom_com((const float*)ref, sizeof(md_float4_t), (const float*)ref + 3, sizeof(md_float4_t), 0, PBC_NP, 0);

    md_float3x3_t R = {0};
    md_float3_t com = {0};

    // Empty set
    EXPECT_NEAR(md_geom_optimal_rotation_pbc(&R, &com, NULL, ref, ref_com, ref, 0, &cell, 8, 1.0e-6f), 0.0f, 1.0e-6f);

    // No cell at all: still a plain Kabsch fit, identity here since target == reference
    const float residual = md_geom_optimal_rotation_pbc(&R, &com, NULL, ref, ref_com, ref, PBC_NP, NULL, 8, 1.0e-6f);
    EXPECT_LT(residual, 1.0e-3f);
    EXPECT_LT(pbc_angle_err(R, md_float3x3_ident()), 0.5f);

    // A set with no coordinates to speak of: its center is the origin
    md_float3_t c = md_float3_set(1, 2, 3);
    md_geom_deperiodize_self(ref, 0, &cell, &c);
    EXPECT_EQ(0.0f, c.x);
    EXPECT_EQ(0.0f, c.y);
    EXPECT_EQ(0.0f, c.z);
}

// An atom that crosses the cell boundary between two frames goes the short way, across the
// boundary, and not back through the box. The orthorhombic branch once discarded its
// minimum image and so did exactly that.
UTEST(geom, interpolate_linear_across_boundary) {
    const md_unitcell_t cell = md_unitcell_from_extent(10, 10, 10);
    md_float3_t a[16] = {0}, b[16] = {0}, out[16] = {0};
    for (int i = 0; i < 16; ++i) {
        a[i] = md_float3_set(9.5f, 5.0f, 0.25f);
        b[i] = md_float3_set(0.5f, 5.0f, 9.75f);
    }
    const md_float3_t* const in[2] = { a, b };
    md_geom_interpolate_linear(out, in, 16, &cell, 0.5f);
    for (int i = 0; i < 16; ++i) {
        // Halfway is on the boundary: 10 or its image 0, never 5
        EXPECT_NEAR(0.0f, fabsf(fmodf(out[i].x + 5.0f, 10.0f) - 5.0f), 1.0e-4f);
        EXPECT_NEAR(5.0f, out[i].y, 1.0e-4f);
        EXPECT_NEAR(0.0f, fabsf(fmodf(out[i].z + 5.0f, 10.0f) - 5.0f), 1.0e-4f);
    }

    // Without a cell it is a plain blend
    md_geom_interpolate_linear(out, in, 16, NULL, 0.5f);
    EXPECT_NEAR(5.0f, out[3].x, 1.0e-4f);
}

// Wrapping a subset moves those atoms and no others; the orthorhombic branch once wrote each result
// to the i:th atom instead of the one indexed.
UTEST(geom, pbc_indexed) {
    const md_unitcell_t cell = md_unitcell_from_extent(10, 10, 10);
    md_float3_t xyz[20];
    for (int i = 0; i < 20; ++i) xyz[i] = md_float3_set(12.0f + i, -3.0f, 5.0f);

    // Without a cell, or one without periodic axes, there is nothing to wrap into
    const md_unitcell_t none = {0};
    md_geom_pbc((float*)xyz, 0, NULL, 20, NULL);
    md_geom_pbc((float*)xyz, 0, NULL, 20, &none);
    for (int i = 0; i < 20; ++i) EXPECT_EQ(12.0f + i, xyz[i].x);

    const int32_t idx[3] = { 17, 4, 11 };
    md_geom_pbc((float*)xyz, 0, idx, 3, &cell);
    for (int i = 0; i < 20; ++i) {
        const bool wrapped = (i == 17 || i == 4 || i == 11);
        const float ex = wrapped ? fmodf(12.0f + i, 10.0f) : 12.0f + i;
        EXPECT_NEAR(ex, xyz[i].x, 1.0e-4f);
        EXPECT_NEAR(wrapped ? 7.0f : -3.0f, xyz[i].y, 1.0e-4f);
    }

    // All of them, through the eight wide path and its tail. An atom on the boundary may come out
    // at either end, which is the same place.
    md_geom_pbc((float*)xyz, 0, NULL, 20, &cell);
    for (int i = 0; i < 20; ++i) {
        EXPECT_GE(xyz[i].x, 0.0f);
        EXPECT_LE(xyz[i].x, 10.0f);
        EXPECT_NEAR(0.0f, remainderf(xyz[i].x - (12.0f + i), 10.0f), 1.0e-4f);
        EXPECT_NEAR(7.0f, xyz[i].y, 1.0e-4f);
        EXPECT_NEAR(5.0f, xyz[i].z, 1.0e-4f);
    }

    // Triclinic: every wrapped atom lands inside the cell
    const md_unitcell_t tri = md_unitcell_from_basis_parameters(10, 10, 10, 2, 1, 3);
    for (int i = 0; i < 20; ++i) xyz[i] = md_float3_set(12.0f + i, -3.0f - i, 25.0f);
    md_geom_pbc((float*)xyz, 0, NULL, 20, &tri);
    md_float3x3_t I = {0};
    md_unitcell_I_extract_float(I.elem, &tri);
    for (int i = 0; i < 20; ++i) {
        const md_float3_t f = md_float3x3_mul_float3(I, xyz[i]);
        for (int k = 0; k < 3; ++k) {
            EXPECT_GE(f.elem[k], -1.0e-4f);
            EXPECT_LE(f.elem[k], 1.0f + 1.0e-4f);
        }
    }
}

// Points of md_float4_t, by their stride: the same positions as packed, in both kinds of cell, all of them or some,
// and their w left alone.
UTEST(geom, pbc_strided) {
    const md_unitcell_t cells[2] = {
        md_unitcell_from_extent(10, 10, 10),
        md_unitcell_from_basis_parameters(10, 10, 10, 2, 1, 3),
    };
    const int32_t idx[3] = { 17, 4, 11 };
    for (int ci = 0; ci < 2; ++ci) {
        for (int with_idx = 0; with_idx < 2; ++with_idx) {
            md_float3_t xyz[20];
            md_float4_t xyzw[20];
            for (int i = 0; i < 20; ++i) {
                xyz[i]  = md_float3_set(12.0f + i, -3.0f - i, 25.0f - 0.5f * i);
                xyzw[i] = md_float4_from_float3(xyz[i], (float)i);
            }
            const int32_t* ii = with_idx ? idx : NULL;
            const size_t   n  = with_idx ? 3 : 20;
            md_geom_pbc((float*)xyz, 0, ii, n, &cells[ci]);
            md_geom_pbc((float*)xyzw, sizeof(md_float4_t), ii, n, &cells[ci]);
            for (int i = 0; i < 20; ++i) {
                for (int k = 0; k < 3; ++k) {
                    EXPECT_NEAR(xyz[i].elem[k], xyzw[i].elem[k], 1.0e-4f);
                }
                EXPECT_EQ((float)i, xyzw[i].w);
            }
        }
    }
}

// Each of the four paths (packed, packed with an index, strided, strided with an index) against a scalar reference,
// on a count that leaves a tail after the eight wide part.
UTEST(geom, aabb_paths) {
    enum { N = 37 };
    md_float3_t xyz[N];
    md_float4_t xyzw[N];
    int32_t idx[N];
    for (int i = 0; i < N; ++i) {
        const float s = (float)((i * 2654435761u) % 1000) / 100.0f - 5.0f;
        xyz[i]  = md_float3_set(s, -2.0f * s + i, 0.5f * i - s);
        xyzw[i] = md_float4_from_float3(xyz[i], 1.0e9f);
        idx[i]  = (i * 11) % N;
    }
    for (int variant = 0; variant < 4; ++variant) {
        const bool     strided = (variant & 1) != 0;
        const int32_t* ii = (variant & 2) ? idx : NULL;
        const size_t   n  = ii ? 29 : N;
        float ref_min[3] = { FLT_MAX,  FLT_MAX,  FLT_MAX};
        float ref_max[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
        for (size_t k = 0; k < n; ++k) {
            const int32_t j = ii ? ii[k] : (int32_t)k;
            for (int c = 0; c < 3; ++c) {
                ref_min[c] = MD_MIN(ref_min[c], xyz[j].elem[c]);
                ref_max[c] = MD_MAX(ref_max[c], xyz[j].elem[c]);
            }
        }
        float got_min[3], got_max[3];
        if (strided) {
            md_geom_aabb(got_min, got_max, (const float*)xyzw, sizeof(md_float4_t), ii, n);
        } else {
            md_geom_aabb(got_min, got_max, (const float*)xyz, 0, ii, n);
        }
        for (int c = 0; c < 3; ++c) {
            EXPECT_EQ(ref_min[c], got_min[c]);
            EXPECT_EQ(ref_max[c], got_max[c]);
        }
    }
}

// The spline passes through the two middle frames at t = 0 and t = 1, and across a boundary it
// takes the neighbouring frames' images rather than their wrapped positions.
UTEST(geom, interpolate_cubic) {
    md_float3_t f[4][16], out[16];
    for (int i = 0; i < 16; ++i) {
        for (int k = 0; k < 4; ++k) {
            f[k][i] = md_float3_set(1.0f * k + i, 2.0f * k, -0.5f * k * k);
        }
    }
    const md_float3_t* const in[4] = { f[0], f[1], f[2], f[3] };
    float w[4];
    md_geom_cubic_spline_weights(w, 0.0f, 0.5f);
    md_geom_interpolate_cubic_weighted(out, in, 16, NULL, w);
    for (int i = 0; i < 16; ++i) {
        EXPECT_NEAR(f[1][i].x, out[i].x, 1.0e-4f);
        EXPECT_NEAR(f[1][i].z, out[i].z, 1.0e-4f);
    }
    md_geom_cubic_spline_weights(w, 1.0f, 0.5f);
    md_geom_interpolate_cubic_weighted(out, in, 16, NULL, w);
    for (int i = 0; i < 16; ++i) {
        EXPECT_NEAR(f[2][i].y, out[i].y, 1.0e-4f);
    }

    // Moving +1 per frame along x, wrapped into a box of 10
    const md_unitcell_t cell = md_unitcell_from_extent(10, 10, 10);
    for (int i = 0; i < 16; ++i) {
        for (int k = 0; k < 4; ++k) {
            f[k][i] = md_float3_set(fmodf(8.0f + k, 10.0f), 5.0f, 5.0f);   // 8, 9, 0, 1
        }
    }
    md_geom_cubic_spline_weights(w, 0.5f, 0.5f);
    md_geom_interpolate_cubic_weighted(out, in, 16, &cell, w);
    for (int i = 0; i < 16; ++i) {
        EXPECT_NEAR(9.5f, out[i].x, 1.0e-3f);
    }
}

// Both schemes' weights sum to one and reproduce uniform motion. The cardinal spline lands on the middle frames exactly,
// through the weighted sum too; the B-spline is short of them by design, at (p0 + 4 p1 + p2) / 6. Across a periodic
// boundary, orthorhombic or triclinic, the weighted sum takes the frames' images.
UTEST(geom, interpolate_cubic_weights) {
    for (int i = 0; i <= 16; ++i) {
        const float t = i / 16.0f;
        float wc[4], wb[4];
        md_geom_cubic_spline_weights(wc, t, 0.5f);
        md_geom_cubic_bspline_weights(wb, t);
        EXPECT_NEAR(1.0f, wc[0] + wc[1] + wc[2] + wc[3], 1.0e-6f);
        EXPECT_NEAR(1.0f, wb[0] + wb[1] + wb[2] + wb[3], 1.0e-6f);
        // Frames at 0, 1, 2 and 3: the interval between the middle two is at 1 + t
        EXPECT_NEAR(1.0f + t, wc[1] + 2.0f * wc[2] + 3.0f * wc[3], 1.0e-6f);
        EXPECT_NEAR(1.0f + t, wb[1] + 2.0f * wb[2] + 3.0f * wb[3], 1.0e-6f);
    }

    float w[4];
    md_geom_cubic_spline_weights(w, 0.0f, 0.5f);
    EXPECT_EQ(0.0f, fabsf(w[0]));
    EXPECT_EQ(1.0f, w[1]);
    EXPECT_EQ(0.0f, fabsf(w[2]));
    EXPECT_EQ(0.0f, fabsf(w[3]));
    md_geom_cubic_spline_weights(w, 1.0f, 0.5f);
    EXPECT_EQ(0.0f, fabsf(w[0]));
    EXPECT_EQ(0.0f, fabsf(w[1]));
    EXPECT_EQ(1.0f, w[2]);
    EXPECT_EQ(0.0f, fabsf(w[3]));
    md_geom_cubic_bspline_weights(w, 0.0f);
    EXPECT_NEAR(1.0f / 6.0f, w[0], 1.0e-7f);
    EXPECT_NEAR(4.0f / 6.0f, w[1], 1.0e-7f);
    EXPECT_NEAR(1.0f / 6.0f, w[2], 1.0e-7f);
    EXPECT_NEAR(0.0f,        w[3], 1.0e-7f);

    md_float3_t f[4][16], out[16];
    for (int i = 0; i < 16; ++i) {
        for (int k = 0; k < 4; ++k) {
            f[k][i] = md_float3_set(1.0f * k + i, 2.0f * k * k, -0.5f * k + 0.1f * i);
        }
    }
    const md_float3_t* const in[4] = { f[0], f[1], f[2], f[3] };

    md_geom_cubic_spline_weights(w, 0.0f, 0.5f);
    md_geom_interpolate_cubic_weighted(out, in, 16, NULL, w);
    for (int i = 0; i < 16; ++i) {
        EXPECT_EQ(f[1][i].x, out[i].x);
        EXPECT_EQ(f[1][i].y, out[i].y);
        EXPECT_EQ(f[1][i].z, out[i].z);
    }

    md_geom_cubic_bspline_weights(w, 0.0f);
    md_geom_interpolate_cubic_weighted(out, in, 16, NULL, w);
    for (int i = 0; i < 16; ++i) {
        for (int c = 0; c < 3; ++c) {
            const float ref = (f[0][i].elem[c] + 4.0f * f[1][i].elem[c] + f[2][i].elem[c]) / 6.0f;
            EXPECT_NEAR(ref, out[i].elem[c], 1.0e-4f);
        }
    }

    // Moving +1 per frame along x, wrapped into a box of 10: 8, 9, 0, 1
    for (int i = 0; i < 16; ++i) {
        for (int k = 0; k < 4; ++k) {
            f[k][i] = md_float3_set(fmodf(8.0f + k, 10.0f), 5.0f, 5.0f);
        }
    }
    const md_unitcell_t cells[2] = {
        md_unitcell_from_extent(10, 10, 10),
        md_unitcell_from_basis_parameters(10, 10, 10, 2, 0, 0),
    };
    for (int c = 0; c < 2; ++c) {
        md_geom_cubic_bspline_weights(w, 0.5f);
        md_geom_interpolate_cubic_weighted(out, in, 16, &cells[c], w);
        for (int i = 0; i < 16; ++i) {
            EXPECT_NEAR(9.5f, out[i].x, 1.0e-3f);
            EXPECT_NEAR(5.0f, out[i].y, 1.0e-3f);
        }
    }
}

// The component wise kernels step over a few atoms to align their loads, so where an array begins
// decides how they split the work. Every start, against scalar references.
UTEST(geom, packed_kernels_any_alignment) {
    enum { N = 37, PAD = 64 };
    const md_unitcell_t cell = md_unitcell_from_extent(10, 12, 14);
    const md_float3_t ext = {10, 12, 14};
    md_float3_t base_a[N + PAD], base_b[N + PAD], base_o[N + PAD];
    float  base_w[N + PAD];
    for (int off = 0; off < 8; ++off) {
        md_float3_t* a = base_a + off;
        md_float3_t* b = base_b + off;
        md_float3_t* o = base_o + off;
        float*  w = base_w + off;
        for (int i = 0; i < N; ++i) {
            const float s = (float)((i * 2654435761u + off) % 1000) / 1000.0f;
            a[i] = md_float3_set(9.0f * s + 0.5f, 11.0f * (1.0f - s) + 0.5f, 13.0f * s * s + 0.5f);
            b[i] = md_float3_add(a[i], md_float3_set(3.0f - 6.0f * s, 5.0f * s - 2.5f, 1.5f));
            w[i] = 0.5f + s;
        }

        // Plain and weighted COM without a cell
        double ref[3] = {0}, ref_w[3] = {0}, sw = 0;
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < 3; ++k) {
                ref[k]   += a[i].elem[k];
                ref_w[k] += a[i].elem[k] * w[i];
            }
            sw += w[i];
        }
        const md_float3_t c  = md_geom_com((const float*)a, 0, NULL, 0, NULL, N, NULL);
        const md_float3_t cw = md_geom_com((const float*)a, 0, w, 0, NULL, N, NULL);
        for (int k = 0; k < 3; ++k) {
            EXPECT_NEAR(ref[k] / N, c.elem[k], 1.0e-4);
            EXPECT_NEAR(ref_w[k] / sw, cw.elem[k], 1.0e-4);
        }

        // Periodic COM of a cloud well inside the cell is its plain mean, weighted or not
        md_float3_t near_center[N];
        for (int i = 0; i < N; ++i) near_center[i] = md_float3_add(md_float3_set(5, 6, 7), md_float3_mul1(md_float3_sub(a[i], md_float3_set(5, 6, 7)), 0.1f));
        md_memcpy(o, near_center, sizeof(near_center));
        const md_float3_t pc  = md_geom_com((const float*)o, 0, NULL, 0, NULL, N, &cell);
        const md_float3_t pcw = md_geom_com((const float*)o, 0, w, 0, NULL, N, &cell);
        const md_float3_t mc  = md_geom_com((const float*)o, 0, NULL, 0, NULL, N, NULL);
        const md_float3_t mcw = md_geom_com((const float*)o, 0, w, 0, NULL, N, NULL);
        for (int k = 0; k < 3; ++k) {
            EXPECT_NEAR(mc.elem[k],  pc.elem[k],  2.0e-2);
            EXPECT_NEAR(mcw.elem[k], pcw.elem[k], 2.0e-2);
        }

        // AABB
        float mn[3], mx[3];
        md_geom_aabb(mn, mx, (const float*)a, 0, NULL, N);
        for (int k = 0; k < 3; ++k) {
            float lo = FLT_MAX, hi = -FLT_MAX;
            for (int i = 0; i < N; ++i) { lo = MD_MIN(lo, a[i].elem[k]); hi = MD_MAX(hi, a[i].elem[k]); }
            EXPECT_EQ(lo, mn[k]);
            EXPECT_EQ(hi, mx[k]);
        }

        // Linear interpolation, no cell and orthorhombic, and nothing written past the end
        base_o[off + N] = md_float3_set(-1, -1, -1);
        const md_float3_t* const in[2] = { a, b };
        md_geom_interpolate_linear(o, in, N, NULL, 0.25f);
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < 3; ++k) {
                EXPECT_NEAR(a[i].elem[k] + 0.25f * (b[i].elem[k] - a[i].elem[k]), o[i].elem[k], 1.0e-4f);
            }
        }
        md_geom_interpolate_linear(o, in, N, &cell, 0.25f);
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < 3; ++k) {
                float d = b[i].elem[k] - a[i].elem[k];
                d -= ext.elem[k] * roundf(d / ext.elem[k]);
                EXPECT_NEAR(a[i].elem[k] + 0.25f * d, o[i].elem[k], 1.0e-4f);
            }
        }
        EXPECT_EQ(-1.0f, base_o[off + N].x);

        // Cubic through the middle frames at the ends of the interval
        const md_float3_t* const in4[4] = { a, a, b, b };
        float w4[4];
        md_geom_cubic_spline_weights(w4, 1.0f, 0.5f);
        md_geom_interpolate_cubic_weighted(o, in4, N, NULL, w4);
        for (int i = 0; i < N; ++i) {
            EXPECT_NEAR(b[i].y, o[i].y, 1.0e-4f);
        }
        EXPECT_EQ(-1.0f, base_o[off + N].x);

        // Wrapping
        md_memcpy(o, b, N * sizeof(md_float3_t));
        md_geom_pbc((float*)o, 0, NULL, N, &cell);
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < 3; ++k) {
                EXPECT_GE(o[i].elem[k], -1.0e-4f);
                EXPECT_LE(o[i].elem[k], ext.elem[k] + 1.0e-4f);
                EXPECT_NEAR(0.0f, remainderf(o[i].elem[k] - b[i].elem[k], ext.elem[k]), 1.0e-4f);
            }
        }
    }
}

// ---- minimum and maximum distance between sets ----------------------------------------------

// The minimum image by exhaustion, in double: reduced into the brick, then every combination of up to three lattice
// vectors either way along the periodic axes. Three is more than any of the cells below needs: a vector A n which
// shortens a reduced d has |n_i| <= 2 |d| |r_i| (r_i a row of the inverse basis), at most 2 for these.
static double dist_ref(md_float3_t a, md_float3_t b, const md_unitcell_t* cell) {
    double d[3] = { (double)a.x - b.x, (double)a.y - b.y, (double)a.z - b.z };
    double A[3][3] = {{0}};
    md_unitcell_A_extract_double(A, cell);
    int per[3] = { (cell->flags & MD_UNITCELL_PBC_X) != 0, (cell->flags & MD_UNITCELL_PBC_Y) != 0, (cell->flags & MD_UNITCELL_PBC_Z) != 0 };
    for (int i = 0; i < 3; ++i) {
        if (!per[i] || A[i][i] == 0.0) { A[i][0] = A[i][1] = A[i][2] = 0.0; per[i] = 0; }
    }
    for (int i = 2; i >= 0; --i) {
        if (per[i]) {
            const double n = nearbyint(d[i] / A[i][i]);
            for (int j = 0; j < 3; ++j) d[j] -= n * A[i][j];
        }
    }
    const int R = 3;
    double best = DBL_MAX;
    for (int nz = per[2] ? -R : 0; nz <= (per[2] ? R : 0); ++nz) {
        for (int ny = per[1] ? -R : 0; ny <= (per[1] ? R : 0); ++ny) {
            for (int nx = per[0] ? -R : 0; nx <= (per[0] ? R : 0); ++nx) {
                double v[3];
                for (int j = 0; j < 3; ++j) v[j] = d[j] + nx * A[0][j] + ny * A[1][j] + nz * A[2][j];
                best = MD_MIN(best, v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            }
        }
    }
    return sqrt(best);
}

static uint64_t dist_rng_state;

static double dist_rand(void) {
    uint64_t x = dist_rng_state;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    dist_rng_state = x;
    return (double)(x >> 11) * (1.0 / 9007199254740992.0);
}

// A point at fractional coordinates in [lo, hi) along the periodic axes - several images of the cell, as in an
// unwrapped trajectory - and within span along the others
static md_float3_t dist_rand_point(const md_unitcell_t* cell, double span, double lo, double hi) {
    double A[3][3] = {{0}};
    md_unitcell_A_extract_double(A, cell);
    double p[3] = {0};
    for (int i = 0; i < 3; ++i) {
        const double f = lo + (hi - lo) * dist_rand();
        if (A[i][i] != 0.0) {
            for (int j = 0; j < 3; ++j) p[j] += f * A[i][j];
        } else {
            p[i] += f * span;
        }
    }
    return md_float3_set((float)p[0], (float)p[1], (float)p[2]);
}

// Groups which are empty, compact, spread over several images, or share points with b, in every kind of cell -
// including cells small enough that most distances are beyond half of them, where reducing a displacement into the
// cell is not the minimum image, and a skewed cell which is not reduced - against the minimum image by exhaustion
UTEST(geom, min_max_distance_groups) {
    const double s2 = sqrt(2.0), s6 = sqrt(6.0), d = 30.0, ds = 11.0;
    const md_unitcell_t cells[] = {
        md_unitcell_none(),
        md_unitcell_from_extent(30, 30, 30),
        md_unitcell_from_extent(9, 12, 10),
        md_unitcell_from_basis_parameters(30, 40, 0, 0, 0, 0),                              // periodic in x and y only
        md_unitcell_from_extent_and_angles(30, 30, 30, 70, 70, 70),
        md_unitcell_from_basis_parameters(d, d, d / s2, 0, d / 2, d / 2),                   // rhombic dodecahedron
        md_unitcell_from_basis_parameters(d, 2 * s2 * d / 3, s6 * d / 3, d / 3, -d / 3, s2 * d / 3),  // truncated octahedron
        md_unitcell_from_basis_parameters(ds, ds, ds / s2, 0, ds / 2, ds / 2),              // small dodecahedron
        md_unitcell_from_basis_parameters(30, 20, 25, 27, 20, 15),                          // skewed, not reduced
        md_unitcell_from_basis_parameters(30, 30, 0, 10, 0, 0),                             // triclinic, periodic in x and y only
    };

    dist_rng_state = 0x243F6A8885A308D3ULL;
    md_float3_t  a[6 * 12];
    md_float3_t  b[60];
    size_t  off[7];
    float   dist[6];
    int64_t ia[6], ib[6];

    for (size_t ci = 0; ci < MD_COUNTOF(cells); ++ci) {
        const md_unitcell_t* cell = &cells[ci];
        for (int trial = 0; trial < 25; ++trial) {
            const size_t nb = trial == 0 ? 0 : (trial % 5 == 0 ? 1 + (size_t)(dist_rand() * 3) : 1 + (size_t)(dist_rand() * 59));
            const bool blob = trial % 3 == 0;
            const md_float3_t bc = dist_rand_point(cell, 40, 0, 1);
            for (size_t j = 0; j < nb; ++j) {
                b[j] = blob ? md_float3_add(bc, md_float3_set((float)(dist_rand() * 8 - 4), (float)(dist_rand() * 8 - 4), (float)(dist_rand() * 8 - 4)))
                            : dist_rand_point(cell, 40, -1.5, 2.5);
            }

            const size_t ng = 1 + (size_t)(dist_rand() * 6);
            size_t na = 0;
            for (size_t g = 0; g < ng; ++g) {
                off[g] = na;
                const size_t n    = dist_rand() < 0.2 ? 0 : 1 + (size_t)(dist_rand() * 12);
                const int    kind = (int)(dist_rand() * 4);
                const md_float3_t cen  = dist_rand_point(cell, 40, -1, 2);
                for (size_t k = 0; k < n; ++k) {
                    if (kind == 0) {
                        a[na++] = dist_rand_point(cell, 40, -1.5, 2.5);
                    } else if (kind == 1 && nb > 0) {
                        a[na++] = b[(size_t)(dist_rand() * nb)];
                    } else {
                        a[na++] = md_float3_add(cen, md_float3_set((float)(dist_rand() * 5 - 2.5), (float)(dist_rand() * 5 - 2.5), (float)(dist_rand() * 5 - 2.5)));
                    }
                }
            }
            off[ng] = na;

            for (int largest = 0; largest < 2; ++largest) {
                if (largest) {
                    md_geom_max_distance_groups(dist, ia, ib, a, off, ng, b, nb, cell);
                } else {
                    md_geom_min_distance_groups(dist, ia, ib, a, off, ng, b, nb, cell);
                }
                for (size_t g = 0; g < ng; ++g) {
                    if (off[g] == off[g + 1] || nb == 0) {
                        EXPECT_EQ(0.0f, dist[g]);
                        EXPECT_EQ(-1, ia[g]);
                        EXPECT_EQ(-1, ib[g]);
                        continue;
                    }
                    double ref = largest ? -1.0 : DBL_MAX;
                    for (size_t i = off[g]; i < off[g + 1]; ++i) {
                        for (size_t j = 0; j < nb; ++j) {
                            const double r = dist_ref(a[i], b[j], cell);
                            ref = largest ? MD_MAX(ref, r) : MD_MIN(ref, r);
                        }
                    }
                    EXPECT_NEAR(ref, dist[g], 2.0e-3);
                    // The pair is one which has that distance
                    ASSERT_GE(ia[g], (int64_t)off[g]);
                    ASSERT_LT(ia[g], (int64_t)off[g + 1]);
                    ASSERT_GE(ib[g], 0);
                    ASSERT_LT(ib[g], (int64_t)nb);
                    EXPECT_NEAR(ref, dist_ref(a[ia[g]], b[ib[g]], cell), 2.0e-3);
                }
            }
        }
    }

    // A single group: 0 and no pair (-1) for an empty set
    const md_float3_t p = md_float3_set(1, 2, 3);
    const md_unitcell_t cell = md_unitcell_from_extent(10, 10, 10);
    const size_t none[2] = { 0, 0 }, one[2] = { 0, 1 };
    float dd = -1.0f;
    int64_t i0 = 7, i1 = 7;
    md_geom_min_distance_groups(&dd, &i0, &i1, &p, none, 1, &p, 1, &cell);
    EXPECT_EQ(0.0f, dd);
    EXPECT_EQ(-1, i0);
    EXPECT_EQ(-1, i1);
    md_geom_max_distance_groups(&dd, &i0, &i1, &p, one, 1, &p, 0, &cell);
    EXPECT_EQ(0.0f, dd);
    EXPECT_EQ(-1, i0);
    EXPECT_EQ(-1, i1);
    const md_float3_t q = md_float3_set(9.5f, 2, 3);
    md_geom_min_distance_groups(&dd, &i0, &i1, &p, one, 1, &q, 1, &cell);
    EXPECT_NEAR(1.5f, dd, 1.0e-5f);
    md_geom_min_distance_groups(&dd, &i0, &i1, &p, one, 1, &q, 1, NULL);
    EXPECT_NEAR(8.5f, dd, 1.0e-5f);
    EXPECT_EQ(0, i0);
    EXPECT_EQ(0, i1);
}

// The pairs near the extreme: every pair within tol of it and no other, against the minimum image by exhaustion, in
// every kind of cell. Points on a coarse grid make many pairs tie. Too small a cap says so.
UTEST(geom, min_max_distance_near) {
    const double s2 = sqrt(2.0), d = 30.0;
    const md_unitcell_t cells[] = {
        md_unitcell_none(),
        md_unitcell_from_extent(30, 30, 30),
        md_unitcell_from_extent(9, 12, 10),
        md_unitcell_from_basis_parameters(30, 40, 0, 0, 0, 0),
        md_unitcell_from_extent_and_angles(30, 30, 30, 70, 70, 70),
        md_unitcell_from_basis_parameters(d, d, d / s2, 0, d / 2, d / 2),
        md_unitcell_from_basis_parameters(30, 20, 25, 27, 20, 15),
    };
    dist_rng_state = 0x13198A2E03707344ULL;
    enum { NA = 70, NB = 90, CAP = NA * NB };
    md_float3_t a[NA], b[NB];
    int64_t* na_idx = malloc(CAP * sizeof(int64_t));
    int64_t* nb_idx = malloc(CAP * sizeof(int64_t));
    uint8_t* in_set = malloc(CAP);
    for (size_t ci = 0; ci < MD_COUNTOF(cells); ++ci) {
        const md_unitcell_t* cell = &cells[ci];
        for (int trial = 0; trial < 12; ++trial) {
            const size_t na = 1 + (size_t)(dist_rand() * (NA - 1)), nb = 1 + (size_t)(dist_rand() * (NB - 1));
            const float grid = trial % 3 == 0 ? 0.0f : trial % 3 == 1 ? 0.5f : 2.0f;
            for (size_t i = 0; i < na; ++i) a[i] = dist_rand_point(cell, 40, -1.5, 2.5);
            for (size_t j = 0; j < nb; ++j) b[j] = trial % 4 == 3 && j < na ? a[j] : dist_rand_point(cell, 40, -1.5, 2.5);
            if (grid > 0) {
                for (size_t i = 0; i < na; ++i) for (int k = 0; k < 3; ++k) a[i].elem[k] = roundf(a[i].elem[k] / grid) * grid;
                for (size_t j = 0; j < nb; ++j) for (int k = 0; k < 3; ++k) b[j].elem[k] = roundf(b[j].elem[k] / grid) * grid;
            }
            const float tol = trial % 2 ? 0.0f : 0.3f;
            for (int largest = 0; largest < 2; ++largest) {
                double ref = largest ? -1.0 : DBL_MAX;
                for (size_t i = 0; i < na; ++i) {
                    for (size_t j = 0; j < nb; ++j) {
                        const double r = dist_ref(a[i], b[j], cell);
                        ref = largest ? MD_MAX(ref, r) : MD_MIN(ref, r);
                    }
                }
                const size_t cnt = largest ? md_geom_max_distance_near(na_idx, nb_idx, CAP, tol, a, na, b, nb, cell)
                                           : md_geom_min_distance_near(na_idx, nb_idx, CAP, tol, a, na, b, nb, cell);
                ASSERT_LE(cnt, (size_t)CAP);
                ASSERT_GE(cnt, (size_t)1);
                memset(in_set, 0, CAP);
                for (size_t c = 0; c < cnt; ++c) {
                    ASSERT_GE(na_idx[c], 0);
                    ASSERT_LT(na_idx[c], (int64_t)na);
                    ASSERT_GE(nb_idx[c], 0);
                    ASSERT_LT(nb_idx[c], (int64_t)nb);
                    const size_t e = (size_t)na_idx[c] * nb + (size_t)nb_idx[c];
                    EXPECT_EQ(0, in_set[e]);      // each pair once
                    in_set[e] = 1;
                }
                // Every pair clearly within tol is there, every pair there is within it, to the kernel's rounding
                for (size_t i = 0; i < na; ++i) {
                    for (size_t j = 0; j < nb; ++j) {
                        const double off = largest ? ref - dist_ref(a[i], b[j], cell) : dist_ref(a[i], b[j], cell) - ref;
                        if (off <= tol - 1.0e-3) EXPECT_EQ(1, in_set[i * nb + j]);
                        if (in_set[i * nb + j]) EXPECT_LE(off, tol + 1.0e-3);
                    }
                }
                // Too small a cap: more than it
                if (cnt > 1) {
                    const size_t c2 = largest ? md_geom_max_distance_near(na_idx, nb_idx, cnt - 1, tol, a, na, b, nb, cell)
                                              : md_geom_min_distance_near(na_idx, nb_idx, cnt - 1, tol, a, na, b, nb, cell);
                    EXPECT_GT(c2, cnt - 1);
                }
            }
        }
    }
    int64_t i0 = 0, i1 = 0;
    EXPECT_EQ((size_t)0, md_geom_min_distance_near(&i0, &i1, 1, 0.0f, a, 0, b, 1, NULL));
    free(na_idx);
    free(nb_idx);
    free(in_set);
}

// The kernels that take packed positions give the same bits wherever the positions start in memory. Each runs on the
// same values copied to every offset of 4 bytes within 64, the alignments a gather or an allocator can hand it. A
// sum whose grouping followed the address (scalar up to an alignment, vector lanes after) fails here.
static double addr_rnd(double lo, double hi) { return lo + (hi - lo) * (double)rand() / (double)RAND_MAX; }

UTEST(geom, packed_kernels_independent_of_address) {
    enum { N = 103, OFFSETS = 16 };
    srand(4711);
    md_float3_t src[N];
    float w[N];
    for (int i = 0; i < N; ++i) {
        src[i] = (md_float3_t){{(float)addr_rnd(-5.0, 45.0), (float)addr_rnd(-5.0, 45.0), (float)addr_rnd(-5.0, 45.0)}};
        w[i] = (float)addr_rnd(1.0, 16.0);
    }
    const md_unitcell_t cell = md_unitcell_from_extent(40.0, 41.0, 42.0);
    const md_unitcell_t none = {0};
    // 64-byte aligned storage with room to slide the positions along it
    float* raw = (float*)md_alloc(md_get_heap_allocator(), (6 * N + 64) * sizeof(float) + 64);
    float* base = (float*)(((uintptr_t)raw + 63) & ~(uintptr_t)63);

    md_float3_t ref_com[4] = {0};
    float ref_min[3], ref_max[3];
    md_float3_t ref_pbc[N], ref_lin[N];
    for (int o = 0; o < OFFSETS; ++o) {
        md_float3_t* p = (md_float3_t*)(base + o);
        memcpy(p, src, sizeof(src));
        const md_float3_t com[4] = {
            md_geom_com((const float*)p, 0, NULL, 0, NULL, N, &none),
            md_geom_com((const float*)p, 0, w, 0, NULL, N, &none),
            md_geom_com((const float*)p, 0, NULL, 0, NULL, N, &cell),
            md_geom_com((const float*)p, 0, w, 0, NULL, N, &cell),
        };
        float mn[3], mx[3];
        md_geom_aabb(mn, mx, (const float*)p, 0, NULL, N);
        // pbc and interpolation write in place or to out: the result at this offset against the first offset's
        memcpy(p, src, sizeof(src));
        md_geom_pbc((float*)p, 0, NULL, N, &cell);
        md_float3_t pbc[N];
        memcpy(pbc, p, sizeof(pbc));
        memcpy(p, src, sizeof(src));
        md_float3_t* q = (md_float3_t*)(base + 3 * N + 32 - o);   // a second frame, at another alignment
        for (int i = 0; i < N; ++i) q[i] = md_float3_add(src[i], (md_float3_t){{0.5f, -0.25f, 0.125f}});
        const md_float3_t* frames[2] = {p, q};
        md_float3_t lin[N];
        md_geom_interpolate_linear(lin, frames, N, &cell, 0.37f);
        if (o == 0) {
            memcpy(ref_com, com, sizeof(com));
            memcpy(ref_min, mn, sizeof(mn));
            memcpy(ref_max, mx, sizeof(mx));
            memcpy(ref_pbc, pbc, sizeof(pbc));
            memcpy(ref_lin, lin, sizeof(lin));
            continue;
        }
        for (int k = 0; k < 4; ++k) {
            if (memcmp(&ref_com[k], &com[k], sizeof(md_float3_t))) printf("com %d differs at offset %d bytes: %.9g %.9g %.9g against %.9g %.9g %.9g\n",
                k, 4 * o, com[k].x, com[k].y, com[k].z, ref_com[k].x, ref_com[k].y, ref_com[k].z);
            EXPECT_EQ(0, memcmp(&ref_com[k], &com[k], sizeof(md_float3_t)));
        }
        EXPECT_EQ(0, memcmp(ref_min, mn, sizeof(mn)));
        EXPECT_EQ(0, memcmp(ref_max, mx, sizeof(mx)));
        EXPECT_EQ(0, memcmp(ref_pbc, pbc, sizeof(pbc)));
        EXPECT_EQ(0, memcmp(ref_lin, lin, sizeof(lin)));
    }
    md_free(md_get_heap_allocator(), raw, (6 * N + 64) * sizeof(float) + 64);
}

// Lengths and angles in, lengths and angles out: a triclinic cell whose alpha is not 90 among them
UTEST(geom, unitcell_extent_angles_round_trip) {
    const double cells[][6] = {
        { 35.446, 35.062, 34.159, 91.328, 61.735, 44.407 },
        { 40.0, 40.0, 40.0, 109.4712206, 109.4712206, 109.4712206 },   // truncated octahedron
        { 30.0, 31.0, 32.0, 60.0, 70.0, 80.0 },
        { 30.0, 31.0, 32.0, 90.0, 90.0, 90.0 },
    };
    for (size_t i = 0; i < MD_COUNTOF(cells); ++i) {
        const double* p = cells[i];
        const md_unitcell_t uc = md_unitcell_from_extent_and_angles(p[0], p[1], p[2], p[3], p[4], p[5]);
        double v[6];
        md_unitcell_extract_extent_angles(&v[0], &v[1], &v[2], &v[3], &v[4], &v[5], &uc);
        for (int k = 0; k < 6; ++k) EXPECT_NEAR(p[k], v[k], 1e-6);
    }
    // No cell: no angles, and no NaN
    const md_unitcell_t none = md_unitcell_none();
    double v[6];
    md_unitcell_extract_extent_angles(&v[0], &v[1], &v[2], &v[3], &v[4], &v[5], &none);
    EXPECT_EQ(0.0, v[3]);
}
