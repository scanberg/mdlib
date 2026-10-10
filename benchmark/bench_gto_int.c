#include "ubench.h"

#include <md_vlx.h>
#include <md_system.h>
#include <md_gto.h>
#include <md_gto_int.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_vec_math.h>
#if MD_ENABLE_GPU
#include <core/md_gpu.h>
#endif

#include <float.h>
#include <math.h>
#include <stdio.h>

// Electrostatic potential of test_data/vlx/mol.h5 (26 atoms, def2-SVP, B3LYP ground state): the
// total density and the nuclei, on a 0.4 bohr grid reaching 6 bohr past the atoms - the resolution
// an isodensity surface needs for its colouring (trilinear error ~1% of the colour range).

#define BENCH_ANGSTROM_TO_BOHR 1.8897261246257702

typedef struct {
    md_allocator_i* arena;
    md_gto_basis_t  basis;
    float*          atom_xyz;   // bohr
    double*         atom_z;
    size_t          num_atoms;
    double*         D;
    md_grid_t       grid;
} bench_esp_t;

static bool bench_esp_load(bench_esp_t* b) {
    MEMSET(b, 0, sizeof(*b));
    b->arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(64));
    md_system_t sys = { .alloc = b->arena };
    md_system_state_t state = { .alloc = b->arena };
    if (!md_vlx_system_init_from_file(&sys, &state, STR_LIT(MD_BENCHMARK_DATA_DIR "/vlx/mol.h5"))) {
        MD_LOG_ERROR("Could not load benchmark vlx file");
        return false;
    }
    const md_attribute_t* c = md_attributes_find(&sys.attributes, STR_LIT("qm/atom/coordinate"));
    const md_attribute_t* z = md_attributes_find(&sys.attributes, STR_LIT("qm/atom/nuclear_charge"));
    const md_attribute_t* d = md_attributes_find(&sys.attributes, STR_LIT("orbital/total/density"));
    if (!d) d = md_attributes_find(&sys.attributes, STR_LIT("orbital/alpha/density"));
    if (!c || !z || !d || !md_gto_basis_extract_attributes(&b->basis, &sys.attributes, b->arena)) {
        MD_LOG_ERROR("Benchmark file lacks a basis, a geometry or a density");
        return false;
    }
    b->num_atoms = md_attribute_value_count(&c->format);
    double* xyz = md_arena_allocator_push(b->arena, sizeof(double) * 3 * b->num_atoms);
    md_attribute_extract_f64(xyz, 3 * b->num_atoms, c, md_attribute_slice_all(), md_unit_none());
    b->atom_xyz = md_arena_allocator_push(b->arena, sizeof(float) * 3 * b->num_atoms);
    b->atom_z   = md_arena_allocator_push(b->arena, sizeof(double) * b->num_atoms);
    md_attribute_extract_f64(b->atom_z, b->num_atoms, z, md_attribute_slice_all(), md_unit_none());

    vec3_t lo = vec3_set1(FLT_MAX), hi = vec3_set1(-FLT_MAX);
    for (size_t i = 0; i < b->num_atoms; ++i) {
        const vec3_t p = vec3_set((float)(xyz[3 * i + 0] * BENCH_ANGSTROM_TO_BOHR), (float)(xyz[3 * i + 1] * BENCH_ANGSTROM_TO_BOHR), (float)(xyz[3 * i + 2] * BENCH_ANGSTROM_TO_BOHR));
        MEMCPY(b->atom_xyz + 3 * i, p.elem, sizeof(float) * 3);
        lo = vec3_min(lo, p);
        hi = vec3_max(hi, p);
    }
    const size_t n = md_gto_basis_num_ao(&b->basis);
    b->D = md_arena_allocator_push(b->arena, sizeof(double) * n * n);
    md_attribute_extract_f64(b->D, n * n, d, md_attribute_slice_all(), md_unit_none());

    const float h = 0.4f, margin = 6.0f;
    lo = vec3_sub1(lo, margin);
    hi = vec3_add1(hi, margin);
    b->grid = (md_grid_t){
        .orientation = mat3_ident(),
        .origin  = lo,
        .spacing = vec3_set1(h),
        .dim     = { (int)ceilf((hi.x - lo.x) / h), (int)ceilf((hi.y - lo.y) / h), (int)ceilf((hi.z - lo.z) / h) },
    };
    return true;
}

static bool bench_esp_charges(md_gto_int_charges_t* q, const bench_esp_t* b, double threshold, md_allocator_i* alloc) {
    md_gto_int_charges_desc_t desc = {
        .basis = &b->basis, .atom_xyz = b->atom_xyz, .density_matrix = b->D, .density_scale = -1.0,
        .point_xyz = b->atom_xyz, .point_charge = b->atom_z, .num_points = b->num_atoms,
        .threshold = threshold,
    };
    return md_gto_int_charges_init(q, &desc, alloc);
}

UBENCH_EX(gto_int, charges_build_mol) {
    bench_esp_t b;
    if (!bench_esp_load(&b)) return;
    md_gto_int_charges_t q = {0};
    UBENCH_DO_BENCHMARK() {
        bench_esp_charges(&q, &b, 1e-8, md_get_heap_allocator());
        md_gto_int_charges_free(&q, md_get_heap_allocator());
    }
    md_arena_allocator_destroy(b.arena);
}

// The CPU reference on 1024 voxels of the grid's middle - a slice of what a worker thread gets
// when a surface band is evaluated on the CPU. Per voxel it is one evaluation of every gaussian.
UBENCH_EX(gto_int, potential_cpu_mol_1024) {
    bench_esp_t b;
    if (!bench_esp_load(&b)) return;
    md_gto_int_charges_t q = {0};
    bench_esp_charges(&q, &b, 1e-8, b.arena);
    const int len[3] = { 16, 16, 4 };
    const int off[3] = { (b.grid.dim[0] - len[0]) / 2, (b.grid.dim[1] - len[1]) / 2, (b.grid.dim[2] - len[2]) / 2 };
    float* out = md_arena_allocator_push(b.arena, sizeof(float) * md_grid_num_points(&b.grid));
    const float so[3] = { 0.5f, 0.5f, 0.5f };
    UBENCH_DO_BENCHMARK() {
        md_gto_int_potential_grid_sub(out, &b.grid, so, off, len, &q);
    }
    md_arena_allocator_destroy(b.arena);
}

#if MD_ENABLE_GPU
UBENCH_EX(gto_int, potential_gpu_mol) {
    bench_esp_t b;
    if (!bench_esp_load(&b)) return;
    md_gpu_device_t device = md_gpu_device_create(NULL);
    if (!device) {
        md_arena_allocator_destroy(b.arena);
        return;
    }
    md_gpu_stream_t stream = md_gpu_stream_default(device, MD_GPU_STREAM_COMPUTE);
    md_gto_int_gpu_initialize(device);

    md_gto_int_charges_t q = {0};
    bench_esp_charges(&q, &b, 1e-8, b.arena);
    {
        uint64_t per_voxel = 0;
        for (uint32_t L = 0; L <= q.max_order; ++L) {
            const uint64_t n = q.order_offset[L + 1] - q.order_offset[L];
            per_voxel += n * (((L + 1) * (L + 2) * (L + 3) * (L + 4)) / 24 + md_gto_int_num_hermite(L));
            printf("  L%u: %llu", L, (unsigned long long)n);
        }
        printf("\n  mol.h5: %zu atoms, %zu AOs, %u gaussians, grid %dx%dx%d (%zu voxels), %.2e recursion steps\n",
               b.num_atoms, md_gto_basis_num_ao(&b.basis), q.num_gaussians, b.grid.dim[0], b.grid.dim[1], b.grid.dim[2],
               md_grid_num_points(&b.grid), (double)per_voxel * (double)md_grid_num_points(&b.grid));
    }

    md_gto_int_gpu_charges_t gq = md_gto_int_gpu_charges_create(stream, &q);
    md_gpu_texture_t tex = md_gpu_texture_create(stream, &(md_gpu_texture_desc_t){
        .type = MD_GPU_TEX_3D, .format = MD_GPU_FORMAT_R32_FLOAT, .usage = MD_GPU_TEX_STORAGE,
        .width = (uint32_t)b.grid.dim[0], .height = (uint32_t)b.grid.dim[1], .depth_or_layers = (uint32_t)b.grid.dim[2],
    });
    md_gto_int_gpu_potential_desc_t desc = {
        .charges = gq, .out_tex = tex, .grid = &b.grid, .sample_offset = {0.5f, 0.5f, 0.5f}, .op = MD_GTO_OP_SET,
    };
    // Warm up: kernel creation, first-use compilation
    md_gto_int_gpu_potential_launch(stream, &desc);
    md_gpu_stream_sync(stream);

    UBENCH_DO_BENCHMARK() {
        md_gto_int_gpu_potential_launch(stream, &desc);
        md_gpu_stream_sync(stream);
    }

    md_gpu_texture_destroy(tex);
    md_gto_int_gpu_charges_destroy(stream, gq);
    md_gpu_stream_sync(stream);
    md_gto_int_gpu_shutdown();
    md_gpu_device_destroy(device);
    md_arena_allocator_destroy(b.arena);
}
#endif
