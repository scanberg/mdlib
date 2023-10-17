#include "ubench.h"

#include <core/md_os.h>
#include <core/md_arena_allocator.h>
#include <md_gro.h>
#include <md_util.h>

UBENCH_EX(util, com_soa) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(4));
    str_t path = STR(MD_BENCHMARK_DATA_DIR "/centered.gro");

    md_molecule_t mol = {0};
    md_gro_molecule_api()->init_from_file(&mol, path, alloc);

    vec3_t com = {0};
    vec3_t box = mat3_mul_vec3(mol.unit_cell.basis, vec3_set(1.0f, 1.0f, 1.0f));
    
    UBENCH_DO_BENCHMARK() {
        com = md_util_compute_com_ortho(mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.mass, 0, mol.atom.count, box);
    }

    printf("com: %.3f %.3f %.3f\n", com.x, com.y, com.z);

    md_arena_allocator_destroy(alloc);
}

UBENCH_EX(util, com_soa_indices) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(4));
    str_t path = STR(MD_BENCHMARK_DATA_DIR "/centered.gro");

    md_molecule_t mol = {0};
    md_gro_molecule_api()->init_from_file(&mol, path, alloc);

    int* indices = md_alloc(alloc, sizeof(int) * mol.atom.count);

    // Scramble indices to some extent
    int i = 0;
    int count = (int)ALIGN_TO(mol.atom.count, md_simd_width_f32) - md_simd_width_f32;
    for (; i < (int)count; i += 8) {
        indices[i + 0] = i + 4;
        indices[i + 1] = i + 0;
        indices[i + 2] = i + 7;
        indices[i + 3] = i + 5;
        indices[i + 4] = i + 6;
        indices[i + 5] = i + 2;
        indices[i + 6] = i + 3;
        indices[i + 7] = i + 1;
    }
    for (; i < mol.atom.count; ++i) {
        indices[i] = (int)i;
    }

    vec3_t com = {0};
    vec3_t box = mat3_mul_vec3(mol.unit_cell.basis, vec3_set(1.0f, 1.0f, 1.0f));

    UBENCH_DO_BENCHMARK() {
        com = md_util_compute_com_ortho(mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.mass, indices, mol.atom.count, box);        
    }

    printf("com: %.3f %.3f %.3f\n", com.x, com.y, com.z);

    md_arena_allocator_destroy(alloc);
}

UBENCH_EX(util, com_vec4) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(4));
    str_t path = STR(MD_BENCHMARK_DATA_DIR "/centered.gro");

    md_molecule_t mol = {0};
    md_gro_molecule_api()->init_from_file(&mol, path, alloc);

    vec4_t* xyzw = md_alloc(alloc, sizeof(vec4_t) * mol.atom.count);
    for (int64_t i = 0; i < mol.atom.count; ++i) {
        xyzw[i] = vec4_set(mol.atom.x[i], mol.atom.y[i], mol.atom.z[i], mol.atom.mass ? mol.atom.mass[i] : 1.0f);
    }

    vec3_t com = {0};
    vec3_t box = mat3_mul_vec3(mol.unit_cell.basis, vec3_set(1.0f, 1.0f, 1.0f));

    UBENCH_DO_BENCHMARK() {
        com = md_util_compute_com_vec4_ortho(xyzw, 0, mol.atom.count, box);
    }

    printf("com: %.3f %.3f %.3f\n", com.x, com.y, com.z);

    md_arena_allocator_destroy(alloc);
}

UBENCH_EX(util, com_vec4_indices) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(4));
    str_t path = STR(MD_BENCHMARK_DATA_DIR "/centered.gro");

    md_molecule_t mol = {0};
    md_gro_molecule_api()->init_from_file(&mol, path, alloc);

    int* indices = md_alloc(alloc, sizeof(int) * mol.atom.count);

    // Scramble indices to some extent
    int i = 0;
    int count = (int)ALIGN_TO(mol.atom.count, md_simd_width_f32) - md_simd_width_f32;
    for (; i < count; i += 8) {
        indices[i + 0] = i + 4;
        indices[i + 1] = i + 0;
        indices[i + 2] = i + 7;
        indices[i + 3] = i + 5;
        indices[i + 4] = i + 6;
        indices[i + 5] = i + 2;
        indices[i + 6] = i + 3;
        indices[i + 7] = i + 1;
    }
    for (; i < mol.atom.count; ++i) {
        indices[i] = (int)i;
    }

    vec4_t* xyzw = md_alloc(alloc, sizeof(vec3_t) * mol.atom.count);
    for (i = 0; i < mol.atom.count; ++i) {
        xyzw[i] = vec4_set(mol.atom.x[i], mol.atom.y[i], mol.atom.z[i], mol.atom.mass ? mol.atom.mass[i] : 1.0f);
    }

    vec3_t com = {0};
    vec3_t box = mat3_mul_vec3(mol.unit_cell.basis, vec3_set(1.0f, 1.0f, 1.0f));

    UBENCH_DO_BENCHMARK() {
        com = md_util_compute_com_vec4_ortho(xyzw, indices, mol.atom.count, box);
    }

    printf("com: %.3f %.3f %.3f\n", com.x, com.y, com.z);

    md_arena_allocator_destroy(alloc);
}