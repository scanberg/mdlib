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
    
    UBENCH_DO_BENCHMARK() {
        vec3_t com = md_util_compute_com(mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.mass, 0, mol.atom.count);
        UBENCH_DO_NOTHING(&com);
    }

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
    for (i = 0; i < (int)mol.atom.count; i += 8) {
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

    UBENCH_DO_BENCHMARK() {
        vec3_t com = md_util_compute_com(mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.mass, indices, mol.atom.count);
        UBENCH_DO_NOTHING(&com);
    }

    md_arena_allocator_destroy(alloc);
}

UBENCH_EX(util, com_vec3) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(4));
    str_t path = STR(MD_BENCHMARK_DATA_DIR "/centered.gro");

    md_molecule_t mol = {0};
    md_gro_molecule_api()->init_from_file(&mol, path, alloc);

    vec3_t* pos = md_alloc(alloc, sizeof(vec3_t) * mol.atom.count);
    for (int64_t i = 0; i < mol.atom.count; ++i) {
        pos[i] = vec3_set(mol.atom.x[i], mol.atom.y[i], mol.atom.z[i]);
    }

    UBENCH_DO_BENCHMARK() {
        vec3_t com = md_util_compute_com_vec3(pos, mol.atom.mass, 0, mol.atom.count);
        UBENCH_DO_NOTHING(&com);
    }

    md_arena_allocator_destroy(alloc);
}

UBENCH_EX(util, com_vec3_indices) {
    md_allocator_i* alloc = md_arena_allocator_create(md_heap_allocator, MEGABYTES(4));
    str_t path = STR(MD_BENCHMARK_DATA_DIR "/centered.gro");

    md_molecule_t mol = {0};
    md_gro_molecule_api()->init_from_file(&mol, path, alloc);

    int* indices = md_alloc(alloc, sizeof(int) * mol.atom.count);

    // Scramble indices to some extent
    int i = 0;
    for (i = 0; i < (int)mol.atom.count; i += 8) {
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

    vec3_t* pos = md_alloc(alloc, sizeof(vec3_t) * mol.atom.count);
    for (i = 0; i < mol.atom.count; ++i) {
        pos[i] = vec3_set(mol.atom.x[i], mol.atom.y[i], mol.atom.z[i]);
    }

    UBENCH_DO_BENCHMARK() {
        vec3_t com = md_util_compute_com_vec3(pos, mol.atom.mass, indices, mol.atom.count);
        UBENCH_DO_NOTHING(&com);
    }

    md_arena_allocator_destroy(alloc);
}