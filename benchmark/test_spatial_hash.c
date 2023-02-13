#include "ubench.h"

#include <md_gro.h>
#include <md_molecule.h>
#include <core/md_spatial_hash.h>
#include <core/md_arena_allocator.h>

struct spatial_hash {
    md_vm_arena_t arena;
    md_allocator_i alloc;
    md_molecule_t mol;
    vec3_t pbc_ext;
};

UBENCH_F_SETUP(spatial_hash) {
    md_vm_arena_init(&ubench_fixture->arena, GIGABYTES(4));

    ubench_fixture->alloc = md_vm_arena_create_interface(&ubench_fixture->arena);

    md_gro_data_t gro_data = {0};
    md_gro_data_parse_file(&gro_data, STR(MD_BENCHMARK_DATA_DIR "/centered.gro"), &ubench_fixture->alloc);
    md_gro_molecule_init(&ubench_fixture->mol, &gro_data, &ubench_fixture->alloc);
    ubench_fixture->pbc_ext = mat3_mul_vec3(ubench_fixture->mol.cell.basis, vec3_set1(1.f));
}

UBENCH_F_TEARDOWN(spatial_hash) {
    md_vm_arena_free(&ubench_fixture->arena);
}

UBENCH_EX_F(spatial_hash, perf_test_init) {
    md_molecule_t* mol = &ubench_fixture->mol;
    md_allocator_i* alloc = &ubench_fixture->alloc;
    vec3_t pbc_ext = ubench_fixture->pbc_ext;

    md_spatial_hash_t spatial_hash = {0};

    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(&ubench_fixture->arena);
    UBENCH_DO_BENCHMARK() {    
        md_spatial_hash_init_soa(&spatial_hash, mol->atom.x, mol->atom.y, mol->atom.z, mol->atom.count, pbc_ext, alloc);
    }

    md_vm_arena_temp_end(temp);
}