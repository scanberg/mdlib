#include "utest.h"

#include "core/md_allocator.h"
#include "core/md_arena_allocator.h"
#include "core/md_spatial_hash.h"
#include "core/md_str.h"
#include <md_pdb.h>
#include <md_molecule.h>

static bool func(uint32_t idx, vec3_t pos, void* param) {
    (void)idx;
    (void)pos;
    uint32_t* count = param;
    (*count)++;
    return true;
}

UTEST(spatial_hash, small) {
    float x[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    float y[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    float z[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    md_spatial_hash_t spatial_hash = {0};
    md_spatial_hash_args_t args = {
        .coords = {
            .count = 10,
            .x = x,
            .y = y,
            .z = z,
        },
        .alloc = default_allocator
    };
    ASSERT_TRUE(md_spatial_hash_init(&spatial_hash, &args));

    vec3_t pos = {5, 0, 0};
    float  rad = 1.5f;
    uint32_t count = 0;
    EXPECT_TRUE(md_spatial_hash_query(&spatial_hash, pos, rad, func, &count));

    EXPECT_EQ(3, count);
}

UTEST(spatial_hash, big) {
    md_allocator_i* alloc = md_arena_allocator_create(default_allocator, KILOBYTES(64));
    const str_t pdb_file = make_cstr(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");

    md_pdb_data_t pdb_data = {0};
    ASSERT_TRUE(md_pdb_data_parse_file(pdb_file, &pdb_data, alloc));

    md_molecule_t mol = {0};
    ASSERT_TRUE(md_pdb_molecule_init(&mol, &pdb_data, alloc));

    md_spatial_hash_t spatial_hash = {0};
    md_spatial_hash_args_t args = {
        .coords = {
            .count = mol.atom.count,
            .x = mol.atom.x,
            .y = mol.atom.y,
            .z = mol.atom.z,
        },
        .alloc = default_allocator
    };
    ASSERT_TRUE(md_spatial_hash_init(&spatial_hash, &args));

    vec3_t pos = {24, 48, 24};
    float  rad = 10.0f;
    uint32_t count = 0;
    EXPECT_TRUE(md_spatial_hash_query(&spatial_hash, pos, rad, func, &count));

    EXPECT_LT(0, count);

    md_arena_allocator_destroy(alloc);
}
