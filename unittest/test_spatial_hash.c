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
    ASSERT_TRUE(md_spatial_hash_init(&spatial_hash, x, y, z, 10, 0, default_allocator));

    vec4_t pos_rad = {5, 0, 0, 1.5f};
    uint32_t count = 0;
    EXPECT_TRUE(md_spatial_hash_query(&spatial_hash, pos_rad, func, &count));

    EXPECT_EQ(3, count);
}

UTEST(spatial_hash, big) {
    md_allocator_i* alloc = md_arena_allocator_create(default_allocator, KILOBYTES(128));
    const str_t pdb_file = make_cstr(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");

    md_pdb_data_t pdb_data = {0};
    ASSERT_TRUE(md_pdb_data_parse_file(pdb_file, &pdb_data, alloc));

    md_molecule_t mol = {0};
    ASSERT_TRUE(md_pdb_molecule_init(&mol, &pdb_data, alloc));

    md_spatial_hash_t spatial_hash = {0};
    ASSERT_TRUE(md_spatial_hash_init(&spatial_hash, mol.atom.x, mol.atom.y, mol.atom.z, mol.atom.count, 0, alloc));

    vec4_t pos_rad = {24, 48, 24, 10.0f};
    uint32_t count = 0;
    EXPECT_TRUE(md_spatial_hash_query(&spatial_hash, pos_rad, func, &count));

    EXPECT_LT(0, count);

    md_arena_allocator_destroy(alloc);
}
