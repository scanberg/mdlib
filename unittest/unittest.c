#include "utest.h"
#include "md_allocator.h"
#include "md_log.h"
#include "md_molecule.h"
#include "md_filter.h"
#include "md_script.h"
#include "md_pdb.h"
#include "md_gro.h"
#include "md_xtc.h"

#include "core/intrinsics.h"
#include "core/array.inl"
#include "core/file.inl"
#include "core/common.h"
#include "core/bitop.h"
#include "core/str_util.h"
#include "core/pool_allocator.h"
#include "core/arena_allocator.h"

extern void filter_func_all     (uint64_t* bits, const md_molecule* mol);
extern void filter_func_none    (uint64_t* bits, const md_molecule* mol);
extern void filter_func_protein (uint64_t* bits, const md_molecule* mol);
extern void filter_func_water   (uint64_t* bits, const md_molecule* mol);
extern void filter_func_name    (const char* str,                uint64_t* bits, const md_molecule* mol);
extern void filter_func_resname (const char* str,                uint64_t* bits, const md_molecule* mol);
extern void filter_func_resid   (int min_range, int max_range,   uint64_t* bits, const md_molecule* mol);
extern void filter_func_residue (int min_range, int max_range,   uint64_t* bits, const md_molecule* mol);
extern void filter_func_chain   (const char* str,                uint64_t* bits, const md_molecule* mol);
extern void filter_func_within  (float min_range, float max_range, const uint64_t* in_bits, uint64_t* out_bits, const md_molecule* mol);

uint32_t str_len(const char* cstr) {
    if (!cstr) return 0;
    uint32_t len = 0;
    while(cstr[len] != '\0') len++;
    return len;
}

uint64_t make_bits(const char* str) {
    uint64_t bits = {0};
    uint64_t len = str_len(str) < 64 ? str_len(str) : 64;
    for (uint64_t i = 0; i < len; i++) {
        if (str[i] != '0') bits |= (1LLU << i);
    }
    return bits;
}

void print_bits(const uint64_t* bits, uint64_t bit_count) {
    for (uint64_t i = 0; i < bit_count; ++i) {
        printf("%i", (int)bit_test(bits, i));
    }
}

typedef struct AllocatorTestData {
    void* mem;
    uint64_t size;
} AllocatorTestData;

typedef struct linear_allocator {
    void* mem;
    uint64_t offset;
    uint64_t size;
} linear_allocator;

static inline void linear_allocator_init(linear_allocator* alloc, void* mem, uint64_t size) {
    alloc->mem = mem;
    alloc->offset = 0;
    alloc->size = size;
}

static inline void* linear_allocator_alloc(linear_allocator* alloc, uint64_t size) {
    ASSERT(alloc->offset + size <= alloc->size);
    void* res = (char*)alloc->mem + alloc->offset;
    alloc->offset += size;
    return res;
}

static inline void linear_allocator_reset(linear_allocator* alloc) {
    alloc->offset = 0;
}

static AllocatorTestData allocator_test_data;
static linear_allocator linear_alloc;

UTEST_STATE();

int main(int argc, const char *const argv[]) {
    allocator_test_data.mem = md_alloc(default_allocator, MEGABYTES(4));
    allocator_test_data.size = MEGABYTES(4);

    linear_allocator_init(&linear_alloc, allocator_test_data.mem, allocator_test_data.size);

    // do your own thing
    return utest_main(argc, argv);
}

UTEST(core, intrinsics) {
    {
        static const uint32_t mask[5] = {0, 0xFFU, 0xFFFFU, 0xFFFFFFU, 0xFFFFFFFFU};

        EXPECT_EQ(clz(mask[1]), 24);
        EXPECT_EQ(clz(mask[2]), 16);
        EXPECT_EQ(clz(mask[3]), 8);
        EXPECT_EQ(clz(mask[4]), 0);

        EXPECT_EQ(popcnt(mask[0]), 0);
        EXPECT_EQ(popcnt(mask[1]), 8);
        EXPECT_EQ(popcnt(mask[2]), 16);
        EXPECT_EQ(popcnt(mask[3]), 24);
        EXPECT_EQ(popcnt(mask[4]), 32);
    }

    {
        static const uint32_t mask[5] = {0, (1U << 8), (1U << 16), (1U << 24), (1U << 31)};

        EXPECT_EQ(bit_scan_forward(mask[0]), 0);
        EXPECT_EQ(bit_scan_forward(mask[1]), 9);
        EXPECT_EQ(bit_scan_forward(mask[2]), 17);
        EXPECT_EQ(bit_scan_forward(mask[3]), 25);
        EXPECT_EQ(bit_scan_forward(mask[4]), 32);
    }

    EXPECT_EQ(next_power_of_two(3), 4);
    EXPECT_EQ(next_power_of_two(31), 32);
    EXPECT_EQ(next_power_of_two(32), 32);
    EXPECT_EQ(next_power_of_two(5), 8);
    EXPECT_EQ(next_power_of_two(8000), 8192);

    EXPECT_EQ(next_power_of_two64(sizeof(void*)), sizeof(void*));

    // @TODO: Implement more tests for example find_first_zero_byte
}

UTEST(allocator_perf, temp) {
    void* ptr[1000];
    for (int j = 0; j < 1000; ++j) {
        for (int i = 0; i < 1000; ++i) {
            ptr[i] = md_alloc(default_temp_allocator, 4);
        }
        /*
        for (int i = 0; i < 1000; ++i) {
            md_free(default_temp_allocator, ptr[i], 4);
        }
        */
    }
}

UTEST(allocator_perf, linear_inline) {
    void* ptr[1000];
    for (int j = 0; j < 1000; ++j) {
        for (int i = 0; i < 1000; ++i) {
            ptr[i] = linear_allocator_alloc(&linear_alloc, 4);
        }
    }
}

UTEST(allocator_perf, default) {
    void* ptr[1000];
    for (int j = 0; j < 1000; ++j) {
        for (int i = 0; i < 1000; ++i) {
            ptr[i] = md_alloc(default_allocator, 4);
        }
        for (int i = 0; i < 1000; ++i) {
            md_free(default_allocator, ptr[i], 4);
        }
    }
}

UTEST(allocator_perf, arena) {
    md_allocator_i* arena = md_arena_allocator_create(default_allocator, MEGABYTES(1));
    void* ptr[1000];
    for (int j = 0; j < 1000; ++j) {
        for (int i = 0; i < 1000; ++i) {
            ptr[i] = md_alloc(arena, 4);
        }
        /*
        for (int i = 0; i < 1000; ++i) {
            md_free(arena, ptr[i], 4);
        }
        */
    }
    md_arena_allocator_destroy(arena);
}

UTEST(array, test) {
    uint64_t arr_data[] = {0, 5, 1, 2, 3, 4, 5};
    uint64_t *arr = arr_data + 2;

    EXPECT_EQ(md_array_size(arr), 5);
    EXPECT_EQ(md_array_capacity(arr), 0);
    EXPECT_EQ(md_array_bytes(arr), 5 * sizeof(uint64_t));
    EXPECT_EQ(md_array_end(arr), arr + 5);
    EXPECT_EQ(md_array_last(arr), arr + 4);

    // TODO: FILL IN
}

UTEST(pool_alloc, test) {
    md_allocator_i* pool = md_pool_allocator_create(default_allocator, sizeof(uint64_t));

    uint64_t **mem = {0};

    for (int j = 0; j < 1000; ++j) {
        for (int i = 0; i < 1000; ++i) {
            uint64_t *item = *md_array_push(mem, md_alloc(pool, sizeof(uint64_t)), default_allocator);
            *item = i;
        }

        for (int i = 100; i < 1000; ++i) {
            md_free(pool, mem[i], sizeof(uint64_t));
        }
        md_array_shrink(mem, 100);

        int indices[10] = {1, 2, 5, 6, 70, 90, 18, 16, 12, 10};

        for (int i = 0; i < ARRAY_SIZE(indices); ++i) {
            int idx = indices[i];
            md_free(pool, mem[idx], sizeof(uint64_t));
        }

        for (int i = 0; i < ARRAY_SIZE(indices); ++i) {
            int idx = indices[i];
            mem[idx] = md_alloc(pool, sizeof(uint64_t));
        }
    }

    md_array_free(mem, default_allocator);
}

UTEST(arena_alloc, test) {
    md_allocator_i* arena = md_arena_allocator_create(default_allocator, MD_ARENA_ALLOCATOR_DEFAULT_PAGE_SIZE);

    for (int j = 0; j < 1000; ++j) {
        EXPECT_EQ((uint64_t)md_alloc(arena, 16) % 16, 0); // Expect to be aligned to 16 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 4) % 16, 0);  // Expect to be aligned to 16 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 4) % 16, 0);  // Expect to be aligned to 16 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 4) % 16, 0);  // Expect to be aligned to 16 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 4) % 16, 0);  // Expect to be aligned to 16 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 4) % 16, 0);  // Expect to be aligned to 16 bytes
        EXPECT_EQ((uint64_t)md_alloc(arena, 1) % 1, 0);  // Expect to be aligned to 1 byte
        EXPECT_EQ((uint64_t)md_alloc(arena, 1) % 1, 0);  // Expect to be aligned to 1 byte
        EXPECT_EQ((uint64_t)md_alloc(arena, 1) % 1, 0);  // Expect to be aligned to 1 byte
        EXPECT_EQ((uint64_t)md_alloc(arena, 1) % 1, 0);  // Expect to be aligned to 1 byte
        EXPECT_EQ((uint64_t)md_alloc(arena, 1) % 1, 0);  // Expect to be aligned to 1 byte
        EXPECT_EQ((uint64_t)md_alloc(arena, 2) % 2, 0);  // Expect to be aligned to 2 bytes

        EXPECT_NE(md_alloc(arena, 9000), NULL); // Exceeds page size, should still be good

        // Make sure we get some internal pages going.
        for (int i = 0; i < 20; ++i) {
            EXPECT_EQ((uint64_t)md_alloc(arena, 1024) % 16, 0);
        }

        md_arena_allocator_reset(arena);
    }

    md_arena_allocator_destroy(arena);
}

UTEST(script, test) {

    // Create molecule for evaulation
    uint32_t atom_count = 16;
    float x[] = {1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8};
    float y[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    float z[] = {2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2};
    float r[] = {1,2,3,4,4,4,5,1,1,2,3,4,4,4,5,1};
    float m[] = {1,2,2,2,2,4,4,4,1,2,2,2,2,4,4,4};
    uint8_t e[] = {1,8,1,8,5,8,8,8,1,8,1,8,5,8,8,8};
    md_label n[] = {"H", "O", "H", "He", "C", "N", "CA", "O"};
    float b[] = {0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0};
    float o[] = {1,1,1,1,1,2,2,2,1,1,1,1,1,2,2,2};
    uint32_t r_idx[] = {0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,3};
    uint32_t c_idx[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    uint32_t res_count = 4;
    md_label r_name[] = {"SOL", "LYS", "PFT", "PFT"};
    uint32_t r_id[] = {1, 2, 3, 4};
    md_range r_range[] = {{0, 3}, {3, 8}, {8,12}, {12,16}};

    uint32_t chain_count = 1;
    md_label c_id[] = {"A"};
    md_range c_arange[] = {0,16};
    md_range c_rrange[] = {0,4};

    md_molecule mol = {
        .atom = {
            .count = atom_count,
            .x = x,
            .y = y,
            .z = z,
            .radius = r,
            .mass = m,
            .element = e,
            .name = n,
            .bfactor = b,
            .occupancy = o,
            .residue_idx = r_idx,
            .chain_idx = c_idx
        },
        .residue = {
            .count = res_count,
            .name = r_name,
            .id = r_id,
            .atom_range = r_range
        },
        .chain = {
            .count = chain_count,
            .id = c_id,
            .atom_range = c_arange,
            .residue_range = c_rrange
        }
    };

    const char script_file[] = MD_UNITTEST_DATA_DIR "/script.txt";
    str_t script_text = load_textfile(script_file, ARRAY_SIZE(script_file), default_allocator);

    struct md_script_ir* ir = md_script_compile(script_text.ptr, script_text.len, &mol, default_allocator);

    free_str(script_text, default_allocator);
}