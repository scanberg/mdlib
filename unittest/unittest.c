#include "utest.h"
#include "md_allocator.h"
#include "md_log.h"
#include "md_molecule.h"
#include "md_filter.h"
#include "md_script.h"
#include "md_pdb.h"
#include "md_gro.h"
#include "md_xtc.h"

#include "core/common.h"
#include "core/bitop.h"
#include "core/strpool.h"
#include "core/file.h"
#include "core/str_util.h"

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

UTEST_MAIN();

UTEST(filter, test) {
    uint32_t atom_count = 8;
    float x[] = {1,2,3,4,5,6,7,8};
    float y[] = {1,1,1,1,1,1,1,1};
    float z[] = {2,2,2,2,2,2,2,2};
    float r[] = {1,2,3,4,4,4,5,1};
    float m[] = {1,2,2,2,2,4,4,4};
    uint8_t e[] = {1,8,1,8,5,8,8,8};
    char* n[] = {"H", "O", "H", "He", "C", "N", "CA", "O"};
    float b[] = {0,0,0,0,0,1,0,0};
    float o[] = {1,1,1,1,1,2,2,2};
    uint32_t r_idx[] = {0,0,0,1,1,1,1,1};
    uint32_t c_idx[] = {0,0,0,0,0,0,0,0};

    uint32_t res_count = 2;
    char* r_name[] = {"SOL", "LYS"};
    uint32_t r_id[] = {1, 2};
    md_range r_range[] = {{0, 3}, {3, 8}};

    uint32_t chain_count = 1;
    char* c_id[] = {"A"};
    md_range c_range[] = {0,8};

    md_molecule mol = {0};

    mol.atom.count = atom_count;
    mol.atom.x = x;
    mol.atom.y = y;
    mol.atom.z = z;
    mol.atom.radius = r;
    mol.atom.mass = m;
    mol.atom.element = e;
    mol.atom.name = n;
    mol.atom.bfactor = b;
    mol.atom.occupancy = o;
    mol.atom.residue_idx = r_idx;
    mol.atom.chain_idx = c_idx;

    mol.residue.count = res_count;
    mol.residue.name = r_name;
    mol.residue.id = r_id;
    mol.residue.atom_range = r_range;

    mol.chain.count = chain_count;
    mol.chain.id = c_id;
    mol.chain.atom_range = c_range;

    uint64_t bits = 0;
    uint64_t ref = 0;

    uint64_t stored_bits = make_bits("00000001");
    md_filter_stored_selection stored;
    stored.ident = "cool";
    stored.bits = &stored_bits;

    md_filter_context ctx = {
        &mol,
        mol.atom.count,
        &bits,
        1,
        &stored
    };

    {
        bits = 0;
        filter_func_all(&bits, &mol);
        ref = make_bits("11111111");
        EXPECT_EQ(bits, ref);
    }

    {
        bits = 0;
        filter_func_none(&bits, &mol);
        ref = make_bits("00000000");
        EXPECT_EQ(bits, ref);
    }

    {
        bits = 0;
        filter_func_protein(&bits, &mol);
        ref = make_bits("00011111");
        EXPECT_EQ(bits, ref);
    }

    {
        bits = 0;
        filter_func_water(&bits, &mol);
        ref = make_bits("11100000");
        EXPECT_EQ(bits, ref);
    }

    {
        bits = 0;
        filter_func_name("H", &bits, &mol);
        ref = make_bits("10100000");
        EXPECT_EQ(bits, ref);

        bits = 0;
        filter_func_name("CA", &bits, &mol);
        ref = make_bits("00000010");
        EXPECT_EQ(bits, ref);
    }

    {
        bits = 0;
        filter_func_resname("SOL", &bits, &mol);
        ref = make_bits("11100000");
        EXPECT_EQ(bits, ref);

        bits = 0;
        filter_func_resname("LYS", &bits, &mol);
        ref = make_bits("00011111");
        EXPECT_EQ(bits, ref);
    }

    {
        bits = 0;
        filter_func_resid(1, 1, &bits, &mol);
        ref = make_bits("11100000");
        EXPECT_EQ(bits, ref);

        bits = 0;
        filter_func_resid(2, 2, &bits, &mol);
        ref = make_bits("00011111");
        EXPECT_EQ(bits, ref);
    }

    {
        bits = 0;
        filter_func_residue(0, 0, &bits, &mol);
        ref = make_bits("11100000");
        EXPECT_EQ(bits, ref);

        bits = 0;
        filter_func_residue(1, 1, &bits, &mol);
        ref = make_bits("00011111");
        EXPECT_EQ(bits, ref);
    }

    {
        bits = 0;
        filter_func_chain("A", &bits, &mol);
        ref = make_bits("11111111");
        EXPECT_EQ(bits, ref);
    }

    {
        uint64_t mask = make_bits("10000000");
        bits = 0;
        filter_func_within(0, 3, &mask, &bits, &mol);
        ref = make_bits("11110000");
        EXPECT_EQ(bits, ref);
    }



    {
        bits = 0;
        ref = make_bits("00011111");
        EXPECT_TRUE(md_filter_apply("protein", &ctx));
        EXPECT_EQ(bits, ref);
    }

    {
        bits = 0;
        ref = make_bits("00011111");
        EXPECT_TRUE(md_filter_apply("resname(LYS)", &ctx));
        EXPECT_EQ(bits, ref);
    }

    {
        bits = 0;
        ref = make_bits("11100000");
        EXPECT_TRUE(md_filter_apply("water", &ctx));
        EXPECT_EQ(bits, ref);
    }

    {
        bits = 0;
        ref = make_bits("11100000");
        EXPECT_TRUE(md_filter_apply("resname(SOL)", &ctx));
        EXPECT_EQ(bits, ref);
    }


    {
        bits = 0;
        ref = make_bits("00000001");
        EXPECT_TRUE(md_filter_apply("@cool", &ctx));
        EXPECT_EQ(bits, ref);
    }

    {
        bits = 0;
        EXPECT_TRUE(md_filter_apply("resname(ALA LYS) and within(2.5, name(C H O) or residue(1:18)) and protein or @cool", &ctx));
        ref = make_bits("00011111");
        EXPECT_EQ(bits, ref);
    }
    {
#if 0
        const char* expr = "all and residue(1:48) or not name(H C O)";
        EXPECT_TRUE(md_filter_apply(expr, NULL));
#endif
    }
    {
        EXPECT_FALSE(md_filter_apply("all and not resjizzle(1 2 3 4)", NULL));
    }
    {
        EXPECT_FALSE(md_filter_apply("residue(1 ALA 8 1:8)", NULL));
    }
}

UTEST(strpool, test) {
    md_strpool pool;

    md_strpool_init(&pool, NULL);
    EXPECT_EQ(pool.page.len, 1);
    EXPECT_EQ(pool.entry.len, 0);
    EXPECT_GT(pool.entry.cap, 0);

    md_strpool_insert(&pool, "A", (int)strlen("A"));
    md_strpool_insert(&pool, "B", (int)strlen("B"));
    md_strpool_insert(&pool, "C", (int)strlen("C"));
    md_strpool_insert(&pool, "D", (int)strlen("D"));
    md_strpool_insert(&pool, "E", (int)strlen("E"));
    md_strpool_insert(&pool, "F", (int)strlen("F"));
    md_strpool_insert(&pool, "G", (int)strlen("G"));
    md_strpool_insert(&pool, "H", (int)strlen("H"));
    md_strpool_insert(&pool, "I", (int)strlen("I"));
    md_strpool_insert(&pool, "J", (int)strlen("J"));
    md_strpool_insert(&pool, "K", (int)strlen("K"));
    md_strpool_insert(&pool, "L", (int)strlen("L"));
    md_strpool_insert(&pool, "M", (int)strlen("M"));
    md_strpool_insert(&pool, "N", (int)strlen("N"));
    md_strpool_insert(&pool, "O", (int)strlen("O"));
    md_strpool_insert(&pool, "P", (int)strlen("P"));   // 16 -> entry grow
    md_strpool_insert(&pool, "Q", (int)strlen("Q"));
    md_strpool_insert(&pool, "R", (int)strlen("R"));
    md_strpool_insert(&pool, "S", (int)strlen("S"));
    md_strpool_insert(&pool, "T", (int)strlen("T"));
    EXPECT_EQ(pool.entry.len, 20);

    md_strpool_insert(&pool, "A", (int)strlen("A"));
    md_strpool_insert(&pool, "B", (int)strlen("B"));
    md_strpool_insert(&pool, "C", (int)strlen("C"));
    md_strpool_insert(&pool, "D", (int)strlen("D"));
    md_strpool_insert(&pool, "Q", (int)strlen("Q"));
    EXPECT_EQ(pool.entry.len, 20);

    md_strpool_insert(&pool, "AA", (int)strlen("AA"));
    md_strpool_insert(&pool, "AB", (int)strlen("AB"));
    md_strpool_insert(&pool, "AC", (int)strlen("AC"));
    md_strpool_insert(&pool, "AD", (int)strlen("AD"));
    md_strpool_insert(&pool, "AE", (int)strlen("AE"));
    EXPECT_EQ(pool.entry.len, 25);

    EXPECT_EQ(pool.page.len, 1);

    char str[] = "this is a very long string which is not meant for such a system but will work since it is less than the page limit of 1024 characters but will work to populate the pool nevertheless so we can trigger page allocations. It's triggered when the page is full.";
    md_strpool_insert(&pool, str, (int)strlen(str));
    md_strpool_insert(&pool, str, (int)strlen(str)-1);
    md_strpool_insert(&pool, str, (int)strlen(str)-2);
    md_strpool_insert(&pool, str, (int)strlen(str)-3);
    EXPECT_EQ(pool.entry.len, 29)
    EXPECT_EQ(pool.page.len, 2);

    md_strpool_clear(&pool);
    EXPECT_EQ(pool.entry.len, 0);
    EXPECT_GT(pool.entry.cap, 0);
    EXPECT_EQ(pool.page.len, 1);

    const char filename[] = MD_UNITTEST_DATA_DIR "/atom_res_chain.txt";
    md_file* file = md_file_open(filename, ARRAY_SIZE(filename), "rb");
    EXPECT_TRUE(file);

    char buf[64];
    int lines_read = 0;
    while (fgets(buf, 64, (FILE*)file)) {
        const char* atom = buf;
        int atom_len = 0;
        while (buf[atom_len] != ' ' && atom_len < 64) ++atom_len;

        const char* res = buf + 4;
        int res_len = 3;

        const char* chain = buf + 8;
        int chain_len = 1;

        md_strpool_insert(&pool, atom,  atom_len);
        md_strpool_insert(&pool, res,   res_len);
        md_strpool_insert(&pool, chain, chain_len);
        lines_read += 1;
    }

    md_file_close(file);

    EXPECT_LT((int)pool.entry.len, lines_read * 3);
}

UTEST(script, test) {
    const char script_file[] = "script.txt";
    str_t script_text = load_textfile(script_file, ARRAY_SIZE(script_file), default_allocator);
    // @TODO: Don't leak here, fix a proper free function for str_t
    struct md_script_ir* ir = md_script_compile(script_text.str, (uint32_t)script_text.len, NULL);

    md_print(MD_LOG_TYPE_INFO, "cool");
}

UTEST(allocator, test) {
    void* ptr;
    for (int i = 0; i < 1000000; ++i) {
        ptr = md_alloc(default_temp_allocator, 1000);
    }
}