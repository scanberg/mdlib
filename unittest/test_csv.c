#include "utest.h"

#include <md_csv.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_str.h>
#include <core/md_writer.h>

#include <stdio.h>
#include <math.h>

static const float col_time[]     = { 0.0f, 0.5f, 1.0f };
static const float col_distance[] = { 1.25f, -2.5f, 3.0f };

UTEST(csv, write_and_parse) {
    md_allocator_i* heap  = md_get_heap_allocator();
    md_allocator_i* arena = md_arena_allocator_create(heap, MD_KILOBYTES(64));

    const float* const values[] = { col_time, col_distance };
    const md_str_t names[] = { MD_STR_LIT("time"), MD_STR_LIT("Distance (nm)") };

    // To any writer
    char buf[256];
    md_buf_writer_t bw = md_buf_writer_create(buf, sizeof(buf));
    md_csv_write(&bw.w, values, names, 2, 3);
    EXPECT_FALSE(md_buf_writer_truncated(&bw));
    EXPECT_STREQ(
        "time,Distance (nm)\n"
        "0.000000,1.250000\n"
        "0.500000,-2.500000\n"
        "1.000000,3.000000\n", md_buf_writer_to_str(&bw).ptr);

    // To a string, the same text
    md_str_t str = md_csv_write_to_str(values, names, 2, 3, arena);
    EXPECT_TRUE(md_str_eq(md_buf_writer_to_str(&bw), str));

    md_csv_t csv = {0};
    ASSERT_TRUE(md_csv_parse_str(&csv, str, arena));
    EXPECT_EQ((size_t)2, csv.num_fields);
    EXPECT_EQ((size_t)3, csv.num_values);
    if (csv.num_fields == 2 && csv.num_values == 3) {
        EXPECT_TRUE(md_str_eq(names[1], csv.field_names[1]));
        for (size_t i = 0; i < 3; ++i) {
            EXPECT_EQ(col_time[i],     csv.field_values[0][i]);
            EXPECT_EQ(col_distance[i], csv.field_values[1][i]);
        }
    }

    // To a file, the same text
    const md_str_t path = MD_STR_LIT("md_csv_test.tmp");
    ASSERT_TRUE(md_csv_write_to_file(values, names, 2, 3, path));
    md_str_t file_str = md_load_textfile(path, arena);
    EXPECT_TRUE(md_str_eq(str, file_str));
    remove(path.ptr);

    // Nothing to write, nothing written
    EXPECT_TRUE(md_str_empty(md_csv_write_to_str(values, names, 2, 0, arena)));

    md_arena_allocator_destroy(arena);
}

// The CSV that spreadsheets, pandas and hand edits write: no header row (the first row is data, not names), ';' or tab
// between the fields, quoted names, CRLF, blank and '#' lines, and empty fields (NaN, where the value is missing).
UTEST(csv, the_csv_that_is_written) {
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_KILOBYTES(64));
    {
        md_csv_t csv = {0};
        ASSERT_TRUE(md_csv_parse_str(&csv, MD_STR_LIT("0,1.5\n1,1.6\n2,1.7\n"), arena));
        EXPECT_EQ((size_t)2, csv.num_fields);
        ASSERT_EQ((size_t)3, csv.num_values);
        EXPECT_EQ(0, (int)md_array_size(csv.field_names));
        EXPECT_EQ(1.5f, csv.field_values[1][0]);
    }
    {
        md_csv_t csv = {0};
        ASSERT_TRUE(md_csv_parse_str(&csv, MD_STR_LIT("\"time\";\"rmsd (A)\"\r\n# a comment\r\n0;1.5\r\n\r\n1;1.6\r\n"), arena));
        ASSERT_EQ((size_t)2, csv.num_fields);
        ASSERT_EQ((size_t)2, csv.num_values);
        EXPECT_TRUE(md_str_eq(csv.field_names[1], MD_STR_LIT("rmsd (A)")));
        EXPECT_EQ(1.6f, csv.field_values[1][1]);
    }
    {
        md_csv_t csv = {0};
        ASSERT_TRUE(md_csv_parse_str(&csv, MD_STR_LIT("time\trmsd\n0\t1.5\n1\t1.6\n"), arena));
        EXPECT_EQ((size_t)2, csv.num_fields);
        EXPECT_EQ((size_t)2, csv.num_values);
    }
    {
        md_csv_t csv = {0};
        ASSERT_TRUE(md_csv_parse_str(&csv, MD_STR_LIT("time,rmsd,rg\n0,1.5,10\n1,,10.1\n2,1.7,\n"), arena));
        ASSERT_EQ((size_t)3, csv.num_values);
        EXPECT_TRUE(isnan(csv.field_values[1][1]));
        EXPECT_TRUE(isnan(csv.field_values[2][2]));
        EXPECT_EQ(1.7f, csv.field_values[1][2]);
    }
    {
        // A row of more fields than the table has is not this table
        md_csv_t csv = {0};
        EXPECT_FALSE(md_csv_parse_str(&csv, MD_STR_LIT("time,rmsd\n0,1.5\n1,1.6,7\n"), arena));
    }
    md_arena_allocator_destroy(arena);
}
