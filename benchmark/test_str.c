#include "ubench.h"

#include <core/md_str.h>
#include <core/md_parse.h>
#include <core/md_os.h>
#include <core/md_allocator.h>

UBENCH_EX(str, read_lines) {
    md_file_o* file = md_file_open(STR(MD_BENCHMARK_DATA_DIR "/centered.gro"), MD_FILE_READ | MD_FILE_BINARY);

    UBENCH_SET_BYTES(md_file_size(file));

    const int64_t cap = MEGABYTES(1);
    char* buf = md_alloc(default_allocator, cap);

    UBENCH_DO_BENCHMARK() {
        md_file_seek(file, 0, MD_FILE_BEG);
        while (md_file_read_lines(file, buf, cap)) {
            // do nothing
        }
    }

    md_free(default_allocator, buf, cap);
}

UBENCH_EX(str, buffered_reader) {
    md_file_o* file = md_file_open(STR(MD_BENCHMARK_DATA_DIR "/centered.gro"), MD_FILE_READ | MD_FILE_BINARY);
    const int64_t cap = MEGABYTES(1);
    char* buf = md_alloc(default_allocator, cap);
    
    md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);

    UBENCH_SET_BYTES(md_file_size(file));

    UBENCH_DO_BENCHMARK() {
        md_file_seek(file, 0, MD_FILE_BEG);
        str_t line;
        while (md_buffered_line_reader_extract_line(&line, &reader)) {
            // do nothing
        }
    }

    md_free(default_allocator, buf, cap);
    md_file_close(file);
}

UBENCH_EX(str, parse_int) {
    str_t str[] = {
        STR("8123"),
        STR("1123"),
        STR("19228123"),
        STR("1921238123"),
        STR("1928123"),
        STR("19"),
        STR("11232323"),
        STR("-20"),
        STR("1928123123123"),
        STR("-11239"),
        STR("11232323"),
        STR("220"),
        STR("1928123"),
        STR("1231419"),
        STR("123123223"),
        STR("02"),
    };

    int64_t num_bytes = 0;
    for (int i = 0; i < ARRAY_SIZE(str); ++i) {
        num_bytes += str[i].len;
    }
    UBENCH_SET_BYTES(num_bytes);


    int64_t acc = 0;
    int i = 0;
    UBENCH_DO_BENCHMARK() {
		acc += str_parse_int_simd(str[0]);
        acc += str_parse_int_simd(str[1]);
        acc += str_parse_int_simd(str[2]);
        acc += str_parse_int_simd(str[3]);
        
        acc += str_parse_int_simd(str[4]);
        acc += str_parse_int_simd(str[5]);
        acc += str_parse_int_simd(str[6]);
        acc += str_parse_int_simd(str[7]);
        
        acc += str_parse_int_simd(str[8]);
        acc += str_parse_int_simd(str[9]);
        acc += str_parse_int_simd(str[10]);
        acc += str_parse_int_simd(str[11]);

        acc += str_parse_int_simd(str[12]);
        acc += str_parse_int_simd(str[13]);
        acc += str_parse_int_simd(str[14]);
        acc += str_parse_int_simd(str[15]);
    }

    printf("acc: %i\n", (int)acc);
}

UBENCH_EX(str, parse_float) {
    str_t str[] = {
        STR("1928123.2767"),
        STR("19.2"),
        STR("12323"),
        STR("0.000000"),
    };

    int64_t num_bytes = 0;
    for (int i = 0; i < (int)ARRAY_SIZE(str); ++i) {
        num_bytes += str[i].len;
    }
    UBENCH_SET_BYTES(num_bytes);

    double acc = 0;
    
    int i = 0;
    UBENCH_DO_BENCHMARK() {
        acc += parse_float(str[i++]);
        acc += parse_float(str[i++]);
        acc += parse_float(str[i++]);
        acc += parse_float(str[i++]);
        
        i = 0;
        
        acc += parse_float(str[i++]);
        acc += parse_float(str[i++]);
        acc += parse_float(str[i++]);
        acc += parse_float(str[i++]);

        i = 0;
    }

    //printf("acc: %f\n", acc);
}