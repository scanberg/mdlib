#ifndef _MD_SCRIPT_H_
#define _MD_SCRIPT_H_

#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_molecule;
struct md_trajectory;

struct md_script_ir;


// COMPILE
struct md_script_error {
    uint32_t line;
    uint32_t offset;
    uint32_t length;
    const char* error_str;
    uint32_t    error_len;
};

/*
struct md_script_visualization_primitives {
    struct md_allocator_i* _alloc;

    uint32_t num_points;
    float (*points)[3];

    uint32_t num_lines;
    float (*lines_beg)[3];
    float (*lines_end)[3];

    uint64_t num_bits;
    uint64_t* bits;
};

struct md_script_visualization_primitives* md_script_create_visualization_primitives(const char* str, uint32_t str_len, const struct md_script_context* ctx, struct md_allocator_i* alloc);
void md_script_free_visualization_primitives(struct md_script_visualization_primitives* primitives);
*/

struct md_script_ir* md_script_compile(const char* str, uint64_t len);
void md_script_free(struct md_script_ir* ir);

uint32_t md_script_get_num_errors(const struct md_script_ir* ir);
void md_script_get_errors(const struct md_script_ir* ir);

// execute
//struct md_script_res* md_script_execute(const struct md_script_ir* ir, struct md_allocator_i* alloc);

//void md_script_get_num_properties(const struct md_script_res* res);

#ifdef __cplusplus
}
#endif

#endif