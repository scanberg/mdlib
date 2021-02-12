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
typedef struct md_script_ir md_script_ir;


// COMPILE
typedef struct md_script_error {
    // Data for indicating where the error occured within the script string
    uint32_t line;      // Line number
    uint32_t offset;    // Offset in characters
    uint32_t length;    // Length in characters

    // Human readable error message
    const char* str_ptr;
    uint32_t    str_len;
} md_script_error;

typedef struct md_script_type_info {
    int cool;
} md_script_type_info;

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

// ### COMPILE THE SCRIPT ###

// Should always return something, even if the compilation is not successful, so you can get the error messages
md_script_ir* md_script_compile(const char* str, uint64_t len, struct md_allocator_i* alloc);
void md_script_free(md_script_ir* ir);

bool md_script_compile_success(const md_script_ir* ir);

// ### ERROR MESSAGES ###
uint32_t md_script_get_num_errors(const md_script_ir* ir);
const md_script_error* md_script_get_errors(const md_script_ir* ir); // Retreives all messages

// ### TYPE INFO ###
// For providing syntax highlighting in text editor and providing human readable context data.
uint32_t md_script_get_num_type_infos(const md_script_ir* ir);
const md_script_type_info* md_script_get_type_info(const md_script_ir* ir); // Retreives all type infos


// ### EXECUTE SCRIPT ###
//struct md_script_res* md_script_execute(const struct md_script_ir* ir, struct md_allocator_i* alloc);

//void md_script_get_num_properties(const struct md_script_res* res);

#ifdef __cplusplus
}
#endif

#endif