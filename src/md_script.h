#ifndef _MD_SCRIPT_H_
#define _MD_SCRIPT_H_

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_molecule;
struct md_trajectory_i;
struct md_allocator_i;

typedef struct md_script_o md_script_o;
typedef struct md_script_error_t md_script_error_t;
typedef struct md_script_data_t md_script_data_t;
typedef struct md_script_property_t md_script_property_t;
typedef enum md_script_unit_t md_script_unit_t;

enum md_script_unit_t {
    MD_SCRIPT_UNIT_NONE,
    MD_SCRIPT_UNIT_ANGSTROM,
    MD_SCRIPT_UNIT_DEGREES,
};

// Opaque script object
struct md_script_o;

struct md_script_error_t {
    // Data for indicating where the error occured within the script string
    int32_t line;      // Line number
    int32_t offset;    // String pointer
    int32_t length;    // Length in characters

    // Human readable error message
    str_t error;
};

struct md_script_data_t {
    // This gives you the length of the data in each dimension
    int32_t dims[4];

    // Flat access of the data
    int64_t num_values;
    float*  values;
    float   min_value;
    float   max_value;

    float*  variance; // optional, only computed if the values are computed as an aggregate

    md_script_unit_t unit;
};

struct md_script_property_t {
    str_t ident;
    struct md_script_data_t data;
};

// ### COMPILE ###
// Should always return an object, even if the compilation is not successful, so you can extract the error messages
md_script_o* md_script_compile(str_t src, const struct md_molecule* mol, struct md_allocator_i* alloc);

bool md_script_compile_success(const md_script_o* ir);

// ### ERROR MESSAGES ###
int64_t md_script_get_num_errors(const md_script_o* ir);
const md_script_error_t* md_script_get_errors(const md_script_o* ir); // Retreives all messages

// ### EVALUATE ###
// This evaluates all frames
bool md_script_evaluate(md_script_o* ir, const struct md_molecule* mol, const struct md_trajectory_i* traj, struct md_allocator_i* alloc);

// This only evaluates the supplied subset of frames, represented in the bitmask
bool md_script_evaluate_subset(md_script_o* ir, const struct md_molecule* mol, const struct md_trajectory_i* traj, uint64_t num_frame_mask_bits, uint64_t* frame_mask, struct md_allocator_i* alloc);

// ### PROPERTIES ###
int64_t md_script_get_num_properties(const md_script_o* ir);
const md_script_property_t* md_script_get_properties(const md_script_o* ir);

void md_script_free(md_script_o* ir);

#ifdef __cplusplus
}
#endif

#endif