#pragma once

// ### FILTERS ###
//
// A filter is an expression of the script language that picks atoms: "protein and not resname('PRO')",
// "within(5, residue(1))", "s1 or chain('B')". It is read against a system, may use the names a compiled script
// defines (s1 above), and is evaluated in a state of the system. What it picks may depend on the coordinates of the
// state (within() does), and then it has to be evaluated again for another state (is_dynamic).
//
// Where something goes wrong, the result says what, and where in the expression, so that an input can point at it.
// While a filter is typed, md_filter_complete says what can be written at the cursor.

#include <core/md_array.h>
#include <core/md_str.h>
#include <md_script.h>

#include <stdint.h>
#include <stdbool.h>

struct md_allocator_i;
struct md_bitset_t;
struct md_system_t;
struct md_system_state_t;
struct md_script_ir_t;

#ifdef __cplusplus
    extern "C" {
#endif

typedef struct md_filter_context_t {
    const struct md_system_t*       sys;    // The system the atoms are picked from. Required.
    const struct md_system_state_t* state;  // The coordinates it is evaluated in. NULL for the system's reference state.
    const struct md_script_ir_t*    ir;     // A compiled script whose names the expression may use. May be NULL.
} md_filter_context_t;

typedef struct md_filter_result_t {
    bool success;                           // The expression picked atoms (possibly none)
    bool is_dynamic;                        // What it picks depends on the coordinates: another state may give others
    md_script_range_marker_t error_range;   // The bytes of the expression the first error is about. Empty when there is
                                            // none, or it is about the expression as a whole.
    char error[256];                        // Why it failed, one line per error. Empty on success.
} md_filter_result_t;

// The atoms that expr picks, all in one set, written to dst. Its previous payload is released with alloc, so it must
// have been produced with alloc (or be empty), and the new one is allocated with it. Empty when it fails.
md_filter_result_t md_filter(struct md_bitset_t* dst, str_t expr, const md_filter_context_t* ctx, struct md_allocator_i* alloc);

// The same, one set for each structure the expression gives rather than all in one: 'residue(within(5, residue(1)))'
// gives a set per residue. The sets are pushed to dst, the array and their payloads allocated with alloc.
md_filter_result_t md_filter_sets(md_array(struct md_bitset_t)* dst, str_t expr, const md_filter_context_t* ctx, struct md_allocator_i* alloc);

// What can be written at the byte offset cursor of expr, while it is typed (see md_script_complete). The variables
// are the names of ctx->ir. Outside of a call only what can be a selection is offered (and the keywords that combine
// them); in the argument of a call everything is, and what has the type of the parameter is marked (fits). The
// system gives the values that selectors take, such as the residue names in resname('|'); ctx->state is not read.
md_script_completions_t md_filter_complete(str_t expr, int cursor, const md_filter_context_t* ctx, struct md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
