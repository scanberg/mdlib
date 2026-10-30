#pragma once

// The checker: from the syntax tree to a checked program. It grows with the steps of phase 2; this is step 2, names
// and scopes. Every name a script reads is resolved to what bound it: a built-in, a prelude def, a top-level name, a
// def, a parameter, a local, a loop or comprehension variable.

#include "ast.h"
#include "builtins.h"

#include <core/md_hash.h>

typedef enum vial_sym_kind_t {
    VIAL_SYM_BUILTIN = 1,   // builtin points into the built-in table
    VIAL_SYM_PRELUDE,       // a def of the prelude; node is its DEF in the prelude's tree
    VIAL_SYM_GLOBAL,        // a top-level name, bound by assignment
    VIAL_SYM_DEF,           // a def; node is its DEF
    VIAL_SYM_PARAM,         // a parameter of a def
    VIAL_SYM_LOCAL,         // bound inside a def, or inside a top-level loop
    VIAL_SYM_LOOP,          // the target of a top-level for: it belongs to the loop
    VIAL_SYM_COMP,          // the target of a comprehension's for
} vial_sym_kind_t;

typedef struct vial_symbol_t {
    str_t                   name;
    uint8_t                 kind;       // vial_sym_kind_t
    bool                    in_prelude; // bound in the prelude's tree
    uint32_t                scope;
    vial_node_id            node;       // what binds it: a NAME, or a DEF
    uint32_t                stmt;       // VIAL_SYM_GLOBAL and VIAL_SYM_DEF: the top-level statement that binds it
    const vial_builtin_t*   builtin;
} vial_symbol_t;

typedef enum vial_scope_kind_t {
    VIAL_SCOPE_MODULE = 0,
    VIAL_SCOPE_DEF,
    VIAL_SCOPE_COMP,
    VIAL_SCOPE_LOOP,        // a for at top level: its target and what its body binds
} vial_scope_kind_t;

typedef struct vial_scope_t {
    uint8_t         kind;
    uint32_t        parent;     // scope index; the module is 0 and is its own parent
    vial_node_id    node;       // the DEF, LISTCOMP or FOR that opens it
    md_hashmap32_t  names;      // hash of name -> symbol index
} vial_scope_t;

typedef struct vial_names_t {
    vial_symbol_t*  symbols;    // md_array; 0 is unused
    vial_scope_t*   scopes;     // md_array; 0 is the module
    uint32_t*       ref;        // md_array by node id: the symbol a NAME reads or binds, 0 when none
    uint32_t*       scope_of;   // md_array by node id: the scope a DEF, LISTCOMP or FOR opens, 0 when none
    vial_node_id*   parent;     // md_array by node id: the node it is a child of
} vial_names_t;

// Resolves the names of a tree. prelude is the prelude's tree, checked and resolved alike, or NULL; its defs are the
// prelude's names. Diagnostics go to diags, the prelude's to prelude_diags.
typedef struct vial_check_t {
    const vial_tree_t*  tree;
    str_t               src;
    vial_diags_t*       diags;
    vial_names_t        names;

    const vial_tree_t*  prelude;
    str_t               prelude_src;
    vial_diags_t*       prelude_diags;
    vial_names_t        prelude_names;

    const struct md_system_t* sys;  // the system selectors are evaluated against, or NULL

    md_allocator_i*     alloc;
} vial_check_t;

void vial_check_names(vial_check_t* c);

// The symbol a NAME node reads or binds, or NULL
const vial_symbol_t* vial_names_symbol(const vial_names_t* n, vial_node_id id);
