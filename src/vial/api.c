#include "ast.h"

#include <core/md_arena_allocator.h>
#include <core/md_common.h>

#include <string.h>

struct md_vial_program_t {
    md_allocator_i*  arena;
    str_t            src;
    vial_token_t*    tokens;
    md_vial_token_t* pub_tokens;
    vial_diags_t     diags;
    vial_tree_t      tree;
};

md_vial_program_t* md_vial_compile(str_t src, const md_vial_compile_desc_t* desc, md_allocator_i* alloc) {
    (void)desc;
    md_allocator_i* arena = md_arena_allocator_create(alloc, KILOBYTES(64));
    if (!arena) return NULL;
    md_vial_program_t* prog = md_alloc(arena, sizeof(md_vial_program_t));
    memset(prog, 0, sizeof(*prog));
    prog->arena = arena;
    prog->src = str_copy(src, arena);

    vial_diags_init(&prog->diags, prog->src, arena);
    prog->tokens = vial_lex(prog->src, &prog->diags, arena);
    prog->tree = vial_parse(prog->tokens, prog->src, &prog->diags, arena);

    const size_t n = md_array_size(prog->tokens);
    md_array_resize(prog->pub_tokens, n, arena);
    for (size_t i = 0; i < n; ++i) {
        prog->pub_tokens[i] = (md_vial_token_t){.kind = prog->tokens[i].kind, .span = prog->tokens[i].span};
    }
    return prog;
}

void md_vial_program_free(md_vial_program_t* prog) {
    if (prog) md_arena_allocator_destroy(prog->arena);
}

size_t md_vial_diagnostics(const md_vial_program_t* prog, const md_vial_diagnostic_t** out) {
    if (!prog) { if (out) *out = NULL; return 0; }
    if (out) *out = prog->diags.items;
    return md_array_size(prog->diags.items);
}

size_t md_vial_tokens(const md_vial_program_t* prog, const md_vial_token_t** out) {
    if (!prog) { if (out) *out = NULL; return 0; }
    if (out) *out = prog->pub_tokens;
    return md_array_size(prog->pub_tokens);
}
