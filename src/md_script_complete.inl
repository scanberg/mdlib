// ### COMPLETION ###
//
// md_script_complete tells an editor what can be written at a position of a script. A script is incomplete while it
// is being typed, so nothing is parsed or compiled here: the tokenizer runs up to the cursor, and the brackets of the
// statement the cursor is in are followed. That is enough to know
//   - what is being typed: an identifier, the contents of a string, or nothing yet,
//   - the procedure call the cursor is in, and which parameter the argument is for, by position or by name,
//   - whether the cursor is left of the '=' of an assignment, where the name of a variable is written.
// What string parameters take is declared in the 'param_domains' table (md_script_functions.inl) and read from the
// system here, with the same accessors as the procedures use to match their arguments.

#include <stdlib.h>

#define COMPLETION_MAX_DEPTH 32

typedef struct completion_frame_t {
    int               open;         // '(', '[' or '{'
    str_t             callee;       // For the '(' of a procedure call: the name of the procedure
    const proc_sig_t* sig;          // Its signature, if it has one
    int               arg;          // Position of the current argument (or array element)
    str_t             arg_name;     // The name the current argument is given by, empty if it is given by position
    int               arg_tokens;   // Number of tokens of the current argument so far
    uint64_t          given;        // The parameters that earlier arguments were given to, one bit per position
} completion_frame_t;

typedef struct completion_context_t {
    completion_frame_t frame[COMPLETION_MAX_DEPTH];
    int  depth;                     // Can exceed COMPLETION_MAX_DEPTH, frames beyond it are not recorded
    bool lhs;                       // Left of the '=' of an assignment
} completion_context_t;

static completion_frame_t* completion_top(completion_context_t* ctx) {
    return (0 < ctx->depth && ctx->depth <= COMPLETION_MAX_DEPTH) ? &ctx->frame[ctx->depth - 1] : NULL;
}

// The parameter that the current argument of a call is for, -1 if unknown
static int completion_param(const completion_frame_t* frame) {
    if (!str_empty(frame->arg_name)) {
        return frame->sig ? find_param_index(frame->sig, frame->arg_name) : -1;
    }
    return frame->arg;
}

// The procedure call whose argument the cursor is in, looking through array literals: resname({"ALA", "GLY"})
static const completion_frame_t* completion_call(const completion_context_t* ctx) {
    if (ctx->depth > COMPLETION_MAX_DEPTH) return NULL;
    for (int d = ctx->depth - 1; d >= 0; --d) {
        const completion_frame_t* frame = &ctx->frame[d];
        if (frame->open == '(') return str_empty(frame->callee) ? NULL : frame;
        if (frame->open != '{') return NULL;
    }
    return NULL;
}

static value_domain_t find_param_domain(str_t proc, int param) {
    if (param < 0) return DOMAIN_NONE;
    for (size_t i = 0; i < ARRAY_SIZE(param_domains); ++i) {
        if (param_domains[i].param == (uint32_t)param && str_eq(param_domains[i].proc, proc)) {
            return param_domains[i].domain;
        }
    }
    return DOMAIN_NONE;
}

static bool is_keyword_token(int type) {
    switch (type) {
    case TOKEN_AND: case TOKEN_OR: case TOKEN_XOR: case TOKEN_NOT: case TOKEN_IN: case TOKEN_OF: case TOKEN_OUT:
        return true;
    default:
        return false;
    }
}

// The tokens of str
static md_array(token_t) completion_tokenize(str_t str, md_allocator_i* alloc) {
    md_array(token_t) tokens = 0;
    tokenizer_t tokenizer = tokenizer_init(str);
    for (;;) {
        const token_t token = tokenizer_consume_next(&tokenizer);
        if (token.type == TOKEN_END || token.str.len == 0) break;
        md_array_push(tokens, token, alloc);
    }
    return tokens;
}

// Follows the tokens of the statement that ends at the last of them
static void completion_scan(completion_context_t* ctx, const token_t* tokens, size_t count) {
    MEMSET(ctx, 0, sizeof(*ctx));
    ctx->lhs = true;

    size_t beg = count;
    while (beg > 0 && tokens[beg - 1].type != ';') --beg;

    for (size_t i = beg; i < count; ++i) {
        const int type = tokens[i].type;
        completion_frame_t* top = completion_top(ctx);

        // Only names, in braces and separated by commas, can be left of the '=' of an assignment
        if (!(type == TOKEN_IDENT || type == '{' || type == '}' || type == ',')) {
            ctx->lhs = false;
        }

        if (type == '(' || type == '[' || type == '{') {
            if (top) top->arg_tokens += 1;
            ctx->depth += 1;
            completion_frame_t* frame = completion_top(ctx);
            if (frame) {
                MEMSET(frame, 0, sizeof(*frame));
                frame->open = type;
                if (type == '(' && i > beg && tokens[i - 1].type == TOKEN_IDENT) {
                    frame->callee = tokens[i - 1].str;
                    frame->sig = find_proc_signature(frame->callee);
                }
            }
        } else if (type == ')' || type == ']' || type == '}') {
            if (ctx->depth > 0) ctx->depth -= 1;
        } else if (type == ',') {
            if (top) {
                const int param = completion_param(top);
                if (0 <= param && param < 64) top->given |= 1ULL << param;
                top->arg += 1;
                top->arg_name = (str_t){0};
                top->arg_tokens = 0;
            }
        } else if (type == '=' && top && top->open == '(' && top->arg_tokens == 1 && tokens[i - 1].type == TOKEN_IDENT) {
            // 'name=' starts an argument given by name
            top->arg_name = tokens[i - 1].str;
            top->arg_tokens = 0;
        } else if (top) {
            top->arg_tokens += 1;
        }
    }
}

// The names the script assigns to, in the order they first appear: the identifiers left of the '=' in 'x = ...' and
// '{a, b} = ...'. Not those of the statement the cursor is in, which are being written: a variable is not available
// in its own definition.
static md_array(str_t) completion_variables(const token_t* tokens, size_t count, int cursor, md_allocator_i* alloc) {
    md_array(str_t) names = 0;
    size_t beg = 0;
    while (beg < count) {
        size_t end = beg;
        while (end < count && tokens[end].type != ';') ++end;

        size_t eq = beg;
        bool names_only = true;
        while (eq < end && tokens[eq].type != '=') {
            const int type = tokens[eq].type;
            if (!(type == TOKEN_IDENT || type == '{' || type == '}' || type == ',')) {
                names_only = false;
                break;
            }
            ++eq;
        }

        const int stmt_beg = tokens[beg].beg;
        const int stmt_end = end < count ? tokens[end].end : INT_MAX;
        const bool at_cursor = stmt_beg <= cursor && cursor < stmt_end;

        if (names_only && eq < end && !at_cursor) {
            for (size_t i = beg; i < eq; ++i) {
                if (tokens[i].type != TOKEN_IDENT) continue;
                bool found = false;
                for (size_t j = 0; j < md_array_size(names) && !found; ++j) found = str_eq(names[j], tokens[i].str);
                if (!found) md_array_push(names, tokens[i].str, alloc);
            }
        }
        beg = end + 1;
    }
    return names;
}

static void completion_add_value(md_array(str_t)* values, md_hashmap32_t* seen, str_t value, md_allocator_i* alloc) {
    if (str_empty(value)) return;
    const uint64_t key = md_hash64_str(value, 0) & 0x7FFFFFFFFFFFFFFFULL;  // Keep clear of the reserved keys
    const uint32_t* idx = md_hashmap_get(seen, key);
    if (idx) {
        if (str_eq((*values)[*idx], value)) return;
        // A different value with the same hash, keep it (but it is not found by its hash)
    } else {
        md_hashmap_add(seen, key, (uint32_t)md_array_size(*values));
    }
    md_array_push(*values, value, alloc);
}

static int completion_cmp_str(const void* a, const void* b) {
    return str_cmp_lex(*(const str_t*)a, *(const str_t*)b);
}

// The paths of the attributes that attr() can read, in the form it resolves them (see static_check_attribute):
// relative to their run where that names exactly one attribute, in full otherwise. In full also when what has been
// typed starts with 'run', which is how a full path starts.
static void completion_attr_paths(md_array(str_t)* values, md_hashmap32_t* seen, const md_attributes_t* attributes, str_t typed, md_allocator_i* alloc) {
    const size_t num_ids = md_attributes_query_flags(NULL, 0, attributes, STR_LIT("run"), MD_ATTRIBUTE_FLAG_TEMPORAL, MD_ATTRIBUTE_FLAG_TEMPORAL);
    if (num_ids == 0) return;
    md_attribute_id_t* ids = md_alloc(alloc, num_ids * sizeof(md_attribute_id_t));
    md_attributes_query_flags(ids, num_ids, attributes, STR_LIT("run"), MD_ATTRIBUTE_FLAG_TEMPORAL, MD_ATTRIBUTE_FLAG_TEMPORAL);

    str_t runs[64];
    const size_t num_runs = MIN(md_attributes_query_children(runs, ARRAY_SIZE(runs), attributes, STR_LIT("run")), ARRAY_SIZE(runs));
    const bool full_form = str_begins_with(typed, STR_LIT("run"));

    for (size_t i = 0; i < num_ids; ++i) {
        const md_attribute_t* attr = md_attributes_get(attributes, ids[i]);
        attr_read_info_t info;
        if (!attr || attr_readable(attributes, attr, &info) != ATTR_READABLE) continue;

        const str_t relative = str_substr(attr->path, info.run.len + 1, SIZE_MAX);
        bool unique = !full_form && !md_attributes_find(attributes, relative);  // A full path is looked for first
        for (size_t r = 0; r < num_runs && unique; ++r) {
            char buf[512];
            const int len = snprintf(buf, sizeof(buf), "run/"STR_FMT"/"STR_FMT, STR_ARG(runs[r]), STR_ARG(relative));
            if (len <= 0 || (size_t)len >= sizeof(buf)) continue;
            const md_attribute_t* other = md_attributes_find(attributes, (str_t){buf, (size_t)len});
            if (other && other != attr) unique = false;
        }
        completion_add_value(values, seen, unique ? relative : attr->path, alloc);
    }
}

// The distinct values of a domain, sorted. The strings point into the system or static memory.
// typed is what has been typed of the value so far, which only decides the form of attribute paths.
static md_array(str_t) completion_domain_values(value_domain_t domain, const md_system_t* sys, str_t typed, md_allocator_i* alloc) {
    md_array(str_t) values = 0;
    md_hashmap32_t seen = {.allocator = alloc};

    switch (domain) {
    case DOMAIN_ATOM_NAME:
        // Atom names are the names of their types (md_atom_name), as matched by _name
        if (sys && sys->atom.type.name) {
            for (size_t i = 0; i < sys->atom.type.count; ++i) {
                completion_add_value(&values, &seen, LBL_TO_STR(sys->atom.type.name[i]), alloc);
            }
        }
        break;
    case DOMAIN_ELEMENT:
        if (sys) {
            bool present[256] = {0};
            for (size_t i = 0; i < sys->atom.count; ++i) {
                present[md_atom_atomic_number(&sys->atom, i)] = true;
            }
            for (int z = 1; z < 256; ++z) {
                if (present[z]) completion_add_value(&values, &seen, md_util_element_symbol((md_element_t)z), alloc);
            }
        }
        break;
    case DOMAIN_COMP_NAME:
        if (sys && sys->component.name) {
            str_t prev = {0};
            for (size_t i = 0; i < sys->component.count; ++i) {
                const str_t name = LBL_TO_STR(sys->component.name[i]);
                if (str_eq(name, prev)) continue;   // Runs of the same component (water) are common
                completion_add_value(&values, &seen, name, alloc);
                prev = name;
            }
        }
        break;
    case DOMAIN_CHAIN_ID:
    case DOMAIN_INST_ID:
    case DOMAIN_INST_AUTH_ID:
        if (sys) {
            for (size_t i = 0; i < md_system_instance_count(sys); ++i) {
                if (domain == DOMAIN_CHAIN_ID && !is_instance_chain(sys, i)) continue;
                const str_t id = domain == DOMAIN_INST_AUTH_ID ? md_system_instance_auth_id(sys, i) : md_system_instance_id(sys, i);
                completion_add_value(&values, &seen, id, alloc);
            }
        }
        break;
    case DOMAIN_COUNT_UNIT:
        for (size_t i = COUNT_TYPE_UNKNOWN + 1; i < COUNT_TYPE_COUNT; ++i) {
            completion_add_value(&values, &seen, count_type_str[i], alloc);
        }
        break;
    case DOMAIN_ATTR_PATH:
        if (sys) {
            completion_attr_paths(&values, &seen, &sys->attributes, typed, alloc);
        }
        break;
    default:
        break;
    }

    md_hashmap_free(&seen);
    if (md_array_size(values) > 1) {
        qsort(values, md_array_size(values), sizeof(str_t), completion_cmp_str);
    }
    return values;
}

static void completion_push(md_array(md_script_completion_t)* items, md_script_completion_kind_t kind, str_t label, str_t pre, str_t post, md_allocator_i* alloc) {
    md_script_completion_t item = {.kind = kind};
    item.label = str_copy(label, alloc);
    if (pre.len == 0 && post.len == 0) {
        item.text = item.label;
    } else {
        const size_t len = pre.len + label.len + post.len;
        char* buf = md_alloc(alloc, len + 1);
        if (pre.len)   MEMCPY(buf, pre.ptr, pre.len);
        if (label.len) MEMCPY(buf + pre.len, label.ptr, label.len);
        if (post.len)  MEMCPY(buf + pre.len + label.len, post.ptr, post.len);
        buf[len] = '\0';
        item.text = (str_t){buf, len};
    }
    md_array_push(*items, item, alloc);
}

static void completion_push_values(md_array(md_script_completion_t)* items, value_domain_t domain, const md_system_t* sys, str_t typed, str_t pre, str_t post, md_allocator_i* alloc, md_allocator_i* temp_alloc) {
    md_array(str_t) values = completion_domain_values(domain, sys, typed, temp_alloc);
    for (size_t i = 0; i < md_array_size(values); ++i) {
        completion_push(items, MD_SCRIPT_COMPLETION_VALUE, values[i], pre, post, alloc);
    }
}

md_script_completions_t md_script_complete(str_t src, int cursor, const md_system_t* sys, md_allocator_i* alloc) {
    ASSERT(alloc);
    md_script_completions_t result = {0};
    if (!src.ptr) src.len = 0;
    cursor = CLAMP(cursor, 0, (int)src.len);
    result.range.beg = cursor;
    result.range.end = cursor;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    md_array(md_script_completion_t) items = 0;

    // What is being typed, from the last token before the cursor
    const str_t head = str_substr(src, 0, (size_t)cursor);
    const token_t* tokens = completion_tokenize(head, temp_alloc);
    size_t num_tokens = md_array_size(tokens);
    const token_t* last = num_tokens ? &tokens[num_tokens - 1] : NULL;

    enum { PARTIAL_NONE, PARTIAL_IDENT, PARTIAL_STRING, PARTIAL_NOTHING_TO_OFFER } partial = PARTIAL_NONE;
    char quote = 0;
    int beg = cursor;

    if (last && last->end == cursor) {
        if (last->type == TOKEN_IDENT || is_keyword_token(last->type)) {
            partial = PARTIAL_IDENT;
            beg = last->beg;
            num_tokens -= 1;
        } else if (last->type == TOKEN_UNDEF && (last->str.ptr[0] == '"' || last->str.ptr[0] == '\'')) {
            // A string that is not terminated (yet)
            partial = PARTIAL_STRING;
            quote = last->str.ptr[0];
            beg = last->beg + 1;
            num_tokens -= 1;
        } else if (last->type == TOKEN_INT || last->type == TOKEN_FLOAT || last->type == TOKEN_STRING) {
            partial = PARTIAL_NOTHING_TO_OFFER;  // In a number, or right after a string
        }
    } else {
        // The tokenizer skips comments, so one that the cursor is in starts after the last token
        for (int i = last ? last->end : 0; i < cursor; ++i) {
            if (src.ptr[i] == '#') {
                partial = PARTIAL_NOTHING_TO_OFFER;
                break;
            }
        }
    }

    int end = cursor;
    if (partial == PARTIAL_IDENT) {
        while (end < (int)src.len && (is_alpha(src.ptr[end]) || is_digit(src.ptr[end]) || src.ptr[end] == '_')) ++end;
    } else if (partial == PARTIAL_STRING) {
        while (end < (int)src.len && src.ptr[end] != quote && src.ptr[end] != ';' && src.ptr[end] != '\n' && src.ptr[end] != '\r') ++end;
    }
    result.range.beg = beg;
    result.range.end = end;
    result.prefix = str_copy(str_substr(src, (size_t)beg, (size_t)(cursor - beg)), alloc);

    if (partial != PARTIAL_NOTHING_TO_OFFER) {
        completion_context_t ctx;
        completion_scan(&ctx, tokens, num_tokens);
        completion_frame_t* top = completion_top(&ctx);
        const completion_frame_t* call = completion_call(&ctx);
        const value_domain_t domain = call ? find_param_domain(call->callee, completion_param(call)) : DOMAIN_NONE;
        const str_t quote_str = {partial == PARTIAL_STRING ? &quote : "\"", 1};

        if (partial == PARTIAL_STRING) {
            // Only the values the argument takes, and the closing quote if the string does not have one yet
            if (domain != DOMAIN_NONE) {
                const bool closed = end < (int)src.len && src.ptr[end] == quote;
                completion_push_values(&items, domain, sys, result.prefix, (str_t){0}, closed ? (str_t){0} : quote_str, alloc, temp_alloc);
            }
        } else {
            const token_t* all_tokens = completion_tokenize(src, temp_alloc);
            md_array(str_t) variables = completion_variables(all_tokens, md_array_size(all_tokens), cursor, temp_alloc);

            if (ctx.lhs) {
                // Where a variable is named, only the names that are taken
                for (size_t i = 0; i < md_array_size(variables); ++i) {
                    completion_push(&items, MD_SCRIPT_COMPLETION_VARIABLE, variables[i], (str_t){0}, (str_t){0}, alloc);
                }
            } else {
                const bool arg_start = top && top->arg_tokens == 0;
                if (domain != DOMAIN_NONE && arg_start) {
                    completion_push_values(&items, domain, sys, result.prefix, quote_str, quote_str, alloc, temp_alloc);
                }
                if (arg_start && top == call && str_empty(top->arg_name) && top->sig) {
                    for (size_t i = 0; i < top->sig->num_params; ++i) {
                        if (i < 64 && (top->given & (1ULL << i))) continue;
                        completion_push(&items, MD_SCRIPT_COMPLETION_PARAMETER, top->sig->param[i].name, (str_t){0}, STR_LIT("="), alloc);
                    }
                }
                for (size_t i = 0; i < md_array_size(variables); ++i) {
                    completion_push(&items, MD_SCRIPT_COMPLETION_VARIABLE, variables[i], (str_t){0}, (str_t){0}, alloc);
                }
                for (size_t i = 0; i < ARRAY_SIZE(procedures); ++i) {
                    if (!procedure_name_before(i, procedures[i].name)) {
                        completion_push(&items, MD_SCRIPT_COMPLETION_PROCEDURE, procedures[i].name, (str_t){0}, (str_t){0}, alloc);
                    }
                }
                for (size_t i = 0; i < ARRAY_SIZE(script_intrinsics); ++i) {
                    if (!procedure_name_before(ARRAY_SIZE(procedures), script_intrinsics[i])) {
                        completion_push(&items, MD_SCRIPT_COMPLETION_PROCEDURE, script_intrinsics[i], (str_t){0}, (str_t){0}, alloc);
                    }
                }
                for (size_t i = 0; i < ARRAY_SIZE(constants); ++i) {
                    completion_push(&items, MD_SCRIPT_COMPLETION_CONSTANT, constants[i].name, (str_t){0}, (str_t){0}, alloc);
                }
                for (size_t i = 0; i < ARRAY_SIZE(script_keywords); ++i) {
                    completion_push(&items, MD_SCRIPT_COMPLETION_KEYWORD, script_keywords[i], (str_t){0}, (str_t){0}, alloc);
                }
            }
        }
    }

    md_temp_end(temp);
    result.items = items;
    result.count = md_array_size(items);
    return result;
}
