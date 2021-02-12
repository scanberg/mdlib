#include "md_script.h"
#include "md_allocator.h"
#include "md_log.h"
#include "core/common.h"
#include "core/str_util.h"
#include "core/array.inl"
#include "core/string_builder.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

typedef enum token_type {
    TOKEN_TYPE_UNDEF = 0,
    // Reserve the first indices for character literals
    TOKEN_TYPE_IDENT = 128, // idenfitier
    TOKEN_TYPE_LT, // '<='
    TOKEN_TYPE_GT, // '>='
    TOKEN_TYPE_EQ, // '=='
    TOKEN_TYPE_AND,
    TOKEN_TYPE_OR,
    TOKEN_TYPE_NOT,
    TOKEN_TYPE_OF,  // Reserved
    TOKEN_TYPE_IN,  // Reserved
    TOKEN_TYPE_NUM, // Numerical literal or Range
    TOKEN_TYPE_RNG, // Range literal
    TOKEN_TYPE_STR, // String literal, defined between two citation marks "LIKE THIS"
    TOKEN_TYPE_END  // End of tokenstream
} token_type;

typedef struct token_t {
    token_type type;
    // Reference into the original string where the token was parsed from
    uint32_t offset;
    uint32_t length;
} token_t;

typedef enum ast_node_type {
    AST_NODE_VALUE,
    AST_NODE_IDENTIFIER,
    AST_NODE_ARRAY_SUBSCRIPT,
    AST_NODE_FUNC_CALL,
    AST_NODE_ASSIGNMENT,
} ast_node_type;

typedef enum value_type {
    TYPE_INT,
    TYPE_FLOAT,
    TYPE_FLOAT2,
    TYPE_FLOAT3,
    TYPE_FLOAT4,
    TYPE_BOOL,
    TYPE_IRANGE,
    TYPE_FRANGE
} value_type;

typedef struct ast_node_t ast_node_t;
typedef struct ast_node_t {
    ast_node_type type;
    value_type  return_type;
    str_t       str;
    ast_node_t* child[2];
} ast_node_t;

typedef struct tokenizer_t {
    const char* str;
    uint32_t len;
    uint32_t cur;
    uint32_t line;
} tokenizer_t;

bool compare_n(const char* a, const char* b, int n) {
    for (int i = 0; i < n; ++i) {
        if (a[i] != b[i]) return false;
    }
    return true;
}

token_t tokenizer_get_next_from_buffer(tokenizer_t* tokenizer) {
    token_t token = {0};

    if (tokenizer->cur >= tokenizer->len) {
        tokenizer->cur = tokenizer->len;
        token.type = TOKEN_TYPE_END;
        return token;
    }

    bool skip_to_eol = false;
    const char* buf = tokenizer->str;
    const uint32_t len = tokenizer->len;
    for (int i = tokenizer->cur; i != len; ++i) {
        int j = 0;

        if (buf[i] == '#') {
            // Skip to end of line
            while (++i != len) {
                if (buf[i] == '\n') {
                    ++tokenizer->line;
                    break;
                }
            }
        }

        else if (buf[i] == '\n') {
            ++tokenizer->line;
        }

        // Alpha numeric block
        else if (is_alpha(buf[i]) || buf[i] == '_') {
            for (j = i+1; j != len; ++j) {
                if (!is_alpha(buf[j]) && buf[j] != '_' && !is_digit(buf[j])) {
                    break;
                }
            }

            // Keywords
            const int n = j - i;
            if (n == 2) {
                if (compare_n(&buf[i], "or", 2)) {
                    token.type = TOKEN_TYPE_OR;
                }
                else if (compare_n(&buf[i], "in", 2)) {
                    token.type = TOKEN_TYPE_IN;
                }
                else if (compare_n(&buf[i], "of", 2)) {
                    token.type = TOKEN_TYPE_OF;
                }
            }
            else if (n == 3) {
                if (compare_n(&buf[i], "and", 3)) {
                    token.type = TOKEN_TYPE_AND;
                }
                else if (compare_n(&buf[i], "not", 3)) {
                    token.type = TOKEN_TYPE_NOT;
                }
            }

            if (token.type == TOKEN_TYPE_UNDEF) {
                token.type = TOKEN_TYPE_IDENT;
            }
        }

        // Numerical or Range literal
        else if (is_digit(buf[i]) || buf[i] == ':') {
            bool is_range = (buf[i] == ':');
            for (j = i+1; j != len; ++j) {
                if (buf[j] == ':') {
                    is_range = true;
                }
                else if (!is_digit(buf[j]) && buf[j] != '.') {
                    break;
                }
            }
            if (is_range) {
                token.type = TOKEN_TYPE_RNG;
            } else {
                token.type = TOKEN_TYPE_NUM;
            }
        }

        // String literal
        else if (buf[i] == '\'' || buf[i] == '"') {
            char terminator = buf[i]; // Store the termination character to search for
            for (j = i+1; j != len; ++j) {
                if (buf[j] == terminator) {
                    token.type = TOKEN_TYPE_STR;
                    ++j;
                    break;
                }
                else if (buf[j] == '\r' && buf[j+1] == '\n') { // We do not want to leak return carry '\r' into the text since that will malform any output
                    break;
                }
            }
        }

        else if (is_symbol(buf[i])) {
            typedef struct Symbol {
                const char* str;
                token_type type;
            } Symbol;

            static Symbol accepted_symbols_2[] = {
                {"<=", TOKEN_TYPE_LT},
                {">=", TOKEN_TYPE_GT},
                {"==", TOKEN_TYPE_EQ}
            };

            static char accepted_symbols_1[] = {
                '(',
                ')',
                '[',
                ']',
                '{',
                '}',
                '<',
                '>',
                '*',
                '+',
                '-',
                '/',
                '=',
                ';',
                ',',
            };

            // Potential symbol count in buf to match (only care about 1 or 2)
            const uint32_t n = (i + 1 < len && is_symbol(buf[i + 1])) ? 2 : 1;
                
            // Match buf + i against the longest accepted symbol that we can find.
            if (n == 2) {
                for (int k = 0; k < ARRAY_SIZE(accepted_symbols_2); ++k) {
                    if (compare_n(accepted_symbols_2[k].str, &buf[i], 2)) {
                        token.type = accepted_symbols_2[k].type;
                        j = i + 2;
                        break;
                    }
                }
            }
            if (!token.type) {
                for (int k = 0; k < ARRAY_SIZE(accepted_symbols_1); ++k) {
                    if (accepted_symbols_1[k] == buf[i]) {
                        // The character is the type
                        token.type = accepted_symbols_1[k];
                        j = i + 1;
                        break;
                    }
                }
            }
        }

        if (j != 0) {
            token.offset = i;
            token.length = j - i;
            break;
        }
    }

    return token;
}

token_t tokenizer_get_next(tokenizer_t* tokenizer) {
    token_t token = tokenizer_get_next_from_buffer(tokenizer);
    if (token.offset && token.length) {
        // Advance current of tokenizer
        tokenizer->cur = token.offset + token.length;
    }
    return token;
}

token_t tokenizer_peek_next(tokenizer_t* tokenizer) {
    token_t token = tokenizer_get_next_from_buffer(tokenizer);
    return token;
}

tokenizer_t tokenizer_init(const char* str, uint32_t len) {
    tokenizer_t tokenizer = {0};
    tokenizer.str = str;
    tokenizer.len = len;
    tokenizer.line = 1; // We start on line 1 not 0
    return tokenizer;
}

const char* get_token_type_str(token_type type) {
    switch (type) {

        case TOKEN_TYPE_UNDEF: return "undefined";
        case TOKEN_TYPE_IDENT: return "identifier";
        case TOKEN_TYPE_LT: return "less than";
        case TOKEN_TYPE_GT: return "greater than";
        case TOKEN_TYPE_EQ: return "equals";
        case TOKEN_TYPE_AND: return "and";
        case TOKEN_TYPE_OR: return "or";
        case TOKEN_TYPE_NOT: return "not";
        case TOKEN_TYPE_OF: return "of";
        case TOKEN_TYPE_IN: return "in";
        case TOKEN_TYPE_NUM: return "number";
        case TOKEN_TYPE_RNG: return "range";
        case TOKEN_TYPE_STR: return "string";
        case TOKEN_TYPE_END: return "end of token stream";
        default: return "symbol";
    }
}

struct md_script_ir {
    struct md_allocator_i* alloc;
    ast_node_t* nodes;
    md_script_error* errors;
    char* error_str; // One long string which contains all generated error string data.
};

void create_error(md_script_ir* ir, uint32_t line, uint32_t offset, uint32_t length, const char* format, ...) {
    char buffer[512] = {0};
    va_list args;
    va_start(args, format);
    int len = vsnprintf(buffer, ARRAY_SIZE(buffer), format, args);
    va_end(args);

    assert(len > 0);
    assert(len < ARRAY_SIZE(buffer));

    uint64_t ptr_offset = md_array_size(ir->error_str);
    md_array_push_array(ir->error_str, buffer, len + 1, ir->alloc); // len + 1 for the zero terminator character 

    md_script_error error;
    error.line = line;
    error.offset = offset;
    error.length = length;
    // We are snieky and store an offset in the pointer instead of the actual data as error_str will grow and be reallocated
    // This will be patched up and converted into a real pointer once we are done with all error messages.
    error.str_ptr = (const char*)ptr_offset;
    error.str_len = len;
    md_array_push(ir->errors, error, ir->alloc);
}

ast_node_type* create_node(md_script_ir* ir, ast_node_type type) {
    ast_node_t node;
    node.type = type;
    return md_array_push(ir->nodes, node, ir->alloc);
}



// Once all allocations are done, we can convert the all the temporary pointers into actual pointers
void finalize_pointers(md_script_ir* ir) {
    for (uint64_t i = 0; i < md_array_size(ir->nodes); ++i) {
        // @TODO: FIX NODE POINTERS
    }

    for (uint64_t i = 0; i < md_array_size(ir->errors); ++i) {
        uint64_t offset = (uint64_t)ir->errors[i].str_ptr;
        ir->errors[i].str_ptr = ir->error_str + offset;
    }
}

typedef struct parse_context {
    md_script_ir* ir;
    ast_node_t* node; // This represents the current top level node in the tree
    tokenizer_t tokenizer;
} parse_context;

bool parse_expression(parse_context* ctx) {

}

bool parse_logical(parse_context* ctx) {

}

bool parse_expression(parse_context* ctx) {
    token_t token;
    while (true) {
        token = tokenizer_get_next(&ctx->tokenizer);
        printf("'%.*s' -> %s\n", token.length, ctx->tokenizer.str + token.offset, get_token_type_str(token.type));
        
        switch (token.type) {
        case TOKEN_TYPE_END:
            return false;
        case TOKEN_TYPE_UNDEF:
            create_error(ctx->ir, ctx->tokenizer.line, token.offset, token.length, "Undefined token!");
            return false;
        case '(':
        {
            token = tokenizer_get_next(&ctx->tokenizer);
        }
        case TOKEN_TYPE_AND:
        case TOKEN_TYPE_OR:
        case TOKEN_TYPE_NOT:
            return parse_logical(ctx, token);
        case '=':
            return parse_assignment(ctx, token);
        case TOKEN_TYPE_GT:
        case TOKEN_TYPE_LT:
        case TOKEN_TYPE_EQ:
            return parse_comparison(ctx, token);
        case ';':
            return true;
        default:
            create_error(ctx->ir, ctx->tokenizer.line, token.offset, token.length, "Unexpected token!");
            return false;
        }
    }
    return false;
}

struct md_script_ir* md_script_compile(const char* str, uint64_t len, md_allocator_i* alloc) {
    md_script_ir* ir = md_alloc(alloc, sizeof(md_script_ir));
    memset(ir, 0, sizeof(md_script_ir));
    ir->alloc = alloc;
    
    parse_context ctx = {0};
    ctx.tokenizer = tokenizer_init(str, len);
    ctx.ir = ir;

    // Parse expressions until we have exhausted the tokenizer
    while (tokenizer_peek_next(&ctx.tokenizer).type != TOKEN_TYPE_END) {
        parse_expression(&ctx);
    }

    finalize_pointers(ir);

    return ir;
}


bool md_script_compile_success(const md_script_ir* ir) {
    // What do we return here?
    return true;
}

uint32_t md_script_get_num_errors(const md_script_ir* ir) {
    if (ir) {
        return md_array_size(ir->errors);
    }
    return 0;
}
const md_script_error* md_script_get_errors(const md_script_ir* ir) {
    if (ir) {
        return ir->errors;
    }
    return NULL;
}