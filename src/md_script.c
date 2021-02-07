#include "md_script.h"
#include "md_allocator.h"

typedef enum token_type {
    TOKEN_TYPE_UNDEFINED = 0,
    // Reserve the first indices for character literals
    TOKEN_TYPE_IDENT = 128, // idenfitier
    TOKEN_TYPE_FUNC,    // built in function
    TOKEN_TYPE_NOT,
    TOKEN_TYPE_AND,
    TOKEN_TYPE_OR,
    TOKEN_TYPE_OF,  // reserved
    TOKEN_TYPE_IN,  // reserved
    TOKEN_TYPE_INT, // Numerical literal
    TOKEN_TYPE_FLT,
    TOKEN_TYPE_IRNG,
    TOKEN_TYPE_FRNG,
    TOKEN_TYPE_STR, // String literal
    TOKEN_TYPE_END  // Marks end of tokenstream
} token_type;

typedef struct token {
    const char* beg;
    const char* end;
} token;

typedef enum ast_node_type {
    AST_NODE_INT,
    AST_NODE_FLOAT,
} ast_node_type;

struct ast_node {
    ast_node_type type;
};


struct md_script_ir* md_script_compile(const char* str, uint64_t len) {

}