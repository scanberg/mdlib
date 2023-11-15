#include "md_smiles.h"

#include <core/md_str.h>
#include <core/md_log.h>

#include <stdbool.h>

typedef struct smiles_parse_state_t {
    const char* c;
    const char* end;
    md_smiles_node_t* nodes;
    int64_t node_len;
    int64_t node_cap;
    bool dot;
    bool failed;
} smiles_parse_state_t;

static bool is_valid(smiles_parse_state_t* s) {
    return s->c != s->end;
}

static void terminate(smiles_parse_state_t* s) {
    s->failed = true;
    s->c = s->end;
}

static char peek_char(smiles_parse_state_t* s) {
    return s->c != s->end ? *s->c : '\0';
}

static char consume_char(smiles_parse_state_t* s) {
    return s->c != s->end ? *s->c++ : '\0';
}

static bool parse_integer(int* result, smiles_parse_state_t* s) {
    if (!is_digit(peek_char(s))) return false;
    char c;

    int val = 0;
    while ((c = peek_char(s)) && is_digit(c)) {
        val = val * 10 + (*s->c - '0');
        consume_char(s);
    }

    *result = val;

    return true;
}

static bool parse_symbol(char symbol[3], smiles_parse_state_t* s) {
    if (!is_alpha(peek_char(s))) return false;

    const char symbol1[] = "bcnopsBCFHIKNOPSUVWY";
    const char* symbol2[] = {
        "Ac", "Ag", "Al", "Am", "Ar", "As", "At", "Au",
        "Ba", "Be", "Bh", "Bi", "Bk", "Br",
        "Ca", "Cd", "Ce", "Cf", "Cl", "Cm", "Co", "Cr", "Cs", "Cu",
        "Db", "Ds", "Dy",
        "Er", "Es", "Eu",
        "Fe", "Fl", "Fm", "Fr",
        "Ga", "Gd", "Ge",
        "He", "Hf", "Hg", "Ho", "Hs",
        "In", "Ir",
        "Kr",
        "La", "Li", "Lr", "Lu", "Lv",
        "Mc", "Mg", "Mn", "Mo", "Mt",
        "Na", "Nb", "Nd", "Ne", "Nh", "Ni", "No", "Np",
        "Og", "Os",
        "Pa", "Pb", "Pd", "Pm", "Po", "Pr", "Pt", "Pu",
        "Ra", "Rb", "Re", "Rf", "Rg", "Rh", "Rn", "Ru",
        "Sb", "Sc", "Se", "Sg", "Si", "Sm", "Sn", "Sr",
        "Ta", "Tb", "Tc", "Te", "Th", "Ti", "Tl", "Tm", "Ts",
        "Xe", "Yb", "Zn", "Zr",
    };

    char c[2] = {consume_char(s), 0};
    for (int i = 0; i < ARRAY_SIZE(symbol1); ++i) {
        if (c[0] == symbol1[i]) {
            symbol[0] = c[0];
            return true;
        }
    }

    if (is_alpha(peek_char(s))) {
        for (int i = 0; i < ARRAY_SIZE(symbol2); ++i) {
            if (c[0] == symbol2[i][0] && c[1] == symbol2[i][1]) {
                symbol[0] = c[0];
                symbol[1] = c[1];
                return true;
            }
        }
    }

    return false;
}

static bool parse_organic_symbol(char symbol[3], smiles_parse_state_t* s) {
    if (!is_alpha(peek_char(s))) return false;

    const char  symbol1[] = "bcnopsBCNOFIPS";
    const char* symbol2[] = { "Br", "Cl", "At", "Ts" };

    char c[2] = {consume_char(s), 0};
    for (int i = 0; i < ARRAY_SIZE(symbol1); ++i) {
        if (c[0] == symbol1[i]) {
            symbol[0] = c[0];
            return true;
        }
    }

    if (is_alpha(peek_char(s))) {
        c[1] = consume_char(s);
        for (int i = 0; i < ARRAY_SIZE(symbol2); ++i) {
            if (c[0] == symbol2[i][0] && c[1] == symbol2[i][1]) {
                symbol[0] = c[0];
                symbol[1] = c[1];
                return true;
            }
        }
    }

    return false;
}

static bool parse_bracket_atom(md_smiles_node_t* node, smiles_parse_state_t* s) {
    if (peek_char(s) != '[') {
        return false;
    }
    consume_char(s);

    // Parse isotope
    int iso;
    if (parse_integer(&iso, s)) {
        node->atom.isotope = (uint8_t)iso;
    }

    // Parse symbol
    if (!parse_symbol(node->atom.symbol, s)) {
        return false;
    }

    // Parse chirality
    if (peek_char(s) == '@') {
        node->atom.flags |= MD_SMILES_FLAG_CHIRAL_1;
        consume_char(s);
        if (peek_char(s) == '@') {
            node->atom.flags |= MD_SMILES_FLAG_CHIRAL_2;
            consume_char(s);
        }
    }

    // Parse h_count
    if (peek_char(s) == 'H') {
        consume_char(s);
        if (is_digit(peek_char(s))) {
            int h_count;
            if (parse_integer(&h_count, s)) {
                node->atom.h_count = (uint8_t)h_count;
            }
        } else {
            node->atom.h_count = 1;
        }
    }

    // Parse charge
    if (peek_char(s) == '+' || peek_char(s) == '-') {
        node->atom.charge = consume_char(s) == '+' ? 1 : -1;
        if (is_digit(peek_char(s))) {
            int magnitude;
            parse_integer(&magnitude, s);
            if (magnitude > 15) {
                return false;
            }
            node->atom.charge = (int8_t)(node->atom.charge * magnitude);
        }
    }

    // Parse atom_class
    if (peek_char(s) == ':') {
        consume_char(s);
        if (is_digit(peek_char(s))) {
            int atom_class;
            if (parse_integer(&atom_class, s)) {
                node->atom.atom_class = (uint8_t)atom_class;
            }
        }
    }

    if (consume_char(s) != ']') {
        // Hard error
        MD_LOG_ERROR("Expected termination character ']' in atom entry\n");
        terminate(s);
        return false;
    }

    return true;
}

static bool push_node(smiles_parse_state_t* s, md_smiles_node_t node) {
    if (s->node_len < s->node_cap) {
		s->nodes[s->node_len++] = node;
		return true;
	}
	return false;
}

static bool parse_bond(smiles_parse_state_t* s) {
    char c = peek_char(s);

    const char bond_symbol[] = "-=#$/\\";
    for (int i = 0; i < ARRAY_SIZE(bond_symbol); ++i) {
        if (c == bond_symbol[i]) {
            consume_char(s);
            md_smiles_node_t node = { .type = MD_SMILES_NODE_BOND };
            node.bond.symbol = c;
            return push_node(s, node);
        }
    }

    return false;
}

static bool parse_dot(smiles_parse_state_t* s) {
    char c = peek_char(s);

    if (c == '.') {
        consume_char(s);
        md_smiles_node_t node = { .type = MD_SMILES_NODE_BOND };
        s->dot = true;
        return true;
    }
    return false;
}

static bool parse_rnum(smiles_parse_state_t* s) {
    if (is_digit(peek_char(s))) {
        char c = consume_char(s);
        int num = c - '0';
        md_smiles_node_t node = { .type = MD_SMILES_NODE_RING_CLOSURE };
        node.ring.index = num;
        return push_node(s, node);
    }

    if (peek_char(s) == '%') {
        consume_char(s);
        char c[2] = { consume_char(s), consume_char(s) };
        if (!is_digit(c[0]) || !is_digit(c[1])) {
            MD_LOG_ERROR("Expected ring number after '%%'\n");
            terminate(s);
            return false;
        }
        int num = (c[0] - '0') * 10 + (c[1] - '0');
        md_smiles_node_t node = { .type = MD_SMILES_NODE_RING_CLOSURE };
        node.ring.index = num;
        return push_node(s, node);
    }

    return false;
}

static bool parse_line(smiles_parse_state_t*);

static bool parse_branch(smiles_parse_state_t* s) {
    if (peek_char(s) != '(') {
        return false;
    }
    consume_char(s);

    md_smiles_node_t node = { .type = MD_SMILES_NODE_BRANCH_OPEN };
    if (!push_node(s, node)) {
        return false;
    }

    int num = 0;
    while (peek_char(s) != ')') {
        parse_bond(s) || parse_dot(s);
        if (!parse_line(s)) {
            MD_LOG_ERROR("Failed to parse line in branch\n");
            terminate(s);
            return false;
        }
        num += 1;
    }

    if (num < 1) {
        MD_LOG_ERROR("Expected at least one entry in branch\n");
        terminate(s);
        return false;
    }

    if (is_valid(s) && consume_char(s) != ')') {
        MD_LOG_ERROR("Expected termination character ')' in branch entry\n");
        terminate(s);
        return false;
    }

    node.type = MD_SMILES_NODE_BRANCH_CLOSE;
    return push_node(s, node);
}

static bool parse_joker_symbol(char symbol[3], smiles_parse_state_t* s) {
    if (peek_char(s) != '*') {
        return false;
    }
    consume_char(s);
    symbol[0] = '*';
    return true;
}

static bool parse_atom(smiles_parse_state_t* s) {
    md_smiles_node_t node = { .type = MD_SMILES_NODE_ATOM };
    if (parse_organic_symbol(node.atom.symbol, s) || parse_bracket_atom(&node, s) || parse_joker_symbol(node.atom.symbol, s)) {
        s->dot = false;
        return push_node(s, node);
    }
    return false;
}

static bool parse_chain(smiles_parse_state_t* s) {
    int num = 0;
    while (is_valid(s)) {
        if (parse_dot(s) && parse_atom(s)) {
            num += 1;
            continue;
        }
        parse_bond(s);
        if (parse_atom(s) || parse_rnum(s)) {
            num += 1;
            continue;
        }
        break;
    }
    if (num < 1) {
        return false;
    }
    return true;
}

static bool parse_line(smiles_parse_state_t* s) {
    if (!parse_atom(s)) {
        return false;
    }
    while (parse_chain(s) || parse_branch(s));
    return true;
}

int64_t md_smiles_parse(md_smiles_node_t* out_nodes, int64_t in_cap, const char* in_str, int64_t in_len) {
    if (!in_str || in_len <= 0) {
        MD_LOG_ERROR("Smiles: Invalid input string");
        return 0;
    }

    str_t str = str_trim((str_t){in_str, in_len});

    smiles_parse_state_t state = {
        .c = str_beg(str),
        .end = str_end(str),
        .nodes = out_nodes,
        .node_len = 0,
        .node_cap = in_cap
    };

    while (parse_chain(&state) || parse_branch(&state));

    return state.node_len;
}