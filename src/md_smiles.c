#include "md_smiles.h"

#if defined(_MSC_VER)
#pragma warning(disable: 4295)
#endif

#include <core/md_log.h>

#include <stdbool.h>

typedef struct state_t {
    const char* c;
    const char* end;
    md_smiles_node_t* nodes;
    int64_t node_len;
    int64_t node_cap;
    bool dot;
    bool abort;
} state_t;

typedef struct symbol_t {
    char symbol[2];
    uint8_t element;
    uint8_t flags;
} symbol_t;


// The matching table is segmented to first accomodate shortcuts (B, C, N, O, P, S, F, Cl, Br, and I)
// Then the rest of the elements follow
static const symbol_t table[] = {
    { "b", 5, MD_SMILES_FLAG_AROMATIC }, { "c", 6, MD_SMILES_FLAG_AROMATIC }, { "n", 7, MD_SMILES_FLAG_AROMATIC }, { "o", 8, MD_SMILES_FLAG_AROMATIC }, { "p", 15, MD_SMILES_FLAG_AROMATIC }, { "s", 16, MD_SMILES_FLAG_AROMATIC },
    { "B", 5 }, { "C", 6 }, { "N", 7 }, { "O", 8 }, { "P", 15 }, { "S", 16 }, { "F", 9 }, { "I", 53 },
    { "Cl", 17 }, { "Br", 35 },
    { "H", 1 }, { "K", 19 }, { "P", 15 }, { "U", 92 }, { "V", 23 }, { "W", 74 }, { "Y", 39 },
	{ "Ac", 89 }, { "Ag", 47 }, { "Al", 13 }, { "Am", 95 }, { "Ar", 18 }, { "As", 33 }, { "At", 85 }, { "Au", 79 },
	{ "Ba", 56 }, { "Be", 4 }, { "Bh", 107 }, { "Bi", 83 }, { "Bk", 97 },
    { "Ca", 20 }, { "Cd", 48 }, { "Ce", 58 }, { "Cf", 98 }, {"Cm", 96 }, { "Co", 27 }, { "Cr", 24 }, { "Cs", 55 }, { "Cu", 29 },
	{ "Db", 105 }, { "Ds", 110 }, { "Dy", 66 },
	{ "Er", 68 }, { "Es", 99 }, { "Eu", 63 },
	{ "Fe", 26 }, { "Fl", 114 }, { "Fm", 100 }, { "Fr", 87 },
	{ "Ga", 31 }, { "Gd", 64 }, { "Ge", 32 },
	{ "He", 2 }, { "Hf", 72 }, { "Hg", 80 }, { "Ho", 67 }, { "Hs", 108 },
	{ "In", 49 }, { "Ir", 77 },
	{ "Kr", 36 },
	{ "La", 57 }, { "Li", 3 }, { "Lr", 103 }, { "Lu", 71 }, { "Lv", 116 },
	{ "Mc", 115 }, { "Mg", 12 }, { "Mn", 25 }, { "Mo", 42 }, { "Mt", 109 },
	{ "Na", 11 }, { "Nb", 41 }, { "Nd", 60 }, { "Ne", 10 }, { "Nh", 113 }, { "Ni", 28 }, { "No", 102 }, { "Np", 93 },
	{ "Og", 118 }, { "Os", 76 },
	{ "Pa", 91 }, { "Pb", 82 }, { "Pd", 46 }, { "Pm", 61 }, { "Po", 84 }, { "Pr", 59 }, { "Pt", 78 }, { "Pu", 94 },
	{ "Ra", 88 }, { "Rb", 37 }, { "Re", 75 }, { "Rf", 104 }, { "Rg", 111 }, { "Rh", 45 }, { "Rn", 86 }, { "Ru", 44 },
	{ "Sb", 51 }, { "Sc", 21 }, { "Se", 34 }, { "Sg", 106 }, { "Si", 14 }, { "Sm", 62 }, { "Sn", 50 }, { "Sr", 38 },
	{ "Ta", 73 }, { "Tb", 65 }, { "Tc", 43 }, { "Te", 52 }, { "Th", 90 }, { "Ti", 22 }, { "Tl", 81 }, { "Tm", 69 }, { "Ts", 117 },
	{ "Xe", 54 }, { "Yb", 70 }, { "Zn", 30 }, { "Zr", 40 },
};

static const symbol_t* shortcut_end = table + 16;
static const symbol_t* single_char_end = table + 23;

static inline bool is_digit(int c) { return '0' <= c && c <= '9'; }
static inline bool is_alpha(int c) { return ('a' <= c && c <= 'z') || ('A' <= c && c <= 'Z'); }
static inline bool is_lower(int c) { return 'a' <= c && c <= 'z'; }
static inline bool is_upper(int c) { return 'A' <= c && c <= 'Z'; }
static inline bool is_whitespace(int c) { return c == ' ' || c == '\n' || c == '\r' || c == '\t'; }

static bool is_valid(state_t* s) {
    return s->c != s->end;
}

static void terminate(state_t* s) {
    s->abort = true;
    s->c = s->end;
}

static char peek_char(state_t* s) {
    return s->c != s->end ? *s->c : 0;
}

static char consume_char(state_t* s) {
    return s->c != s->end ? *s->c++ : 0;
}

static bool parse_integer(int* result, state_t* s) {
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

static const symbol_t* match_symbol(bool shortcut, state_t* s) {
    if (!is_alpha(peek_char(s))) {
        return false;
    }
    char c0 = consume_char(s);
    char c1 = peek_char(s);
    const char str[2] = {c0, (is_upper(c0) && is_lower(c1)) ? c1 : '\0'};

    const symbol_t* sym = table;
    const symbol_t* end = shortcut ? shortcut_end : (str[1] ? table + ARRAY_SIZE(table) : single_char_end);
    const symbol_t* match1 = NULL;

	while (sym != end) {
        if (sym->symbol[1] == '\0' && sym->symbol[0] == str[0]) {
        	match1 = sym;
            if (str[1] == '\0') {
                return match1;
            }
        }
		if (sym->symbol[0] == str[0] && sym->symbol[1] == str[1]) {
			consume_char(s);
            return sym;
		}
		++sym;
	}
	return match1;
}

static bool parse_shortcut(md_smiles_node_t* node, state_t* s) {
    const symbol_t* sym = match_symbol(true, s);
    if (!sym) {
		return false;
	}
    node->atom.element = sym->element;
	node->atom.flags  |= sym->flags;

	return true;
}

static bool parse_bracket(md_smiles_node_t* node, state_t* s) {
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
    const symbol_t* sym = match_symbol(false, s);
    if (!sym) {
        return false;
    }

    node->atom.element = sym->element;
    node->atom.flags  |= sym->flags;

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
                node->atom.hydrogen_count = (uint8_t)h_count;
            }
        } else {
            node->atom.hydrogen_count = 1;
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

static bool push_node(state_t* s, md_smiles_node_t node) {
    if (s->node_len < s->node_cap) {
		s->nodes[s->node_len++] = node;
		return true;
	}
	return false;
}

static bool parse_bond(state_t* s) {
    char c = peek_char(s);

    const char bond_symbol[] = "-=:#$/\\";
    for (int i = 0; i < (int)ARRAY_SIZE(bond_symbol); ++i) {
        if (c == bond_symbol[i]) {
            consume_char(s);
            md_smiles_node_t node = { .type = MD_SMILES_NODE_BOND };
            node.bond.symbol = c;
            node.bond.flags |= (c == ':')  ? MD_SMILES_FLAG_AROMATIC : 0;
            node.bond.flags |= (c == '/')  ? MD_SMILES_FLAG_UP : 0;
            node.bond.flags |= (c == '\\') ? MD_SMILES_FLAG_DOWN : 0;
            return push_node(s, node);
        }
    }

    return false;
}

static bool parse_dot(state_t* s) {
    char c = peek_char(s);

    if (c == '.') {
        consume_char(s);
        //md_smiles_node_t node = { .type = MD_SMILES_NODE_BOND };
        // Dot is not supported yet, so now we just ignore it and set the state flag
        s->dot = true;
        return true;
    }
    return false;
}

static bool parse_rnum(state_t* s) {
    if (is_digit(peek_char(s))) {
        char c = consume_char(s);
        int num = c - '0';
        md_smiles_node_t node = { .type = MD_SMILES_NODE_BRIDGE };
        node.bridge.index = num;
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
        md_smiles_node_t node = { .type = MD_SMILES_NODE_BRIDGE };
        node.bridge.index = num;
        return push_node(s, node);
    }

    return false;
}

static bool parse_line(state_t*);

static bool parse_branch(state_t* s) {
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
        parse_bond(s);
        parse_dot(s);
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

static bool parse_star(md_smiles_node_t* node, state_t* s) {
    if (peek_char(s) != '*') {
        return false;
    }
    consume_char(s);
    node->atom.element = 0;
    return true;
}

static bool parse_atom(state_t* s) {
    md_smiles_node_t node = { .type = MD_SMILES_NODE_ATOM };
    if (parse_star(&node, s) || parse_shortcut(&node, s) || parse_bracket(&node, s)) {
        s->dot = false;
        return push_node(s, node);
    }
    return false;
}

static bool parse_chain(state_t* s) {
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

static bool parse_line(state_t* s) {
    if (!parse_atom(s)) {
        return false;
    }
    while (parse_chain(s) || parse_branch(s));
    return true;
}

int64_t md_smiles_parse(md_smiles_node_t* out_nodes, size_t in_cap, const char* in_str, size_t in_len) {
    if (!in_str || in_len <= 0) {
        MD_LOG_ERROR("Smiles: Invalid input string");
        return 0;
    }

    const char* beg = in_str;
    const char* end = in_str + in_len;

    // Do some trimming
    while (beg < end && is_whitespace(*beg)) {
        beg++;
    }

    while (beg < end && (is_whitespace(end[-1]) || end[-1] == '\0')) {
		end--;
	}

    state_t state = {
        .c = beg,
        .end = end,
        .nodes = out_nodes,
        .node_len = 0,
        .node_cap = in_cap
    };

    while (parse_chain(&state) || parse_branch(&state));

    return state.node_len;
}
