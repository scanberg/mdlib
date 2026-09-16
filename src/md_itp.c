#include <md_itp.h>

#include <md_system.h>
#include <md_util.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_parse.h>
#include <core/md_str.h>
#include <core/md_str_builder.h>
#include <core/md_unit.h>

#include <stdlib.h>
#include <math.h>

#define ITP_MAX_INCLUDE_DEPTH 16
#define ITP_MAX_COND_DEPTH    64
#define ITP_MAX_TOKENS        32

typedef enum {
    SECTION_NONE = 0,
    SECTION_IGNORE,
    SECTION_ATOMTYPES,
    SECTION_MOLECULETYPE,
    SECTION_ATOMS,
    SECTION_BONDS,          // [ bonds ], [ constraints ]: ai aj ...
    SECTION_SETTLES,        // ow funct doh dhh -> ow-(ow+1), ow-(ow+2)
    SECTION_VSITE,          // site ai ... -> site-ai
    SECTION_VSITEN,         // site funct ai ... -> site-ai
    SECTION_MOLECULES,
    SECTION_SYSTEM,
} section_t;

typedef struct {
    md_itp_data_t*  data;
    md_allocator_i* alloc;
    md_allocator_i* temp;

    str_t* defines;         // md_array, temp

    // Conditional stack. value: what the #ifdef/#ifndef tested. parent: whether the enclosing block was active.
    bool cond_value [ITP_MAX_COND_DEPTH];
    bool cond_else  [ITP_MAX_COND_DEPTH];
    bool cond_parent[ITP_MAX_COND_DEPTH];
    int  cond_depth;

    section_t section;
    int64_t   mol;          // Current moleculetype, -1 when none (or after [ intermolecular_interactions ])
    int       include_depth;
    bool      error;
} itp_parser_t;

static bool parser_active(const itp_parser_t* p) {
    if (p->cond_depth == 0) return true;
    const int d = p->cond_depth - 1;
    return p->cond_parent[d] && (p->cond_else[d] ? !p->cond_value[d] : p->cond_value[d]);
}

static bool is_defined(const itp_parser_t* p, str_t name) {
    for (size_t i = 0; i < md_array_size(p->defines); ++i) {
        if (str_eq(p->defines[i], name)) return true;
    }
    return false;
}

static section_t section_from_name(str_t name) {
    static const struct { const char* name; section_t section; } table[] = {
        {"atomtypes",       SECTION_ATOMTYPES},
        {"moleculetype",    SECTION_MOLECULETYPE},
        {"atoms",           SECTION_ATOMS},
        {"bonds",           SECTION_BONDS},
        {"constraints",     SECTION_BONDS},
        {"settles",         SECTION_SETTLES},
        {"virtual_sites1",  SECTION_VSITE},
        {"virtual_sites2",  SECTION_VSITE},
        {"virtual_sites3",  SECTION_VSITE},
        {"virtual_sites4",  SECTION_VSITE},
        {"dummies2",        SECTION_VSITE},
        {"dummies3",        SECTION_VSITE},
        {"dummies4",        SECTION_VSITE},
        {"virtual_sitesn",  SECTION_VSITEN},
        {"dummiesn",        SECTION_VSITEN},
        {"molecules",       SECTION_MOLECULES},
        {"system",          SECTION_SYSTEM},
    };
    for (size_t i = 0; i < ARRAY_SIZE(table); ++i) {
        if (str_eq_cstr_ignore_case(name, table[i].name)) return table[i].section;
    }
    return SECTION_IGNORE;
}

static md_itp_moleculetype_t* current_mol(itp_parser_t* p) {
    if (p->mol < 0) return NULL;
    return &p->data->moleculetypes[p->mol];
}

static void add_local_bond(itp_parser_t* p, int64_t a1, int64_t b1) {
    md_itp_moleculetype_t* mol = current_mol(p);
    if (!mol) return;
    const int64_t n = (int64_t)md_array_size(mol->atoms);
    if (a1 < 1 || b1 < 1 || a1 > n || b1 > n || a1 == b1) {
        MD_LOG_INFO("itp: skipping connection %lld-%lld outside of moleculetype '" STR_FMT "' (%lld atoms)", (long long)a1, (long long)b1, STR_ARG(mol->name), (long long)n);
        return;
    }
    md_atom_pair_t pair = {{ (md_atom_idx_t)(MIN(a1, b1) - 1), (md_atom_idx_t)(MAX(a1, b1) - 1) }};
    md_array_push(mol->bonds, pair, p->alloc);
}

static bool parse_i64(int64_t* out, str_t tok) {
    if (!is_int(tok)) return false;
    *out = parse_int(tok);
    return true;
}

// parse_float does not take a leading '+', which is_float accepts and topologies do write
static float parse_f32(str_t tok) {
    if (tok.len > 1 && tok.ptr[0] == '+') {
        tok.ptr += 1;
        tok.len -= 1;
    }
    return (float)parse_float(tok);
}

static void parse_data_line(itp_parser_t* p, str_t line) {
    str_t tok[ITP_MAX_TOKENS];
    str_t rest = line;
    const size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &rest);
    if (num_tok == 0) return;

    switch (p->section) {
    case SECTION_ATOMTYPES: {
        // Column layouts differ between force fields:
        //   name mass charge ptype ...
        //   name at.num mass charge ptype ...
        //   name bond_type at.num mass charge ptype ...
        // The particle type (A, S, V, D, B) is the anchor; mass and charge precede it.
        size_t pt = 0;
        for (size_t i = 3; i < num_tok; ++i) {
            if (tok[i].len == 1 && (tok[i].ptr[0] == 'A' || tok[i].ptr[0] == 'S' || tok[i].ptr[0] == 'V' || tok[i].ptr[0] == 'D' || tok[i].ptr[0] == 'B')) {
                pt = i;
                break;
            }
        }
        if (pt == 0 || !is_float(tok[pt - 1]) || !is_float(tok[pt - 2])) {
            MD_LOG_INFO("itp: could not interpret atomtype line '" STR_FMT "'", STR_ARG(line));
            return;
        }
        md_itp_atomtype_t type = {
            .name   = str_copy(tok[0], p->alloc),
            .atomic_number = -1,
            .mass   = parse_f32(tok[pt - 2]),
            .charge = parse_f32(tok[pt - 1]),
        };
        if (pt >= 4 && is_int(tok[pt - 3])) {
            type.atomic_number = (int32_t)parse_int(tok[pt - 3]);
        }
        md_array_push(p->data->atomtypes, type, p->alloc);
        break;
    }
    case SECTION_MOLECULETYPE: {
        md_itp_moleculetype_t mol = {
            .name   = str_copy(tok[0], p->alloc),
            .nrexcl = num_tok > 1 ? (int32_t)parse_int(tok[1]) : 0,
        };
        md_array_push(p->data->moleculetypes, mol, p->alloc);
        p->mol = (int64_t)md_array_size(p->data->moleculetypes) - 1;
        // A moleculetype has exactly one line, what follows belongs to its subsections
        p->section = SECTION_IGNORE;
        break;
    }
    case SECTION_ATOMS: {
        md_itp_moleculetype_t* mol = current_mol(p);
        if (!mol) return;
        int64_t nr = 0;
        if (num_tok < 5 || !parse_i64(&nr, tok[0])) {
            MD_LOG_INFO("itp: could not interpret atom line '" STR_FMT "'", STR_ARG(line));
            return;
        }
        if (nr != (int64_t)md_array_size(mol->atoms) + 1) {
            MD_LOG_ERROR("itp: atoms of moleculetype '" STR_FMT "' are not numbered consecutively (got %lld, expected %zu)", STR_ARG(mol->name), (long long)nr, md_array_size(mol->atoms) + 1);
            p->error = true;
            return;
        }
        md_itp_atom_t atom = {
            .type     = str_copy(tok[1], p->alloc),
            .res_nr   = (int32_t)parse_int(tok[2]),
            .res_name = str_copy(tok[3], p->alloc),
            .name     = str_copy(tok[4], p->alloc),
        };
        if (num_tok > 6 && is_float(tok[6])) {
            atom.charge = parse_f32(tok[6]);
            atom.has_charge = true;
        }
        if (num_tok > 7 && is_float(tok[7])) {
            atom.mass = parse_f32(tok[7]);
            atom.has_mass = true;
        }
        md_array_push(mol->atoms, atom, p->alloc);
        break;
    }
    case SECTION_BONDS: {
        int64_t a, b;
        if (num_tok >= 2 && parse_i64(&a, tok[0]) && parse_i64(&b, tok[1])) {
            add_local_bond(p, a, b);
        }
        break;
    }
    case SECTION_SETTLES: {
        int64_t o;
        if (parse_i64(&o, tok[0])) {
            add_local_bond(p, o, o + 1);
            add_local_bond(p, o, o + 2);
        }
        break;
    }
    case SECTION_VSITE: {
        int64_t s, a;
        if (num_tok >= 2 && parse_i64(&s, tok[0]) && parse_i64(&a, tok[1])) {
            add_local_bond(p, s, a);
        }
        break;
    }
    case SECTION_VSITEN: {
        int64_t s, a;
        if (num_tok >= 3 && parse_i64(&s, tok[0]) && parse_i64(&a, tok[2])) {
            add_local_bond(p, s, a);
        }
        break;
    }
    case SECTION_MOLECULES: {
        int64_t count = 0;
        if (num_tok >= 2 && parse_i64(&count, tok[1])) {
            md_itp_molecules_t entry = { .name = str_copy(tok[0], p->alloc), .count = (int32_t)count };
            md_array_push(p->data->molecules, entry, p->alloc);
        }
        break;
    }
    case SECTION_SYSTEM:
        if (str_empty(p->data->system_name)) {
            p->data->system_name = str_copy(line, p->alloc);
        }
        break;
    default:
        break;
    }
}

static bool parse_text(itp_parser_t* p, str_t text, str_t folder);

static void parse_include(itp_parser_t* p, str_t arg, str_t folder) {
    arg = str_trim(arg);
    if (arg.len < 2) return;
    const char open = arg.ptr[0];
    const char close = open == '"' ? '"' : (open == '<' ? '>' : 0);
    if (!close) {
        MD_LOG_INFO("itp: malformed #include '" STR_FMT "'", STR_ARG(arg));
        return;
    }
    size_t end = 0;
    str_t inner = str_substr(arg, 1, SIZE_MAX);
    if (!str_find_char(&end, inner, close)) {
        MD_LOG_INFO("itp: malformed #include '" STR_FMT "'", STR_ARG(arg));
        return;
    }
    str_t file = str_substr(inner, 0, end);

    if (p->include_depth >= ITP_MAX_INCLUDE_DEPTH) {
        MD_LOG_ERROR("itp: #include nested deeper than %d, skipping '" STR_FMT "'", ITP_MAX_INCLUDE_DEPTH, STR_ARG(file));
        return;
    }

    md_strb_t sb = md_strb_create(p->temp);
    str_t resolved = {0};
    if (md_path_is_absolute(file)) {
        if (md_path_is_valid(file)) resolved = file;
    } else {
        if (!str_empty(folder)) {
            md_strb_push_str(&sb, folder);
            md_strb_push_str(&sb, file);
            str_t cand = md_strb_to_str(sb);
            if (md_path_is_valid(cand)) resolved = cand;
        }
        const char* gmxlib = getenv("GMXLIB");
        if (str_empty(resolved) && gmxlib) {
            md_strb_reset(&sb);
            md_strb_push_cstr(&sb, gmxlib);
            md_strb_push_char(&sb, '/');
            md_strb_push_str(&sb, file);
            str_t cand = md_strb_to_str(sb);
            if (md_path_is_valid(cand)) resolved = cand;
        }
    }

    if (str_empty(resolved)) {
        // Common in a .top: the force field lives in the GROMACS share folder. Its atomtypes are only a
        // fallback for masses and charges, so this is worth a note and not a failure.
        MD_LOG_INFO("itp: could not resolve #include '" STR_FMT "', skipping it", STR_ARG(file));
        return;
    }

    str_t text = load_textfile(resolved, p->temp);
    if (str_empty(text)) {
        MD_LOG_INFO("itp: #include '" STR_FMT "' is empty or unreadable", STR_ARG(resolved));
        return;
    }
    str_t sub_folder = {0};
    extract_folder_path(&sub_folder, resolved);
    // The included file's sections do not leak into the including one's position, but its current
    // moleculetype does: a molecule .itp is commonly followed by position restraints in the parent.
    p->include_depth += 1;
    parse_text(p, text, str_copy(sub_folder, p->temp));
    p->include_depth -= 1;
}

static void parse_directive(itp_parser_t* p, str_t line, str_t folder) {
    str_t rest = str_trim(str_substr(line, 1, SIZE_MAX));
    str_t word = {0};
    str_t args = rest;
    extract_token(&word, &args);
    str_t arg0 = {0};
    {
        str_t tmp = args;
        extract_token(&arg0, &tmp);
    }

    if (str_eq_cstr(word, "ifdef") || str_eq_cstr(word, "ifndef")) {
        if (p->cond_depth >= ITP_MAX_COND_DEPTH) {
            MD_LOG_ERROR("itp: conditionals nested deeper than %d", ITP_MAX_COND_DEPTH);
            p->error = true;
            return;
        }
        const bool defined = is_defined(p, arg0);
        const int d = p->cond_depth;
        p->cond_parent[d] = parser_active(p);
        p->cond_value[d]  = str_eq_cstr(word, "ifdef") ? defined : !defined;
        p->cond_else[d]   = false;
        p->cond_depth += 1;
    } else if (str_eq_cstr(word, "else")) {
        if (p->cond_depth == 0) {
            MD_LOG_ERROR("itp: #else without #ifdef");
            p->error = true;
            return;
        }
        p->cond_else[p->cond_depth - 1] = true;
    } else if (str_eq_cstr(word, "endif")) {
        if (p->cond_depth == 0) {
            MD_LOG_ERROR("itp: #endif without #ifdef");
            p->error = true;
            return;
        }
        p->cond_depth -= 1;
    } else if (!parser_active(p)) {
        return;
    } else if (str_eq_cstr(word, "define")) {
        if (!str_empty(arg0) && !is_defined(p, arg0)) {
            md_array_push(p->defines, str_copy(arg0, p->temp), p->temp);
        }
    } else if (str_eq_cstr(word, "undef")) {
        for (size_t i = 0; i < md_array_size(p->defines); ++i) {
            if (str_eq(p->defines[i], arg0)) {
                p->defines[i] = *md_array_last(p->defines);
                md_array_pop(p->defines);
                break;
            }
        }
    } else if (str_eq_cstr(word, "include")) {
        parse_include(p, args, folder);
    } else {
        MD_LOG_INFO("itp: ignoring unsupported directive '" STR_FMT "'", STR_ARG(line));
    }
}

static bool parse_text(itp_parser_t* p, str_t text, str_t folder) {
    md_strb_t joined = md_strb_create(p->temp);
    bool continuing = false;

    str_t raw;
    while (str_extract_line(&raw, &text)) {
        str_t line = str_trim_end(raw);

        // Line continuation, before comments like the C preprocessor it imitates
        if (line.len > 0 && line.ptr[line.len - 1] == '\\') {
            md_strb_push_str(&joined, str_substr(line, 0, line.len - 1));
            md_strb_push_char(&joined, ' ');
            continuing = true;
            continue;
        }
        if (continuing) {
            md_strb_push_str(&joined, line);
            line = md_strb_to_str(joined);
            continuing = false;
        }

        size_t loc = 0;
        if (str_find_char(&loc, line, ';')) {
            line = str_substr(line, 0, loc);
        }
        line = str_trim(line);

        if (line.len > 0) {
            if (line.ptr[0] == '#') {
                parse_directive(p, line, folder);
            } else if (parser_active(p)) {
                if (line.ptr[0] == '[') {
                    size_t close = 0;
                    if (str_find_char(&close, line, ']')) {
                        str_t name = str_trim(str_substr(line, 1, close - 1));
                        p->section = section_from_name(name);
                        if (str_eq_cstr_ignore_case(name, "intermolecular_interactions")) {
                            // Global atom indices from here on, nothing that belongs to a moleculetype
                            p->mol = -1;
                        }
                        if (p->section == SECTION_MOLECULES || p->section == SECTION_SYSTEM) {
                            p->mol = -1;
                        }
                    } else {
                        MD_LOG_INFO("itp: malformed section header '" STR_FMT "'", STR_ARG(line));
                        p->section = SECTION_IGNORE;
                    }
                } else {
                    parse_data_line(p, line);
                }
            }
        }

        if (!continuing) {
            md_strb_reset(&joined);
        }
        if (p->error) return false;
    }
    return true;
}

static int compare_pair(const void* a, const void* b) {
    const md_atom_pair_t* pa = (const md_atom_pair_t*)a;
    const md_atom_pair_t* pb = (const md_atom_pair_t*)b;
    if (pa->idx[0] != pb->idx[0]) return pa->idx[0] < pb->idx[0] ? -1 : 1;
    if (pa->idx[1] != pb->idx[1]) return pa->idx[1] < pb->idx[1] ? -1 : 1;
    return 0;
}

bool md_itp_data_parse_str(md_itp_data_t* data, str_t str, str_t base_folder, md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(alloc);

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    itp_parser_t p = {
        .data  = data,
        .alloc = alloc,
        .temp  = md_temp_allocator(temp),
        .mol   = -1,
    };

    bool ok = parse_text(&p, str, base_folder);
    if (ok && p.cond_depth != 0) {
        MD_LOG_ERROR("itp: missing #endif");
        ok = false;
    }

    // Bonds, constraints and settles may name the same pair more than once
    for (size_t i = 0; i < md_array_size(data->moleculetypes); ++i) {
        md_itp_moleculetype_t* mol = &data->moleculetypes[i];
        size_t n = md_array_size(mol->bonds);
        if (n > 1) {
            qsort(mol->bonds, n, sizeof(md_atom_pair_t), compare_pair);
            size_t w = 1;
            for (size_t r = 1; r < n; ++r) {
                if (compare_pair(&mol->bonds[r], &mol->bonds[w - 1]) != 0) {
                    mol->bonds[w++] = mol->bonds[r];
                }
            }
            md_array_shrink(mol->bonds, w);
        }
    }

    md_temp_end(temp);
    return ok;
}

bool md_itp_data_parse_file(md_itp_data_t* data, str_t filename, md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(alloc);

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    bool ok = false;
    str_t text = load_textfile(filename, md_temp_allocator(temp));
    if (str_empty(text)) {
        MD_LOG_ERROR("Could not open topology file '" STR_FMT "'", STR_ARG(filename));
    } else {
        str_t folder = {0};
        extract_folder_path(&folder, filename);
        ok = md_itp_data_parse_str(data, text, folder, alloc);
    }
    md_temp_end(temp);
    return ok;
}

void md_itp_data_free(md_itp_data_t* data, md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(alloc);
    for (size_t i = 0; i < md_array_size(data->atomtypes); ++i) {
        str_free(data->atomtypes[i].name, alloc);
    }
    for (size_t i = 0; i < md_array_size(data->moleculetypes); ++i) {
        md_itp_moleculetype_t* mol = &data->moleculetypes[i];
        for (size_t j = 0; j < md_array_size(mol->atoms); ++j) {
            str_free(mol->atoms[j].type, alloc);
            str_free(mol->atoms[j].name, alloc);
            str_free(mol->atoms[j].res_name, alloc);
        }
        str_free(mol->name, alloc);
        md_array_free(mol->atoms, alloc);
        md_array_free(mol->bonds, alloc);
    }
    for (size_t i = 0; i < md_array_size(data->molecules); ++i) {
        str_free(data->molecules[i].name, alloc);
    }
    if (!str_empty(data->system_name)) str_free(data->system_name, alloc);
    md_array_free(data->atomtypes, alloc);
    md_array_free(data->moleculetypes, alloc);
    md_array_free(data->molecules, alloc);
    MEMSET(data, 0, sizeof(*data));
}

// ### MATCHING ###

static inline bool names_match(str_t sys_name, str_t top_name) {
    if (str_eq(sys_name, top_name)) return true;
    // Structure formats truncate: gro keeps 5 characters, pdb 4 (3 for residues), the system label 6
    return sys_name.len >= 3 && top_name.len > sys_name.len && str_begins_with(top_name, sys_name);
}

typedef struct {
    const md_system_t* sys;
    int32_t* atom_comp;     // [atom] -> component, NULL when the system has no components
} match_ctx_t;

// Atom names alone are ambiguous for repeating polymers: two consecutive 50-mers carry exactly the names
// of one 100-mer. So a molecule must also agree with the system's residues: it starts and ends on a
// residue boundary, and its own residue changes (resnr or residue name) fall exactly where the system's do.
static bool moleculetype_matches_at(const match_ctx_t* ctx, const md_itp_moleculetype_t* mol, size_t offset) {
    const md_system_t* sys = ctx->sys;
    const size_t n = md_array_size(mol->atoms);
    if (n == 0 || offset + n > sys->atom.count) return false;

    for (size_t i = 0; i < n; ++i) {
        if (!names_match(md_atom_name(&sys->atom, offset + i), mol->atoms[i].name)) return false;
    }

    if (ctx->atom_comp) {
        const int32_t* ac = ctx->atom_comp;
        if (offset > 0 && ac[offset - 1] == ac[offset]) return false;
        if (offset + n < sys->atom.count && ac[offset + n - 1] == ac[offset + n]) return false;
        for (size_t i = 0; i < n; ++i) {
            const bool top_new = i > 0 && (mol->atoms[i].res_nr != mol->atoms[i - 1].res_nr || !str_eq(mol->atoms[i].res_name, mol->atoms[i - 1].res_name));
            const bool sys_new = i > 0 && ac[offset + i] != ac[offset + i - 1];
            if (top_new != sys_new) return false;
            if (i == 0 || sys_new) {
                if (!names_match(md_component_name(&sys->component, (size_t)ac[offset + i]), mol->atoms[i].res_name)) return false;
            }
        }
    }
    return true;
}

static int64_t find_moleculetype(const md_itp_data_t* data, str_t name) {
    for (size_t i = 0; i < md_array_size(data->moleculetypes); ++i) {
        if (str_eq(data->moleculetypes[i].name, name)) return (int64_t)i;
    }
    for (size_t i = 0; i < md_array_size(data->moleculetypes); ++i) {
        if (str_eq_ignore_case(data->moleculetypes[i].name, name)) return (int64_t)i;
    }
    return -1;
}

static bool match_sequential(md_itp_instance_t** out, const match_ctx_t* ctx, const md_itp_data_t* data, md_allocator_i* alloc) {
    size_t offset = 0;
    for (size_t e = 0; e < md_array_size(data->molecules); ++e) {
        const md_itp_molecules_t* entry = &data->molecules[e];
        const int64_t type = find_moleculetype(data, entry->name);
        if (type < 0) {
            MD_LOG_INFO("itp: [ molecules ] names '" STR_FMT "', which has no [ moleculetype ]", STR_ARG(entry->name));
            return false;
        }
        const md_itp_moleculetype_t* mol = &data->moleculetypes[type];
        for (int32_t c = 0; c < entry->count; ++c) {
            if (!moleculetype_matches_at(ctx, mol, offset)) {
                MD_LOG_INFO("itp: molecule %d of '" STR_FMT "' does not match the system at atom %zu", c + 1, STR_ARG(entry->name), offset + 1);
                return false;
            }
            md_itp_instance_t inst = { (uint32_t)type, (uint32_t)offset };
            md_array_push(*out, inst, alloc);
            offset += md_array_size(mol->atoms);
        }
    }
    if (offset != ctx->sys->atom.count) {
        MD_LOG_INFO("itp: [ molecules ] covers %zu atoms, the system has %zu", offset, ctx->sys->atom.count);
    }
    return true;
}

static void match_scan(md_itp_instance_t** out, const match_ctx_t* ctx, const md_itp_data_t* data, md_allocator_i* alloc) {
    const md_system_t* sys = ctx->sys;
    const size_t num_types = md_array_size(data->moleculetypes);
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);

    // Longest first, so that a molecule is not claimed by a smaller moleculetype that happens to fit its start
    size_t* order = md_temp_alloc_array(temp, size_t, num_types);
    for (size_t i = 0; i < num_types; ++i) order[i] = i;
    for (size_t i = 1; i < num_types; ++i) {
        size_t k = order[i];
        size_t j = i;
        while (j > 0 && md_array_size(data->moleculetypes[order[j - 1]].atoms) < md_array_size(data->moleculetypes[k].atoms)) {
            order[j] = order[j - 1];
            --j;
        }
        order[j] = k;
    }

    size_t offset = 0;
    while (offset < sys->atom.count) {
        size_t step = 1;
        for (size_t o = 0; o < num_types; ++o) {
            const md_itp_moleculetype_t* mol = &data->moleculetypes[order[o]];
            if (moleculetype_matches_at(ctx, mol, offset)) {
                md_itp_instance_t inst = { (uint32_t)order[o], (uint32_t)offset };
                md_array_push(*out, inst, alloc);
                step = md_array_size(mol->atoms);
                break;
            }
        }
        if (step == 1 && ctx->atom_comp) {
            // Nothing starts here: a molecule can only start on a residue boundary, so skip the rest of this residue
            const int32_t c = ctx->atom_comp[offset];
            while (offset + step < sys->atom.count && ctx->atom_comp[offset + step] == c) ++step;
        }
        offset += step;
    }

    md_temp_end(temp);
}

size_t md_itp_match_system(md_itp_instance_t** out_instances, const md_itp_data_t* data, const md_system_t* sys, md_allocator_i* alloc) {
    ASSERT(out_instances);
    ASSERT(data);
    ASSERT(sys);
    ASSERT(alloc);

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    match_ctx_t ctx = { .sys = sys };
    if (sys->component.count > 0 && sys->component.atom_offset) {
        ctx.atom_comp = md_temp_alloc_array(temp, int32_t, sys->atom.count);
        for (size_t i = 0; i < sys->atom.count; ++i) ctx.atom_comp[i] = -1;
        for (size_t c = 0; c < sys->component.count; ++c) {
            md_urange_t range = md_component_atom_range(&sys->component, c);
            for (uint32_t i = range.beg; i < range.end && i < sys->atom.count; ++i) ctx.atom_comp[i] = (int32_t)c;
        }
    }

    const size_t beg = md_array_size(*out_instances);
    bool done = false;
    if (md_array_size(data->molecules) > 0) {
        done = match_sequential(out_instances, &ctx, data, alloc);
        if (!done) {
            MD_LOG_INFO("itp: [ molecules ] does not describe the system in order, placing molecules by atom and residue names instead");
            md_array_shrink(*out_instances, beg);
        }
    }
    if (!done) {
        match_scan(out_instances, &ctx, data, alloc);
    }

    md_temp_end(temp);
    return md_array_size(*out_instances) - beg;
}

// ### APPLICATION ###

static const md_itp_atomtype_t* find_atomtype(const md_itp_data_t* data, str_t name) {
    for (size_t i = 0; i < md_array_size(data->atomtypes); ++i) {
        if (str_eq(data->atomtypes[i].name, name)) return &data->atomtypes[i];
    }
    return NULL;
}

bool md_itp_system_supplement(md_system_t* sys, const md_itp_data_t* data) {
    ASSERT(sys);
    ASSERT(data);

    if (!sys->alloc || sys->atom.count == 0) {
        MD_LOG_ERROR("A topology can only supplement a loaded system");
        return false;
    }

    const size_t atom_count = sys->atom.count;
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp);
    bool result = false;

    md_itp_instance_t* instances = 0;
    const size_t num_instances = md_itp_match_system(&instances, data, sys, temp_arena);
    if (num_instances == 0) {
        MD_LOG_ERROR("None of the topology's molecules match the loaded system");
        goto done;
    }

    // Charge and mass per moleculetype atom, resolved against the atomtypes once rather than per instance
    const size_t num_types = md_array_size(data->moleculetypes);
    float** type_charge = md_temp_alloc_array(temp, float*, num_types);
    float** type_mass   = md_temp_alloc_array(temp, float*, num_types);
    for (size_t t = 0; t < num_types; ++t) {
        const md_itp_moleculetype_t* mol = &data->moleculetypes[t];
        const size_t n = md_array_size(mol->atoms);
        type_charge[t] = md_temp_alloc_array(temp, float, n);
        type_mass[t]   = md_temp_alloc_array(temp, float, n);
        for (size_t i = 0; i < n; ++i) {
            const md_itp_atom_t* atom = &mol->atoms[i];
            const md_itp_atomtype_t* at = (!atom->has_charge || !atom->has_mass) ? find_atomtype(data, atom->type) : NULL;
            type_charge[t][i] = atom->has_charge ? atom->charge : (at ? at->charge : NAN);
            type_mass[t][i]   = atom->has_mass   ? atom->mass   : (at ? at->mass   : NAN);
        }
    }

    uint8_t* covered = md_temp_alloc_array(temp, uint8_t, atom_count);
    float*   charge  = md_temp_alloc_array(temp, float, atom_count);
    float*   mass    = md_temp_alloc_array(temp, float, atom_count);
    MEMSET(covered, 0, atom_count);
    for (size_t i = 0; i < atom_count; ++i) {
        charge[i] = NAN;
        mass[i]   = NAN;
    }

    size_t num_topology_bonds = 0;
    size_t num_covered = 0;
    for (size_t k = 0; k < num_instances; ++k) {
        const uint32_t t = instances[k].type;
        const md_itp_moleculetype_t* mol = &data->moleculetypes[t];
        const size_t off = instances[k].atom_offset;
        const size_t n = md_array_size(mol->atoms);
        for (size_t i = 0; i < n; ++i) {
            covered[off + i] = 1;
            charge[off + i]  = type_charge[t][i];
            mass[off + i]    = type_mass[t][i];
        }
        num_covered += n;
        num_topology_bonds += md_array_size(mol->bonds);
    }

    // ## Bonds
    // Laid out as [kept] [kept topology] [new topology] [user defined]. Re-inference keeps every bond without
    // MD_BOND_FLAG_INFERRED wherever it sits, so the order is for readability, not a contract.
    {
        const size_t old_count = sys->bond.count;
        const size_t cap = old_count + num_topology_bonds;
        md_atom_pair_t*  pairs = md_temp_alloc_array(temp, md_atom_pair_t,  cap);
        md_bond_flags_t* flags = md_temp_alloc_array(temp, md_bond_flags_t, cap);
        size_t count = 0;

        for (int pass = 0; pass < 2; ++pass) {
            for (size_t i = 0; i < old_count; ++i) {
                const md_atom_pair_t  pair = sys->bond.pairs[i];
                const md_bond_flags_t f    = sys->bond.flags[i];
                if (f & MD_BOND_FLAG_USER_DEFINED) continue;
                if (((f & MD_BOND_FLAG_TOPOLOGY) != 0) != (pass == 1)) continue;
                const bool valid = pair.idx[0] >= 0 && pair.idx[1] >= 0 && (size_t)pair.idx[0] < atom_count && (size_t)pair.idx[1] < atom_count;
                if (valid && covered[pair.idx[0]] && covered[pair.idx[1]]) continue;  // Replaced by the topology
                pairs[count] = pair;
                flags[count] = f;
                count += 1;
            }
        }

        for (size_t k = 0; k < num_instances; ++k) {
            const md_itp_moleculetype_t* mol = &data->moleculetypes[instances[k].type];
            const md_atom_idx_t off = (md_atom_idx_t)instances[k].atom_offset;
            for (size_t b = 0; b < md_array_size(mol->bonds); ++b) {
                pairs[count].idx[0] = off + mol->bonds[b].idx[0];
                pairs[count].idx[1] = off + mol->bonds[b].idx[1];
                flags[count] = MD_BOND_FLAG_COVALENT | MD_BOND_FLAG_TOPOLOGY;
                count += 1;
            }
        }

        size_t num_user = 0;
        for (size_t i = 0; i < old_count; ++i) {
            if (sys->bond.flags[i] & MD_BOND_FLAG_USER_DEFINED) num_user += 1;
        }
        md_atom_pair_t*  user_pairs = md_temp_alloc_array(temp, md_atom_pair_t,  num_user + 1);
        md_bond_flags_t* user_flags = md_temp_alloc_array(temp, md_bond_flags_t, num_user + 1);
        for (size_t i = 0, u = 0; i < old_count; ++i) {
            if (sys->bond.flags[i] & MD_BOND_FLAG_USER_DEFINED) {
                user_pairs[u] = sys->bond.pairs[i];
                user_flags[u] = sys->bond.flags[i];
                u += 1;
            }
        }

        md_bond_data_clear(&sys->bond);
        md_array_resize(sys->bond.pairs, count + num_user, sys->alloc);
        md_array_resize(sys->bond.flags, count + num_user, sys->alloc);
        if (count) {
            MEMCPY(sys->bond.pairs, pairs, sizeof(md_atom_pair_t) * count);
            MEMCPY(sys->bond.flags, flags, sizeof(md_bond_flags_t) * count);
        }
        if (num_user) {
            MEMCPY(sys->bond.pairs + count, user_pairs, sizeof(md_atom_pair_t) * num_user);
            MEMCPY(sys->bond.flags + count, user_flags, sizeof(md_bond_flags_t) * num_user);
        }
        sys->bond.count = count + num_user;
        md_bond_build_connectivity(&sys->bond, atom_count, sys->alloc);
    }

    // ## Charge and mass
    if (sys->attributes.alloc) {
        md_attributes_publish_atom_column(&sys->attributes, STR_LIT("atom/charge"), md_unit_elementary_charge(), 1, charge, atom_count);
        md_attributes_publish_atom_column(&sys->attributes, STR_LIT("atom/mass"),   md_unit_dalton(),            1, mass,   atom_count);
    }

    // ## Blank type masses
    // Only where the type has none, and only when every matched atom of the type agrees: a type is shared
    // by atoms the topology may give different masses (the end beads of a chain), and picking one of them
    // would be a guess dressed up as data. 'atom/mass' carries the per atom truth either way.
    if (sys->atom.type_idx && sys->atom.type.count > 0 && sys->atom.type.mass) {
        const size_t num_atom_types = sys->atom.type.count;
        float*   fill  = md_temp_alloc_array(temp, float,   num_atom_types);
        uint8_t* state = md_temp_alloc_array(temp, uint8_t, num_atom_types);  // 0 unseen, 1 agreeing, 2 conflicting
        MEMSET(state, 0, num_atom_types);
        for (size_t i = 0; i < atom_count; ++i) {
            if (!covered[i]) continue;
            const size_t t = sys->atom.type_idx[i];
            if (t >= num_atom_types || sys->atom.type.mass[t] != 0.0f) continue;
            const float m = mass[i];
            if (!(m > 0.0f)) {
                state[t] = 2;
            } else if (state[t] == 0) {
                fill[t] = m;
                state[t] = 1;
            } else if (state[t] == 1 && fill[t] != m) {
                state[t] = 2;
            }
        }
        for (size_t t = 0; t < num_atom_types; ++t) {
            if (state[t] == 1) sys->atom.type.mass[t] = fill[t];
        }
    }

    // ## Derived topology
    md_util_system_infer_structures(sys);
    if (sys->bond.count) {
        md_util_system_infer_rings(sys);
    }

    MD_LOG_INFO("Topology applied: %zu molecules covering %zu of %zu atoms, %zu bonds, %zu structures",
        num_instances, num_covered, atom_count, num_topology_bonds, md_structure_count(&sys->structure));
    result = true;

done:
    md_temp_end(temp);
    return result;
}

bool md_itp_system_supplement_from_file(md_system_t* sys, str_t filename) {
    ASSERT(sys);
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    md_itp_data_t data = {0};
    bool result = md_itp_data_parse_file(&data, filename, temp_arena) && md_itp_system_supplement(sys, &data);

    md_temp_end(temp);
    return result;
}
