#include <md_amber_prmtop.h>

#include <md_system.h>
#include <md_util.h>
#include <md_topology_build.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_hash.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_parse.h>
#include <core/md_str.h>
#include <core/md_unit.h>

#include <math.h>
#include <string.h>

// ## THE FORMAT
//
//     %VERSION  VERSION_STAMP = V0001.000  DATE = ...
//     %FLAG TITLE
//     %FORMAT(20a4)
//     ...
//     %FLAG POINTERS
//     %FORMAT(10I8)
//        NATOM  NTYPES   NBONH   MBONA  NTHETH  MTHETA   NPHIH   MPHIA  NHPARM   NPARM
//          NNB    NRES   NBONA  NTHETA   NPHIA  NUMBND  NUMANG  NPTRA   NATYP    NPHB
//       IFPERT   NBPER   NGPER   NDPER   MBPER   MGPER   MDPER   IFBOX   NMXRS   IFCAP
//     NUMEXTRA   NCOPY
//     %FLAG ATOM_NAME
//     ...
//
// Each block is a %FLAG, a %FORMAT in Fortran notation (count per line, kind, width) and its values in fixed width
// fields, as many lines as they take. A %COMMENT may come anywhere between. The fields are read by their width, as
// Fortran does: a name may hold any character, and nothing separates two numbers that fill their fields.

enum {
    POINTER_NATOM = 0,
    POINTER_NRES  = 11,
    POINTER_IFBOX = 27,
};

typedef struct prmtop_block_t {
    md_str_t    name;
    char     kind;          // 'a' text, 'i' integer, 'e' real (E, D and F alike)
    uint32_t per_line;
    uint32_t width;
    md_str_t    text;          // the value lines
} prmtop_block_t;

typedef struct prmtop_t {
    md_array(prmtop_block_t) blocks;
    bool chamber;
} prmtop_t;

static bool prmtop_parse_format(prmtop_block_t* b, md_str_t line) {
    // %FORMAT(20a4) %FORMAT(5E16.8) %FORMAT(10I8) %FORMAT(1a80)
    size_t open = 0;
    if (!md_str_find_char(&open, line, '(')) return false;
    md_str_t s = md_str_substr(line, open + 1, SIZE_MAX);
    int64_t count = 0, width = 0;
    const size_t n = md_parse_i64(&count, s);
    s = md_str_substr(s, n, SIZE_MAX);     // no count: one per line
    if (n == 0) count = 1;
    if (s.len == 0) return false;
    const char k = (char)(s.ptr[0] | 0x20);
    s = md_str_substr(s, 1, SIZE_MAX);
    const size_t w = md_parse_i64(&width, s);
    if (!w || count <= 0 || count > 1000 || width <= 0 || width > 1000) return false;
    switch (k) {
    case 'a': b->kind = 'a'; break;
    case 'i': b->kind = 'i'; break;
    case 'e': case 'd': case 'f': case 'g': b->kind = 'e'; break;
    default: return false;
    }
    b->per_line = (uint32_t)count;
    b->width = (uint32_t)width;
    return true;
}

static md_str_t prmtop_next_line(md_str_t* text) {
    size_t nl = 0;
    md_str_t line;
    if (!md_str_find_char(&nl, *text, '\n')) {
        line = *text;
        *text = (md_str_t){0};
    } else {
        line = md_str_substr(*text, 0, nl);
        *text = md_str_substr(*text, nl + 1, SIZE_MAX);
    }
    if (line.len && line.ptr[line.len - 1] == '\r') line.len -= 1;
    return line;
}

// Finds the blocks of the file. The values are left where they are, read when asked for.
static bool prmtop_parse(prmtop_t* top, md_str_t text, md_allocator_i* alloc, md_str_t path) {
    MD_MEMSET(top, 0, sizeof(*top));
    const char* text_end = text.ptr + text.len;
    prmtop_block_t* cur = NULL;
    const char* data_beg = NULL;
    bool seen_flag = false;
    while (text.len) {
        const char* line_ptr = text.ptr;
        const md_str_t line = prmtop_next_line(&text);
        if (line.len && line.ptr[0] == '%') {
            // A % line ends the values of the block before it
            if (cur && data_beg) {
                cur->text = (md_str_t){data_beg, (size_t)(line_ptr - data_beg)};
                data_beg = NULL;
            }
            if (md_str_begins_with(line, MD_STR_LIT("%FLAG"))) {
                seen_flag = true;
                prmtop_block_t b = { .name = md_str_trim(md_str_substr(line, 5, SIZE_MAX)) };
                md_array_push(top->blocks, b, alloc);
                cur = md_array_last(top->blocks);
                if (md_str_eq_cstr(cur->name, "CTITLE") || md_str_begins_with(cur->name, MD_STR_LIT("CHARMM_"))) {
                    top->chamber = true;
                }
            } else if (md_str_begins_with(line, MD_STR_LIT("%FORMAT"))) {
                if (!cur || !prmtop_parse_format(cur, line)) {
                    MD_LOG_ERROR("Amber topology: '" MD_STR_FMT "' has a %%FORMAT that cannot be read: '" MD_STR_FMT "'", MD_STR_ARG(path), MD_STR_ARG(line));
                    return false;
                }
                data_beg = text.ptr;
            }
            // %VERSION, %COMMENT: nothing to keep
        } else if (!seen_flag && md_str_trim(line).len) {
            break;  // not this format
        }
    }
    if (cur && data_beg) {
        cur->text = (md_str_t){data_beg, (size_t)(text_end - data_beg)};
    }
    if (!seen_flag) {
        MD_LOG_ERROR("Amber topology: '" MD_STR_FMT "' is not a %%FLAG/%%FORMAT topology (Amber 7 or later)", MD_STR_ARG(path));
        return false;
    }
    return true;
}

static const prmtop_block_t* prmtop_block(const prmtop_t* top, const char* name) {
    for (size_t i = 0; i < md_array_size(top->blocks); ++i) {
        if (md_str_eq_cstr(top->blocks[i].name, name) && top->blocks[i].width) return &top->blocks[i];
    }
    return NULL;
}

// Calls back per field, in order. Blank fields are skipped: they pad a last line, or stand for nothing.
typedef bool (*prmtop_field_fn)(md_str_t field, size_t index, void* user);

static size_t prmtop_for_each(const prmtop_block_t* b, prmtop_field_fn fn, void* user, size_t max) {
    md_str_t text = b->text;
    size_t index = 0;
    while (text.len && index < max) {
        const md_str_t line = prmtop_next_line(&text);
        for (uint32_t k = 0; k < b->per_line && index < max; ++k) {
            const size_t beg = (size_t)k * b->width;
            if (beg >= line.len) break;
            const md_str_t field = md_str_substr(line, beg, b->width);
            if (b->kind != 'a' && md_str_empty(md_str_trim(field))) continue;
            if (!fn(field, index, user)) return SIZE_MAX;
            index += 1;
        }
    }
    return index;
}

static bool count_field(md_str_t field, size_t index, void* user) { (void)field; (void)index; (void)user; return true; }

static size_t prmtop_count(const prmtop_block_t* b) {
    return b ? prmtop_for_each(b, count_field, NULL, SIZE_MAX) : 0;
}

static bool int_field(md_str_t field, size_t index, void* user) {
    int64_t* out = (int64_t*)user;
    field = md_str_trim(field);
    return md_parse_i64(&out[index], field) == field.len;
}

static bool real_field(md_str_t field, size_t index, void* user) {
    double* out = (double*)user;
    field = md_str_trim(field);
    char buf[64];
    if (field.len >= sizeof(buf)) return false;
    // Fortran's D exponent
    for (size_t i = 0; i < field.len; ++i) buf[i] = (field.ptr[i] == 'D' || field.ptr[i] == 'd') ? 'E' : field.ptr[i];
    return md_parse_f64(&out[index], (md_str_t){buf, field.len}) == field.len && isfinite(out[index]);
}

static bool label_field(md_str_t field, size_t index, void* user) {
    md_str_t* out = (md_str_t*)user;
    out[index] = md_str_trim(field);
    return true;
}

// Exactly count values of the block, or false (and logged when the block is required)
static bool prmtop_ints(int64_t* out, const prmtop_t* top, const char* name, size_t count, md_str_t path) {
    const prmtop_block_t* b = prmtop_block(top, name);
    if (!b) return false;
    const size_t n = (b->kind == 'i') ? prmtop_for_each(b, int_field, out, count) : SIZE_MAX;
    if (n != count) {
        MD_LOG_ERROR("Amber topology: %s in '" MD_STR_FMT "' does not hold %zu integers", name, MD_STR_ARG(path), count);
        return false;
    }
    return true;
}

static bool prmtop_reals(double* out, const prmtop_t* top, const char* name, size_t count, md_str_t path) {
    const prmtop_block_t* b = prmtop_block(top, name);
    if (!b) return false;
    const size_t n = (b->kind != 'a') ? prmtop_for_each(b, real_field, out, count) : SIZE_MAX;
    if (n != count) {
        MD_LOG_ERROR("Amber topology: %s in '" MD_STR_FMT "' does not hold %zu numbers", name, MD_STR_ARG(path), count);
        return false;
    }
    return true;
}

static bool prmtop_labels(md_str_t* out, const prmtop_t* top, const char* name, size_t count, md_str_t path) {
    const prmtop_block_t* b = prmtop_block(top, name);
    if (!b) return false;
    const size_t n = (b->kind == 'a') ? prmtop_for_each(b, label_field, out, count) : SIZE_MAX;
    if (n != count) {
        MD_LOG_ERROR("Amber topology: %s in '" MD_STR_FMT "' does not hold %zu names", name, MD_STR_ARG(path), count);
        return false;
    }
    return true;
}

// ## WHAT A SYSTEM IS MADE OF, read and checked once for both building and supplementing

typedef struct prmtop_data_t {
    size_t num_atoms;
    size_t num_residues;
    size_t num_bonds;
    int64_t ifbox;
    md_str_t*   atom_name;
    md_str_t*   atom_type;         // NULL without AMBER_ATOM_TYPE
    double*  charge;            // e, NULL without CHARGE
    double*  mass;              // NULL without MASS
    int64_t* atomic_number;     // NULL without ATOMIC_NUMBER
    md_str_t*   res_name;
    uint32_t* res_offset;       // num_residues + 1
    int64_t* res_number;        // NULL without RESIDUE_NUMBER
    md_str_t*   res_chain;         // NULL without RESIDUE_CHAINID
    md_atom_pair_t* bonds;
    double   box[4];            // beta, a, b, c
    bool     has_box;
    md_str_t    title;
    bool     chamber;
} prmtop_data_t;

static bool prmtop_read(prmtop_data_t* d, md_str_t path, md_allocator_i* alloc) {
    MD_MEMSET(d, 0, sizeof(*d));
    const md_str_t text = md_load_textfile(path, alloc);
    if (md_str_empty(text)) {
        MD_LOG_ERROR("Amber topology: cannot read '" MD_STR_FMT "'", MD_STR_ARG(path));
        return false;
    }
    prmtop_t top;
    if (!prmtop_parse(&top, text, alloc, path)) return false;
    d->chamber = top.chamber;

    const prmtop_block_t* pb = prmtop_block(&top, "POINTERS");
    const size_t num_pointers = prmtop_count(pb);
    if (!pb || num_pointers < POINTER_NRES + 1) {
        MD_LOG_ERROR("Amber topology: '" MD_STR_FMT "' has no POINTERS", MD_STR_ARG(path));
        return false;
    }
    int64_t* pointers = md_alloc(alloc, num_pointers * sizeof(int64_t));
    if (!prmtop_ints(pointers, &top, "POINTERS", num_pointers, path)) return false;
    const int64_t natom = pointers[POINTER_NATOM];
    const int64_t nres  = pointers[POINTER_NRES];
    d->ifbox = num_pointers > POINTER_IFBOX ? pointers[POINTER_IFBOX] : 0;
    if (natom <= 0 || natom > INT32_MAX || nres <= 0 || nres > natom) {
        MD_LOG_ERROR("Amber topology: '" MD_STR_FMT "' has %lld atoms in %lld residues", MD_STR_ARG(path), (long long)natom, (long long)nres);
        return false;
    }
    const size_t N = (size_t)natom, R = (size_t)nres;
    d->num_atoms = N;
    d->num_residues = R;

    d->atom_name = md_alloc(alloc, N * sizeof(md_str_t));
    d->res_name  = md_alloc(alloc, R * sizeof(md_str_t));
    int64_t* res_ptr = md_alloc(alloc, R * sizeof(int64_t));
    if (!prmtop_labels(d->atom_name, &top, "ATOM_NAME", N, path) ||
        !prmtop_labels(d->res_name, &top, "RESIDUE_LABEL", R, path) ||
        !prmtop_ints(res_ptr, &top, "RESIDUE_POINTER", R, path)) {
        if (!prmtop_block(&top, "ATOM_NAME") || !prmtop_block(&top, "RESIDUE_LABEL") || !prmtop_block(&top, "RESIDUE_POINTER")) {
            MD_LOG_ERROR("Amber topology: '" MD_STR_FMT "' lacks ATOM_NAME, RESIDUE_LABEL or RESIDUE_POINTER", MD_STR_ARG(path));
        }
        return false;
    }
    d->res_offset = md_alloc(alloc, (R + 1) * sizeof(uint32_t));
    for (size_t r = 0; r < R; ++r) {
        const int64_t first = res_ptr[r] - 1;
        if (first < 0 || first >= natom || (r == 0 && first != 0) || (r > 0 && first <= (int64_t)d->res_offset[r - 1])) {
            MD_LOG_ERROR("Amber topology: RESIDUE_POINTER of residue %zu in '" MD_STR_FMT "' is out of order", r + 1, MD_STR_ARG(path));
            return false;
        }
        d->res_offset[r] = (uint32_t)first;
    }
    d->res_offset[R] = (uint32_t)N;

    // Optional, each checked for its count when there
    if (prmtop_block(&top, "CHARGE")) {
        d->charge = md_alloc(alloc, N * sizeof(double));
        if (!prmtop_reals(d->charge, &top, "CHARGE", N, path)) return false;
        for (size_t i = 0; i < N; ++i) d->charge[i] /= 18.2223;     // Amber's charge unit: e * sqrt(332.0522173)
    }
    if (prmtop_block(&top, "MASS")) {
        d->mass = md_alloc(alloc, N * sizeof(double));
        if (!prmtop_reals(d->mass, &top, "MASS", N, path)) return false;
    }
    if (prmtop_block(&top, "ATOMIC_NUMBER")) {
        d->atomic_number = md_alloc(alloc, N * sizeof(int64_t));
        if (!prmtop_ints(d->atomic_number, &top, "ATOMIC_NUMBER", N, path)) return false;
    }
    if (prmtop_block(&top, "AMBER_ATOM_TYPE")) {
        d->atom_type = md_alloc(alloc, N * sizeof(md_str_t));
        if (!prmtop_labels(d->atom_type, &top, "AMBER_ATOM_TYPE", N, path)) return false;
    }
    if (prmtop_block(&top, "RESIDUE_NUMBER")) {
        d->res_number = md_alloc(alloc, R * sizeof(int64_t));
        if (!prmtop_ints(d->res_number, &top, "RESIDUE_NUMBER", R, path)) d->res_number = NULL;
    }
    if (prmtop_block(&top, "RESIDUE_CHAINID")) {
        d->res_chain = md_alloc(alloc, R * sizeof(md_str_t));
        if (!prmtop_labels(d->res_chain, &top, "RESIDUE_CHAINID", R, path)) d->res_chain = NULL;
    }
    if (d->ifbox > 0 && prmtop_block(&top, "BOX_DIMENSIONS")) {
        d->has_box = prmtop_reals(d->box, &top, "BOX_DIMENSIONS", 4, path) && d->box[1] > 0 && d->box[2] > 0 && d->box[3] > 0;
    }
    {
        const prmtop_block_t* t = prmtop_block(&top, top.chamber ? "CTITLE" : "TITLE");
        if (!t) t = prmtop_block(&top, "TITLE");
        if (t) d->title = md_str_trim(t->text);
    }

    // ## Bonds: triples of (3 * (i - 1), 3 * (j - 1), parameter index), the coordinate offsets of the atoms
    const char* bond_blocks[] = { "BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN" };
    for (size_t k = 0; k < MD_COUNTOF(bond_blocks); ++k) {
        const prmtop_block_t* b = prmtop_block(&top, bond_blocks[k]);
        if (!b) continue;
        const size_t n = prmtop_count(b);
        if (n % 3) {
            MD_LOG_ERROR("Amber topology: %s in '" MD_STR_FMT "' is not made of triples", bond_blocks[k], MD_STR_ARG(path));
            return false;
        }
        if (!n) continue;
        int64_t* v = md_alloc(alloc, n * sizeof(int64_t));
        if (!prmtop_ints(v, &top, bond_blocks[k], n, path)) return false;
        for (size_t i = 0; i < n; i += 3) {
            const int64_t a = v[i], b2 = v[i + 1];
            if (a < 0 || b2 < 0 || a % 3 || b2 % 3 || a / 3 >= natom || b2 / 3 >= natom || a == b2) {
                MD_LOG_ERROR("Amber topology: %s in '" MD_STR_FMT "' has a bond between atoms that are not there", bond_blocks[k], MD_STR_ARG(path));
                return false;
            }
            const md_atom_pair_t pair = {{ (md_atom_idx_t)(a / 3), (md_atom_idx_t)(b2 / 3) }};
            md_array_push(d->bonds, pair, alloc);
        }
    }
    d->num_bonds = md_array_size(d->bonds);
    return true;
}

bool md_amber_prmtop_file_is_topology(md_str_t filename) {
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) return false;
    char buf[256];
    const size_t len = md_file_read(file, buf, sizeof(buf));
    md_file_close(&file);
    const md_str_t head = md_str_trim_beg((md_str_t){buf, len});
    return md_str_begins_with(head, MD_STR_LIT("%VERSION")) || md_str_begins_with(head, MD_STR_LIT("%FLAG"));
}

size_t md_amber_prmtop_file_atom_count(md_str_t filename) {
    // POINTERS comes after TITLE, near the top: the first block of the file is enough
    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) return 0;
    const size_t cap = 64 * 1024;
    char* buf = md_alloc(md_get_heap_allocator(), cap + 1);
    const size_t len = md_file_read(file, buf, cap);
    md_file_close(&file);
    buf[len] = '\0';

    size_t natom = 0;
    md_str_t text = {buf, len};
    bool in_pointers = false;
    while (text.len) {
        const md_str_t line = prmtop_next_line(&text);
        if (md_str_begins_with(line, MD_STR_LIT("%FLAG"))) {
            in_pointers = md_str_eq_cstr(md_str_trim(md_str_substr(line, 5, SIZE_MAX)), "POINTERS");
        } else if (in_pointers && line.len && line.ptr[0] != '%') {
            int64_t v = 0;
            const md_str_t field = md_str_trim(md_str_substr(line, 0, 8));
            if (md_parse_i64(&v, field) == field.len && v > 0) natom = (size_t)v;
            break;
        }
    }
    md_free(md_get_heap_allocator(), buf, cap + 1);
    return natom;
}

// The unit cell of BOX_DIMENSIONS: beta and the three lengths. A truncated octahedron (IFBOX 2) has all three angles
// beta; otherwise beta is the only angle that is not 90.
static md_unitcell_t prmtop_cell(const prmtop_data_t* d) {
    if (!d->has_box) return md_unitcell_none();
    const double beta = d->box[0];
    if (d->ifbox == 2 || fabs(beta - 109.4712206) < 0.01) {
        return md_unitcell_from_extent_and_angles(d->box[1], d->box[2], d->box[3], beta, beta, beta);
    }
    return md_unitcell_from_extent_and_angles(d->box[1], d->box[2], d->box[3], 90.0, beta > 0.0 && beta < 180.0 ? beta : 90.0, 90.0);
}

static md_topology_desc_t prmtop_desc(const prmtop_data_t* d) {
    return (md_topology_desc_t){
        .format        = "Amber topology",
        .num_atoms     = d->num_atoms,
        .atom_name     = d->atom_name,
        .atom_type     = d->atom_type,
        .charge        = d->charge,
        .mass          = d->mass,
        .atomic_number = d->atomic_number,
        .num_residues  = d->num_residues,
        .res_name      = d->res_name,
        .res_offset    = d->res_offset,
        .res_number    = d->res_number,
        .res_chain     = d->res_chain,
        .num_bonds     = d->num_bonds,
        .bonds         = d->bonds,
        .cell          = prmtop_cell(d),
        .title         = d->title,
    };
}

bool md_amber_prmtop_system_init_from_file(md_system_t* sys, md_system_state_t* state, md_str_t filename) {
    MD_ASSERT(sys);
    MD_ASSERT(state);
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    prmtop_data_t d;
    bool result = prmtop_read(&d, filename, md_temp_allocator(temp));
    if (result) {
        const md_topology_desc_t desc = prmtop_desc(&d);
        result = md_topology_system_init(sys, state, &desc);
        if (result && d.chamber) {
            MD_LOG_INFO("Amber topology: '" MD_STR_FMT "' is a CHAMBER topology (a CHARMM force field); its CHARMM terms are not read", MD_STR_ARG(filename));
        }
    }
    md_temp_end(temp);
    return result;
}

bool md_amber_prmtop_system_supplement_from_file(md_system_t* sys, md_str_t filename) {
    MD_ASSERT(sys);
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    prmtop_data_t d;
    bool result = prmtop_read(&d, filename, md_temp_allocator(temp));
    if (result) {
        const md_topology_desc_t desc = prmtop_desc(&d);
        result = md_topology_system_supplement(sys, &desc, filename);
    }
    md_temp_end(temp);
    return result;
}
