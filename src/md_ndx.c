#include <md_ndx.h>
#include <md_system.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_bitset.h>
#include <core/md_log.h>
#include <core/md_parse.h>
#include <core/md_reader.h>
#include <core/md_str.h>
#include <core/md_file.h>

static bool ndx_parse(md_ndx_data_t* data, md_reader_i* r, md_allocator_i* alloc) {
    md_memset(data, 0, sizeof(md_ndx_data_t));
    md_array(md_ndx_group_t) groups = 0;
    md_array(uint32_t) atoms = 0;
    size_t line_num = 0;
    bool ok = true;

    md_str_t line;
    while (ok && md_reader_extract_line(r, &line)) {
        ++line_num;
        line = md_str_substr(line, 0, md_str_find_char(line, ';'));
        line = md_str_trim(line);
        if (md_str_empty(line)) continue;

        if (line.ptr[0] == '[') {
            size_t close;
            if ((close = md_str_find_char(line, ']')) == MD_STR_NPOS) {
                MD_LOG_ERROR("NDX: Line %zu: a group header without ']'", line_num);
                ok = false;
                break;
            }
            // The atoms of the group before go with it
            if (md_array_size(groups)) {
                md_ndx_group_t* g = md_array_last(groups);
                g->count = md_array_size(atoms);
                if (g->count) md_array_push_array(g->atom_idx, atoms, g->count, alloc);
                md_array_shrink(atoms, 0);
            }
            const md_str_t name = md_str_trim(md_str_substr(line, 1, close - 1));
            md_ndx_group_t g = { .name = md_str_empty(name) ? (md_str_t){0} : md_alloc_str(alloc, name) };
            md_array_push(groups, g, alloc);
            continue;
        }

        if (!md_array_size(groups)) {
            MD_LOG_ERROR("NDX: Line %zu: atom numbers before the first group header", line_num);
            ok = false;
            break;
        }
        md_str_t tok;
        while (md_parse_next_token(&tok, &line)) {
            int64_t v = 0;
            if (md_parse_i64(&v, tok) != tok.len || v < 1 || v > (int64_t)MD_BITSET_MAX_INDEX) {
                MD_LOG_ERROR("NDX: Line %zu: '" MD_STR_FMT "' is not an atom number", line_num, MD_STR_ARG(tok));
                ok = false;
                break;
            }
            const uint32_t idx = (uint32_t)(v - 1);
            md_array_push(atoms, idx, md_get_heap_allocator());
            if (idx + 1 > data->max_atom) data->max_atom = idx + 1;
        }
    }

    if (ok && md_array_size(groups)) {
        md_ndx_group_t* g = md_array_last(groups);
        g->count = md_array_size(atoms);
        if (g->count) md_array_push_array(g->atom_idx, atoms, g->count, alloc);
    }
    md_array_free(atoms, md_get_heap_allocator());

    data->groups = groups;
    data->num_groups = md_array_size(groups);
    if (!ok) {
        md_ndx_data_free(data, alloc);
        return false;
    }
    if (!data->num_groups) {
        MD_LOG_ERROR("NDX: No groups");
        return false;
    }
    return true;
}

bool md_ndx_data_parse_str(md_ndx_data_t* data, md_str_t str, md_allocator_i* alloc) {
    MD_ASSERT(data);
    MD_ASSERT(alloc);
    md_reader_i r = md_reader_from_str(str);
    return ndx_parse(data, &r, alloc);
}

bool md_ndx_data_parse_file(md_ndx_data_t* data, md_str_t filename, md_allocator_i* alloc) {
    MD_ASSERT(data);
    MD_ASSERT(alloc);
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_file_reader_t fr;
    bool result = false;
    if (md_file_reader_open(&fr, filename, MD_KILOBYTES(64), md_temp_allocator(temp))) {
        result = ndx_parse(data, &fr.r, alloc);
        md_file_reader_close(&fr);
    } else {
        MD_LOG_ERROR("NDX: Could not open '" MD_STR_FMT "'", MD_STR_ARG(filename));
    }
    md_temp_end(temp);
    return result;
}

void md_ndx_data_free(md_ndx_data_t* data, md_allocator_i* alloc) {
    MD_ASSERT(data);
    MD_ASSERT(alloc);
    for (size_t i = 0; i < md_array_size(data->groups); ++i) {
        if (!md_str_empty(data->groups[i].name)) md_free_str(alloc, data->groups[i].name);
        md_array_free(data->groups[i].atom_idx, alloc);
    }
    md_array_free(data->groups, alloc);
    md_memset(data, 0, sizeof(md_ndx_data_t));
}

bool md_ndx_system_supplement(md_system_t* sys, const md_ndx_data_t* data) {
    MD_ASSERT(sys);
    MD_ASSERT(data);
    if (data->max_atom > sys->atom.count) {
        MD_LOG_ERROR("NDX: The index file names atom %u, the system has %zu atoms: it is of another system", data->max_atom, sys->atom.count);
        return false;
    }
    md_allocator_i* heap = md_get_heap_allocator();
    for (size_t i = 0; i < data->num_groups; ++i) {
        const md_ndx_group_t* g = &data->groups[i];
        md_bitset_t atoms = md_bitset_from_indices(g->atom_idx, g->count, heap);
        md_system_group_add(sys, g->name, atoms, 0);
        md_bitset_free(&atoms, heap);
    }
    MD_LOG_INFO("NDX: %zu groups", data->num_groups);
    return true;
}

bool md_ndx_system_supplement_from_file(md_system_t* sys, md_str_t filename) {
    MD_ASSERT(sys);
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_ndx_data_t data = {0};
    bool result = md_ndx_data_parse_file(&data, filename, arena) && md_ndx_system_supplement(sys, &data);
    md_arena_allocator_destroy(arena);
    return result;
}

// ### GROMACS' DEFAULT GROUPS (gmx analyse(), src/gromacs/topology/index.cpp) ###

enum {
    GMX_RESTYPE_OTHER = 0,
    GMX_RESTYPE_PROTEIN,
    GMX_RESTYPE_DNA,
    GMX_RESTYPE_RNA,
    GMX_RESTYPE_ION,
    GMX_RESTYPE_WATER,
};

#include "md_gmx_restypes.inl"

static const char* const gmx_restype_names[] = { "Other", "Protein", "DNA", "RNA", "Ion", "Water" };

static uint8_t gmx_restype(md_str_t resname) {
    for (size_t i = 0; i < MD_COUNTOF(gmx_restypes); ++i) {
        if (md_str_eq_cstr_ignore_case(resname, gmx_restypes[i].name)) return gmx_restypes[i].type;
    }
    return GMX_RESTYPE_OTHER;
}

// What GROMACS makes the groups of protein by: names of atoms, of which the first 'whole' are compared whole and the
// rest as a prefix (-1: all whole), and whether the group is the atoms that match or the others
typedef struct gmx_prot_group_t {
    const char* name;
    const char* const* atoms;
    int num_atoms;
    bool complement;
    int whole;
} gmx_prot_group_t;

static const char* const gmx_pnoh[]   = { "H", "HN" };
static const char* const gmx_pnodum[] = { "MN1", "MN2", "MCB1", "MCB2", "MCG1", "MCG2", "MCD1", "MCD2", "MCE1", "MCE2", "MNZ1", "MNZ2" };
static const char* const gmx_calpha[] = { "CA" };
static const char* const gmx_bb[]     = { "N", "CA", "C" };
static const char* const gmx_mc[]     = { "N", "CA", "C", "O", "O1", "O2", "OC1", "OC2", "OT", "OXT" };
static const char* const gmx_mcb[]    = { "N", "CA", "CB", "C", "O", "O1", "O2", "OC1", "OC2", "OT", "OXT" };
static const char* const gmx_mch[]    = { "N", "CA", "C", "O", "O1", "O2", "OC1", "OC2", "OT", "OXT", "H1", "H2", "H3", "H", "HN" };

static const gmx_prot_group_t gmx_prot_groups[] = {
    { "Protein",      NULL,       0,                         true,  -1 },
    { "Protein-H",    gmx_pnoh,   MD_COUNTOF(gmx_pnoh),      true,   0 },
    { "C-alpha",      gmx_calpha, MD_COUNTOF(gmx_calpha),    false, -1 },
    { "Backbone",     gmx_bb,     MD_COUNTOF(gmx_bb),        false, -1 },
    { "MainChain",    gmx_mc,     MD_COUNTOF(gmx_mc),        false, -1 },
    { "MainChain+Cb", gmx_mcb,    MD_COUNTOF(gmx_mcb),       false, -1 },
    { "MainChain+H",  gmx_mch,    MD_COUNTOF(gmx_mch),       false, -1 },
    { "SideChain",    gmx_mch,    MD_COUNTOF(gmx_mch),       true,  -1 },
    { "SideChain-H",  gmx_mch,    MD_COUNTOF(gmx_mch),       true,  11 },
    { "Prot-Masses",  gmx_pnodum, MD_COUNTOF(gmx_pnodum),    true,   0 },
};

static bool gmx_prot_atom_in(const gmx_prot_group_t* g, md_str_t atom_name) {
    // Digits at the beginning of a name are skipped (1H)
    while (atom_name.len && atom_name.ptr[0] >= '0' && atom_name.ptr[0] <= '9') { atom_name.ptr++; atom_name.len--; }
    bool match = false;
    for (int j = 0; j < g->num_atoms && !match; ++j) {
        const md_str_t def = md_str_from_cstr(g->atoms[j]);
        if (g->whole == -1 || j < g->whole) {
            match = md_str_eq_ignore_case(def, atom_name);
        } else {
            match = atom_name.len >= def.len && md_str_eq_ignore_case(def, md_str_substr(atom_name, 0, def.len));
        }
    }
    return g->complement != match;
}

typedef struct gmx_default_ctx_t {
    const md_system_t* sys;
    const uint8_t* type;        // [residue] gmx restype
    md_str_t want;
    md_bitmask_t mask;
    bool found;
} gmx_default_ctx_t;

// The residues of a type (or of every type but it), as atoms
static void gmx_mask_type(gmx_default_ctx_t* c, uint8_t type, bool other_than) {
    for (size_t r = 0; r < md_system_residue_count(c->sys); ++r) {
        if ((c->type[r] == type) == other_than) continue;
        const md_urange_t ar = md_system_residue_atom_range(c->sys, r);
        md_bitmask_set_range(&c->mask, ar.beg, ar.end);
    }
}

static bool gmx_is(gmx_default_ctx_t* c, const char* name) {
    if (c->found || !md_str_eq_cstr_ignore_case(c->want, name)) return false;
    c->found = true;
    return true;
}

bool md_ndx_default_group(md_bitset_t* out, const md_system_t* sys, md_str_t name, md_allocator_i* alloc) {
    MD_ASSERT(out);
    MD_ASSERT(sys);
    MD_ASSERT(alloc);
    if (md_str_empty(name) || !sys->atom.count) return false;
    md_allocator_i* heap = md_get_heap_allocator();
    const size_t R = md_system_residue_count(sys);
    const uint32_t N = (uint32_t)sys->atom.count;
    gmx_default_ctx_t c = { .sys = sys, .want = name };
    md_bitmask_init(&c.mask, 0, N, heap);
    uint8_t* type = R ? md_alloc(heap, R) : NULL;
    // By the residue types where the system has them: a few names to look up rather than one per residue
    const size_t T = sys->residue.type.count;
    uint8_t* of_type = (T && sys->residue.type_idx) ? md_alloc(heap, T) : NULL;
    for (size_t t = 0; of_type && t < T; ++t) of_type[t] = 0xFF;
    for (size_t r = 0; r < R; ++r) {
        const size_t t = of_type ? (size_t)sys->residue.type_idx[r] : SIZE_MAX;
        if (t < T) {
            if (of_type[t] == 0xFF) of_type[t] = gmx_restype(md_system_residue_name(sys, r));
            type[r] = of_type[t];
        } else {
            type[r] = gmx_restype(md_system_residue_name(sys, r));
        }
    }
    if (of_type) md_free(heap, of_type, T);
    c.type = type;

    if (gmx_is(&c, "System")) {
        md_bitmask_set_range(&c.mask, 0, N);
    }
    // The residue types in the order they first appear, each with what GROMACS makes of it
    bool seen[6] = {0};
    bool analysed_other = false, have_water = false, have_ion = false;
    for (size_t r = 0; r < R && !c.found; ++r) {
        const uint8_t t = type[r];
        if (seen[t]) continue;
        seen[t] = true;
        if (t == GMX_RESTYPE_PROTEIN) {
            for (size_t g = 0; g < MD_COUNTOF(gmx_prot_groups) && !c.found; ++g) {
                if (!gmx_is(&c, gmx_prot_groups[g].name)) continue;
                for (size_t q = 0; q < R; ++q) {
                    if (type[q] != GMX_RESTYPE_PROTEIN) continue;
                    const md_urange_t ar = md_system_residue_atom_range(sys, q);
                    for (uint32_t a = ar.beg; a < ar.end; ++a) {
                        if (gmx_prot_atom_in(&gmx_prot_groups[g], md_atom_name(&sys->atom, a))) md_bitmask_set_bit(&c.mask, a);
                    }
                }
            }
            if (gmx_is(&c, "non-Protein")) gmx_mask_type(&c, GMX_RESTYPE_PROTEIN, true);
        } else if (t == GMX_RESTYPE_WATER) {
            have_water = true;
            if (gmx_is(&c, "Water") || gmx_is(&c, "SOL")) gmx_mask_type(&c, GMX_RESTYPE_WATER, false);
            if (gmx_is(&c, "non-Water")) gmx_mask_type(&c, GMX_RESTYPE_WATER, true);
        } else if (t == GMX_RESTYPE_ION) {
            have_ion = true;
            if (gmx_is(&c, "Ion")) gmx_mask_type(&c, GMX_RESTYPE_ION, false);
        } else if (!analysed_other) {
            // The first of DNA, RNA and Other only, as GROMACS has it; then a group per residue name of all that is
            // not protein, DNA, RNA or water, ions included
            analysed_other = true;
            if (gmx_is(&c, gmx_restype_names[t])) gmx_mask_type(&c, t, false);
            // GROMACS makes them in the order the names first appear, and grompp takes the first that matches whatever
            // the case: so the first residue whose name matches names it (in its own case), and every residue of that
            // name is in it. One pass, however many such residues (a Martini system's are all of them).
            md_str_t rn = {0};
            for (size_t q = 0; q < R && (!c.found || !md_str_empty(rn)); ++q) {
                const uint8_t tq = type[q];
                if (tq == GMX_RESTYPE_PROTEIN || tq == GMX_RESTYPE_DNA || tq == GMX_RESTYPE_RNA || tq == GMX_RESTYPE_WATER) continue;
                const md_str_t name_q = md_system_residue_name(sys, q);
                if (md_str_empty(rn)) {
                    if (!md_str_eq_ignore_case(c.want, name_q)) continue;
                    rn = name_q;
                    c.found = true;
                } else if (!md_str_eq(name_q, rn)) {
                    continue;
                }
                const md_urange_t ar = md_system_residue_atom_range(sys, q);
                md_bitmask_set_range(&c.mask, ar.beg, ar.end);
            }
        }
    }
    if (!c.found && R) {
        for (size_t r = 0; r < R; ++r) { have_water |= type[r] == GMX_RESTYPE_WATER; have_ion |= type[r] == GMX_RESTYPE_ION; }
        if (have_water && have_ion && gmx_is(&c, "Water_and_ions")) {
            gmx_mask_type(&c, GMX_RESTYPE_WATER, false);
            gmx_mask_type(&c, GMX_RESTYPE_ION, false);
        }
    }

    // As GROMACS makes them only where they are some atoms and not all
    if (c.found && (md_str_eq_cstr_ignore_case(name, "non-Protein") || md_str_eq_cstr_ignore_case(name, "non-Water"))) {
        const uint64_t n = md_bitmask_count(&c.mask);
        c.found = n > 0 && n < N;
    }
    if (c.found) *out = md_bitmask_extract(&c.mask, alloc);
    md_bitmask_free(&c.mask);
    if (type) md_free(heap, type, R);
    return c.found;
}
