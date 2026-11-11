#include "md_infer_internal.h"

// Particle type inference: the types of the atoms and beads from their names and residues (md_util.h).

// The type of the name, element, mass, size and particle kind, added when there is none. md_particle_type_find_or_add goes by
// name and element alone, which gives every bead of a name the mass and size of the first: a Martini SC1 is 36 in Ala
// and 54 in Leu (Martini 3), and Trp's SC3 is a massless virtual site.
md_particle_type_idx_t mdi_find_or_add_particle_type(md_particle_type_data_t* type, md_str_t name, md_atomic_number_t z, float mass, float radius, uint32_t color, md_particle_type_flags_t flags, md_allocator_i* alloc) {
    for (size_t t = 1; t < type->count; ++t) {
        if (type->z[t] == z && type->mass[t] == mass && type->radius[t] == radius && md_particle_type_flags(type, t) == flags && md_str_eq(md_particle_type_name(type, t), name)) {
            return (md_particle_type_idx_t)t;
        }
    }
    return md_particle_type_add(type, name, (md_str_t){0}, z, mass, radius, color, flags, alloc);
}

md_particle_type_idx_t md_util_particle_type_find_or_add(md_system_t* sys, md_hashmap32_t* map, md_str_t name, md_atomic_number_t z) {
    MD_ASSERT(sys && map);
    // The name as the type keeps it, its first 8 characters
    const uint64_t kv[2] = { md_label_key(md_label_make(name)), (uint64_t)z };
    const uint64_t key = md_hash64(kv, sizeof(kv), 0);
    const uint32_t* found = md_hashmap_get(map, key);
    if (found) return (md_particle_type_idx_t)*found;
    const md_particle_type_idx_t type = md_particle_type_find_or_add(&sys->atom.type, name, z, md_atomic_number_mass(z), md_atomic_number_vdw_radius(z),
        md_atomic_number_cpk_color(z), 0, sys->alloc);
    md_hashmap_add(map, key, (uint32_t)type);
    return type;
}

typedef struct {
    const char* res;                // NULL for any component named as its res_kind: an amino acid or a nucleotide
    const char* atom;
    int   atomic_nr;
    float mass;
    float radius;
    md_particle_kind_t  particle;   // What the particle is
    md_atom_flags_t     role;       // Its role in its residue (backbone, side chain)
    md_residue_kind_t res_kind;     // What its residue is, when the particle tells (OTHER when it does not)
} particle_type_t;


// Predefined particle types (This include coarse grained types)
static const particle_type_t predefined_particle_types[] = {
    // Martini proteins: the backbone bead and the side chain beads of any residue named as an amino acid, which covers
    // the names martinize2 gives protonation states and lipidated residues. No atomistic force field names an atom BB
    // or SC#. The mass and radius are those of a regular bead; the residues of martini_db have those of each version.
    // Martini molecules other than amino acids and nucleotides (lipids, water, ions, ...) are told by their residue
    // name and all their bead names at once (mdi_martini_db_match_molecule), not here.
    {NULL, "BB",      0, 72.0f, MARTINI_R, MD_PARTICLE_BEAD, MD_ATOM_FLAG_BACKBONE,   MD_RESIDUE_KIND_AMINO_ACID},
    {NULL, "SC#",     0, 72.0f, MARTINI_R, MD_PARTICLE_BEAD, MD_ATOM_FLAG_SIDE_CHAIN, MD_RESIDUE_KIND_AMINO_ACID},
    // The charge dummies of the polar side chains of polarizable Martini 2 (martini22p, elnedyn22p): no Lennard-Jones,
    // so the radius of element 0, as md_tpr gives a particle without
    {NULL, "SCP|SCN", 0, 36.0f, 1.0f,      MD_PARTICLE_BEAD, MD_ATOM_FLAG_SIDE_CHAIN, MD_RESIDUE_KIND_AMINO_ACID},
    // The lipid tails of the lipidated residues of Martini 3: farnesyl, geranylgeranyl and palmitoyl Cys, myristoyl Gly
    {"CYSF|CYSG|CYSP|GLYM", "C#", 0, 72.0f, MARTINI_R, MD_PARTICLE_BEAD, MD_ATOM_FLAG_SIDE_CHAIN, MD_RESIDUE_KIND_AMINO_ACID},

    // Martini nucleotides (Martini 2 DNA, Uusitalo et al. 2015): the phosphate (BB1) and the sugar (BB2, BB3) of the
    // backbone, and the base (SC1-SC4), of any residue named as a nucleotide
    {NULL, "BB1",     0, 72.0f, MARTINI_R, MD_PARTICLE_BEAD, MD_ATOM_FLAG_BACKBONE, MD_RESIDUE_KIND_NUCLEOTIDE},
    {NULL, "BB#",     0, 72.0f, MARTINI_R, MD_PARTICLE_BEAD, MD_ATOM_FLAG_BACKBONE | MD_ATOM_FLAG_NUCLEOSIDE, MD_RESIDUE_KIND_NUCLEOTIDE},
    {NULL, "SC#",     0, 72.0f, MARTINI_R, MD_PARTICLE_BEAD, MD_ATOM_FLAG_NUCLEOBASE | MD_ATOM_FLAG_NUCLEOSIDE, MD_RESIDUE_KIND_NUCLEOTIDE},

    {"RAMP", "PO1", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GM1", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GM2", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GM3", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GM4", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GM5", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GM6", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "PO2", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GL1", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GL2", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C1A", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C2A", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C3A", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C1B", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C2B", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C3B", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GL3", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GL4", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C1C", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C2C", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C3C", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C1D", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C2D", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C3D", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GL5", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GL6", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C1E", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C2E", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GL7", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "GL8", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C1F", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "C2F", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S01", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S02", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S03", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S04", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S05", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S06", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S07", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S08", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S09", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S10", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S11", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S12", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S13", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S14", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S15", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S16", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S17", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S18", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S19", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S20", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S21", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S22", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S23", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S24", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S25", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S26", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S27", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S28", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S29", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S30", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S31", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S32", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S33", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S34", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S35", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S36", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S37", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S38", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"RAMP", "S39", 0, 50.f, 4.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},

    // The massless charge sites of 4 and 5 site water models: M of TIP4P and its variants and OPC (MW in GROMACS, OM
    // in CHARMM, EPW in AMBER, M in OpenMM) and the two lone pairs of TIP5P (LP1 LP2, EP1 EP2, M1 M2). These are all
    // atom models of a single molecule, NOT coarse grained: the coarse grained flag would turn off covalent bond and
    // hydrogen bond inference for the whole system. No element, no mass, no radius, and no bonds.
    { "SOL|WAT|HOH|H2O|TIP4|TIP5|TIP4P|TIP5P|TP4|TP5|T4P|T5P|T4E|OPC", "MW|M|OM|EPW|EP|EP#|LP#|M#", 0, 0.0f, 0.0f, MD_PARTICLE_VIRTUAL_SITE, 0, MD_RESIDUE_KIND_OTHER },

    { "ION", "CL",  17, 35.45f,     1.8f, MD_PARTICLE_ATOM, 0, MD_RESIDUE_KIND_ION },
    { "ION", "NA",  11, 22.99f,     2.3f, MD_PARTICLE_ATOM, 0, MD_RESIDUE_KIND_ION },
    { "ION", "K",   19, 39.10f,     2.7f, MD_PARTICLE_ATOM, 0, MD_RESIDUE_KIND_ION },
    { "ION", "CA",  20, 40.08f,     2.7f, MD_PARTICLE_ATOM, 0, MD_RESIDUE_KIND_ION },
	{ "ION", "MG",  12, 24.31f,     2.0f, MD_PARTICLE_ATOM, 0, MD_RESIDUE_KIND_ION },

    { "CLA", "CLA", 17, 35.45f,     1.8f, MD_PARTICLE_ATOM, 0, MD_RESIDUE_KIND_ION },

    {"*", "IC", 0,   3240.0f,   12.5f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"*", "OC", 0,   3240.0f,   12.5f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"*", "CC", 0,   3900.0f,   14.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},

    {"P10", "X1", 0, 3.31E6f,   50.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"P20", "X2", 0, 2.65E7f,   100.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"P30", "X3", 0, 8.93E7f,   150.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
    {"P40", "X4", 0, 2.12E8f,   200.0f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},

    {"C##", "C#?", 0, 1727.7f,   13.5f, MD_PARTICLE_BEAD, 0, MD_RESIDUE_KIND_OTHER},
};

static inline uint64_t gen_key_from_names(md_str_t res_name, md_str_t atom_name) {
    char buf[64];
    int len = md_snprintf(buf, sizeof(buf), MD_STR_FMT "_" MD_STR_FMT, MD_STR_ARG(res_name), MD_STR_ARG(atom_name));
    md_str_t str = {buf, (size_t)len};
    return md_hash64_str(str, 0);
}

// Searches the predefined particle types for a matching atom name within a residue
// We support simple pattern matching here like wildcards and numerical substitutes
// * is supported to match any sequence of characters
// ? is supported to match any single character
// # is supported to match any single digit
// | separates alternatives

static inline bool pattern_match(const char* pattern, const char* str) {
    while (*pattern && *str) {
        if (*pattern == '*') {
            pattern++;
            if (!*pattern) return true;
            while (*str) {
                if (pattern_match(pattern, str)) return true;
                str++;
            }
            return false;
        } else if (*pattern == '?') {
            pattern++;
            str++;
        } else if (*pattern == '#') {
            if (!md_ascii_is_digit((unsigned char)*str)) return false;
            pattern++;
            str++;
        } else {
            if (*pattern != *str) return false;
            pattern++;
            str++;
        }
    }
    while (*pattern == '*') pattern++;
    return !*pattern && !*str;
}

// Alternatives separated by '|': "SOL|WAT" matches either
static bool pattern_match_any(const char* patterns, const char* str) {
    char buf[128];
    const char* beg = patterns;
    for (;;) {
        const char* end = beg;
        while (*end && *end != '|') ++end;
        const size_t len = (size_t)(end - beg);
        if (len < sizeof(buf)) {
            md_memcpy(buf, beg, len);
            buf[len] = '\0';
            if (pattern_match(buf, str)) return true;
        }
        if (!*end) return false;
        beg = end + 1;
    }
}

// Whether a residue name is that of the kind: what a predefined type without a residue name matches
static bool resname_of_kind(md_str_t res_name, md_residue_kind_t kind) {
    switch (kind) {
    case MD_RESIDUE_KIND_AMINO_ACID: return md_util_resname_amino_acid(res_name);
    case MD_RESIDUE_KIND_NUCLEOTIDE: return md_util_resname_nucleotide(res_name);
    default: return false;
    }
}

static particle_type_t* find_predefined_particle_type(md_str_t res_name, md_str_t atom_name) {
    // The patterns are matched against C strings, and the names need not end in a zero (a label of 8 characters, a
    // name read from a file)
    char res_buf[32], atom_buf[32];
    md_str_copy_to_char_buf(res_buf, sizeof(res_buf), res_name);
    md_str_copy_to_char_buf(atom_buf, sizeof(atom_buf), atom_name);
    for (size_t i = 0; i < MD_COUNTOF(predefined_particle_types); ++i) {
        particle_type_t type = predefined_particle_types[i];
        if (!pattern_match_any(type.atom, atom_buf)) continue;
        const bool res_match = type.res ? pattern_match_any(type.res, res_buf) : resname_of_kind(res_name, type.res_kind);
        if (res_match) {
            return (particle_type_t*)&predefined_particle_types[i];
        }
    }
    return NULL;
}

// The residue kind a predefined type tells, applied when the residue has none yet
static inline void apply_predefined_res_kind(md_system_t* sys, size_t res_idx, md_residue_kind_t kind) {
    if (kind != MD_RESIDUE_KIND_OTHER && md_system_residue_kind(sys, res_idx) == MD_RESIDUE_KIND_OTHER) {
        md_system_set_residue_kind(sys, res_idx, kind);
    }
}

typedef struct residue_template_t {
    size_t            res_idx;          // The first residue of the key: its atoms' types are the template's
    size_t            role_offset;      // In the roles: the atom flags of the residue before it was typed
    md_residue_kind_t kind;             // What the typing said the residue is
    int32_t           unknown_variant;  // The Martini molecule of no topology it was, -1 if none
} residue_template_t;

// The roles the typing gave the template residue's atoms: their flags now, without what they had before
static void add_residue_template(md_hashmap32_t* cache, md_array(residue_template_t)* templates, md_atom_flags_t* roles, residue_template_t* tmpl, uint64_t key, const md_system_t* sys, md_allocator_i* alloc) {
    const md_urange_t range = md_system_residue_atom_range(sys, tmpl->res_idx);
    for (uint32_t i = range.beg; i < range.end; ++i) {
        md_atom_flags_t* role = &roles[tmpl->role_offset + (i - range.beg)];
        const md_atom_flags_t now = sys->atom.flags ? sys->atom.flags[i] : 0;
        *role = (md_atom_flags_t)((uint32_t)now & ~(uint32_t)*role);
    }
    md_hashmap_add(cache, key, (uint32_t)md_array_size(*templates));
    md_array_push(*templates, *tmpl, alloc);
}

// The key of each residue and whether its atoms have no types yet (untyped), a chunk of INFER_KEY_CHUNK at a time: what
// the typing of a residue looks up the typing of the first of its key by. Of each residue alone, so in parallel.
#define INFER_KEY_CHUNK 4096

typedef struct infer_keys_t {
    const md_system_t* sys;
    md_str_t (*label)(const void*, size_t);
    const void* label_data;
    const md_label_t* words;    // The labels as md_label_t, when they are (else NULL)
    uint64_t* key;
    uint8_t*  untyped;
} infer_keys_t;

static md_str_t label_of_labels(const void* data, size_t i) {
    return md_label_str(&((const md_label_t*)data)[i]);
}

static void infer_residue_keys(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    const infer_keys_t* t = (const infer_keys_t*)arg;
    const md_system_t* sys = t->sys;
    const size_t r_end = MD_MIN((size_t)end * INFER_KEY_CHUNK, sys->residue.count);
    for (size_t res_idx = (size_t)beg * INFER_KEY_CHUNK; res_idx < r_end; ++res_idx) {
        const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
        const size_t res_size = range.end - range.beg;
        // The labels are part of the key, as which of them a name stands for can depend on the others; and
        // everything the types of the residue's atoms come from is in it, so the atoms of a residue of a key seen
        // before (a water among thousands) take what the first one's did, without a look at each atom
        const md_str_t res_name = md_system_residue_name(sys, res_idx);
        uint64_t key = XXH3_64bits_withSeed(res_name.ptr, res_name.len, res_size);
        if (t->words) {
            // Labels as words (zero padded, of no zero byte): equal bytes are equal labels, so the residue's are hashed
            // in one go
            key = XXH3_64bits_withSeed(t->words + range.beg, res_size * sizeof(md_label_t), key);
        } else {
            for (size_t i = range.beg; i < range.end; ++i) {
                const md_str_t lbl = t->label(t->label_data, i);
                key = XXH3_64bits_withSeed(lbl.ptr, lbl.len, key);
            }
        }
        bool untyped = true;
        for (size_t i = range.beg; i < range.end && untyped; ++i) untyped = sys->atom.type_idx[i] == 0;
        t->key[res_idx] = key;
        t->untyped[res_idx] = untyped;
    }
}

// The residues of a key seen before: the types and roles of the first of it (its template), and what it is
typedef struct infer_apply_t {
    md_system_t* sys;
    const uint32_t* tmpl_of;    // [residues] The template of each, UINT32_MAX for none
    const residue_template_t* templates;
    const md_atom_flags_t* roles;
} infer_apply_t;

static void infer_apply_templates(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    const infer_apply_t* t = (const infer_apply_t*)arg;
    md_system_t* sys = t->sys;
    const size_t r_end = MD_MIN((size_t)end * INFER_KEY_CHUNK, sys->residue.count);
    for (size_t res_idx = (size_t)beg * INFER_KEY_CHUNK; res_idx < r_end; ++res_idx) {
        if (t->tmpl_of[res_idx] == UINT32_MAX) continue;
        const residue_template_t* tmpl = &t->templates[t->tmpl_of[res_idx]];
        const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
        const md_urange_t from = md_system_residue_atom_range(sys, tmpl->res_idx);
        const size_t res_size = range.end - range.beg;
        md_memcpy(sys->atom.type_idx + range.beg, sys->atom.type_idx + from.beg, res_size * sizeof(md_particle_type_idx_t));
        if (sys->atom.flags) {
            for (size_t k = 0; k < res_size; ++k) sys->atom.flags[range.beg + k] |= t->roles[tmpl->role_offset + k];
        }
        apply_predefined_res_kind(sys, res_idx, tmpl->kind);
    }
}

// The labels of the atoms through an accessor: an array of strings, or of labels. The residues are typed in order, each
// the first of its key: the others take what it got, over the threads of ex (NULL: the calling thread), as do the keys.
static void infer_particle_types(md_system_t* sys, md_str_t (*label)(const void*, size_t), const void* label_data, const md_executor_t* ex) {
    if (!sys) {
        MD_LOG_ERROR("system was null");
        return;
    }
    if (!sys->alloc) {
        MD_LOG_ERROR("System allocator was null");
        return;
    }
    md_allocator_i* alloc = sys->alloc;
    md_system_topology_changed(sys);

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    // Maps (residue name and size, atom label) to the type in the low 16 bits and the predefined type it came from
    // (its index + 1, 0 for none) in the high 16 bits: the predefined type also gives the atom its role.
    md_hashmap32_t particle_type_cache = {.allocator = temp_alloc };
    md_hashmap_reserve(&particle_type_cache, 256);
    MD_STATIC_ASSERT(sizeof(md_particle_type_idx_t) == 2, "The cache packs a type index in 16 bits");
    MD_STATIC_ASSERT(MD_COUNTOF(predefined_particle_types) < 0xFFFF, "The cache packs a predefined type in 16 bits");

    if (sys->residue.count > 0) {
        // The masses and sizes of Martini beads depend on the version of Martini
        const martini_system_t martini = mdi_martini_system_ex(sys, label, label_data, ex);
        // The molecules of the force field of beads of none of its topologies, by name and size: the first and how many
        typedef struct { md_str_t name; size_t size; size_t first; size_t count; } unknown_variant_t;
        md_array(unknown_variant_t) unknown_variants = 0;
        md_array(md_str_t) res_labels = 0;
        // The residues whose atoms were typed by a look at each: by their key, the types and roles their atoms got
        md_hashmap32_t residue_template_cache = {.allocator = temp_alloc };
        md_array(residue_template_t) templates = 0;
        md_array(md_atom_flags_t) template_roles = 0;

        const size_t num_res = sys->residue.count;
        const uint32_t num_chunks = (uint32_t)MD_DIV_UP(num_res, INFER_KEY_CHUNK);
        uint64_t* res_keys = md_temp_alloc(temp, num_res * sizeof(uint64_t));
        uint8_t*  res_untyped = md_temp_alloc(temp, num_res);
        uint32_t* tmpl_of = md_temp_alloc(temp, num_res * sizeof(uint32_t));
        {
            const md_label_t* words = label == label_of_labels ? (const md_label_t*)label_data : NULL;
            infer_keys_t task = { sys, label, label_data, words, res_keys, res_untyped };
            md_parallel_for(ex, num_chunks, 1, infer_residue_keys, &task);
        }
        // The template a residue of the key before was given: a run of one key (the waters, the lipids) looks it up once
        uint64_t hit_key = 0;
        uint32_t hit_tmpl = UINT32_MAX;
        size_t num_hits = 0;

        for (size_t res_idx = 0; res_idx < num_res; ++res_idx) {
            tmpl_of[res_idx] = UINT32_MAX;
            const uint64_t res_key = res_keys[res_idx];
            const bool untyped = res_untyped[res_idx] != 0;
            if (untyped) {
                uint32_t t = UINT32_MAX;
                if (hit_tmpl != UINT32_MAX && res_key == hit_key) {
                    t = hit_tmpl;
                } else {
                    const uint32_t* found = md_hashmap_get(&residue_template_cache, res_key);
                    if (found) {
                        t = *found;
                        hit_key = res_key;
                        hit_tmpl = t;
                    }
                }
                if (t != UINT32_MAX) {
                    // Given its template's types once all are typed
                    tmpl_of[res_idx] = t;
                    num_hits += 1;
                    if (templates[t].unknown_variant >= 0) unknown_variants[templates[t].unknown_variant].count += 1;
                    continue;
                }
            }

            md_str_t res_name         = md_system_residue_name(sys, res_idx);
            md_urange_t res_range = md_system_residue_atom_range(sys, res_idx);
            size_t res_size        = res_range.end - res_range.beg;

            // The labels of the residue, by its atoms from 0
            md_array_resize(res_labels, res_size, temp_alloc);
            for (size_t k = 0; k < res_size; ++k) res_labels[k] = label(label_data, res_range.beg + k);
            const md_str_t* labels = res_labels;

            // What the residue's atoms got, for those to follow of its key
            residue_template_t tmpl = { .res_idx = res_idx, .role_offset = md_array_size(template_roles), .kind = MD_RESIDUE_KIND_OTHER, .unknown_variant = -1 };
            if (untyped) {
                for (size_t i = res_range.beg; i < res_range.end; ++i) {
                    md_array_push(template_roles, sys->atom.flags ? sys->atom.flags[i] : 0, temp_alloc);
                }
            }

            // A Martini molecule (a lipid, water, an ion, ...): its residue name and all its bead names tell it, and give
            // every bead its mass, size and particle kind. A residue of a single bead is one only in a Martini system:
            // an ion of Martini 3 (ION NA) is named like an atomistic one.
            uint8_t map[256];
            const martini_db_residue_t* mol = mdi_martini_db_match_molecule(martini.version, res_name, labels, res_size, map);
            if (mol && res_size == 1 && !martini.martini && mol->kind != MD_RESIDUE_KIND_WATER) mol = NULL;
            // A molecule renamed (CHL for CHOL, the three letters of a PDB file) is told by its beads, in a Martini system
            // only: atomistic residues are named as their atoms, which may be the names of beads
            if (!mol && martini.martini && !md_util_resname_amino_acid(res_name) && !md_util_resname_nucleotide(res_name)) {
                mol = mdi_martini_db_match_beads(martini.version, labels, NULL, res_size, map);
            }
            if (mol) {
                for (size_t i = res_range.beg; i < res_range.end; ++i) {
                    if (sys->atom.type_idx[i] != 0) continue;
                    const martini_db_bead_t* bead = mdi_martini_db_bead(mol, map[i - res_range.beg]);
                    const md_particle_type_flags_t flags = md_particle_type_flags_set_particle_kind(MD_PARTICLE_TYPE_FLAG_NONE, bead->particle);
                    sys->atom.type_idx[i] = mdi_find_or_add_particle_type(&sys->atom.type, labels[i - res_range.beg], 0, bead->mass, bead->radius, 0, flags, alloc);
                }
                apply_predefined_res_kind(sys, res_idx, mol->kind);
                if (untyped) {
                    tmpl.kind = mol->kind;
                    add_residue_template(&residue_template_cache, &templates, template_roles, &tmpl, res_key, sys, temp_alloc);
                }
                continue;
            }
            // A molecule of the force field whose beads are those of none of its topologies (a variant of it, a residue
            // edited), in a Martini system: beads still, rather than the elements their names would read as (the dummy V1
            // of a sterol vanadium, its C1 carbon). Regular beads, as nothing tells their masses and sizes, of the kind of
            // the molecule, and without a topology of it no bonds. Not of a polymer, whose residues are bonded by name.
            if (martini.martini && !md_util_resname_amino_acid(res_name) && !md_util_resname_nucleotide(res_name) && !mdi_martini_db_linked(res_name)) {
                const martini_db_residue_t* like = mdi_martini_db_unknown_variant(res_name, labels, res_size);
                if (like) {
                    const md_particle_type_flags_t flags = md_particle_type_flags_set_particle_kind(MD_PARTICLE_TYPE_FLAG_NONE, MD_PARTICLE_BEAD);
                    for (size_t i = res_range.beg; i < res_range.end; ++i) {
                        if (sys->atom.type_idx[i] != 0) continue;
                        sys->atom.type_idx[i] = mdi_find_or_add_particle_type(&sys->atom.type, labels[i - res_range.beg], 0, 72.0f, MARTINI_R, 0, flags, alloc);
                    }
                    apply_predefined_res_kind(sys, res_idx, like->kind);
                    // Told once per name and size
                    size_t u = 0;
                    while (u < md_array_size(unknown_variants) && !(md_str_eq(unknown_variants[u].name, res_name) && unknown_variants[u].size == res_size)) ++u;
                    if (u == md_array_size(unknown_variants)) {
                        unknown_variant_t v = { res_name, res_size, res_idx, 0 };
                        md_array_push(unknown_variants, v, temp_alloc);
                    }
                    unknown_variants[u].count += 1;
                    if (untyped) {
                        tmpl.kind = like->kind;
                        tmpl.unknown_variant = (int32_t)u;
                        add_residue_template(&residue_template_cache, &templates, template_roles, &tmpl, res_key, sys, temp_alloc);
                    }
                    continue;
                }
            }

            // What the predefined types of its atoms say the residue is
            md_residue_kind_t res_kind = MD_RESIDUE_KIND_OTHER;
            for (size_t i = res_range.beg; i < res_range.end; ++i) {
                if (sys->atom.type_idx[i] != 0) continue;

                const particle_type_t* predef = NULL;
                uint64_t key = md_hash64_str(labels[i - res_range.beg], res_key);
                uint32_t* cached = md_hashmap_get(&particle_type_cache, key);
                if (cached) {
                    sys->atom.type_idx[i] = (md_particle_type_idx_t)(*cached & 0xFFFF);
                    const uint32_t entry = *cached >> 16;
                    predef = entry ? &predefined_particle_types[entry - 1] : NULL;
                } else {
                    md_str_t atom_name = labels[i - res_range.beg];
                    md_atomic_number_t z = 0;
                    float mass = 0;
                    float radius = 0;
                    uint32_t color = 0;
                    md_particle_type_flags_t flags = MD_PARTICLE_TYPE_FLAG_NONE;

                    // Try to find in predefined set
                    predef = find_predefined_particle_type(res_name, atom_name);
                    if (predef) {
                        z       = (md_atomic_number_t)predef->atomic_nr;
                        mass    = predef->mass;
                        radius  = predef->radius;
                        color   = z == 0 ? 0 : md_atomic_number_cpk_color(z);
                        flags   = md_particle_type_flags_set_particle_kind(flags, predef->particle);
                        if (predef->res_kind == MD_RESIDUE_KIND_AMINO_ACID || predef->res_kind == MD_RESIDUE_KIND_NUCLEOTIDE) {
                            // A Martini bead of an amino acid or a nucleotide: its version tells its mass and size
                            const martini_db_bead_t* bead = mdi_martini_db_biopolymer_bead(martini.version, res_name, labels, res_size, i - res_range.beg);
                            if (bead) {
                                mass   = bead->mass;
                                radius = bead->radius;
                                flags  = md_particle_type_flags_set_particle_kind(MD_PARTICLE_TYPE_FLAG_NONE, bead->particle);
                            }
                        }
                    } else {
                        z       = md_atomic_number_infer_from_label(atom_name, res_name, res_size);
                        mass    = md_atomic_number_mass(z);
                        radius  = md_atomic_number_vdw_radius(z);
                        color   = md_atomic_number_cpk_color(z);
                    }

                    md_particle_type_idx_t type = mdi_find_or_add_particle_type(&sys->atom.type, atom_name, z, mass, radius, color, flags, alloc);
                    sys->atom.type_idx[i] = type;
                    const uint32_t entry = predef ? (uint32_t)(predef - predefined_particle_types) + 1 : 0;
                    md_hashmap_add(&particle_type_cache, key, (uint32_t)type | (entry << 16));
                }
                if (predef) {
                    if (sys->atom.flags) sys->atom.flags[i] |= predef->role;
                    if (predef->res_kind != MD_RESIDUE_KIND_OTHER) res_kind = predef->res_kind;
                }
            }
            apply_predefined_res_kind(sys, res_idx, res_kind);
            if (untyped) {
                tmpl.kind = res_kind;
                add_residue_template(&residue_template_cache, &templates, template_roles, &tmpl, res_key, sys, temp_alloc);
            }
        }
        if (num_hits > 0) {
            infer_apply_t task = { sys, tmpl_of, templates, template_roles };
            md_parallel_for(ex, num_chunks, 1, infer_apply_templates, &task);
        }
        for (size_t u = 0; u < md_array_size(unknown_variants); ++u) {
            const unknown_variant_t* v = &unknown_variants[u];
            const md_urange_t range = md_system_residue_atom_range(sys, v->first);
            char beads[256];
            size_t len = 0;
            beads[0] = '\0';
            for (uint32_t i = range.beg; i < range.end && len + 8 < sizeof(beads); ++i) {
                len += (size_t)snprintf(beads + len, sizeof(beads) - len, "%s" MD_STR_FMT, len ? " " : "", MD_STR_ARG(label(label_data, i)));
            }
            MD_LOG_INFO(MD_STR_FMT " of %zu beads (%s), %zu of them: no Martini topology of the name has these beads; shown as regular beads, without bonds",
                MD_STR_ARG(v->name), v->size, beads, v->count);
        }
    } else {
        for (size_t i = 0; i < sys->atom.count; ++i) {
            if (sys->atom.type_idx[i] != 0) continue;

            const md_str_t atom_label = label(label_data, i);
            uint64_t key = md_hash64_str(atom_label, 0);
            uint32_t* cached_type = md_hashmap_get(&particle_type_cache, key);
            if (cached_type) {
                sys->atom.type_idx[i] = (md_particle_type_idx_t)*cached_type;
            } else {
                md_atomic_number_t z = md_atomic_number_infer_from_label(atom_label, (md_str_t){0}, 0);
                float mass   = md_atomic_number_mass(z);
                float radius = md_atomic_number_vdw_radius(z);
                uint32_t color = md_atomic_number_cpk_color(z);
                md_particle_type_idx_t type = md_particle_type_find_or_add(&sys->atom.type, atom_label, z, mass, radius, color, MD_PARTICLE_TYPE_FLAG_NONE, alloc);
                sys->atom.type_idx[i] = type;
                md_hashmap_add(&particle_type_cache, key, (uint32_t)type);
            }
        }
    }

    md_temp_end(temp);
}

void md_system_infer_particle_types(md_system_t* sys, const md_str_t atom_labels[], const md_executor_t* ex) {
    if (!atom_labels) {
        MD_LOG_ERROR("Must provide atom_labels");
        return;
    }
    infer_particle_types(sys, label_of_array, atom_labels, ex);
}

void md_system_infer_particle_types_from_labels(md_system_t* sys, const md_label_t atom_labels[], const md_executor_t* ex) {
    if (!atom_labels) {
        MD_LOG_ERROR("Must provide atom_labels");
        return;
    }
    infer_particle_types(sys, label_of_labels, atom_labels, ex);
}

void md_system_augment_particle_types(md_system_t* sys) {
    if (!sys || !sys->atom.type_idx || sys->residue.count == 0) {
        return;
    }
    md_system_topology_changed(sys);

    const size_t num_types = sys->atom.type.count;

    // The entry of an atom depends on its residue's name and its type alone, of which a system has few: each pair is
    // matched against the tables once (found, by a hash of the two), and while residues of the same name follow one
    // another, an atom takes its type's entry from the slot of its type, without a lookup. A slot holds the entry of the
    // run of residue names it was filled in (by the run's number: 0 for none), so a new run needs no clearing.
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_hashmap32_t found = { .allocator = md_temp_allocator(temp) };
    uint32_t* slot_run   = md_temp_alloc_array(temp, uint32_t, num_types + 1);
    uint32_t* slot_entry = md_temp_alloc_array(temp, uint32_t, num_types + 1);     // The entry + 1, 0 for none
    md_memset(slot_run, 0, (num_types + 1) * sizeof(uint32_t));
    uint32_t run = 0;
    md_str_t run_name = {0};
    uint64_t run_key = 0;

    for (size_t res_idx = 0; res_idx < sys->residue.count; ++res_idx) {
        const md_str_t res_name = md_system_residue_name(sys, res_idx);
        if (run == 0 || !md_str_eq(res_name, run_name)) {
            run += 1;
            run_name = res_name;
            run_key = md_hash64_str(res_name, 0);
        }
        const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
        md_residue_kind_t res_kind = MD_RESIDUE_KIND_OTHER;
        for (uint32_t i = range.beg; i < range.end; ++i) {
            const md_particle_type_idx_t t = sys->atom.type_idx[i];
            // Only types without an element: an atom with one is not what these tables describe
            if (t <= 0 || (size_t)t >= num_types || sys->atom.type.z[t] != 0) continue;

            if (slot_run[t] != run) {
                const uint64_t key = md_hash64(&t, sizeof(t), run_key);
                const uint32_t* known = md_hashmap_get(&found, key);
                uint32_t e;
                if (known) {
                    e = *known;
                } else {
                    const particle_type_t* match = find_predefined_particle_type(res_name, md_particle_type_name(&sys->atom.type, t));
                    e = match ? (uint32_t)(match - predefined_particle_types) + 1 : 0;
                    md_hashmap_add(&found, key, e);
                }
                slot_run[t] = run;
                slot_entry[t] = e;
            }
            if (slot_entry[t] == 0) continue;
            const particle_type_t* entry = &predefined_particle_types[slot_entry[t] - 1];

            // The loader's word on what the particle is (a tpr knows its virtual sites) stands over the table's
            if (sys->atom.type.flags && md_particle_type_particle_kind(&sys->atom.type, t) == MD_PARTICLE_ATOM) {
                sys->atom.type.flags[t] = md_particle_type_flags_set_particle_kind(sys->atom.type.flags[t], entry->particle);
            }
            if (sys->atom.flags) sys->atom.flags[i] |= entry->role;
            if (entry->res_kind != MD_RESIDUE_KIND_OTHER) res_kind = entry->res_kind;
        }
        apply_predefined_res_kind(sys, res_idx, res_kind);
    }
    md_temp_end(temp);
}
