#include <md_fragment.h>

#include <md_system.h>
#include <md_util.h>
#include <md_match.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_bitset.h>
#include <core/md_common.h>
#include <core/md_log.h>

#include <stdio.h>
#include <string.h>

// ---------------------------------------------------------------------------
// Elements
// ---------------------------------------------------------------------------

static inline bool frag_is_metal(md_atomic_number_t z) {
    return (z == 3 || z == 4) ||                // Li, Be
           (z >= 11 && z <= 13) ||              // Na, Mg, Al
           (z >= 19 && z <= 31) ||              // K .. Ga
           (z >= 37 && z <= 50) ||              // Rb .. Sn
           (z >= 55 && z <= 84) ||              // Cs .. Po
           (z >= 87);
}

static inline bool frag_is_halogen(md_atomic_number_t z) {
    return z == MD_Z_F || z == MD_Z_Cl || z == MD_Z_Br || z == MD_Z_I;
}

static inline bool frag_is_noble(md_atomic_number_t z) {
    return z == MD_Z_He || z == MD_Z_Ne || z == MD_Z_Ar || z == MD_Z_Kr || z == MD_Z_Xe || z == MD_Z_Rn;
}

// Ertl's heteroatom: anything but carbon and hydrogen - here also not a halogen (a single one goes with its
// atom, see the header), a metal (a unit of its own) or a noble gas
static inline bool frag_is_hetero(md_atomic_number_t z) {
    return z > 1 && z != MD_Z_C && !frag_is_halogen(z) && !frag_is_metal(z) && !frag_is_noble(z);
}

// A bond that is more than single for the purpose of marking: a double or triple, or one in resonance
static inline bool frag_bond_multiple(md_bond_flags_t f) {
    return md_bond_order(f) >= 2 || (f & (MD_BOND_FLAG_AROMATIC | MD_BOND_FLAG_DELOCALIZED)) != 0;
}

const char* md_fragment_kind_str(md_fragment_kind_t kind) {
    switch (kind) {
    case MD_FRAGMENT_MOLECULE_SET: return "molecules";
    case MD_FRAGMENT_MOLECULE:     return "molecule";
    case MD_FRAGMENT_RESIDUE:      return "residue";
    case MD_FRAGMENT_METAL:        return "metal";
    case MD_FRAGMENT_RING_SYSTEM:  return "ring system";
    case MD_FRAGMENT_GROUP:        return "functional group";
    case MD_FRAGMENT_LINKER:       return "linker";
    case MD_FRAGMENT_ALKYL:        return "alkyl";
    case MD_FRAGMENT_ATOM:         return "atom";
    default:                       return "";
    }
}

// ---------------------------------------------------------------------------
// The library the cores are named from, most specific first: identification takes the first entry which
// matches the core as a whole. Lowercase atoms must be aromatic; uppercase ones in a ring of the pattern are
// not held to be aliphatic, which is why the aromatic systems come before the saturated rings.
// ---------------------------------------------------------------------------

typedef struct frag_lib_entry_t {
    const char* smiles;
    const char* name;
} frag_lib_entry_t;

static const frag_lib_entry_t frag_library[] = {
    // fused aromatic systems
    { "c1ccc2cc3ccccc3cc2c1",       "anthracenyl" },
    { "c1ccc2c(c1)ccc1ccccc12",     "phenanthrenyl" },
    { "c1ccc2c(c1)c1ccccc1n2",      "carbazolyl" },
    { "c1ccc2c(c1)oc1ccccc12",      "dibenzofuranyl" },
    { "c1ccc2c(c1)sc1ccccc12",      "dibenzothiophenyl" },
    { "c1ccc2c(c1)c1ccccc1c2",      "fluorenyl" },
    { "c1ccc2ccccc2c1",             "naphthyl" },
    { "c1ccc2ncccc2c1",             "quinolinyl" },
    { "c1ccc2cnccc2c1",             "isoquinolinyl" },
    { "c1ccc2c(c1)ccn2",            "indolyl" },
    { "c1ccc2c(c1)ncn2",            "benzimidazolyl" },
    { "c1ccc2c(c1)ocn2",            "benzoxazolyl" },
    { "c1ccc2c(c1)scn2",            "benzothiazolyl" },
    { "c1ccc2c(c1)nsn2",            "benzothiadiazolyl" },
    { "c1ncc2ncnc2n1",              "purinyl" },
    // aromatic rings
    { "c1ccccc1",                   "phenyl" },
    { "c1ccncc1",                   "pyridyl" },
    { "c1cncnc1",                   "pyrimidinyl" },
    { "c1cnccn1",                   "pyrazinyl" },
    { "c1ccnnc1",                   "pyridazinyl" },
    { "c1ncncn1",                   "triazinyl" },
    { "c1ccsc1",                    "thienyl" },
    { "c1ccoc1",                    "furyl" },
    { "c1ccnc1",                    "pyrrolyl" },
    { "c1cncn1",                    "imidazolyl" },
    { "c1ccnn1",                    "pyrazolyl" },
    { "c1cscn1",                    "thiazolyl" },
    { "c1cocn1",                    "oxazolyl" },
    { "c1ncnn1",                    "triazolyl" },
    { "c1nnnn1",                    "tetrazolyl" },
    // saturated and partly saturated rings
    { "O=C1C=CC(=O)C=C1",           "benzoquinone" },
    { "O=C1CCCCC1",                 "cyclohexanone" },
    { "C1CCCCC1",                   "cyclohexyl" },
    { "C1CCCC1",                    "cyclopentyl" },
    { "C1CCC1",                     "cyclobutyl" },
    { "C1CC1",                      "cyclopropyl" },
    { "C1CNCCN1",                   "piperazinyl" },
    { "C1COCCN1",                   "morpholinyl" },
    { "C1CCNCC1",                   "piperidyl" },
    { "C1CCOCC1",                   "oxanyl" },
    { "C1CCNC1",                    "pyrrolidinyl" },
    { "C1CCOC1",                    "oxolanyl" },
    { "C1CO1",                      "epoxide" },
    // functional groups
    { "N(~O)~O",                    "nitro" },
    { "S(=O)(=O)N",                 "sulfonamide" },
    { "S(=O)(=O)[O-]",              "sulfonate" },
    { "S(=O)(=O)O",                 "sulfonic ester" },
    { "S(=O)=O",                    "sulfonyl" },
    { "S=O",                        "sulfinyl" },
    { "P(~O)(~O)(~O)~O",            "phosphate" },
    { "P(~O)(~O)~O",                "phosphonate" },
    { "P=O",                        "phosphoryl" },
    { "NC(=O)N",                    "urea" },
    { "NC(=O)O",                    "carbamate" },
    { "OC(=O)O",                    "carbonate" },
    { "NC(=N)N",                    "guanidine" },
    { "C(=O)N",                     "amide" },
    { "C(=S)N",                     "thioamide" },
    { "C(=N)N",                     "amidine" },
    { "C(=O)[O-]",                  "carboxylate" },
    { "C(=O)[OH]",                  "carboxylic acid" },
    { "C(=O)O",                     "ester" },
    { "C(=O)S",                     "thioester" },
    { "C(=O)OC=O",                  "anhydride" },
    { "N=C=O",                      "isocyanate" },
    { "N=C=S",                      "isothiocyanate" },
    { "N~N~N",                      "azide" },
    { "N=N",                        "azo" },
    { "C=N[OH]",                    "oxime" },
    { "C=NN",                       "hydrazone" },
    { "NN",                         "hydrazine" },
    { "C=N",                        "imine" },
    { "C#N",                        "nitrile" },
    { "C=CC=O",                     "enone" },
    { "C=O",                        "carbonyl" },
    { "C=S",                        "thiocarbonyl" },
    { "OO",                         "peroxide" },
    { "SS",                         "disulfide" },
    { "C(O)O",                      "acetal" },
    { "C(N)O",                      "hemiaminal" },
    { "C(F)(F)F",                   "trifluoromethyl" },
    { "C(F)F",                      "difluoromethyl" },
    { "C(Cl)(Cl)Cl",                "trichloromethyl" },
    { "C(Cl)Cl",                    "dichloromethyl" },
    { "C=CC=C",                     "diene" },
    { "C=C",                        "alkene" },
    { "C#C",                        "alkyne" },
    { "C=C=C",                      "allene" },
    { "NO",                         "hydroxylamine" },
};

// ---------------------------------------------------------------------------
// Work state
// ---------------------------------------------------------------------------

typedef struct frag_ctx_t {
    const md_system_t* sys;
    size_t n;
    md_atomic_number_t* z;
    uint32_t* nb_off;               // [n + 1] covalent neighbours
    uint32_t* nb;
    md_bond_flags_t* nb_flags;
    int32_t* unit_of;               // [n]
    int32_t* core_of;               // [n] the unit whose core the atom is in, -1 for none
    md_array(uint8_t) unit_kind;
    md_array(int32_t) unit_residue; // residue of a RESIDUE unit, -1 otherwise
    md_allocator_i* temp;
} frag_ctx_t;

static int32_t frag_new_unit(frag_ctx_t* c, md_fragment_kind_t kind, int32_t residue) {
    md_array_push(c->unit_kind, (uint8_t)kind, c->temp);
    md_array_push(c->unit_residue, residue, c->temp);
    return (int32_t)md_array_size(c->unit_kind) - 1;
}

static inline void frag_assign(frag_ctx_t* c, uint32_t atom, int32_t unit, bool core) {
    c->unit_of[atom] = unit;
    if (core) c->core_of[atom] = unit;
}

static inline bool frag_heavy(const frag_ctx_t* c, uint32_t i) {
    return c->z[i] > 1;
}

static uint32_t frag_heavy_degree(const frag_ctx_t* c, uint32_t i) {
    uint32_t d = 0;
    for (uint32_t k = c->nb_off[i]; k < c->nb_off[i + 1]; ++k) d += frag_heavy(c, c->nb[k]) ? 1 : 0;
    return d;
}

static bool frag_bonded(const frag_ctx_t* c, uint32_t a, uint32_t b) {
    for (uint32_t k = c->nb_off[a]; k < c->nb_off[a + 1]; ++k) if (c->nb[k] == b) return true;
    return false;
}

static uint32_t uf_find(uint32_t* p, uint32_t x) {
    while (p[x] != x) { p[x] = p[p[x]]; x = p[x]; }
    return x;
}
static void uf_union(uint32_t* p, uint32_t a, uint32_t b) {
    a = uf_find(p, a); b = uf_find(p, b);
    if (a != b) { if (a < b) p[b] = a; else p[a] = b; }
}

// Composition of a set of atoms, carbon first and the rest alphabetically (Hill), without hydrogens unless
// 'with_h'
static void frag_formula(char* buf, size_t cap, const frag_ctx_t* c, const uint32_t* atoms, size_t count, bool with_h) {
    uint32_t num[128] = {0};
    for (size_t i = 0; i < count; ++i) {
        const md_atomic_number_t z = c->z[atoms[i]];
        if (z < 128) num[z]++;
    }
    size_t len = 0;
    buf[0] = '\0';
    #define FRAG_APPEND(Z) do { \
        const md_str_t sym = md_util_element_symbol((md_element_t)(Z)); \
        int w = (num[Z] > 1) ? snprintf(buf + len, cap - len, "%.*s%u", (int)sym.len, sym.ptr, num[Z]) \
                             : snprintf(buf + len, cap - len, "%.*s", (int)sym.len, sym.ptr); \
        if (w > 0) len = MD_MIN(cap - 1, len + (size_t)w); \
    } while (0)
    if (num[MD_Z_C]) FRAG_APPEND(MD_Z_C);
    if (with_h && num[MD_Z_H]) FRAG_APPEND(MD_Z_H);
    // the rest alphabetically by symbol
    bool done[128] = {0};
    done[MD_Z_C] = done[MD_Z_H] = done[0] = true;
    for (;;) {
        int best = -1;
        md_str_t best_sym = {0};
        for (int z = 1; z < 128; ++z) {
            if (done[z] || !num[z]) continue;
            const md_str_t sym = md_util_element_symbol((md_element_t)z);
            if (best < 0 || md_str_cmp_lex(sym, best_sym) < 0) { best = z; best_sym = sym; }
        }
        if (best < 0) break;
        done[best] = true;
        FRAG_APPEND(best);
    }
    #undef FRAG_APPEND
}

// Hydrogens on an atom: the system's count where it has one, the bonded hydrogen atoms otherwise
static int frag_hydrogens(const frag_ctx_t* c, uint32_t i) {
    if (c->sys->atom.hydrogen_count) return md_atom_hydrogen_count(&c->sys->atom, i);
    int h = 0;
    for (uint32_t k = c->nb_off[i]; k < c->nb_off[i + 1]; ++k) h += c->z[c->nb[k]] == MD_Z_H ? 1 : 0;
    return h;
}

static void frag_name_single(char* buf, size_t cap, const frag_ctx_t* c, uint32_t i) {
    const md_atomic_number_t z = c->z[i];
    const uint32_t deg = frag_heavy_degree(c, i);
    const int charge = c->sys->atom.formal_charge ? md_atom_formal_charge(&c->sys->atom, i) : 0;
    const char* name = NULL;
    switch (z) {
    case MD_Z_N: name = (charge > 0 || deg + (uint32_t)frag_hydrogens(c, i) >= 4) ? "ammonium" : "amine"; break;
    case MD_Z_O: name = deg == 0 ? "water" : deg == 1 ? (charge < 0 ? "oxide" : "hydroxyl") : "ether"; break;
    case MD_Z_S: name = deg == 0 ? "sulfide" : deg == 1 ? (charge < 0 ? "thiolate" : "thiol") : "thioether"; break;
    case 15:     name = "phosphine"; break;
    case 5:      name = "boryl"; break;
    case 14:     name = "silyl"; break;
    case 34:     name = "selenyl"; break;
    default: break;
    }
    if (name) {
        snprintf(buf, cap, "%s", name);
    } else {
        const md_str_t sym = md_util_element_symbol((md_element_t)z);
        snprintf(buf, cap, "%.*s", (int)sym.len, sym.ptr);
    }
}

// ---------------------------------------------------------------------------

bool md_fragments_compute(md_fragments_t* out, const md_system_t* sys, md_allocator_i* alloc) {
    if (!out || !sys || !alloc) return false;
    md_memset(out, 0, sizeof(*out));
    out->alloc = alloc;

    const size_t n = md_atom_count(&sys->atom);
    out->num_atoms = n;
    if (n == 0) return true;

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* talloc = md_temp_allocator(temp);

    frag_ctx_t c = {0};
    c.sys = sys;
    c.n = n;
    c.temp = talloc;
    c.z = md_temp_alloc_array(temp, md_atomic_number_t, n);
    c.unit_of = md_temp_alloc_array(temp, int32_t, n);
    c.core_of = md_temp_alloc_array(temp, int32_t, n);
    for (size_t i = 0; i < n; ++i) {
        c.z[i] = md_atom_atomic_number(&sys->atom, i);
        c.unit_of[i] = -1;
        c.core_of[i] = -1;
    }

    // Covalent neighbours. Coordination bonds join nothing here: a metal and its ligands are units, and
    // molecules, of their own.
    c.nb_off = md_temp_alloc_zero_array(temp, uint32_t, n + 1);
    for (size_t i = 0; i < n; ++i) {
        uint32_t k = 0;
        for (md_bond_iter_t it = md_bond_iter(&sys->bond, i); md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
            if (md_bond_kind(md_bond_iter_bond_flags(&it)) == MD_BOND_KIND_COVALENT) k++;
        }
        c.nb_off[i + 1] = c.nb_off[i] + k;
    }
    c.nb = md_temp_alloc_array(temp, uint32_t, c.nb_off[n] ? c.nb_off[n] : 1);
    c.nb_flags = md_temp_alloc_array(temp, md_bond_flags_t, c.nb_off[n] ? c.nb_off[n] : 1);
    for (size_t i = 0; i < n; ++i) {
        uint32_t k = c.nb_off[i];
        for (md_bond_iter_t it = md_bond_iter(&sys->bond, i); md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
            const md_bond_flags_t f = (md_bond_flags_t)md_bond_iter_bond_flags(&it);
            if (md_bond_kind(f) != MD_BOND_KIND_COVALENT) continue;
            c.nb[k] = (uint32_t)md_bond_iter_atom_index(&it);
            c.nb_flags[k] = f;
            k++;
        }
    }

    // 1. Residues of polymers
    for (size_t r = 0; r < md_system_residue_count(sys); ++r) {
        const md_residue_kind_t kind = md_system_residue_kind(sys, r);
        if (kind != MD_RESIDUE_KIND_AMINO_ACID && kind != MD_RESIDUE_KIND_NUCLEOTIDE) continue;
        const md_urange_t range = md_system_residue_atom_range(sys, r);
        int32_t u = -1;
        for (uint32_t i = range.beg; i < range.end && i < n; ++i) {
            if (c.z[i] == 0 || c.unit_of[i] >= 0) continue;
            if (u < 0) u = frag_new_unit(&c, MD_FRAGMENT_RESIDUE, (int32_t)r);
            frag_assign(&c, i, u, true);
        }
    }

    // 2. Metals
    for (uint32_t i = 0; i < n; ++i) {
        if (c.unit_of[i] < 0 && frag_is_metal(c.z[i])) {
            frag_assign(&c, i, frag_new_unit(&c, MD_FRAGMENT_METAL, -1), true);
        }
    }

    // 3. Ring systems
    {
        md_ring_data_t local = {0};
        const md_ring_data_t* rings = &sys->ring;
        if (md_ring_data_count(rings) == 0 && sys->bond.count > 0 && sys->bond.conn.offset) {
            md_util_infer_rings(&local, &sys->bond, n, talloc);
            rings = &local;
        }

        uint32_t* uf = md_temp_alloc_array(temp, uint32_t, n);
        bool* in_ring = md_temp_alloc_zero_array(temp, bool, n);
        for (uint32_t i = 0; i < n; ++i) uf[i] = i;

        for (size_t r = 0; r < md_ring_data_count(rings); ++r) {
            const int32_t* atoms = md_ring_atoms(rings, r);
            const size_t size = md_ring_size(rings, r);
            bool ok = size >= 3;
            for (size_t k = 0; k < size && ok; ++k) {
                const uint32_t a = (uint32_t)atoms[k];
                const uint32_t b = (uint32_t)atoms[(k + 1) % size];
                // A ring through a metal is a chelate ring, and a ring through a residue already placed is
                // that residue's
                ok = a < n && b < n && c.unit_of[a] < 0 && frag_heavy(&c, a) && frag_bonded(&c, a, b);
            }
            if (!ok) continue;
            for (size_t k = 0; k < size; ++k) {
                in_ring[atoms[k]] = true;
                uf_union(uf, (uint32_t)atoms[k], (uint32_t)atoms[(k + 1) % size]);
            }
        }

        int32_t* root_unit = md_temp_alloc_array(temp, int32_t, n);
        for (uint32_t i = 0; i < n; ++i) root_unit[i] = -1;
        for (uint32_t i = 0; i < n; ++i) {
            if (!in_ring[i]) continue;
            const uint32_t r = uf_find(uf, i);
            if (root_unit[r] < 0) root_unit[r] = frag_new_unit(&c, MD_FRAGMENT_RING_SYSTEM, -1);
            frag_assign(&c, i, root_unit[r], true);
        }

        // Terminal atoms double bonded to a ring atom: the O of a quinone or a lactam
        for (uint32_t i = 0; i < n; ++i) {
            if (!in_ring[i]) continue;
            for (uint32_t k = c.nb_off[i]; k < c.nb_off[i + 1]; ++k) {
                const uint32_t j = c.nb[k];
                if (in_ring[j] || c.unit_of[j] >= 0 || !frag_heavy(&c, j) || frag_is_halogen(c.z[j])) continue;
                if (!frag_bond_multiple(c.nb_flags[k]) || frag_heavy_degree(&c, j) != 1) continue;
                frag_assign(&c, j, c.unit_of[i], true);
            }
        }

        if (rings == &local) md_ring_data_free(&local, talloc);
    }

    // 4. Ertl's marking, over what is left
    bool* marked = md_temp_alloc_zero_array(temp, bool, n);
    for (uint32_t i = 0; i < n; ++i) {
        if (c.unit_of[i] >= 0 || !frag_heavy(&c, i)) continue;
        if (frag_is_hetero(c.z[i])) {
            marked[i] = true;
            continue;
        }
        if (c.z[i] != MD_Z_C) continue;

        uint32_t num_halogen = 0, num_single_hetero = 0;
        bool all_single = true;
        for (uint32_t k = c.nb_off[i]; k < c.nb_off[i + 1]; ++k) {
            const uint32_t j = c.nb[k];
            const md_atomic_number_t zj = c.z[j];
            if (frag_bond_multiple(c.nb_flags[k])) {
                all_single = false;
                if (frag_is_hetero(zj) || zj == MD_Z_C) marked[i] = true;   // C=O, C=N, C#N, C=S; C=C, C#C
            } else if (zj == MD_Z_O || zj == MD_Z_N || zj == MD_Z_S) {
                num_single_hetero++;
            }
            if (frag_is_halogen(zj)) num_halogen++;
        }
        if (all_single && num_single_hetero >= 2) marked[i] = true;        // acetal
        if (num_halogen >= 2) {                                             // CX2, CX3
            marked[i] = true;
            for (uint32_t k = c.nb_off[i]; k < c.nb_off[i + 1]; ++k) {
                if (frag_is_halogen(c.z[c.nb[k]]) && c.unit_of[c.nb[k]] < 0) marked[c.nb[k]] = true;
            }
        }
    }

    uint32_t* stack = md_temp_alloc_array(temp, uint32_t, n);

    // Each connected set of marked atoms is a group
    for (uint32_t s = 0; s < n; ++s) {
        if (!marked[s] || c.unit_of[s] >= 0) continue;
        const int32_t u = frag_new_unit(&c, MD_FRAGMENT_GROUP, -1);
        size_t top = 0;
        stack[top++] = s;
        frag_assign(&c, s, u, true);
        while (top) {
            const uint32_t a = stack[--top];
            for (uint32_t k = c.nb_off[a]; k < c.nb_off[a + 1]; ++k) {
                const uint32_t b = c.nb[k];
                if (!marked[b] || c.unit_of[b] >= 0) continue;
                frag_assign(&c, b, u, true);
                stack[top++] = b;
            }
        }
    }

    // 5. The rest of the skeleton, in connected pieces (single halogens and hydrogens come after)
    {
        int32_t* piece = md_temp_alloc_array(temp, int32_t, n);
        uint32_t* members = md_temp_alloc_array(temp, uint32_t, n);
        for (uint32_t i = 0; i < n; ++i) piece[i] = -1;

        for (uint32_t s = 0; s < n; ++s) {
            if (c.unit_of[s] >= 0 || piece[s] >= 0 || !frag_heavy(&c, s) || frag_is_halogen(c.z[s]) || frag_is_noble(c.z[s])) continue;

            size_t num = 0, top = 0;
            stack[top++] = s;
            piece[s] = (int32_t)s;
            while (top) {
                const uint32_t a = stack[--top];
                members[num++] = a;
                for (uint32_t k = c.nb_off[a]; k < c.nb_off[a + 1]; ++k) {
                    const uint32_t b = c.nb[k];
                    if (c.unit_of[b] >= 0 || piece[b] >= 0 || !frag_heavy(&c, b) || frag_is_halogen(c.z[b]) || frag_is_noble(c.z[b])) continue;
                    piece[b] = (int32_t)s;
                    stack[top++] = b;
                }
            }

            // The units it touches
            int32_t touch[2] = { -1, -1 };
            uint32_t num_touch = 0;
            for (size_t m = 0; m < num; ++m) {
                const uint32_t a = members[m];
                for (uint32_t k = c.nb_off[a]; k < c.nb_off[a + 1]; ++k) {
                    const int32_t u = c.unit_of[c.nb[k]];
                    if (u < 0 || u == touch[0] || u == touch[1]) continue;
                    if (num_touch < 2) touch[num_touch] = u;
                    num_touch++;
                }
            }

            if (num <= 2 && num_touch == 1) {
                for (size_t m = 0; m < num; ++m) frag_assign(&c, members[m], touch[0], false);
            } else {
                const int32_t u = frag_new_unit(&c, num_touch >= 2 ? MD_FRAGMENT_LINKER : MD_FRAGMENT_ALKYL, -1);
                for (size_t m = 0; m < num; ++m) frag_assign(&c, members[m], u, true);
            }
        }
    }

    // 6. Single halogens with their atom, hydrogens with theirs, anything else on its own
    for (int pass = 0; pass < 2; ++pass) {
        for (uint32_t i = 0; i < n; ++i) {
            if (c.unit_of[i] >= 0 || c.z[i] == 0) continue;
            const bool hydrogen = c.z[i] == MD_Z_H;
            if (hydrogen != (pass == 1)) continue;
            int32_t u = -1;
            for (uint32_t k = c.nb_off[i]; k < c.nb_off[i + 1] && u < 0; ++k) {
                const uint32_t j = c.nb[k];
                if (c.unit_of[j] >= 0 && (frag_heavy(&c, j) || hydrogen)) u = c.unit_of[j];
            }
            if (u >= 0) {
                frag_assign(&c, i, u, false);
            } else {
                frag_assign(&c, i, frag_new_unit(&c, MD_FRAGMENT_ATOM, -1), true);
            }
        }
    }

    const size_t num_units = md_array_size(c.unit_kind);

    // Atoms of each unit, in index order
    uint32_t* unit_off = md_temp_alloc_zero_array(temp, uint32_t, num_units + 1);
    for (uint32_t i = 0; i < n; ++i) if (c.unit_of[i] >= 0) unit_off[c.unit_of[i] + 1]++;
    for (size_t u = 0; u < num_units; ++u) unit_off[u + 1] += unit_off[u];
    uint32_t* unit_atoms = md_temp_alloc_array(temp, uint32_t, unit_off[num_units] ? unit_off[num_units] : 1);
    {
        uint32_t* fill = md_temp_alloc_zero_array(temp, uint32_t, num_units + 1);
        for (uint32_t i = 0; i < n; ++i) {
            const int32_t u = c.unit_of[i];
            if (u >= 0) unit_atoms[unit_off[u] + fill[u]++] = i;
        }
    }

    // Molecules: covalent components, a unit always within one
    uint32_t* mol_uf = md_temp_alloc_array(temp, uint32_t, n);
    for (uint32_t i = 0; i < n; ++i) mol_uf[i] = i;
    for (uint32_t i = 0; i < n; ++i) {
        if (c.unit_of[i] < 0) continue;
        for (uint32_t k = c.nb_off[i]; k < c.nb_off[i + 1]; ++k) {
            if (c.unit_of[c.nb[k]] >= 0) uf_union(mol_uf, i, c.nb[k]);
        }
    }
    for (size_t u = 0; u < num_units; ++u) {
        for (uint32_t k = unit_off[u] + 1; k < unit_off[u + 1]; ++k) uf_union(mol_uf, unit_atoms[unit_off[u]], unit_atoms[k]);
    }

    // Units ordered by first atom are already in that order within a molecule if visited by first atom
    int32_t* mol_of_root = md_temp_alloc_array(temp, int32_t, n);
    for (uint32_t i = 0; i < n; ++i) mol_of_root[i] = -1;
    md_array(uint32_t) mol_first_unit_order = 0;    // units in order of first atom
    int32_t* unit_mol = md_temp_alloc_array(temp, int32_t, num_units ? num_units : 1);
    bool* unit_seen = md_temp_alloc_zero_array(temp, bool, num_units ? num_units : 1);
    size_t num_mols = 0;
    for (uint32_t i = 0; i < n; ++i) {
        const int32_t u = c.unit_of[i];
        if (u < 0 || unit_seen[u]) continue;
        unit_seen[u] = true;
        const uint32_t r = uf_find(mol_uf, i);
        if (mol_of_root[r] < 0) mol_of_root[r] = (int32_t)num_mols++;
        unit_mol[u] = mol_of_root[r];
        md_array_push(mol_first_unit_order, (uint32_t)u, talloc);
    }

    // ---- names of the units ----
    char* unit_name = md_temp_alloc_zero_array(temp, char, (num_units ? num_units : 1) * MD_FRAGMENT_NAME_CAP);
    {
        // Cores of two or more heavy atoms are identified against the library
        md_array(md_bitset_t) units = 0;
        md_array(uint32_t) unit_of_query = 0;
        uint32_t* core = md_temp_alloc_array(temp, uint32_t, n);

        for (size_t u = 0; u < num_units; ++u) {
            char* name = unit_name + u * MD_FRAGMENT_NAME_CAP;
            const md_fragment_kind_t kind = (md_fragment_kind_t)c.unit_kind[u];

            size_t num_core = 0;
            for (uint32_t k = unit_off[u]; k < unit_off[u + 1]; ++k) {
                const uint32_t a = unit_atoms[k];
                if (c.core_of[a] == (int32_t)u && frag_heavy(&c, a)) core[num_core++] = a;
            }

            switch (kind) {
            case MD_FRAGMENT_RESIDUE: {
                const int32_t r = c.unit_residue[u];
                const md_str_t rn = md_system_residue_name(sys, (size_t)r);
                snprintf(name, MD_FRAGMENT_NAME_CAP, "%.*s %d", (int)rn.len, rn.ptr, (int)md_system_residue_seq_id(sys, (size_t)r));
                break;
            }
            case MD_FRAGMENT_METAL:
            case MD_FRAGMENT_ATOM: {
                const md_str_t sym = md_util_element_symbol((md_element_t)c.z[unit_atoms[unit_off[u]]]);
                snprintf(name, MD_FRAGMENT_NAME_CAP, "%.*s", (int)sym.len, sym.ptr);
                break;
            }
            case MD_FRAGMENT_ALKYL: {
                static const char* alkyl[] = { "", "methyl", "ethyl", "propyl", "butyl", "pentyl", "hexyl" };
                if (num_core < MD_COUNTOF(alkyl)) snprintf(name, MD_FRAGMENT_NAME_CAP, "%s", alkyl[num_core]);
                else                              snprintf(name, MD_FRAGMENT_NAME_CAP, "C%zu alkyl", num_core);
                break;
            }
            case MD_FRAGMENT_LINKER: {
                if      (num_core == 1) snprintf(name, MD_FRAGMENT_NAME_CAP, "methylene");
                else if (num_core == 2) snprintf(name, MD_FRAGMENT_NAME_CAP, "ethylene");
                else                    snprintf(name, MD_FRAGMENT_NAME_CAP, "C%zu linker", num_core);
                break;
            }
            case MD_FRAGMENT_GROUP:
            case MD_FRAGMENT_RING_SYSTEM:
                if (num_core == 1) {
                    frag_name_single(name, MD_FRAGMENT_NAME_CAP, &c, core[0]);
                } else if (num_core > 1) {
                    // sorted, as md_bitset_from_indices wants them
                    for (size_t a = 1; a < num_core; ++a) {
                        const uint32_t v = core[a];
                        size_t b = a;
                        while (b > 0 && core[b - 1] > v) { core[b] = core[b - 1]; b--; }
                        core[b] = v;
                    }
                    md_array_push(units, md_bitset_from_indices(core, num_core, talloc), talloc);
                    md_array_push(unit_of_query, (uint32_t)u, talloc);
                    // the fallback, overwritten when identified
                    char formula[24];
                    frag_formula(formula, sizeof(formula), &c, core, num_core, false);
                    bool aromatic = false;
                    for (size_t a = 0; a < num_core; ++a) aromatic |= (md_atom_flags(&sys->atom, core[a]) & MD_ATOM_FLAG_AROMATIC) != 0;
                    bool carbon_only = true;
                    for (size_t a = 0; a < num_core; ++a) carbon_only &= c.z[core[a]] == MD_Z_C;
                    if (kind == MD_FRAGMENT_GROUP && carbon_only && num_core >= 4) {
                        snprintf(name, MD_FRAGMENT_NAME_CAP, "polyene");    // conjugated C=C beyond a diene
                    } else {
                        snprintf(name, MD_FRAGMENT_NAME_CAP, "%s %s",
                                 kind == MD_FRAGMENT_GROUP ? "group" : aromatic ? "aryl" : "ring", formula);
                    }
                }
                break;
            default:
                break;
            }
        }

        if (md_array_size(units) > 0) {
            md_match_library_t* lib = md_match_library_create(talloc);
            for (size_t e = 0; e < MD_COUNTOF(frag_library); ++e) {
                const int32_t idx = md_match_library_add_smiles(lib, md_str_from_cstr(frag_library[e].smiles),
                                                                md_str_from_cstr(frag_library[e].name), NULL);
                if (idx < 0) MD_LOG_ERROR("Fragments: library entry '%s' does not parse", frag_library[e].smiles);
            }
            md_match_identify_desc_t desc = {
                .library = lib,
                .level = MD_MATCH_LEVEL_MOLECULE,
                .flags = MD_MATCH_FLAG_WHOLE,
                .units = units,
                .num_units = md_array_size(units),
            };
            md_match_identify_result_t res = {0};
            if (md_match_identify(&res, &desc, sys, talloc)) {
                for (size_t k = 0; k < res.count; ++k) {
                    const uint32_t u = unit_of_query[res.unit[k]];
                    const md_str_t nm = md_match_library_name(lib, res.entry[k]);
                    snprintf(unit_name + u * MD_FRAGMENT_NAME_CAP, MD_FRAGMENT_NAME_CAP, "%.*s", (int)nm.len, nm.ptr);
                }
                md_match_identify_result_free(&res);
            }
            md_match_library_destroy(lib);
        }
    }

    // ---- molecules: their units, names, and sets of identical small ones ----
    uint32_t* mol_unit_off = md_temp_alloc_zero_array(temp, uint32_t, num_mols + 1);
    for (size_t u = 0; u < num_units; ++u) mol_unit_off[unit_mol[u] + 1]++;
    for (size_t m = 0; m < num_mols; ++m) mol_unit_off[m + 1] += mol_unit_off[m];
    uint32_t* mol_units = md_temp_alloc_array(temp, uint32_t, num_units ? num_units : 1);
    {
        uint32_t* fill = md_temp_alloc_zero_array(temp, uint32_t, num_mols + 1);
        for (size_t k = 0; k < md_array_size(mol_first_unit_order); ++k) {
            const uint32_t u = mol_first_unit_order[k];
            const int32_t m = unit_mol[u];
            mol_units[mol_unit_off[m] + fill[m]++] = u;
        }
    }

    // Repeated unit names within a molecule are numbered
    for (size_t m = 0; m < num_mols; ++m) {
        const uint32_t beg = mol_unit_off[m], end = mol_unit_off[m + 1];
        for (uint32_t a = beg; a < end; ++a) {
            char* na = unit_name + mol_units[a] * MD_FRAGMENT_NAME_CAP;
            if (c.unit_kind[mol_units[a]] == MD_FRAGMENT_RESIDUE) continue;
            uint32_t total = 0;
            for (uint32_t b = beg; b < end; ++b) total += strcmp(na, unit_name + mol_units[b] * MD_FRAGMENT_NAME_CAP) == 0 ? 1 : 0;
            if (total < 2) continue;
            // number all of them, in order; the copies compared against are still unnumbered
            char base[MD_FRAGMENT_NAME_CAP];
            snprintf(base, sizeof(base), "%s", na);
            uint32_t idx = 0;
            for (uint32_t b = beg; b < end; ++b) {
                char* nb_ = unit_name + mol_units[b] * MD_FRAGMENT_NAME_CAP;
                if (strcmp(nb_, base) != 0) continue;
                char tmp[MD_FRAGMENT_NAME_CAP + 8];
                snprintf(tmp, sizeof(tmp), "%s %u", base, ++idx);
                snprintf(nb_, MD_FRAGMENT_NAME_CAP, "%s", tmp);
            }
        }
    }

    char* mol_name = md_temp_alloc_zero_array(temp, char, (num_mols ? num_mols : 1) * MD_FRAGMENT_NAME_CAP);
    uint32_t* mol_atom_count = md_temp_alloc_zero_array(temp, uint32_t, num_mols + 1);
    for (size_t m = 0; m < num_mols; ++m) {
        const uint32_t beg = mol_unit_off[m], end = mol_unit_off[m + 1];
        char* name = mol_name + m * MD_FRAGMENT_NAME_CAP;
        md_array(uint32_t) atoms = 0;
        bool protein = false, nucleic = false;
        for (uint32_t k = beg; k < end; ++k) {
            const uint32_t u = mol_units[k];
            for (uint32_t a = unit_off[u]; a < unit_off[u + 1]; ++a) md_array_push(atoms, unit_atoms[a], talloc);
            if (c.unit_kind[u] == MD_FRAGMENT_RESIDUE) {
                const md_residue_kind_t rk = md_system_residue_kind(sys, (size_t)c.unit_residue[u]);
                protein |= rk == MD_RESIDUE_KIND_AMINO_ACID;
                nucleic |= rk == MD_RESIDUE_KIND_NUCLEOTIDE;
            }
        }
        mol_atom_count[m] = (uint32_t)md_array_size(atoms);
        if (end - beg == 1) {
            snprintf(name, MD_FRAGMENT_NAME_CAP, "%s", unit_name + mol_units[beg] * MD_FRAGMENT_NAME_CAP);
        } else if (protein) {
            snprintf(name, MD_FRAGMENT_NAME_CAP, "protein");
        } else if (nucleic) {
            snprintf(name, MD_FRAGMENT_NAME_CAP, "nucleic acid");
        } else {
            frag_formula(name, MD_FRAGMENT_NAME_CAP, &c, atoms, md_array_size(atoms), true);
        }
    }

    // Sets: single unit molecules of one name and size, two or more of them. set_of[m] is the first member.
    int32_t* set_of = md_temp_alloc_array(temp, int32_t, num_mols ? num_mols : 1);
    uint32_t* set_size = md_temp_alloc_zero_array(temp, uint32_t, num_mols ? num_mols : 1);
    for (size_t m = 0; m < num_mols; ++m) {
        set_of[m] = (int32_t)m;
        if (mol_unit_off[m + 1] - mol_unit_off[m] != 1) {
            set_size[m] = 1;
            continue;
        }
        set_of[m] = -1;
        for (size_t p = 0; p < m; ++p) {
            if (set_of[p] != (int32_t)p) continue;
            if (mol_atom_count[p] == mol_atom_count[m] &&
                strcmp(mol_name + p * MD_FRAGMENT_NAME_CAP, mol_name + m * MD_FRAGMENT_NAME_CAP) == 0) {
                set_of[m] = (int32_t)p;
                break;
            }
        }
        if (set_of[m] < 0) set_of[m] = (int32_t)m;
        set_size[set_of[m]]++;
    }

    // Repeated names among what stays a root are numbered too
    {
        for (size_t m = 0; m < num_mols; ++m) {
            if (set_of[m] != (int32_t)m || set_size[m] > 1) continue;
            const char* nm = mol_name + m * MD_FRAGMENT_NAME_CAP;
            uint32_t total = 0;
            for (size_t p = 0; p < num_mols; ++p) {
                if (set_of[p] == (int32_t)p && set_size[p] == 1 && strcmp(nm, mol_name + p * MD_FRAGMENT_NAME_CAP) == 0) total++;
            }
            if (total < 2) continue;
            char base[MD_FRAGMENT_NAME_CAP];
            snprintf(base, sizeof(base), "%s", nm);
            uint32_t idx = 0;
            for (size_t p = 0; p < num_mols; ++p) {
                char* np = mol_name + p * MD_FRAGMENT_NAME_CAP;
                if (set_of[p] != (int32_t)p || set_size[p] != 1 || strcmp(np, base) != 0) continue;
                char tmp[MD_FRAGMENT_NAME_CAP + 8];
                snprintf(tmp, sizeof(tmp), "%s %u", base, ++idx);
                snprintf(np, MD_FRAGMENT_NAME_CAP, "%s", tmp);
            }
        }
    }

    // ---- the tree, in pre-order ----
    md_array(uint32_t) parent = 0;
    md_array(uint8_t)  kind   = 0;
    md_array(char)     names  = 0;
    md_array(uint32_t) abeg   = 0;
    md_array(uint32_t) aend   = 0;
    md_array(int32_t)  aidx   = 0;
    out->atom_fragment = md_array_create(uint32_t, n, alloc);
    for (size_t i = 0; i < n; ++i) out->atom_fragment[i] = MD_FRAGMENT_NONE;

    #define FRAG_PUSH(PARENT, KIND, NAME) do { \
        md_array_push(parent, (PARENT), alloc); \
        md_array_push(kind, (uint8_t)(KIND), alloc); \
        const size_t base_ = md_array_size(names); \
        md_array_resize(names, base_ + MD_FRAGMENT_NAME_CAP, alloc); \
        snprintf(names + base_, MD_FRAGMENT_NAME_CAP, "%s", (NAME)); \
        md_array_push(abeg, (uint32_t)md_array_size(aidx), alloc); \
        md_array_push(aend, (uint32_t)md_array_size(aidx), alloc); \
    } while (0)

    #define FRAG_ATOMS_OF_UNIT(FRAG, U) do { \
        for (uint32_t k_ = unit_off[(U)]; k_ < unit_off[(U) + 1]; ++k_) { \
            md_array_push(aidx, (int32_t)unit_atoms[k_], alloc); \
            out->atom_fragment[unit_atoms[k_]] = (FRAG); \
        } \
    } while (0)

    // ends the ranges of a fragment and all of its ancestors at the current end
    #define FRAG_CLOSE(FRAG) do { \
        for (uint32_t f_ = (FRAG); f_ != MD_FRAGMENT_NONE; f_ = parent[f_]) aend[f_] = (uint32_t)md_array_size(aidx); \
    } while (0)

    for (size_t m = 0; m < num_mols; ++m) {
        const bool in_set = set_size[set_of[m]] > 1;
        if (in_set && set_of[m] != (int32_t)m) continue;    // emitted with the first of its set

        if (in_set) {
            char set_name[MD_FRAGMENT_NAME_CAP];
            snprintf(set_name, sizeof(set_name), "%s \xc3\x97%u", mol_name + m * MD_FRAGMENT_NAME_CAP, set_size[m]);
            const uint32_t set_idx = (uint32_t)md_array_size(parent);
            FRAG_PUSH(MD_FRAGMENT_NONE, MD_FRAGMENT_MOLECULE_SET, set_name);
            uint32_t member = 0;
            for (size_t p = m; p < num_mols; ++p) {
                if (set_of[p] != (int32_t)m) continue;
                const uint32_t u = mol_units[mol_unit_off[p]];
                char member_name[MD_FRAGMENT_NAME_CAP + 12];
                snprintf(member_name, sizeof(member_name), "%s %u", mol_name + m * MD_FRAGMENT_NAME_CAP, ++member);
                char clipped[MD_FRAGMENT_NAME_CAP];
                snprintf(clipped, sizeof(clipped), "%s", member_name);
                const uint32_t f = (uint32_t)md_array_size(parent);
                FRAG_PUSH(set_idx, c.unit_kind[u], clipped);
                FRAG_ATOMS_OF_UNIT(f, u);
                FRAG_CLOSE(f);
            }
            continue;
        }

        const uint32_t beg = mol_unit_off[m], end = mol_unit_off[m + 1];
        if (end - beg == 1) {
            const uint32_t u = mol_units[beg];
            const uint32_t f = (uint32_t)md_array_size(parent);
            FRAG_PUSH(MD_FRAGMENT_NONE, c.unit_kind[u], mol_name + m * MD_FRAGMENT_NAME_CAP);
            FRAG_ATOMS_OF_UNIT(f, u);
            FRAG_CLOSE(f);
        } else {
            const uint32_t mf = (uint32_t)md_array_size(parent);
            FRAG_PUSH(MD_FRAGMENT_NONE, MD_FRAGMENT_MOLECULE, mol_name + m * MD_FRAGMENT_NAME_CAP);
            for (uint32_t k = beg; k < end; ++k) {
                const uint32_t u = mol_units[k];
                const uint32_t f = (uint32_t)md_array_size(parent);
                FRAG_PUSH(mf, c.unit_kind[u], unit_name + u * MD_FRAGMENT_NAME_CAP);
                FRAG_ATOMS_OF_UNIT(f, u);
                FRAG_CLOSE(f);
            }
        }
    }
    #undef FRAG_PUSH
    #undef FRAG_ATOMS_OF_UNIT
    #undef FRAG_CLOSE

    out->count        = md_array_size(parent);
    out->parent       = parent;
    out->kind         = kind;
    out->name         = names;
    out->atom_beg     = abeg;
    out->atom_end     = aend;
    out->atom_idx     = aidx;
    out->num_atom_idx = md_array_size(aidx);

    md_temp_end(temp);
    return true;
}

void md_fragments_free(md_fragments_t* f) {
    if (!f) return;
    if (f->alloc) {
        md_array_free(f->parent, f->alloc);
        md_array_free(f->kind, f->alloc);
        md_array_free(f->name, f->alloc);
        md_array_free(f->atom_beg, f->alloc);
        md_array_free(f->atom_end, f->alloc);
        md_array_free(f->atom_idx, f->alloc);
        md_array_free(f->atom_fragment, f->alloc);
    }
    md_memset(f, 0, sizeof(*f));
}
