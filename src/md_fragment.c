#include <md_fragment.h>

#include <md_system.h>
#include <md_match.h>
#include <md_element.h>

#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_bitset.h>
#include <core/md_common.h>
#include <core/md_hash.h>
#include <core/md_log.h>
#include <core/md_os.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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
// atom, see the header), a metal (a fragment of its own) or a noble gas
static inline bool frag_is_hetero(md_atomic_number_t z) {
    return z > 1 && z != MD_Z_C && !frag_is_halogen(z) && !frag_is_metal(z) && !frag_is_noble(z);
}

// A bond that is more than single for the purpose of marking: a double or triple, or one in resonance
static inline bool frag_bond_multiple(md_bond_flags_t f) {
    return md_bond_order(f) >= 2 || (f & (MD_BOND_FLAG_AROMATIC | MD_BOND_FLAG_DELOCALIZED)) != 0;
}

const char* md_fragment_kind_str(md_fragment_kind_t kind) {
    switch (kind) {
    case MD_FRAGMENT_BACKBONE:    return "backbone";
    case MD_FRAGMENT_LIPID_TAIL:  return "lipid tail";
    case MD_FRAGMENT_RING_SYSTEM: return "ring system";
    case MD_FRAGMENT_GROUP:       return "functional group";
    case MD_FRAGMENT_LINKER:      return "linker";
    case MD_FRAGMENT_ALKYL:       return "alkyl";
    case MD_FRAGMENT_ATOM:        return "atom";
    default:                      return "";
    }
}

static int frag_cmp_u64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
    return (x > y) - (x < y);
}

static int frag_cmp_u32(const void* a, const void* b) {
    const uint32_t x = *(const uint32_t*)a, y = *(const uint32_t*)b;
    return (x > y) - (x < y);
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
    { "C1~C~C2~C~C3~C~C~C(~N3)~C~C3~C~C~C(~N3)~C~C3~C~C~C(~N3)~C~C1~N2", "porphyrin" },
    { "B1~N2~C~C~C~C~2~C~C2~C~C~C~N~1~2", "BODIPY" },
    { "C1~C~C~C2~C(~C1)~C~C~C1~C2~C~C~C2~C~C~C~C~1~2", "steroid" },
    { "O~C1~N~C(~O)~C2~N~C~N~C~2~N1", "xanthine" },
    { "O~C1~C~C~N~C2~C~C~C~C~C~1~2", "quinolone" },
    { "O~C1~C~C2~S~C~C~N~1~2",      "penam" },
    { "c1ncc2ncnc2n1",              "purinyl" },
    // the ring systems of the nucleobases, with their oxygens (an amine on them is a group of its own): of guanine,
    // of thymine and uracil, of cytosine. Any bonds, as the tautomer perceived varies.
    { "O~C1~N~C~N~C2~C1~N~C~N2",    "purinone" },
    { "O~C1~N~C(~O)~N~C~C1",        "pyrimidinedione" },
    { "O~C1~N~C~C~C~N1",            "pyrimidinone" },
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
    { "c1ccon1",                    "isoxazolyl" },
    { "c1ccsn1",                    "isothiazolyl" },
    { "c1ncnn1",                    "triazolyl" },
    { "c1nnnn1",                    "tetrazolyl" },
    // saturated and partly saturated rings
    { "O=C1C=CC(=O)C=C1",           "benzoquinone" },
    { "O=C1CCCCC1",                 "cyclohexanone" },
    { "C1CCCCC1",                   "cyclohexyl" },
    { "C1CCCC1",                    "cyclopentyl" },
    { "C1CCC1",                     "cyclobutyl" },
    { "C1CC1",                      "cyclopropyl" },
    { "C1=CCCCC1",                  "cyclohexenyl" },
    { "C1=CCCC1",                   "cyclopentenyl" },
    { "C1CNCCN1",                   "piperazinyl" },
    { "C1COCCN1",                   "morpholinyl" },
    { "C1CCNCC1",                   "piperidyl" },
    { "C1CCOCC1",                   "oxanyl" },
    { "C1~C~C~O~C~C~1",             "pyranyl" },     // unsaturated, after oxanyl
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
    { "O~P(~O)(~O)~O~P(~O)(~O)~O~P(~O)(~O)~O", "triphosphate" },
    { "O~P(~O)(~O)~O~P(~O)(~O)~O",  "diphosphate" },
    { "P(~O)(~O)(~O)~O",            "phosphate" },
    { "P(~S)(~O)(~O)~O",            "thiophosphate" },
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
    { "C(~C#N)~C#N",                "dicyanomethylene" },
    { "C#N",                        "nitrile" },
    { "C=CC=O",                     "enone" },
    { "C=O",                        "carbonyl" },
    { "C=S",                        "thiocarbonyl" },
    { "OO",                         "peroxide" },
    { "SS",                         "disulfide" },
    { "C(O)O",                      "acetal" },
    { "C(N)O",                      "hemiaminal" },
    { "C(F)(F)(F)C(~O)~O",          "trifluoroacetyl" },
    { "C(F)(F)F",                   "trifluoromethyl" },
    { "C(F)F",                      "difluoromethyl" },
    { "C(Cl)(Cl)Cl",                "trichloromethyl" },
    { "C(Cl)Cl",                    "dichloromethyl" },
    { "C(Br)(Br)Br",                "tribromomethyl" },
    { "C=CC=C",                     "diene" },
    { "C=C",                        "alkene" },
    { "C#C",                        "alkyne" },
    { "C=C=C",                      "allene" },
    { "NO",                         "hydroxylamine" },
    { "N=O",                        "nitroso" },
};

// ---------------------------------------------------------------------------
// One molecule, in local indices (its slots)
// ---------------------------------------------------------------------------

#define FRAG_NAME_CAP 48

// The atom flags a template depends on
#define FRAG_ROLE_MASK (MD_ATOM_FLAG_BACKBONE | MD_ATOM_FLAG_SIDE_CHAIN | MD_ATOM_FLAG_NUCLEOSIDE | MD_ATOM_FLAG_NUCLEOBASE | \
                        MD_ATOM_FLAG_TERMINAL_BEG | MD_ATOM_FLAG_TERMINAL_END | \
                        MD_ATOM_FLAG_AROMATIC)

typedef struct frag_mol_t {
    const md_system_t* sys;
    size_t           n;
    const int32_t*   atom;          // [n] global atom of each slot
    const int32_t*   atom_res;      // [atom count] residue of each atom, -1 for none (shared)
    const uint8_t*   lipid_tail;    // [atom count] 1 for an atom of a tail of a lipid, NULL for none (shared)
    md_atomic_number_t* z;          // [n]
    uint32_t*        nb_off;        // [n + 1] covalent neighbours, local
    uint32_t*        nb;
    md_bond_flags_t* nb_flags;
    int32_t*         unit_of;       // [n] the fragment (local) of each atom, -1 for none yet
    int32_t*         core_of;       // [n] the fragment whose core the atom is in, -1 for none
    md_array(uint8_t) unit_kind;
    md_array(uint8_t) unit_role;    // frag_role_t: what names a fragment other than its core
    md_allocator_i*  temp;
} frag_mol_t;

// Fragments named by what they are in a polymer rather than by their core
typedef enum frag_role_t {
    FRAG_ROLE_NONE = 0,
    FRAG_ROLE_ALPHA_CARBON,         // CA of an amino acid, with its hydrogen and a side chain of one carbon
    FRAG_ROLE_SUGAR,                // The sugar of a nucleotide: its ring and C5'
} frag_role_t;

// The covalent graph of molecule m in local indices. Coordination bonds join nothing here: a metal and its ligands
// are fragments of their own.
// lipid_tail: the atoms of the tails of the lipids (md_lipid_tail_atom_mask), NULL for none
static void frag_mol_init(frag_mol_t* c, const md_system_t* sys, size_t m, const int32_t* atom_res, const uint8_t* lipid_tail, md_allocator_i* temp) {
    const md_molecule_data_t* mol = &sys->molecule;
    const uint32_t off = mol->offset[m];
    const size_t n = mol->offset[m + 1] - off;
    md_memset(c, 0, sizeof(*c));
    c->sys = sys;
    c->n = n;
    c->atom = mol->atom_idx + off;
    c->atom_res = atom_res;
    c->lipid_tail = lipid_tail;
    c->temp = temp;
    c->z = md_alloc(temp, sizeof(md_atomic_number_t) * (n ? n : 1));
    c->nb_off = md_alloc(temp, sizeof(uint32_t) * (n + 1));
    c->unit_of = md_alloc(temp, sizeof(int32_t) * (n ? n : 1));
    c->core_of = md_alloc(temp, sizeof(int32_t) * (n ? n : 1));
    c->nb_off[0] = 0;
    for (size_t i = 0; i < n; ++i) {
        c->z[i] = md_atom_atomic_number(&sys->atom, (size_t)c->atom[i]);
        c->unit_of[i] = -1;
        c->core_of[i] = -1;
        uint32_t k = 0;
        for (md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)c->atom[i]); md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
            if (md_bond_kind(md_bond_iter_bond_flags(&it)) != MD_BOND_KIND_COVALENT) continue;
            const int32_t s = mol->atom_slot[md_bond_iter_atom_index(&it)] - (int32_t)off;
            k += (s >= 0 && (size_t)s < n);
        }
        c->nb_off[i + 1] = c->nb_off[i] + k;
    }
    const size_t nn = c->nb_off[n];
    c->nb = md_alloc(temp, sizeof(uint32_t) * (nn ? nn : 1));
    c->nb_flags = md_alloc(temp, sizeof(md_bond_flags_t) * (nn ? nn : 1));
    for (size_t i = 0; i < n; ++i) {
        uint32_t k = c->nb_off[i];
        for (md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)c->atom[i]); md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
            const md_bond_flags_t f = (md_bond_flags_t)md_bond_iter_bond_flags(&it);
            if (md_bond_kind(f) != MD_BOND_KIND_COVALENT) continue;
            const int32_t s = mol->atom_slot[md_bond_iter_atom_index(&it)] - (int32_t)off;
            if (s < 0 || (size_t)s >= n) continue;
            c->nb[k] = (uint32_t)s;
            c->nb_flags[k] = f;
            k++;
        }
    }
}

static inline md_atom_flags_t frag_aflags(const frag_mol_t* c, uint32_t i) {
    return md_atom_flags(&c->sys->atom, (size_t)c->atom[i]);
}

static inline bool frag_lipid_tail(const frag_mol_t* c, uint32_t i) {
    return c->lipid_tail && c->lipid_tail[c->atom[i]];
}

static inline int32_t frag_res(const frag_mol_t* c, uint32_t i) {
    return c->atom_res ? c->atom_res[c->atom[i]] : -1;
}

static int32_t frag_new_unit_role(frag_mol_t* c, md_fragment_kind_t kind, frag_role_t role) {
    md_array_push(c->unit_kind, (uint8_t)kind, c->temp);
    md_array_push(c->unit_role, (uint8_t)role, c->temp);
    return (int32_t)md_array_size(c->unit_kind) - 1;
}

static int32_t frag_new_unit(frag_mol_t* c, md_fragment_kind_t kind) {
    return frag_new_unit_role(c, kind, FRAG_ROLE_NONE);
}

static inline void frag_assign(frag_mol_t* c, uint32_t i, int32_t unit, bool core) {
    c->unit_of[i] = unit;
    if (core) c->core_of[i] = unit;
}

static inline bool frag_heavy(const frag_mol_t* c, uint32_t i) {
    return c->z[i] > 1;
}

// The kind of the residue of an atom where it is a resolved amino acid or nucleotide (its atoms have roles),
// MD_RESIDUE_KIND_OTHER otherwise
static md_residue_kind_t frag_resolved_kind(const frag_mol_t* c, uint32_t i) {
    const int32_t r = c->atom_res ? c->atom_res[c->atom[i]] : -1;
    if (r < 0) return MD_RESIDUE_KIND_OTHER;
    const md_residue_flags_t rf = md_system_residue_flags(c->sys, (size_t)r);
    if (!(rf & MD_RESIDUE_FLAG_RESOLVED)) return MD_RESIDUE_KIND_OTHER;
    const md_residue_kind_t rk = md_system_residue_kind(c->sys, (size_t)r);
    return (rk == MD_RESIDUE_KIND_AMINO_ACID || rk == MD_RESIDUE_KIND_NUCLEOTIDE) ? rk : MD_RESIDUE_KIND_OTHER;
}

static inline bool frag_on_backbone(const frag_mol_t* c, uint32_t i) {
    return (md_atom_flags(&c->sys->atom, (size_t)c->atom[i]) & MD_ATOM_FLAG_BACKBONE) && frag_resolved_kind(c, i) != MD_RESIDUE_KIND_OTHER;
}

static uint32_t frag_heavy_degree(const frag_mol_t* c, uint32_t i) {
    uint32_t d = 0;
    for (uint32_t k = c->nb_off[i]; k < c->nb_off[i + 1]; ++k) d += frag_heavy(c, c->nb[k]) ? 1 : 0;
    return d;
}

static bool frag_bonded(const frag_mol_t* c, uint32_t a, uint32_t b) {
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

// Composition of a set of atoms without hydrogens, carbon first and the rest alphabetically (Hill)
static void frag_formula(char* buf, size_t cap, const frag_mol_t* c, const uint32_t* atoms, size_t count) {
    // The elements present and their counts: a handful, in a short list kept in the order written
    md_atomic_number_t el[16];
    uint32_t num[16];
    size_t num_el = 0;
    for (size_t i = 0; i < count; ++i) {
        const md_atomic_number_t z = c->z[atoms[i]];
        if (z == MD_Z_H || z == 0) continue;
        size_t k = 0;
        while (k < num_el && el[k] != z) ++k;
        if (k == num_el) {
            if (num_el == MD_COUNTOF(el)) continue;
            el[num_el] = z;
            num[num_el++] = 0;
        }
        num[k]++;
    }
    // carbon first, the rest by symbol
    for (size_t i = 1; i < num_el; ++i) {
        const md_atomic_number_t z = el[i];
        const uint32_t v = num[i];
        size_t j = i;
        for (; j > 0; --j) {
            const md_atomic_number_t y = el[j - 1];
            const bool before = z == MD_Z_C ? y != MD_Z_C : (y != MD_Z_C && md_str_cmp_lex(md_atomic_number_symbol((md_atomic_number_t)z), md_atomic_number_symbol((md_atomic_number_t)y)) < 0);
            if (!before) break;
            el[j] = el[j - 1];
            num[j] = num[j - 1];
        }
        el[j] = z;
        num[j] = v;
    }
    size_t len = 0;
    for (size_t k = 0; k < num_el; ++k) {
        const md_str_t sym = md_atomic_number_symbol((md_atomic_number_t)el[k]);
        for (size_t s_ = 0; s_ < sym.len && len + 1 < cap; ++s_) buf[len++] = sym.ptr[s_];
        if (num[k] > 1) {
            char digits[12];
            int nd = 0;
            for (uint32_t v = num[k]; v; v /= 10) digits[nd++] = (char)('0' + v % 10);
            while (nd && len + 1 < cap) buf[len++] = digits[--nd];
        }
    }
    buf[len] = '\0';
}

// Hydrogens on an atom: the system's count where it has one, the bonded hydrogen atoms otherwise
static int frag_hydrogens(const frag_mol_t* c, uint32_t i) {
    if (c->sys->atom.hydrogen_count) return md_atom_hydrogen_count(&c->sys->atom, (size_t)c->atom[i]);
    int h = 0;
    for (uint32_t k = c->nb_off[i]; k < c->nb_off[i + 1]; ++k) h += c->z[c->nb[k]] == MD_Z_H ? 1 : 0;
    return h;
}

static void frag_name_single(char* buf, size_t cap, const frag_mol_t* c, uint32_t i) {
    const md_atomic_number_t z = c->z[i];
    const uint32_t deg = frag_heavy_degree(c, i);
    const int charge = c->sys->atom.formal_charge ? md_atom_formal_charge(&c->sys->atom, (size_t)c->atom[i]) : 0;
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
        const md_str_t sym = md_atomic_number_symbol((md_atomic_number_t)z);
        snprintf(buf, cap, "%.*s", (int)sym.len, sym.ptr);
    }
}

// The bonds between heavy atoms which are on no cycle (Tarjan's bridges, iteratively): bridge[k] for c->nb[k]
static void frag_find_bridges(bool* bridge, const frag_mol_t* c, md_allocator_i* temp) {
    const size_t n = c->n;
    const size_t nn = c->nb_off[n];
    uint32_t* disc = md_alloc(temp, sizeof(uint32_t) * (n ? n : 1));
    uint32_t* low  = md_alloc(temp, sizeof(uint32_t) * (n ? n : 1));
    uint32_t* next = md_alloc(temp, sizeof(uint32_t) * (n ? n : 1));     // the adjacency entry to look at next
    uint32_t* from = md_alloc(temp, sizeof(uint32_t) * (n ? n : 1));     // the parent's adjacency entry it was reached by
    uint32_t* stack = md_alloc(temp, sizeof(uint32_t) * (n ? n : 1));
    for (size_t i = 0; i < n; ++i) disc[i] = UINT32_MAX;
    for (size_t k = 0; k < nn; ++k) bridge[k] = false;
    uint32_t time = 0;
    for (uint32_t s = 0; s < n; ++s) {
        if (disc[s] != UINT32_MAX || !frag_heavy(c, s)) continue;
        size_t top = 0;
        stack[top++] = s;
        disc[s] = low[s] = time++;
        next[s] = c->nb_off[s];
        from[s] = UINT32_MAX;
        while (top) {
            const uint32_t v = stack[top - 1];
            if (next[v] < c->nb_off[v + 1]) {
                const uint32_t k = next[v]++;
                const uint32_t w = c->nb[k];
                if (!frag_heavy(c, w)) continue;
                if (disc[w] == UINT32_MAX) {
                    disc[w] = low[w] = time++;
                    next[w] = c->nb_off[w];
                    from[w] = k;
                    stack[top++] = w;
                } else {
                    // A back edge, unless it is the edge back to the parent v was reached from
                    const uint32_t parent = top >= 2 ? stack[top - 2] : UINT32_MAX;
                    if (w != parent) low[v] = MD_MIN(low[v], disc[w]);
                }
            } else {
                --top;
                if (top) {
                    const uint32_t parent = stack[top - 1];
                    low[parent] = MD_MIN(low[parent], low[v]);
                    if (low[v] > disc[parent]) {
                        // parent - v is a bridge: both of its adjacency entries
                        bridge[from[v]] = true;
                        for (uint32_t k = c->nb_off[v]; k < c->nb_off[v + 1]; ++k) if (c->nb[k] == parent) bridge[k] = true;
                    }
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// The fragments of one molecule: c->unit_of and c->unit_kind. rings: the rings of the molecule (global atoms).
// ---------------------------------------------------------------------------

static void frag_analyse(frag_mol_t* c, const md_ring_data_t* rings, const uint32_t* ring_idx, size_t num_rings,
                         int32_t* res_unit /* [residue count], all -1, left so */) {
    const size_t n = c->n;
    const md_system_t* sys = c->sys;
    md_allocator_i* temp = c->temp;
    uint32_t* stack = md_alloc(temp, sizeof(uint32_t) * (n ? n : 1));
    md_array(int32_t) touched_res = 0;

    // 0a. Of the backbone of a polymer whose residues are resolved, what Ertl's marking would not make sense of: the
    // alpha carbon of each amino acid, and the sugar of each nucleotide (its ring with C5'), which would otherwise be
    // a ring system and a methylene linker. The rest of the backbone is found as groups below: the peptide groups
    // and the phosphates, each between two residues, and the termini.
    for (uint32_t i = 0; i < n; ++i) {
        if (c->unit_of[i] >= 0 || !frag_heavy(c, i)) continue;
        const md_residue_kind_t rk = frag_resolved_kind(c, i);
        if (rk == MD_RESIDUE_KIND_OTHER) continue;
        const md_atom_flags_t af = frag_aflags(c, i);
        if (rk == MD_RESIDUE_KIND_AMINO_ACID) {
            // CA: the carbon of the backbone bonded to no oxygen (C is bonded to its O)
            if (c->z[i] != MD_Z_C || !(af & MD_ATOM_FLAG_BACKBONE)) continue;
            bool to_o = false;
            for (uint32_t k = c->nb_off[i]; k < c->nb_off[i + 1]; ++k) to_o |= c->z[c->nb[k]] == MD_Z_O;
            if (to_o) continue;
            frag_assign(c, i, frag_new_unit_role(c, MD_FRAGMENT_BACKBONE, FRAG_ROLE_ALPHA_CARBON), true);
        } else {
            // The ring of the sugar (C1', C2', C3', C4', O4' are of the nucleoside, not of the base) and C5' (a carbon
            // of the backbone), one fragment per residue
            const bool ring = (af & MD_ATOM_FLAG_NUCLEOSIDE) && !(af & MD_ATOM_FLAG_NUCLEOBASE);
            const bool c5 = c->z[i] == MD_Z_C && (af & MD_ATOM_FLAG_BACKBONE);
            if (!ring && !c5) continue;
            const int32_t r = frag_res(c, i);
            if (res_unit[r] < 0) {
                res_unit[r] = frag_new_unit_role(c, MD_FRAGMENT_BACKBONE, FRAG_ROLE_SUGAR);
                md_array_push(touched_res, r, temp);
            }
            frag_assign(c, i, res_unit[r], true);
        }
    }
    for (size_t k = 0; k < md_array_size(touched_res); ++k) res_unit[touched_res[k]] = -1;

    // 0b. Lipid tails: each connected set of their heavy atoms
    for (uint32_t s = 0; s < n; ++s) {
        if (c->unit_of[s] >= 0 || !frag_heavy(c, s) || !frag_lipid_tail(c, s)) continue;
        const int32_t u = frag_new_unit(c, MD_FRAGMENT_LIPID_TAIL);
        size_t top = 0;
        stack[top++] = s;
        frag_assign(c, s, u, true);
        while (top) {
            const uint32_t a = stack[--top];
            for (uint32_t k = c->nb_off[a]; k < c->nb_off[a + 1]; ++k) {
                const uint32_t b = c->nb[k];
                if (c->unit_of[b] >= 0 || !frag_heavy(c, b) || !frag_lipid_tail(c, b)) continue;
                frag_assign(c, b, u, true);
                stack[top++] = b;
            }
        }
    }

    // 1. Metals
    for (uint32_t i = 0; i < n; ++i) {
        if (c->unit_of[i] < 0 && frag_is_metal(c->z[i])) frag_assign(c, i, frag_new_unit(c, MD_FRAGMENT_ATOM), true);
    }

    // 2. Ring systems
    {
        uint32_t* uf = md_alloc(temp, sizeof(uint32_t) * (n ? n : 1));
        bool* in_ring = md_alloc(temp, sizeof(bool) * (n ? n : 1));
        for (uint32_t i = 0; i < n; ++i) { uf[i] = i; in_ring[i] = false; }
        const md_molecule_data_t* mol = &sys->molecule;
        int32_t local[MD_RING_MAX_SIZE];
        for (size_t q = 0; q < num_rings; ++q) {
            const size_t r = ring_idx[q];
            const int32_t* atoms = md_ring_atoms(rings, r);
            const size_t size = md_ring_size(rings, r);
            bool ok = size >= 3 && size <= MD_RING_MAX_SIZE;
            for (size_t k = 0; k < size && ok; ++k) {
                const int32_t s = mol->atom_slot[atoms[k]] - (int32_t)(c->atom - mol->atom_idx);
                ok = s >= 0 && (size_t)s < n;
                if (ok) local[k] = s;
            }
            for (size_t k = 0; k < size && ok; ++k) {
                const uint32_t a = (uint32_t)local[k], b = (uint32_t)local[(k + 1) % size];
                // A ring through a metal is a chelate ring, and one through a backbone or a tail is theirs
                ok = c->unit_of[a] < 0 && frag_heavy(c, a) && frag_bonded(c, a, b);
            }
            if (!ok) continue;
            for (size_t k = 0; k < size; ++k) {
                in_ring[local[k]] = true;
                uf_union(uf, (uint32_t)local[k], (uint32_t)local[(k + 1) % size]);
            }
        }
        int32_t* root_unit = md_alloc(temp, sizeof(int32_t) * (n ? n : 1));
        for (uint32_t i = 0; i < n; ++i) root_unit[i] = -1;
        for (uint32_t i = 0; i < n; ++i) {
            if (!in_ring[i]) continue;
            const uint32_t r = uf_find(uf, i);
            if (root_unit[r] < 0) root_unit[r] = frag_new_unit(c, MD_FRAGMENT_RING_SYSTEM);
            frag_assign(c, i, root_unit[r], true);
        }
        // The conjugated bridges of a macrocycle, which is no ring here (rings are of at most MD_RING_MAX_SIZE atoms):
        // an atom outside the rings with a multiple or aromatic bond, on a cycle (by no bridge bond) between two
        // atoms of ring systems, joins them and they one another. The meso carbons of a porphyrin make it one ring
        // system, as the aza bridges of a phthalocyanine; a C=C between two rings on no cycle (a stilbene) does not.
        {
            bool* bridge = md_alloc(temp, sizeof(bool) * (c->nb_off[n] ? c->nb_off[n] : 1));
            frag_find_bridges(bridge, c, temp);
            // The ring systems joined, as a union-find over the units so far: a framework of thousands of bridges
            // costs a pass, not a pass per bridge
            const uint32_t num_units = (uint32_t)md_array_size(c->unit_kind);
            uint32_t* join = md_alloc(temp, sizeof(uint32_t) * (num_units ? num_units : 1));
            for (uint32_t u = 0; u < num_units; ++u) join[u] = u;
            bool joined = false;
            for (uint32_t i = 0; i < n; ++i) {
                if (c->unit_of[i] >= 0 || in_ring[i] || !frag_heavy(c, i)) continue;
                bool conjugated = false;
                int32_t sys_a = -1;
                uint32_t num_ring_nb = 0;
                for (uint32_t k = c->nb_off[i]; k < c->nb_off[i + 1]; ++k) {
                    conjugated |= frag_bond_multiple(c->nb_flags[k]);
                    const uint32_t j = c->nb[k];
                    if (bridge[k] || !in_ring[j]) continue;
                    const int32_t u = c->unit_of[j];
                    if (u < 0 || c->unit_kind[u] != MD_FRAGMENT_RING_SYSTEM) continue;
                    num_ring_nb++;
                    if (sys_a < 0) sys_a = u;
                    else uf_union(join, (uint32_t)sys_a, (uint32_t)u);
                }
                if (!conjugated || num_ring_nb < 2) continue;
                frag_assign(c, i, sys_a, true);
                in_ring[i] = true;
                joined = true;
            }
            if (joined) {
                for (uint32_t a = 0; a < n; ++a) {
                    if (c->unit_of[a] >= 0 && (uint32_t)c->unit_of[a] < num_units) c->unit_of[a] = (int32_t)uf_find(join, (uint32_t)c->unit_of[a]);
                    if (c->core_of[a] >= 0 && (uint32_t)c->core_of[a] < num_units) c->core_of[a] = (int32_t)uf_find(join, (uint32_t)c->core_of[a]);
                }
            }
        }
        // Terminal atoms double bonded to a ring atom: the O of a quinone or a lactam
        for (uint32_t i = 0; i < n; ++i) {
            if (!in_ring[i]) continue;
            for (uint32_t k = c->nb_off[i]; k < c->nb_off[i + 1]; ++k) {
                const uint32_t j = c->nb[k];
                if (in_ring[j] || c->unit_of[j] >= 0 || !frag_heavy(c, j) || frag_is_halogen(c->z[j])) continue;
                if (!frag_bond_multiple(c->nb_flags[k]) || frag_heavy_degree(c, j) != 1) continue;
                frag_assign(c, j, c->unit_of[i], true);
            }
        }
    }

    // 3. Ertl's marking, over what is left
    bool* marked = md_alloc(temp, sizeof(bool) * (n ? n : 1));
    for (uint32_t i = 0; i < n; ++i) marked[i] = false;
    for (uint32_t i = 0; i < n; ++i) {
        if (c->unit_of[i] >= 0 || !frag_heavy(c, i)) continue;
        if (frag_is_hetero(c->z[i])) {
            marked[i] = true;
            continue;
        }
        if (c->z[i] != MD_Z_C) continue;
        uint32_t num_halogen = 0, num_single_hetero = 0;
        bool all_single = true;
        for (uint32_t k = c->nb_off[i]; k < c->nb_off[i + 1]; ++k) {
            const md_atomic_number_t zj = c->z[c->nb[k]];
            if (frag_bond_multiple(c->nb_flags[k])) {
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
            for (uint32_t k = c->nb_off[i]; k < c->nb_off[i + 1]; ++k) {
                if (frag_is_halogen(c->z[c->nb[k]]) && c->unit_of[c->nb[k]] < 0) marked[c->nb[k]] = true;
            }
        }
    }
    // each connected set of marked atoms is a group
    for (uint32_t s = 0; s < n; ++s) {
        if (!marked[s] || c->unit_of[s] >= 0) continue;
        const int32_t u = frag_new_unit(c, MD_FRAGMENT_GROUP);
        size_t top = 0;
        stack[top++] = s;
        frag_assign(c, s, u, true);
        while (top) {
            const uint32_t a = stack[--top];
            for (uint32_t k = c->nb_off[a]; k < c->nb_off[a + 1]; ++k) {
                const uint32_t b = c->nb[k];
                if (!marked[b] || c->unit_of[b] >= 0) continue;
                frag_assign(c, b, u, true);
                stack[top++] = b;
            }
        }
    }

    // 4. The rest of the skeleton, in connected pieces (single halogens and hydrogens come after)
    {
        int32_t* piece = md_alloc(temp, sizeof(int32_t) * (n ? n : 1));
        uint32_t* members = md_alloc(temp, sizeof(uint32_t) * (n ? n : 1));
        for (uint32_t i = 0; i < n; ++i) piece[i] = -1;
        #define FRAG_SKELETON(A) (c->unit_of[A] < 0 && piece[A] < 0 && frag_heavy(c, A) && !frag_is_halogen(c->z[A]) && !frag_is_noble(c->z[A]))
        for (uint32_t s = 0; s < n; ++s) {
            if (!FRAG_SKELETON(s)) continue;
            size_t num = 0, top = 0;
            stack[top++] = s;
            piece[s] = (int32_t)s;
            while (top) {
                const uint32_t a = stack[--top];
                members[num++] = a;
                for (uint32_t k = c->nb_off[a]; k < c->nb_off[a + 1]; ++k) {
                    const uint32_t b = c->nb[k];
                    if (!FRAG_SKELETON(b)) continue;
                    piece[b] = (int32_t)s;
                    stack[top++] = b;
                }
            }
            // The fragments it touches: the backbone apart, which takes no side chain (a CB goes with its side chain);
            // and of them the functional groups
            int32_t touch[2] = { -1, -1 };
            uint32_t num_touch = 0;
            int32_t group[2] = { -1, -1 };
            uint32_t num_group = 0;
            int32_t backbone = -1;
            uint32_t num_backbone = 0;
            for (size_t m = 0; m < num; ++m) {
                const uint32_t a = members[m];
                for (uint32_t k = c->nb_off[a]; k < c->nb_off[a + 1]; ++k) {
                    const uint32_t j = c->nb[k];
                    const int32_t u = c->unit_of[j];
                    if (u < 0) continue;
                    if (c->unit_kind[u] == MD_FRAGMENT_BACKBONE) {
                        if (u != backbone) { backbone = u; num_backbone++; }
                        continue;
                    }
                    if (c->unit_kind[u] == MD_FRAGMENT_GROUP && u != group[0] && u != group[1]) {
                        if (num_group < 2) group[num_group] = u;
                        num_group++;
                    }
                    if (u == touch[0] || u == touch[1]) continue;
                    if (num_touch < 2) touch[num_touch] = u;
                    num_touch++;
                }
            }
            // A methyl or ethyl on a functional group and on nothing else is of it, as a chemist names them
            // together: the methyls of a dimethylamino, the methyl of a methoxy or an acetyl, the ethyl of an ethyl
            // ester. A piece between the group and something else (a ring, the backbone, another group) is a linker:
            // the carboxylate of a glutamate is its COO-, the CH2CH2 between it and the alpha carbon a linker, as the
            // longer ones of a lysine or an arginine are. Nor does a ring or the backbone take a piece (toluene is a
            // phenyl and a methyl, the CB of an alanine a methyl).
            if (num <= 2 && num_group == 1 && num_touch == 1 && num_backbone == 0) {
                for (size_t m = 0; m < num; ++m) frag_assign(c, members[m], group[0], false);
            } else {
                const int32_t u = frag_new_unit(c, num_touch + num_backbone >= 2 ? MD_FRAGMENT_LINKER : MD_FRAGMENT_ALKYL);
                for (size_t m = 0; m < num; ++m) frag_assign(c, members[m], u, true);
            }
        }
        #undef FRAG_SKELETON
    }

    // 5. Single halogens with their atom, hydrogens with theirs, anything else on its own
    for (int pass = 0; pass < 2; ++pass) {
        for (uint32_t i = 0; i < n; ++i) {
            if (c->unit_of[i] >= 0 || c->z[i] == 0) continue;
            const bool hydrogen = c->z[i] == MD_Z_H;
            if (hydrogen != (pass == 1)) continue;
            int32_t u = -1;
            for (uint32_t k = c->nb_off[i]; k < c->nb_off[i + 1] && u < 0; ++k) {
                const uint32_t j = c->nb[k];
                if (c->unit_of[j] >= 0 && (frag_heavy(c, j) || hydrogen)) u = c->unit_of[j];
            }
            if (u >= 0) frag_assign(c, i, u, false);
            else        frag_assign(c, i, frag_new_unit(c, MD_FRAGMENT_ATOM), true);
        }
    }

    // 6. A group of the backbone (a peptide group, a phosphate, a terminus) is of the backbone
    for (uint32_t i = 0; i < n; ++i) {
        const int32_t u = c->unit_of[i];
        if (u >= 0 && c->unit_kind[u] == MD_FRAGMENT_GROUP && frag_on_backbone(c, i)) c->unit_kind[u] = MD_FRAGMENT_BACKBONE;
    }
}

// ---------------------------------------------------------------------------
// Encoding: what makes two molecules alike, in slot order
// ---------------------------------------------------------------------------

// metallic: the atoms of the metallic phases (md_bond_metallic_atoms), NULL for none
static void frag_encode(md_array(uint64_t)* enc, const md_system_t* sys, size_t m, const int32_t* atom_res, const uint8_t* metallic, const uint8_t* lipid_tail, md_allocator_i* temp) {
    const md_molecule_data_t* mol = &sys->molecule;
    const uint32_t off = mol->offset[m];
    const size_t n = mol->offset[m + 1] - off;
    md_array_shrink(*enc, 0);
    md_array_push(*enc, (uint64_t)n, temp);
    for (size_t i = 0; i < n; ++i) {
        const size_t a = (size_t)mol->atom_idx[off + i];
        const uint64_t z = md_atom_atomic_number(&sys->atom, a);
        const uint64_t q = sys->atom.formal_charge ? (uint8_t)(int8_t)md_atom_formal_charge(&sys->atom, a) : 0;
        const uint64_t h = sys->atom.hydrogen_count ? (uint8_t)md_atom_hydrogen_count(&sys->atom, a) : 0;
        const uint64_t f = (uint64_t)(md_atom_flags(&sys->atom, a) & FRAG_ROLE_MASK) & 0xFFFF;
        uint64_t rk = 0;
        const int32_t r = atom_res ? atom_res[a] : -1;
        if (r >= 0) {
            const md_residue_flags_t rf = md_system_residue_flags(sys, (size_t)r);
            rk = (uint64_t)md_system_residue_kind(sys, (size_t)r) | ((rf & MD_RESIDUE_FLAG_RESOLVED) ? 8u : 0u);
        }
        const uint64_t mt = metallic && metallic[a] ? 1u : 0u;
        const uint64_t lt = lipid_tail && lipid_tail[a] ? 1u : 0u;
        md_array_push(*enc, z | (q << 8) | (h << 16) | (f << 24) | (rk << 40) | (mt << 44) | (lt << 45), temp);
    }
    for (size_t i = 0; i < n; ++i) {
        // The bonds of slot i to later slots, in order of those: made in slot order, and each atom's few sorted
        // in place, which orders them all, whatever the order of the connectivity
        const size_t beg = md_array_size(*enc);
        for (md_bond_iter_t it = md_bond_iter(&sys->bond, (size_t)mol->atom_idx[off + i]); md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
            const md_bond_flags_t bf = (md_bond_flags_t)md_bond_iter_bond_flags(&it);
            if (md_bond_kind(bf) != MD_BOND_KIND_COVALENT) continue;
            const int32_t j = mol->atom_slot[md_bond_iter_atom_index(&it)] - (int32_t)off;
            if (j <= (int32_t)i || (size_t)j >= n) continue;
            const uint64_t b = (uint64_t)md_bond_order(bf) | ((bf & MD_BOND_FLAG_AROMATIC) ? 4u : 0u) | ((bf & MD_BOND_FLAG_DELOCALIZED) ? 8u : 0u);
            md_array_push(*enc, ((uint64_t)i << 40) | ((uint64_t)j << 16) | b, temp);
        }
        uint64_t* e = *enc;
        for (size_t k = beg + 1; k < md_array_size(*enc); ++k) {
            const uint64_t v = e[k];
            size_t q = k;
            for (; q > beg && e[q - 1] > v; --q) e[q] = e[q - 1];
            e[q] = v;
        }
    }
}

// ---------------------------------------------------------------------------
// Types
// ---------------------------------------------------------------------------

static inline uint64_t frag_type_key(md_str_t name, md_fragment_kind_t kind) {
    uint64_t k = md_hash64_str(name, (uint64_t)kind + 1);
    return k >= MD_HASH_TOMBSTONE ? k - 2 : k;
}

static uint32_t frag_default_color(md_str_t name) {
    // The golden ratio spreads the hues of names which hash alike; channels within [0.45, 0.9]
    const uint32_t h = (uint32_t)md_hash64_str(name, 0) * 0x9E3779B9u;
    const float hue = (float)(h >> 8) * (1.0f / 16777216.0f) * 6.0f;
    const float s = 0.5f, v = 0.9f;
    const int sector = (int)hue % 6;
    const float f = hue - floorf(hue);
    const float p = v * (1 - s), q = v * (1 - s * f), t = v * (1 - s * (1 - f));
    float r, g, b;
    switch (sector) {
    case 0:  r = v; g = t; b = p; break;
    case 1:  r = q; g = v; b = p; break;
    case 2:  r = p; g = v; b = t; break;
    case 3:  r = p; g = q; b = v; break;
    case 4:  r = t; g = p; b = v; break;
    default: r = v; g = p; b = q; break;
    }
    return (0xFFu << 24) | ((uint32_t)(b * 255.0f + 0.5f) << 16) | ((uint32_t)(g * 255.0f + 0.5f) << 8) | (uint32_t)(r * 255.0f + 0.5f);
}

void md_fragment_data_free(md_fragment_data_t* d, md_allocator_i* alloc) {
    if (!d || !alloc) return;
    md_array_free(d->mol_template, alloc);
    md_array_free(d->frag_offset, alloc);
    md_array_free(d->template_molecule, alloc);
    md_array_free(d->template_slot_offset, alloc);
    md_array_free(d->slot_frag, alloc);
    md_array_free(d->template_frag_offset, alloc);
    md_array_free(d->frag_type, alloc);
    md_array_free(d->frag_slot_offset, alloc);
    md_array_free(d->frag_slot, alloc);
    md_array_free(d->type.kind, alloc);
    for (size_t i = 0; i < md_array_size(d->type.name); ++i) {
        if (!md_str_empty(d->type.name[i])) md_free_str(alloc, d->type.name[i]);
    }
    md_array_free(d->type.name, alloc);
    md_array_free(d->type.color, alloc);
    md_memset(d, 0, sizeof(*d));
}

// ---------------------------------------------------------------------------
// Inference
// ---------------------------------------------------------------------------

bool md_util_system_infer_fragments(md_system_t* sys) {
    if (!sys || !sys->alloc) return false;
    md_allocator_i* alloc = sys->alloc;

    // Scratch on a temporary arena other than the system's: the system may live on one (a scope begun on it would
    // take back what is given to the system when it ends). One scope, for all of it: a second, nested one would have
    // no third arena to go to.
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* talloc = md_temp_allocator(temp);
    md_tick_t tick = md_tick_now();
    double ms_encode = 0, ms_analyse = 0, ms_library = 0;

    // The colours of the types there were, by kind and name

    md_hashmap32_t old_color = { .allocator = talloc };
    for (size_t t = 0; t < sys->fragment.type.count; ++t) {
        md_hashmap_add(&old_color, frag_type_key(sys->fragment.type.name[t], (md_fragment_kind_t)sys->fragment.type.kind[t]), sys->fragment.type.color[t]);
    }
    md_fragment_data_free(&sys->fragment, alloc);
    md_fragment_data_t* out = &sys->fragment;
    md_system_topology_changed(sys);

    const md_molecule_data_t* mol = &sys->molecule;
    const size_t num_mols = mol->count;
    const size_t num_atoms = sys->atom.count;
    if (num_mols == 0 || !mol->offset || !mol->atom_slot || num_atoms == 0) {
        md_temp_end(temp);
        return true;
    }

    // The residue of each atom
    int32_t* atom_res = NULL;
    if (sys->residue.count && sys->residue.atom_offset) {
        atom_res = md_temp_alloc_array(temp, int32_t, num_atoms);
        for (size_t i = 0; i < num_atoms; ++i) atom_res[i] = -1;
        for (size_t r = 0; r < sys->residue.count; ++r) {
            const md_urange_t range = md_system_residue_atom_range(sys, r);
            for (uint32_t a = range.beg; a < range.end && a < num_atoms; ++a) atom_res[a] = (int32_t)r;
        }
    }

    // The atoms of the metallic phases, which a template tells apart from the same metal outside of one
    uint8_t* metallic = md_temp_alloc_array(temp, uint8_t, num_atoms);
    if (!md_bond_metallic_atoms(metallic, &sys->bond, num_atoms)) metallic = NULL;
    // The atoms of the tails of the lipids, of a fragment of their own
    uint8_t* lipid_tail = md_temp_alloc_array(temp, uint8_t, num_atoms);
    if (!md_lipid_tail_atom_mask(lipid_tail, &sys->lipid, num_atoms)) lipid_tail = NULL;

    // ---- templates: the molecules by their encoding ----
    md_array(uint32_t) tmpl_enc_off = 0;    // [template + 1] its encoding in enc_pool
    md_array(uint64_t) enc_pool = 0;
    md_array(uint32_t) tmpl_next = 0;       // the next template of the same hash
    md_array(uint64_t) enc = 0;
    md_hashmap32_t by_hash = { .allocator = talloc };
    md_array_push(tmpl_enc_off, 0, talloc);
    md_array_resize(out->mol_template, num_mols, alloc);

    for (size_t m = 0; m < num_mols; ++m) {
        frag_encode(&enc, sys, m, atom_res, metallic, lipid_tail, talloc);
        const size_t len = md_array_size(enc);
        uint64_t h = md_hash64(enc, len * sizeof(uint64_t), 0);
        if (h >= MD_HASH_TOMBSTONE) h -= 2;
        uint32_t t = UINT32_MAX;
        const uint32_t* first = md_hashmap_get(&by_hash, h);
        for (uint32_t c = first ? *first : UINT32_MAX; c != UINT32_MAX; c = tmpl_next[c]) {
            const uint32_t beg = tmpl_enc_off[c], end = tmpl_enc_off[c + 1];
            if (end - beg == len && memcmp(enc_pool + beg, enc, len * sizeof(uint64_t)) == 0) { t = c; break; }
        }
        if (t == UINT32_MAX) {
            t = (uint32_t)md_array_size(out->template_molecule);
            md_array_push(out->template_molecule, (uint32_t)m, alloc);
            md_array_push_array(enc_pool, enc, len, talloc);
            md_array_push(tmpl_enc_off, (uint32_t)md_array_size(enc_pool), talloc);
            md_array_push(tmpl_next, first ? *first : UINT32_MAX, talloc);
            if (first) *(uint32_t*)md_hashmap_get(&by_hash, h) = t;
            else       md_hashmap_add(&by_hash, h, t);
        }
        out->mol_template[m] = t;
    }
    const size_t num_tmpl = md_array_size(out->template_molecule);
    out->template_count = num_tmpl;
    ms_encode = md_tick_to_milliseconds(md_tick_now() - tick);
    tick = md_tick_now();

    // The rings of each template's molecule, by the molecule of their first atom
    const size_t num_rings = md_ring_data_count(&sys->ring);
    uint32_t* ring_off = md_temp_alloc_zero_array(temp, uint32_t, num_mols + 1);
    uint32_t* ring_list = md_temp_alloc_array(temp, uint32_t, num_rings ? num_rings : 1);
    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t m = md_system_molecule_find_by_atom_idx(sys, (size_t)md_ring_atoms(&sys->ring, r)[0]);
        if (m >= 0) ring_off[m + 1]++;
    }
    for (size_t m = 0; m < num_mols; ++m) ring_off[m + 1] += ring_off[m];
    {
        uint32_t* fill = md_temp_alloc_zero_array(temp, uint32_t, num_mols + 1);
        for (size_t r = 0; r < num_rings; ++r) {
            const int32_t m = md_system_molecule_find_by_atom_idx(sys, (size_t)md_ring_atoms(&sys->ring, r)[0]);
            if (m >= 0) ring_list[ring_off[m] + fill[m]++] = (uint32_t)r;
        }
    }

    int32_t* res_unit = md_temp_alloc_array(temp, int32_t, sys->residue.count ? sys->residue.count : 1);
    for (size_t r = 0; r < sys->residue.count; ++r) res_unit[r] = -1;

    // ---- the analysis, per template; names of the multi atom cores identified together after ----
    md_array(char) names = 0;                   // FRAG_NAME_CAP per template fragment
    md_array(uint8_t) kinds = 0;
    md_array(md_bitset_t) queries = 0;
    md_array(uint32_t) query_frag = 0;          // the template fragment of each query (the first of its cores)
    // Cores alike are identified once: a core by its atoms (element, charge, hydrogens, aromaticity) in slot order
    // and its bonds, compared in full where the hashes agree. Alike is exact, so equal cores are equal graphs and
    // name alike; the same core in another order is only identified again.
    md_array(uint32_t) frag_query = 0;          // the query of each template fragment, UINT32_MAX for none
    md_array(uint64_t) core_pool = 0;           // the encodings of the queries
    md_array(uint32_t) core_off = 0;            // [query + 1]
    md_array(uint32_t) core_next = 0;           // the next query of the same hash
    md_array(uint64_t) core_enc = 0;
    md_hashmap32_t core_by_hash = { .allocator = talloc };
    md_array_push(core_off, 0, talloc);

    md_array_push(out->template_slot_offset, 0, alloc);
    md_array_push(out->template_frag_offset, 0, alloc);
    md_array_push(out->frag_slot_offset, 0, alloc);

    for (size_t t = 0; t < num_tmpl; ++t) {
        const size_t m = out->template_molecule[t];
        md_allocator_i* salloc = talloc;     // a template's scratch: the unique molecules', in all

        frag_mol_t c;
        frag_mol_init(&c, sys, m, atom_res, lipid_tail, salloc);
        frag_analyse(&c, &sys->ring, ring_list + ring_off[m], ring_off[m + 1] - ring_off[m], res_unit);
        const size_t n = c.n;
        const size_t num_units = md_array_size(c.unit_kind);

        // Fragments in the order of their first slot
        int32_t* order = md_alloc(salloc, sizeof(int32_t) * (num_units ? num_units : 1));
        for (size_t u = 0; u < num_units; ++u) order[u] = -1;
        uint32_t num_frags = 0;
        for (size_t i = 0; i < n; ++i) {
            const int32_t u = c.unit_of[i];
            if (u >= 0 && order[u] < 0) order[u] = (int32_t)num_frags++;
        }

        const uint32_t frag_base = (uint32_t)md_array_size(kinds);
        // Grown by doubling (md_array_resize alone would grow them to the size asked, a copy per template)
        md_array_ensure(kinds, frag_base + num_frags, talloc);
        md_array_resize(kinds, frag_base + num_frags, talloc);
        md_array_ensure(frag_query, frag_base + num_frags, talloc);
        md_array_resize(frag_query, frag_base + num_frags, talloc);
        for (uint32_t f = frag_base; f < frag_base + num_frags; ++f) frag_query[f] = UINT32_MAX;
        md_array_ensure(names, (size_t)(frag_base + num_frags) * FRAG_NAME_CAP, talloc);
        md_array_resize(names, (size_t)(frag_base + num_frags) * FRAG_NAME_CAP, talloc);
        for (size_t u = 0; u < num_units; ++u) if (order[u] >= 0) kinds[frag_base + order[u]] = c.unit_kind[u];

        // slot -> fragment
        for (size_t i = 0; i < n; ++i) {
            const int32_t u = c.unit_of[i];
            md_array_push(out->slot_frag, u >= 0 ? (uint32_t)order[u] : MD_FRAGMENT_NONE, alloc);
        }
        md_array_push(out->template_slot_offset, (uint32_t)md_array_size(out->slot_frag), alloc);

        // fragment -> slots, ascending
        uint32_t* cnt = md_alloc(salloc, sizeof(uint32_t) * (num_frags + 1));
        md_memset(cnt, 0, sizeof(uint32_t) * (num_frags + 1));
        for (size_t i = 0; i < n; ++i) if (c.unit_of[i] >= 0) cnt[order[c.unit_of[i]] + 1]++;
        for (uint32_t f = 0; f < num_frags; ++f) cnt[f + 1] += cnt[f];
        const uint32_t slot_base = (uint32_t)md_array_size(out->frag_slot);
        md_array_ensure(out->frag_slot, slot_base + cnt[num_frags], alloc);
        md_array_resize(out->frag_slot, slot_base + cnt[num_frags], alloc);
        for (uint32_t f = 0; f < num_frags; ++f) md_array_push(out->frag_slot_offset, slot_base + cnt[f + 1], alloc);
        {
            uint32_t* fill = md_alloc(salloc, sizeof(uint32_t) * (num_frags ? num_frags : 1));
            md_memset(fill, 0, sizeof(uint32_t) * (num_frags ? num_frags : 1));
            for (size_t i = 0; i < n; ++i) {
                if (c.unit_of[i] < 0) continue;
                const uint32_t f = (uint32_t)order[c.unit_of[i]];
                out->frag_slot[slot_base + cnt[f] + fill[f]++] = (uint32_t)i;
            }
        }
        md_array_push(out->template_frag_offset, frag_base + num_frags, alloc);

        // names
        uint32_t* core = md_alloc(salloc, sizeof(uint32_t) * (n ? n : 1));
        uint32_t* core_global = md_alloc(salloc, sizeof(uint32_t) * (n ? n : 1));
        for (size_t u = 0; u < num_units; ++u) {
            if (order[u] < 0) continue;
            const uint32_t f = frag_base + (uint32_t)order[u];
            char* name = names + (size_t)f * FRAG_NAME_CAP;
            const md_fragment_kind_t kind = (md_fragment_kind_t)c.unit_kind[u];
            // its slots, ascending, from the fragment -> slots table above
            size_t num_core = 0;
            const uint32_t* fs = out->frag_slot + slot_base + cnt[order[u]];
            const uint32_t num_fs = cnt[order[u] + 1] - cnt[order[u]];
            const uint32_t first = num_fs ? fs[0] : UINT32_MAX;
            for (uint32_t q = 0; q < num_fs; ++q) {
                const uint32_t i = fs[q];
                if (c.core_of[i] == (int32_t)u && frag_heavy(&c, i)) core[num_core++] = i;
            }
            const frag_role_t role = (frag_role_t)c.unit_role[u];
            if (role == FRAG_ROLE_ALPHA_CARBON) {
                snprintf(name, FRAG_NAME_CAP, "alpha carbon");
                continue;
            }
            if (role == FRAG_ROLE_SUGAR) {
                // Ribose has an oxygen on C2' as well as on C3' and C5'
                uint32_t num_o = 0;
                for (size_t k = 0; k < num_core; ++k) {
                    if (c.z[core[k]] != MD_Z_C) continue;
                    for (uint32_t q = c.nb_off[core[k]]; q < c.nb_off[core[k] + 1]; ++q) {
                        const uint32_t j = c.nb[q];
                        num_o += c.z[j] == MD_Z_O && c.unit_of[j] != (int32_t)u;
                    }
                }
                snprintf(name, FRAG_NAME_CAP, "%s", num_o >= 3 ? "ribose" : "deoxyribose");
                continue;
            }
            switch (kind) {
            case MD_FRAGMENT_LIPID_TAIL: {
                uint32_t num_c = 0;
                for (size_t k = 0; k < num_core; ++k) num_c += c.z[core[k]] == MD_Z_C;
                snprintf(name, FRAG_NAME_CAP, "C%u tail", num_c);
                break;
            }
            case MD_FRAGMENT_ATOM: {
                const md_str_t sym = md_atomic_number_symbol((md_atomic_number_t)c.z[first]);
                snprintf(name, FRAG_NAME_CAP, "%.*s", (int)sym.len, sym.ptr);
                break;
            }
            case MD_FRAGMENT_ALKYL: {
                static const char* alkyl[] = { "", "methyl", "ethyl", "propyl", "butyl", "pentyl", "hexyl" };
                if (num_core < MD_COUNTOF(alkyl)) snprintf(name, FRAG_NAME_CAP, "%s", alkyl[num_core]);
                else                              snprintf(name, FRAG_NAME_CAP, "C%zu alkyl", num_core);
                if (num_core == 3 || num_core == 4) {
                    // Branched by the carbons' degrees within it and the one it is attached by: the side chains of
                    // valine (isopropyl), leucine (isobutyl) and isoleucine (sec-butyl)
                    uint32_t attach_deg = 0, max_deg = 0;
                    bool attached = false;
                    for (size_t q = 0; q < num_core; ++q) {
                        const uint32_t a = core[q];
                        uint32_t deg = 0;
                        bool out = false;
                        for (uint32_t e = c.nb_off[a]; e < c.nb_off[a + 1]; ++e) {
                            const uint32_t b = c.nb[e];
                            if (!frag_heavy(&c, b)) continue;
                            if (c.unit_of[b] == (int32_t)u) deg++;
                            else out = true;
                        }
                        max_deg = MD_MAX(max_deg, deg);
                        if (out && !attached) { attached = true; attach_deg = deg; }
                    }
                    if (attached) {
                        const char* branched = NULL;
                        if (num_core == 3 && attach_deg == 2)                    branched = "isopropyl";
                        else if (num_core == 4 && attach_deg == 3)               branched = "tert-butyl";
                        else if (num_core == 4 && max_deg == 3)                  branched = "isobutyl";
                        else if (num_core == 4 && attach_deg == 2)               branched = "sec-butyl";
                        if (branched) snprintf(name, FRAG_NAME_CAP, "%s", branched);
                    }
                }
                break;
            }
            case MD_FRAGMENT_LINKER:
                if (num_core == 1) {
                    snprintf(name, FRAG_NAME_CAP, "methylene");
                } else if (num_core == 2) {
                    // -CH2CH2- between its ends, or -CH(CH3)- with both ends on one carbon (threonine, ibuprofen)
                    uint32_t ends[2] = { 0, 0 };
                    for (size_t q = 0; q < 2; ++q) {
                        for (uint32_t e = c.nb_off[core[q]]; e < c.nb_off[core[q] + 1]; ++e) {
                            const uint32_t j = c.nb[e];
                            if (frag_heavy(&c, j) && c.unit_of[j] != (int32_t)u) ends[q]++;
                        }
                    }
                    snprintf(name, FRAG_NAME_CAP, "%s", (ends[0] == 0 || ends[1] == 0) ? "ethylidene" : "ethylene");
                } else {
                    snprintf(name, FRAG_NAME_CAP, "C%zu linker", num_core);
                }
                break;
            case MD_FRAGMENT_BACKBONE:     // a group of the backbone: a peptide group, a phosphate, a terminus
            case MD_FRAGMENT_GROUP:
            case MD_FRAGMENT_RING_SYSTEM:
                if (num_core == 1) {
                    frag_name_single(name, FRAG_NAME_CAP, &c, core[0]);
                } else if (num_core > 1) {
                    char formula[32];
                    frag_formula(formula, sizeof(formula), &c, core, num_core);
                    bool aromatic = false, carbon_only = true;
                    for (size_t k = 0; k < num_core; ++k) {
                        aromatic |= (frag_aflags(&c, core[k]) & MD_ATOM_FLAG_AROMATIC) != 0;
                        carbon_only &= c.z[core[k]] == MD_Z_C;
                    }
                    // A carbonyl of carbons only, and of a conjugated chain: by its carbon's hydrogens an aldehyde or
                    // a ketone (the library's "carbonyl" and "enone" say less)
                    uint32_t num_o = 0, carbonyl_c = UINT32_MAX;
                    for (size_t k = 0; k < num_core; ++k) {
                        if (c.z[core[k]] == MD_Z_O) num_o++;
                        else if (c.z[core[k]] != MD_Z_C) num_o = 99;
                    }
                    if (num_o == 1) {
                        for (size_t k = 0; k < num_core && carbonyl_c == UINT32_MAX; ++k) {
                            if (c.z[core[k]] != MD_Z_O || frag_heavy_degree(&c, core[k]) != 1) continue;
                            for (uint32_t e = c.nb_off[core[k]]; e < c.nb_off[core[k] + 1]; ++e) {
                                if (c.z[c.nb[e]] == MD_Z_C && md_bond_order(c.nb_flags[e]) == 2) carbonyl_c = c.nb[e];
                            }
                        }
                    }
                    const bool carbonyl = kind == MD_FRAGMENT_GROUP && num_o == 1 && carbonyl_c != UINT32_MAX;
                    const bool aldehyde = carbonyl && frag_hydrogens(&c, carbonyl_c) > 0;
                    if (kind == MD_FRAGMENT_GROUP && carbon_only && num_core >= 4) {
                        snprintf(name, FRAG_NAME_CAP, "polyene");      // conjugated C=C beyond a diene
                    } else if (carbonyl && num_core == 2) {
                        snprintf(name, FRAG_NAME_CAP, "%s", aldehyde ? "aldehyde" : "ketone");
                    } else if (carbonyl && num_core >= 6) {
                        snprintf(name, FRAG_NAME_CAP, "%s", aldehyde ? "polyenal" : "polyenone");
                    } else {
                        snprintf(name, FRAG_NAME_CAP, "%s %s", kind != MD_FRAGMENT_RING_SYSTEM ? "group" : aromatic ? "aryl" : "ring", formula);
                    }
                    if (carbonyl && num_core == 2) break;   // named: not to the library, whose "carbonyl" would win
                    // identified as a whole against the library, below, once per core alike
                    md_array_shrink(core_enc, 0);
                    for (size_t k = 0; k < num_core; ++k) {
                        const uint32_t i = core[k];
                        const size_t a = (size_t)c.atom[i];
                        const uint64_t q = sys->atom.formal_charge ? (uint8_t)(int8_t)md_atom_formal_charge(&sys->atom, a) : 0;
                        const uint64_t ar = (frag_aflags(&c, i) & MD_ATOM_FLAG_AROMATIC) ? 1 : 0;
                        md_array_push(core_enc, (uint64_t)c.z[i] | (q << 8) | ((uint64_t)frag_hydrogens(&c, i) << 16) | (ar << 24), talloc);
                    }
                    for (size_t k = 0; k < num_core; ++k) {
                        const uint32_t i = core[k];
                        for (uint32_t e = c.nb_off[i]; e < c.nb_off[i + 1]; ++e) {
                            const uint32_t j = c.nb[e];
                            if (j <= i || c.core_of[j] != (int32_t)u || !frag_heavy(&c, j)) continue;
                            // j's position among the core, which is in ascending slot order
                            size_t lo = k + 1, hi = num_core;
                            while (lo < hi) { const size_t mid = (lo + hi) / 2; if (core[mid] < j) lo = mid + 1; else hi = mid; }
                            const md_bond_flags_t bf = c.nb_flags[e];
                            const uint64_t b = (uint64_t)md_bond_order(bf) | ((bf & MD_BOND_FLAG_AROMATIC) ? 4u : 0u) | ((bf & MD_BOND_FLAG_DELOCALIZED) ? 8u : 0u);
                            md_array_push(core_enc, ((uint64_t)1 << 63) | ((uint64_t)k << 40) | ((uint64_t)lo << 16) | b, talloc);
                        }
                    }
                    const size_t len = md_array_size(core_enc);
                    uint64_t h = md_hash64(core_enc, len * sizeof(uint64_t), 0);
                    if (h >= MD_HASH_TOMBSTONE) h -= 2;
                    const uint32_t* first_q = md_hashmap_get(&core_by_hash, h);
                    uint32_t qi = UINT32_MAX;
                    for (uint32_t cq = first_q ? *first_q : UINT32_MAX; cq != UINT32_MAX; cq = core_next[cq]) {
                        if (core_off[cq + 1] - core_off[cq] == len && memcmp(core_pool + core_off[cq], core_enc, len * sizeof(uint64_t)) == 0) { qi = cq; break; }
                    }
                    if (qi == UINT32_MAX) {
                        qi = (uint32_t)md_array_size(queries);
                        for (size_t k = 0; k < num_core; ++k) core_global[k] = (uint32_t)c.atom[core[k]];
                        qsort(core_global, num_core, sizeof(uint32_t), frag_cmp_u32);
                        md_array_push(queries, md_bitset_from_indices(core_global, num_core, talloc), talloc);
                        md_array_push(query_frag, f, talloc);
                        md_array_push_array(core_pool, core_enc, len, talloc);
                        md_array_push(core_off, (uint32_t)md_array_size(core_pool), talloc);
                        md_array_push(core_next, first_q ? *first_q : UINT32_MAX, talloc);
                        if (first_q) *(uint32_t*)md_hashmap_get(&core_by_hash, h) = qi;
                        else         md_hashmap_add(&core_by_hash, h, qi);
                    }
                    frag_query[f] = qi;
                }
                break;
            default:
                snprintf(name, FRAG_NAME_CAP, "%s", md_fragment_kind_str(kind));
                break;
            }
        }
        md_array_free(c.unit_kind, salloc);
        md_array_free(c.unit_role, salloc);
    }

    ms_analyse = md_tick_to_milliseconds(md_tick_now() - tick);
    tick = md_tick_now();

    // ---- the library ----
    if (md_array_size(queries) > 0) {
        md_match_library_t* lib = md_match_library_create(talloc);
        for (size_t e = 0; e < MD_COUNTOF(frag_library); ++e) {
            if (md_match_library_add_smiles(lib, md_str_from_cstr(frag_library[e].smiles), md_str_from_cstr(frag_library[e].name), NULL) < 0) {
                MD_LOG_ERROR("Fragments: library entry '%s' does not parse", frag_library[e].smiles);
            }
        }
        md_match_identify_desc_t desc = {
            .library = lib,
            .level = MD_MATCH_LEVEL_MOLECULE,
            .flags = MD_MATCH_FLAG_WHOLE,
            .units = queries,
            .num_units = md_array_size(queries),
        };
        md_match_identify_result_t res = {0};
        if (md_match_identify(&res, &desc, sys, talloc)) {
            for (size_t k = 0; k < res.count; ++k) {
                const uint32_t f = query_frag[res.unit[k]];
                const md_str_t nm = md_match_library_name(lib, res.entry[k]);
                snprintf(names + (size_t)f * FRAG_NAME_CAP, FRAG_NAME_CAP, "%.*s", (int)nm.len, nm.ptr);
            }
            md_match_identify_result_free(&res);
        }
        md_match_library_destroy(lib);
        // and the cores alike named as theirs was
        for (size_t f = 0; f < md_array_size(frag_query); ++f) {
            const uint32_t q = frag_query[f];
            if (q == UINT32_MAX || query_frag[q] == (uint32_t)f) continue;
            memcpy(names + f * FRAG_NAME_CAP, names + (size_t)query_frag[q] * FRAG_NAME_CAP, FRAG_NAME_CAP);
        }
    }
    // The amide of a backbone is the peptide group, apart from the amides of side chains (Asn, Gln)
    for (size_t f = 0; f < md_array_size(kinds); ++f) {
        char* name = names + f * FRAG_NAME_CAP;
        if (kinds[f] == MD_FRAGMENT_BACKBONE && strcmp(name, "amide") == 0) snprintf(name, FRAG_NAME_CAP, "peptide");
    }

    ms_library = md_tick_to_milliseconds(md_tick_now() - tick);
    tick = md_tick_now();

    // ---- types ----
    const size_t num_tfrags = md_array_size(kinds);
    md_array_resize(out->frag_type, num_tfrags, alloc);
    md_hashmap32_t type_of = { .allocator = talloc };
    for (size_t f = 0; f < num_tfrags; ++f) {
        const md_str_t name = md_str_from_cstr(names + f * FRAG_NAME_CAP);
        const md_fragment_kind_t kind = (md_fragment_kind_t)kinds[f];
        const uint64_t key = frag_type_key(name, kind);
        const uint32_t* found = md_hashmap_get(&type_of, key);
        if (found) {
            out->frag_type[f] = *found;
            continue;
        }
        const uint32_t t = (uint32_t)out->type.count++;
        md_array_push(out->type.kind, (uint8_t)kind, alloc);
        md_array_push(out->type.name, md_alloc_str(alloc, name), alloc);
        const uint32_t* old = md_hashmap_get(&old_color, key);
        md_array_push(out->type.color, old ? *old : frag_default_color(name), alloc);
        md_hashmap_add(&type_of, key, t);
        out->frag_type[f] = t;
    }

    // ---- the molecules' fragments ----
    md_array_resize(out->frag_offset, num_mols + 1, alloc);
    out->frag_offset[0] = 0;
    for (size_t m = 0; m < num_mols; ++m) {
        const uint32_t t = out->mol_template[m];
        out->frag_offset[m + 1] = out->frag_offset[m] + (out->template_frag_offset[t + 1] - out->template_frag_offset[t]);
    }
    out->count = out->frag_offset[num_mols];

    MD_LOG_DEBUG("Fragments: %zu molecules, %zu templates, %zu fragments, %zu cores identified, in ms: templates %.1f | analysis %.1f | "
        "library %.1f | types %.1f", num_mols, num_tmpl, out->count, md_array_size(queries), ms_encode, ms_analyse, ms_library,
        md_tick_to_milliseconds(md_tick_now() - tick));
    md_temp_end(temp);
    md_system_topology_changed(sys);
    return true;
}

// ---------------------------------------------------------------------------
// Queries
// ---------------------------------------------------------------------------

size_t md_system_fragment_count(const md_system_t* sys) {
    return sys ? sys->fragment.count : 0;
}

size_t md_system_fragment_type_count(const md_system_t* sys) {
    return sys ? sys->fragment.type.count : 0;
}

// The molecule of a fragment and the fragment's index within its template's
static inline bool frag_locate(const md_system_t* sys, md_fragment_idx_t frag, size_t* mol, uint32_t* tfrag) {
    const md_fragment_data_t* d = &sys->fragment;
    if (frag >= d->count || !d->frag_offset) return false;
    const int32_t m = md_offset_range_find(d->frag_offset, md_array_size(d->mol_template), frag);
    if (m < 0) return false;
    *mol = (size_t)m;
    *tfrag = d->template_frag_offset[d->mol_template[m]] + (frag - d->frag_offset[m]);
    return true;
}

md_fragment_idx_t md_system_fragment_find_by_atom_idx(const md_system_t* sys, size_t atom_idx) {
    if (!sys || !sys->fragment.frag_offset || atom_idx >= sys->atom.count) return MD_FRAGMENT_NONE;
    const md_fragment_data_t* d = &sys->fragment;
    const int32_t m = md_system_molecule_find_by_atom_idx(sys, atom_idx);
    if (m < 0 || (size_t)m >= md_array_size(d->mol_template)) return MD_FRAGMENT_NONE;
    const uint32_t t = d->mol_template[m];
    const uint32_t local = (uint32_t)sys->molecule.atom_slot[atom_idx] - sys->molecule.offset[m];
    const uint32_t f = d->slot_frag[d->template_slot_offset[t] + local];
    return f == MD_FRAGMENT_NONE ? MD_FRAGMENT_NONE : d->frag_offset[m] + f;
}

md_fragment_idx_t md_system_fragment_find_by_ring(const md_system_t* sys, size_t ring_idx) {
    if (!sys || ring_idx >= md_ring_data_count(&sys->ring)) return MD_FRAGMENT_NONE;
    return md_system_fragment_find_by_atom_idx(sys, (size_t)md_ring_atoms(&sys->ring, ring_idx)[0]);
}

size_t md_system_fragment_molecule(const md_system_t* sys, md_fragment_idx_t frag) {
    size_t m = 0; uint32_t tf = 0;
    return frag_locate(sys, frag, &m, &tf) ? m : SIZE_MAX;
}

md_fragment_type_idx_t md_system_fragment_type(const md_system_t* sys, md_fragment_idx_t frag) {
    size_t m = 0; uint32_t tf = 0;
    return frag_locate(sys, frag, &m, &tf) ? sys->fragment.frag_type[tf] : MD_FRAGMENT_NONE;
}

md_fragment_kind_t md_system_fragment_kind(const md_system_t* sys, md_fragment_idx_t frag) {
    return md_fragment_type_kind(sys, md_system_fragment_type(sys, frag));
}

md_str_t md_system_fragment_name(const md_system_t* sys, md_fragment_idx_t frag) {
    return md_fragment_type_name(sys, md_system_fragment_type(sys, frag));
}

size_t md_system_fragment_atom_count(const md_system_t* sys, md_fragment_idx_t frag) {
    size_t m = 0; uint32_t tf = 0;
    if (!frag_locate(sys, frag, &m, &tf)) return 0;
    return sys->fragment.frag_slot_offset[tf + 1] - sys->fragment.frag_slot_offset[tf];
}

md_fragment_atom_iter_t md_fragment_atom_iter(const md_system_t* sys, md_fragment_idx_t frag) {
    md_fragment_atom_iter_t it = {0};
    it.atom = -1;
    size_t m = 0; uint32_t tf = 0;
    if (!frag_locate(sys, frag, &m, &tf)) return it;
    const md_fragment_data_t* d = &sys->fragment;
    it.slot = d->frag_slot + d->frag_slot_offset[tf];
    it.end  = d->frag_slot + d->frag_slot_offset[tf + 1];
    it.mol_atom = sys->molecule.atom_idx + sys->molecule.offset[m];
    return it;
}

bool md_fragment_atom_iter_next(md_fragment_atom_iter_t* it) {
    if (!it || !it->slot || it->slot >= it->end) return false;
    it->atom = it->mol_atom[*it->slot++];
    return true;
}

// Sorts and removes repeats of a small list, returns the new count
// A few distinct values, ascending: inserted into a short list, which is all a residue's fragments or a fragment's
// residues ever are; to a sorted array past it
typedef struct frag_set_t {
    uint32_t  local[32];
    uint32_t* v;
    size_t    n;
    size_t    cap;
    md_array(uint32_t) spill;
    md_allocator_i* alloc;
} frag_set_t;

static void frag_set_add(frag_set_t* s, uint32_t x) {
    size_t lo = 0, hi = s->n;
    while (lo < hi) { const size_t mid = (lo + hi) / 2; if (s->v[mid] < x) lo = mid + 1; else hi = mid; }
    if (lo < s->n && s->v[lo] == x) return;
    if (s->n == s->cap) {
        md_array_resize(s->spill, s->cap * 2, s->alloc);
        if (s->v == s->local) md_memcpy(s->spill, s->local, sizeof(uint32_t) * s->n);
        s->v = s->spill;
        s->cap *= 2;
    }
    md_memmove(s->v + lo + 1, s->v + lo, sizeof(uint32_t) * (s->n - lo));
    s->v[lo] = x;
    s->n++;
}

size_t md_system_residue_fragments(md_fragment_idx_t* out, size_t cap, const md_system_t* sys, size_t res_idx) {
    if (!sys || res_idx >= sys->residue.count) return 0;
    const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
    frag_set_t set = { .cap = 32, .alloc = md_get_heap_allocator() };
    set.v = set.local;
    md_fragment_idx_t last = MD_FRAGMENT_NONE;
    for (uint32_t a = range.beg; a < range.end; ++a) {
        const md_fragment_idx_t x = md_system_fragment_find_by_atom_idx(sys, a);
        if (x != MD_FRAGMENT_NONE && x != last) frag_set_add(&set, x);
        last = x;
    }
    if (out) md_memcpy(out, set.v, sizeof(uint32_t) * MD_MIN(set.n, cap));
    md_array_free(set.spill, set.alloc);
    return set.n;
}

size_t md_system_fragment_residues(md_residue_idx_t* out, size_t cap, const md_system_t* sys, md_fragment_idx_t frag) {
    if (!sys || !sys->residue.count) return 0;
    frag_set_t set = { .cap = 32, .alloc = md_get_heap_allocator() };
    set.v = set.local;
    md_residue_idx_t last = -1;
    for (md_fragment_atom_iter_t it = md_fragment_atom_iter(sys, frag); md_fragment_atom_iter_next(&it);) {
        const md_residue_idx_t x = md_system_residue_find_by_atom_idx(sys, (size_t)it.atom);
        if (x >= 0 && x != last) frag_set_add(&set, (uint32_t)x);
        last = x;
    }
    if (out) md_memcpy(out, set.v, sizeof(uint32_t) * MD_MIN(set.n, cap));
    md_array_free(set.spill, set.alloc);
    return set.n;
}

size_t md_system_fragment_find_by_type(md_fragment_idx_t* out, size_t cap, const md_system_t* sys, md_fragment_type_idx_t type) {
    if (!sys || type >= sys->fragment.type.count) return 0;
    const md_fragment_data_t* d = &sys->fragment;
    size_t n = 0;
    for (size_t m = 0; m < md_array_size(d->mol_template); ++m) {
        const uint32_t t = d->mol_template[m];
        const uint32_t beg = d->template_frag_offset[t], end = d->template_frag_offset[t + 1];
        for (uint32_t f = beg; f < end; ++f) {
            if (d->frag_type[f] != type) continue;
            if (out && n < cap) out[n] = d->frag_offset[m] + (f - beg);
            n++;
        }
    }
    return n;
}

md_fragment_kind_t md_fragment_type_kind(const md_system_t* sys, md_fragment_type_idx_t type) {
    return (sys && type < sys->fragment.type.count) ? (md_fragment_kind_t)sys->fragment.type.kind[type] : MD_FRAGMENT_KIND_COUNT;
}

md_str_t md_fragment_type_name(const md_system_t* sys, md_fragment_type_idx_t type) {
    return (sys && type < sys->fragment.type.count) ? sys->fragment.type.name[type] : (md_str_t){0};
}

uint32_t md_fragment_type_color(const md_system_t* sys, md_fragment_type_idx_t type) {
    return (sys && type < sys->fragment.type.count) ? sys->fragment.type.color[type] : 0xFFFFFFFFu;
}

md_fragment_type_idx_t md_fragment_type_find(const md_system_t* sys, md_str_t name, md_fragment_kind_t kind) {
    if (!sys) return MD_FRAGMENT_NONE;
    for (size_t t = 0; t < sys->fragment.type.count; ++t) {
        if (kind != MD_FRAGMENT_KIND_COUNT && sys->fragment.type.kind[t] != (uint8_t)kind) continue;
        if (md_str_eq(sys->fragment.type.name[t], name)) return (md_fragment_type_idx_t)t;
    }
    return MD_FRAGMENT_NONE;
}
