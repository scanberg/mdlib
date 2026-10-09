#include <md_match.h>

#include <md_smiles.h>
#include <md_system.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_bitfield.h>
#include <core/md_common.h>
#include <core/md_hash.h>
#include <core/md_intrinsics.h>
#include <core/md_log.h>

#include <stdlib.h>
#include <string.h>

// See md_match.h for what a match is and how the search runs. In short:
//   - The system is read through a TARGET: the graph which is searched (THE GRAPH in md_match.h), its molecules,
//     what each molecule can tell and its resonance groups. It is built once per search or identification, in time
//     linear in the size of the system.
//   - The query is turned into a PLAN once per search: the order in which its atoms are matched, for each atom the
//     earlier atom it is reached from (its parent) and the further bonds back into the earlier atoms, its constraints
//     with the hydrogens folded into it, and the hydrogens that are mapped after the rest.
//   - The search walks the graph unit by unit. The only state of the size of the system is a few bit arrays (atoms in
//     the current mapping, atoms taken by earlier matches), set and cleared per step.

#define H_OPEN 255
#define MAX_HYDROGENS_PER_ATOM 16
#define MAX_LEAVES 16

// ### BITS ###

static inline bool bit_test(const uint64_t* bits, size_t i) {
    return (bits[i >> 6] >> (i & 63)) & 1;
}

static inline void bit_set(uint64_t* bits, size_t i) {
    bits[i >> 6] |= 1ULL << (i & 63);
}

static inline void bit_clear(uint64_t* bits, size_t i) {
    bits[i >> 6] &= ~(1ULL << (i & 63));
}

static uint64_t* bits_create(size_t num_bits, md_allocator_i* alloc) {
    const size_t bytes = MAX(1, DIV_UP(num_bits, 64)) * sizeof(uint64_t);
    uint64_t* bits = md_alloc(alloc, bytes);
    MEMSET(bits, 0, bytes);
    return bits;
}

static inline bool label_eq(md_label_t a, md_label_t b) {
    return a.len == b.len && MEMCMP(a.buf, b.buf, a.len) == 0;
}

static int cmp_u64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a;
    const uint64_t y = *(const uint64_t*)b;
    return (x > y) - (x < y);
}

// Index i of the range off[i] <= value < off[i+1] among count ranges (off holds count + 1 entries), -1 if none
static int64_t range_find(const uint32_t* off, size_t count, uint32_t value) {
    if (!off || count == 0 || value < off[0] || value >= off[count]) return -1;
    size_t lo = 0, hi = count;
    while (hi - lo > 1) {
        const size_t mid = (lo + hi) / 2;
        if (off[mid] <= value) lo = mid;
        else hi = mid;
    }
    return (int64_t)lo;
}

// ### TARGET ###
// Read only access to what the search needs of the system, see THE GRAPH and WHAT THE SYSTEM CAN TELL in md_match.h

typedef struct resolution_t {
    uint64_t h_elem[2];     // Elements (z < 128) of which some atom of the molecule has a hydrogen bonded to it
    bool     any_h;         // The molecule has hydrogen atoms: its protonation is known
} resolution_t;

typedef struct target_t {
    const md_system_t* sys;
    size_t num_atoms;

    const md_atom_type_idx_t* type_idx;
    const md_atomic_number_t* type_z;
    const md_label_t*         type_name;
    size_t                    num_types;

    // Chemistry, NULL when the system has none (see md_chem.h)
    const uint8_t*    h_count;      // Hydrogens of each atom, explicit and implicit
    const int8_t*     charge;       // Formal charges
    const md_atom_flags_t* atom_flags;

    // The graph in compressed rows, the system's own where nothing is left out of it. NULL without bonds.
    const uint32_t*        conn_off;    // [num_atoms + 1]
    const md_atom_idx_t*   conn_atom;
    const md_bond_idx_t*   conn_bond;
    const md_bond_flags_t* bond_flags;  // By bond index

    // Atoms outside of the graph (virtual sites), NULL for none
    const uint8_t* skip;

    // Molecules: the connected parts of the graph, numbered in the order of their first atom. Where nothing is left
    // out of the graph they are the structures of the system, read as they are: the atoms of molecule m are
    // mol_atom[mol_off[m] .. mol_off[m + 1]) either way, and the molecule of an atom is mol_of[a], or the structure
    // whose slots hold mol_slot[a].
    size_t          num_mol;
    const uint32_t* mol_off;        // [num_mol + 1]
    const int32_t*  mol_atom;
    const int32_t*  mol_of;         // [num_atoms] NULL where the molecules are the structures
    const int32_t*  mol_slot;       // [num_atoms] Where the molecules are the structures
    resolution_t*   mol_res;        // [num_mol] What each can tell, worked out when first asked (TARGET_RESOLUTION)
    uint8_t*        mol_res_done;   // [num_mol]

    // Resonance groups: atoms joined by aromatic or delocalized bonds. NULL unless built (TARGET_GROUPS).
    int32_t*  group_of;         // [num_atoms] -1 for an atom in none
    int32_t*  group_charge;     // [num_groups] The formal charge of the group
    uint32_t* group_size;
    size_t    num_groups;
} target_t;

enum {
    TARGET_RESOLUTION = 0x1,    // What each molecule can tell (mol_res)
    TARGET_GROUPS     = 0x2,    // Resonance groups
};

static inline md_atomic_number_t target_z(const target_t* tg, int32_t a) {
    if (!tg->type_idx || !tg->type_z) return 0;
    const md_atom_type_idx_t t = tg->type_idx[a];
    return t < tg->num_types ? tg->type_z[t] : 0;
}

static inline md_label_t target_name(const target_t* tg, int32_t a) {
    md_label_t lbl = {0};
    if (tg->type_idx && tg->type_name) {
        const md_atom_type_idx_t t = tg->type_idx[a];
        if (t < tg->num_types) lbl = tg->type_name[t];
    }
    return lbl;
}

static inline uint32_t conn_beg(const target_t* tg, int32_t a) {
    return tg->conn_off ? tg->conn_off[a] : 0;
}

static inline uint32_t conn_end(const target_t* tg, int32_t a) {
    return tg->conn_off ? tg->conn_off[a + 1] : 0;
}

static inline md_bond_flags_t target_bond_flags(const target_t* tg, md_bond_idx_t b) {
    return tg->bond_flags ? tg->bond_flags[b] : MD_BOND_FLAG_NONE;
}

// Hydrogen atoms bonded to the atom
static inline uint32_t count_hydrogens(const target_t* tg, int32_t a) {
    uint32_t count = 0;
    for (uint32_t c = conn_beg(tg, a); c < conn_end(tg, a); ++c) {
        count += target_z(tg, tg->conn_atom[c]) == MD_Z_H;
    }
    return count;
}

// Hydrogens of the atom: the count of the system where it has one (explicit and implicit), the bonded ones otherwise.
// Atoms outside of what md_chem_perceive handles have a count of 0, whatever is bonded to them.
static inline uint32_t target_hydrogens(const target_t* tg, int32_t a) {
    const uint32_t bonded = count_hydrogens(tg, a);
    return tg->h_count ? MAX((uint32_t)tg->h_count[a], bonded) : bonded;
}

static inline int target_charge(const target_t* tg, int32_t a) {
    return tg->charge ? tg->charge[a] : 0;
}

static inline bool bond_order_known(md_bond_flags_t f) {
    return md_bond_order(f) != MD_BOND_ORDER_UNKNOWN || (f & (MD_BOND_FLAG_AROMATIC | MD_BOND_FLAG_DELOCALIZED));
}

// 1 aromatic, 0 not aromatic, -1 unknown (no bond of known order)
static inline int target_aromatic(const target_t* tg, int32_t a) {
    if (tg->atom_flags && (tg->atom_flags[a] & MD_ATOM_FLAG_AROMATIC)) return 1;
    bool known = false;
    for (uint32_t c = conn_beg(tg, a); c < conn_end(tg, a); ++c) {
        const md_bond_flags_t f = target_bond_flags(tg, tg->conn_bond[c]);
        if (f & MD_BOND_FLAG_AROMATIC) return 1;
        known |= bond_order_known(f);
    }
    return known ? 0 : -1;
}

static inline bool target_skipped(const target_t* tg, int32_t a) {
    return tg->skip && tg->skip[a];
}

static inline int64_t target_molecule(const target_t* tg, int32_t a) {
    if (tg->mol_of) return tg->mol_of[a];
    return range_find(tg->mol_off, tg->num_mol, (uint32_t)tg->mol_slot[a]);
}

// Elements (z < 128) of which some atom of the molecule has a hydrogen bonded to it, and whether it has hydrogens
static const resolution_t* molecule_resolution(const target_t* tg, int64_t m) {
    ASSERT(tg->mol_res && 0 <= m && (size_t)m < tg->num_mol);
    resolution_t* res = &tg->mol_res[m];
    if (!tg->mol_res_done[m]) {
        MEMSET(res, 0, sizeof(resolution_t));
        for (uint32_t i = tg->mol_off[m]; i < tg->mol_off[m + 1]; ++i) {
            const int32_t a = tg->mol_atom[i];
            if (target_z(tg, a) != MD_Z_H) continue;
            res->any_h = true;
            for (uint32_t c = conn_beg(tg, a); c < conn_end(tg, a); ++c) {
                const md_atomic_number_t z = target_z(tg, tg->conn_atom[c]);
                if (z != MD_Z_H && z < 128) res->h_elem[z >> 6] |= 1ULL << (z & 63);
            }
        }
        tg->mol_res_done[m] = 1;
    }
    return res;
}

static inline bool resolution_h(const resolution_t* res, md_atomic_number_t z) {
    return z < 128 && ((res->h_elem[z >> 6] >> (z & 63)) & 1);
}

// A coordination bond (a metal to a non-metal) between residues, or any in a system without residues
static bool is_coordination(const md_system_t* sys, md_bond_flags_t f, int32_t a, int32_t b) {
    if (!(f & MD_BOND_FLAG_COORDINATE)) return false;
    if (!(sys->component.count > 0 && sys->component.atom_offset)) return true;
    return range_find(sys->component.atom_offset, sys->component.count, (uint32_t)a) != range_find(sys->component.atom_offset, sys->component.count, (uint32_t)b);
}

static int32_t uf_find(int32_t* parent, int32_t x) {
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

static void uf_union(int32_t* parent, int32_t a, int32_t b) {
    a = uf_find(parent, a);
    b = uf_find(parent, b);
    if (a != b) parent[MAX(a, b)] = MIN(a, b);
}

// what: TARGET_RESOLUTION, TARGET_GROUPS. Everything is allocated from arena.
static void target_build(target_t* tg, const md_system_t* sys, uint32_t what, md_allocator_i* arena) {
    MEMSET(tg, 0, sizeof(target_t));
    const size_t N = sys->atom.count;
    tg->sys        = sys;
    tg->num_atoms  = N;
    tg->type_idx   = sys->atom.type_idx;
    tg->type_z     = sys->atom.type.z;
    tg->type_name  = sys->atom.type.name;
    tg->num_types  = sys->atom.type.count;
    tg->h_count    = sys->atom.hydrogen_count;
    tg->charge     = sys->atom.formal_charge;
    tg->atom_flags = sys->atom.flags;
    if (N == 0) return;

    // Virtual sites are not atoms of the graph
    uint8_t* skip = NULL;
    for (size_t i = 0; i < N; ++i) {
        if (md_atom_particle_kind(&sys->atom, i) == MD_PARTICLE_VIRTUAL_SITE) {
            if (!skip) {
                skip = md_alloc(arena, N);
                MEMSET(skip, 0, N);
            }
            skip[i] = 1;
        }
    }

    const md_bond_conn_data_t* conn = &sys->bond.conn;
    const bool has_bonds = sys->bond.count > 0 && conn->offset && conn->offset_count >= N + 1 && conn->atom_idx && conn->bond_idx;
    uint8_t* drop = NULL;
    if (has_bonds) {
        tg->bond_flags = sys->bond.flags;

        // Bonds left out: coordination, and any bond of a virtual site
        for (size_t b = 0; b < sys->bond.count; ++b) {
            const int32_t i = sys->bond.pairs[b].idx[0];
            const int32_t j = sys->bond.pairs[b].idx[1];
            const md_bond_flags_t f = sys->bond.flags ? sys->bond.flags[b] : MD_BOND_FLAG_NONE;
            if (is_coordination(sys, f, i, j) || (skip && (skip[i] || skip[j]))) {
                if (!drop) {
                    drop = md_alloc(arena, sys->bond.count);
                    MEMSET(drop, 0, sys->bond.count);
                }
                drop[b] = 1;
            }
        }

        if (!drop) {
            tg->conn_off  = conn->offset;
            tg->conn_atom = conn->atom_idx;
            tg->conn_bond = conn->bond_idx;
        } else {
            const uint32_t total = conn->offset[N];
            uint32_t*      off   = md_alloc(arena, sizeof(uint32_t) * (N + 1));
            md_atom_idx_t* atom  = md_alloc(arena, sizeof(md_atom_idx_t) * MAX(1, total));
            md_bond_idx_t* bond  = md_alloc(arena, sizeof(md_bond_idx_t) * MAX(1, total));
            uint32_t n = 0;
            for (size_t a = 0; a < N; ++a) {
                off[a] = n;
                for (uint32_t c = conn->offset[a]; c < conn->offset[a + 1]; ++c) {
                    if (drop[conn->bond_idx[c]]) continue;
                    atom[n] = conn->atom_idx[c];
                    bond[n] = conn->bond_idx[c];
                    n++;
                }
            }
            off[N] = n;
            tg->conn_off  = off;
            tg->conn_atom = atom;
            tg->conn_bond = bond;
        }
    }

    tg->skip = skip;

    // Molecules. The structures of the system are its bond graph, cut into connected parts, with links of their own
    // where atoms have no bonds to hold them (coarse grained beads, virtual sites): without those, and without
    // coordination bonds to leave out, they are the molecules.
    const bool coarse_grained = md_system_is_coarse_grained(sys);
    const md_structure_data_t* st = &sys->structure;
    if (!drop && !skip && !coarse_grained && st->count > 0 && st->offset && st->atom_idx && st->atom_slot && st->offset[st->count] == N) {
        tg->num_mol  = st->count;
        tg->mol_off  = st->offset;
        tg->mol_atom = st->atom_idx;
        tg->mol_slot = st->atom_slot;
    } else {
        // Union find over the bonds of the graph, each part rooted at its first atom, numbered in that order
        int32_t* mol_of = md_alloc(arena, sizeof(int32_t) * N);
        int32_t* root   = md_alloc(arena, sizeof(int32_t) * N);
        for (size_t i = 0; i < N; ++i) root[i] = (int32_t)i;
        if (has_bonds) {
            for (size_t b = 0; b < sys->bond.count; ++b) {
                if (drop && drop[b]) continue;
                uf_union(root, sys->bond.pairs[b].idx[0], sys->bond.pairs[b].idx[1]);
            }
        }
        size_t num_mol = 0;
        for (size_t i = 0; i < N; ++i) {
            if (skip && skip[i]) {
                mol_of[i] = -1;
                continue;
            }
            const int32_t r = uf_find(root, (int32_t)i);
            mol_of[i] = r == (int32_t)i ? (int32_t)num_mol++ : mol_of[r];
        }
        // The atoms of each, ascending
        uint32_t* mol_off = md_alloc(arena, sizeof(uint32_t) * (num_mol + 1));
        uint32_t* fill    = md_alloc(arena, sizeof(uint32_t) * MAX(1, num_mol));
        MEMSET(mol_off, 0, sizeof(uint32_t) * (num_mol + 1));
        for (size_t i = 0; i < N; ++i) {
            if (mol_of[i] >= 0) mol_off[mol_of[i] + 1]++;
        }
        for (size_t m = 0; m < num_mol; ++m) {
            mol_off[m + 1] += mol_off[m];
            fill[m] = mol_off[m];
        }
        int32_t* mol_atom = root;
        for (size_t i = 0; i < N; ++i) {
            if (mol_of[i] >= 0) mol_atom[fill[mol_of[i]]++] = (int32_t)i;
        }
        tg->num_mol  = num_mol;
        tg->mol_off  = mol_off;
        tg->mol_atom = mol_atom;
        tg->mol_of   = mol_of;
    }

    if (what & TARGET_RESOLUTION) {
        tg->mol_res      = md_alloc(arena, sizeof(resolution_t) * MAX(1, tg->num_mol));
        tg->mol_res_done = md_alloc(arena, MAX(1, tg->num_mol));
        MEMSET(tg->mol_res_done, 0, MAX(1, tg->num_mol));
    }

    if ((what & TARGET_GROUPS) && tg->conn_off) {
        int32_t* group_of = md_alloc(arena, sizeof(int32_t) * N);
        int32_t* parent   = md_alloc(arena, sizeof(int32_t) * N);
        for (size_t i = 0; i < N; ++i) {
            group_of[i] = -1;
            parent[i]   = (int32_t)i;
        }
        for (size_t a = 0; a < N; ++a) {
            for (uint32_t c = conn_beg(tg, (int32_t)a); c < conn_end(tg, (int32_t)a); ++c) {
                const int32_t b = tg->conn_atom[c];
                if (b <= (int32_t)a) continue;
                if (!(target_bond_flags(tg, tg->conn_bond[c]) & (MD_BOND_FLAG_AROMATIC | MD_BOND_FLAG_DELOCALIZED))) continue;
                uf_union(parent, (int32_t)a, b);
                group_of[a] = -2;
                group_of[b] = -2;
            }
        }
        // A root is the lowest atom of its group, and the groups are numbered in the order of their first atom
        size_t num_groups = 0;
        for (size_t a = 0; a < N; ++a) {
            if (group_of[a] == -1) continue;
            const int32_t r = uf_find(parent, (int32_t)a);
            group_of[a] = (r == (int32_t)a) ? (int32_t)num_groups++ : group_of[r];
        }
        tg->group_charge = md_alloc(arena, sizeof(int32_t)  * MAX(1, num_groups));
        tg->group_size   = md_alloc(arena, sizeof(uint32_t) * MAX(1, num_groups));
        MEMSET(tg->group_charge, 0, sizeof(int32_t)  * MAX(1, num_groups));
        MEMSET(tg->group_size,   0, sizeof(uint32_t) * MAX(1, num_groups));
        for (size_t a = 0; a < N; ++a) {
            const int32_t g = group_of[a];
            if (g < 0) continue;
            tg->group_size[g]   += 1;
            tg->group_charge[g] += target_charge(tg, (int32_t)a);
        }
        tg->group_of   = group_of;
        tg->num_groups = num_groups;
    }
}

// ### PLAN ###

typedef struct pos_atom_t {
    uint32_t           flags;   // MD_MATCH_ATOM_ELEMENT, _NAME, _CHARGE, _AROMATIC, _ALIPHATIC (HCOUNT is in h_lo, h_hi)
    md_atomic_number_t z;
    bool               is_h;
    bool               h_test;  // h_lo > 0 or h_hi < H_OPEN
    bool               h_stated;// The query gives the count (MD_MATCH_ATOM_HCOUNT), rather than only folded hydrogens
    uint8_t            h_lo;
    uint8_t            h_hi;
    int8_t             charge;
    uint16_t           degree;  // Bonds to searched atoms which are not hydrogens
    md_label_t         name;
} pos_atom_t;

typedef struct plan_t {
    size_t    num_atoms;        // Atoms of the query
    uint32_t  num;              // Atoms searched for, one position each
    bool      impossible;       // Contradicting constraints: nothing can match
    bool      has_charge;       // Some atom states a charge

    uint32_t*   order;          // [num]       Query atom at each position
    int32_t*    parent;         // [num]       Position of the parent, -1 for the first atom of a part
    uint8_t*    parent_order;   // [num]       Order of the bond to the parent
    uint32_t*   back_off;       // [num + 1]   Bonds back to earlier positions other than the parent
    uint32_t*   back_pos;
    uint8_t*    back_order;
    pos_atom_t* atom;           // [num]

    // Hydrogens mapped after the search, grouped by the position of the atom carrying them
    uint32_t    num_folded;
    uint32_t*   folded_atom;    // Query atom
    uint32_t*   folded_host;    // Position of the atom it is bonded to
    pos_atom_t* folded;         // Its own constraints

    // Positions [0, num_core) are searched atom by atom, the leaves [num_core, num) after them in groups: group g
    // is [leaf_beg[g], leaf_beg[g + 1]), the leaves which hang from the atom at position leaf_host[g]
    uint32_t  num_core;
    uint32_t  num_leaf_groups;
    uint32_t* leaf_beg;
    uint32_t* leaf_host;

    // Symmetry breaking (see plan_symmetry), outside of ALL: the atom at position k is below (sym_less) or above the
    // atom at position sym_other[i], for i in [sym_off[k], sym_off[k + 1]). NULL without conditions.
    uint32_t* sym_off;
    uint32_t* sym_other;
    bool*     sym_less;
    uint32_t  num_sym;

    // For MD_MATCH_FLAG_WHOLE: searched atoms which are not hydrogens and the bonds between them
    uint32_t  num_heavy;
    uint32_t  num_heavy_bonds;
} plan_t;

static pos_atom_t make_pos_atom(const md_match_atom_t* a, uint32_t num_folded) {
    pos_atom_t p = {0};
    p.flags  = a->flags & (MD_MATCH_ATOM_ELEMENT | MD_MATCH_ATOM_NAME | MD_MATCH_ATOM_CHARGE | MD_MATCH_ATOM_AROMATIC | MD_MATCH_ATOM_ALIPHATIC);
    p.z      = a->z;
    p.charge = a->charge;
    p.name   = a->name;
    p.is_h   = (a->flags & MD_MATCH_ATOM_ELEMENT) && a->z == MD_Z_H;
    p.h_stated = (a->flags & MD_MATCH_ATOM_HCOUNT) != 0;
    uint32_t lo = p.h_stated ? a->h_min : 0;
    uint32_t hi = p.h_stated ? a->h_max : H_OPEN;
    // The count includes the hydrogens given as atoms (see md_match_atom_t)
    lo = MAX(lo, num_folded);
    p.h_lo   = (uint8_t)MIN(lo, 255);
    p.h_hi   = (uint8_t)hi;
    p.h_test = lo > 0 || hi < H_OPEN;
    return p;
}

static bool query_validate(const md_match_query_t* q) {
    if (!q) {
        MD_LOG_ERROR("md_match: no query");
        return false;
    }
    if (q->num_atoms == 0 || !q->atoms) {
        MD_LOG_ERROR("md_match: the query has no atoms");
        return false;
    }
    if (q->num_atoms >= INT32_MAX) {
        MD_LOG_ERROR("md_match: the query has too many atoms");
        return false;
    }
    if (q->num_bonds > 0 && !q->bonds) {
        MD_LOG_ERROR("md_match: the query has bonds but no bond array");
        return false;
    }
    for (size_t i = 0; i < q->num_atoms; ++i) {
        const md_match_atom_t* a = &q->atoms[i];
        if ((a->flags & MD_MATCH_ATOM_HCOUNT) && a->h_min > a->h_max) {
            MD_LOG_ERROR("md_match: query atom %zu has h_min > h_max", i);
            return false;
        }
    }
    for (size_t i = 0; i < q->num_bonds; ++i) {
        const md_match_bond_t* b = &q->bonds[i];
        if (b->a >= q->num_atoms || b->b >= q->num_atoms) {
            MD_LOG_ERROR("md_match: query bond %zu refers to an atom out of range", i);
            return false;
        }
        if (b->a == b->b) {
            MD_LOG_ERROR("md_match: query bond %zu bonds atom %u to itself", i, b->a);
            return false;
        }
        if (b->order > MD_MATCH_BOND_AROMATIC) {
            MD_LOG_ERROR("md_match: query bond %zu has an unknown order", i);
            return false;
        }
    }
    if (q->num_bonds > 1) {
        md_temp_scope_t temp = md_temp_begin();
        uint64_t* keys = md_temp_alloc_array(temp, uint64_t, q->num_bonds);
        for (size_t i = 0; i < q->num_bonds; ++i) {
            const uint64_t a = MIN(q->bonds[i].a, q->bonds[i].b);
            const uint64_t b = MAX(q->bonds[i].a, q->bonds[i].b);
            keys[i] = (a << 32) | b;
        }
        qsort(keys, q->num_bonds, sizeof(uint64_t), cmp_u64);
        bool dup = false;
        for (size_t i = 1; i < q->num_bonds; ++i) {
            if (keys[i] == keys[i - 1]) {
                MD_LOG_ERROR("md_match: the query bonds atoms %u and %u twice", (uint32_t)(keys[i] >> 32), (uint32_t)keys[i]);
                dup = true;
                break;
            }
        }
        md_temp_end(temp);
        if (dup) return false;
    }
    return true;
}

static bool query_states_charge(const md_match_query_t* q) {
    for (size_t i = 0; i < q->num_atoms; ++i) {
        if (q->atoms[i].flags & MD_MATCH_ATOM_CHARGE) return true;
    }
    return false;
}

// ### SYMMETRY ###
// A query with symmetries matches the same atoms once for each of them: a phenyl ring flipped over is the same ring,
// and a chain with 30 of them would be found 2^30 times over before UNIQUE saw that it is one match. Outside of ALL the
// search only takes the mappings which satisfy a set of conditions map[v] < map[u], which every set of atoms satisfies
// in one of its symmetric mappings at least (Grochow and Kellis 2007): map[v] < map[u] for an atom v and every u a
// symmetry of the query takes v to, then the same among the symmetries which leave v where it is, and so on. The atom
// v chosen is the first in query order, so that a reference given in the order of its atoms still maps onto itself.
// The symmetries are those of the core of the plan (the leaves are mapped as sets already), with every atom labelled
// by all that the search tests of it, its leaves and hydrogens included. They are found by partition refinement and
// a search for each candidate, as nauty does, and each is checked before it is used. Where finding them takes too long
// the search goes on with the conditions found so far, which is slower, never wrong.

typedef struct sym_part_t {
    uint32_t* elem;         // The atoms (core positions), cell by cell
    uint32_t* pos;          // Where each atom is in elem
    uint32_t* cell;         // The start of the cell of each atom
    uint32_t* end;          // For the start of a cell, its end
    uint32_t  num_cells;
} sym_part_t;

typedef struct sym_count_t {
    uint64_t count;
    uint32_t atom;
} sym_count_t;

typedef struct sym_t {
    uint32_t  n;
    uint32_t* off;          // [n + 1] The graph of the core, bond orders in ord
    uint32_t* nbr;
    uint8_t*  ord;
    uint32_t* label;        // [n]

    // Scratch of sym_refine
    uint64_t*    count;     // [n] Bonds into the cell split by, by order
    uint32_t*    touched;   // [n]
    uint32_t*    cells;     // [n]
    uint8_t*     cell_mark; // [n]
    uint32_t*    runs;      // [n + 1]
    sym_count_t* sorted;    // [n]
    uint32_t*    stack;     // [n] Cells to split by
    uint8_t*     in_stack;  // [n]
    uint32_t     num_stack;

    size_t work;
    size_t budget;
    md_allocator_i* arena;
} sym_t;

static void sym_part_alloc(sym_part_t* p, uint32_t n, md_allocator_i* arena) {
    p->elem = md_alloc(arena, sizeof(uint32_t) * n);
    p->pos  = md_alloc(arena, sizeof(uint32_t) * n);
    p->cell = md_alloc(arena, sizeof(uint32_t) * n);
    p->end  = md_alloc(arena, sizeof(uint32_t) * n);
    p->num_cells = 0;
}

static void sym_part_copy(sym_part_t* dst, const sym_part_t* src, uint32_t n) {
    MEMCPY(dst->elem, src->elem, sizeof(uint32_t) * n);
    MEMCPY(dst->pos,  src->pos,  sizeof(uint32_t) * n);
    MEMCPY(dst->cell, src->cell, sizeof(uint32_t) * n);
    MEMCPY(dst->end,  src->end,  sizeof(uint32_t) * n);
    dst->num_cells = src->num_cells;
}

static inline void sym_push(sym_t* y, uint32_t cell) {
    if (y->in_stack[cell]) return;
    y->in_stack[cell] = 1;
    y->stack[y->num_stack++] = cell;
}

static int cmp_u32(const void* a, const void* b) {
    const uint32_t x = *(const uint32_t*)a;
    const uint32_t y = *(const uint32_t*)b;
    return (x > y) - (x < y);
}

static int cmp_sym_count(const void* a, const void* b) {
    const sym_count_t* x = (const sym_count_t*)a;
    const sym_count_t* y = (const sym_count_t*)b;
    if (x->count != y->count) return x->count < y->count ? -1 : 1;
    return (x->atom > y->atom) - (x->atom < y->atom);
}

// Splits the cells of p until every atom of a cell has as many bonds of each order into each cell (an equitable
// partition), by the cells on the stack. Cells split in an order which depends on the graph alone, so that partitions
// which correspond refine to partitions which correspond. Returns false when out of budget.
static bool sym_refine(sym_t* y, sym_part_t* p) {
    while (y->num_stack > 0) {
        if (y->work > y->budget) {
            while (y->num_stack > 0) y->in_stack[y->stack[--y->num_stack]] = 0;
            return false;
        }
        const uint32_t w = y->stack[--y->num_stack];
        y->in_stack[w] = 0;

        // Bonds of each atom into cell w, counted by order in 10 bits each
        uint32_t num_touched = 0;
        for (uint32_t i = w; i < p->end[w]; ++i) {
            const uint32_t v = p->elem[i];
            for (uint32_t e = y->off[v]; e < y->off[v + 1]; ++e) {
                const uint32_t u = y->nbr[e];
                if (y->count[u] == 0) y->touched[num_touched++] = u;
                y->count[u] += 1ULL << (10 * y->ord[e]);
            }
            y->work += 1 + y->off[v + 1] - y->off[v];
        }

        uint32_t num_cells = 0;
        for (uint32_t i = 0; i < num_touched; ++i) {
            const uint32_t c = p->cell[y->touched[i]];
            if (!y->cell_mark[c]) {
                y->cell_mark[c] = 1;
                y->cells[num_cells++] = c;
            }
        }
        qsort(y->cells, num_cells, sizeof(uint32_t), cmp_u32);

        for (uint32_t ci = 0; ci < num_cells; ++ci) {
            const uint32_t c = y->cells[ci];
            const uint32_t e = p->end[c];
            y->cell_mark[c] = 0;
            if (e - c == 1) continue;
            for (uint32_t i = c; i < e; ++i) {
                y->sorted[i - c] = (sym_count_t){.count = y->count[p->elem[i]], .atom = p->elem[i]};
            }
            qsort(y->sorted, e - c, sizeof(sym_count_t), cmp_sym_count);
            y->work += e - c;
            if (y->sorted[0].count == y->sorted[e - c - 1].count) continue;

            // Split into runs of equal counts, in ascending order of the counts
            uint32_t num_runs = 0;
            for (uint32_t i = c; i < e; ++i) {
                const sym_count_t* x = &y->sorted[i - c];
                p->elem[i] = x->atom;
                p->pos[x->atom] = i;
                if (i == c || x->count != y->sorted[i - c - 1].count) y->runs[num_runs++] = i;
            }
            y->runs[num_runs] = e;
            uint32_t largest = 0;
            for (uint32_t r = 0; r < num_runs; ++r) {
                const uint32_t rb = y->runs[r], re = y->runs[r + 1];
                p->end[rb] = re;
                for (uint32_t i = rb; i < re; ++i) p->cell[p->elem[i]] = rb;
                if (re - rb > y->runs[largest + 1] - y->runs[largest]) largest = r;
            }
            p->num_cells += num_runs - 1;
            // Each of them splits further, but the largest only by what the others imply, unless it was to anyway
            const bool was_in = y->in_stack[c] != 0;
            for (uint32_t r = 0; r < num_runs; ++r) {
                if (was_in || r != largest) sym_push(y, y->runs[r]);
            }
        }
        for (uint32_t i = 0; i < num_touched; ++i) y->count[y->touched[i]] = 0;
    }
    return true;
}

// Makes v a cell of its own, at the start of its cell
static void sym_individualize(sym_t* y, sym_part_t* p, uint32_t v) {
    const uint32_t c = p->cell[v];
    const uint32_t e = p->end[c];
    if (e - c == 1) return;
    const uint32_t i = p->pos[v];
    const uint32_t first = p->elem[c];
    p->elem[i] = first;
    p->pos[first] = i;
    p->elem[c] = v;
    p->pos[v] = c;
    p->end[c] = c + 1;
    p->end[c + 1] = e;
    for (uint32_t j = c + 1; j < e; ++j) p->cell[p->elem[j]] = c + 1;
    p->num_cells += 1;
    if (y->in_stack[c]) sym_push(y, c + 1);
    else sym_push(y, c);
}

// Is gamma a symmetry: does it keep the labels, and take bonds to bonds of the same order
static bool sym_verify(const sym_t* y, const uint32_t* gamma) {
    for (uint32_t v = 0; v < y->n; ++v) {
        const uint32_t g = gamma[v];
        if (y->label[v] != y->label[g]) return false;
        if (y->off[v + 1] - y->off[v] != y->off[g + 1] - y->off[g]) return false;
        for (uint32_t e = y->off[v]; e < y->off[v + 1]; ++e) {
            const uint32_t gu = gamma[y->nbr[e]];
            bool found = false;
            for (uint32_t f = y->off[g]; f < y->off[g + 1] && !found; ++f) {
                found = y->nbr[f] == gu && y->ord[f] == y->ord[e];
            }
            if (!found) return false;
        }
    }
    return true;
}

// A symmetry which takes partition a onto partition b, both refined: 1 found (in gamma), 0 none, -1 out of budget
static int sym_search(sym_t* y, const sym_part_t* a, const sym_part_t* b, uint32_t* gamma) {
    const uint32_t n = y->n;
    if (y->work > y->budget) return -1;
    if (a->num_cells != b->num_cells) return 0;
    for (uint32_t i = 0; i < n; ++i) {
        if (a->cell[a->elem[i]] != b->cell[b->elem[i]]) return 0;
    }
    y->work += n;
    if (a->num_cells == n) {
        for (uint32_t i = 0; i < n; ++i) gamma[a->elem[i]] = b->elem[i];
        y->work += n;
        return sym_verify(y, gamma) ? 1 : 0;
    }

    // The first cell of more than one atom: its first atom in a, onto each of the cell in b
    uint32_t c = 0;
    while (a->end[c] - c == 1) c = a->end[c];
    md_vm_arena_temp_t temp = md_vm_arena_temp_begin(y->arena);
    sym_part_t a2, b2;
    sym_part_alloc(&a2, n, y->arena);
    sym_part_alloc(&b2, n, y->arena);
    sym_part_copy(&a2, a, n);
    sym_individualize(y, &a2, a->elem[c]);
    int result = sym_refine(y, &a2) ? 0 : -1;
    for (uint32_t j = c; j < b->end[c] && result == 0; ++j) {
        sym_part_copy(&b2, b, n);
        sym_individualize(y, &b2, b->elem[j]);
        result = sym_refine(y, &b2) ? sym_search(y, &a2, &b2, gamma) : -1;
    }
    md_vm_arena_temp_end(temp);
    return result;
}

typedef struct sym_key_t {
    const uint8_t* ptr;
    uint32_t len;
    uint32_t idx;
} sym_key_t;

static int cmp_sym_key_bytes(const sym_key_t* x, const sym_key_t* y) {
    if (x->len != y->len) return x->len < y->len ? -1 : 1;
    return MEMCMP(x->ptr, y->ptr, x->len);
}

static int cmp_sym_key(const void* a, const void* b) {
    const sym_key_t* x = (const sym_key_t*)a;
    const sym_key_t* y = (const sym_key_t*)b;
    const int c = cmp_sym_key_bytes(x, y);
    if (c) return c;
    return (x->idx > y->idx) - (x->idx < y->idx);
}

// What the search tests of an atom, as bytes
static void sym_key_atom(md_array(uint8_t)* key, const pos_atom_t* a, md_allocator_i* arena) {
    uint8_t b[20] = {0};
    MEMCPY(b, &a->flags, 4);
    b[4] = a->z;
    b[5] = a->is_h;
    b[6] = a->h_lo;
    b[7] = a->h_hi;
    b[8] = a->h_stated;
    b[9] = (uint8_t)a->charge;
    MEMCPY(b + 10, &a->degree, 2);
    MEMCPY(b + 12, a->name.buf, 7);
    b[19] = a->name.len;
    md_array_push_array(*key, b, 20, arena);
}

// Appends the keys given, sorted, after their number
static void sym_key_append_sorted(md_array(uint8_t)* key, sym_key_t* parts, uint32_t num, md_allocator_i* arena) {
    qsort(parts, num, sizeof(sym_key_t), cmp_sym_key);
    md_array_push_array(*key, (const uint8_t*)&num, sizeof(num), arena);
    for (uint32_t i = 0; i < num; ++i) {
        md_array_push_array(*key, (const uint8_t*)&parts[i].len, sizeof(parts[i].len), arena);
        md_array_push_array(*key, parts[i].ptr, parts[i].len, arena);
    }
}

static uint32_t log2_floor(uint32_t x) {
    uint32_t r = 0;
    while (x >>= 1) ++r;
    return r;
}

// pos_of: the position of each query atom, -1 for a folded hydrogen
static void plan_symmetry(plan_t* plan, const md_match_query_t* q, const int32_t* pos_of, md_allocator_i* arena) {
    const uint32_t n = plan->num_core;
    if (n < 2) return;

    sym_t y = {0};
    y.n = n;
    y.arena = arena;

    // The graph of the core
    y.off = md_alloc(arena, sizeof(uint32_t) * (n + 1));
    MEMSET(y.off, 0, sizeof(uint32_t) * (n + 1));
    for (size_t i = 0; i < q->num_bonds; ++i) {
        const int32_t a = pos_of[q->bonds[i].a], b = pos_of[q->bonds[i].b];
        if (a < 0 || b < 0 || (uint32_t)a >= n || (uint32_t)b >= n) continue;
        y.off[a + 1]++;
        y.off[b + 1]++;
    }
    for (uint32_t i = 0; i < n; ++i) y.off[i + 1] += y.off[i];
    const uint32_t num_edges = y.off[n];
    y.nbr = md_alloc(arena, sizeof(uint32_t) * MAX(1, num_edges));
    y.ord = md_alloc(arena, sizeof(uint8_t)  * MAX(1, num_edges));
    {
        uint32_t* fill = md_alloc(arena, sizeof(uint32_t) * n);
        MEMCPY(fill, y.off, sizeof(uint32_t) * n);
        for (size_t i = 0; i < q->num_bonds; ++i) {
            const int32_t a = pos_of[q->bonds[i].a], b = pos_of[q->bonds[i].b];
            if (a < 0 || b < 0 || (uint32_t)a >= n || (uint32_t)b >= n) continue;
            y.nbr[fill[a]] = (uint32_t)b; y.ord[fill[a]++] = (uint8_t)q->bonds[i].order;
            y.nbr[fill[b]] = (uint32_t)a; y.ord[fill[b]++] = (uint8_t)q->bonds[i].order;
        }
    }

    // The label of an atom: itself, its hydrogens, and its leaves with their hydrogens
    md_array(uint8_t)* key_of = md_alloc(arena, sizeof(md_array(uint8_t)) * plan->num);
    for (uint32_t p = 0; p < plan->num; ++p) {
        key_of[p] = 0;
        sym_key_atom(&key_of[p], &plan->atom[p], arena);
    }
    {
        md_array(uint8_t) h_keys = 0;
        for (uint32_t f = 0; f < plan->num_folded; ++f) sym_key_atom(&h_keys, &plan->folded[f], arena);
        sym_key_t* parts = md_alloc(arena, sizeof(sym_key_t) * MAX(1, MAX(plan->num_folded, plan->num - n)));
        uint32_t f = 0;
        for (uint32_t p = 0; p < plan->num; ++p) {
            // The folded hydrogens are grouped by the position of their atom, ascending
            while (f < plan->num_folded && plan->folded_host[f] < p) ++f;
            uint32_t num = 0;
            for (uint32_t g = f; g < plan->num_folded && plan->folded_host[g] == p; ++g) {
                parts[num++] = (sym_key_t){.ptr = h_keys + 20 * g, .len = 20, .idx = 0};
            }
            sym_key_append_sorted(&key_of[p], parts, num, arena);
        }
        bool* has_leaves = md_alloc(arena, sizeof(bool) * n);
        MEMSET(has_leaves, 0, sizeof(bool) * n);
        for (uint32_t g = 0; g < plan->num_leaf_groups; ++g) {
            uint32_t num = 0;
            for (uint32_t p = plan->leaf_beg[g]; p < plan->leaf_beg[g + 1]; ++p) {
                md_array_push(key_of[p], plan->parent_order[p], arena);
                parts[num++] = (sym_key_t){.ptr = key_of[p], .len = (uint32_t)md_array_size(key_of[p]), .idx = 0};
            }
            sym_key_append_sorted(&key_of[plan->leaf_host[g]], parts, num, arena);
            has_leaves[plan->leaf_host[g]] = true;
        }
        for (uint32_t p = 0; p < n; ++p) {
            if (!has_leaves[p]) sym_key_append_sorted(&key_of[p], parts, 0, arena);
        }
    }

    // The first partition: by label, in the order of the labels
    sym_key_t* order = md_alloc(arena, sizeof(sym_key_t) * n);
    for (uint32_t p = 0; p < n; ++p) order[p] = (sym_key_t){.ptr = key_of[p], .len = (uint32_t)md_array_size(key_of[p]), .idx = p};
    qsort(order, n, sizeof(sym_key_t), cmp_sym_key);

    y.label     = md_alloc(arena, sizeof(uint32_t) * n);
    y.count     = md_alloc(arena, sizeof(uint64_t) * n);
    y.touched   = md_alloc(arena, sizeof(uint32_t) * n);
    y.cells     = md_alloc(arena, sizeof(uint32_t) * n);
    y.cell_mark = md_alloc(arena, sizeof(uint8_t)  * n);
    y.runs      = md_alloc(arena, sizeof(uint32_t) * (n + 1));
    y.sorted    = md_alloc(arena, sizeof(sym_count_t) * n);
    y.stack     = md_alloc(arena, sizeof(uint32_t) * n);
    y.in_stack  = md_alloc(arena, sizeof(uint8_t)  * n);
    MEMSET(y.count, 0, sizeof(uint64_t) * n);
    MEMSET(y.cell_mark, 0, n);
    MEMSET(y.in_stack, 0, n);
    // Generous for molecules: refinement is close to linear in the size of the graph, and a search for a symmetry
    // mostly takes one branch
    y.budget = 200000 + (size_t)64 * (n + num_edges) * (1 + log2_floor(n));

    sym_part_t cur;
    sym_part_alloc(&cur, n, arena);
    uint32_t label = 0;
    for (uint32_t i = 0; i < n; ++i) {
        if (i > 0 && cmp_sym_key_bytes(&order[i], &order[i - 1]) != 0) label++;
        y.label[order[i].idx] = label;
        cur.elem[i] = order[i].idx;
        cur.pos[order[i].idx] = i;
    }
    for (uint32_t i = 0; i < n;) {
        uint32_t e = i + 1;
        while (e < n && y.label[cur.elem[e]] == y.label[cur.elem[i]]) ++e;
        cur.end[i] = e;
        for (uint32_t j = i; j < e; ++j) cur.cell[cur.elem[j]] = i;
        cur.num_cells++;
        sym_push(&y, i);
        i = e;
    }

    md_array(uint32_t) cond_v = 0;      // Conditions map[cond_v] < map[cond_u], by position
    md_array(uint32_t) cond_u = 0;
    if (sym_refine(&y, &cur)) {
        md_array(uint32_t*) gens = 0;   // Symmetries found which leave the atoms fixed so far where they are
        int32_t*  orbit = md_alloc(arena, sizeof(int32_t) * n);
        uint32_t* gamma = md_alloc(arena, sizeof(uint32_t) * n);
        sym_part_t a, b;
        sym_part_alloc(&a, n, arena);
        sym_part_alloc(&b, n, arena);
        while (cur.num_cells < n) {
            // v: the first atom in query order of those in cells of more than one
            uint32_t v = UINT32_MAX;
            for (uint32_t i = 0; i < n; i = cur.end[i]) {
                if (cur.end[i] - i == 1) continue;
                for (uint32_t j = i; j < cur.end[i]; ++j) {
                    const uint32_t x = cur.elem[j];
                    if (v == UINT32_MAX || plan->order[x] < plan->order[v]) v = x;
                }
            }
            for (uint32_t i = 0; i < n; ++i) orbit[i] = (int32_t)i;
            for (size_t k = 0; k < md_array_size(gens); ++k) {
                for (uint32_t i = 0; i < n; ++i) uf_union(orbit, (int32_t)i, (int32_t)gens[k][i]);
            }
            // The atoms of its cell which a symmetry takes it to: known by those found, or found now
            bool out_of_budget = false;
            const uint32_t c = cur.cell[v];
            for (uint32_t j = c; j < cur.end[c]; ++j) {
                const uint32_t u = cur.elem[j];
                if (u == v) continue;
                if (uf_find(orbit, (int32_t)u) != uf_find(orbit, (int32_t)v)) {
                    sym_part_copy(&a, &cur, n);
                    sym_individualize(&y, &a, v);
                    int found = sym_refine(&y, &a) ? 0 : -1;
                    if (found == 0) {
                        sym_part_copy(&b, &cur, n);
                        sym_individualize(&y, &b, u);
                        found = sym_refine(&y, &b) ? sym_search(&y, &a, &b, gamma) : -1;
                    }
                    if (found < 0) {
                        out_of_budget = true;
                        break;
                    }
                    if (found > 0) {
                        uint32_t* gen = md_alloc(arena, sizeof(uint32_t) * n);
                        MEMCPY(gen, gamma, sizeof(uint32_t) * n);
                        md_array_push(gens, gen, arena);
                        for (uint32_t i = 0; i < n; ++i) uf_union(orbit, (int32_t)i, (int32_t)gen[i]);
                    }
                }
                if (uf_find(orbit, (int32_t)u) == uf_find(orbit, (int32_t)v)) {
                    md_array_push(cond_v, v, arena);
                    md_array_push(cond_u, u, arena);
                }
            }
            if (out_of_budget) break;
            sym_individualize(&y, &cur, v);
            if (!sym_refine(&y, &cur)) break;
            size_t num_gens = 0;
            for (size_t k = 0; k < md_array_size(gens); ++k) {
                if (gens[k][v] == v) gens[num_gens++] = gens[k];
            }
            md_array_shrink(gens, num_gens);
        }
    }

    // The conditions, by the later of their two positions
    const uint32_t num_cond = (uint32_t)md_array_size(cond_v);
    if (num_cond == 0) return;
    plan->num_sym   = num_cond;
    plan->sym_off   = md_alloc(arena, sizeof(uint32_t) * (plan->num + 1));
    plan->sym_other = md_alloc(arena, sizeof(uint32_t) * num_cond);
    plan->sym_less  = md_alloc(arena, sizeof(bool) * num_cond);
    MEMSET(plan->sym_off, 0, sizeof(uint32_t) * (plan->num + 1));
    for (uint32_t i = 0; i < num_cond; ++i) plan->sym_off[MAX(cond_v[i], cond_u[i]) + 1]++;
    for (uint32_t p = 0; p < plan->num; ++p) plan->sym_off[p + 1] += plan->sym_off[p];
    uint32_t* fill = md_alloc(arena, sizeof(uint32_t) * plan->num);
    MEMCPY(fill, plan->sym_off, sizeof(uint32_t) * plan->num);
    for (uint32_t i = 0; i < num_cond; ++i) {
        const uint32_t v = cond_v[i], u = cond_u[i];
        const uint32_t k = MAX(v, u);
        plan->sym_other[fill[k]] = MIN(v, u);
        plan->sym_less[fill[k]++] = k == v;
    }
}

// rarity[z]: candidate atoms of the system with that element, num_candidates: all of them. symmetry: work out the
// conditions which break the symmetries of the query (see plan_symmetry), which all modes but ALL use.
static void plan_build(plan_t* plan, const md_match_query_t* q, const uint32_t* rarity, uint32_t num_candidates, bool symmetry, md_allocator_i* arena) {
    MEMSET(plan, 0, sizeof(plan_t));
    const uint32_t n = (uint32_t)q->num_atoms;
    const uint32_t m = (uint32_t)q->num_bonds;
    plan->num_atoms  = n;
    plan->has_charge = query_states_charge(q);

    // Adjacency of the whole query in compressed rows
    uint32_t* off   = md_alloc(arena, sizeof(uint32_t) * (n + 1));
    uint32_t* nbr   = md_alloc(arena, sizeof(uint32_t) * MAX(1, 2 * m));
    uint8_t*  nbr_o = md_alloc(arena, sizeof(uint8_t)  * MAX(1, 2 * m));
    MEMSET(off, 0, sizeof(uint32_t) * (n + 1));
    for (uint32_t i = 0; i < m; ++i) {
        off[q->bonds[i].a + 1]++;
        off[q->bonds[i].b + 1]++;
    }
    for (uint32_t i = 0; i < n; ++i) off[i + 1] += off[i];
    {
        uint32_t* fill = md_alloc(arena, sizeof(uint32_t) * n);
        MEMCPY(fill, off, sizeof(uint32_t) * n);
        for (uint32_t i = 0; i < m; ++i) {
            const uint32_t a = q->bonds[i].a, b = q->bonds[i].b;
            nbr[fill[a]] = b; nbr_o[fill[a]++] = (uint8_t)q->bonds[i].order;
            nbr[fill[b]] = a; nbr_o[fill[b]++] = (uint8_t)q->bonds[i].order;
        }
    }

    // Hydrogens with exactly one bond, to an atom which is not a hydrogen, are folded into that atom
    bool*     is_h        = md_alloc(arena, sizeof(bool) * n);
    int32_t*  host        = md_alloc(arena, sizeof(int32_t) * n);
    uint32_t* num_folded  = md_alloc(arena, sizeof(uint32_t) * n);
    for (uint32_t i = 0; i < n; ++i) {
        is_h[i] = (q->atoms[i].flags & MD_MATCH_ATOM_ELEMENT) && q->atoms[i].z == MD_Z_H;
        num_folded[i] = 0;
    }
    for (uint32_t i = 0; i < n; ++i) {
        host[i] = -1;
        if (is_h[i] && off[i + 1] - off[i] == 1 && !is_h[nbr[off[i]]]) {
            host[i] = (int32_t)nbr[off[i]];
            num_folded[host[i]]++;
            plan->num_folded++;
        }
    }

    const uint32_t num_search = n - plan->num_folded;
    plan->num          = num_search;
    plan->order        = md_alloc(arena, sizeof(uint32_t)   * num_search);
    plan->parent       = md_alloc(arena, sizeof(int32_t)    * num_search);
    plan->parent_order = md_alloc(arena, sizeof(uint8_t)    * num_search);
    plan->back_off     = md_alloc(arena, sizeof(uint32_t)   * (num_search + 1));
    plan->back_pos     = md_alloc(arena, sizeof(uint32_t)   * MAX(1, m));
    plan->back_order   = md_alloc(arena, sizeof(uint8_t)    * MAX(1, m));
    plan->atom         = md_alloc(arena, sizeof(pos_atom_t) * num_search);

    // Number of bonds to searched atoms, and the rarity of each atom among the candidates of the system
    uint32_t* degree = md_alloc(arena, sizeof(uint32_t) * n);
    uint32_t* rank   = md_alloc(arena, sizeof(uint32_t) * n);
    for (uint32_t i = 0; i < n; ++i) {
        degree[i] = 0;
        for (uint32_t e = off[i]; e < off[i + 1]; ++e) {
            degree[i] += host[nbr[e]] < 0;
        }
        rank[i] = (q->atoms[i].flags & MD_MATCH_ATOM_ELEMENT) ? rarity[q->atoms[i].z] : num_candidates;
    }

    int32_t*  pos_of   = md_alloc(arena, sizeof(int32_t)  * n);   // -1 until ordered
    uint32_t* conn_cnt = md_alloc(arena, sizeof(uint32_t) * n);   // Bonds into the ordered atoms
    uint32_t* seq      = md_alloc(arena, sizeof(uint32_t) * n);   // When the atom joined the frontier
    uint32_t  num_seq  = 0;
    for (uint32_t i = 0; i < n; ++i) {
        pos_of[i]   = -1;
        conn_cnt[i] = 0;
        seq[i]      = 0;
    }
    // Leaves: atoms with one bond to a searched atom, which has more (hydrogens aside). They are ordered after all
    // the others, grouped by the atom they hang from, and outside of ALL each group is mapped as a set onto the
    // neighbours of that atom rather than atom by atom (see leaves_search): the three O of a sulfonate fit its S in
    // 3! ways, which are one set of atoms, and a molecule with ten of them would be tried 6^10 times over. The first
    // atom of a part is never a leaf, so that a rare one (the Cl of a chloride) still starts the search.
    bool*     leaf      = md_alloc(arena, sizeof(bool) * n);
    int32_t*  leaf_host = md_alloc(arena, sizeof(int32_t) * n);
    uint32_t* num_leaves = md_alloc(arena, sizeof(uint32_t) * n);
    int32_t*  part      = md_alloc(arena, sizeof(int32_t) * n);
    bool*     part_started = md_alloc(arena, sizeof(bool) * n);
    for (uint32_t i = 0; i < n; ++i) {
        leaf[i] = false;
        leaf_host[i] = -1;
        num_leaves[i] = 0;
        part[i] = (int32_t)i;
        part_started[i] = false;
    }
    for (uint32_t i = 0; i < m; ++i) {
        const uint32_t a = q->bonds[i].a, b = q->bonds[i].b;
        if (host[a] < 0 && host[b] < 0) uf_union(part, (int32_t)a, (int32_t)b);
    }
    for (uint32_t i = 0; i < n; ++i) {
        if (host[i] >= 0 || is_h[i] || degree[i] != 1) continue;
        for (uint32_t e = off[i]; e < off[i + 1]; ++e) {
            const uint32_t y = nbr[e];
            if (host[y] >= 0) continue;
            if (!is_h[y] && degree[y] >= 2) {
                leaf[i] = true;
                leaf_host[i] = (int32_t)y;
                num_leaves[y]++;
            }
            break;
        }
    }
    // An atom with very many of them (they are matched as subsets of its neighbours, MAX_LEAVES at most) keeps them
    for (uint32_t i = 0; i < n; ++i) {
        if (leaf[i] && num_leaves[leaf_host[i]] > MAX_LEAVES) leaf[i] = false;
    }

    md_array(uint32_t) frontier = 0;
    md_array_ensure(frontier, 16, arena);

    uint32_t num_back = 0;
    uint32_t k = 0;
    for (;;) {
        uint32_t x = UINT32_MAX;
        if (md_array_size(frontier) > 0) {
            // The atom with the most bonds back into the ordered atoms, so that rings close as early as possible.
            // Among equals the one longest in the frontier: breadth first, which keeps the atoms around a matched
            // atom close behind it in the order, and a wrong choice among its neighbours is found out a few steps
            // later rather than at the far end of a chain.
            size_t best = 0;
            for (size_t f = 1; f < md_array_size(frontier); ++f) {
                const uint32_t a = frontier[f], b = frontier[best];
                if (conn_cnt[a] != conn_cnt[b]) { if (conn_cnt[a] > conn_cnt[b]) best = f; continue; }
                if (seq[a]      != seq[b])      { if (seq[a]      < seq[b])      best = f; continue; }
                if (rank[a]     != rank[b])     { if (rank[a]     < rank[b])     best = f; continue; }
                if (degree[a]   != degree[b])   { if (degree[a]   > degree[b])   best = f; continue; }
                if (a < b) best = f;
            }
            x = frontier[best];
            md_array_swap_back_and_pop(frontier, best);
        } else {
            // First atom of a part not yet begun: the rarest, then the one with most bonds, then a named one
            for (uint32_t i = 0; i < n; ++i) {
                if (host[i] >= 0 || pos_of[i] >= 0 || part_started[uf_find(part, (int32_t)i)]) continue;
                if (x == UINT32_MAX) { x = i; continue; }
                if (rank[i]   != rank[x])   { if (rank[i]   < rank[x])   x = i; continue; }
                if (degree[i] != degree[x]) { if (degree[i] > degree[x]) x = i; continue; }
                const bool named_i = (q->atoms[i].flags & MD_MATCH_ATOM_NAME) != 0;
                const bool named_x = (q->atoms[x].flags & MD_MATCH_ATOM_NAME) != 0;
                if (named_i && !named_x) x = i;
            }
            if (x == UINT32_MAX) break;
            part_started[uf_find(part, (int32_t)x)] = true;
            if (leaf[x]) {
                leaf[x] = false;
                num_leaves[leaf_host[x]]--;
            }
        }

        pos_of[x] = (int32_t)k;
        plan->order[k]  = x;
        plan->atom[k]   = make_pos_atom(&q->atoms[x], num_folded[x]);
        for (uint32_t e = off[x]; e < off[x + 1]; ++e) {
            const uint32_t y = nbr[e];
            plan->atom[k].degree += host[y] < 0 && !is_h[y];
        }
        plan->parent[k] = -1;
        plan->parent_order[k] = MD_MATCH_BOND_ANY;
        if (plan->atom[k].h_lo > plan->atom[k].h_hi) {
            plan->impossible = true;
        }

        // The earliest ordered neighbour is the parent, the others are bonds back
        uint32_t parent_e = UINT32_MAX;
        for (uint32_t e = off[x]; e < off[x + 1]; ++e) {
            const uint32_t y = nbr[e];
            if (host[y] >= 0 || pos_of[y] < 0 || y == x) continue;
            if (parent_e == UINT32_MAX || pos_of[y] < pos_of[nbr[parent_e]]) parent_e = e;
        }
        plan->back_off[k] = num_back;
        if (parent_e != UINT32_MAX) {
            plan->parent[k]       = pos_of[nbr[parent_e]];
            plan->parent_order[k] = nbr_o[parent_e];
            for (uint32_t e = off[x]; e < off[x + 1]; ++e) {
                const uint32_t y = nbr[e];
                if (e == parent_e || host[y] >= 0 || pos_of[y] < 0 || y == x) continue;
                plan->back_pos[num_back]   = (uint32_t)pos_of[y];
                plan->back_order[num_back] = nbr_o[e];
                num_back++;
            }
        }

        // Its unordered neighbours join the frontier, the leaves aside
        for (uint32_t e = off[x]; e < off[x + 1]; ++e) {
            const uint32_t y = nbr[e];
            if (host[y] >= 0 || pos_of[y] >= 0 || leaf[y]) continue;
            if (conn_cnt[y]++ == 0) {
                seq[y] = num_seq++;
                md_array_push(frontier, y, arena);
            }
        }
        k++;
    }
    plan->num_core = k;

    // The leaves, grouped by the position of the atom they hang from, in query order within a group
    const uint32_t num_leaf = num_search - k;
    if (num_leaf > 0) {
        uint64_t* keys = md_alloc(arena, sizeof(uint64_t) * num_leaf);
        uint32_t f = 0;
        for (uint32_t i = 0; i < n; ++i) {
            if (leaf[i]) keys[f++] = ((uint64_t)pos_of[leaf_host[i]] << 32) | i;
        }
        ASSERT(f == num_leaf);
        qsort(keys, num_leaf, sizeof(uint64_t), cmp_u64);
        plan->leaf_beg  = md_alloc(arena, sizeof(uint32_t) * (num_leaf + 1));
        plan->leaf_host = md_alloc(arena, sizeof(uint32_t) * num_leaf);
        for (f = 0; f < num_leaf; ++f) {
            const uint32_t x = (uint32_t)keys[f];
            const uint32_t h = (uint32_t)(keys[f] >> 32);
            if (f == 0 || plan->leaf_host[plan->num_leaf_groups - 1] != h) {
                plan->leaf_beg[plan->num_leaf_groups]  = k;
                plan->leaf_host[plan->num_leaf_groups] = h;
                plan->num_leaf_groups++;
            }
            pos_of[x] = (int32_t)k;
            plan->order[k]  = x;
            plan->atom[k]   = make_pos_atom(&q->atoms[x], num_folded[x]);
            plan->atom[k].degree = 1;
            if (plan->atom[k].h_lo > plan->atom[k].h_hi) plan->impossible = true;
            plan->parent[k] = (int32_t)h;
            plan->parent_order[k] = MD_MATCH_BOND_ANY;
            for (uint32_t e = off[x]; e < off[x + 1]; ++e) {
                if (nbr[e] == (uint32_t)leaf_host[x]) plan->parent_order[k] = nbr_o[e];
            }
            plan->back_off[k] = num_back;
            k++;
        }
        plan->leaf_beg[plan->num_leaf_groups] = k;
    }
    ASSERT(k == num_search);
    plan->back_off[num_search] = num_back;

    // Folded hydrogens, grouped by the position of their atom, in query order within a group
    if (plan->num_folded > 0) {
        uint64_t* keys = md_alloc(arena, sizeof(uint64_t) * plan->num_folded);
        uint32_t f = 0;
        for (uint32_t i = 0; i < n; ++i) {
            if (host[i] >= 0) {
                keys[f++] = ((uint64_t)pos_of[host[i]] << 32) | i;
            }
        }
        qsort(keys, plan->num_folded, sizeof(uint64_t), cmp_u64);
        plan->folded_atom = md_alloc(arena, sizeof(uint32_t)   * plan->num_folded);
        plan->folded_host = md_alloc(arena, sizeof(uint32_t)   * plan->num_folded);
        plan->folded      = md_alloc(arena, sizeof(pos_atom_t) * plan->num_folded);
        for (f = 0; f < plan->num_folded; ++f) {
            const uint32_t i = (uint32_t)keys[f];
            plan->folded_atom[f] = i;
            plan->folded_host[f] = (uint32_t)(keys[f] >> 32);
            plan->folded[f]      = make_pos_atom(&q->atoms[i], 0);
            if (plan->folded[f].h_lo > plan->folded[f].h_hi) plan->impossible = true;
        }
    }

    // Size of the query as MD_MATCH_FLAG_WHOLE compares it
    for (uint32_t p = 0; p < num_search; ++p) {
        plan->num_heavy += !plan->atom[p].is_h;
    }
    for (uint32_t i = 0; i < m; ++i) {
        const uint32_t a = q->bonds[i].a, b = q->bonds[i].b;
        if (host[a] < 0 && host[b] < 0 && !is_h[a] && !is_h[b]) plan->num_heavy_bonds++;
    }

    if (symmetry && !plan->impossible) plan_symmetry(plan, q, pos_of, arena);
}

// ### UNIQUE ###
// The atom sets seen so far: sorted atom indices, hashed

typedef struct uniq_t {
    uint32_t  width;
    md_array(int32_t) keys;
    uint64_t* slot_hash;    // 0 for an empty slot
    uint32_t* slot_entry;
    uint32_t  cap;          // Power of two
    uint32_t  count;
    md_allocator_i* alloc;
} uniq_t;

static void uniq_grow(uniq_t* u) {
    const uint32_t new_cap = u->cap ? u->cap * 2 : 1024;
    uint64_t* hash  = md_alloc(u->alloc, sizeof(uint64_t) * new_cap);
    uint32_t* entry = md_alloc(u->alloc, sizeof(uint32_t) * new_cap);
    MEMSET(hash, 0, sizeof(uint64_t) * new_cap);
    for (uint32_t i = 0; i < u->cap; ++i) {
        if (!u->slot_hash[i]) continue;
        uint32_t s = (uint32_t)(u->slot_hash[i] & (new_cap - 1));
        while (hash[s]) s = (s + 1) & (new_cap - 1);
        hash[s]  = u->slot_hash[i];
        entry[s] = u->slot_entry[i];
    }
    if (u->cap) {
        md_free(u->alloc, u->slot_hash,  sizeof(uint64_t) * u->cap);
        md_free(u->alloc, u->slot_entry, sizeof(uint32_t) * u->cap);
    }
    u->slot_hash  = hash;
    u->slot_entry = entry;
    u->cap = new_cap;
}

// Returns true if the set was not seen before
static bool uniq_insert(uniq_t* u, const int32_t* key) {
    if ((u->count + 1) * 2 > u->cap) uniq_grow(u);
    const size_t bytes = sizeof(int32_t) * u->width;
    const uint64_t h = md_hash64(key, bytes, 0) | 1;
    uint32_t s = (uint32_t)(h & (u->cap - 1));
    while (u->slot_hash[s]) {
        if (u->slot_hash[s] == h && MEMCMP(u->keys + (size_t)u->slot_entry[s] * u->width, key, bytes) == 0) {
            return false;
        }
        s = (s + 1) & (u->cap - 1);
    }
    u->slot_hash[s]  = h;
    u->slot_entry[s] = u->count++;
    md_array_push_array(u->keys, key, u->width, u->alloc);
    return true;
}

static int cmp_i32(const void* a, const void* b) {
    const int32_t x = *(const int32_t*)a;
    const int32_t y = *(const int32_t*)b;
    return (x > y) - (x < y);
}

static void sort_i32(int32_t* arr, size_t n) {
    if (n <= 32) {
        for (size_t i = 1; i < n; ++i) {
            const int32_t v = arr[i];
            size_t j = i;
            while (j > 0 && arr[j - 1] > v) { arr[j] = arr[j - 1]; --j; }
            arr[j] = v;
        }
    } else {
        qsort(arr, n, sizeof(int32_t), cmp_i32);
    }
}

// ### SEARCH ###

// What the charges of a resonance group add up to over the atoms of one match, see groups_ok
typedef struct group_acc_t {
    uint32_t stamp;
    uint32_t count;         // Atoms of the group matched by query atoms which state a charge
    int32_t  q_query;       // Charges, stated and of the system
    int32_t  q_target;
    int32_t  v_query;       // Hydrogens minus charge, stated and of the system
    int32_t  v_target;
    bool     proton;        // Some of the atoms are of unknown protonation
    bool     open;          // ... and do not state their hydrogens exactly: the group cannot be compared
} group_acc_t;

typedef struct search_t {
    const target_t* tg;
    const plan_t*   plan;

    md_match_level_t   level;
    md_match_mode_t    mode;
    bool               whole;
    md_match_resolve_t hydrogens;
    md_match_resolve_t bond_orders;
    md_match_resolve_t charges;
    size_t             limit;

    uint64_t* mask;         // Atoms which may be matched, NULL for all
    uint64_t* mapped;       // Atoms of the current partial mapping
    uint64_t* used;         // DISJOINT: atoms of the matches reported

    // The current unit: its atoms are unit_atoms[0 .. unit_count) or, without unit_atoms, the range
    // [unit_beg, unit_end) of which the graph is cut at the edges (ranged)
    int64_t        unit;
    bool           ranged;
    const int32_t* unit_atoms;
    uint32_t       unit_beg;
    uint32_t       unit_end;
    uint32_t       unit_count;
    bool           unit_done;

    int32_t*  map;          // [plan.num] Atom at each position
    uint32_t* cursor;       // [plan.num] Where the search for the candidates of each position is
    int32_t*  row;          // [plan.num_atoms] The match in query order
    int32_t*  key;          // [plan.num] UNIQUE: sorted atoms
    uniq_t    uniq;

    // Leaves, see leaves_search
    md_array(struct leaf_cand_t) leaf_cand;     // The candidates of the groups entered, one group after the other
    uint32_t* leaf_base;    // [plan.num_leaf_groups + 1] Where the candidates of each group begin
    uint32_t* leaf_comb;    // [plan.num] The candidates chosen for a group, ascending indices (by leaf position)
    bool*     leaf_fresh;   // [plan.num_leaf_groups] The group was just entered: its first choice is still to be tried

    // The molecule asked about last, see search_resolution
    struct res_cache_t* res_cache;

    group_acc_t* group_acc;     // [tg.num_groups]
    int32_t*     group_touched; // [plan.num]
    uint32_t     group_stamp;

    md_match_callback_t callback;
    void*    user_param;
    size_t   count;
    md_allocator_i* arena;
    bool     stop;          // The search is over
    bool     root_done;     // The search from the current start atom is over
} search_t;

typedef struct res_cache_t {
    int64_t  mol;
    uint32_t slot_beg;      // Where the molecules are the structures, its slots
    uint32_t slot_end;
    const resolution_t* res;
} res_cache_t;

// What the molecule of atom t can tell. The atoms of a search are mostly of one molecule.
static inline const resolution_t* search_resolution(const search_t* s, int32_t t) {
    res_cache_t* rc = s->res_cache;
    const target_t* tg = s->tg;
    if (tg->mol_of) {
        const int32_t m = tg->mol_of[t];
        if (m != rc->mol) {
            rc->mol = m;
            rc->res = molecule_resolution(tg, m);
        }
    } else {
        const uint32_t slot = (uint32_t)tg->mol_slot[t];
        if (rc->mol < 0 || slot < rc->slot_beg || slot >= rc->slot_end) {
            const int64_t m = range_find(tg->mol_off, tg->num_mol, slot);
            rc->mol      = m;
            rc->slot_beg = tg->mol_off[m];
            rc->slot_end = tg->mol_off[m + 1];
            rc->res      = molecule_resolution(tg, m);
        }
    }
    return rc->res;
}

enum {
    CHEM_NONE   = 0,        // Not tested
    CHEM_EXACT  = 1,        // Tested as it is
    CHEM_PROTON = 2,        // Known up to a proton: tested as hydrogens minus charge
};

// The atoms whose hydrogens and charge are a guess at their protonation where the molecule has no hydrogen atoms
static inline bool protonation_site(md_atomic_number_t z) {
    return z == MD_Z_N || z == MD_Z_O || z == MD_Z_P || z == MD_Z_S || z == 34;
}

static inline int hydrogens_mode(const search_t* s, int32_t t, md_atomic_number_t z) {
    switch (s->hydrogens) {
    case MD_MATCH_RESOLVE_ALWAYS: return CHEM_EXACT;
    case MD_MATCH_RESOLVE_NEVER:  return CHEM_NONE;
    default: break;
    }
    const resolution_t* res = search_resolution(s, t);
    if (s->tg->h_count) {
        return (protonation_site(z) && !res->any_h) ? CHEM_PROTON : CHEM_EXACT;
    }
    return resolution_h(res, z) ? CHEM_EXACT : CHEM_NONE;
}

static inline int charge_mode(const search_t* s, int32_t t, md_atomic_number_t z) {
    switch (s->charges) {
    case MD_MATCH_RESOLVE_ALWAYS: return CHEM_EXACT;
    case MD_MATCH_RESOLVE_NEVER:  return CHEM_NONE;
    default: break;
    }
    if (!s->tg->charge) return CHEM_NONE;
    return (protonation_site(z) && !search_resolution(s, t)->any_h) ? CHEM_PROTON : CHEM_EXACT;
}

static inline bool in_group(const search_t* s, int32_t t) {
    return s->tg->group_of && s->tg->group_of[t] >= 0;
}

// The constraints of a query atom on its own
static inline bool atom_props_ok(const search_t* s, const pos_atom_t* pa, int32_t t) {
    const md_atomic_number_t z = target_z(s->tg, t);
    if (pa->is_h != (z == MD_Z_H)) return false;
    if ((pa->flags & MD_MATCH_ATOM_ELEMENT) && pa->z != z) return false;
    if ((pa->flags & MD_MATCH_ATOM_NAME) && !label_eq(pa->name, target_name(s->tg, t))) return false;
    if ((pa->flags & (MD_MATCH_ATOM_AROMATIC | MD_MATCH_ATOM_ALIPHATIC)) && s->bond_orders != MD_MATCH_RESOLVE_NEVER) {
        const int aromatic = target_aromatic(s->tg, t);
        if (pa->flags & MD_MATCH_ATOM_AROMATIC) {
            if (aromatic == 0 || (aromatic < 0 && s->bond_orders == MD_MATCH_RESOLVE_ALWAYS)) return false;
        }
        if ((pa->flags & MD_MATCH_ATOM_ALIPHATIC) && aromatic == 1) return false;
    }
    if (pa->h_test || (pa->flags & MD_MATCH_ATOM_CHARGE)) {
        const int h_mode = hydrogens_mode(s, t, z);
        if (pa->h_test && h_mode == CHEM_EXACT) {
            const uint32_t h = target_hydrogens(s->tg, t);
            if (h < pa->h_lo || h > pa->h_hi) return false;
        }
        if (pa->flags & MD_MATCH_ATOM_CHARGE) {
            const bool grouped = in_group(s, t);
            const int q_mode = charge_mode(s, t, z);
            if (q_mode == CHEM_EXACT) {
                // The charges of a resonance group are compared once the match is complete (groups_ok)
                if (!grouped && target_charge(s->tg, t) != pa->charge) return false;
            } else if (q_mode == CHEM_PROTON && h_mode == CHEM_PROTON && pa->h_stated) {
                // Hydrogens minus charge, which a proton leaves as they are. Within a resonance group they also move
                // by one with a double bond, which md_chem may have put on another atom of the group.
                const int slack = grouped ? 1 : 0;
                const int v = (int)target_hydrogens(s->tg, t) - target_charge(s->tg, t);
                if (v < (int)pa->h_lo - pa->charge - slack || v > (int)pa->h_hi - pa->charge + slack) return false;
            }
        }
    }
    return true;
}

// And where it may be
static inline bool in_graph(const search_t* s, int32_t t) {
    if (s->mask && !bit_test(s->mask, t)) return false;
    if (s->ranged && ((uint32_t)t < s->unit_beg || (uint32_t)t >= s->unit_end)) return false;
    return true;
}

// Bonds to atoms which are not hydrogens, within the unit and the mask
static inline uint32_t heavy_degree(const search_t* s, int32_t t) {
    uint32_t d = 0;
    for (uint32_t c = conn_beg(s->tg, t); c < conn_end(s->tg, t); ++c) {
        const int32_t b = s->tg->conn_atom[c];
        d += target_z(s->tg, b) != MD_Z_H && in_graph(s, b);
    }
    return d;
}

static inline bool atom_ok(const search_t* s, const pos_atom_t* pa, int32_t t) {
    if (bit_test(s->mapped, t)) return false;
    if (!in_graph(s, t)) return false;
    if (s->used && bit_test(s->used, t)) return false;
    if (!atom_props_ok(s, pa, t)) return false;
    // An atom with fewer bonds than the query atom cannot take its place; with WHOLE not with more either
    if (pa->degree > 1 || s->whole) {
        const uint32_t d = heavy_degree(s, t);
        if (d < pa->degree || (s->whole && d != pa->degree)) return false;
    }
    return true;
}

// See md_match_bond_order_t
static inline bool bond_ok(const search_t* s, uint8_t q_order, md_bond_idx_t b) {
    if (q_order == MD_MATCH_BOND_ANY || s->bond_orders == MD_MATCH_RESOLVE_NEVER) return true;
    const md_bond_flags_t f = target_bond_flags(s->tg, b);
    int order = md_bond_order(f);
    const bool resonance = (f & (MD_BOND_FLAG_AROMATIC | MD_BOND_FLAG_DELOCALIZED)) != 0;
    if (order == MD_BOND_ORDER_UNKNOWN && !resonance) {
        if (s->bond_orders != MD_MATCH_RESOLVE_ALWAYS) return true;
        order = MD_BOND_ORDER_SINGLE;
    }

    switch (q_order) {
    case MD_MATCH_BOND_SINGLE:    return order == MD_BOND_ORDER_SINGLE || resonance;
    case MD_MATCH_BOND_DOUBLE:    return order == MD_BOND_ORDER_DOUBLE || resonance;
    case MD_MATCH_BOND_AROMATIC:  return order == MD_BOND_ORDER_SINGLE || order == MD_BOND_ORDER_DOUBLE || resonance;
    case MD_MATCH_BOND_TRIPLE:    return order == MD_BOND_ORDER_TRIPLE;
    case MD_MATCH_BOND_QUADRUPLE: return order == MD_BOND_ORDER_QUADRUPLE;
    default:                      return true;
    }
}

static inline bool back_bonds_ok(const search_t* s, uint32_t k, int32_t t) {
    const plan_t* pl = s->plan;
    for (uint32_t e = pl->back_off[k]; e < pl->back_off[k + 1]; ++e) {
        const int32_t other = s->map[pl->back_pos[e]];
        bool found = false;
        for (uint32_t c = conn_beg(s->tg, t); c < conn_end(s->tg, t); ++c) {
            if (s->tg->conn_atom[c] == other) {
                if (!bond_ok(s, pl->back_order[e], s->tg->conn_bond[c])) return false;
                found = true;
                break;
            }
        }
        if (!found) return false;
    }
    return true;
}

// The conditions which break the symmetries of the query, see plan_symmetry
static inline bool symmetry_ok(const search_t* s, uint32_t k, int32_t t) {
    const plan_t* pl = s->plan;
    if (!pl->sym_off || s->mode == MD_MATCH_MODE_ALL) return true;
    for (uint32_t i = pl->sym_off[k]; i < pl->sym_off[k + 1]; ++i) {
        const int32_t other = s->map[pl->sym_other[i]];
        if (pl->sym_less[i] ? !(t < other) : !(t > other)) return false;
    }
    return true;
}

static inline int32_t unit_atom(const search_t* s, uint32_t i) {
    return s->unit_atoms ? s->unit_atoms[i] : (int32_t)(s->unit_beg + i);
}

static inline void cursor_begin(search_t* s, uint32_t k) {
    const int32_t p = s->plan->parent[k];
    s->cursor[k] = p >= 0 ? conn_beg(s->tg, s->map[p]) : 0;
}

static bool next_candidate(search_t* s, uint32_t k, int32_t* out) {
    const plan_t* pl = s->plan;
    const pos_atom_t* pa = &pl->atom[k];
    const int32_t p = pl->parent[k];
    if (p >= 0) {
        const uint32_t end = conn_end(s->tg, s->map[p]);
        while (s->cursor[k] < end) {
            const uint32_t c = s->cursor[k]++;
            const int32_t t = s->tg->conn_atom[c];
            if (!atom_ok(s, pa, t)) continue;
            if (!bond_ok(s, pl->parent_order[k], s->tg->conn_bond[c])) continue;
            if (!back_bonds_ok(s, k, t)) continue;
            if (!symmetry_ok(s, k, t)) continue;
            *out = t;
            return true;
        }
    } else {
        // The first atom of a later part of the query: anywhere in the unit
        while (s->cursor[k] < s->unit_count) {
            const int32_t t = unit_atom(s, s->cursor[k]++);
            if (target_skipped(s->tg, t)) continue;
            if (!atom_ok(s, pa, t)) continue;
            if (!symmetry_ok(s, k, t)) continue;
            *out = t;
            return true;
        }
    }
    return false;
}

// Maps the folded hydrogens onto hydrogens of the atoms they are bonded to. Returns false if the match fails on them.
static bool map_hydrogens(search_t* s) {
    const plan_t* pl = s->plan;
    uint32_t i = 0;
    while (i < pl->num_folded) {
        const uint32_t host_pos = pl->folded_host[i];
        uint32_t j = i;
        while (j < pl->num_folded && pl->folded_host[j] == host_pos) ++j;

        const int32_t host = s->map[host_pos];
        int32_t  hs[MAX_HYDROGENS_PER_ATOM];
        uint32_t nh = 0;
        for (uint32_t c = conn_beg(s->tg, host); c < conn_end(s->tg, host) && nh < MAX_HYDROGENS_PER_ATOM; ++c) {
            const int32_t b = s->tg->conn_atom[c];
            if (target_z(s->tg, b) != MD_Z_H) continue;
            uint32_t x = nh++;
            while (x > 0 && hs[x - 1] > b) { hs[x] = hs[x - 1]; --x; }
            hs[x] = b;
        }

        uint32_t taken = 0;
        uint32_t num_taken = 0;
        for (uint32_t f = i; f < j; ++f) {
            int32_t found = -1;
            for (uint32_t x = 0; x < nh; ++x) {
                if (taken & (1u << x)) continue;
                if (atom_props_ok(s, &pl->folded[f], hs[x])) {
                    found = (int32_t)x;
                    break;
                }
            }
            if (found >= 0) {
                taken |= 1u << found;
                num_taken++;
                s->row[pl->folded_atom[f]] = hs[found];
            } else if (num_taken < nh) {
                // There is a hydrogen, but not the one asked for
                return false;
            } else {
                s->row[pl->folded_atom[f]] = -1;
            }
        }
        i = j;
    }
    return true;
}

// The charges of resonance groups, see CHARGES in md_match.h. Over the matched atoms of a group whose query atoms
// state a charge: where the protonation is known, their charges add up to those of the system on the same atoms, or,
// as the charge of the group may sit on any of its atoms, to a share of the group's charge (0 up to all of it), and to
// all of it when they are the whole group. Where it is not known, hydrogens minus charge add up as on the system, and
// only the whole group can be compared: a proton may sit on any of its atoms as well.
static bool groups_ok(search_t* s) {
    const plan_t* pl = s->plan;
    const target_t* tg = s->tg;
    if (!pl->has_charge || !tg->group_of) return true;

    if (++s->group_stamp == 0) {
        for (size_t g = 0; g < tg->num_groups; ++g) s->group_acc[g].stamp = 0;
        s->group_stamp = 1;
    }
    uint32_t num_touched = 0;
    for (uint32_t p = 0; p < pl->num; ++p) {
        const pos_atom_t* pa = &pl->atom[p];
        if (!(pa->flags & MD_MATCH_ATOM_CHARGE)) continue;
        const int32_t t = s->map[p];
        const int32_t g = tg->group_of[t];
        if (g < 0) continue;
        const md_atomic_number_t z = target_z(tg, t);
        const int q_mode = charge_mode(s, t, z);
        if (q_mode == CHEM_NONE) continue;

        group_acc_t* acc = &s->group_acc[g];
        if (acc->stamp != s->group_stamp) {
            MEMSET(acc, 0, sizeof(group_acc_t));
            acc->stamp = s->group_stamp;
            s->group_touched[num_touched++] = g;
        }
        const int h_t = (int)target_hydrogens(tg, t);
        const int q_t = target_charge(tg, t);
        acc->count += 1;
        if (q_mode == CHEM_EXACT) {
            acc->q_query  += pa->charge;
            acc->q_target += q_t;
            acc->v_query  += (pa->h_stated ? (int)pa->h_lo : h_t) - pa->charge;
            acc->v_target += h_t - q_t;
        } else {
            acc->proton = true;
            if (hydrogens_mode(s, t, z) == CHEM_PROTON && pa->h_stated && pa->h_lo == pa->h_hi) {
                acc->v_query  += (int)pa->h_lo - pa->charge;
                acc->v_target += h_t - q_t;
            } else {
                acc->open = true;
            }
        }
    }

    for (uint32_t i = 0; i < num_touched; ++i) {
        const int32_t g = s->group_touched[i];
        const group_acc_t* acc = &s->group_acc[g];
        const bool whole_group = acc->count == tg->group_size[g];
        if (acc->proton) {
            if (acc->open || !whole_group) continue;
            if (acc->v_query != acc->v_target) return false;
        } else {
            if (acc->q_query == acc->q_target) continue;
            if (whole_group) return false;
            const int32_t q = tg->group_charge[g];
            if (acc->q_query < MIN(0, q) || acc->q_query > MAX(0, q)) return false;
        }
    }
    return true;
}

static void report(search_t* s) {
    const plan_t* pl = s->plan;
    for (uint32_t p = 0; p < pl->num; ++p) {
        s->row[pl->order[p]] = s->map[p];
    }
    if (pl->num_folded > 0 && !map_hydrogens(s)) return;
    if (!groups_ok(s)) return;

    switch (s->mode) {
    case MD_MATCH_MODE_UNIQUE:
        if (pl->num > 1) {
            MEMCPY(s->key, s->map, sizeof(int32_t) * pl->num);
            sort_i32(s->key, pl->num);
            if (!uniq_insert(&s->uniq, s->key)) return;
        }
        break;
    case MD_MATCH_MODE_ONE_PER_UNIT:
        s->unit_done = true;
        s->root_done = true;
        break;
    case MD_MATCH_MODE_DISJOINT:
        for (uint32_t p = 0; p < pl->num; ++p) {
            bit_set(s->used, s->map[p]);
        }
        s->root_done = true;
        break;
    default:
        break;
    }

    s->count++;
    if (s->callback && !s->callback(s->row, pl->num_atoms, (uint32_t)s->unit, s->user_param)) {
        s->stop = true;
    }
    if (s->limit && s->count >= s->limit) {
        s->stop = true;
    }
}

// ### LEAVES ###
// Outside of ALL, the leaves of a group (see plan_build) are mapped as a set: for each set of as many neighbours of the
// atom they hang from as there are leaves, one assignment of the leaves to them which satisfies their constraints, if
// there is one (a bipartite matching). The search is the same backtracking as the one over atoms, with a group of
// leaves as a step and a set of neighbours as a choice.

typedef struct leaf_cand_t {
    int32_t  atom;
    uint32_t compat;        // The leaves of the group it may be: bit i for the i-th
} leaf_cand_t;

static bool leaves_augment(uint32_t j, const uint32_t* compat, uint32_t k, int8_t* cand_of_leaf, uint32_t* visited) {
    for (uint32_t i = 0; i < k; ++i) {
        if (!((compat[j] >> i) & 1) || ((*visited >> i) & 1)) continue;
        *visited |= 1u << i;
        if (cand_of_leaf[i] < 0 || leaves_augment((uint32_t)cand_of_leaf[i], compat, k, cand_of_leaf, visited)) {
            cand_of_leaf[i] = (int8_t)j;
            return true;
        }
    }
    return false;
}

// The k chosen candidates onto the k leaves (Kuhn): cand_of_leaf[i] is the candidate (0 .. k-1) of the i-th leaf.
// Each candidate takes the first leaf free for it where there is one, so that leaves which are alike take the
// candidates in order.
static bool leaves_assign(const uint32_t* compat, uint32_t k, int8_t* cand_of_leaf) {
    for (uint32_t i = 0; i < k; ++i) cand_of_leaf[i] = -1;
    for (uint32_t j = 0; j < k; ++j) {
        uint32_t i = 0;
        while (i < k && !(((compat[j] >> i) & 1) && cand_of_leaf[i] < 0)) ++i;
        if (i < k) {
            cand_of_leaf[i] = (int8_t)j;
            continue;
        }
        uint32_t visited = 0;
        if (!leaves_augment(j, compat, k, cand_of_leaf, &visited)) return false;
    }
    return true;
}

static void leaves_enter(search_t* s, uint32_t g) {
    const plan_t* pl = s->plan;
    const uint32_t beg = pl->leaf_beg[g];
    const uint32_t k   = pl->leaf_beg[g + 1] - beg;
    const int32_t host = s->map[pl->leaf_host[g]];
    md_array_shrink(s->leaf_cand, s->leaf_base[g]);
    for (uint32_t c = conn_beg(s->tg, host); c < conn_end(s->tg, host); ++c) {
        const int32_t t = s->tg->conn_atom[c];
        uint32_t compat = 0;
        for (uint32_t i = 0; i < k; ++i) {
            const uint32_t p = beg + i;
            if (atom_ok(s, &pl->atom[p], t) && bond_ok(s, pl->parent_order[p], s->tg->conn_bond[c])) compat |= 1u << i;
        }
        if (compat) {
            // In ascending order of the atoms, so that equal leaves take the atoms in their own order
            leaf_cand_t cand = {.atom = t, .compat = compat};
            md_array_push(s->leaf_cand, cand, s->arena);
            size_t x = md_array_size(s->leaf_cand) - 1;
            while (x > s->leaf_base[g] && s->leaf_cand[x - 1].atom > t) {
                s->leaf_cand[x] = s->leaf_cand[x - 1];
                --x;
            }
            s->leaf_cand[x] = cand;
        }
    }
    s->leaf_base[g + 1] = (uint32_t)md_array_size(s->leaf_cand);
    for (uint32_t i = 0; i < k; ++i) s->leaf_comb[beg + i] = i;
    s->leaf_fresh[g] = true;
}

// The next set of candidates for group g which its leaves can be assigned to, mapped. False when there is none left.
static bool leaves_next(search_t* s, uint32_t g) {
    const plan_t* pl = s->plan;
    const uint32_t beg = pl->leaf_beg[g];
    const uint32_t k   = pl->leaf_beg[g + 1] - beg;
    const uint32_t m   = s->leaf_base[g + 1] - s->leaf_base[g];
    const leaf_cand_t* cand = s->leaf_cand + s->leaf_base[g];
    uint32_t* comb = s->leaf_comb + beg;
    if (m < k) return false;
    for (;;) {
        if (s->leaf_fresh[g]) {
            s->leaf_fresh[g] = false;
        } else {
            // The next combination, in lexicographic order
            int32_t i = (int32_t)k - 1;
            while (i >= 0 && comb[i] == m - k + (uint32_t)i) --i;
            if (i < 0) return false;
            comb[i]++;
            for (uint32_t j = (uint32_t)i + 1; j < k; ++j) comb[j] = comb[j - 1] + 1;
        }
        uint32_t compat[MAX_LEAVES];
        uint32_t all = 0;
        for (uint32_t i = 0; i < k; ++i) {
            compat[i] = cand[comb[i]].compat;
            all |= compat[i];
        }
        if (all != (1u << k) - 1) continue;
        int8_t cand_of_leaf[MAX_LEAVES];
        if (!leaves_assign(compat, k, cand_of_leaf)) continue;
        for (uint32_t i = 0; i < k; ++i) {
            s->map[beg + i] = cand[comb[cand_of_leaf[i]]].atom;
        }
        return true;
    }
}

static void leaves_mark(search_t* s, uint32_t g, bool set) {
    const plan_t* pl = s->plan;
    for (uint32_t p = pl->leaf_beg[g]; p < pl->leaf_beg[g + 1]; ++p) {
        if (set) bit_set(s->mapped, s->map[p]);
        else bit_clear(s->mapped, s->map[p]);
    }
}

// With the atoms before them mapped (and marked), the leaves
static void leaves_search(search_t* s) {
    const plan_t* pl = s->plan;
    const uint32_t num_groups = pl->num_leaf_groups;
    s->leaf_base[0] = 0;
    md_array_shrink(s->leaf_cand, 0);

    uint32_t g = 0;
    leaves_enter(s, 0);
    for (;;) {
        if (leaves_next(s, g)) {
            if (g + 1 == num_groups) {
                report(s);
                if (s->stop || s->root_done) break;
                continue;
            }
            leaves_mark(s, g, true);
            ++g;
            leaves_enter(s, g);
        } else {
            if (g == 0) break;
            --g;
            leaves_mark(s, g, false);
        }
    }
    for (uint32_t i = 0; i < g; ++i) leaves_mark(s, i, false);
}

// The atoms before the leaves are mapped: on to the leaves, or the match is complete
static void core_mapped(search_t* s, uint32_t n) {
    if (n == s->plan->num) {
        report(s);
        return;
    }
    const int32_t last = s->map[n - 1];
    bit_set(s->mapped, last);
    leaves_search(s);
    bit_clear(s->mapped, last);
}

static void search_from_root(search_t* s) {
    // ALL reports every mapping, the leaves among them: they are searched atom by atom like the others
    const uint32_t n = s->mode == MD_MATCH_MODE_ALL ? s->plan->num : s->plan->num_core;
    if (n == 1) {
        core_mapped(s, n);
        return;
    }

    uint32_t k = 1;
    cursor_begin(s, k);
    for (;;) {
        int32_t t;
        if (next_candidate(s, k, &t)) {
            s->map[k] = t;
            if (k + 1 == n) {
                core_mapped(s, n);
                if (s->stop || s->root_done) break;
                continue;
            }
            bit_set(s->mapped, t);
            ++k;
            cursor_begin(s, k);
        } else {
            --k;
            if (k == 0) break;
            bit_clear(s->mapped, s->map[k]);
        }
    }

    // Positions 1 .. k-1 are still mapped if the search stopped early
    for (uint32_t i = 1; i < k; ++i) {
        bit_clear(s->mapped, s->map[i]);
    }
}

static size_t level_num_units(const target_t* tg, md_match_level_t level) {
    switch (level) {
    case MD_MATCH_LEVEL_STRUCTURE: return tg->num_mol;
    case MD_MATCH_LEVEL_COMPONENT: return tg->sys->component.count;
    case MD_MATCH_LEVEL_INSTANCE:  return tg->sys->instance.count;
    default: return 0;
    }
}

static void unit_begin(search_t* s, int64_t unit) {
    const target_t* tg = s->tg;
    const md_system_t* sys = tg->sys;
    s->unit = unit;
    s->unit_atoms = NULL;
    uint32_t beg = 0, end = 0;
    switch (s->level) {
    case MD_MATCH_LEVEL_STRUCTURE:
        s->unit_atoms = tg->mol_atom + tg->mol_off[unit];
        s->unit_count = tg->mol_off[unit + 1] - tg->mol_off[unit];
        return;
    case MD_MATCH_LEVEL_COMPONENT:
        beg = sys->component.atom_offset[unit];
        end = sys->component.atom_offset[unit + 1];
        break;
    case MD_MATCH_LEVEL_INSTANCE:
        beg = sys->component.atom_offset[sys->instance.comp_offset[unit]];
        end = sys->component.atom_offset[sys->instance.comp_offset[unit + 1]];
        break;
    default:
        ASSERT(false);
        break;
    }
    end = MIN(end, (uint32_t)tg->num_atoms);
    beg = MIN(beg, end);
    s->unit_beg   = beg;
    s->unit_end   = end;
    s->unit_count = end - beg;
}

// The heavy atoms of the current unit and the bonds between them, within the mask
static void unit_size(const search_t* s, uint32_t* out_heavy, uint32_t* out_bonds) {
    uint32_t num_heavy = 0;
    uint32_t num_bonds = 0;
    for (uint32_t i = 0; i < s->unit_count; ++i) {
        const int32_t a = unit_atom(s, i);
        if (target_skipped(s->tg, a)) continue;
        if (s->mask && !bit_test(s->mask, a)) continue;
        if (target_z(s->tg, a) == MD_Z_H) continue;
        num_heavy++;
        for (uint32_t c = conn_beg(s->tg, a); c < conn_end(s->tg, a); ++c) {
            const int32_t b = s->tg->conn_atom[c];
            if (b <= a || target_z(s->tg, b) == MD_Z_H || !in_graph(s, b)) continue;
            num_bonds++;
        }
    }
    *out_heavy = num_heavy;
    *out_bonds = num_bonds;
}

// Searches the current unit (unit_begin) for the current plan. whole_checked: the caller has made sure the unit has the
// size of the query, for MD_MATCH_FLAG_WHOLE.
static void unit_search(search_t* s, bool whole_checked) {
    const plan_t* pl = s->plan;
    s->unit_done = false;
    if (s->unit_count < pl->num) return;
    if (s->whole && !whole_checked) {
        uint32_t num_heavy, num_bonds;
        unit_size(s, &num_heavy, &num_bonds);
        if (num_heavy != pl->num_heavy || num_bonds != pl->num_heavy_bonds) return;
    }

    const pos_atom_t* p0 = &pl->atom[0];
    for (uint32_t i = 0; i < s->unit_count && !s->stop && !s->unit_done; ++i) {
        const int32_t t = unit_atom(s, i);
        // Cheap rejections first
        const md_atomic_number_t z = target_z(s->tg, t);
        if (p0->is_h != (z == MD_Z_H)) continue;
        if ((p0->flags & MD_MATCH_ATOM_ELEMENT) && p0->z != z) continue;
        if (target_skipped(s->tg, t)) continue;
        if (!atom_ok(s, p0, t)) continue;

        s->map[0] = t;
        bit_set(s->mapped, t);
        s->root_done = false;
        search_from_root(s);
        bit_clear(s->mapped, t);
    }
}

static void search_init(search_t* s, const target_t* tg, md_match_level_t level, md_match_mode_t mode, uint32_t flags, md_match_resolve_t hydrogens,
                        md_match_resolve_t bond_orders, md_match_resolve_t charges, const md_bitfield_t* mask, md_allocator_i* arena) {
    MEMSET(s, 0, sizeof(search_t));
    s->tg          = tg;
    s->level       = level;
    s->mode        = mode;
    s->whole       = (flags & MD_MATCH_FLAG_WHOLE) != 0;
    s->hydrogens   = hydrogens;
    s->bond_orders = bond_orders;
    s->charges     = charges;
    s->ranged      = level != MD_MATCH_LEVEL_STRUCTURE;
    s->unit        = -1;
    s->mapped      = bits_create(tg->num_atoms, arena);
    s->res_cache   = md_alloc(arena, sizeof(res_cache_t));
    MEMSET(s->res_cache, 0, sizeof(res_cache_t));
    s->res_cache->mol = -1;
    if (mask) {
        s->mask = bits_create(tg->num_atoms, arena);
        md_bitfield_iter_t it = md_bitfield_iter_create(mask);
        while (md_bitfield_iter_next(&it)) {
            const uint64_t i = md_bitfield_iter_idx(&it);
            if (i >= tg->num_atoms) break;
            bit_set(s->mask, i);
        }
    }
    if (mode == MD_MATCH_MODE_DISJOINT) s->used = bits_create(tg->num_atoms, arena);
    if (tg->group_of) {
        s->group_acc = md_alloc(arena, sizeof(group_acc_t) * MAX(1, tg->num_groups));
        MEMSET(s->group_acc, 0, sizeof(group_acc_t) * MAX(1, tg->num_groups));
    }
}

// Room for plans of up to num positions and num_atoms atoms
static void search_reserve(search_t* s, uint32_t num, size_t num_atoms, md_allocator_i* arena) {
    s->map           = md_alloc(arena, sizeof(int32_t)  * MAX(1, num));
    s->cursor        = md_alloc(arena, sizeof(uint32_t) * MAX(1, num));
    s->key           = md_alloc(arena, sizeof(int32_t)  * MAX(1, num));
    s->group_touched = md_alloc(arena, sizeof(int32_t)  * MAX(1, num));
    s->row           = md_alloc(arena, sizeof(int32_t)  * MAX(1, num_atoms));
    s->leaf_base     = md_alloc(arena, sizeof(uint32_t) * (num + 1));
    s->leaf_comb     = md_alloc(arena, sizeof(uint32_t) * MAX(1, num));
    s->leaf_fresh    = md_alloc(arena, sizeof(bool)     * MAX(1, num));
    s->leaf_cand     = 0;
    md_array_ensure(s->leaf_cand, 64, arena);
    s->arena         = arena;
}

// How common each element is among the candidate atoms, to start from the rarest. Returns the number of candidates.
static uint32_t count_candidates(uint32_t rarity[256], const search_t* s) {
    MEMSET(rarity, 0, sizeof(uint32_t) * 256);
    uint32_t num_candidates = 0;
    for (size_t i = 0; i < s->tg->num_atoms; ++i) {
        if (target_skipped(s->tg, (int32_t)i)) continue;
        if (s->mask && !bit_test(s->mask, i)) continue;
        rarity[target_z(s->tg, (int32_t)i)]++;
        num_candidates++;
    }
    return num_candidates;
}

static bool resolve_valid(md_match_resolve_t r) {
    return (uint32_t)r <= MD_MATCH_RESOLVE_NEVER;
}

static bool level_valid(md_match_level_t level, const md_system_t* sys) {
    if ((uint32_t)level > MD_MATCH_LEVEL_INSTANCE) {
        MD_LOG_ERROR("md_match: unknown level");
        return false;
    }
    if (level == MD_MATCH_LEVEL_COMPONENT && !(sys->component.count > 0 && sys->component.atom_offset)) {
        MD_LOG_ERROR("md_match: matching within components, but the system has none");
        return false;
    }
    if (level == MD_MATCH_LEVEL_INSTANCE && !(sys->instance.count > 0 && sys->instance.comp_offset && sys->component.atom_offset)) {
        MD_LOG_ERROR("md_match: matching within instances, but the system has none");
        return false;
    }
    return true;
}

static bool desc_validate(const md_match_desc_t* desc, const md_system_t* sys) {
    if (!desc) {
        MD_LOG_ERROR("md_match: no description");
        return false;
    }
    if (!sys) {
        MD_LOG_ERROR("md_match: no system");
        return false;
    }
    if (!level_valid(desc->level, sys)) return false;
    if ((uint32_t)desc->mode > MD_MATCH_MODE_DISJOINT) {
        MD_LOG_ERROR("md_match: unknown mode");
        return false;
    }
    if (!resolve_valid(desc->hydrogens) || !resolve_valid(desc->bond_orders) || !resolve_valid(desc->charges)) {
        MD_LOG_ERROR("md_match: unknown resolution");
        return false;
    }
    return query_validate(desc->query);
}

static uint32_t target_needs(md_match_resolve_t hydrogens, md_match_resolve_t charges, bool has_charge) {
    uint32_t what = 0;
    if (hydrogens == MD_MATCH_RESOLVE_AUTO || charges == MD_MATCH_RESOLVE_AUTO) what |= TARGET_RESOLUTION;
    if (has_charge && charges != MD_MATCH_RESOLVE_NEVER) what |= TARGET_GROUPS;
    return what;
}

// Does the system have as many atoms of each element as the query asks for (hydrogens folded into their atoms aside)
static bool elements_suffice(const md_match_query_t* q, const md_system_t* sys) {
    uint32_t need[256] = {0};
    md_temp_scope_t temp = md_temp_begin();
    uint32_t* degree = md_temp_alloc_array(temp, uint32_t, q->num_atoms);
    int32_t*  other  = md_temp_alloc_array(temp, int32_t, q->num_atoms);
    MEMSET(degree, 0, sizeof(uint32_t) * q->num_atoms);
    for (size_t i = 0; i < q->num_bonds; ++i) {
        degree[q->bonds[i].a]++;
        degree[q->bonds[i].b]++;
        other[q->bonds[i].a] = (int32_t)q->bonds[i].b;
        other[q->bonds[i].b] = (int32_t)q->bonds[i].a;
    }
    for (size_t i = 0; i < q->num_atoms; ++i) {
        const md_match_atom_t* a = &q->atoms[i];
        if (!(a->flags & MD_MATCH_ATOM_ELEMENT)) continue;
        if (a->z == MD_Z_H && degree[i] == 1) {
            const md_match_atom_t* b = &q->atoms[other[i]];
            if (!((b->flags & MD_MATCH_ATOM_ELEMENT) && b->z == MD_Z_H)) continue;
        }
        need[a->z]++;
    }
    md_temp_end(temp);

    uint32_t have[256] = {0};
    const md_atom_type_idx_t* type_idx = sys->atom.type_idx;
    const md_atomic_number_t* type_z   = sys->atom.type.z;
    if (type_idx && type_z) {
        for (size_t i = 0; i < sys->atom.count; ++i) {
            const md_atom_type_idx_t t = type_idx[i];
            have[t < sys->atom.type.count ? type_z[t] : 0]++;
        }
    } else {
        have[0] = (uint32_t)sys->atom.count;
    }
    for (int z = 0; z < 256; ++z) {
        if (need[z] > have[z]) return false;
    }
    return true;
}

static bool match_run(size_t* out_count, const md_match_desc_t* desc, const md_system_t* sys, md_match_callback_t callback, void* user_param) {
    if (out_count) *out_count = 0;
    if (!desc_validate(desc, sys)) return false;
    if (sys->atom.count == 0 || !elements_suffice(desc->query, sys)) return true;

    // The scratch memory lives in an arena of its own: the callback may allocate from any allocator, temp arenas
    // included, without interleaving with it.
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));

    target_t tg;
    target_build(&tg, sys, target_needs(desc->hydrogens, desc->charges, query_states_charge(desc->query)), arena);

    search_t s;
    search_init(&s, &tg, desc->level, desc->mode, desc->flags, desc->hydrogens, desc->bond_orders, desc->charges, desc->mask, arena);
    s.limit      = desc->limit;
    s.callback   = callback;
    s.user_param = user_param;

    uint32_t rarity[256];
    const uint32_t num_candidates = count_candidates(rarity, &s);

    plan_t plan;
    plan_build(&plan, desc->query, rarity, num_candidates, desc->mode != MD_MATCH_MODE_ALL, arena);

    if (!plan.impossible && plan.num <= num_candidates) {
        s.plan = &plan;
        search_reserve(&s, plan.num, plan.num_atoms, arena);
        s.uniq.width = plan.num;
        s.uniq.alloc = arena;

        const size_t num_units = level_num_units(&tg, desc->level);
        for (size_t u = 0; u < num_units && !s.stop; ++u) {
            unit_begin(&s, (int64_t)u);
            unit_search(&s, false);
        }
    }

    if (out_count) *out_count = s.count;
    md_vm_arena_destroy(arena);
    return true;
}

// ### PUBLIC ###

// The atoms of a SMILES graph which lie on a ring: the paths through the tree of its chain and branch bonds which
// each ring closure bond closes
static void smiles_ring_atoms(bool* in_ring, const md_smiles_t* g, md_allocator_i* alloc) {
    const size_t n = g->num_atoms;
    int32_t*  parent = md_alloc(alloc, sizeof(int32_t)  * n);
    uint32_t* depth  = md_alloc(alloc, sizeof(uint32_t) * n);
    for (size_t i = 0; i < n; ++i) {
        parent[i] = -1;
        in_ring[i] = false;
    }
    for (size_t i = 0; i < g->num_bonds; ++i) {
        const md_smiles_bond_t* b = &g->bonds[i];
        if (b->flags & MD_SMILES_BOND_RING) continue;
        // A chain or branch bond joins an atom to the one it follows, written before it
        const uint32_t lo = MIN(b->a, b->b), hi = MAX(b->a, b->b);
        parent[hi] = (int32_t)lo;
    }
    for (size_t i = 0; i < n; ++i) {
        depth[i] = parent[i] >= 0 ? depth[parent[i]] + 1 : 0;
    }
    for (size_t i = 0; i < g->num_bonds; ++i) {
        const md_smiles_bond_t* b = &g->bonds[i];
        if (!(b->flags & MD_SMILES_BOND_RING)) continue;
        // The lowest common ancestor, if the two are in the same tree (C1.C1 is not a ring)
        int32_t x = (int32_t)b->a, y = (int32_t)b->b;
        while (depth[x] > depth[y]) x = parent[x];
        while (depth[y] > depth[x]) y = parent[y];
        while (x != y && x >= 0 && y >= 0) {
            x = parent[x];
            y = parent[y];
        }
        if (x < 0 || y < 0 || x != y) continue;
        const int32_t top = x;
        for (int32_t v = (int32_t)b->a; v != top; v = parent[v]) in_ring[v] = true;
        for (int32_t v = (int32_t)b->b; v != top; v = parent[v]) in_ring[v] = true;
        in_ring[top] = true;
    }
}

bool md_match_query_init_smiles(md_match_query_t* query, str_t smiles, struct md_allocator_i* alloc, struct md_smiles_error_t* err) {
    ASSERT(query);
    ASSERT(alloc);
    MEMSET(query, 0, sizeof(md_match_query_t));

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    md_smiles_t g = {0};
    if (!md_smiles_parse(&g, smiles, temp_alloc, err)) {
        md_temp_end(temp);
        return false;
    }
    const size_t n = g.num_atoms;

    // Hydrogens written as atoms, which the search folds into the atom they are bonded to, add to the count of a
    // bracket atom: [C@@]([H]) has one hydrogen, like [C@@H]
    uint32_t* degree   = md_alloc(temp_alloc, sizeof(uint32_t) * n);
    uint32_t* h_atoms  = md_alloc(temp_alloc, sizeof(uint32_t) * n);
    bool any_aromatic = false;
    for (size_t i = 0; i < n; ++i) {
        degree[i]  = 0;
        h_atoms[i] = 0;
        any_aromatic |= (g.atoms[i].flags & MD_SMILES_ATOM_AROMATIC) != 0;
    }
    for (size_t i = 0; i < g.num_bonds; ++i) {
        degree[g.bonds[i].a]++;
        degree[g.bonds[i].b]++;
    }
    for (size_t i = 0; i < g.num_bonds; ++i) {
        const uint32_t a = g.bonds[i].a, b = g.bonds[i].b;
        if (g.atoms[b].z == MD_Z_H && degree[b] == 1 && g.atoms[a].z != MD_Z_H) h_atoms[a]++;
        if (g.atoms[a].z == MD_Z_H && degree[a] == 1 && g.atoms[b].z != MD_Z_H) h_atoms[b]++;
    }

    // Uppercase atoms outside of rings are aliphatic in a SMILES which writes aromatic atoms lowercase
    bool* in_ring = NULL;
    if (any_aromatic) {
        in_ring = md_alloc(temp_alloc, sizeof(bool) * n);
        smiles_ring_atoms(in_ring, &g, temp_alloc);
    }

    query->alloc     = alloc;
    query->num_atoms = n;
    query->atoms     = md_alloc(alloc, sizeof(md_match_atom_t) * n);
    MEMSET(query->atoms, 0, sizeof(md_match_atom_t) * n);
    for (size_t i = 0; i < n; ++i) {
        const md_smiles_atom_t* src = &g.atoms[i];
        md_match_atom_t* dst = &query->atoms[i];
        if (src->z != 0) {
            dst->flags |= MD_MATCH_ATOM_ELEMENT;
            dst->z = src->z;
        }
        if (src->flags & MD_SMILES_ATOM_AROMATIC) {
            dst->flags |= MD_MATCH_ATOM_AROMATIC;
        } else if (in_ring && !in_ring[i] && src->z != 0 && src->z != MD_Z_H) {
            dst->flags |= MD_MATCH_ATOM_ALIPHATIC;
        }
        // Bracket atoms state their hydrogens and charge, organic subset atoms leave them open
        if (src->flags & MD_SMILES_ATOM_BRACKET) {
            const uint32_t h = MIN((uint32_t)src->h_count + h_atoms[i], 254u);
            dst->flags |= MD_MATCH_ATOM_HCOUNT | MD_MATCH_ATOM_CHARGE;
            dst->h_min  = (uint8_t)h;
            dst->h_max  = (uint8_t)h;
            dst->charge = src->charge;
        }
        dst->source = (int32_t)src->offset;
        dst->tag    = src->atom_class;
    }

    if (g.num_bonds > 0) {
        query->num_bonds = g.num_bonds;
        query->bonds     = md_alloc(alloc, sizeof(md_match_bond_t) * g.num_bonds);
        for (size_t i = 0; i < g.num_bonds; ++i) {
            const md_smiles_bond_t* src = &g.bonds[i];
            uint32_t order = MD_MATCH_BOND_ANY;
            switch (src->order) {
            case MD_SMILES_BOND_IMPLICIT: {
                const bool aromatic = (g.atoms[src->a].flags & MD_SMILES_ATOM_AROMATIC) && (g.atoms[src->b].flags & MD_SMILES_ATOM_AROMATIC);
                order = aromatic ? MD_MATCH_BOND_AROMATIC : MD_MATCH_BOND_SINGLE;
                break;
            }
            case MD_SMILES_BOND_SINGLE:    order = MD_MATCH_BOND_SINGLE;    break;
            case MD_SMILES_BOND_DOUBLE:    order = MD_MATCH_BOND_DOUBLE;    break;
            case MD_SMILES_BOND_TRIPLE:    order = MD_MATCH_BOND_TRIPLE;    break;
            case MD_SMILES_BOND_QUADRUPLE: order = MD_MATCH_BOND_QUADRUPLE; break;
            case MD_SMILES_BOND_AROMATIC:  order = MD_MATCH_BOND_AROMATIC;  break;
            default: break;
            }
            query->bonds[i] = (md_match_bond_t){.a = src->a, .b = src->b, .order = order};
        }
    }

    md_temp_end(temp);
    return true;
}

bool md_match_query_init_atoms(md_match_query_t* query, const int32_t* atom_idx, size_t count, md_match_label_t label, const struct md_system_t* sys, struct md_allocator_i* alloc) {
    ASSERT(query);
    ASSERT(alloc);
    MEMSET(query, 0, sizeof(md_match_query_t));

    if (!sys) {
        MD_LOG_ERROR("md_match: no system");
        return false;
    }
    if (!atom_idx || count == 0) {
        MD_LOG_ERROR("md_match: no atoms given for the query");
        return false;
    }
    if (count >= INT32_MAX) {
        MD_LOG_ERROR("md_match: too many atoms given for the query");
        return false;
    }

    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));
    target_t tg;
    target_build(&tg, sys, TARGET_RESOLUTION, arena);

    // (atom << 32 | query index), sorted, to find the query index of an atom
    uint64_t* lookup = md_alloc(arena, sizeof(uint64_t) * count);
    bool has_h = false;
    for (size_t i = 0; i < count; ++i) {
        const int32_t a = atom_idx[i];
        if (a < 0 || (size_t)a >= tg.num_atoms) {
            MD_LOG_ERROR("md_match: atom index %d out of range", a);
            md_vm_arena_destroy(arena);
            return false;
        }
        lookup[i] = ((uint64_t)(uint32_t)a << 32) | (uint64_t)i;
        has_h |= target_z(&tg, a) == MD_Z_H;
    }
    qsort(lookup, count, sizeof(uint64_t), cmp_u64);
    for (size_t i = 1; i < count; ++i) {
        if ((lookup[i] >> 32) == (lookup[i - 1] >> 32)) {
            MD_LOG_ERROR("md_match: atom %u given twice", (uint32_t)(lookup[i] >> 32));
            md_vm_arena_destroy(arena);
            return false;
        }
    }

    #define FIND_QUERY_IDX(out, atom) do {                                      \
        size_t lo_ = 0, hi_ = count;                                            \
        (out) = -1;                                                             \
        while (lo_ < hi_) {                                                     \
            const size_t mid_ = (lo_ + hi_) / 2;                                \
            const uint32_t v_ = (uint32_t)(lookup[mid_] >> 32);                 \
            if (v_ < (uint32_t)(atom)) lo_ = mid_ + 1;                          \
            else if (v_ > (uint32_t)(atom)) hi_ = mid_;                         \
            else { (out) = (int64_t)(uint32_t)lookup[mid_]; break; }            \
        }                                                                       \
    } while (0)

    query->alloc     = alloc;
    query->num_atoms = count;
    query->atoms     = md_alloc(alloc, sizeof(md_match_atom_t) * count);
    MEMSET(query->atoms, 0, sizeof(md_match_atom_t) * count);

    md_array(md_match_bond_t) bonds = 0;
    for (size_t i = 0; i < count; ++i) {
        const int32_t a = atom_idx[i];
        const md_atomic_number_t z = target_z(&tg, a);
        md_match_atom_t* dst = &query->atoms[i];
        dst->flags  = MD_MATCH_ATOM_ELEMENT;
        dst->z      = z;
        dst->source = a;
        if (label == MD_MATCH_LABEL_NAME) {
            dst->flags |= MD_MATCH_ATOM_NAME;
            dst->name   = target_name(&tg, a);
        }

        if (target_aromatic(&tg, a) == 1) {
            dst->flags |= MD_MATCH_ATOM_AROMATIC;
        }

        // An atom whose hydrogens are all given has its protonation fixed, see md_match.h. With hydrogen counts in
        // the system every atom resolves its hydrogens. Virtual sites are in no molecule and have none.
        const bool resolved = !target_skipped(&tg, a) && (tg.h_count || resolution_h(molecule_resolution(&tg, target_molecule(&tg, a)), z));
        if (has_h && z != MD_Z_H && resolved) {
            uint32_t h_given = 0;
            for (uint32_t c = conn_beg(&tg, a); c < conn_end(&tg, a); ++c) {
                const int32_t b = tg.conn_atom[c];
                if (target_z(&tg, b) != MD_Z_H) continue;
                int64_t j;
                FIND_QUERY_IDX(j, b);
                h_given += j >= 0;
            }
            const uint32_t h_total = target_hydrogens(&tg, a);
            if (h_given == h_total) {
                dst->flags |= MD_MATCH_ATOM_HCOUNT;
                dst->h_min = (uint8_t)MIN(h_total, 254);
                dst->h_max = (uint8_t)MIN(h_total, 254);
                if (tg.charge) {
                    dst->flags |= MD_MATCH_ATOM_CHARGE;
                    dst->charge = tg.charge[a];
                }
            }
        }

        for (uint32_t c = conn_beg(&tg, a); c < conn_end(&tg, a); ++c) {
            const int32_t b = tg.conn_atom[c];
            if (b <= a) continue;
            int64_t j;
            FIND_QUERY_IDX(j, b);
            if (j < 0) continue;
            const md_bond_flags_t f = target_bond_flags(&tg, tg.conn_bond[c]);
            uint32_t order = MD_MATCH_BOND_ANY;
            if (f & (MD_BOND_FLAG_AROMATIC | MD_BOND_FLAG_DELOCALIZED)) {
                order = MD_MATCH_BOND_AROMATIC;
            } else {
                switch (md_bond_order(f)) {
                case MD_BOND_ORDER_SINGLE:    order = MD_MATCH_BOND_SINGLE;    break;
                case MD_BOND_ORDER_DOUBLE:    order = MD_MATCH_BOND_DOUBLE;    break;
                case MD_BOND_ORDER_TRIPLE:    order = MD_MATCH_BOND_TRIPLE;    break;
                case MD_BOND_ORDER_QUADRUPLE: order = MD_MATCH_BOND_QUADRUPLE; break;
                default: break;
                }
            }
            md_match_bond_t bond = {.a = (uint32_t)i, .b = (uint32_t)j, .order = order};
            md_array_push(bonds, bond, arena);
        }
    }
    #undef FIND_QUERY_IDX

    const size_t num_bonds = md_array_size(bonds);
    if (num_bonds > 0) {
        query->num_bonds = num_bonds;
        query->bonds = md_alloc(alloc, sizeof(md_match_bond_t) * num_bonds);
        MEMCPY(query->bonds, bonds, sizeof(md_match_bond_t) * num_bonds);
    }

    // A reference in several parts is matched part by part within one unit (a chain with a gap in it)
    if (count > 1) {
        int32_t* parent = md_alloc(arena, sizeof(int32_t) * count);
        for (size_t i = 0; i < count; ++i) parent[i] = (int32_t)i;
        for (size_t i = 0; i < num_bonds; ++i) uf_union(parent, (int32_t)bonds[i].a, (int32_t)bonds[i].b);
        size_t num_parts = 0;
        for (size_t i = 0; i < count; ++i) num_parts += parent[i] == (int32_t)i;
        if (num_parts > 1) {
            MD_LOG_DEBUG("md_match: the reference atoms are in %zu parts, which are matched within one unit", num_parts);
        }
    }

    md_vm_arena_destroy(arena);
    return true;
}

void md_match_query_free(md_match_query_t* query) {
    ASSERT(query);
    if (query->alloc) {
        if (query->atoms) md_free(query->alloc, query->atoms, sizeof(md_match_atom_t) * query->num_atoms);
        if (query->bonds) md_free(query->alloc, query->bonds, sizeof(md_match_bond_t) * query->num_bonds);
    }
    MEMSET(query, 0, sizeof(md_match_query_t));
}

bool md_match_for_each(size_t* out_count, const md_match_desc_t* desc, const struct md_system_t* sys, md_match_callback_t callback, void* user_param) {
    return match_run(out_count, desc, sys, callback, user_param);
}

typedef struct collect_t {
    md_array(int32_t)  atom_idx;
    md_array(uint32_t) unit;
    md_allocator_i*    alloc;
} collect_t;

static bool collect_callback(const int32_t* atom_idx, size_t width, uint32_t unit, void* user_param) {
    collect_t* c = (collect_t*)user_param;
    md_array_push_array(c->atom_idx, atom_idx, width, c->alloc);
    md_array_push(c->unit, unit, c->alloc);
    return true;
}

bool md_match_find(md_match_result_t* out, const md_match_desc_t* desc, const struct md_system_t* sys, struct md_allocator_i* alloc) {
    ASSERT(out);
    ASSERT(alloc);
    MEMSET(out, 0, sizeof(md_match_result_t));

    collect_t c = {.alloc = alloc};
    size_t count = 0;
    if (!match_run(&count, desc, sys, collect_callback, &c)) {
        md_array_free(c.atom_idx, alloc);
        md_array_free(c.unit, alloc);
        return false;
    }

    out->count    = count;
    out->width    = desc->query->num_atoms;
    out->atom_idx = c.atom_idx;
    out->unit     = c.unit;
    out->alloc    = alloc;
    return true;
}

void md_match_result_free(md_match_result_t* result) {
    ASSERT(result);
    if (result->alloc) {
        md_array_free(result->atom_idx, result->alloc);
        md_array_free(result->unit, result->alloc);
    }
    MEMSET(result, 0, sizeof(md_match_result_t));
}

// The atoms of one match, with the hydrogens of its atoms when asked for
static void select_row(md_bitfield_t* bf, const int32_t* row, size_t width, const md_system_t* sys, bool hydrogens) {
    const md_bond_conn_data_t* conn = sys ? &sys->bond.conn : NULL;
    hydrogens = hydrogens && conn && conn->offset && conn->offset_count >= sys->atom.count + 1;
    for (size_t j = 0; j < width; ++j) {
        const int32_t a = row[j];
        if (a < 0) continue;
        md_bitfield_set_bit(bf, (uint64_t)a);
        if (hydrogens && md_atom_atomic_number(&sys->atom, a) != MD_Z_H) {
            for (uint32_t c = conn->offset[a]; c < conn->offset[a + 1]; ++c) {
                const int32_t b = conn->atom_idx[c];
                if (md_atom_atomic_number(&sys->atom, b) == MD_Z_H) md_bitfield_set_bit(bf, (uint64_t)b);
            }
        }
    }
}

size_t md_match_result_select(struct md_bitfield_t* out, size_t cap, const md_match_result_t* result, const struct md_system_t* sys, uint32_t flags) {
    ASSERT(result);
    if (!out || cap == 0 || result->count == 0) return 0;
    const bool hydrogens = (flags & MD_MATCH_SELECT_HYDROGENS) != 0;
    const size_t n = MIN(cap, result->count);
    for (size_t i = 0; i < n; ++i) {
        md_bitfield_clear(&out[i]);
        select_row(&out[i], result->atom_idx + i * result->width, result->width, sys, hydrogens);
    }
    return n;
}

void md_match_result_select_all(struct md_bitfield_t* out, const md_match_result_t* result, const struct md_system_t* sys, uint32_t flags) {
    ASSERT(out);
    ASSERT(result);
    md_bitfield_clear(out);
    const bool hydrogens = (flags & MD_MATCH_SELECT_HYDROGENS) != 0;
    for (size_t i = 0; i < result->count; ++i) {
        select_row(out, result->atom_idx + i * result->width, result->width, sys, hydrogens);
    }
}

// ### LIBRARY ###

typedef struct lib_entry_t {
    md_match_query_t query;
    str_t    name;
    // What the prefilter of md_match_identify compares: the heavy atoms of the query, their elements and bonds
    uint64_t key;               // Sum of the element hashes of its heavy atoms, and of its bonds
    uint32_t num_heavy;
    uint32_t num_heavy_bonds;
    uint32_t comp_off;          // Its (z, count) pairs in md_match_library_t.comp
    uint32_t comp_len;
    bool     wild;              // A heavy atom without an element
} lib_entry_t;

struct md_match_library_t {
    md_allocator_i* alloc;
    md_array(lib_entry_t) entries;
    md_array(uint32_t)    comp;     // z << 16 | count
};

static inline uint64_t element_hash(md_atomic_number_t z) {
    // splitmix64 of the element: a multiset of elements hashes to the sum of its members
    uint64_t x = (uint64_t)z * 0x9E3779B97F4A7C15ULL + 0x632BE59BD9B4E019ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

#define BOND_HASH 0xD6E8FEB86659FD93ULL

md_match_library_t* md_match_library_create(struct md_allocator_i* alloc) {
    ASSERT(alloc);
    md_match_library_t* lib = md_alloc(alloc, sizeof(md_match_library_t));
    MEMSET(lib, 0, sizeof(md_match_library_t));
    lib->alloc = alloc;
    return lib;
}

void md_match_library_destroy(md_match_library_t* lib) {
    if (!lib) return;
    md_allocator_i* alloc = lib->alloc;
    for (size_t i = 0; i < md_array_size(lib->entries); ++i) {
        md_match_query_free(&lib->entries[i].query);
        if (lib->entries[i].name.ptr) str_free(lib->entries[i].name, alloc);
    }
    md_array_free(lib->entries, alloc);
    md_array_free(lib->comp, alloc);
    md_free(alloc, lib, sizeof(md_match_library_t));
}

// Takes the query over
static int32_t library_push(md_match_library_t* lib, md_match_query_t query, str_t name) {
    lib_entry_t e = {0};
    e.query = query;
    e.name  = str_empty(name) ? (str_t){0} : str_copy(name, lib->alloc);

    const size_t n = query.num_atoms;
    md_temp_scope_t temp = md_temp_begin_avoid(lib->alloc);
    uint32_t* degree = md_temp_alloc_array(temp, uint32_t, n);
    uint32_t  count[256] = {0};
    MEMSET(degree, 0, sizeof(uint32_t) * n);
    for (size_t i = 0; i < query.num_bonds; ++i) {
        degree[query.bonds[i].a]++;
        degree[query.bonds[i].b]++;
    }
    for (size_t i = 0; i < n; ++i) {
        const md_match_atom_t* a = &query.atoms[i];
        const bool is_h = (a->flags & MD_MATCH_ATOM_ELEMENT) && a->z == MD_Z_H;
        if (is_h) continue;
        e.num_heavy += 1;
        if (a->flags & MD_MATCH_ATOM_ELEMENT) {
            count[a->z] += 1;
            e.key += element_hash(a->z);
        } else {
            e.wild = true;
        }
    }
    for (size_t i = 0; i < query.num_bonds; ++i) {
        const md_match_atom_t* a = &query.atoms[query.bonds[i].a];
        const md_match_atom_t* b = &query.atoms[query.bonds[i].b];
        const bool a_h = (a->flags & MD_MATCH_ATOM_ELEMENT) && a->z == MD_Z_H;
        const bool b_h = (b->flags & MD_MATCH_ATOM_ELEMENT) && b->z == MD_Z_H;
        if (!a_h && !b_h) e.num_heavy_bonds += 1;
    }
    md_temp_end(temp);
    e.key += (uint64_t)e.num_heavy_bonds * BOND_HASH;

    e.comp_off = (uint32_t)md_array_size(lib->comp);
    for (uint32_t z = 0; z < 256; ++z) {
        if (count[z]) md_array_push(lib->comp, z << 16 | MIN(count[z], 0xFFFFu), lib->alloc);
    }
    e.comp_len = (uint32_t)md_array_size(lib->comp) - e.comp_off;

    md_array_push(lib->entries, e, lib->alloc);
    return (int32_t)md_array_size(lib->entries) - 1;
}

int32_t md_match_library_add(md_match_library_t* lib, const md_match_query_t* query, str_t name) {
    ASSERT(lib);
    if (!query_validate(query)) return -1;
    md_match_query_t copy = {0};
    copy.alloc     = lib->alloc;
    copy.num_atoms = query->num_atoms;
    copy.atoms     = md_alloc(lib->alloc, sizeof(md_match_atom_t) * query->num_atoms);
    MEMCPY(copy.atoms, query->atoms, sizeof(md_match_atom_t) * query->num_atoms);
    if (query->num_bonds > 0) {
        copy.num_bonds = query->num_bonds;
        copy.bonds     = md_alloc(lib->alloc, sizeof(md_match_bond_t) * query->num_bonds);
        MEMCPY(copy.bonds, query->bonds, sizeof(md_match_bond_t) * query->num_bonds);
    }
    return library_push(lib, copy, name);
}

int32_t md_match_library_add_smiles(md_match_library_t* lib, str_t smiles, str_t name, struct md_smiles_error_t* err) {
    ASSERT(lib);
    md_match_query_t query = {0};
    if (!md_match_query_init_smiles(&query, smiles, lib->alloc, err)) return -1;
    return library_push(lib, query, name);
}

size_t md_match_library_count(const md_match_library_t* lib) {
    return lib ? md_array_size(lib->entries) : 0;
}

const md_match_query_t* md_match_library_query(const md_match_library_t* lib, size_t entry) {
    if (!lib || entry >= md_array_size(lib->entries)) return NULL;
    return &lib->entries[entry].query;
}

str_t md_match_library_name(const md_match_library_t* lib, size_t entry) {
    if (!lib || entry >= md_array_size(lib->entries)) return (str_t){0};
    return lib->entries[entry].name;
}

typedef struct keyed_entry_t {
    uint64_t key;
    uint32_t entry;
} keyed_entry_t;

static int cmp_keyed_entry(const void* a, const void* b) {
    const keyed_entry_t* x = (const keyed_entry_t*)a;
    const keyed_entry_t* y = (const keyed_entry_t*)b;
    if (x->key != y->key) return x->key < y->key ? -1 : 1;
    return (x->entry > y->entry) - (x->entry < y->entry);
}

typedef struct identify_row_t {
    int32_t* row;
    size_t   width;
    bool     found;
} identify_row_t;

static bool identify_callback(const int32_t* atom_idx, size_t width, uint32_t unit, void* user_param) {
    (void)unit;
    identify_row_t* r = (identify_row_t*)user_param;
    MEMCPY(r->row, atom_idx, sizeof(int32_t) * width);
    r->width = width;
    r->found = true;
    return true;
}

bool md_match_identify(md_match_identify_result_t* out, const md_match_identify_desc_t* desc, const struct md_system_t* sys, struct md_allocator_i* alloc) {
    ASSERT(out);
    ASSERT(alloc);
    MEMSET(out, 0, sizeof(md_match_identify_result_t));
    out->alloc = alloc;

    if (!desc || !desc->library) {
        MD_LOG_ERROR("md_match: no library to identify with");
        return false;
    }
    if (!sys) {
        MD_LOG_ERROR("md_match: no system");
        return false;
    }
    if (!level_valid(desc->level, sys)) return false;
    if (!resolve_valid(desc->hydrogens) || !resolve_valid(desc->bond_orders) || !resolve_valid(desc->charges)) {
        MD_LOG_ERROR("md_match: unknown resolution");
        return false;
    }

    md_array_push(out->offset, 0, alloc);
    const md_match_library_t* lib = desc->library;
    const size_t num_entries = md_array_size(lib->entries);
    if (sys->atom.count == 0 || num_entries == 0) return true;

    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));

    bool has_charge = false;
    for (size_t e = 0; e < num_entries; ++e) has_charge |= query_states_charge(&lib->entries[e].query);

    target_t tg;
    target_build(&tg, sys, target_needs(desc->hydrogens, desc->charges, has_charge), arena);

    const bool whole = (desc->flags & MD_MATCH_FLAG_WHOLE) != 0;
    search_t s;
    search_init(&s, &tg, desc->level, MD_MATCH_MODE_ONE_PER_UNIT, desc->flags, desc->hydrogens, desc->bond_orders, desc->charges, desc->mask, arena);

    uint32_t rarity[256];
    const uint32_t num_candidates = count_candidates(rarity, &s);

    // A plan per entry
    plan_t* plans = md_alloc(arena, sizeof(plan_t) * num_entries);
    uint32_t max_num = 1;
    size_t max_atoms = 1;
    for (size_t e = 0; e < num_entries; ++e) {
        plan_build(&plans[e], &lib->entries[e].query, rarity, num_candidates, true, arena);
        max_num   = MAX(max_num, plans[e].num);
        max_atoms = MAX(max_atoms, plans[e].num_atoms);
    }
    search_reserve(&s, max_num, max_atoms, arena);
    identify_row_t found = {.row = md_alloc(arena, sizeof(int32_t) * max_atoms)};
    s.callback   = identify_callback;
    s.user_param = &found;

    // The prefilter. With WHOLE an entry without wildcards fits a unit of exactly its composition: those are looked up
    // by key. Every other entry is compared element by element.
    keyed_entry_t* keyed = md_alloc(arena, sizeof(keyed_entry_t) * num_entries);
    uint32_t*      other = md_alloc(arena, sizeof(uint32_t) * num_entries);
    size_t num_keyed = 0, num_other = 0;
    for (size_t e = 0; e < num_entries; ++e) {
        if (plans[e].impossible || plans[e].num > num_candidates) continue;
        if (whole && !lib->entries[e].wild) {
            keyed[num_keyed++] = (keyed_entry_t){.key = lib->entries[e].key, .entry = (uint32_t)e};
        } else {
            other[num_other++] = (uint32_t)e;
        }
    }
    qsort(keyed, num_keyed, sizeof(keyed_entry_t), cmp_keyed_entry);

    uint32_t  count[256] = {0};
    uint8_t   touched[256];
    uint32_t* cand = md_alloc(arena, sizeof(uint32_t) * num_entries);

    const size_t num_units = level_num_units(&tg, desc->level);
    for (size_t u = 0; u < num_units; ++u) {
        unit_begin(&s, (int64_t)u);
        if (s.unit_count == 0) continue;

        // Composition of the unit
        uint32_t num_heavy = 0, num_bonds = 0, num_touched = 0;
        uint64_t key = 0;
        for (uint32_t i = 0; i < s.unit_count; ++i) {
            const int32_t a = unit_atom(&s, i);
            if (target_skipped(&tg, a)) continue;
            if (s.mask && !bit_test(s.mask, a)) continue;
            const md_atomic_number_t z = target_z(&tg, a);
            if (z == MD_Z_H) continue;
            num_heavy++;
            key += element_hash(z);
            if (count[z]++ == 0) touched[num_touched++] = z;
            for (uint32_t c = conn_beg(&tg, a); c < conn_end(&tg, a); ++c) {
                const int32_t b = tg.conn_atom[c];
                if (b <= a || target_z(&tg, b) == MD_Z_H || !in_graph(&s, b)) continue;
                num_bonds++;
            }
        }
        key += (uint64_t)num_bonds * BOND_HASH;

        // The entries which fit, in the order of the library
        size_t num_cand = 0;
        {
            // The keyed entries of this key, merged with the others
            size_t lo = 0, hi = num_keyed;
            while (lo < hi) {
                const size_t mid = (lo + hi) / 2;
                if (keyed[mid].key < key) lo = mid + 1;
                else hi = mid;
            }
            size_t k = lo;
            size_t o = 0;
            for (;;) {
                const bool has_k = k < num_keyed && keyed[k].key == key;
                const bool has_o = o < num_other;
                if (!has_k && !has_o) break;
                uint32_t e;
                if (has_k && (!has_o || keyed[k].entry < other[o])) {
                    e = keyed[k].entry;
                    k++;
                    const lib_entry_t* le = &lib->entries[e];
                    if (le->num_heavy != num_heavy || le->num_heavy_bonds != num_bonds) continue;
                } else {
                    e = other[o++];
                    const lib_entry_t* le = &lib->entries[e];
                    if (whole ? (le->num_heavy != num_heavy || le->num_heavy_bonds != num_bonds) : (le->num_heavy > num_heavy || le->num_heavy_bonds > num_bonds)) continue;
                    bool fits = true;
                    for (uint32_t i = 0; i < le->comp_len && fits; ++i) {
                        const uint32_t pair = lib->comp[le->comp_off + i];
                        fits = count[pair >> 16] >= (pair & 0xFFFF);
                    }
                    if (!fits) continue;
                }
                cand[num_cand++] = e;
            }
        }
        for (uint32_t i = 0; i < num_touched; ++i) count[touched[i]] = 0;

        for (size_t i = 0; i < num_cand; ++i) {
            const uint32_t e = cand[i];
            s.plan = &plans[e];
            found.found = false;
            unit_search(&s, true);
            if (found.found) {
                md_array_push(out->unit,  (uint32_t)u, alloc);
                md_array_push(out->entry, e, alloc);
                md_array_push_array(out->atom_idx, found.row, found.width, alloc);
                md_array_push(out->offset, (uint32_t)md_array_size(out->atom_idx), alloc);
                out->count++;
                break;
            }
        }
    }

    md_vm_arena_destroy(arena);
    return true;
}

void md_match_identify_result_free(md_match_identify_result_t* result) {
    ASSERT(result);
    if (result->alloc) {
        md_array_free(result->unit, result->alloc);
        md_array_free(result->entry, result->alloc);
        md_array_free(result->offset, result->alloc);
        md_array_free(result->atom_idx, result->alloc);
    }
    MEMSET(result, 0, sizeof(md_match_identify_result_t));
}
