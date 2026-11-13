#include "md_infer_internal.h"
#include <core/md_arena_allocator.h>

// Lipid inference: the lipids of a system, their heads and tails (md_util.h).

// ### LIPIDS ###
//
// md_lipid_infer: the heads and tails of the lipids (md_lipid_data_t). A molecule is a set of residues
// which are not amino acids, nucleotides, water or ions, bonded to each other, wherever they are in the system. Its apolar atoms
// are the carbons bonded to carbons (and hydrogens) only, or, for coarse grained beads, the beads named as the beads of
// the tails of the Martini lipidomes (lipid_bead_kind). Each connected part of them which hangs off the rest of the
// molecule at one atom (its anchor) is a tail: a chain, or the skeleton of a sterol when it holds fused rings (or the
// ring beads of a coarse grained sterol). A chain which hangs between two anchors spans them (a bolalipid). The rest of the molecule, in connected parts, is its heads.
//
// A molecule whose residues are lipids by name is a lipid with the tails found. One of no lipid is a lipid by its
// structure: a polar head (an atom other than carbon outside of the tails) and a chain without rings of at least
// LIPID_MIN_CHAIN_UNNAMED carbons, or the skeleton of a sterol. Coarse grained molecules are lipids by the force field
// only: by their names, by their beads where the name is not Martini's (md_system_infer_residue_flags), or by beads
// named as those of the tails of the Martini lipidomes where the molecule is none of the force field's.
//
// The lipid of a molecule of several residues is the part of it which holds the tails, with what is bonded to that
// part (lipid_extent): a sugar chain or a polymer attached to a lipid, each sugar or monomer a residue of its own, is no
// part of the lipid.

#define LIPID_MIN_CHAIN_NAMED       4   // Carbons on the spine of a chain of a lipid by name (hexanoyl: 5)
#define LIPID_MIN_CHAIN_UNNAMED     8   // Of a molecule which is a lipid by its structure alone (octanoyl: 7, decanoyl: 9)
#define LIPID_MIN_CHAIN_BEADS       1   // Beads on the spine of a chain of a coarse grained lipid (the C1B of a lyso PC)
#define LIPID_MIN_CHAIN_LETTERLESS  2   // Of a chain of beads named without the letter of their tail (COO C1 C2 C3 C4)
#define LIPID_MIN_SKELETON_ATOMS    15  // Carbons of the skeleton of a sterol (cholesterol: 26 besides C3)
#define LIPID_MIN_SKELETON_BEADS    4   // Beads of it (Martini cholesterol: 7 and 8)
#define LIPID_MAX_ATOMS             16384   // Larger molecules (polymers) are no lipids

// What a bead of a lipid of the Martini lipidomes is to its tails, by its name:
//   CHAIN       A bead of a tail: C (saturated), D (unsaturated) or T (trans), the number of the bead and the letter of
//               its tail, with the number of the acyl chain after it in cardiolipins (C1A, D2A, C4B, T1A, C1A1); the ring
//               skeleton of a sterol R and a number (R1 .. R6), and in a sterol its side chain C and a number (C1, C2)
//               and its dummies V (VOH V1 V2 V3): the optimal LINCS cholesterol of Martini 2 carries the mass of its
//               rings on dummies which do not interact, its ring beads virtual sites of them, so its skeleton holds
//               together through them, and its head is the hydroxyl bead ROH alone, as that of the cholesterol of 8 beads.
//   LETTERLESS  C, D or T and a number alone elsewhere: the chain of a fatty acid or a detergent of Martini 2 (COO C1 C2
//               C3 C4, SO3 C1 C2 C4, the B1 B2 B3 C1 C2 of an octyl glucoside), but also the inositol of a
//               phosphatidylinositol (C1 .. C4) and the sugars of glycolipids. Such beads are a tail only where they
//               make a chain: an unbranched run of two beads or more, without rings, which the inositol and sugars are.
//   POLAR       Any other: the beads of the heads (NC3, NH3, CNO, PO4, GL1, GL2, AM1, ROH, COO, ...).
typedef enum {
    LIPID_BEAD_POLAR = 0,
    LIPID_BEAD_CHAIN,
    LIPID_BEAD_LETTERLESS,
} lipid_bead_t;

static lipid_bead_t lipid_bead_kind(md_str_t name, bool sterol) {
    if (name.len < 2) return LIPID_BEAD_POLAR;
    const char c = name.ptr[0];
    if (c == 'V' && sterol) return LIPID_BEAD_CHAIN;
    if (c != 'C' && c != 'D' && c != 'T' && c != 'R') return LIPID_BEAD_POLAR;
    size_t i = 1;
    if (!md_ascii_is_digit(name.ptr[i])) return LIPID_BEAD_POLAR;
    while (i < name.len && md_ascii_is_digit(name.ptr[i])) ++i;
    if (i == name.len) {
        if (c == 'R' || (c == 'C' && sterol)) return LIPID_BEAD_CHAIN;
        return LIPID_BEAD_LETTERLESS;
    }
    if (name.ptr[i] < 'A' || name.ptr[i] > 'H') return LIPID_BEAD_POLAR;
    ++i;
    while (i < name.len && md_ascii_is_digit(name.ptr[i])) ++i;
    return i == name.len ? LIPID_BEAD_CHAIN : LIPID_BEAD_POLAR;
}

static inline bool lipid_atom_is_bead(const md_system_t* sys, size_t i) {
    return md_atom_particle_kind(&sys->atom, i) != MD_PARTICLE_ATOM && md_atom_atomic_number(&sys->atom, i) == 0;
}

// A heavy atom, or a bead: what a head or the spine of a tail is made of
static inline bool lipid_atom_is_heavy(const md_system_t* sys, size_t i) {
    return lipid_atom_is_bead(sys, i) || md_atom_atomic_number(&sys->atom, i) > MD_Z_H;
}

// A ring bead of a coarse grained sterol (R1 .. R6)
static inline bool lipid_bead_sterol_ring(md_str_t name) {
    return name.len >= 2 && name.ptr[0] == 'R' && md_ascii_is_digit(name.ptr[1]);
}

// What an atom is to the tails: a bead by its name, an atom a CHAIN when it is a carbon bonded to carbons and hydrogens only
static lipid_bead_t lipid_atom_kind(const md_system_t* sys, size_t i, bool sterol) {
    if (lipid_atom_is_bead(sys, i)) {
        return lipid_bead_kind(md_atom_name(&sys->atom, i), sterol);
    }
    if (md_atom_atomic_number(&sys->atom, i) != MD_Z_C) return LIPID_BEAD_POLAR;
    md_bond_iter_t it = md_bond_iter(&sys->bond, i);
    while (md_bond_iter_has_next(&it)) {
        const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, md_bond_iter_atom_index(&it));
        if (z != MD_Z_C && z != MD_Z_H) return LIPID_BEAD_POLAR;
        md_bond_iter_next(&it);
    }
    return LIPID_BEAD_CHAIN;
}

static inline bool lipid_candidate_res(const md_system_t* sys, size_t res_idx) {
    const md_residue_kind_t kind = md_system_residue_kind(sys, res_idx);
    return kind == MD_RESIDUE_KIND_OTHER || kind == MD_RESIDUE_KIND_LIPID;
}

// A bead named as one of a tail with the letter of its tail (C1A, D2B, T1A, C1A1), as only the lipids of the Martini
// lipidomes have them
static inline bool lipid_bead_lettered(md_str_t name) {
    return name.len >= 3 && name.ptr[0] != 'R' && lipid_bead_kind(name, false) == LIPID_BEAD_CHAIN;
}

// The residues of the molecules are joined by the bonds between them (union-find, each set rooted at its first
// residue)
static int32_t lipid_find(int32_t* root, int32_t c) {
    while (root[c] != c) {
        root[c] = root[root[c]];
        c = root[c];
    }
    return c;
}

static void lipid_union(int32_t* root, int32_t a, int32_t b) {
    a = lipid_find(root, a);
    b = lipid_find(root, b);
    if (a < b) root[b] = a;
    else if (b < a) root[a] = b;
}

// A molecule: residues bonded to each other, ascending, and their atoms. Its residues need not be consecutive (the
// lipid A of a PDB entry, whose acyl chains are ligands of their own listed after its sugars).
typedef struct {
    const int32_t* residues; // [num_residues] ascending
    int32_t  num_residues;
    int32_t* atoms;         // [n] ascending: the atoms of its residues one after the other
    int32_t* member;        // [n] the residue (index into residues) of each atom
    int32_t  n;
    bool     contiguous;    // atoms[k] == atoms[0] + k
} lipid_mol_t;

// The index within the molecule of an atom, -1 for one not of it
static inline int32_t lipid_mol_local(const lipid_mol_t* mol, int32_t atom) {
    if (mol->contiguous) {
        const int32_t k = atom - mol->atoms[0];
        return (0 <= k && k < mol->n) ? k : -1;
    }
    int32_t lo = 0, hi = mol->n;
    while (lo < hi) {
        const int32_t mid = (lo + hi) / 2;
        if (mol->atoms[mid] < atom) lo = mid + 1;
        else hi = mid;
    }
    return (lo < mol->n && mol->atoms[lo] == atom) ? lo : -1;
}

typedef struct {
    int32_t  beg, end;      // its atoms (local indices) in part_atoms
    int32_t  anchor[2];     // atoms (global) of the rest of the molecule it hangs off
    int32_t  num_anchors;   // distinct anchors, more than 2 counts as 3
    bool     outside;       // hangs off an atom outside of the molecule
    bool     skeleton;      // holds an atom of a fused ring, or a ring bead of a sterol (R1 .. R6)
    bool     letterless;    // all of its beads are named without the letter of a tail (LIPID_BEAD_LETTERLESS)
    bool     lettered;      // all of its atoms are beads named with the letter of a tail (lipid_bead_lettered)
    bool     branched;      // an atom of it is bonded to more than two others of it
    int32_t  bonds;         // bonds within it, each counted from both ends
} lipid_part_t;

static int lipid_cmp_i32(const void* a, const void* b) {
    const int32_t x = *(const int32_t*)a, y = *(const int32_t*)b;
    return (x > y) - (x < y);
}

static int lipid_cmp_u64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
    return (x > y) - (x < y);
}

// Appends a range of indices to a ragged array of offsets (the first of them 0) and indices
static void push_index_range(uint32_t** offset, int32_t** idx, const int32_t* src, size_t n, md_allocator_i* alloc) {
    if (md_array_size(*offset) == 0) md_array_push(*offset, 0, alloc);
    if (n > 0) md_array_push_array(*idx, src, n, alloc);
    md_array_push(*offset, (uint32_t)md_array_size(*idx), alloc);
}

// The empty lipid data, whose offsets start at 0
static void init_lipid_data(md_lipid_data_t* lipid, md_allocator_i* alloc) {
    md_array_push(lipid->residue_offset, 0, alloc);
    md_array_push(lipid->head_offset, 0, alloc);
    md_array_push(lipid->tail_offset, 0, alloc);
    md_array_push(lipid->head.offset, 0, alloc);
    md_array_push(lipid->tail.spine_offset, 0, alloc);
    md_array_push(lipid->tail.atom_offset, 0, alloc);
}

// Records a lipid of the residues residues (ascending), with the heads and tails pushed for it since the last one (none
// when its anatomy is not known)
static void record_lipid(md_lipid_data_t* lipid, const int32_t* residues, size_t num_residues, md_allocator_i* alloc) {
    push_index_range(&lipid->residue_offset, &lipid->residue_idx, residues, num_residues, alloc);
    md_array_push(lipid->head_offset, (uint32_t)lipid->head.count, alloc);
    md_array_push(lipid->tail_offset, (uint32_t)lipid->tail.count, alloc);
    lipid->count += 1;
}

// The lipid of a molecule of several residues: the residues which hold its tails and those between them (the core),
// and the residues bonded to the core: a head group which is a residue of its own (the PC of a lyso PC of the AMBER
// lipids), the first sugar of a glycolipid, the first Kdo of a lipid A. What hangs off that (the further sugars of a
// glycolipid, the core and O antigen of a lipopolysaccharide, the polymer of a PEG lipid of several residues) is attached
// to the lipid but not of it. The core is what is left of the residues when those without tails are taken off the ends
// one after the other. in_lipid[c] for each residue of the molecule.
static void lipid_extent(uint8_t* in_lipid, const md_system_t* sys, const lipid_mol_t* mol, const uint8_t* has_tail, md_allocator_i* temp_arena) {
    const int32_t nc = mol->num_residues;
    md_memset(in_lipid, 1, (size_t)nc);
    if (nc < 2) return;

    // The bonds between the residues, each pair once
    md_array(uint64_t) edges = 0;
    for (int32_t k = 0; k < mol->n; ++k) {
        md_bond_iter_t it = md_bond_iter(&sys->bond, mol->atoms[k]);
        while (md_bond_iter_has_next(&it)) {
            const int32_t l = lipid_mol_local(mol, md_bond_iter_atom_index(&it));
            md_bond_iter_next(&it);
            if (l < 0) continue;
            const uint32_t a = (uint32_t)mol->member[k], b = (uint32_t)mol->member[l];
            if (a < b) md_array_push(edges, ((uint64_t)a << 32) | b, temp_arena);
        }
    }
    size_t num_edges = md_array_size(edges);
    if (num_edges == 0) return;
    qsort(edges, num_edges, sizeof(uint64_t), lipid_cmp_u64);
    size_t u = 0;
    for (size_t e = 0; e < num_edges; ++e) {
        if (u == 0 || edges[e] != edges[u - 1]) edges[u++] = edges[e];
    }
    num_edges = u;

    int32_t* deg = md_alloc(temp_arena, (size_t)nc * sizeof(int32_t));
    int32_t* off = md_alloc(temp_arena, ((size_t)nc + 1) * sizeof(int32_t));
    int32_t* adj = md_alloc(temp_arena, 2 * num_edges * sizeof(int32_t));
    int32_t* queue = md_alloc(temp_arena, (size_t)nc * sizeof(int32_t));
    uint8_t* alive = md_alloc(temp_arena, (size_t)nc);
    md_memset(deg, 0, (size_t)nc * sizeof(int32_t));
    for (size_t e = 0; e < num_edges; ++e) {
        deg[edges[e] >> 32] += 1;
        deg[edges[e] & 0xFFFFFFFF] += 1;
    }
    off[0] = 0;
    for (int32_t c = 0; c < nc; ++c) off[c + 1] = off[c] + deg[c];
    md_memset(deg, 0, (size_t)nc * sizeof(int32_t));
    for (size_t e = 0; e < num_edges; ++e) {
        const int32_t a = (int32_t)(edges[e] >> 32), b = (int32_t)(edges[e] & 0xFFFFFFFF);
        adj[off[a] + deg[a]++] = b;
        adj[off[b] + deg[b]++] = a;
    }

    // The core: the ends without tails taken off until none is left
    int32_t qe = 0;
    for (int32_t c = 0; c < nc; ++c) {
        alive[c] = 1;
        if (!has_tail[c] && deg[c] <= 1) {
            alive[c] = 0;
            queue[qe++] = c;
        }
    }
    for (int32_t qb = 0; qb < qe; ++qb) {
        const int32_t c = queue[qb];
        for (int32_t e = off[c]; e < off[c + 1]; ++e) {
            const int32_t d = adj[e];
            if (!alive[d]) continue;
            deg[d] -= 1;
            if (!has_tail[d] && deg[d] <= 1) {
                alive[d] = 0;
                queue[qe++] = d;
            }
        }
    }

    // And what is bonded to it
    for (int32_t c = 0; c < nc; ++c) {
        bool in = alive[c];
        for (int32_t e = off[c]; e < off[c + 1] && !in; ++e) in = alive[adj[e]];
        in_lipid[c] = (uint8_t)in;
    }
}

// The lipid of a molecule, if it is or holds one
static void infer_lipid_molecule(md_lipid_data_t* lipid, const md_system_t* sys, const lipid_mol_t* mol, const uint8_t* in_fused_ring, md_allocator_i* temp_arena, md_allocator_i* alloc) {
    const int32_t n = mol->n;
    const int32_t* atoms = mol->atoms;
    const int32_t nc = mol->num_residues;
    if (n <= 0) return;

    // Coarse grained when it holds a bead (a molecule too large to look into: when the first atom of a residue is one)
    bool cg = false;
    bool sterol_beads = false;
    if (atoms) {
        for (int32_t k = 0; k < n; ++k) {
            if (lipid_atom_is_bead(sys, atoms[k])) {
                cg = true;
                sterol_beads |= lipid_bead_sterol_ring(md_atom_name(&sys->atom, atoms[k]));
            }
        }
    } else {
        for (int32_t c = 0; c < nc && !cg; ++c) {
            const md_urange_t r = md_system_residue_atom_range(sys, (size_t)mol->residues[c]);
            cg = r.beg < r.end && lipid_atom_is_bead(sys, r.beg);
        }
    }
    // Lipids by name. A coarse grained lipid is one of the force field, by its name or, where that is not Martini's, by
    // its beads (the kind md_system_infer_residue_flags gave it), or one whose beads are named as those of the
    // lipids of the force field (below): never by its structure alone. A residue which is a lipid by its structure alone
    // (the kind an earlier inference gave it) is none by name: it is a lipid again only if it is found to be one
    // (md_system_infer_lipid_kinds), so that inferring again gives what inferring once did.
    bool named = false;
    for (int32_t c = 0; c < nc; ++c) {
        const size_t res = (size_t)mol->residues[c];
        if (md_resname_lipid(md_system_residue_name(sys, res)) || (cg && md_system_residue_kind(sys, res) == MD_RESIDUE_KIND_LIPID)) {
            named = true;
        }
    }
    // A lipid by name too small or too large to look into is one without heads and tails
    if (n <= 1 || n > LIPID_MAX_ATOMS || !atoms) {
        if (named) record_lipid(lipid, mol->residues, (size_t)nc, alloc);
        return;
    }

    const size_t pos = md_vm_arena_get_pos(temp_arena);

    uint8_t* apolar     = md_alloc(temp_arena, n);
    uint8_t* letterless = md_alloc(temp_arena, n);                      // apolar, if it makes a chain (LIPID_BEAD_LETTERLESS)
    int32_t* part       = md_alloc(temp_arena, n * sizeof(int32_t));   // the apolar part of an atom, -1
    int32_t* part_atoms = md_alloc(temp_arena, n * sizeof(int32_t));   // the atoms of each part, contiguous
    int32_t* dist       = md_alloc(temp_arena, n * sizeof(int32_t));
    int32_t* prev       = md_alloc(temp_arena, n * sizeof(int32_t));
    int32_t* queue      = md_alloc(temp_arena, n * sizeof(int32_t));
    int32_t* head       = md_alloc(temp_arena, n * sizeof(int32_t));   // the head of a heavy atom outside of the tails, -1
    uint8_t* tail       = md_alloc(temp_arena, n);                      // in a tail
    uint8_t* has_tail   = md_alloc(temp_arena, nc);                     // a residue holds an atom of a tail
    uint8_t* in_lipid   = md_alloc(temp_arena, nc);                     // a residue is of the lipid (lipid_extent)
    int32_t* lipid_residues = md_alloc(temp_arena, nc * sizeof(int32_t));
    md_array(lipid_part_t) parts = 0;

    for (int32_t k = 0; k < n; ++k) {
        const lipid_bead_t kind = lipid_atom_is_heavy(sys, atoms[k]) ? lipid_atom_kind(sys, atoms[k], sterol_beads) : LIPID_BEAD_POLAR;
        apolar[k] = (uint8_t)(kind != LIPID_BEAD_POLAR);
        letterless[k] = (uint8_t)(kind == LIPID_BEAD_LETTERLESS);
        part[k] = -1;
        head[k] = -1;
        tail[k] = 0;
    }
    md_memset(has_tail, 0, nc);

    // The connected apolar parts and what they hang off. The atoms of a part are part_atoms[beg, end), in the order
    // they were reached.
    int32_t num_part_atoms = 0;
    for (int32_t s = 0; s < n; ++s) {
        if (!apolar[s] || part[s] != -1) continue;
        const int32_t id = (int32_t)md_array_size(parts);
        lipid_part_t p = { .beg = num_part_atoms, .anchor = {-1, -1}, .letterless = true, .lettered = true };
        int32_t qb = num_part_atoms;
        part_atoms[num_part_atoms++] = s;
        part[s] = id;
        while (qb < num_part_atoms) {
            const int32_t k = part_atoms[qb++];
            int32_t degree = 0;
            p.letterless &= letterless[k] != 0;
            const bool bead = lipid_atom_is_bead(sys, atoms[k]);
            p.lettered &= bead && lipid_bead_lettered(md_atom_name(&sys->atom, atoms[k]));
            if (in_fused_ring && in_fused_ring[atoms[k]]) p.skeleton = true;
            // The ring beads of a sterol, whose bonds (a star of virtual sites in Martini 3) need not form rings
            if (bead && lipid_bead_sterol_ring(md_atom_name(&sys->atom, atoms[k]))) p.skeleton = true;
            md_bond_iter_t it = md_bond_iter(&sys->bond, atoms[k]);
            while (md_bond_iter_has_next(&it)) {
                const int32_t j = md_bond_iter_atom_index(&it);
                md_bond_iter_next(&it);
                if (!lipid_atom_is_heavy(sys, j)) continue;
                const int32_t l = lipid_mol_local(mol, j);
                if (l < 0) {
                    p.outside = true;
                    continue;
                }
                if (apolar[l]) {
                    p.bonds += 1;
                    degree += 1;
                    if (part[l] == -1) {
                        part[l] = id;
                        part_atoms[num_part_atoms++] = l;
                    }
                } else if (j != p.anchor[0] && j != p.anchor[1]) {
                    if (p.num_anchors < 2) p.anchor[p.num_anchors] = j;
                    p.num_anchors = MD_MIN(p.num_anchors + 1, 3);
                }
            }
            p.branched |= degree > 2;
        }
        p.end = num_part_atoms;
        md_array_push(parts, p, temp_arena);
    }

    // A coarse grained molecule the force field does not know, from a topology of its own (a ganglioside, a sphingomyelin
    // named otherwise): a lipid when a chain of beads named as those of the tails of the Martini lipidomes (C1A C2A C3A,
    // the letter of the tail after the number) hangs off the rest of the molecule, which has a bead of its own. Then it is
    // looked into as a lipid by name, so that inferring again (the kind lipid now its name) gives the same.
    if (cg && !named) {
        bool lettered_tail = false;
        const size_t num_parts = md_array_size(parts);
        for (size_t pi = 0; pi < num_parts && !lettered_tail; ++pi) {
            const lipid_part_t* p = &parts[pi];
            const int32_t size = p->end - p->beg;
            lettered_tail = p->lettered && !p->skeleton && !p->outside && !p->branched && p->num_anchors >= 1 && p->num_anchors <= 2 &&
                            size >= 2 && p->bonds / 2 == size - 1;
        }
        bool polar = false;
        for (int32_t k = 0; k < n && !polar; ++k) polar = !apolar[k] && lipid_atom_is_heavy(sys, atoms[k]);
        if (!lettered_tail || !polar) goto done;
        named = true;
    }

    // The tails among the parts, and their spines. The spine is found by a breadth first search within the part: from
    // every atom bonded to the anchor out to the farthest atom, or for a spanning chain to the nearest atom bonded to its
    // other anchor.
    typedef struct { int32_t beg, end; } spine_range_t;
    md_array(int32_t) spines = 0;
    md_array(spine_range_t) spine_ranges = 0;
    md_array(int32_t) tail_part = 0;
    md_array(md_lipid_tail_kind_t) tail_kind = 0;
    md_array(int32_t) tail_anchor_beg = 0;
    md_array(int32_t) tail_anchor_end = 0;

    const size_t num_parts = md_array_size(parts);
    for (size_t pi = 0; pi < num_parts; ++pi) {
        const lipid_part_t* p = &parts[pi];
        if (p->outside || p->num_anchors == 0 || p->num_anchors > 2) continue;
        const int32_t size = p->end - p->beg;
        const bool skeleton = p->skeleton;
        if (skeleton) {
            // The skeleton of a sterol hangs off its hydroxyl (or ester) alone
            if (p->num_anchors != 1 || size < (cg ? LIPID_MIN_SKELETON_BEADS : LIPID_MIN_SKELETON_ATOMS)) continue;
        } else if (p->letterless && (p->bonds / 2 != size - 1 || p->branched || size < LIPID_MIN_CHAIN_LETTERLESS)) {
            // Beads named without the letter of a tail are one only as a chain: the inositol of a phosphatidylinositol and
            // the sugars of a glycolipid, named so too, are rings
            continue;
        } else if (!named && p->bonds / 2 != size - 1) {
            // The chains of a lipid of no name hold no rings: a ring with alkyl groups (a drug, a dye) is no tail. Those of
            // lipids by name may (the cyclopropane of a fatty acid, the cyclopentanes of the chains of archaeal lipids).
            continue;
        }

        // Anchors in order of their index, so that a spanning chain runs from the head of the lower one
        int32_t anchor_beg = p->anchor[0];
        int32_t anchor_end = p->num_anchors == 2 ? p->anchor[1] : -1;
        if (anchor_end != -1 && anchor_end < anchor_beg) {
            const int32_t t = anchor_beg; anchor_beg = anchor_end; anchor_end = t;
        }

        int32_t qb = 0, qe = 0;
        for (int32_t m = p->beg; m < p->end; ++m) {
            const int32_t k = part_atoms[m];
            prev[k] = -1;
            dist[k] = -1;
            if (md_bond_find(&sys->bond, atoms[k], anchor_beg) >= 0) {
                dist[k] = 0;
                queue[qe++] = k;
            }
        }
        int32_t tip = -1;
        while (qb < qe) {
            const int32_t k = queue[qb++];
            if (anchor_end != -1) {
                if (md_bond_find(&sys->bond, atoms[k], anchor_end) >= 0) {
                    tip = k;
                    break;
                }
            } else if (tip == -1 || dist[k] > dist[tip] || (dist[k] == dist[tip] && k > tip)) {
                // The farthest atom, the last of them on a tie (the end of the side chain of a Martini cholesterol)
                tip = k;
            }
            md_bond_iter_t it = md_bond_iter(&sys->bond, atoms[k]);
            while (md_bond_iter_has_next(&it)) {
                const int32_t l = lipid_mol_local(mol, md_bond_iter_atom_index(&it));
                md_bond_iter_next(&it);
                if (l >= 0 && part[l] == (int32_t)pi && dist[l] == -1) {
                    dist[l] = dist[k] + 1;
                    prev[l] = k;
                    queue[qe++] = l;
                }
            }
        }
        if (tip == -1) continue;

        const int32_t len = dist[tip] + 1;
        const int32_t min_len = skeleton ? 1 : (cg ? (p->letterless ? LIPID_MIN_CHAIN_LETTERLESS : LIPID_MIN_CHAIN_BEADS) : (named ? LIPID_MIN_CHAIN_NAMED : LIPID_MIN_CHAIN_UNNAMED));
        if (len < min_len) continue;

        const int32_t base = (int32_t)md_array_size(spines);
        spine_range_t r = { base, base + len };
        md_array_resize(spines, (size_t)r.end, temp_arena);
        int32_t at = r.end - 1;
        for (int32_t k = tip; k != -1; k = prev[k]) {
            spines[at--] = atoms[k];
        }
        MD_ASSERT(at == base - 1);
        md_array_push(spine_ranges, r, temp_arena);
        md_array_push(tail_part, (int32_t)pi, temp_arena);
        md_array_push(tail_kind, skeleton ? MD_LIPID_TAIL_SKELETON : MD_LIPID_TAIL_CHAIN, temp_arena);
        md_array_push(tail_anchor_beg, anchor_beg, temp_arena);
        md_array_push(tail_anchor_end, anchor_end, temp_arena);
        for (int32_t m = p->beg; m < p->end; ++m) {
            tail[part_atoms[m]] = 1;
            has_tail[mol->member[part_atoms[m]]] = 1;
        }
    }

    // A lipid by name in which no tails are found (no bonds, a naming the beads do not follow) is one without heads and
    // tails: it is a lipid molecule still
    const size_t num_tails = md_array_size(spine_ranges);
    if (num_tails == 0) {
        if (named) record_lipid(lipid, mol->residues, (size_t)nc, alloc);
        goto done;
    }

    // The residues of the lipid, and of what is attached to it
    lipid_extent(in_lipid, sys, mol, has_tail, temp_arena);

    // A molecule of no lipid needs a polar head: an atom other than carbon outside of its tails. Its head is of the
    // elements of lipids (C, N, O, P, S) and is not aromatic: a porphyrin with a phytyl or farnesyl tail (chlorophyll,
    // pheophytin, heme A) and the like are cofactors to be shown as such.
    if (!named) {
        bool polar = false;
        for (int32_t k = 0; k < n; ++k) {
            if (tail[k] || !in_lipid[mol->member[k]]) continue;
            const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, atoms[k]);
            if (z <= MD_Z_H) continue;
            if (z != MD_Z_C && z != MD_Z_N && z != MD_Z_O && z != MD_Z_P && z != MD_Z_S) goto done;
            if (md_atom_flags(&sys->atom, atoms[k]) & MD_ATOM_FLAG_AROMATIC) goto done;
            polar |= z != MD_Z_C;
        }
        if (!polar) goto done;
    }

    // The heads: the connected parts of the heavy atoms of the lipid outside of the tails, each in ascending order
    const uint32_t head_base = (uint32_t)lipid->head.count;
    uint32_t num_heads = 0;
    for (int32_t s = 0; s < n; ++s) {
        if (tail[s] || head[s] != -1 || !in_lipid[mol->member[s]] || !lipid_atom_is_heavy(sys, atoms[s])) continue;
        int32_t qb = 0, qe = 0;
        queue[qe++] = s;
        head[s] = (int32_t)num_heads;
        while (qb < qe) {
            const int32_t k = queue[qb++];
            md_bond_iter_t it = md_bond_iter(&sys->bond, atoms[k]);
            while (md_bond_iter_has_next(&it)) {
                const int32_t j = md_bond_iter_atom_index(&it);
                md_bond_iter_next(&it);
                const int32_t l = lipid_mol_local(mol, j);
                if (l < 0) continue;
                if (!tail[l] && head[l] == -1 && in_lipid[mol->member[l]] && lipid_atom_is_heavy(sys, j)) {
                    head[l] = (int32_t)num_heads;
                    queue[qe++] = l;
                }
            }
        }
        for (int32_t m = 0; m < qe; ++m) queue[m] = atoms[queue[m]];
        qsort(queue, (size_t)qe, sizeof(int32_t), lipid_cmp_i32);
        push_index_range(&lipid->head.offset, &lipid->head.atom_idx, queue, (size_t)qe, alloc);
        num_heads += 1;
    }

    for (size_t t = 0; t < num_tails; ++t) {
        const spine_range_t r = spine_ranges[t];
        push_index_range(&lipid->tail.spine_offset, &lipid->tail.spine_idx, spines + r.beg, (size_t)(r.end - r.beg), alloc);

        // Its atoms and their hydrogens, ascending
        const lipid_part_t* p = &parts[tail_part[t]];
        int32_t m = 0;
        for (int32_t q = p->beg; q < p->end; ++q) {
            const int32_t i = atoms[part_atoms[q]];
            queue[m++] = i;
            md_bond_iter_t it = md_bond_iter(&sys->bond, i);
            while (md_bond_iter_has_next(&it)) {
                const int32_t j = md_bond_iter_atom_index(&it);
                md_bond_iter_next(&it);
                if (lipid_mol_local(mol, j) >= 0 && md_atom_atomic_number(&sys->atom, j) == MD_Z_H && m < n) {
                    queue[m++] = j;
                }
            }
        }
        qsort(queue, (size_t)m, sizeof(int32_t), lipid_cmp_i32);
        {
            // A hydrogen bonded to two atoms of the tail (an inferred bond too many) is in it once
            int32_t u = 0;
            for (int32_t q = 0; q < m; ++q) {
                if (u == 0 || queue[q] != queue[u - 1]) queue[u++] = queue[q];
            }
            m = u;
        }
        push_index_range(&lipid->tail.atom_offset, &lipid->tail.atom_idx, queue, (size_t)m, alloc);

        const int32_t ab = tail_anchor_beg[t];
        const int32_t ae = tail_anchor_end[t];
        const uint32_t hb = head_base + (uint32_t)head[lipid_mol_local(mol, ab)];
        const uint32_t he = ae != -1 ? head_base + (uint32_t)head[lipid_mol_local(mol, ae)] : UINT32_MAX;
        md_array_push(lipid->tail.anchor_beg, ab, alloc);
        md_array_push(lipid->tail.anchor_end, ae, alloc);
        md_array_push(lipid->tail.head_beg, hb, alloc);
        md_array_push(lipid->tail.head_end, he, alloc);
        md_array_push(lipid->tail.kind, tail_kind[t], alloc);
    }

    lipid->head.count += num_heads;
    lipid->tail.count += num_tails;

    {
        // What is attached to the lipid is no lipid, whatever its name
        int32_t num_lipid_residues = 0;
        for (int32_t c = 0; c < nc; ++c) {
            if (in_lipid[c]) {
                lipid_residues[num_lipid_residues++] = mol->residues[c];
            }
        }
        record_lipid(lipid, lipid_residues, (size_t)num_lipid_residues, alloc);
    }

done:
    md_vm_arena_set_pos_back(temp_arena, pos);
}

bool md_lipid_infer(md_lipid_data_t* out, const md_system_t* sys, md_allocator_i* alloc) {
    MD_ASSERT(out);
    MD_ASSERT(sys);
    MD_ASSERT(alloc);
    MD_ASSERT(out != &sys->lipid);
    md_lipid_data_free(out, alloc);
    md_lipid_data_t* lipid = out;
    init_lipid_data(lipid, alloc);
    // Without bonds every lipid by name is a molecule of its own residue, without heads and tails
    if (sys->atom.count == 0 || sys->residue.count == 0) {
        return true;
    }

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    // The atoms of fused rings: a skeleton (a sterol) rather than a chain
    uint8_t* in_fused_ring = NULL;
    const size_t num_rings = md_ring_data_count(&sys->ring);
    for (size_t r = 0; r < num_rings; ++r) {
        if (!(md_system_ring_flags(sys, r) & MD_RING_FLAG_FUSED)) continue;
        if (!in_fused_ring) {
            in_fused_ring = md_alloc(temp_arena, sys->atom.count);
            md_memset(in_fused_ring, 0, sys->atom.count);
        }
        const int32_t* atoms = md_ring_atoms(&sys->ring, r);
        const size_t n = md_ring_size(&sys->ring, r);
        for (size_t k = 0; k < n; ++k) in_fused_ring[atoms[k]] = 1;
    }

    // The molecules: the residues which may be of a lipid (not amino acids, nucleotides, water or ions), joined by the
    // bonds between them, wherever they are in the system
    const size_t num_residues = sys->residue.count;
    int32_t* root = md_alloc(temp_arena, num_residues * sizeof(int32_t));
    for (size_t c = 0; c < num_residues; ++c) {
        root[c] = lipid_candidate_res(sys, c) ? (int32_t)c : -1;
    }
    for (size_t c = 0; c < num_residues; ++c) {
        if (root[c] < 0) continue;
        const md_urange_t range = md_system_residue_atom_range(sys, c);
        for (uint32_t i = range.beg; i < range.end; ++i) {
            md_bond_iter_t it = md_bond_iter(&sys->bond, i);
            while (md_bond_iter_has_next(&it)) {
                const uint32_t j = (uint32_t)md_bond_iter_atom_index(&it);
                md_bond_iter_next(&it);
                if (range.beg <= j && j < range.end) continue;
                const md_residue_idx_t d = md_system_residue_find_by_atom_idx(sys, j);
                if (d < 0 || root[d] < 0) continue;
                lipid_union(root, (int32_t)c, (int32_t)d);
            }
        }
    }

    // The residues of each molecule, ascending, after one another
    int32_t* count = md_alloc(temp_arena, num_residues * sizeof(int32_t));
    int32_t* start = md_alloc(temp_arena, num_residues * sizeof(int32_t));
    md_memset(count, 0, num_residues * sizeof(int32_t));
    size_t num_members = 0;
    for (size_t c = 0; c < num_residues; ++c) {
        if (root[c] < 0) continue;
        root[c] = lipid_find(root, (int32_t)c);
        count[root[c]] += 1;
        num_members += 1;
    }
    int32_t* members = md_alloc(temp_arena, MD_MAX(num_members, 1) * sizeof(int32_t));
    int32_t pos = 0;
    for (size_t c = 0; c < num_residues; ++c) {
        start[c] = pos;
        if (root[c] == (int32_t)c) pos += count[c];
    }
    md_memset(count, 0, num_residues * sizeof(int32_t));
    for (size_t c = 0; c < num_residues; ++c) {
        if (root[c] < 0) continue;
        const int32_t r = root[c];
        members[start[r] + count[r]++] = (int32_t)c;
    }

    for (size_t c = 0; c < num_residues; ++c) {
        if (root[c] != (int32_t)c) continue;
        const int32_t* residues = members + start[c];
        const int32_t nc = count[c];

        size_t n = 0;
        for (int32_t m = 0; m < nc; ++m) {
            const md_urange_t r = md_system_residue_atom_range(sys, (size_t)residues[m]);
            n += r.end - r.beg;
        }
        lipid_mol_t mol = { .residues = residues, .num_residues = nc, .n = (int32_t)MD_MIN(n, (size_t)INT32_MAX) };
        const size_t mark = md_vm_arena_get_pos(temp_arena);
        if (n <= LIPID_MAX_ATOMS) {
            mol.atoms  = md_alloc(temp_arena, MD_MAX(n, 1) * sizeof(int32_t));
            mol.member = md_alloc(temp_arena, MD_MAX(n, 1) * sizeof(int32_t));
            int32_t k = 0;
            for (int32_t m = 0; m < nc; ++m) {
                const md_urange_t r = md_system_residue_atom_range(sys, (size_t)residues[m]);
                for (uint32_t i = r.beg; i < r.end; ++i) {
                    mol.atoms[k] = (int32_t)i;
                    mol.member[k] = m;
                    k += 1;
                }
            }
            mol.contiguous = n > 0 && mol.atoms[n - 1] - mol.atoms[0] == (int32_t)n - 1;
        }
        infer_lipid_molecule(lipid, sys, &mol, in_fused_ring, temp_arena, alloc);
        md_vm_arena_set_pos_back(temp_arena, mark);
    }

    md_temp_end(temp);
    return true;
}

void md_system_infer_lipid_kinds(md_system_t* sys) {
    MD_ASSERT(sys);
    const md_lipid_data_t* lipid = &sys->lipid;
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    uint8_t* in_lipid = md_temp_alloc_zero_array(temp, uint8_t, MD_MAX(sys->residue.count, 1));
    for (size_t l = 0; l < lipid->count; ++l) {
        const int32_t* res = md_lipid_residues(lipid, l);
        const size_t n = md_lipid_residue_count(lipid, l);
        for (size_t k = 0; k < n; ++k) {
            if (res[k] >= 0 && (size_t)res[k] < sys->residue.count) in_lipid[res[k]] = 1;
        }
    }
    for (size_t r = 0; r < sys->residue.count; ++r) {
        if (in_lipid[r]) {
            if (md_system_residue_kind(sys, r) != MD_RESIDUE_KIND_LIPID) md_system_set_residue_kind(sys, r, MD_RESIDUE_KIND_LIPID);
            continue;
        }
        // A lipid by its structure (an earlier inference) which is none now; one by name or by its beads stays
        if (md_system_residue_kind(sys, r) != MD_RESIDUE_KIND_LIPID) continue;
        if (md_resname_lipid(md_system_residue_name(sys, r))) continue;
        const md_urange_t range = md_system_residue_atom_range(sys, r);
        if (range.beg < range.end && lipid_atom_is_bead(sys, range.beg)) continue;
        md_system_set_residue_kind(sys, r, MD_RESIDUE_KIND_OTHER);
    }
    md_temp_end(temp);
}

size_t md_lipid_tail_atom_mask(uint8_t* out, const md_lipid_data_t* lipid, size_t atom_count) {
    MD_ASSERT(out);
    MD_ASSERT(lipid);
    md_memset(out, 0, atom_count);
    size_t count = 0;
    const size_t n = md_array_size(lipid->tail.atom_idx);
    for (size_t k = 0; k < n; ++k) {
        const md_atom_idx_t a = lipid->tail.atom_idx[k];
        if (a >= 0 && (size_t)a < atom_count && !out[a]) { out[a] = 1; count += 1; }
    }
    return count;
}
