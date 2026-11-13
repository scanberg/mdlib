#include "md_infer_internal.h"

// The Martini force fields (md_martini_db.inl): which version a system is from, the residues and beads of the database
// a residue matches, and the bonds of a coarse grained system.

// The residues of the Martini force fields: proteins (vermouth), nucleic acids, the lipidomes, sterols, solvents, ions,
// small molecules, sugars and polymers (polyply) of Martini 2 and 3, read from their topology files by
// tools/gen_martini_db.py, and the links between the residues of the polymers.
#include "md_martini_db.inl"
#include <core/md_executor.h>
#include <core/md_geom.h>

const martini_db_bead_t* mdi_martini_db_bead(const martini_db_residue_t* res, size_t k) {
    return &martini_db_beads[martini_db_residue_beads[res->bead_offset + k]];
}

// An amino acid or a nucleotide, of a chain whose residues vary (termini, modifications, missing beads)
static inline bool martini_db_biopolymer(const martini_db_residue_t* res) {
    return res->kind == MD_RESIDUE_KIND_AMINO_ACID || res->kind == MD_RESIDUE_KIND_NUCLEOTIDE;
}

MD_STATIC_ASSERT(MD_COUNTOF(martini_db_residues) <= UINT16_MAX, "martini_votes keeps the ranges of names as 16 bits");

// The residues of the name, [*beg, *end)
static void martini_db_find(md_str_t name, size_t* beg, size_t* end) {
    size_t lo = 0, hi = MD_COUNTOF(martini_db_residues);
    while (lo < hi) {
        const size_t mid = (lo + hi) / 2;
        if (md_str_cmp_lex(md_str_from_cstr(martini_db_residues[mid].name), name) < 0) lo = mid + 1; else hi = mid;
    }
    size_t e = lo;
    while (e < MD_COUNTOF(martini_db_residues) && md_str_eq_cstr(name, martini_db_residues[e].name)) ++e;
    *beg = lo;
    *end = e;
}

// Whether the labels are the beads of the residue, each once (in any order). map[i] is then the bead of label i.
static bool martini_db_match(const martini_db_residue_t* res, const md_str_t labels[], size_t num_labels, uint8_t map[]) {
    if (num_labels != res->num_beads) return false;
    bool used[256] = {0};
    for (size_t i = 0; i < num_labels; ++i) {
        size_t k = 0;
        while (k < res->num_beads && (used[k] || !md_str_eq_cstr(labels[i], mdi_martini_db_bead(res, k)->name))) ++k;
        if (k == res->num_beads) return false;
        used[k] = true;
        map[i] = (uint8_t)k;
    }
    return true;
}

// The molecule (not an amino acid or a nucleotide) a residue of the labels is, of the version if there is one of it,
// else of any version. NULL when it is none.
const martini_db_residue_t* mdi_martini_db_match_molecule(uint32_t version, md_str_t res_name, const md_str_t labels[], size_t num_labels, uint8_t map[]) {
    if (num_labels == 0 || num_labels > 255) return NULL;
    size_t beg, end;
    martini_db_find(res_name, &beg, &end);
    const martini_db_residue_t* other = NULL;
    for (size_t k = beg; k < end; ++k) {
        const martini_db_residue_t* res = &martini_db_residues[k];
        if (martini_db_biopolymer(res) || !martini_db_match(res, labels, num_labels, map)) continue;
        if (res->versions & version) return res;
        if (!other) other = res;
    }
    if (other) martini_db_match(other, labels, num_labels, map);
    return other;
}

// The molecule (not an amino acid or a nucleotide) whose beads the labels are, whatever its name: a residue renamed, as the
// three letters of a PDB file make CHOL CHL, PSM SPM and PCER CER. Of the version where there is one of it, and with the
// masses (optional) of its beads. NULL for a single bead, which says too little, and where molecules of different kinds
// have those beads. The names of a molecule of the Martini force fields are those of its chemistry, so molecules of the
// same beads are the same molecule as far as anything here goes (its kind, its bonds, its anatomy as a lipid).
const martini_db_residue_t* mdi_martini_db_match_beads(uint32_t version, const md_str_t labels[], const float masses[], size_t num_labels, uint8_t map[]) {
    if (num_labels < 2 || num_labels > 255) return NULL;
    const martini_db_residue_t* match = NULL;   // Of the version
    const martini_db_residue_t* other = NULL;   // Of any
    for (size_t r = 0; r < MD_COUNTOF(martini_db_residues); ++r) {
        const martini_db_residue_t* res = &martini_db_residues[r];
        if (res->num_beads != num_labels || martini_db_biopolymer(res) || !martini_db_match(res, labels, num_labels, map)) continue;
        bool same = true;
        for (size_t i = 0; i < num_labels && masses && same; ++i) same = masses[i] == mdi_martini_db_bead(res, map[i])->mass;
        if (!same) continue;
        if (other && other->kind != res->kind) return NULL;
        if (!other) other = res;
        if (!match && (res->versions & version)) match = res;
    }
    const martini_db_residue_t* res = match ? match : other;
    if (res) martini_db_match(res, labels, num_labels, map);
    return res;
}

// A molecule of the force field (not an amino acid or a nucleotide) of the name whose beads are those of none of its
// topologies, as most of the labels (more than half) are beads of the name: a variant of the molecule whose topology is
// not in the database, or a residue edited. NULL for a name of no such molecule, for a single bead (which says too
// little), and where most labels are not beads of it (atoms named as an atomistic force field names them, of which C1 and
// C2 are beads of a Martini cholesterol too). Of several molecules of the name, the first.
const martini_db_residue_t* mdi_martini_db_unknown_variant(md_str_t res_name, const md_str_t labels[], size_t num_labels) {
    if (num_labels < 2) return NULL;
    size_t beg, end;
    martini_db_find(res_name, &beg, &end);
    const martini_db_residue_t* first = NULL;
    for (size_t k = beg; k < end && !first; ++k) {
        if (!martini_db_biopolymer(&martini_db_residues[k])) first = &martini_db_residues[k];
    }
    if (!first) return NULL;
    size_t known = 0;
    for (size_t i = 0; i < num_labels; ++i) {
        bool found = false;
        for (size_t k = beg; k < end && !found; ++k) {
            const martini_db_residue_t* res = &martini_db_residues[k];
            if (martini_db_biopolymer(res)) continue;
            for (size_t b = 0; b < res->num_beads && !found; ++b) found = md_str_eq_cstr(labels[i], mdi_martini_db_bead(res, b)->name);
        }
        known += found;
    }
    return 2 * known > num_labels ? first : NULL;
}

// Whether the residue is a monomer of a polymer: a link names it
bool mdi_martini_db_linked(md_str_t name) {
    for (size_t k = 0; k < MD_COUNTOF(martini_db_links); ++k) {
        if (md_str_eq_cstr(name, martini_db_links[k].name) || md_str_eq_cstr(name, martini_db_links[k].next_name)) return true;
    }
    return false;
}

// Bead idx of an amino acid or a nucleotide of the labels: of the residue of the version whose beads they are, else the
// first bead of its name in a residue of the version (the residue may lack beads), else in one of any version
const martini_db_bead_t* mdi_martini_db_biopolymer_bead(uint32_t version, md_str_t res_name, const md_str_t labels[], size_t num_labels, size_t idx) {
    size_t beg, end;
    martini_db_find(res_name, &beg, &end);
    if (num_labels <= 255) {
        uint8_t map[256];
        for (size_t r = beg; r < end; ++r) {
            const martini_db_residue_t* res = &martini_db_residues[r];
            if (martini_db_biopolymer(res) && (res->versions & version) && martini_db_match(res, labels, num_labels, map)) return mdi_martini_db_bead(res, map[idx]);
        }
    }
    const md_str_t label = labels[idx];
    const martini_db_bead_t* other = NULL;
    for (size_t r = beg; r < end; ++r) {
        const martini_db_residue_t* res = &martini_db_residues[r];
        if (!martini_db_biopolymer(res)) continue;
        for (size_t k = 0; k < res->num_beads; ++k) {
            const martini_db_bead_t* bead = mdi_martini_db_bead(res, k);
            if (!md_str_eq_cstr(label, bead->name)) continue;
            if (res->versions & version) return bead;
            if (!other) other = bead;
        }
    }
    return other;
}

// Which version of Martini a system is from. Each residue whose beads are those of a residue of the database in some
// versions only votes for those: Ala with SC1, Tyr with SC4 and Trp with SC5 for Martini 3, the charge dummies and the
// polarizable water for polarizable Martini 2, cholesterol of 8 beads for Martini 2, ... Martini 3 when nothing tells,
// the current one and martinize2's default, and Martini 2 over polarizable Martini 2 on a tie.
// The votes and the word of the residues [beg, end[
typedef struct martini_votes_t {
    size_t votes[3];    // MARTINI_2, MARTINI_2P, MARTINI_3
    bool   martini;
} martini_votes_t;

static martini_votes_t martini_votes(const md_system_t* sys, md_str_t (*label)(const void*, size_t), const void* data, size_t res_beg, size_t res_end) {
    martini_votes_t v = {0};
    md_str_t labels[256];
    uint8_t map[256];
    // The residues of the database of a name, by the name's word: a system has few names (the amino acids, a few
    // lipids), whose binary searches through the database are done once
    struct { uint64_t key; uint16_t beg, end; } found[256];
    md_memset(found, 0, sizeof(found));
    for (size_t res_idx = res_beg; res_idx < res_end; ++res_idx) {
        const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
        const size_t n = range.end - range.beg;
        if (n == 0 || n > 255) continue;
        size_t beg, end;
        const uint64_t key = sys->residue.name ? md_label_key(sys->residue.name[res_idx]) : 0;
        const size_t slot = (size_t)((key * 0x9E3779B97F4A7C15ull) >> 56);
        if (key != 0 && found[slot].key == key) {
            beg = found[slot].beg;
            end = found[slot].end;
        } else {
            martini_db_find(md_system_residue_name(sys, res_idx), &beg, &end);
            if (key != 0) {
                found[slot].key = key;
                found[slot].beg = (uint16_t)beg;
                found[slot].end = (uint16_t)end;
            }
        }
        // Only residues of as many beads as it has atoms can be it: most atomistic residues are none
        bool fits = false;
        for (size_t k = beg; k < end && !fits; ++k) fits = martini_db_residues[k].num_beads == n;
        if (!fits) continue;
        for (size_t i = 0; i < n; ++i) labels[i] = label(data, range.beg + i);
        uint32_t versions = 0;
        for (size_t k = beg; k < end; ++k) {
            const martini_db_residue_t* res = &martini_db_residues[k];
            if (!martini_db_match(res, labels, n, map)) continue;
            versions |= res->versions;
            // No atomistic residue has the name and the atom names of a Martini residue of several beads, of Martini water
            // or of a monomer of a Martini polymer
            if (n > 1 || res->kind == MD_RESIDUE_KIND_WATER || mdi_martini_db_linked(md_str_from_cstr(res->name))) v.martini = true;
        }
        if (versions == 0 || versions == MARTINI_ANY) continue;
        if (versions & MARTINI_2)  v.votes[0] += 1;
        if (versions & MARTINI_2P) v.votes[1] += 1;
        if (versions & MARTINI_3)  v.votes[2] += 1;
    }
    return v;
}

#define MARTINI_VOTE_CHUNK 4096

typedef struct martini_vote_task_t {
    const md_system_t* sys;
    md_str_t (*label)(const void*, size_t);
    const void* data;
    martini_votes_t* chunk;
} martini_vote_task_t;

static void martini_vote_chunks(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    martini_vote_task_t* t = (martini_vote_task_t*)arg;
    for (uint32_t c = beg; c < end; ++c) {
        const size_t r_beg = (size_t)c * MARTINI_VOTE_CHUNK;
        t->chunk[c] = martini_votes(t->sys, t->label, t->data, r_beg, MD_MIN(r_beg + MARTINI_VOTE_CHUNK, t->sys->residue.count));
    }
}

martini_system_t mdi_martini_system_ex(const md_system_t* sys, md_str_t (*label)(const void*, size_t), const void* data, const md_executor_t* ex) {
    martini_votes_t v = {0};
    const size_t num_chunks = MD_DIV_UP(sys->residue.count, MARTINI_VOTE_CHUNK);
    if (md_executor_num_threads(ex) > 1 && num_chunks > 1) {
        // A count of each chunk, summed: exact, whichever thread had which
        md_temp_scope_t temp = md_temp_begin();
        martini_vote_task_t task = { sys, label, data, md_temp_alloc(temp, num_chunks * sizeof(martini_votes_t)) };
        md_memset(task.chunk, 0, num_chunks * sizeof(martini_votes_t));
        md_parallel_for(ex, (uint32_t)num_chunks, 1, martini_vote_chunks, &task);
        for (size_t c = 0; c < num_chunks; ++c) {
            for (int k = 0; k < 3; ++k) v.votes[k] += task.chunk[c].votes[k];
            v.martini |= task.chunk[c].martini;
        }
        md_temp_end(temp);
    } else {
        v = martini_votes(sys, label, data, 0, sys->residue.count);
    }
    martini_system_t result = { .version = MARTINI_3, .martini = v.martini };
    if (v.votes[0] > v.votes[2] && v.votes[0] >= v.votes[1]) result.version = MARTINI_2;
    else if (v.votes[1] > v.votes[2] && v.votes[1] > v.votes[0]) result.version = MARTINI_2P;
    return result;
}

martini_system_t mdi_martini_system(const md_system_t* sys, md_str_t (*label)(const void*, size_t), const void* data) {
    return mdi_martini_system_ex(sys, label, data, NULL);
}

// Lipids by name: those of the Martini lipidomes (whose names the atomistic force fields mostly share: POPC, DPPC, CHOL,
// ...) and the cholesterol of CHARMM
bool md_resname_lipid(md_str_t str) {
    size_t beg, end;
    martini_db_find(str, &beg, &end);
    for (size_t k = beg; k < end; ++k) {
        if (martini_db_residues[k].kind == MD_RESIDUE_KIND_LIPID) return true;
    }
    return md_str_eq_cstr(str, "CHL1");
}

// Coarse grained (Martini) bond inference. A residue whose beads are those of a residue of the Martini force fields
// (martini_db), with their masses, has the bonds of its topology: those within the residue for an amino acid or a
// nucleotide, all of them for a molecule (a lipid, a sterol, a solvent, a sugar, ...). Amino acids and nucleotides which
// are not (missing beads, modified residues, other force fields) are bonded by their geometry:
//
// Within a residue
// - A minimum spanning tree of its beads, of bonds up to 6 A: the beads hang together the way their topologies have
//   them, and a stray bead stays unbonded.
// - The rings: side chain beads closer than 3.6 A (base beads closer than 3.0 A) besides those the tree bonds. Longer
//   ring bonds (Trp SC2-SC4 of Martini 3, 4.7 A; purine SC2-SC4 of Martini 2 DNA) are left out, as beads of the rings
//   that are not bonded come as close as 3.3 A.
// - A massless bead is a virtual site, which the topologies do not bond (Trp SC3 of Martini 3), and stays out, unless
//   it is all of the side chain there is (the polar SC1 of polarizable Martini 2, bonded to BB in its topology).
// - The charge dummies of polarizable Martini 2 (SCP, SCN) are bonded to the bead they sit by, and to each other.
// Between residues (of any amino acid or nucleotide)
// - The last backbone bead of a residue to the first of the next (BB to BB; BB3 to BB1 in DNA and RNA, and BB2 to BB1
//   as well in Martini 3 RNA), up to 5.5 A when the next is numbered next (bonds stretch in a simulation), else up to
//   4.3 A for amino acids (bonded pairs of structures reach 4.05 A) and 5.5 A for nucleotides (4.4 A), unless the
//   numbering starts over: a new chain. The numbering may jump forward over a bond (an engineered insert), so a gap is
//   told by the distance alone.
// - The SC1 beads of two cysteines closer than 3 A: a disulfide bridge (2.3-2.5 A; unbonded pairs are 4.5 A apart).
//
// The residues of a polymer (polyply's library: PEO, PS, PE, PMMA, dextran, ..., a PEGylated lipid) are bonded to the
// residue before them by the links of the force field.
// A virtual site is bonded to the first bead it is constructed from in the topologies, as md_tpr bonds it.
// Validated against the topologies martinize2 (vermouth 0.15) writes for 16 structures (dipeptides to 300 residues,
// several chains, disulfides, a chain with gaps) in Martini 3, 2 and polarizable 2, against Martini 2 DNA (polyply)
// mapped from a nucleosome, and against simulated systems of proteins, DNA, RNA, membranes and solvents in Martini 2
// and 3 (their .tpr): no bond which the topologies do not have, besides elastic networks and the restraints of
// secondary structure, which are not chemical bonds.
static inline float cg_bead_distance(const md_system_state_t* state, int i, int j) {
    md_float3_t d = md_float3_sub(state->xyz[i], state->xyz[j]);
    md_geom_min_image(d.elem, 0, 1, &state->unitcell);
    return md_float3_length(d);
}

static inline void cg_push_bond(md_array(bond_pair_t)* candidates, int i, int j, float dist, md_allocator_i* alloc) {
    bond_pair_t cp = { .atom_i = i, .atom_j = j, .dist = dist };
    md_array_push(*candidates, cp, alloc);
}

static inline bool martini_charge_dummy(md_str_t atom_name) {
    return md_str_eq_cstr(atom_name, "SCP") || md_str_eq_cstr(atom_name, "SCN");
}

// The residue of the database a residue is: the one of its name whose beads are its beads, each with the mass of its
// atom (its type was given by it), and of the version where there are several. map[i] is then the bead of atom i.
static const martini_db_residue_t* martini_db_match_residue(const md_system_t* sys, size_t res_idx, uint32_t version, md_str_t labels[], uint8_t map[]) {
    const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
    const size_t n = range.end - range.beg;
    if (n == 0 || n > 255) return NULL;
    const md_residue_kind_t kind = md_system_residue_kind(sys, res_idx);
    const bool biopolymer = kind == MD_RESIDUE_KIND_AMINO_ACID || kind == MD_RESIDUE_KIND_NUCLEOTIDE;
    for (size_t i = 0; i < n; ++i) labels[i] = md_atom_name(&sys->atom, range.beg + i);
    size_t beg, end;
    martini_db_find(md_system_residue_name(sys, res_idx), &beg, &end);
    const martini_db_residue_t* other = NULL;
    for (size_t r = beg; r < end; ++r) {
        const martini_db_residue_t* res = &martini_db_residues[r];
        if (martini_db_biopolymer(res) != biopolymer || !martini_db_match(res, labels, n, map)) continue;
        bool same = true;
        for (size_t i = 0; i < n && same; ++i) same = md_atom_mass(&sys->atom, range.beg + i) == mdi_martini_db_bead(res, map[i])->mass;
        if (!same) continue;
        if (res->versions & version) return res;
        if (!other) other = res;
    }
    if (other) {
        martini_db_match(other, labels, n, map);
        return other;
    }
    // A molecule renamed: by its beads and their masses
    if (biopolymer) return NULL;
    float masses[256];
    for (size_t i = 0; i < n; ++i) masses[i] = md_atom_mass(&sys->atom, range.beg + i);
    return mdi_martini_db_match_beads(version, labels, masses, n, map);
}

// The bonds of a residue of a polymer (polyply) to the one before it, by the links of the force field, up to 6 A apart
// (the longest are 4.7 A), unless the numbering starts over: a new chain
static void cg_polymer_links(md_array(bond_pair_t)* candidates, const md_system_t* sys, const md_system_state_t* state, size_t ci, md_allocator_i* alloc) {
    const float polymer_cap = 6.0f;
    const md_str_t prev_name = md_system_residue_name(sys, ci - 1);
    const md_str_t name      = md_system_residue_name(sys, ci);
    const md_sequence_id_t prev_seq = md_system_residue_seq_id(sys, ci - 1);
    const md_sequence_id_t seq      = md_system_residue_seq_id(sys, ci);
    if (seq < prev_seq && !(seq <= 1 && prev_seq >= 9999)) return;
    const md_urange_t prev = md_system_residue_atom_range(sys, ci - 1);
    const md_urange_t cur  = md_system_residue_atom_range(sys, ci);
    for (size_t k = 0; k < MD_COUNTOF(martini_db_links); ++k) {
        const martini_db_link_t* link = &martini_db_links[k];
        if (!md_str_eq_cstr(prev_name, link->name) || !md_str_eq_cstr(name, link->next_name)) continue;
        int a = -1, b = -1;
        for (uint32_t i = prev.beg; i < prev.end && a == -1; ++i) if (md_str_eq_cstr(md_atom_name(&sys->atom, i), link->bead)) a = (int)i;
        for (uint32_t i = cur.beg;  i < cur.end  && b == -1; ++i) if (md_str_eq_cstr(md_atom_name(&sys->atom, i), link->next_bead)) b = (int)i;
        if (a == -1 || b == -1) continue;
        const float d = cg_bead_distance(state, a, b);
        if (d < polymer_cap) cg_push_bond(candidates, a, b, d, alloc);
    }
}

void mdi_infer_coarse_grained_bonds(md_array(bond_pair_t)* candidates, const md_system_t* sys, const md_system_state_t* state, md_allocator_i* alloc) {
    const float tree_cap      = 6.0f;
    const float dummy_cap     = 3.0f;   // A charge dummy sits 1.1-1.4 A from its bead
    const float dummy_pair    = 3.2f;   // SCP-SCN are constrained at 2.8 A
    const float disulfide_cap = 3.0f;

    size_t max_len = 0;
    for (size_t ci = 0; ci < sys->residue.count; ++ci) {
        const md_urange_t r = md_system_residue_atom_range(sys, ci);
        max_len = MD_MAX(max_len, (size_t)(r.end - r.beg));
    }
    int*   tree   = md_alloc(alloc, sizeof(int)   * (max_len + 1));
    int*   parent = md_alloc(alloc, sizeof(int)   * (max_len + 1));
    float* best   = md_alloc(alloc, sizeof(float) * (max_len + 1));
    bool*  done   = md_alloc(alloc, sizeof(bool)  * (max_len + 1));
    md_array(int) cys_sc1 = 0;

    int prev_last_bb = -1, prev_also_bb = -1;
    size_t prev_ci = 0;
    md_residue_kind_t prev_kind = MD_RESIDUE_KIND_OTHER;

    md_str_t* labels = md_alloc(alloc, sizeof(md_str_t) * (max_len + 1));
    uint8_t map[256];
    const uint32_t version = mdi_martini_system(sys, label_of_atom, &sys->atom).version;

    for (size_t ci = 0; ci < sys->residue.count; ++ci) {
        const md_residue_kind_t kind = md_system_residue_kind(sys, ci);
        const bool biopolymer = kind == MD_RESIDUE_KIND_AMINO_ACID || kind == MD_RESIDUE_KIND_NUCLEOTIDE;
        const md_urange_t range = md_system_residue_atom_range(sys, ci);

        if (ci > 0) cg_polymer_links(candidates, sys, state, ci, alloc);

        // The bonds of its topology. A charge site (a virtual site without Lennard-Jones) has none, as md_tpr gives it.
        const martini_db_residue_t* res = martini_db_match_residue(sys, ci, version, labels, map);
        if (res) {
            uint8_t local[256];
            for (size_t i = 0; i < range.end - range.beg; ++i) local[map[i]] = (uint8_t)i;
            for (size_t b = 0; b < res->num_bonds; ++b) {
                const int x = (int)range.beg + local[martini_db_residue_bonds[res->bond_offset + b][0]];
                const int y = (int)range.beg + local[martini_db_residue_bonds[res->bond_offset + b][1]];
                if (md_atom_particle_kind(&sys->atom, x) == MD_PARTICLE_VIRTUAL_SITE || md_atom_particle_kind(&sys->atom, y) == MD_PARTICLE_VIRTUAL_SITE) continue;
                cg_push_bond(candidates, x, y, cg_bead_distance(state, x, y), alloc);
            }
        }
        if (!biopolymer) {
            prev_last_bb = prev_also_bb = -1;
            continue;
        }
        const bool  amino_acid = kind == MD_RESIDUE_KIND_AMINO_ACID;
        const float ring_cap   = amino_acid ? 3.6f : 3.0f;
        const float link_cap      = amino_acid ? 4.3f : 5.5f;
        const float link_cap_next = 5.5f;

        // Without a topology, by the geometry. Which beads make up the tree: all but the charge dummies and, where the
        // side chain has beads with a mass, the massless ones
        bool side_with_mass = false;
        for (uint32_t i = range.beg; i < range.end; ++i) {
            if (!(md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_BACKBONE) && !martini_charge_dummy(md_atom_name(&sys->atom, i)) && md_atom_mass(&sys->atom, i) > 0.0f) {
                side_with_mass = true;
            }
        }
        // The next residue is linked to the last bead of the backbone, and, in a nucleotide of Martini 3, to BB2 as well
        const bool m3_nucleotide = !amino_acid && ((res ? res->versions : version) & (MARTINI_2 | MARTINI_2P)) == 0;
        int n = 0;
        int first_bb = -1, last_bb = -1, also_bb = -1;
        for (uint32_t i = range.beg; i < range.end; ++i) {
            const bool backbone = md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_BACKBONE;
            if (backbone) {
                if (first_bb == -1) first_bb = (int)i;
                last_bb = (int)i;
                if (m3_nucleotide && md_str_eq_cstr(md_atom_name(&sys->atom, i), "BB2")) also_bb = (int)i;
            }
            if (res || martini_charge_dummy(md_atom_name(&sys->atom, i))) continue;
            if (!backbone && side_with_mass && md_atom_mass(&sys->atom, i) == 0.0f) continue;
            tree[n++] = (int)i;
        }

        // Prim's minimum spanning tree
        for (int k = 0; k < n; ++k) {
            best[k] = FLT_MAX;
            parent[k] = -1;
            done[k] = false;
        }
        if (n > 0) best[0] = 0.0f;
        for (int iter = 0; iter < n; ++iter) {
            int k = -1;
            for (int m = 0; m < n; ++m) {
                if (!done[m] && (k == -1 || best[m] < best[k])) k = m;
            }
            done[k] = true;
            if (parent[k] >= 0 && best[k] <= tree_cap) {
                cg_push_bond(candidates, tree[parent[k]], tree[k], best[k], alloc);
            } else {
                parent[k] = -1;
            }
            for (int m = 0; m < n; ++m) {
                if (done[m]) continue;
                const float d = cg_bead_distance(state, tree[k], tree[m]);
                if (d < best[m]) {
                    best[m] = d;
                    parent[m] = k;
                }
            }
        }

        // Ring bonds between side chain beads, besides the tree's
        for (int a = 0; a < n; ++a) {
            if (md_atom_flags(&sys->atom, tree[a]) & MD_ATOM_FLAG_BACKBONE) continue;
            for (int b = a + 1; b < n; ++b) {
                if (md_atom_flags(&sys->atom, tree[b]) & MD_ATOM_FLAG_BACKBONE) continue;
                if (parent[a] == b || parent[b] == a) continue;
                const float d = cg_bead_distance(state, tree[a], tree[b]);
                if (d < ring_cap) cg_push_bond(candidates, tree[a], tree[b], d, alloc);
            }
        }

        // The charge dummies
        for (uint32_t i = range.beg; i < range.end && !res; ++i) {
            if (!martini_charge_dummy(md_atom_name(&sys->atom, i))) continue;
            int nearest = -1;
            float nearest_d = dummy_cap;
            for (uint32_t j = range.beg; j < range.end; ++j) {
                if (j == i) continue;
                const float d = cg_bead_distance(state, (int)i, (int)j);
                if (martini_charge_dummy(md_atom_name(&sys->atom, j))) {
                    if (j > i && d < dummy_pair) cg_push_bond(candidates, (int)i, (int)j, d, alloc);
                } else if (d < nearest_d) {
                    nearest = (int)j;
                    nearest_d = d;
                }
            }
            if (nearest != -1) cg_push_bond(candidates, nearest, (int)i, nearest_d, alloc);
        }

        // To the previous residue of the chain. The next number allows for the stretch of a bond in a simulation
        // (Martini 2's backbone bonds, the softest, pass 4.3 A in one of 200 frames and 5.5 A in one of 10^5).
        if (prev_last_bb != -1 && first_bb != -1 && prev_ci + 1 == ci && prev_kind == kind) {
            const md_sequence_id_t prev_seq = md_system_residue_seq_id(sys, prev_ci);
            const md_sequence_id_t seq      = md_system_residue_seq_id(sys, ci);
            const bool restart = seq < prev_seq && !(seq <= 1 && prev_seq >= 9999);  // GRO and PDB wrap at 99999 and 9999
            const bool next    = seq == prev_seq + 1 || (seq <= 1 && prev_seq >= 9999);
            if (!restart) {
                const float cap = next ? link_cap_next : link_cap;
                const float d = cg_bead_distance(state, prev_last_bb, first_bb);
                if (d < cap) {
                    cg_push_bond(candidates, prev_last_bb, first_bb, d, alloc);
                    if (prev_also_bb != -1 && prev_also_bb != prev_last_bb) {
                        const float d2 = cg_bead_distance(state, prev_also_bb, first_bb);
                        if (d2 < cap) cg_push_bond(candidates, prev_also_bb, first_bb, d2, alloc);
                    }
                }
            }
        }
        prev_last_bb = last_bb;
        prev_also_bb = also_bb;
        prev_ci = ci;
        prev_kind = kind;

        const md_str_t res_name = md_system_residue_name(sys, ci);
        if (amino_acid && (md_str_eq_cstr(res_name, "CYS") || md_str_eq_cstr(res_name, "CYX"))) {
            for (uint32_t i = range.beg; i < range.end; ++i) {
                if (md_str_eq_cstr(md_atom_name(&sys->atom, i), "SC1")) md_array_push(cys_sc1, (int)i, alloc);
            }
        }
    }

    // Disulfide bridges
    const size_t num_cys = md_array_size(cys_sc1);
    for (size_t a = 0; a < num_cys; ++a) {
        for (size_t b = a + 1; b < num_cys; ++b) {
            const float d = cg_bead_distance(state, cys_sc1[a], cys_sc1[b]);
            if (d < disulfide_cap) cg_push_bond(candidates, cys_sc1[a], cys_sc1[b], d, alloc);
        }
    }
}
