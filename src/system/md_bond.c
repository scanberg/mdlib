#include "md_infer_internal.h"
#include <core/md_executor.h>
#include <core/md_fifo.h>
#include <core/md_log.h>
#include <core/md_spatial_acc.h>
#include <core/md_geom.h>

// Bond inference: covalent and metallic bonds from the coordinates, coordination and elastic bonds, and the bonds
// formed and broken in a state (md_util.h).

// ### BOND DATA ###

const char* md_bond_origin_name(md_bond_origin_t origin) {
    switch (origin) {
    case MD_BOND_ORIGIN_FILE:       return "file";
    case MD_BOND_ORIGIN_TOPOLOGY:   return "topology";
    case MD_BOND_ORIGIN_USER:       return "user";
    case MD_BOND_ORIGIN_INFERRED:   return "inferred";
    default:                        return "";
    }
}

static const md_bond_kind_name_t bond_kind_names[] = {
    {"covalent",    MD_BOND_FLAG_KIND_MASK,     (md_bond_flags_t)MD_BOND_KIND_COVALENT},
    {"coordinate",  MD_BOND_FLAG_KIND_MASK,     (md_bond_flags_t)MD_BOND_KIND_COORDINATE},
    {"elastic",     MD_BOND_FLAG_KIND_MASK,     (md_bond_flags_t)MD_BOND_KIND_ELASTIC},
    {"metallic",    MD_BOND_FLAG_KIND_MASK,     (md_bond_flags_t)MD_BOND_KIND_METALLIC},
    {"noncovalent", (md_bond_flags_t)0xC0,      (md_bond_flags_t)0xC0},     // md_bond_is_noncovalent
    {"interaction", MD_BOND_FLAG_KIND_MASK,     (md_bond_flags_t)MD_BOND_KIND_INTERACTION},
    {"hydrogen",    MD_BOND_FLAG_KIND_MASK,     (md_bond_flags_t)MD_BOND_KIND_HYDROGEN},
};

const md_bond_kind_name_t* md_bond_kind_names(size_t* count) {
    if (count) *count = MD_COUNTOF(bond_kind_names);
    return bond_kind_names;
}

// The connectivity is the graph of the molecules: every bond of the first range but the elastic ones
// (MD_BOND_KIND_ELASTIC), springs of the force field which connect nothing. bond_flags NULL: every bond is covalent.
static inline bool conn_has_bond(const md_bond_flags_t* bond_flags, size_t i) {
    return !bond_flags || md_bond_is_graph(bond_flags[i]);
}

static void build_connectivity(md_bond_conn_data_t* conn, const md_atom_pair_t* bond_pairs, const md_bond_flags_t* bond_flags, size_t bond_pair_count, size_t atom_count, md_allocator_i* alloc) {
    MD_ASSERT(conn);
    MD_ASSERT(alloc);

    if (bond_pairs == NULL) return;
    if (bond_pair_count == 0) return;
    if (atom_count == 0) return;

    conn->offset_count = atom_count + 1;
    md_array_resize(conn->offset, conn->offset_count, alloc);
    md_memset(conn->offset, 0, md_array_bytes(conn->offset));

    size_t graph_count = 0;
    for (size_t i = 0; i < bond_pair_count; ++i) graph_count += conn_has_bond(bond_flags, i);

    // This have length of 2 * bond_count (one for each direction of the bond)
    conn->count = 2 * graph_count;
    md_array_resize(conn->atom_idx, conn->count, alloc);
    md_array_resize(conn->bond_idx, conn->count, alloc);

    typedef struct {
        uint16_t off[2];
    } offset_t;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);

    offset_t* local_offset = md_temp_alloc_zero_array(temp, offset_t, bond_pair_count);
    MD_ASSERT(local_offset);

    // Two packed 16-bit local offsets for each of the bond idx
    // Use offsets as accumulators for length
    for (size_t i = 0; i < bond_pair_count; ++i) {
        if (!conn_has_bond(bond_flags, i)) continue;
		MD_ASSERT(bond_pairs[i].idx[0] < (md_atom_idx_t)atom_count);
		MD_ASSERT(bond_pairs[i].idx[1] < (md_atom_idx_t)atom_count);
        local_offset[i].off[0] = (uint16_t)conn->offset[bond_pairs[i].idx[0]]++;
        local_offset[i].off[1] = (uint16_t)conn->offset[bond_pairs[i].idx[1]]++;
    }

    // Compute complete edge offsets (exclusive scan)
    uint32_t off = 0;
    for (size_t i = 0; i < conn->offset_count; ++i) {
        const uint32_t len = conn->offset[i];
        conn->offset[i] = off;
        off += len;
    }

    // Write edge indices to correct location
    for (size_t i = 0; i < bond_pair_count; ++i) {
        if (!conn_has_bond(bond_flags, i)) continue;
        const md_atom_pair_t p = bond_pairs[i];
        const int atom_a = p.idx[0];
        const int atom_b = p.idx[1];
        const int local_a = (int)local_offset[i].off[0];
        const int local_b = (int)local_offset[i].off[1];
        const int off_a = conn->offset[atom_a];
        const int off_b = conn->offset[atom_b];

        const int idx_a = off_a + local_a;
        const int idx_b = off_b + local_b;

        MD_ASSERT(idx_a < (int)conn->count);
        MD_ASSERT(idx_b < (int)conn->count);

        // Store the cross references to the 'other' atom index signified by the bond in the correct location
        conn->atom_idx[idx_a] = atom_b;
        conn->atom_idx[idx_b] = atom_a;

        conn->bond_idx[idx_a] = (md_bond_idx_t)i;
        conn->bond_idx[idx_b] = (md_bond_idx_t)i;
    }

    // Each atom's neighbours ascending by atom index, whatever order the bonds were found in (a spatial search finds
    // them in cell order). What walks the graph - the slot order of the molecules (md_molecule_data_t) above all -
    // then depends on the atoms' order alone, so that molecules alike, written alike, are walked alike.
    for (size_t a = 0; a < atom_count; ++a) {
        const uint32_t beg = conn->offset[a], end = conn->offset[a + 1];
        for (uint32_t k = beg + 1; k < end; ++k) {
            const md_atom_idx_t ai = conn->atom_idx[k];
            const md_bond_idx_t bi = conn->bond_idx[k];
            uint32_t j = k;
            for (; j > beg && conn->atom_idx[j - 1] > ai; --j) {
                conn->atom_idx[j] = conn->atom_idx[j - 1];
                conn->bond_idx[j] = conn->bond_idx[j - 1];
            }
            conn->atom_idx[j] = ai;
            conn->bond_idx[j] = bi;
        }
    }

    md_temp_end(temp);
}

void md_bond_build_connectivity(md_bond_data_t* in_out_bond, size_t atom_count, md_allocator_i* alloc) {
    MD_ASSERT(in_out_bond);
    MD_ASSERT(alloc);
	build_connectivity(&in_out_bond->conn, in_out_bond->pairs, in_out_bond->flags, in_out_bond->count, atom_count, alloc);
}

void md_bond_data_free(md_bond_data_t* bond, md_allocator_i* alloc) {
    MD_ASSERT(bond);
    MD_ASSERT(alloc);
    md_array_free(bond->pairs, alloc);
    md_array_free(bond->flags, alloc);
    md_array_free(bond->conn.atom_idx, alloc);
    md_array_free(bond->conn.bond_idx, alloc);
    md_array_free(bond->conn.offset, alloc);
    md_memset(bond, 0, sizeof(md_bond_data_t));
}

size_t md_bond_metallic_atoms(uint8_t* out, const md_bond_data_t* bond, size_t atom_count) {
    MD_ASSERT(out);
    MD_ASSERT(bond);
    md_memset(out, 0, atom_count);
    size_t count = 0;
    for (size_t i = 0; i < bond->count; ++i) {
        if (md_bond_kind(bond->flags[i]) != MD_BOND_KIND_METALLIC) continue;
        for (int k = 0; k < 2; ++k) {
            const md_atom_idx_t a = bond->pairs[i].idx[k];
            if (a >= 0 && (size_t)a < atom_count && !out[a]) { out[a] = 1; count += 1; }
        }
    }
    return count;
}

// This is bitmasks expressed as combinations of two uint64_t that encodes properties of atomic elements
static const uint64_t element_metal_mask[2] = {
    0xff87ffe0fff83818,
    0x00001fffff9fffff,
};

static const uint64_t element_halogen_mask[2] = {
    0x0020000800020200,
    0x0000000000200000,
};

static const uint64_t element_alkali_mask[2] = {
    0x0080002000080808,
    0x0000000000800000,
};

// Groups 3-12: Sc-Zn, Y-Cd, Hf-Hg, Rf-Cn
static const uint64_t element_transition_metal_mask[2] = {
    0x0001FF807FE00000ULL,
    0x0001FF000001FF00ULL,
};

#if 0
static const uint64_t element_alkaline_earth[2] = {
    0x0100004000101010,
    0x0000000001000000,
};

static const uint64_t element_covalent_mask[2] = {
    0x003E0000001FC7E6ULL, // atomic numbers 1�63
    0x0000000000000000ULL, // atomic numbers 64�127 (unused)
};
#endif

// This is a macro to generate test functions of the bitmasks above
// This seems to generate decent assembly across Clang, GCC and MSVC
// If this is a static inline function, it fails to represent the raw values of the mask in assembly
#define GENERATE_MASK_FUNC(name, mask) \
static inline bool name(uint32_t elem) { \
    uint64_t off = 0; \
    uint64_t pat = mask[0]; \
    if (elem >= 64) { \
        off = 64; \
        pat = mask[1]; \
    } \
    uint64_t msk = (elem < 119) ? 0xFFFFFFFFULL : 0; \
    uint64_t idx = (elem - off) & msk; \
    return pat & (1ULL << idx); \
}

GENERATE_MASK_FUNC(is_metal, element_metal_mask)

GENERATE_MASK_FUNC(is_alkali, element_alkali_mask)

GENERATE_MASK_FUNC(is_transition_metal, element_transition_metal_mask)

GENERATE_MASK_FUNC(is_halogen, element_halogen_mask)

static inline bool is_metal_donor(size_t atomic_nr) {
    switch(atomic_nr) {
    case O:
    case N:
    case S:
    case P:
    case F:
    case Cl:
    case Br:
    case I:
        return true;
    default:
        return false;
    }
}

// The covalent radius of md_element.h, in the bond rules
static inline float element_covalent_radius(md_atomic_number_t element) {
    return md_atomic_number_covalent_radius(element);
}

typedef struct {
	const float* radius; // Bonding radius for each atom
    float k_min;
    float k_max;
    // The metals, which reach further than the other atoms and are searched for on their own (cov_bond_candidates):
    // NULL for none. The pass over all pairs leaves out the pairs with a metal, the pass from the metals ('i') keeps
    // a pair of metals once and no metal with itself.
    const uint8_t* metal;
    bool from_metals;
    md_array(bond_pair_t)* candidates;
    md_allocator_i* alloc;
} cov_bond_callback_param_t;

static void test_cov_bond_pair_callback(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    cov_bond_callback_param_t* data = (cov_bond_callback_param_t*)user_param;

    const md_256  v0   = md_mm256_setzero_ps();
    const md_256i v8   = md_mm256_set1_epi32(8);
    const md_256i vlen = md_mm256_set1_epi32((int)num_pairs);

    const md_256 vk_min = md_mm256_set1_ps(data->k_min);
    const md_256 vk_max = md_mm256_set1_ps(data->k_max);

    md_256i vi = md_mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
    for (size_t i = 0; i < num_pairs; i += 8) {
        md_256 vmask = md_mm256_castsi256_ps(md_mm256_cmplt_epi32(vi, vlen));

        md_256i vi_idx = md_mm256_loadu_epi32(i_idx + i);
        md_256i vj_idx = md_mm256_loadu_epi32(j_idx + i);
        md_256 v_dist2 = md_mm256_loadu_ps(ij_dist2 + i);
        md_256 v_dist  = md_mm256_sqrt_ps(v_dist2);  // There are alot of gathers here, so we can afford latency of sqrt

        float dist[8];
        md_mm256_storeu_ps(dist, v_dist);

        md_256 rad_i   = md_mm256_mask_i32gather_ps(v0, data->radius, vi_idx, vmask, 4);
        md_256 rad_j   = md_mm256_mask_i32gather_ps(v0, data->radius, vj_idx, vmask, 4);
        md_256 rad_sum = md_mm256_add_ps(rad_i, rad_j);

        md_256 r_min_cov = md_mm256_mul_ps(vk_min, rad_sum);
        md_256 r_max_cov = md_mm256_mul_ps(vk_max, rad_sum);

        md_256 cmp_min = md_mm256_cmpgt_ps(v_dist, r_min_cov);
        md_256 cmp_max = md_mm256_cmplt_ps(v_dist, r_max_cov);

        int mask = md_mm256_movemask_ps(md_mm256_and_ps(cmp_min, cmp_max));
        if (mask) {
            size_t new_size = md_array_size(*data->candidates) + md_popcnt32(mask);
            md_array_ensure(*data->candidates, new_size, data->alloc);

            while (mask) {
                int bit_idx = md_ctz32(mask);
                mask = mask & ~(1 << bit_idx);
                int idx = (int)i + bit_idx;
                if (data->metal) {
                    const uint32_t a = i_idx[idx], b = j_idx[idx];
                    if (data->from_metals ? (a == b || (data->metal[b] && b < a)) : (data->metal[a] || data->metal[b])) continue;
                }
                bond_pair_t cp = {
                    .atom_i = i_idx[idx],
                    .atom_j = j_idx[idx],
                    .dist   = dist[bit_idx],
                };
                md_array_push_no_grow(*data->candidates, cp);
            }
        }
        vi = md_mm256_add_epi32(vi, v8);
    }
}

// The candidate pairs of the whole structure, on the threads of an executor: each thread collects the pairs of the
// pieces it is handed into an array of its own, and the pieces are put together in their order, which is the order the
// serial query gives them in. So the candidates, and the bonds made of them, do not depend on the threads.
typedef struct cov_bond_segment_t {
    size_t   piece_beg;
    uint32_t thread;
    size_t   offset;
    size_t   count;
} cov_bond_segment_t;

typedef struct cov_bond_pass_t {
    const md_spatial_acc_t*     acc;
    double                      cutoff;
    // The pass from the metals: their coordinates and atoms, METAL_PIECE of them a piece. NULL for the pass over all pairs.
    const md_float3_t*          xyz;
    const int*                  metals;
    size_t                      num_metals;
    cov_bond_callback_param_t*  param;      // [num_threads], each with its own candidates
    md_array(cov_bond_segment_t)* segments; // [num_threads]
    md_allocator_i*             heap;
} cov_bond_pass_t;

#define METAL_PIECE 64

static void cov_bond_pass_range(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    cov_bond_pass_t* pass = (cov_bond_pass_t*)arg;
    cov_bond_callback_param_t* param = &pass->param[thread_idx];
    const size_t offset = md_array_size(*param->candidates);
    if (pass->metals) {
        const size_t m_beg = (size_t)beg * METAL_PIECE;
        const size_t m_end = MD_MIN((size_t)end * METAL_PIECE, pass->num_metals);
        md_coord_stream_t ext = md_coord_stream_from_aos((const float*)pass->xyz, sizeof(md_float3_t), pass->metals + m_beg, m_end - m_beg);
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(pass->acc, &ext, pass->cutoff, test_cov_bond_pair_callback, param, MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX);
    } else {
        md_spatial_acc_for_each_internal_pair_within_cutoff_in_pieces(pass->acc, pass->cutoff, beg, end, test_cov_bond_pair_callback, param);
    }
    const cov_bond_segment_t seg = { beg, thread_idx, offset, md_array_size(*param->candidates) - offset };
    if (seg.count) md_array_push(pass->segments[thread_idx], seg, pass->heap);
}

static int compare_cov_bond_segment(const void* a, const void* b) {
    const size_t x = ((const cov_bond_segment_t*)a)->piece_beg;
    const size_t y = ((const cov_bond_segment_t*)b)->piece_beg;
    return (x > y) - (x < y);
}

// The pairs within the bonding range of their radii: of all pairs in the order of
// md_spatial_acc_for_each_internal_pair_within_cutoff, or (metals given) from the metals in their order, each metal's
// pairs in the order of md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff
static void cov_bond_candidates(md_array(bond_pair_t)* out, const md_spatial_acc_t* acc, double cutoff, const md_float3_t* xyz, const int* metals, size_t num_metals,
                                const cov_bond_callback_param_t* proto, const md_executor_t* ex, md_allocator_i* alloc) {
    md_allocator_i* heap = md_get_heap_allocator();
    const uint32_t num_threads = md_executor_num_threads(ex);
    const size_t num_pieces = metals ? (num_metals + METAL_PIECE - 1) / METAL_PIECE : md_spatial_acc_num_pieces(acc);
    if (num_pieces == 0) return;

    md_array(bond_pair_t)*        cand = md_alloc(heap, num_threads * sizeof(md_array(bond_pair_t)));
    md_array(cov_bond_segment_t)* segs = md_alloc(heap, num_threads * sizeof(md_array(cov_bond_segment_t)));
    cov_bond_callback_param_t*    param = md_alloc(heap, num_threads * sizeof(cov_bond_callback_param_t));
    for (uint32_t t = 0; t < num_threads; ++t) {
        cand[t] = 0;
        segs[t] = 0;
        param[t] = *proto;
        param[t].candidates = &cand[t];
        param[t].alloc = heap;
        param[t].from_metals = metals != NULL;
    }
    cov_bond_pass_t pass = { .acc = acc, .cutoff = cutoff, .xyz = xyz, .metals = metals, .num_metals = num_metals, .param = param, .segments = segs, .heap = heap };
    const uint32_t grain = metals ? 1 : (uint32_t)MD_MAX(1, num_pieces / ((size_t)num_threads * 32));
    md_parallel_for(ex, (uint32_t)num_pieces, grain, cov_bond_pass_range, &pass);

    md_array(cov_bond_segment_t) all = 0;
    size_t total = 0;
    for (uint32_t t = 0; t < num_threads; ++t) {
        md_array_push_array(all, segs[t], md_array_size(segs[t]), heap);
        total += md_array_size(cand[t]);
    }
    // Not with a NULL array: qsort's pointer is declared non-null, which lets the compiler drop the NULL checks after it
    if (md_array_size(all) > 1) qsort(all, md_array_size(all), sizeof(cov_bond_segment_t), compare_cov_bond_segment);
    md_array_ensure(*out, md_array_size(*out) + total, alloc);
    for (size_t k = 0; k < md_array_size(all); ++k) {
        md_array_push_array(*out, cand[all[k].thread] + all[k].offset, all[k].count, alloc);
    }

    md_array_free(all, heap);
    for (uint32_t t = 0; t < num_threads; ++t) {
        md_array_free(cand[t], heap);
        md_array_free(segs[t], heap);
    }
    md_free(heap, cand,  num_threads * sizeof(md_array(bond_pair_t)));
    md_free(heap, segs,  num_threads * sizeof(md_array(cov_bond_segment_t)));
    md_free(heap, param, num_threads * sizeof(cov_bond_callback_param_t));
}

static int compare_atom_pair(const void* a, const void* b) {
    const md_atom_pair_t* pa = (const md_atom_pair_t*)a;
    const md_atom_pair_t* pb = (const md_atom_pair_t*)b;
    if (pa->idx[0] != pb->idx[0]) return pa->idx[0] < pb->idx[0] ? -1 : 1;
    if (pa->idx[1] != pb->idx[1]) return pa->idx[1] < pb->idx[1] ? -1 : 1;
    return 0;
}

// The most covalent bonds an atom of the element takes, 0 for no limit. Generous: the limit is there to catch
// atoms that clash (a hydrogen within bonding distance of two atoms), not to judge chemistry.
static inline int max_covalent_bonds(md_atomic_number_t z) {
    switch (z) {
    case H: case F:                     return 1;
    case O:                             return 3;   // Oxonium, within one molecule (see below)
    case B: case C: case N: case Si:    return 4;
    case Cl: case Br: case I:           return 4;   // Perchlorate and its kind
    case P: case S: case Se: case As:   return 6;
    default:                            return 0;
    }
}

// The most coordination bonds a metal takes
static inline int max_coordination(md_atomic_number_t z) {
    if (z == Be || z == Mg || z == Al || z == Zn || z == Ga) return 6;
    if ((z >= 57 && z <= 71) || (z >= 89 && z <= 103)) return 9;   // Lanthanides, actinides
    return 8;   // Seven coordinate Ru, eight coordinate Ca, ...
}

// A bond in excess of what its atoms can take. Removed in order: bonds between residues apart (bond_separation)
// first, then bonds between residues next to each other, then bonds within one, and within each the bonds whose
// length is furthest from the sum of the radii (a clash as much as a stretch).
typedef struct bond_excess_t {
    uint32_t inter;
    float    deviation;
    uint32_t idx;
} bond_excess_t;

// How far apart the residues of a bond's atoms are, the first key of the order of bond_excess_t: 0 within one
// residue, 1 between residues next to each other (along a chain, and the head and tails of a lipid that AMBER splits
// into residues), 2 between residues apart or of an atom of none. Of an atom whose bonds are too many, the bond to a
// residue apart is the clash of molecules that overlap (a packed starting structure), and goes before the bond that
// joins the residues of a molecule: the ester of a tail's C12 to its head, 1.54 A, is further from the radii than a
// clash of that C12 with an ester oxygen of a lipid that overlaps it, 1.42 A.
static inline uint32_t bond_separation(int res_i, int res_j) {
    if (res_i == res_j) return 0;
    if (res_i < 0 || res_j < 0) return 2;
    return (res_i - res_j == 1 || res_j - res_i == 1) ? 1 : 2;
}

static int compare_bond_excess(const void* a, const void* b) {
    const bond_excess_t* x = (const bond_excess_t*)a;
    const bond_excess_t* y = (const bond_excess_t*)b;
    if (x->inter != y->inter) return x->inter > y->inter ? -1 : 1;
    if (x->deviation != y->deviation) return x->deviation > y->deviation ? -1 : 1;
    return (x->idx > y->idx) - (x->idx < y->idx);
}

// Factors on the sum of the covalent radii (Cordero et al. 2008, which covers metal-ligand distances too)
static const float bond_k_min    = 0.65f;   // Shorter than this is a clash, not a bond
static const float bond_k_cov    = 1.15f;   // Covalent bond between non-metals
static const float bond_k_coord  = 1.30f;   // Upper end of the candidate search
static const float bond_k_ligand = 1.25f;   // Metal to a donor (N, O, S, P, halogen): 1.0-1.1 in well resolved structures
static const float bond_k_carbon = 1.15f;   // Transition metal to carbon (carbenes, carbonyls, cyclopentadienyl)
static const float bond_k_metal  = 1.00f;   // Metal-metal

// How two atoms are bonded from their distance: within factor times the sum of their covalent radii (0: never, whatever
// the distance), and whether the bond is a coordination bond. same_res: the atoms are of one residue, where a bond
// may be a little longer (large flexible ligands). water: one of them is of a water molecule.
typedef struct bond_rule_t {
    float factor;
    bool  coordinate;
} bond_rule_t;

static bond_rule_t bond_rule(int ei, int ej, bool same_res, bool water) {
    const bond_rule_t none = {0};
    const bool mi = is_metal(ei);
    const bool mj = is_metal(ej);
    if (!mi && !mj) {
        const bond_rule_t r = { same_res ? bond_k_cov * 1.05f : bond_k_cov, false };
        return r;
    }
    if (mi != mj) {
        const int metal_e = mi ? ei : ej;
        const int non_metal_e = mi ? ej : ei;
        // Alkali ions are free ions in solution: a coordination frozen from the reference frame would be
        // wrong in every other. Water around metals is left out for the same reason.
        if (is_alkali(metal_e) || water) return none;
        if (is_metal_donor(non_metal_e)) {
            const bond_rule_t r = { bond_k_ligand, true };
            return r;
        }
        if (non_metal_e == C && is_transition_metal(metal_e)) {
            const bond_rule_t r = { bond_k_carbon, true };
            return r;
        }
        return none;
    }
    if (is_transition_metal(ei) && is_transition_metal(ej)) {
        const bond_rule_t r = { bond_k_metal, false };
        return r;
    }
    return none;
}

// A halogen bonded to a non-metal (an organohalogen, is_halogen) is no ligand of a metal: its contact with one, a
// halogen lying on a metal surface, is no coordination bond. A halide is.

// ### METALLIC PHASES ###
//
// A metal in bulk, at a surface or in a particle is not bonded as a molecule is: each atom has as many neighbours as its
// lattice gives it, and a cutoff on the sum of covalent radii falls within their shell (2.72 Å for gold, against a
// first shell of 2.6-3.1 Å), which bonds a part of it that looks random and changes from frame to frame. The atoms of
// a metallic phase (md_bond_metallic_atoms) are instead bonded to every metal atom of their first coordination shell,
// by METALLIC bonds, which are not bonds of the graph (md_bond.h), and to nothing else: what touches a surface moves
// across it, so its contacts are those of a frame (md_bond_delta_compute) rather than of the topology.
//
// THE FIRST SHELL ends at the first minimum of the distances between the metal atoms, each over the sum of the
// covalent radii of its atoms: the gap between two shells of a crystal, or the first minimum of the radial distribution
// of a liquid or a strained lattice. It is the lowest point of their histogram beyond its first peak and up to 1.5
// times it, which holds the gap after the first shell in the lattices of the metals: their second shell is at 1.41 times
// the first in a close packed lattice (fcc, hcp), and at 1.15 in a body centred cubic one, whose two shells (8 and 6
// neighbours) merge into one when warm. Told from the system rather than given per element, as the covalent radii put
// the first shell of the metals at 0.94 (Fe) to 1.18 (Al) times their sum, and strain and temperature move and widen it.
// METALLIC_SHELL without enough pairs to tell.
//
// A METALLIC ATOM is a metal with at least METALLIC_MIN_NEIGHBORS metal atoms in its first shell, and at least three
// times as many as the non-metals within reach of a bond to it (a donor, N, O, S, P or a halogen, within bond_k_ligand;
// carbon within bond_k_carbon; whether the pair would be bonded or not, as an alkali to its anion is not): a cation of
// an ionic crystal or an oxide is not (NaCl and MgO: 12 metal neighbours to 6 anions, CaF2: 12 to 8, TiO2: 10 to 6),
// nor a metal of a molecule or a cluster (Fe4S4, FeMoco, Mn4CaO5: 3 or 4 metal neighbours at most), whereas an atom at
// the surface of a metal has 7 to 9 metal neighbours and an adsorbate or two. Then so is a metal atom in the first
// shell of a metallic atom with at least METALLIC_MIN_EDGE_NEIGHBORS metal neighbours, again at least three times its
// non-metals: the edges and corners of a slab or a particle, an adatom in a hollow site.
#define METALLIC_MIN_NEIGHBORS      6
#define METALLIC_MIN_EDGE_NEIGHBORS 3
#define METALLIC_SHELL              1.25f       // On the sum of the covalent radii, without enough pairs to tell
#define METALLIC_SHELL_MIN          1.10f
#define METALLIC_SHELL_MAX          1.60f
#define METALLIC_MIN_PAIRS          64          // Pairs within the shell for its end to be told from them
#define METALLIC_HIST_LO            0.80f       // The histogram of the pairs, on the sum of the radii
#define METALLIC_HIST_HI            1.90f
#define METALLIC_HIST_BIN           0.02f
#define METALLIC_HIST_BINS          55          // (HI - LO) / BIN

typedef struct metal_pair_search_t {
    const float* radius;
    float        q_max;                         // On the sum of the radii
    uint32_t*    hist;                          // [METALLIC_HIST_BINS], or NULL
    md_array(bond_pair_t)* pairs;               // Within q_max, or NULL
    md_allocator_i* alloc;
} metal_pair_search_t;

static void metal_pair_callback(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    metal_pair_search_t* s = (metal_pair_search_t*)user_param;
    for (size_t k = 0; k < num_pairs; ++k) {
        const uint32_t i = i_idx[k];
        const uint32_t j = j_idx[k];
        const float sum_r = s->radius[i] + s->radius[j];
        if (i == j || sum_r <= 0.0f) continue;
        const float d = sqrtf(ij_dist2[k]);
        const float q = d / sum_r;
        if (s->hist && q >= METALLIC_HIST_LO && q < METALLIC_HIST_HI) {
            s->hist[(int)((q - METALLIC_HIST_LO) / METALLIC_HIST_BIN)] += 1;
        }
        if (s->pairs && q > bond_k_min && q <= s->q_max) {
            const bond_pair_t p = { (int)MD_MIN(i, j), (int)MD_MAX(i, j), d };
            md_array_push(*s->pairs, p, s->alloc);
        }
    }
}

// The end of the first shell, on the sum of the radii, from the histogram of the metal-metal pairs (THE FIRST SHELL)
static float metallic_shell(const uint32_t hist[METALLIC_HIST_BINS]) {
    float smooth[METALLIC_HIST_BINS];
    for (int b = 0; b < METALLIC_HIST_BINS; ++b) {
        const uint32_t l = b > 0 ? hist[b - 1] : 0;
        const uint32_t r = b + 1 < METALLIC_HIST_BINS ? hist[b + 1] : 0;
        smooth[b] = (float)(l + 2 * hist[b] + r) * 0.25f;
    }
    const float q_peak_max = 1.35f;
    int peak = 0;
    for (int b = 1; b < METALLIC_HIST_BINS && METALLIC_HIST_LO + (b + 0.5f) * METALLIC_HIST_BIN < q_peak_max; ++b) {
        if (smooth[b] > smooth[peak]) peak = b;
    }
    uint32_t in_shell = 0;
    for (int b = 0; b <= peak; ++b) in_shell += hist[b];
    if (smooth[peak] <= 0.0f || in_shell < METALLIC_MIN_PAIRS) return METALLIC_SHELL;

    // The lowest point beyond the peak up to 1.5 times it, the middle of the run of bins as low as it
    const float q_peak = METALLIC_HIST_LO + (peak + 0.5f) * METALLIC_HIST_BIN;
    const int last = MD_MIN(METALLIC_HIST_BINS - 1, (int)((1.5f * q_peak - METALLIC_HIST_LO) / METALLIC_HIST_BIN));
    int lo = -1;
    for (int b = peak + 1; b <= last; ++b) {
        if (lo < 0 || smooth[b] < smooth[lo]) lo = b;
    }
    if (lo < 0) return METALLIC_SHELL;
    int hi = lo;
    while (hi + 1 <= last && smooth[hi + 1] == smooth[lo]) hi += 1;
    const float q = METALLIC_HIST_LO + (0.5f * (lo + hi) + 0.5f) * METALLIC_HIST_BIN;
    return MD_CLAMP(q, METALLIC_SHELL_MIN, METALLIC_SHELL_MAX);
}

// The metallic atoms (out_metallic, one per atom) and the pairs of them within the first shell (out_pairs), see
// METALLIC PHASES. contacts: the number of non-metals within reach of each metal atom.
static void infer_metallic_phases(uint8_t* out_metallic, md_array(bond_pair_t)* out_pairs, const md_system_state_t* state, const md_system_t* sys,
                                  const md_atomic_number_t* atomic_nr, const float* radius, const uint32_t* contacts, md_allocator_i* temp_arena) {
    const size_t num_atoms = state->num_atoms;
    md_memset(out_metallic, 0, num_atoms);

    md_array(int) metals = 0;
    float max_rad = 0.0f;
    for (size_t i = 0; i < num_atoms; ++i) {
        if (!is_metal(atomic_nr[i]) || md_atom_particle_kind(&sys->atom, i) == MD_PARTICLE_VIRTUAL_SITE) continue;
        md_array_push(metals, (int)i, temp_arena);
        max_rad = MD_MAX(max_rad, radius[i]);
    }
    const size_t num_metals = md_array_size(metals);
    if (num_metals < METALLIC_MIN_NEIGHBORS + 1 || max_rad <= 0.0f) return;

    // The histogram of the pairs, then the pairs within the first shell
    uint32_t hist[METALLIC_HIST_BINS] = {0};
    {
        const double cutoff = 2.0 * max_rad * METALLIC_HIST_HI;
        md_spatial_acc_t acc = { .alloc = temp_arena };
        md_spatial_acc_init(&acc, &(md_spatial_acc_desc_t){ .xyz = (const float*)state->xyz, .stride = sizeof(md_float3_t), .idx = metals, .count = num_metals, .cutoff = cutoff, .unitcell = &state->unitcell, .flags = MD_SPATIAL_ACC_FLAG_USE_SOURCE_IDX });
        metal_pair_search_t search = { .radius = radius, .hist = hist };
        md_spatial_acc_for_each_internal_pair_within_cutoff(&acc, cutoff, metal_pair_callback, &search);
    }
    const float shell = metallic_shell(hist);

    md_array(bond_pair_t) pairs = 0;
    {
        const double cutoff = 2.0 * max_rad * shell;
        md_spatial_acc_t acc = { .alloc = temp_arena };
        md_spatial_acc_init(&acc, &(md_spatial_acc_desc_t){ .xyz = (const float*)state->xyz, .stride = sizeof(md_float3_t), .idx = metals, .count = num_metals, .cutoff = cutoff, .unitcell = &state->unitcell, .flags = MD_SPATIAL_ACC_FLAG_USE_SOURCE_IDX });
        metal_pair_search_t search = { .radius = radius, .q_max = shell, .pairs = &pairs, .alloc = temp_arena };
        md_spatial_acc_for_each_internal_pair_within_cutoff(&acc, cutoff, metal_pair_callback, &search);
    }
    const size_t num_pairs = md_array_size(pairs);

    uint32_t* neighbors = md_alloc(temp_arena, sizeof(uint32_t) * num_atoms);
    md_memset(neighbors, 0, sizeof(uint32_t) * num_atoms);
    for (size_t k = 0; k < num_pairs; ++k) {
        neighbors[pairs[k].atom_i] += 1;
        neighbors[pairs[k].atom_j] += 1;
    }
    for (size_t m = 0; m < num_metals; ++m) {
        const int i = metals[m];
        out_metallic[i] = neighbors[i] >= METALLIC_MIN_NEIGHBORS && neighbors[i] >= 3 * contacts[i];
    }
    // The edges, corners and adatoms: in the shell of a metallic atom of the core
    uint8_t* edge = md_alloc(temp_arena, num_atoms);
    md_memset(edge, 0, num_atoms);
    for (size_t k = 0; k < num_pairs; ++k) {
        const int i = pairs[k].atom_i;
        const int j = pairs[k].atom_j;
        if (out_metallic[i] && !out_metallic[j]) edge[j] = 1;
        if (out_metallic[j] && !out_metallic[i]) edge[i] = 1;
    }
    for (size_t m = 0; m < num_metals; ++m) {
        const int i = metals[m];
        if (edge[i] && neighbors[i] >= METALLIC_MIN_EDGE_NEIGHBORS && neighbors[i] >= 3 * contacts[i]) out_metallic[i] = 1;
    }

    for (size_t k = 0; k < num_pairs; ++k) {
        if (out_metallic[pairs[k].atom_i] && out_metallic[pairs[k].atom_j]) md_array_push(*out_pairs, pairs[k], temp_arena);
    }
}

// The bonds inferred from the coordinates, into bond (empty), every one of origin MD_BOND_ORIGIN_INFERRED
static void infer_bonds(md_bond_data_t* bond, const md_system_state_t* state, const md_system_t* sys, md_allocator_i* alloc, const md_executor_t* executor) {
    md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp_scope);
    const md_bond_flags_t inferred = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_INFERRED);

    const bool is_coarse_grained = md_system_is_coarse_grained(sys);

    size_t num_atoms = state->num_atoms;
    md_array(bond_pair_t) candidates = 0;

    if (is_coarse_grained) {
        md_array_ensure(candidates, num_atoms, temp_arena);
        mdi_infer_coarse_grained_bonds(&candidates, sys, state, temp_arena);

        size_t num_candidates = md_array_size(candidates);
		md_array_ensure(bond->pairs, num_candidates, alloc);
        md_array_ensure(bond->flags, num_candidates, alloc);

        for (size_t i = 0; i < num_candidates; ++i) {
            md_atom_pair_t pair = {
                candidates[i].atom_i,
                candidates[i].atom_j,
            };
            md_array_push_no_grow(bond->pairs, pair);
            md_array_push_no_grow(bond->flags, inferred);
            bond->count += 1;
        }
	} else {
        // The largest radius of the metals and of the other atoms: a metal reaches further than the rest (2.0 Å for
        // potassium against 1.1 Å for sulfur), and a few ions are not to set how far every atom of the system is
        // searched from (cov_bond_candidates)
        double max_metal_rad = 0.0, max_other_rad = 0.0;
        size_t num_metals = 0;

        float* atom_radius = md_temp_alloc_array(temp_scope, float, num_atoms);
        md_atomic_number_t* atomic_nr = md_temp_alloc_array(temp_scope, md_atomic_number_t, num_atoms);
        uint8_t* metal = md_temp_alloc_array(temp_scope, uint8_t, num_atoms);
        for (size_t i = 0; i < num_atoms; ++i) {
            atomic_nr[i]  = md_atom_atomic_number(&sys->atom, i);
            // The covalent radius for metals too: the ionic radius of the metal plus the covalent radius of the
            // ligand mixes two scales and comes out 0.4-0.6 Å short of real metal-ligand distances.
            atom_radius[i] = element_covalent_radius(atomic_nr[i]);
            metal[i] = is_metal(atomic_nr[i]);
            num_metals += metal[i];
            if (metal[i]) max_metal_rad = MD_MAX(max_metal_rad, atom_radius[i]);
            else          max_other_rad = MD_MAX(max_other_rad, atom_radius[i]);
        }

        int* comp_id = md_temp_alloc_array(temp_scope, int, num_atoms);
        md_memset(comp_id, -1, sizeof(int) * num_atoms);
        size_t num_res = md_system_residue_count(sys);

        for (size_t i = 0; i < num_res; ++i) {
            md_urange_t range = md_system_residue_atom_range(sys, i);
            for (uint32_t j = range.beg; j < range.end; ++j) {
                comp_id[j] = (int)i;
            }
        }

        // The residues are usually not classified yet (that needs the bonds), so whether a residue is a water
        // molecule is told here from its atoms. Bond perception needs it around metals.
        bool* res_water = md_temp_alloc_array(temp_scope, bool, num_res);
        for (size_t i = 0; i < num_res; ++i) {
            res_water[i] = md_system_residue_kind(sys, i) == MD_RESIDUE_KIND_WATER || atoms_are_water(sys, md_system_residue_atom_range(sys, i));
        }

        // The candidates live on the heap: the spatial structure is given back (to the scratch it came from) before
        // anything else is made, and what is made after it reuses that memory rather than adding to it
        md_allocator_i* heap = md_get_heap_allocator();
        {
            // Find connections for all atoms
            cov_bond_callback_param_t param = {
                .radius = atom_radius,
                .k_min = bond_k_min,
                .k_max = bond_k_coord,
                .metal = num_metals ? metal : NULL,
            };

            // The pairs of other atoms within reach of their largest covalent radius: about 2.7 Å for organic matter.
            // Then the pairs of the metals, from each metal within reach of its radius and the largest of all, which
            // is at most twice as far as the structure is built for.
            const double cutoff       = MD_MAX(1.0, 2.0 * max_other_rad * bond_k_coord);
            const double metal_cutoff = (max_metal_rad + MD_MAX(max_metal_rad, max_other_rad)) * bond_k_coord;
            const double acc_cutoff   = num_metals ? MD_MAX(cutoff, 0.51 * metal_cutoff) : cutoff;

            // Build candidate list
            const md_temp_mark_t before_acc = md_temp_mark(temp_scope);
            md_spatial_acc_t acc = {.alloc = temp_arena};
            md_spatial_acc_init(&acc, &(md_spatial_acc_desc_t){ .xyz = (const float*)state->xyz, .stride = sizeof(md_float3_t), .count = state->num_atoms, .cutoff = cutoff, .unitcell = &state->unitcell });
            cov_bond_candidates(&candidates, &acc, cutoff, NULL, NULL, 0, &param, executor, heap);
            if (num_metals) {
                int* metals = md_temp_alloc_array(temp_scope, int, num_metals);
                size_t k = 0;
                for (size_t i = 0; i < num_atoms; ++i) if (metal[i]) metals[k++] = (int)i;
                cov_bond_candidates(&candidates, &acc, metal_cutoff, state->xyz, metals, num_metals, &param, executor, heap);
            }
            md_temp_reset(temp_scope, before_acc);
        }

        size_t num_candidates = md_array_size(candidates);

        //MD_LOG_DEBUG("Found %zu candidates", num_candidates);

        // The rule of each candidate, 0 for a pair which is not bonded: virtual sites sit within bonding distance of
        // their parent atom (the lone pairs of TIP5P at 0.7 Å) but are not bonded to anything
        bond_rule_t* rules = md_temp_alloc_zero_array(temp_scope, bond_rule_t, MD_MAX(num_candidates, 1));
        uint8_t* organohalogen = md_temp_alloc_zero_array(temp_scope, uint8_t, num_atoms);
        for (size_t i = 0; i < num_candidates; ++i) {
            const int ai = candidates[i].atom_i;
            const int aj = candidates[i].atom_j;
            if (md_atom_particle_kind(&sys->atom, ai) == MD_PARTICLE_VIRTUAL_SITE || md_atom_particle_kind(&sys->atom, aj) == MD_PARTICLE_VIRTUAL_SITE) continue;
            const int ci = comp_id[ai];
            const int cj = comp_id[aj];
            const bool is_water = (ci >= 0 && res_water[ci]) || (cj >= 0 && res_water[cj]);
            const bond_rule_t rule = bond_rule(atomic_nr[ai], atomic_nr[aj], ci == cj, is_water);
            if (rule.factor <= 0.0f || candidates[i].dist >= rule.factor * (atom_radius[ai] + atom_radius[aj])) continue;
            rules[i] = rule;
            if (!rule.coordinate && !is_metal(atomic_nr[ai]) && !is_metal(atomic_nr[aj])) {
                if (is_halogen(atomic_nr[ai])) organohalogen[ai] = 1;
                if (is_halogen(atomic_nr[aj])) organohalogen[aj] = 1;
            }
        }
        // An organohalogen coordinates nothing
        for (size_t i = 0; i < num_candidates; ++i) {
            if (rules[i].coordinate && (organohalogen[candidates[i].atom_i] || organohalogen[candidates[i].atom_j])) {
                rules[i].factor = 0.0f;
            }
        }

        // The metallic phases, from the metal atoms and the non-metals within reach of each (METALLIC PHASES)
        uint32_t* contacts = md_temp_alloc_zero_array(temp_scope, uint32_t, num_atoms);
        for (size_t i = 0; i < num_candidates; ++i) {
            const int ai = candidates[i].atom_i;
            const int aj = candidates[i].atom_j;
            if (is_metal(atomic_nr[ai]) == is_metal(atomic_nr[aj])) continue;
            const int metal = is_metal(atomic_nr[ai]) ? ai : aj;
            const int other = metal == ai ? aj : ai;
            if (md_atom_particle_kind(&sys->atom, other) == MD_PARTICLE_VIRTUAL_SITE) continue;
            const int z = atomic_nr[other];
            const float factor = z == C ? bond_k_carbon : is_metal_donor(z) ? bond_k_ligand : 0.0f;
            if (candidates[i].dist < factor * (atom_radius[ai] + atom_radius[aj])) contacts[metal] += 1;
        }
        uint8_t* metallic = md_temp_alloc_array(temp_scope, uint8_t, num_atoms);
        md_array(bond_pair_t) metallic_pairs = 0;
        infer_metallic_phases(metallic, &metallic_pairs, state, sys, atomic_nr, atom_radius, contacts, temp_arena);

        md_bond_data_t temp_bond = { 0 };
        // We are lazy and propagate the correct periodic distance for each bond, so we can do more informed pruning later (e.g. for coordination bonds)
        md_array(float) temp_bond_dist = 0;
        {
            // Room for all of them at once: arrays which grow in an arena leave their copies behind
            size_t num_bonds = md_array_size(metallic_pairs);
            for (size_t i = 0; i < num_candidates; ++i) num_bonds += rules[i].factor > 0.0f;
            md_array_ensure(temp_bond.pairs, num_bonds, temp_arena);
            md_array_ensure(temp_bond.flags, num_bonds, temp_arena);
            md_array_ensure(temp_bond_dist,  num_bonds, temp_arena);
        }

        for (size_t i = 0; i < num_candidates; ++i) {
            if (rules[i].factor <= 0.0f) continue;
            const int ai = candidates[i].atom_i;
            const int aj = candidates[i].atom_j;
            // An atom of a metallic phase has its metallic bonds, below, and no other: what touches a metal moves
            // across it, and has its contacts in a frame (md_bond_delta_compute)
            if (metallic[ai] || metallic[aj]) continue;
            md_bond_insert(&temp_bond, ai, aj, rules[i].coordinate ? md_bond_flags_set_kind(inferred, MD_BOND_KIND_COORDINATE) : inferred, temp_arena);
            md_array_push(temp_bond_dist, candidates[i].dist, temp_arena);
        }
        md_array_free(candidates, heap);
        {
            const md_bond_flags_t flags = md_bond_flags_set_kind(inferred, MD_BOND_KIND_METALLIC);
            for (size_t k = 0; k < md_array_size(metallic_pairs); ++k) {
                md_bond_insert(&temp_bond, metallic_pairs[k].atom_i, metallic_pairs[k].atom_j, flags, temp_arena);
                md_array_push(temp_bond_dist, metallic_pairs[k].dist, temp_arena);
            }
        }

        // Create bitfield to keep track of valid bonds during pruning (we do not remove bonds to ensure index stability)
        md_array(uint64_t) temp_bond_valid = make_bitfield(temp_bond.count, temp_arena);
        bitfield_set_all(temp_bond_valid, temp_bond.count); // Set all bits to 1 (valid)

        // ### PRUNE ###
        // By the counts of bonds of each atom: no connectivity is needed (the triangle pruning below would want it)

        // Atoms with more bonds than they can take: atoms that clash (overlapping molecules in a packed starting
        // structure, unresolved alternate conformations), most often a hydrogen within reach of two atoms. Bonds go
        // in the order of compare_bond_excess until every atom is within its limit. Covalent bonds and coordination
        // bonds are limited separately.
        {
            const size_t num_temp = temp_bond.count;
            int* num_cov   = md_temp_alloc_zero_array(temp_scope, int, num_atoms);
            int* num_coord = md_temp_alloc_zero_array(temp_scope, int, num_atoms);
            int* max_cov   = md_temp_alloc_zero_array(temp_scope, int, num_atoms);
            for (size_t b = 0; b < num_temp; ++b) {
                // A metallic bond counts towards nothing: its atoms have as many as their lattice gives them
                if (md_bond_kind(temp_bond.flags[b]) == MD_BOND_KIND_METALLIC) continue;
                const int x = temp_bond.pairs[b].idx[0];
                const int y = temp_bond.pairs[b].idx[1];
                int* n = md_bond_kind(temp_bond.flags[b]) == MD_BOND_KIND_COORDINATE ? num_coord : num_cov;
                n[x] += 1;
                n[y] += 1;
                // A third bond on an oxygen only within its own molecule (oxonium); one reaching into another
                // molecule is a clash
                if (md_bond_kind(temp_bond.flags[b]) != MD_BOND_KIND_COORDINATE && comp_id[x] != comp_id[y]) {
                    if (atomic_nr[x] == O) max_cov[x] = 2;
                    if (atomic_nr[y] == O) max_cov[y] = 2;
                }
            }
            for (size_t i = 0; i < num_atoms; ++i) {
                if (!max_cov[i]) max_cov[i] = max_covalent_bonds(atomic_nr[i]);
            }

            md_array(bond_excess_t) over = 0;
            for (size_t b = 0; b < num_temp; ++b) {
                if (md_bond_kind(temp_bond.flags[b]) == MD_BOND_KIND_METALLIC) continue;
                const int x = temp_bond.pairs[b].idx[0];
                const int y = temp_bond.pairs[b].idx[1];
                bool exceeds = false;
                if (md_bond_kind(temp_bond.flags[b]) == MD_BOND_KIND_COORDINATE) {
                    const int metal = is_metal(atomic_nr[x]) ? x : y;
                    exceeds = num_coord[metal] > max_coordination(atomic_nr[metal]);
                } else {
                    const int mx = max_cov[x];
                    const int my = max_cov[y];
                    exceeds = (mx && num_cov[x] > mx) || (my && num_cov[y] > my);
                }
                if (exceeds) {
                    const float sum_r = atom_radius[x] + atom_radius[y];
                    bond_excess_t e = {
                        .inter = bond_separation(comp_id[x], comp_id[y]),
                        .deviation = sum_r > 0 ? fabsf(temp_bond_dist[b] / sum_r - 1.0f) : 0.0f,
                        .idx = (uint32_t)b,
                    };
                    md_array_push(over, e, temp_arena);
                }
            }

            const size_t num_over = md_array_size(over);
            if (num_over) {
                qsort(over, num_over, sizeof(bond_excess_t), compare_bond_excess);
                for (size_t k = 0; k < num_over; ++k) {
                    const uint32_t b = over[k].idx;
                    const int x = temp_bond.pairs[b].idx[0];
                    const int y = temp_bond.pairs[b].idx[1];
                    if (md_bond_kind(temp_bond.flags[b]) == MD_BOND_KIND_COORDINATE) {
                        const int metal = is_metal(atomic_nr[x]) ? x : y;
                        if (num_coord[metal] <= max_coordination(atomic_nr[metal])) continue;
                        num_coord[x] -= 1;
                        num_coord[y] -= 1;
                    } else {
                        const int mx = max_cov[x];
                        const int my = max_cov[y];
                        if (!((mx && num_cov[x] > mx) || (my && num_cov[y] > my))) continue;
                        num_cov[x] -= 1;
                        num_cov[y] -= 1;
                    }
                    bitfield_clear_bit(temp_bond_valid, b);
                }
            }
        }

# if 0
        // @NOTE(Robin): If this is to be re-enabled, it needs to be verified as it did not produce exactly the expected results
        // without side effects in some dataset
        // 
        // Remove *triangles* i.e. if A-B and B-C are covalent, then A-C should not be covalent, (unless all carbon)
        // In such case remove the longest bond
        for (int bij = 0; bij < (int)temp_bond.count; ++bij) {
            if (!bitfield_test_bit(temp_bond_valid, bij)) {
                continue; // Skip invalidated bonds
            }

            int i = temp_bond.pairs[bij].idx[0];
            int j = temp_bond.pairs[bij].idx[1];
            float dij = temp_bond_dist[bij];

            // Check if atoms i and j have a common neighbor k
            md_bond_iter_t iter_i = md_bond_iter(&temp_bond, i);
            while (md_bond_iter_has_next(&iter_i)) {
                int k = md_bond_iter_atom_index(&iter_i);
                bool all_carbon = (atomic_nr[i] == 6 && atomic_nr[j] == 6 && atomic_nr[k] == 6);

                if (k != j && !all_carbon) {
                    int bik = md_bond_iter_bond_index(&iter_i);
                    float dik = temp_bond_dist[bik];

                    // Find other path between i and j, through k
                    int bjk = md_bond_find(&temp_bond, j, k);
                    if (bjk != -1) {
                        float djk = temp_bond_dist[bjk];
                    
                        // Compare all three distances and flag the longest bond for removal
                        if (dij >= dik && dij >= djk) {
                            bitfield_clear_bit(temp_bond_valid, bij);
                        } else if (dik >= dij && dik >= djk) {
                            bitfield_clear_bit(temp_bond_valid, bik);
                        } else {
                            bitfield_clear_bit(temp_bond_valid, bjk);
                        }
                    }
                }
                md_bond_iter_next(&iter_i);
            }
        }
#endif

        size_t num_valid_bonds = bitfield_popcount(temp_bond_valid, temp_bond.count);

        // Populate real bond data
        md_array_ensure(bond->pairs, num_valid_bonds, alloc);
        md_array_ensure(bond->flags, num_valid_bonds, alloc);

        for (size_t i = 0; i < temp_bond.count; ++i) {
            if (bitfield_test_bit(temp_bond_valid, i)) {
                md_bond_insert(bond, temp_bond.pairs[i].idx[0], temp_bond.pairs[i].idx[1], temp_bond.flags[i], alloc);
            }
        }
    }

    md_temp_end(temp_scope);
}

bool md_bond_infer(md_bond_data_t* out, const md_system_t* sys, const md_system_state_t* state, md_allocator_i* alloc, const md_executor_t* executor) {
    MD_ASSERT(out);
    MD_ASSERT(sys);
    MD_ASSERT(state);
    MD_ASSERT(alloc);
    MD_ASSERT(out != &sys->bond);
    md_bond_data_clear(out);
    if (!state->xyz) {
        MD_LOG_ERROR("Missing atom coordinates");
        return false;
    }
    infer_bonds(out, state, sys, alloc, executor);
    md_bond_build_connectivity(out, state->num_atoms, alloc);
    return true;
}

void md_bond_merge(md_bond_data_t* out, const md_bond_data_t* given, const md_bond_data_t* inferred, size_t atom_count, md_allocator_i* alloc) {
    MD_ASSERT(out);
    MD_ASSERT(given);
    MD_ASSERT(inferred);
    MD_ASSERT(alloc);
    MD_ASSERT(out != given && out != inferred);
    md_bond_data_clear(out);

    md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp_scope);

    // Every given bond which was not inferred is kept, in its order, after the inferred ones
    md_array(md_atom_pair_t)  bond_pairs = 0;
    md_array(md_bond_flags_t) bond_flags = 0;
    for (size_t i = 0; i < given->count; ++i) {
        if (md_bond_origin(given->flags[i]) != MD_BOND_ORIGIN_INFERRED) {
            md_array_push(bond_pairs, given->pairs[i], temp_arena);
            md_array_push(bond_flags, given->flags[i], temp_arena);
        }
    }
    // The non-covalent bonds are not inferred: they are kept as they are, after the graph (md_bond_data_t)
    const size_t num_noncovalent = given->noncovalent_count;
    const md_atom_pair_t*  noncovalent_pairs = given->pairs + given->count;
    const md_bond_flags_t* noncovalent_flags = given->flags + given->count;

    md_bond_data_t* bond = out;
    if (inferred->count) {
        md_array_push_array(bond->pairs, inferred->pairs, inferred->count, alloc);
        md_array_push_array(bond->flags, inferred->flags, inferred->count, alloc);
        bond->count = inferred->count;
    }

    // Reconcile the fresh inferred bonds with the kept ones:
    // - an inferred bond naming a pair that is already kept is a duplicate, the kept one wins
    // - where a topology (md_itp) defines the connectivity it is authoritative, so inferred bonds between two
    //   atoms it covers are dropped, or re-inferring would add the distance guesses back on top of it
    // - a pair with a non-covalent bond keeps it: an atom pair has one bond, and that one was given
    {
        const size_t num_kept  = md_array_size(bond_pairs);
        const size_t num_saved = num_kept + num_noncovalent;
        if (num_saved > 0 && bond->count > 0) {
            const size_t num_atoms = atom_count;
            uint64_t* covered = make_bitfield(num_atoms, temp_arena);

            // An elastic bond is no bond of the molecule: it neither stands in for an inferred bond nor says that a
            // topology covers its atoms
            md_atom_pair_t* saved_sorted = md_temp_alloc_array(temp_scope, md_atom_pair_t, num_saved);
            size_t num_sorted = 0;
            for (size_t i = 0; i < num_saved; ++i) {
                const md_bond_flags_t f = i < num_kept ? bond_flags[i] : noncovalent_flags[i - num_kept];
                if (md_bond_kind(f) == MD_BOND_KIND_ELASTIC) continue;
                md_atom_pair_t pair = i < num_kept ? bond_pairs[i] : noncovalent_pairs[i - num_kept];
                if (pair.idx[0] > pair.idx[1]) { md_atom_idx_t t = pair.idx[0]; pair.idx[0] = pair.idx[1]; pair.idx[1] = t; }
                saved_sorted[num_sorted++] = pair;
                if (i < num_kept && md_bond_origin(f) == MD_BOND_ORIGIN_TOPOLOGY) {
                    for (int k = 0; k < 2; ++k) {
                        if (pair.idx[k] >= 0 && (size_t)pair.idx[k] < num_atoms) bitfield_set_bit(covered, pair.idx[k]);
                    }
                }
            }
            qsort(saved_sorted, num_sorted, sizeof(md_atom_pair_t), compare_atom_pair);

            size_t w = 0;
            for (size_t i = 0; i < bond->count; ++i) {
                const md_atom_pair_t pair = bond->pairs[i];
                if (bitfield_test_bit(covered, pair.idx[0]) && bitfield_test_bit(covered, pair.idx[1])) continue;
                md_atom_pair_t key = pair;
                if (key.idx[0] > key.idx[1]) { md_atom_idx_t t = key.idx[0]; key.idx[0] = key.idx[1]; key.idx[1] = t; }
                if (num_sorted && bsearch(&key, saved_sorted, num_sorted, sizeof(md_atom_pair_t), compare_atom_pair)) continue;
                bond->pairs[w] = pair;
                bond->flags[w] = bond->flags[i];
                w += 1;
            }
            bond->count = w;
            md_array_shrink(bond->pairs, w);
            md_array_shrink(bond->flags, w);
        }
    }

    for (size_t i = 0; i < md_array_size(bond_pairs); ++i) {
        md_array_push(bond->pairs, bond_pairs[i], alloc);
        md_array_push(bond->flags, bond_flags[i], alloc);
        bond->count += 1;
    }
    for (size_t i = 0; i < num_noncovalent; ++i) {
        md_array_push(bond->pairs, noncovalent_pairs[i], alloc);
        md_array_push(bond->flags, noncovalent_flags[i], alloc);
    }
    bond->noncovalent_count = num_noncovalent;
    md_bond_build_connectivity(bond, atom_count, alloc);

    md_temp_end(temp_scope);
}

size_t md_system_infer_coordination(md_system_t* sys) {
    MD_ASSERT(sys);
    size_t count = 0;
    if (!sys->bond.flags) return 0;
    for (size_t i = 0; i < sys->bond.count; ++i) {
        const md_atom_pair_t p = sys->bond.pairs[i];
        if (p.idx[0] < 0 || p.idx[1] < 0 || (size_t)p.idx[0] >= sys->atom.count || (size_t)p.idx[1] >= sys->atom.count) continue;
        if (!md_bond_is_graph(sys->bond.flags[i])) continue;   // A spring, of no chemistry
        const bool metal_a = is_metal(md_atom_atomic_number(&sys->atom, p.idx[0]));
        const bool metal_b = is_metal(md_atom_atomic_number(&sys->atom, p.idx[1]));
        // A kind which was given stays as it was given: a metal-ligand bond a user says is covalent
        const md_bond_flags_t f = sys->bond.flags[i];
        if (metal_a != metal_b && md_bond_kind(f) == MD_BOND_KIND_COVALENT && !(f & MD_BOND_FLAG_KIND_GIVEN)) {
            sys->bond.flags[i] = md_bond_flags_set_kind(f, MD_BOND_KIND_COORDINATE);
            count += 1;
        }
    }
    if (count) md_system_topology_changed(sys);
    return count;
}

// ### BONDS FORMED AND BROKEN ###

typedef struct delta_search_t {
    const md_bitset_t* scope;   // NULL when every atom is looked at, and the pairs come once each
    md_array(bond_pair_t) pairs;
    md_allocator_i* alloc;
} delta_search_t;

static void delta_pair_callback(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    delta_search_t* s = (delta_search_t*)user_param;
    for (size_t k = 0; k < num_pairs; ++k) {
        const uint32_t i = i_idx[k];
        const uint32_t j = j_idx[k];
        if (i == j) continue;
        // i is of the scope: a pair with both atoms in it comes twice, and is taken as the one with the lower atom first
        if (s->scope && j < i && md_bitset_test(*s->scope, j)) continue;
        const bond_pair_t p = { (int)MD_MIN(i, j), (int)MD_MAX(i, j), sqrtf(ij_dist2[k]) };
        md_array_push(s->pairs, p, s->alloc);
    }
}

// What the delta needs to know of an atom, looked up as atoms come up rather than for every atom of the system
typedef struct delta_atoms_t {
    const md_system_t* sys;
    int8_t* res_water;          // [residue count] -1 not looked at yet, 0, 1
    const uint8_t* metallic;    // [atom count] Of a metallic phase: an atom of a metallic bond (md_bond_metallic_atoms)
} delta_atoms_t;

static inline int delta_res(const delta_atoms_t* a, uint32_t atom) {
    return a->sys->residue.count ? (int)md_system_residue_find_by_atom_idx(a->sys, atom) : -1;
}

static inline bool delta_water(delta_atoms_t* a, int res) {
    if (res < 0) return false;
    if (a->res_water[res] < 0) {
        a->res_water[res] = (int8_t)(md_system_residue_kind(a->sys, res) == MD_RESIDUE_KIND_WATER || atoms_are_water(a->sys, md_system_residue_atom_range(a->sys, res)));
    }
    return a->res_water[res] != 0;
}

// The cutoff of a pair (factor times the sum of the radii, 0 when it is never bonded) and whether it is coordination
static inline bond_rule_t delta_rule(delta_atoms_t* a, uint32_t i, uint32_t j, float* out_sum_r) {
    const md_atomic_number_t zi = md_atom_atomic_number(&a->sys->atom, i);
    const md_atomic_number_t zj = md_atom_atomic_number(&a->sys->atom, j);
    *out_sum_r = element_covalent_radius(zi) + element_covalent_radius(zj);
    const bond_rule_t none = {0};
    if (zi == 0 || zj == 0) return none;
    if (md_atom_particle_kind(&a->sys->atom, i) == MD_PARTICLE_VIRTUAL_SITE || md_atom_particle_kind(&a->sys->atom, j) == MD_PARTICLE_VIRTUAL_SITE) return none;
    // The bonds within a metallic phase are its metallic bonds, which do not change (METALLIC PHASES)
    if (a->metallic[i] && a->metallic[j]) return none;
    const int ci = delta_res(a, i);
    const int cj = delta_res(a, j);
    return bond_rule(zi, zj, ci == cj, delta_water(a, ci) || delta_water(a, cj));
}

static inline float delta_distance(const md_system_state_t* state, uint32_t i, uint32_t j) {
    md_float3_t dx = md_float3_sub(state->xyz[j], state->xyz[i]);
    md_geom_min_image(dx.elem, 0, 1, &state->unitcell);
    return md_float3_length(dx);
}

static int compare_u64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
    return (x > y) - (x < y);
}

static inline uint64_t delta_pair_key(md_atom_pair_t p) {
    return ((uint64_t)(uint32_t)MD_MIN(p.idx[0], p.idx[1]) << 32) | (uint32_t)MD_MAX(p.idx[0], p.idx[1]);
}

// The position of atom in the sorted atoms, which holds it
static inline size_t delta_atom_slot(const uint32_t* atoms, size_t count, uint32_t atom) {
    size_t lo = 0, hi = count;
    while (lo < hi) {
        const size_t mid = (lo + hi) / 2;
        if (atoms[mid] < atom) lo = mid + 1;
        else hi = mid;
    }
    MD_ASSERT(lo < count && atoms[lo] == atom);
    return lo;
}

static int compare_u32(const void* a, const void* b) {
    const uint32_t x = *(const uint32_t*)a, y = *(const uint32_t*)b;
    return (x > y) - (x < y);
}

void md_bond_delta_free(md_bond_delta_t* delta, md_allocator_i* alloc) {
    MD_ASSERT(delta);
    MD_ASSERT(alloc);
    md_bitset_free(&delta->broken, alloc);
    if (delta->formed_pair)  md_free(alloc, delta->formed_pair,  sizeof(md_atom_pair_t)  * delta->num_formed);
    if (delta->formed_flags) md_free(alloc, delta->formed_flags, sizeof(md_bond_flags_t) * delta->num_formed);
    md_memset(delta, 0, sizeof(md_bond_delta_t));
}

bool md_bond_delta_compute(md_bond_delta_t* out, const md_system_t* sys, const md_system_state_t* state, const md_bond_delta_desc_t* desc, md_allocator_i* alloc) {
    MD_ASSERT(out);
    MD_ASSERT(sys);
    MD_ASSERT(state);
    MD_ASSERT(desc);
    MD_ASSERT(alloc);
    md_memset(out, 0, sizeof(md_bond_delta_t));

    const size_t num_atoms = sys->atom.count;
    if (!state->xyz || state->num_atoms != num_atoms || num_atoms == 0) return false;
    if (md_system_is_coarse_grained(sys)) return false;
    if (sys->bond.count > 0 && !sys->bond.conn.offset) {
        MD_LOG_ERROR("The bonds of the system have no connectivity (md_bond_build_connectivity)");
        return false;
    }
    const float tol = MD_CLAMP(desc->tolerance, 0.0f, 0.49f);
    const md_bitset_t* scope = desc->scope;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    delta_atoms_t atoms = { .sys = sys };
    atoms.res_water = md_temp_alloc_array(temp, int8_t, MD_MAX(sys->residue.count, 1));
    md_memset(atoms.res_water, -1, MD_MAX(sys->residue.count, 1));
    uint8_t* metallic = md_temp_alloc_array(temp, uint8_t, num_atoms);
    md_bond_metallic_atoms(metallic, &sys->bond, num_atoms);
    atoms.metallic = metallic;

    // ## The bonds of the graph which touch the scope: broken beyond the band, and stretched beyond the cutoff
    // (they no longer count against what their atoms can take)
    const size_t num_bonds = sys->bond.count;
    md_bitmask_t broken    = md_bitmask_create(0, (uint32_t)num_bonds, temp_alloc);
    md_bitmask_t stretched = md_bitmask_create(0, (uint32_t)num_bonds, temp_alloc);
    {
        md_bitset_iter_t it = scope ? md_bitset_iter(*scope) : md_bitset_iter(md_bitset_range(0, (uint32_t)num_atoms));
        while (md_bitset_iter_next(&it)) {
            const uint32_t a = (uint32_t)it.idx;
            if (a >= num_atoms) break;
            md_bond_iter_t bit = md_bond_iter(&sys->bond, a);
            while (md_bond_iter_has_next(&bit)) {
                const uint32_t o = (uint32_t)md_bond_iter_atom_index(&bit);
                const uint32_t b = (uint32_t)md_bond_iter_bond_index(&bit);
                md_bond_iter_next(&bit);
                // Each bond once: from its lower atom when both are in the scope
                const bool o_in = scope ? md_bitset_test(*scope, o) : true;
                if (o_in && o < a) continue;
                float sum_r = 0.0f;
                const bond_rule_t rule = delta_rule(&atoms, a, o, &sum_r);
                if (sum_r <= 0.0f) continue;
                const md_atomic_number_t za = md_atom_atomic_number(&sys->atom, a);
                const md_atomic_number_t zo = md_atom_atomic_number(&sys->atom, o);
                if (za == 0 || zo == 0) continue;
                // A pair inference would not bond (an alkali ion and its ligand from a file, ...) breaks beyond the reach
                // of the search for candidates
                const float cutoff = (rule.factor > 0.0f ? rule.factor : bond_k_coord) * sum_r;
                const float d = delta_distance(state, a, o);
                if (d > cutoff) md_bitmask_set_bit(&stretched, b);
                const uint32_t origin_bit = 1u << md_bond_origin(sys->bond.flags[b]);
                if ((desc->breakable & origin_bit) && d > (1.0f + tol) * cutoff) {
                    md_bitmask_set_bit(&broken, b);
                }
            }
        }
    }

    // ## Candidates for the formed bonds: the pairs with an atom of the scope within reach
    md_array(bond_pair_t) formed = 0;
    {
        // The reach of the largest rule, between the largest atoms of the scope and of the system
        float max_rad_scope = 0.0f, max_rad_all = 0.0f;
        for (size_t i = 0; i < num_atoms; ++i) {
            const float r = element_covalent_radius(md_atom_atomic_number(&sys->atom, i));
            max_rad_all = MD_MAX(max_rad_all, r);
            if (!scope || md_bitset_test(*scope, (uint32_t)i)) max_rad_scope = MD_MAX(max_rad_scope, r);
        }
        const double cutoff = MD_MAX(1.0, (double)(max_rad_scope + max_rad_all) * bond_k_coord);

        delta_search_t search = { .scope = scope, .alloc = temp_alloc };
        md_spatial_acc_t acc = { .alloc = temp_alloc };
        md_spatial_acc_init(&acc, &(md_spatial_acc_desc_t){ .xyz = (const float*)state->xyz, .stride = sizeof(md_float3_t), .count = num_atoms, .cutoff = cutoff, .unitcell = &state->unitcell });
        if (!scope) {
            md_spatial_acc_for_each_internal_pair_within_cutoff(&acc, cutoff, delta_pair_callback, &search);
        } else {
            const size_t num_scope = md_bitset_count(*scope);
            if (num_scope > 0) {
                int* idx = md_temp_alloc_array(temp, int, num_scope);
                size_t n = 0;
                md_bitset_iter_t it = md_bitset_iter(*scope);
                while (md_bitset_iter_next(&it) && n < num_scope) {
                    if (it.idx < num_atoms) idx[n++] = (int)it.idx;
                }
                md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&acc, (const float*)state->xyz, sizeof(md_float3_t), idx, n, cutoff, delta_pair_callback, &search, MD_SPATIAL_ACC_FLAG_USE_SOURCE_IDX);
            }
        }

        // The non-covalent bonds of the system, whose pairs form nothing
        const size_t num_nc = sys->bond.noncovalent_count;
        uint64_t* nc_keys = md_temp_alloc_array(temp, uint64_t, MD_MAX(num_nc, 1));
        for (size_t k = 0; k < num_nc; ++k) nc_keys[k] = delta_pair_key(sys->bond.pairs[num_bonds + k]);
        qsort(nc_keys, num_nc, sizeof(uint64_t), compare_u64);

        for (size_t k = 0; k < md_array_size(search.pairs); ++k) {
            const bond_pair_t p = search.pairs[k];
            float sum_r = 0.0f;
            const bond_rule_t rule = delta_rule(&atoms, (uint32_t)p.atom_i, (uint32_t)p.atom_j, &sum_r);
            if (rule.factor <= 0.0f) continue;
            if (p.dist <= bond_k_min * sum_r || p.dist >= (1.0f - tol) * rule.factor * sum_r) continue;
            if (md_bond_find(&sys->bond, p.atom_i, p.atom_j) != -1) continue;
            const md_atom_pair_t pair = {{ p.atom_i, p.atom_j }};
            const uint64_t key = delta_pair_key(pair);
            if (num_nc && bsearch(&key, nc_keys, num_nc, sizeof(uint64_t), compare_u64)) continue;
            bond_pair_t f = p;
            // The kind rides in the sign of the distance until the flags are written
            if (rule.coordinate) f.dist = -f.dist;
            md_array_push(formed, f, temp_alloc);
        }
    }

    // ## An organohalogen coordinates nothing: a halogen bonded to a non-metal in the state, by a bond of the graph which
    // is not broken in it or by one formed in it (is_halogen)
    if (md_array_size(formed)) {
        md_array(uint32_t) halogens = 0;    // With a covalent bond formed
        for (size_t k = 0; k < md_array_size(formed); ++k) {
            if (formed[k].dist < 0.0f) continue;
            const uint32_t i = (uint32_t)formed[k].atom_i;
            const uint32_t j = (uint32_t)formed[k].atom_j;
            const md_atomic_number_t zi = md_atom_atomic_number(&sys->atom, i);
            const md_atomic_number_t zj = md_atom_atomic_number(&sys->atom, j);
            if (is_halogen(zi) && !is_metal(zj)) md_array_push(halogens, i, temp_alloc);
            if (is_halogen(zj) && !is_metal(zi)) md_array_push(halogens, j, temp_alloc);
        }
        const size_t num_halogens = md_array_size(halogens);
        if (num_halogens) qsort(halogens, num_halogens, sizeof(uint32_t), compare_u32);

        size_t w = 0;
        for (size_t k = 0; k < md_array_size(formed); ++k) {
            if (formed[k].dist < 0.0f) {
                const uint32_t i = (uint32_t)formed[k].atom_i;
                const uint32_t j = (uint32_t)formed[k].atom_j;
                const uint32_t h = is_halogen(md_atom_atomic_number(&sys->atom, i)) ? i : j;
                if (is_halogen(md_atom_atomic_number(&sys->atom, h))) {
                    bool organo = num_halogens && bsearch(&h, halogens, num_halogens, sizeof(uint32_t), compare_u32);
                    md_bond_iter_t it = md_bond_iter(&sys->bond, h);
                    for (; md_bond_iter_has_next(&it) && !organo; md_bond_iter_next(&it)) {
                        const uint32_t o = (uint32_t)md_bond_iter_atom_index(&it);
                        const uint32_t b = (uint32_t)md_bond_iter_bond_index(&it);
                        organo = md_bond_kind((md_bond_flags_t)md_bond_iter_bond_flags(&it)) == MD_BOND_KIND_COVALENT &&
                                 !is_metal(md_atom_atomic_number(&sys->atom, o)) && !md_bitmask_test_bit(&broken, b);
                    }
                    if (organo) continue;
                }
            }
            formed[w++] = formed[k];
        }
        md_array_shrink(formed, w);
    }

    // ## What the atoms of the formed bonds can take. The bonds of the graph count unless broken or stretched, and
    // only a formed bond is ever removed.
    const size_t num_candidates = md_array_size(formed);
    md_array(bool) keep = 0;
    if (num_candidates) {
        md_array_resize(keep, num_candidates, temp_alloc);
        for (size_t k = 0; k < num_candidates; ++k) keep[k] = true;

        uint32_t* ends = md_temp_alloc_array(temp, uint32_t, 2 * num_candidates);
        for (size_t k = 0; k < num_candidates; ++k) {
            ends[2 * k + 0] = (uint32_t)formed[k].atom_i;
            ends[2 * k + 1] = (uint32_t)formed[k].atom_j;
        }
        qsort(ends, 2 * num_candidates, sizeof(uint32_t), compare_u32);
        size_t num_ends = 0;
        for (size_t k = 0; k < 2 * num_candidates; ++k) {
            if (num_ends == 0 || ends[num_ends - 1] != ends[k]) ends[num_ends++] = ends[k];
        }
        int* num_cov   = md_temp_alloc_array(temp, int, num_ends);
        int* num_coord = md_temp_alloc_array(temp, int, num_ends);
        md_memset(num_cov,   0, sizeof(int) * num_ends);
        md_memset(num_coord, 0, sizeof(int) * num_ends);

        for (size_t e = 0; e < num_ends; ++e) {
            const uint32_t a = ends[e];
            md_bond_iter_t bit = md_bond_iter(&sys->bond, a);
            while (md_bond_iter_has_next(&bit)) {
                const uint32_t o = (uint32_t)md_bond_iter_atom_index(&bit);
                const uint32_t b = (uint32_t)md_bond_iter_bond_index(&bit);
                const md_bond_flags_t f = (md_bond_flags_t)md_bond_iter_bond_flags(&bit);
                md_bond_iter_next(&bit);
                if (md_bitmask_test_bit(&broken, b)) continue;
                // A bond of an atom outside of the scope has not been measured
                bool stretch = md_bitmask_test_bit(&stretched, b);
                if (!stretch && scope && !md_bitset_test(*scope, a) && !md_bitset_test(*scope, o)) {
                    float sum_r = 0.0f;
                    const bond_rule_t rule = delta_rule(&atoms, a, o, &sum_r);
                    stretch = sum_r > 0.0f && delta_distance(state, a, o) > (rule.factor > 0.0f ? rule.factor : bond_k_coord) * sum_r;
                }
                if (stretch) continue;
                if (md_bond_kind(f) == MD_BOND_KIND_COORDINATE) num_coord[e] += 1;
                else                             num_cov[e]   += 1;
            }
        }
        for (size_t k = 0; k < num_candidates; ++k) {
            const size_t ei = delta_atom_slot(ends, num_ends, (uint32_t)formed[k].atom_i);
            const size_t ej = delta_atom_slot(ends, num_ends, (uint32_t)formed[k].atom_j);
            int* n = formed[k].dist < 0.0f ? num_coord : num_cov;
            n[ei] += 1;
            n[ej] += 1;
        }

        // In the order of md_bond_infer: bonds between residues apart first, then between residues next
        // to each other, then the furthest from the radii (bond_excess_t)
        md_array(bond_excess_t) over = 0;
        for (size_t k = 0; k < num_candidates; ++k) {
            const uint32_t i = (uint32_t)formed[k].atom_i;
            const uint32_t j = (uint32_t)formed[k].atom_j;
            const size_t ei = delta_atom_slot(ends, num_ends, i);
            const size_t ej = delta_atom_slot(ends, num_ends, j);
            const md_atomic_number_t zi = md_atom_atomic_number(&sys->atom, i);
            const md_atomic_number_t zj = md_atom_atomic_number(&sys->atom, j);
            bool exceeds = false;
            if (formed[k].dist < 0.0f) {
                const size_t m = is_metal(zi) ? ei : ej;
                exceeds = num_coord[m] > max_coordination(is_metal(zi) ? zi : zj);
            } else {
                const int mi = max_covalent_bonds(zi);
                const int mj = max_covalent_bonds(zj);
                exceeds = (mi && num_cov[ei] > mi) || (mj && num_cov[ej] > mj);
            }
            if (exceeds) {
                const float sum_r = element_covalent_radius(zi) + element_covalent_radius(zj);
                bond_excess_t x = {
                    .inter = bond_separation(delta_res(&atoms, i), delta_res(&atoms, j)),
                    .deviation = sum_r > 0 ? fabsf(fabsf(formed[k].dist) / sum_r - 1.0f) : 0.0f,
                    .idx = (uint32_t)k,
                };
                md_array_push(over, x, temp_alloc);
            }
        }
        const size_t num_over = md_array_size(over);
        if (num_over) {
            qsort(over, num_over, sizeof(bond_excess_t), compare_bond_excess);
            for (size_t o = 0; o < num_over; ++o) {
                const uint32_t k = over[o].idx;
                const uint32_t i = (uint32_t)formed[k].atom_i;
                const uint32_t j = (uint32_t)formed[k].atom_j;
                const size_t ei = delta_atom_slot(ends, num_ends, i);
                const size_t ej = delta_atom_slot(ends, num_ends, j);
                const md_atomic_number_t zi = md_atom_atomic_number(&sys->atom, i);
                const md_atomic_number_t zj = md_atom_atomic_number(&sys->atom, j);
                if (formed[k].dist < 0.0f) {
                    const size_t m = is_metal(zi) ? ei : ej;
                    if (num_coord[m] <= max_coordination(is_metal(zi) ? zi : zj)) continue;
                    num_coord[ei] -= 1;
                    num_coord[ej] -= 1;
                } else {
                    const int mi = max_covalent_bonds(zi);
                    const int mj = max_covalent_bonds(zj);
                    if (!((mi && num_cov[ei] > mi) || (mj && num_cov[ej] > mj))) continue;
                    num_cov[ei] -= 1;
                    num_cov[ej] -= 1;
                }
                keep[k] = false;
            }
        }
    }

    // ## The delta, in the caller's allocator
    out->broken = md_bitmask_extract(&broken, alloc);
    size_t num_formed = 0;
    for (size_t k = 0; k < num_candidates; ++k) num_formed += keep[k];
    if (num_formed) {
        // Sorted by pair: the candidates came in the order of the search
        uint64_t* keys = md_temp_alloc_array(temp, uint64_t, num_formed);
        size_t n = 0;
        for (size_t k = 0; k < num_candidates; ++k) {
            if (!keep[k]) continue;
            const md_atom_pair_t pair = {{ formed[k].atom_i, formed[k].atom_j }};
            keys[n++] = delta_pair_key(pair);
        }
        qsort(keys, num_formed, sizeof(uint64_t), compare_u64);
        out->formed_pair  = md_alloc(alloc, sizeof(md_atom_pair_t)  * num_formed);
        out->formed_flags = md_alloc(alloc, sizeof(md_bond_flags_t) * num_formed);
        const md_bond_flags_t inferred = md_bond_flags_set_origin(MD_BOND_FLAG_NONE, MD_BOND_ORIGIN_INFERRED);
        for (size_t k = 0; k < num_formed; ++k) {
            const md_atom_pair_t pair = {{ (md_atom_idx_t)(keys[k] >> 32), (md_atom_idx_t)(keys[k] & 0xFFFFFFFFu) }};
            out->formed_pair[k] = pair;
            // A rule gives a coordination bond for a metal and a non-metal, and only for them
            const md_atomic_number_t zi = md_atom_atomic_number(&sys->atom, pair.idx[0]);
            const md_atomic_number_t zj = md_atom_atomic_number(&sys->atom, pair.idx[1]);
            out->formed_flags[k] = (is_metal(zi) != is_metal(zj)) ? md_bond_flags_set_kind(inferred, MD_BOND_KIND_COORDINATE) : inferred;
        }
        out->num_formed = num_formed;
    }

    md_bitmask_free(&broken);
    md_bitmask_free(&stretched);
    md_temp_end(temp);
    return true;
}

// ### ELASTIC BONDS ###
//
// The chemistry of a coarse grained chain is that of mdi_infer_coarse_grained_bonds: between two of its residues there is the
// link from the last backbone bead of one to the first of the next, and a disulfide between the SC1 of two cysteines. A
// chain goes on while the residues are of one kind and their numbering does not start over.
static inline bool cg_chain_restarts(md_sequence_id_t prev_seq, md_sequence_id_t seq) {
    return seq < prev_seq && !(seq <= 1 && prev_seq >= 9999);  // GRO and PDB wrap at 99999 and 9999
}

static inline bool cg_cysteine_sc1(const md_system_t* sys, int32_t res_idx, int32_t atom_idx) {
    const md_str_t res_name = md_system_residue_name(sys, res_idx);
    return (md_str_eq_cstr(res_name, "CYS") || md_str_eq_cstr(res_name, "CYX")) && md_str_eq_cstr(md_atom_name(&sys->atom, atom_idx), "SC1");
}

size_t md_system_infer_elastic_bonds(md_system_t* sys) {
    MD_ASSERT(sys);
    const size_t N = sys->atom.count;
    const size_t C = sys->residue.count;
    if (sys->bond.count == 0 || !sys->bond.flags || N == 0 || C == 0 || !sys->residue.atom_offset) return 0;
    if (!md_system_is_coarse_grained(sys)) return 0;

    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* arena = md_temp_allocator(temp);

    int32_t* res_of = md_alloc(arena, sizeof(int32_t) * N);
    for (size_t i = 0; i < N; ++i) res_of[i] = -1;
    // The beads a residue of a chain links by: the first and the last bead of its backbone, and in a nucleotide of
    // Martini 3 its BB2, which is linked to the next residue as well
    int32_t* first_bb = md_alloc(arena, sizeof(int32_t) * C);
    int32_t* last_bb  = md_alloc(arena, sizeof(int32_t) * C);
    int32_t* also_bb  = md_alloc(arena, sizeof(int32_t) * C);
    const uint32_t version = mdi_martini_system(sys, label_of_atom, &sys->atom).version;
    const bool m3 = (version & (MARTINI_2 | MARTINI_2P)) == 0;
    for (size_t ci = 0; ci < C; ++ci) {
        first_bb[ci] = last_bb[ci] = also_bb[ci] = -1;
        const md_urange_t range = md_system_residue_atom_range(sys, ci);
        const md_residue_kind_t kind = md_system_residue_kind(sys, ci);
        const bool biopolymer = kind == MD_RESIDUE_KIND_AMINO_ACID || kind == MD_RESIDUE_KIND_NUCLEOTIDE;
        for (uint32_t i = range.beg; i < range.end && i < N; ++i) {
            res_of[i] = (int32_t)ci;
            if (!biopolymer || md_atom_particle_kind(&sys->atom, i) != MD_PARTICLE_BEAD) continue;
            if (!(md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_BACKBONE)) continue;
            if (first_bb[ci] == -1) first_bb[ci] = (int32_t)i;
            last_bb[ci] = (int32_t)i;
            if (m3 && kind == MD_RESIDUE_KIND_NUCLEOTIDE && md_str_eq_cstr(md_atom_name(&sys->atom, i), "BB2")) also_bb[ci] = (int32_t)i;
        }
    }

    // The links of the chains first: they say where a chain begins and ends, which the bond closing a cyclic one needs
    uint8_t* has_prev = md_alloc(arena, C);
    uint8_t* has_next = md_alloc(arena, C);
    md_memset(has_prev, 0, C);
    md_memset(has_next, 0, C);
    int32_t* chain = md_alloc(arena, sizeof(int32_t) * C);
    for (size_t ci = 0; ci < C; ++ci) chain[ci] = (int32_t)ci;
    md_array(uint32_t) candidates = 0;

    for (size_t b = 0; b < sys->bond.count; ++b) {
        if (!md_bond_is_graph(sys->bond.flags[b])) continue;
        int32_t x = sys->bond.pairs[b].idx[0];
        int32_t y = sys->bond.pairs[b].idx[1];
        if (x < 0 || y < 0 || (size_t)x >= N || (size_t)y >= N) continue;
        if (res_of[x] > res_of[y]) { const int32_t t = x; x = y; y = t; }
        const int32_t ca = res_of[x];
        const int32_t cb = res_of[y];
        if (ca < 0 || ca == cb) continue;  // Within a residue
        if (md_atom_particle_kind(&sys->atom, x) != MD_PARTICLE_BEAD || md_atom_particle_kind(&sys->atom, y) != MD_PARTICLE_BEAD) continue;
        const md_residue_kind_t ka = md_system_residue_kind(sys, ca);
        const md_residue_kind_t kb = md_system_residue_kind(sys, cb);
        if (!(ka == MD_RESIDUE_KIND_AMINO_ACID || ka == MD_RESIDUE_KIND_NUCLEOTIDE) || ka != kb) continue;

        const bool link = cb == ca + 1 && (x == last_bb[ca] || x == also_bb[ca]) && y == first_bb[cb] &&
                          !cg_chain_restarts(md_system_residue_seq_id(sys, ca), md_system_residue_seq_id(sys, cb));
        if (link) {
            has_next[ca] = 1;
            has_prev[cb] = 1;
            uf_union(chain, ca, cb);
        } else if (ka == MD_RESIDUE_KIND_AMINO_ACID && cg_cysteine_sc1(sys, ca, x) && cg_cysteine_sc1(sys, cb, y)) {
            // A disulfide
        } else {
            md_array_push(candidates, (uint32_t)b, arena);
        }
    }

    size_t count = 0;
    for (size_t k = 0; k < md_array_size(candidates); ++k) {
        const uint32_t b = candidates[k];
        int32_t x = sys->bond.pairs[b].idx[0];
        int32_t y = sys->bond.pairs[b].idx[1];
        if (res_of[x] > res_of[y]) { const int32_t t = x; x = y; y = t; }
        const int32_t ca = res_of[x];
        const int32_t cb = res_of[y];
        // The bond which closes a cyclic chain: from the first bead of the backbone at its start to the last at its end
        const bool closure = cb > ca + 1 && x == first_bb[ca] && y == last_bb[cb] && !has_prev[ca] && !has_next[cb] && uf_find(chain, ca) == uf_find(chain, cb);
        if (closure) continue;
        // A spring has neither an order nor resonance
        sys->bond.flags[b] = md_bond_flags_set_kind(sys->bond.flags[b] & ~(MD_BOND_FLAG_ORDER_MASK | MD_BOND_FLAG_RESONANCE_MASK), MD_BOND_KIND_ELASTIC);
        count += 1;
    }

    md_temp_end(temp);
    if (count) {
        md_bond_build_connectivity(&sys->bond, N, sys->alloc);
        md_system_topology_changed(sys);
    }
    return count;
}

void md_bond_mask_grow(md_bitmask_t* mask, const md_bond_data_t* bond, size_t atom_count, size_t extent, const md_bitset_t* viable_mask) {
    MD_ASSERT(mask);
    MD_ASSERT(bond);

    if (extent >= 255) {
        MD_LOG_DEBUG("Maximum supported growth extent is 255, the extent will be clamped to 255");
        extent = 255;
    }
    // The initial depth is 1, so we need to add one.
    extent += 1;
    
    const size_t mask_size = md_bitmask_count(mask);
    if (!mask_size) return;

    md_temp_scope_t temp_scope = md_temp_begin_avoid(mask->alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);

    int* indices = md_alloc(temp_alloc, mask_size * sizeof(int));
    const size_t num_indices = md_bitset_iter_extract_indices(indices, mask_size, md_bitmask_iter(mask));
    MD_ASSERT(num_indices == mask_size);

    md_fifo_t queue = md_fifo_create(64, temp_alloc);

    // Zero means 'not reached'. Temp memory is not cleared, so this has to be, or the walk depends on what the
    // allocator handed out last: fresh pages happen to be zero, reused ones are not.
    uint8_t* depth = md_temp_alloc_array(temp_scope, uint8_t, atom_count);
    md_memset(depth, 0, atom_count * sizeof(uint8_t));

    {
        md_bitset_iter_t it = md_bitmask_iter(mask);
        while (md_bitset_iter_next(&it)) {
            int idx = (int)it.idx;
            depth[idx] = 1;
        }
    }
    
    for (size_t j = 0; j < num_indices; ++j) {
        int i = indices[j];

        md_fifo_clear(&queue);
        md_fifo_push(&queue, i);
        
        while (!md_fifo_empty(&queue)) {
            int idx = md_fifo_pop(&queue);
            md_bitmask_set_bit(mask, (uint32_t)idx);

            md_bond_iter_t it = md_bond_iter(bond, idx);
            while (md_bond_iter_has_next(&it)) {
                int next = md_bond_iter_atom_index(&it);
                md_bond_iter_next(&it);
                if (viable_mask && md_bitset_test(*viable_mask, (uint32_t)next) == false) {
                    continue;
                }
                // Only explore the path if it is not already marked and it is within the extent.
                // Or if the current path is 'cheaper'
                const uint8_t dc = depth[idx] + 1;
                const uint8_t dn = depth[next];
                if (dn == 0 || dc < dn) {
                    depth[next] = dc;
                    if (dc <= extent) {
                        md_fifo_push(&queue, next);
                    }
                }
            }
        }
    }

    md_temp_end(temp_scope);
}
