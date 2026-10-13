#include <md_dssp.h>

#include <md_system.h>
#include <md_unitcell.h>
#include <md_util.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_coord_stream.h>
#include <core/md_hash.h>
#include <core/md_log.h>
#include <core/md_spatial_acc.h>
#include <core/md_vec_math.h>

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

// ### SECONDARY STRUCTURE ###
//
// DSSP (Kabsch & Sander 1983, doi:10.1002/bip.360221211) assigns secondary structure from the hydrogen bonds of the
// backbone: an n-turn at i is a bond from the C=O of residue i to the N-H of residue i+n, two consecutive n-turns make a
// helix, and bridges between residues are patterns of bonds between two strands. It is done here in two stages. The
// bonds go into dssp_res_hbonds_t, where every donor keeps its two best acceptors and every acceptor its two best donors
// (lowest energy first, a bond being an energy below -0.5). The patterns (dssp_assign) then read the donors' acceptors
// (dssp_test_bond), the same for atomistic ranges and for traces of coarse grained beads.
//
// The bonds of an atomistic range are those of DSSP: the electrostatic energy between N, H, C and O, of each residue
// whole about its CA and of the two in their nearest images, should the periodic boundary split a residue or a pair.
//
// CG-DSSP. A trace (one bead per residue, the BB bead of Martini) has no atoms to bond, but it has peptide groups, and
// those are what the patterns are about. Peptide k is the bond between beads k and k+1: it holds the C=O of residue k
// and the N-H of residue k+1, which point opposite ways in its plane. Its site is the midpoint of the bond and its
// dipole the direction of the C=O, which the beads do not give, as the peptide plane turns about the bond with the
// conformation: a linear model of the trace around the residue places its C and O (cg_dssp_dipole, fit by
// tools/fit_cg_dssp.py). Peptide k bonds peptide l, the C=O of k to the N-H of l+1, when their sites are within R and
// both dipoles point from k towards l, head to tail, by more than c (the cosine of each dipole with the line between the
// sites). The local pairs, which make the turns (|k - l| <= 4 within a chain), and the others, which make the sheets,
// have their own R and c. The pseudo energy -1 - (cos_k + cos_l) / 2 ranks the bonds of a donor. The pairs and their
// distances are those of the spatial structure (minimum image), and each chain is taken in one image, bead by bonded bead,
// for the frames of its residues.
//
// The model is for Martini 3 BB beads, the center of N, CA, C and O as martinize2 maps them (Martini 2 maps BB to the
// center of mass of the same atoms; ELNEDYN to the CA, which it is not fit for). The positions of C and O are fit on
// atomistic structures mapped to beads (9 proteins, 3272 residues), the thresholds on Martini 3 simulations with elastic
// networks of 3 of them (1 ns), every frame labelled by the DSSP (the atomistic one here) of the structure the model was
// built from, as the force field holds the secondary structure it was built with. On simulations of 2 others (ubiquitin and
// lysozyme, 2 ns) helix, sheet and coil agree with that for 83 % of the residues of a frame, 89 to 93 % of the helices
// and 68 to 78 % of the strands found, where the CA distance criteria of TM-align (Zhang & Skolnick 2005) get 68 to 77 %;
// 84 to 91 % on the simulations of the tuning. A trace costs about a microsecond per residue, like the atomistic DSSP.

// The thresholds of CG-DSSP: distance between the sites (A) and alignment of the dipoles with the line between them
#define CG_DSSP_TURN_DISTANCE    5.25f
#define CG_DSSP_TURN_ALIGNMENT   0.55f
#define CG_DSSP_SHEET_DISTANCE   7.0f
#define CG_DSSP_SHEET_ALIGNMENT  0.3f

// Longer than any bond between the backbone beads of a chain (A): a longer vector between two crosses the boundary
#define CG_DSSP_MAX_BOND         6.0f

#include "md_dssp_cg.inl"

typedef struct dssp_hbond_t {
    float energy;
    uint32_t res_idx;
} dssp_hbond_t;

typedef struct dssp_res_coords_t {
    vec4_t N;
    vec4_t O;
    vec4_t C;
    vec4_t H;
    vec4_t CA;
} dssp_res_coords_t;

typedef struct dssp_res_hbonds_t {
    dssp_hbond_t acc[2];
    dssp_hbond_t don[2];
} dssp_res_hbonds_t;

// Keeps a bond of the donor don_idx to the acceptor acc_idx among the two best of each
static inline void dssp_hbond_insert(dssp_res_hbonds_t res_hbonds[], size_t don_idx, size_t acc_idx, float energy) {
    dssp_res_hbonds_t* don_hbonds = &res_hbonds[don_idx];
    dssp_res_hbonds_t* acc_hbonds = &res_hbonds[acc_idx];

    if (energy < don_hbonds->acc[0].energy) {
        don_hbonds->acc[1] = don_hbonds->acc[0];
        don_hbonds->acc[0].res_idx = (uint32_t)acc_idx;
        don_hbonds->acc[0].energy  = energy;
    } else if (energy < don_hbonds->acc[1].energy) {
        don_hbonds->acc[1].res_idx = (uint32_t)acc_idx;
        don_hbonds->acc[1].energy  = energy;
    }

    if (energy < acc_hbonds->don[0].energy) {
        acc_hbonds->don[1] = acc_hbonds->don[0];
        acc_hbonds->don[0].res_idx = (uint32_t)don_idx;
        acc_hbonds->don[0].energy  = energy;
    } else if (energy < acc_hbonds->don[1].energy) {
        acc_hbonds->don[1].res_idx = (uint32_t)don_idx;
        acc_hbonds->don[1].energy  = energy;
    }
}

// The energy of the bond from donor don_idx to acceptor acc_idx, with the acceptor moved by acc_shift: the lattice vector
// which takes it to the periodic image nearest the donor, zero (exactly) when the two are in one image already
static inline float calc_hbond_energy(dssp_res_hbonds_t res_hbonds[], const dssp_res_coords_t res_coords[], size_t don_idx, size_t acc_idx, vec4_t acc_shift) {
    /*
    float d_NO = vec3_distance(donor->N, acceptor->O);
    float d_HC = vec3_distance(donor->H, acceptor->C);
    float d_NC = vec3_distance(donor->N, acceptor->C);
    float d_OH = vec3_distance(donor->O, acceptor->H);
    double e_ref = 0.084 * (1.0 / d_NO + 1.0 / d_HC - 1.0 / d_NC - 1.0 / d_OH) * 332.0;
    */

    const dssp_res_coords_t* don = &res_coords[don_idx];
    const dssp_res_coords_t* acc = &res_coords[acc_idx];

    float energy = 0.0f;

    const float min_distance_sq = 0.5f * 0.5f;
    const float min_bond_energy = -9.9f;
	const float energy_cutoff   = -0.5f;

    const vec4_t acc_O = vec4_add(acc->O, acc_shift);
    const vec4_t acc_C = vec4_add(acc->C, acc_shift);
    const vec4_t d2 = {
        vec4_distance_squared(don->N, acc_O),
        vec4_distance_squared(don->H, acc_C),
        vec4_distance_squared(don->N, acc_C),
        vec4_distance_squared(don->H, acc_O),
    };

    const int mask_min = vec4_move_mask(vec4_cmp_lt(d2, vec4_set1(min_distance_sq)));
    if (mask_min > 0) {
        energy = min_bond_energy;
    } else {
        const vec4_t inv_d = vec4_rsqrt(d2);
		const vec4_t w = vec4_set(1.0f, 1.0f, -1.0f, -1.0f);
        const float Q  = 27.888f; // 0.084f * 332.0f;
        energy = Q * vec4_reduce_add(vec4_mul(inv_d, w));
    }

    if (energy > energy_cutoff)
		return 0.0f;

    energy = MAX(energy, min_bond_energy);
    dssp_hbond_insert(res_hbonds, don_idx, acc_idx, energy);
    return energy;
}

static inline vec4_t estimate_HN(vec4_t N, vec4_t CA, vec4_t C_prev) {
    const float NH_dist = 1.01f;
    vec4_t U = vec4_normalize(vec4_sub(N, CA));
    vec4_t V = vec4_normalize(vec4_sub(N, C_prev));
    return vec4_add(N, vec4_mul1(vec4_normalize(vec4_add(U, V)), NH_dist));
}

static inline vec4_t estimate_HN_term(vec4_t N, vec4_t CA, vec4_t C) {
    // i == range.beg (N-terminus): estimate previous C by extrapolation
    // C_prev_est = N + (CA - C)
    vec4_t CA_minus_C = vec4_sub(CA, C);
    vec4_t Cprev_est  = vec4_add(N, CA_minus_C);

    // Try to construct H using the standard estimator
    return estimate_HN(N, CA, Cprev_est);
}

static inline bool dssp_test_bond(const dssp_res_hbonds_t res_hbonds[], size_t acceptor, size_t donor) {
    const float max_hbond_energy = -0.5f;
    return (res_hbonds[donor].acc[0].res_idx == acceptor && res_hbonds[donor].acc[0].energy < max_hbond_energy) ||
        (res_hbonds[donor].acc[1].res_idx == acceptor && res_hbonds[donor].acc[1].energy < max_hbond_energy);
}

enum {
    SS_FLAG_NONE = 0,
    SS_FLAG_TURN = 1,
    SS_FLAG_BEND = 2,
    SS_FLAG_HELIX_310 = 4,
    SS_FLAG_HELIX_ALPHA = 8,
    SS_FLAG_HELIX_PI = 16,
    SS_FLAG_SHEET = 32,
    SS_FLAG_BRIDGE = 64,
};

typedef struct residue_pair_t {
    uint32_t i;
    uint32_t j;
} residue_pair_t;

// A bend at the middle of five residues: the vectors from its CA (or bead) to those two before and two after it, in one
// image, meet at less than 70 degrees. The angle is atan2(|v1 x v2|, v1 . v2); below 90 degrees it is less than 70 where
// the cross is less than tan(70) times the dot, which spares the arctangent.
static inline bool dssp_bend(vec4_t to_prev, vec4_t to_next) {
    const float dot = vec4_dot(to_prev, to_next);
    const vec4_t cross = vec4_cross(to_prev, to_next);
    const float tan2_70 = 7.5486322f;   // tan(70 degrees)^2
    return dot >= 0.0f && vec4_dot(cross, cross) <= tan2_70 * (dot * dot);
}

// ### THE BONDS OF ATOMISTIC RANGES ###

// The region of the cell further than a cutoff from each of its periodic faces. Two points in it within the cutoff of
// each other (by minimum image) are in their nearest images as they are given: a periodic image of one is never within
// the cutoff of the other, as the fractions of a vector that short along each axis are within those margins.
typedef struct dssp_interior_t {
    float I[3][3];  // The inverse basis, I[col][row] (cartesian to fractional)
    float lo[3];
    float hi[3];
} dssp_interior_t;

static inline dssp_interior_t dssp_interior_init(const md_unitcell_t* cell, float cutoff) {
    dssp_interior_t t = {0};
    md_unitcell_I_extract_float(t.I, cell);
    const uint32_t flags = md_unitcell_flags(cell);
    for (int k = 0; k < 3; ++k) {
        if (flags & (MD_UNITCELL_PBC_X << k)) {
            // The fraction of a vector of length cutoff along k is at most cutoff over the width of the cell across k,
            // the length of row k of I. A cell less than twice the cutoff across k has no interior.
            const float w_inv = sqrtf(t.I[0][k] * t.I[0][k] + t.I[1][k] * t.I[1][k] + t.I[2][k] * t.I[2][k]);
            const float margin = cutoff * w_inv * 1.001f + 1.0e-5f;
            t.lo[k] = margin;
            t.hi[k] = 1.0f - margin;
        } else {
            t.lo[k] = -FLT_MAX;
            t.hi[k] =  FLT_MAX;
        }
    }
    return t;
}

static inline bool dssp_interior_test(const dssp_interior_t* t, vec3_t p) {
    bool inside = true;
    for (int k = 0; k < 3; ++k) {
        const float f = t->I[0][k] * p.x + t->I[1][k] * p.y + t->I[2][k] * p.z;
        inside &= (t->lo[k] <= f) & (f <= t->hi[k]);
    }
    return inside;
}

typedef struct dssp_hbond_energy_user_param_t {
	const uint32_t* res_key;        // [segment] range_idx << 1, | 1 where the CA is outside the interior of the cell (dssp_interior_t)
    const dssp_res_coords_t* res_coords;
    dssp_res_hbonds_t* res_hbonds;
	md_array(residue_pair_t)* sheet_candidates;
    const md_unitcell_t* cell;
	md_allocator_i* alloc;
} dssp_hbond_energy_user_param_t;

static void dssp_hbond_energy_cb(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    dssp_hbond_energy_user_param_t* param = (dssp_hbond_energy_user_param_t*)user_param;
    for (size_t idx = 0; idx < num_pairs; ++idx) {
        uint32_t i = MIN(i_idx[idx], j_idx[idx]);
        uint32_t j = MAX(i_idx[idx], j_idx[idx]);

		const uint32_t ki = param->res_key[i];
		const uint32_t kj = param->res_key[j];
        const bool same_range = (ki >> 1) == (kj >> 1);

        int d = (int)j - (int)i;

        if (same_range && d <= 1) continue;

        // The acceptor in its image nearest the donor. The query measured the pair by minimum image, and their CA as
        // given are that far apart unless the pair crosses the periodic boundary: only then is j moved, by the lattice
        // vector to its image nearest i. Every other pair keeps its coordinates as they are, bit for bit. A pair in
        // the interior of the cell cannot cross.
        vec4_t shift = vec4_zero();
        if ((ki | kj) & 1) {
            const vec4_t raw = vec4_sub(param->res_coords[j].CA, param->res_coords[i].CA);
            if (vec4_dot(raw, raw) > ij_dist2[idx] * 1.0001f + 1.0e-4f) {
                vec4_t img = raw;
                md_util_min_image_vec4(&img, 1, param->cell);
                shift = vec4_sub(img, raw);
            }
        }

        calc_hbond_energy(param->res_hbonds, param->res_coords, i, j, shift);
        calc_hbond_energy(param->res_hbonds, param->res_coords, j, i, vec4_sub(vec4_zero(), shift));

        // DSSP excludes very local sequence contacts on the same chain
        if (same_range && d <= 2) continue;

		// Ensure that pairs are sorted (i < j) to avoid duplicates
		residue_pair_t pair = { i, j };
		md_array_push(*param->sheet_candidates, pair, param->alloc);
    }
}

// p in the image of ref, for atoms bonded within a residue (within 3 A): only a longer vector crosses the periodic boundary
// and is taken to its minimum image, the others stay as they are, bit for bit
static inline vec4_t dssp_whole_about(vec4_t ref, vec4_t p, const md_unitcell_t* cell) {
    vec4_t d = vec4_sub(p, ref);
    if (vec4_dot(d, d) <= 3.0f * 3.0f) return p;
    md_util_min_image_vec4(&d, 1, cell);
    return vec4_add(ref, d);
}

// The bonds and the bends of the atomistic ranges (num_trace residues are in traces)
static void dssp_atomistic_bonds(dssp_res_hbonds_t res_hbonds[], md_array(residue_pair_t)* sheet_candidates, uint8_t res_bend[], const vec3_t* xyz, const md_unitcell_t* cell,
    const md_protein_backbone_data_t* backbone, size_t num_trace, md_allocator_i* temp_alloc)
{
    const size_t backbone_segment_count = backbone->segment.count;
    const md_amino_acid_atoms_t* backbone_atoms = backbone->segment.atoms;
    dssp_res_coords_t* res_coords = md_alloc(temp_alloc, sizeof(dssp_res_coords_t) * backbone_segment_count);
    vec3_t* res_ca = md_alloc(temp_alloc, sizeof(vec3_t) * backbone_segment_count);
    uint32_t* res_key = md_alloc(temp_alloc, sizeof(uint32_t) * backbone_segment_count);
    const float cutoff = 9.0f;
    const dssp_interior_t interior = dssp_interior_init(cell, cutoff);
    // Alongside traces, the query runs over the atomistic residues by their indices
    int* res_idx = num_trace ? md_alloc(temp_alloc, sizeof(int) * (backbone_segment_count - num_trace)) : NULL;
    size_t num_res = 0;

    for (size_t range_idx = 0; range_idx < backbone->range.count; ++range_idx) {
        if (md_protein_backbone_range_is_trace(backbone, range_idx)) continue;
        const md_irange_t range = {(int)backbone->range.offset[range_idx], (int)backbone->range.offset[range_idx + 1]};

        for (int i = range.beg; i < range.end; ++i) {
            md_atom_idx_t ca_idx = backbone_atoms[i].ca;
            md_atom_idx_t n_idx  = backbone_atoms[i].n;
            md_atom_idx_t c_idx  = backbone_atoms[i].c;
            md_atom_idx_t o_idx  = backbone_atoms[i].o;
            md_atom_idx_t hn_idx = backbone_atoms[i].hn;

            res_ca[i] = xyz[ca_idx];
            res_key[i] = (uint32_t)range_idx << 1 | (dssp_interior_test(&interior, xyz[ca_idx]) ? 0 : 1);

            // The residue whole, in the image of its CA, should the periodic boundary split it
            const vec4_t CA = vec4_from_vec3(xyz[ca_idx], 0);
            res_coords[i].CA = CA;
            res_coords[i].N  = dssp_whole_about(CA, vec4_from_vec3(xyz[n_idx], 0), cell);
            res_coords[i].C  = dssp_whole_about(CA, vec4_from_vec3(xyz[c_idx], 0), cell);
            res_coords[i].O  = dssp_whole_about(CA, vec4_from_vec3(xyz[o_idx], 0), cell);

            if (hn_idx >= 0) {
                res_coords[i].H = dssp_whole_about(CA, vec4_from_vec3(xyz[hn_idx], 0), cell);
            } else {
                if (i > range.beg)
                    res_coords[i].H = estimate_HN(res_coords[i].N, res_coords[i].CA, dssp_whole_about(res_coords[i].N, res_coords[i - 1].C, cell));
                else
                    res_coords[i].H = estimate_HN_term(res_coords[i].N, res_coords[i].CA, res_coords[i].C);
            }

            if (res_idx) res_idx[num_res] = i;
            num_res += 1;
        }

        // Bends. Residues two apart are never more than 8 A apart, unless split by the boundary, and in the interior
        // of the cell they cannot be.
        for (int i = range.beg + 2; i + 2 < range.end; ++i) {
            vec4_t v1 = vec4_sub(res_coords[i - 2].CA, res_coords[i].CA);
            vec4_t v2 = vec4_sub(res_coords[i + 2].CA, res_coords[i].CA);
            if ((res_key[i - 2] | res_key[i] | res_key[i + 2]) & 1) {
                if (vec4_dot(v1, v1) > 8.0f * 8.0f) md_util_min_image_vec4(&v1, 1, cell);
                if (vec4_dot(v2, v2) > 8.0f * 8.0f) md_util_min_image_vec4(&v2, 1, cell);
            }
            res_bend[i] = dssp_bend(v1, v2);
        }
    }

    md_coord_stream_t stream = md_coord_stream_from_aos((const float*)res_ca, sizeof(vec3_t), res_idx, num_res);
    md_spatial_acc_t acc = { .alloc = temp_alloc };
    md_spatial_acc_init(&acc, &(md_spatial_acc_desc_t){ .coords = &stream, .cutoff = cutoff, .unitcell = cell, .flags = res_idx ? MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX : MD_SPATIAL_ACC_FLAG_NONE });

    dssp_hbond_energy_user_param_t user_param = {
		.res_key = res_key,
        .res_coords = res_coords,
        .res_hbonds = res_hbonds,
		.sheet_candidates = sheet_candidates,
        .cell = cell,
		.alloc = temp_alloc,
	};
    md_spatial_acc_for_each_internal_pair_within_cutoff(&acc, cutoff, dssp_hbond_energy_cb, &user_param);
}

// ### THE BONDS OF TRACES (CG-DSSP) ###

// The difference b - a by minimum image, as md_hbond takes it: dist2 is the minimum image distance squared of the pair
// as the spatial structure found it (negative when there is none, and max_len bounds the vector instead). Only a vector
// longer than that crosses the periodic boundary and pays for the general minimum image.
static inline vec3_t min_image_diff(vec3_t a, vec3_t b, float dist2, float max_len, const md_unitcell_t* cell) {
    vec3_t d = vec3_sub(b, a);
    const float l2 = vec3_dot(d, d);
    const float ref2 = dist2 >= 0.0f ? dist2 * 1.0001f + 1.0e-4f : max_len * max_len;
    if (l2 > ref2) md_util_min_image_vec3(&d, 1, cell);
    return d;
}

// IUPAC sign: positive in an alpha helix
static inline float cg_dssp_dihedral(vec3_t p0, vec3_t p1, vec3_t p2, vec3_t p3) {
    const vec3_t b0 = vec3_sub(p0, p1);
    const vec3_t b1 = vec3_normalize(vec3_sub(p2, p1));
    const vec3_t b2 = vec3_sub(p3, p2);
    const vec3_t v = vec3_sub(b0, vec3_mul1(b1, vec3_dot(b0, b1)));
    const vec3_t w = vec3_sub(b2, vec3_mul1(b1, vec3_dot(b2, b1)));
    return atan2f(vec3_dot(vec3_cross(b1, v), w), vec3_dot(v, w));
}

static inline float cg_dssp_angle(vec3_t a, vec3_t b, vec3_t c) {
    const vec3_t u = vec3_sub(a, b);
    const vec3_t v = vec3_sub(c, b);
    return acosf(CLAMP(vec3_dot(u, v) / (vec3_length(u) * vec3_length(v)), -1.0f, 1.0f));
}

// The direction of the C=O of residue i of the chain x[beg, end) (beads in one image): the virtual C and O of the residue,
// placed in the frame of its bead by the fit of md_dssp_cg.inl, whose script (tools/fit_cg_dssp.py) defines the
// frame and the features in the order used here. False for a residue without a frame: the ends of the chain, and
// chains of three beads.
static bool cg_dssp_dipole(vec3_t* out_dipole, const vec3_t* x, int beg, int end, int i) {
    if (i < beg + 1 || i + 1 >= end) return false;
    const bool has_a = i - 2 >= beg;
    const bool has_b = i + 2 < end;
    if (!has_a && !has_b) return false;

    const vec3_t t = vec3_normalize(vec3_sub(x[i + 1], x[i - 1]));
    vec3_t m = vec3_sub(vec3_add(x[i - 1], x[i + 1]), vec3_mul1(x[i], 2.0f));
    m = vec3_sub(m, vec3_mul1(t, vec3_dot(m, t)));
    const vec3_t n = vec3_normalize(m);
    const vec3_t b = vec3_cross(t, n);

    const float tau = cg_dssp_angle(x[i - 1], x[i], x[i + 1]);
    const float ta  = has_a ? cg_dssp_dihedral(x[i - 2], x[i - 1], x[i], x[i + 1]) : 0.0f;
    const float tb  = has_b ? cg_dssp_dihedral(x[i - 1], x[i], x[i + 1], x[i + 2]) : 0.0f;
    const float tp  = has_a ? cg_dssp_angle(x[i - 2], x[i - 1], x[i]) : tau;
    const float tn  = has_b ? cg_dssp_angle(x[i], x[i + 1], x[i + 2]) : tau;

    const float ct = cosf(tau), st = sinf(tau);
    const float ca = cosf(ta), sa = sinf(ta);
    const float cb = cosf(tb), sb = sinf(tb);
    const float cp = cosf(tp), cn = cosf(tn);
    const float da[4] = { ca, sa, ca * ca - sa * sa, 2.0f * sa * ca };    // cos, sin of theta_a and of 2 theta_a
    const float db[4] = { cb, sb, cb * cb - sb * sb, 2.0f * sb * cb };

    float f[34];
    int k = 0;
    f[k++] = 1.0f; f[k++] = ct; f[k++] = st; f[k++] = ct * ct;
    f[k++] = cp; f[k++] = cn; f[k++] = cp * cb; f[k++] = cn * cb; f[k++] = cp * ca; f[k++] = cn * ca;
    f[k++] = cp * sb; f[k++] = cn * sb; f[k++] = cp * sa; f[k++] = cn * sa;
    const float (*w)[6];
    if (has_a && has_b) {
        for (int j = 0; j < 4; ++j) f[k++] = da[j];
        for (int j = 0; j < 4; ++j) f[k++] = db[j];
        f[k++] = ca * cb; f[k++] = sa * sb; f[k++] = sa * cb; f[k++] = ca * sb;
        for (int j = 0; j < 4; ++j) f[k++] = ct * da[j];
        for (int j = 0; j < 4; ++j) f[k++] = ct * db[j];
        w = cg_dssp_carbonyl_w1;
        ASSERT(k == (int)ARRAY_SIZE(cg_dssp_carbonyl_w1));
    } else if (has_b) {
        for (int j = 0; j < 4; ++j) f[k++] = db[j];
        for (int j = 0; j < 4; ++j) f[k++] = ct * db[j];
        w = cg_dssp_carbonyl_w2;
        ASSERT(k == (int)ARRAY_SIZE(cg_dssp_carbonyl_w2));
    } else {
        for (int j = 0; j < 4; ++j) f[k++] = da[j];
        for (int j = 0; j < 4; ++j) f[k++] = ct * da[j];
        w = cg_dssp_carbonyl_w3;
        ASSERT(k == (int)ARRAY_SIZE(cg_dssp_carbonyl_w3));
    }

    // {C, O} in (t, n, b), and O - C in the world
    float co[6] = {0};
    for (int j = 0; j < k; ++j) {
        for (int c = 0; c < 6; ++c) co[c] += f[j] * w[j][c];
    }
    const vec3_t d = vec3_add(vec3_add(vec3_mul1(t, co[3] - co[0]), vec3_mul1(n, co[4] - co[1])), vec3_mul1(b, co[5] - co[2]));
    const float len = vec3_length(d);
    if (!(len > 0.0f)) return false;
    *out_dipole = vec3_mul1(d, 1.0f / len);
    return true;
}

typedef struct cg_dssp_param_t {
    const vec3_t* site;             // [peptide]
    const vec3_t* dipole;           // [peptide]
    const uint32_t* peptide_seg;    // [peptide] -> segment, whose C=O it holds
    const uint32_t* res_range_id;   // [segment]
    const md_unitcell_t* cell;
    dssp_res_hbonds_t* res_hbonds;
    md_array(residue_pair_t)* sheet_candidates;
    md_allocator_i* alloc;
} cg_dssp_param_t;

// The bond of the C=O of acc to the N-H of don, and the residue pairs whose bridge pattern it can begin: hb(i-1, j),
// hb(j-1, i), hb(i, j) and hb(i-1, j+1) (see dssp_assign)
static inline void cg_dssp_bond(cg_dssp_param_t* p, uint32_t acc, uint32_t don, float energy) {
    dssp_hbond_insert(p->res_hbonds, don, acc, energy);
    const int64_t a = acc, d = don;
    const int64_t cand[4][2] = { {a + 1, d}, {d, a + 1}, {a, d}, {a + 1, d - 1} };
    for (int c = 0; c < 4; ++c) {
        const int64_t i = MIN(cand[c][0], cand[c][1]);
        const int64_t j = MAX(cand[c][0], cand[c][1]);
        // In range: a bond is never to the last segment's successor, and acc + 1 is a segment as acc holds a peptide
        if (i < 0 || (p->res_range_id[i] == p->res_range_id[j] && j - i <= 2)) continue;
        residue_pair_t pair = { (uint32_t)i, (uint32_t)j };
        md_array_push(*p->sheet_candidates, pair, p->alloc);
    }
}

static void cg_dssp_pair_cb(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    cg_dssp_param_t* p = (cg_dssp_param_t*)user_param;
    for (size_t idx = 0; idx < num_pairs; ++idx) {
        const uint32_t a = i_idx[idx];
        const uint32_t b = j_idx[idx];
        const uint32_t k = p->peptide_seg[a];
        const uint32_t l = p->peptide_seg[b];
        const bool same_chain = p->res_range_id[k] == p->res_range_id[l];
        const uint32_t sep = k > l ? k - l : l - k;
        // Neighbors share a bead
        if (same_chain && sep < 2) continue;

        const bool turn = same_chain && sep <= 4;
        const float max_dist = turn ? CG_DSSP_TURN_DISTANCE  : CG_DSSP_SHEET_DISTANCE;
        const float min_cos  = turn ? CG_DSSP_TURN_ALIGNMENT : CG_DSSP_SHEET_ALIGNMENT;
        if (ij_dist2[idx] >= max_dist * max_dist) continue;

        // The distance is the spatial structure's (minimum image), the direction the vector of that image
        const float d = sqrtf(ij_dist2[idx]);
        if (!(d > 0.0f)) continue;
        const vec3_t r = vec3_mul1(min_image_diff(p->site[a], p->site[b], ij_dist2[idx], 0.0f, p->cell), 1.0f / d);

        const float cos_a = vec3_dot(p->dipole[a], r);
        const float cos_b = vec3_dot(p->dipole[b], r);
        if (cos_a > min_cos && cos_b > min_cos) {
            // Both point from a to b: the C=O of k to the N-H of l+1
            cg_dssp_bond(p, k, l + 1, -1.0f - 0.5f * (cos_a + cos_b));
        } else if (cos_a < -min_cos && cos_b < -min_cos) {
            cg_dssp_bond(p, l, k + 1, -1.0f + 0.5f * (cos_a + cos_b));
        }
    }
}

// The bonds and the bends of the traces
static void dssp_trace_bonds(dssp_res_hbonds_t res_hbonds[], md_array(residue_pair_t)* sheet_candidates, uint8_t res_bend[], const vec3_t* xyz, const md_unitcell_t* cell,
    const md_protein_backbone_data_t* backbone, const uint32_t* res_range_id, md_allocator_i* temp_alloc)
{
    const size_t num_segments = backbone->segment.count;
    vec3_t* x = md_alloc(temp_alloc, sizeof(vec3_t) * num_segments);
    vec3_t* site = md_alloc(temp_alloc, sizeof(vec3_t) * num_segments);
    vec3_t* dipole = md_alloc(temp_alloc, sizeof(vec3_t) * num_segments);
    uint32_t* peptide_seg = md_alloc(temp_alloc, sizeof(uint32_t) * num_segments);
    size_t num_peptides = 0;

    for (size_t range_idx = 0; range_idx < backbone->range.count; ++range_idx) {
        if (!md_protein_backbone_range_is_trace(backbone, range_idx)) continue;
        const int beg = (int)backbone->range.offset[range_idx];
        const int end = (int)backbone->range.offset[range_idx + 1];
        // The chain in one image: the beads are bonded, and so each the minimum image of the one before it
        x[beg] = xyz[backbone->segment.atoms[beg].ca];
        for (int i = beg + 1; i < end; ++i) {
            x[i] = vec3_add(x[i - 1], min_image_diff(xyz[backbone->segment.atoms[i - 1].ca], xyz[backbone->segment.atoms[i].ca], -1.0f, CG_DSSP_MAX_BOND, cell));
        }
        for (int i = beg + 2; i + 2 < end; ++i) {
            res_bend[i] = dssp_bend(vec4_from_vec3(vec3_sub(x[i - 2], x[i]), 0), vec4_from_vec3(vec3_sub(x[i + 2], x[i]), 0));
        }
        for (int k = beg; k + 1 < end; ++k) {
            vec3_t mu;
            if (!cg_dssp_dipole(&mu, x, beg, end, k)) continue;
            site[num_peptides] = vec3_mul1(vec3_add(x[k], x[k + 1]), 0.5f);
            dipole[num_peptides] = mu;
            peptide_seg[num_peptides] = (uint32_t)k;
            num_peptides += 1;
        }
    }
    if (num_peptides < 2) return;

    md_coord_stream_t stream = md_coord_stream_from_aos((const float*)site, sizeof(vec3_t), NULL, num_peptides);
    md_spatial_acc_t acc = { .alloc = temp_alloc };
    const double cutoff = MAX(CG_DSSP_TURN_DISTANCE, CG_DSSP_SHEET_DISTANCE);
    md_spatial_acc_init(&acc, &(md_spatial_acc_desc_t){ .coords = &stream, .cutoff = cutoff, .unitcell = cell });

    cg_dssp_param_t param = {
        .site = site,
        .dipole = dipole,
        .peptide_seg = peptide_seg,
        .res_range_id = res_range_id,
        .cell = cell,
        .res_hbonds = res_hbonds,
        .sheet_candidates = sheet_candidates,
        .alloc = temp_alloc,
    };
    md_spatial_acc_for_each_internal_pair_within_cutoff(&acc, cutoff, cg_dssp_pair_cb, &param);
}

static int residue_pair_cmp(const void* a, const void* b) {
    const residue_pair_t* x = (const residue_pair_t*)a;
    const residue_pair_t* y = (const residue_pair_t*)b;
    if (x->i != y->i) return x->i < y->i ? -1 : 1;
    if (x->j != y->j) return x->j < y->j ? -1 : 1;
    return 0;
}

// ### THE PATTERNS ###

// Turns, helices and bends per range, bridges and ladders between any two residues, from the bonds, the bends (dssp_bend,
// at the middle residue) and the candidate pairs for bridges (i < j, at least 3 apart within a range)
static void dssp_assign(md_secondary_structure_t out_secondary_structure[], const md_protein_backbone_data_t* backbone, const dssp_res_hbonds_t res_hbonds[],
    const uint32_t* res_range_id, const uint8_t* res_bend, const residue_pair_t* sheet_candidates, size_t num_sheet_candidates, md_temp_scope_t temp)
{
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    const size_t backbone_segment_count = backbone->segment.count;
    const size_t backbone_range_count   = backbone->range.count;
    const uint32_t* backbone_range_offsets = backbone->range.offset;

    uint8_t* turn_mask = md_temp_alloc_zero_array(temp, uint8_t, backbone_segment_count);
    uint32_t* ss_flags = md_temp_alloc_zero_array(temp, uint32_t, backbone_segment_count);

    // ---- Per-range passes: helices, turns, bends (set flags only) ----
    // Ground-truth DSSP helix assignment:
    // Two consecutive n-turns define a helix, and then residues i+1..i+n get helix flag.
    for (size_t range_idx = 0; range_idx < backbone_range_count; ++range_idx) {
        const md_urange_t range = {backbone_range_offsets[range_idx], backbone_range_offsets[range_idx + 1]};

        // Bits in turn_mask (local indices): 1 = 3-turn, 2 = 4-turn, 4 = 5-turn
        for (uint32_t i = range.beg; i < range.end; ++i) {
            if (i + 3 < range.end && dssp_test_bond(res_hbonds, i, i + 3)) {
                turn_mask[i] |= 1;
            }
            if (i + 4 < range.end && dssp_test_bond(res_hbonds, i, i + 4)) {
                turn_mask[i] |= 2;
            }
            if (i + 5 < range.end && dssp_test_bond(res_hbonds, i, i + 5)) {
                turn_mask[i] |= 4;
            }
        }

        // DSSP helix assignment from consecutive turns:
        // If turn_n(i) && turn_n(i+1) then mark residues (i+1 ... i+n) as helix_n.
        // n=3 => 3_10 helix, n=4 => alpha helix, n=5 => pi helix.
        for (uint32_t i = range.beg; i + 1 < range.end; ++i) {
            const uint8_t a = turn_mask[i];
            const uint8_t b = turn_mask[i + 1];

            // 3_10 helix (G)
            if ((a & 1) && (b & 1)) {
                const uint32_t beg = i + 1;
                const uint32_t end = MIN(i + 3, range.end - 1);
                for (uint32_t r = beg; r <= end; ++r) {
                    ss_flags[r] |= SS_FLAG_HELIX_310;
                }
            }

            // alpha helix (H)
            if ((a & 2) && (b & 2)) {
                const uint32_t beg = i + 1;
                const uint32_t end = MIN(i + 4, range.end - 1);
                for (uint32_t r = beg; r <= end; ++r) {
                    ss_flags[r] |= SS_FLAG_HELIX_ALPHA;
                }
            }

            // pi helix (I)
            if ((a & 4) && (b & 4)) {
                const uint32_t beg = i + 1;
                const uint32_t end = MIN(i + 5, range.end - 1);
                for (uint32_t r = beg; r <= end; ++r) {
                    ss_flags[r] |= SS_FLAG_HELIX_PI;
                }
            }

            // Just turn
            if (a) {
				const int len = (a & 4) ? 5 : (a & 2) ? 4 : 3;
                const uint32_t beg = i + 1;
                const uint32_t end = MIN(i + len, range.end - 1);
                for (uint32_t r = beg; r <= end; ++r) {
                    ss_flags[r] |= SS_FLAG_TURN;
                }
			}
        }

        // Bends
        for (uint32_t i = range.beg + 2; i + 2 < range.end; ++i) {
            if (res_bend[i]) {
                for (int k = -2; k <= 2; ++k) {
                    ss_flags[i + k] |= SS_FLAG_BEND;
				}
            }
        }
    }

    // ---- Global pass: sheets / bridges ----

    typedef struct {
        uint32_t i;
        uint32_t j;
		uint64_t chain_pair_id_and_type; // Unique id for the chain_pair from (res_range_id[i], res_range_id[j]) and the bridge type (parallel/antiparallel)
    } dssp_bridge_t;

    // At most one bridge per candidate
    dssp_bridge_t* bridges = md_temp_alloc_array(temp, dssp_bridge_t, MAX(num_sheet_candidates, 1));
    size_t num_bridges = 0;

    // Helper lambdas are not allowed in C; keep as local inline-ish macros
#define SAME_RANGE(a,b) (res_range_id[(a)] == res_range_id[(b)])

	md_hashmap32_t parallel_bridge_map     = { .allocator = temp_alloc };
	md_hashmap32_t antiparallel_bridge_map = { .allocator = temp_alloc };

    for (size_t c = 0; c < num_sheet_candidates; ++c) {
        uint32_t i = sheet_candidates[c].i;
        uint32_t j = sheet_candidates[c].j;

        // Neighbor existence within each chain/range
        const bool i_has_im1 = (i > 0) && SAME_RANGE(i - 1, i);
        const bool i_has_ip1 = (i + 1 < (uint32_t)backbone_segment_count) && SAME_RANGE(i + 1, i);
        const bool j_has_jm1 = (j > 0) && SAME_RANGE(j - 1, j);
        const bool j_has_jp1 = (j + 1 < (uint32_t)backbone_segment_count) && SAME_RANGE(j + 1, j);

        // Evaluate DSSP bridge patterns (directed H-bonds)
        bool antiparallel = false;
        bool parallel     = false;

        // Antiparallel:
        // (i -> j && j -> i) OR (i-1 -> j+1 && j-1 -> i+1)
        if (dssp_test_bond(res_hbonds, i, j) && dssp_test_bond(res_hbonds, j, i)) {
            antiparallel = true;
        } else if (i_has_im1 && i_has_ip1 && j_has_jm1 && j_has_jp1) {
            if (dssp_test_bond(res_hbonds, i - 1, j + 1) && dssp_test_bond(res_hbonds, j - 1, i + 1)) {
                antiparallel = true;
            }
        }

        // Parallel:
        // (i-1 -> j && j -> i+1) OR (j-1 -> i && i -> j+1)
        if (!antiparallel) {
            if (i_has_im1 && i_has_ip1) {
                if (dssp_test_bond(res_hbonds, i - 1, j) && dssp_test_bond(res_hbonds, j, i + 1)) {
                    parallel = true;
                }
            }
            if (!parallel && j_has_jm1 && j_has_jp1) {
                if (dssp_test_bond(res_hbonds, j - 1, i) && dssp_test_bond(res_hbonds, i, j + 1)) {
                    parallel = true;
                }
            }
        }

        if (!(antiparallel || parallel)) {
            continue;
        }

		uint32_t ri = res_range_id[i];
		uint32_t rj = res_range_id[j];

		uint64_t res_range_packed = (uint64_t)MIN(ri, rj) << 32 | MAX(ri, rj);

        uint32_t idx = (uint32_t)(num_bridges++);

        // Record bridge (normalize ordering i<j already true)
        bridges[idx] = (dssp_bridge_t){
            .i = i,
            .j = j,
            .chain_pair_id_and_type = (uint64_t)(parallel ? 1 : 0) << 63 | res_range_packed,
        };

        // DSSP marks both residues as "bridge" even before ladder/sheet conversion
        ss_flags[i] |= SS_FLAG_BRIDGE;
        ss_flags[j] |= SS_FLAG_BRIDGE;

		uint64_t key = md_hash64(&bridges[idx], sizeof(dssp_bridge_t), 0);
		if (parallel) {
			md_hashmap_add(&parallel_bridge_map, key, idx);
		} else {
			md_hashmap_add(&antiparallel_bridge_map, key, idx);
		}
    }
#undef SAME_RANGE

	// ---- Build ladders by directly looking up bridges in hashmaps ----
    // We mark as SHEET if a ladder has >= 2 bridges (i.e. spans >= 3 residues on at least one strand).
    bool* used = md_temp_alloc_zero_array(temp, bool, MAX(num_bridges, 1));
    for (size_t b = 0; b < num_bridges; ++b) {
        if (used[b]) continue;
        const dssp_bridge_t seed = bridges[b];
        used[b] = true;
        // Ladder endpoints
        uint32_t i_min = seed.i, i_max = seed.i;
        uint32_t j_min = seed.j, j_max = seed.j;
        // Chain pair and type
        const uint64_t chain_pair_id_and_type = seed.chain_pair_id_and_type;
        const bool is_parallel = (chain_pair_id_and_type >> 63) != 0;
        md_hashmap32_t* bridge_map = is_parallel ? &parallel_bridge_map : &antiparallel_bridge_map;

        bool grew = true;
        while (grew) {
            grew = false;
            // Look for the next bridge in the ladder by checking the expected next positions based on the ladder type
            uint32_t next_i_min, next_i_max, next_j_min, next_j_max;
            if (is_parallel) {
                // Parallel ladder: (i+1, j+1)
                next_i_min = i_max + 1;
                next_j_min = j_max + 1;
                next_i_max = i_min - 1;
                next_j_max = j_min - 1;
            } else {
                // Antiparallel ladder: (i+1, j-1)
                next_i_min = i_max + 1;
                next_j_min = j_min - 1;
                next_i_max = i_min - 1;
                next_j_max = j_max + 1;
            }
            // Check both possible extensions of the ladder
            dssp_bridge_t query_bridge_1 = { .i = next_i_min, .j = next_j_min, .chain_pair_id_and_type = chain_pair_id_and_type };
            dssp_bridge_t query_bridge_2 = { .i = next_i_max, .j = next_j_max, .chain_pair_id_and_type = chain_pair_id_and_type };
			uint64_t key_1 = md_hash64(&query_bridge_1, sizeof(dssp_bridge_t), 0);
			uint64_t key_2 = md_hash64(&query_bridge_2, sizeof(dssp_bridge_t), 0);

            uint32_t* found_idx_1 = md_hashmap_get(bridge_map, key_1);
            uint32_t* found_idx_2 = md_hashmap_get(bridge_map, key_2);
            if (found_idx_1) {
                const dssp_bridge_t next_bridge = bridges[*found_idx_1];
                used[*found_idx_1] = true;
                i_min = MIN(i_min, next_bridge.i);
                i_max = MAX(i_max, next_bridge.i);
                j_min = MIN(j_min, next_bridge.j);
                j_max = MAX(j_max, next_bridge.j);
                grew = true;
            } else if (found_idx_2) {
                const dssp_bridge_t next_bridge = bridges[*found_idx_2];
                used[*found_idx_2] = true;
                i_min = MIN(i_min, next_bridge.i);
                i_max = MAX(i_max, next_bridge.i);
                j_min = MIN(j_min, next_bridge.j);
                j_max = MAX(j_max, next_bridge.j);
                grew = true;
            }
        }
        const uint32_t len_i = i_max - i_min + 1;
        const uint32_t len_j = j_max - j_min + 1;
        if (len_i >= 3 || len_j >= 3) {
            for (uint32_t k = i_min; k <= i_max; ++k) ss_flags[k] |= SS_FLAG_SHEET;
            for (uint32_t k = j_min; k <= j_max; ++k) ss_flags[k] |= SS_FLAG_SHEET;
        }
	}

    // ---- Final resolution: convert flags to enum with priority ----
    // Priority: HELIX (PI > 310 > ALPHA) > SHEET > BRIDGE > TURN > BEND > COIL
    MEMSET(out_secondary_structure, 0, sizeof(md_secondary_structure_t) * backbone_segment_count);

    for (size_t i = 0; i < backbone_segment_count; ++i) {
        uint32_t f = ss_flags[i];

        if (f & SS_FLAG_HELIX_PI) {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_HELIX_PI;
        } else if (f & SS_FLAG_HELIX_ALPHA) {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_HELIX_ALPHA;
        } else if (f & SS_FLAG_HELIX_310) {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_HELIX_310;
        } else if (f & SS_FLAG_SHEET) {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_BETA_SHEET;
        } else if (f & SS_FLAG_BRIDGE) {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_BETA_BRIDGE;
        } else if (f & SS_FLAG_TURN) {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_TURN;
        } else if (f & SS_FLAG_BEND) {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_BEND;
        } else {
            out_secondary_structure[i] = MD_SECONDARY_STRUCTURE_COIL;
        }
    }
}

static void dssp(md_secondary_structure_t out_secondary_structure[], size_t capacity, const vec3_t* xyz, const md_unitcell_t* cell, const md_protein_backbone_data_t* backbone) {
    ASSERT(out_secondary_structure);
    ASSERT(xyz);
    ASSERT(cell);
    ASSERT(backbone);

    const size_t backbone_segment_count = backbone->segment.count;
    const size_t backbone_range_count   = backbone->range.count;
    const uint32_t* backbone_range_offsets = backbone->range.offset;

    if (capacity < backbone_segment_count) {
        MD_LOG_ERROR("Capacity of out_secondary_structure is insufficient to hold the results");
        return;
    }

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    dssp_res_hbonds_t* res_hbonds = md_temp_alloc_array(temp, dssp_res_hbonds_t, backbone_segment_count);
    uint32_t* res_range_id = md_temp_alloc_array(temp, uint32_t, backbone_segment_count);
    uint8_t* res_bend = md_temp_alloc_zero_array(temp, uint8_t, backbone_segment_count);

    md_array(residue_pair_t) sheet_candidates = NULL;
	md_array_ensure(sheet_candidates, 1024, temp_alloc);

    size_t num_trace = 0;
    for (uint32_t range_idx = 0; range_idx < backbone_range_count; ++range_idx) {
        const md_urange_t range = { backbone_range_offsets[range_idx], backbone_range_offsets[range_idx + 1] };
        for (uint32_t i = range.beg; i < range.end; ++i) {
            res_range_id[i] = range_idx;
        }
        num_trace += md_protein_backbone_range_is_trace(backbone, range_idx) ? range.end - range.beg : 0;
    }

	// Safe to set to zero, as we will only be comparing energies (which are negative)
	MEMSET(res_hbonds, 0, sizeof(dssp_res_hbonds_t) * backbone_segment_count);

    if (num_trace < backbone_segment_count) {
        dssp_atomistic_bonds(res_hbonds, &sheet_candidates, res_bend, xyz, cell, backbone, num_trace, temp_alloc);
    }
    if (num_trace > 0) {
        dssp_trace_bonds(res_hbonds, &sheet_candidates, res_bend, xyz, cell, backbone, res_range_id, temp_alloc);
        // A bond of a trace gives each pair it can begin a bridge in, so the same pair comes more than once
        const size_t n = md_array_size(sheet_candidates);
        if (n > 1) {
            qsort(sheet_candidates, n, sizeof(residue_pair_t), residue_pair_cmp);
            size_t u = 1;
            for (size_t c = 1; c < n; ++c) {
                if (residue_pair_cmp(&sheet_candidates[u - 1], &sheet_candidates[c]) != 0) sheet_candidates[u++] = sheet_candidates[c];
            }
            md_array_shrink(sheet_candidates, u);
        }
    }

    dssp_assign(out_secondary_structure, backbone, res_hbonds, res_range_id, res_bend, sheet_candidates, md_array_size(sheet_candidates), temp);

    md_temp_end(temp);
}

bool md_dssp_compute(md_secondary_structure_t secondary_structure[], size_t capacity, const vec3_t* xyz, const md_unitcell_t* cell, const md_protein_backbone_data_t* backbone) {

    if (!backbone) return false;
    if (capacity < backbone->segment.count) {
        MD_LOG_ERROR("Secondary structure capacity is not sufficient to hold the number of segments in the protein backbone");
        return false;
    }

    if (secondary_structure && xyz && backbone) {
        dssp(secondary_structure, capacity, xyz, cell, backbone);
        return true;
    }
    return false;
}
