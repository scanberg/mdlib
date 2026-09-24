#pragma once

// Brute force reference for md_contact, shared by the kernel tests (test_contact.c) and the script tests
// (test_script.c), so both are held to the same oracle. Include after utest.h and the mdlib headers.

#include <md_contact.h>
#include <md_system.h>
#include <md_util.h>
#include <core/md_allocator.h>
#include <core/md_array.h>
#include <core/md_bitfield.h>
#include <core/md_vec_math.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

// Reference: every group pair, every atom pair of it, straight from the definition.

// Distances within REF_EPS of the limit are counted as ambiguous: the result and the reference compute them in a
// different order in single precision, and coordinates on a grid (gro stores 0.01 Å) put some exactly on a round cutoff.
#define REF_EPS 1.0e-4f

typedef struct ref_pair_t {
    uint32_t i, j;
    uint32_t atom_pairs_lo;     // Atom pairs certainly within the limit
    uint32_t atom_pairs_hi;     // ... and those which may be, depending on rounding
    float d_min;
} ref_pair_t;

static inline float ref_distance(const md_system_state_t* state, uint32_t a, uint32_t b) {
    vec3_t dx = { state->x[b] - state->x[a], state->y[b] - state->y[a], state->z[b] - state->z[a] };
    md_util_min_image_vec3(&dx, 1, &state->unitcell);
    return sqrtf(dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
}

// Atoms within max_bonds bonds of a, using md_util_mask_grow_by_bonds
static inline bool ref_excluded(md_bitfield_t* cache, const md_system_t* sys, uint32_t max_bonds, uint32_t a, uint32_t b, md_allocator_i* alloc) {
    if (max_bonds == 0) return false;
    if (!cache[a].alloc) {
        cache[a] = md_bitfield_create(alloc);
        md_bitfield_set_bit(&cache[a], a);
        md_util_mask_grow_by_bonds(&cache[a], sys, max_bonds, NULL);
    }
    return md_bitfield_test_bit(&cache[a], b);
}

static inline md_array(ref_pair_t) ref_contacts(const md_contact_desc_t* desc, const md_system_t* sys, const md_system_state_t* state, md_allocator_i* alloc) {
    const bool self = desc->group_b == NULL;
    const md_bitfield_t* A = desc->group_a;
    const md_bitfield_t* B = self ? desc->group_a : desc->group_b;
    const size_t num_b = self ? desc->num_a : desc->num_b;

    float* radius = NULL;
    if (desc->criterion == MD_CONTACT_CRITERION_RADII) {
        radius = md_alloc(alloc, sizeof(float) * sys->atom.count);
        md_atom_extract_radii(radius, 0, sys->atom.count, &sys->atom);
    }
    md_bitfield_t* excl = md_alloc(alloc, sizeof(md_bitfield_t) * sys->atom.count);
    memset(excl, 0, sizeof(md_bitfield_t) * sys->atom.count);

    md_array(ref_pair_t) result = 0;
    for (uint32_t i = 0; i < desc->num_a; ++i) {
        for (uint32_t j = self ? i + 1 : 0; j < num_b; ++j) {
            if (self && desc->min_separation > 1 && j - i < desc->min_separation && (!desc->group_parent || desc->group_parent[i] == desc->group_parent[j])) continue;
            ref_pair_t p = { i, j, 0, 0, FLT_MAX };
            md_bitfield_iter_t ia = md_bitfield_iter_create(&A[i]);
            while (md_bitfield_iter_next(&ia)) {
                const uint32_t a = (uint32_t)md_bitfield_iter_idx(&ia);
                md_bitfield_iter_t jb = md_bitfield_iter_create(&B[j]);
                while (md_bitfield_iter_next(&jb)) {
                    const uint32_t b = (uint32_t)md_bitfield_iter_idx(&jb);
                    if (a == b) continue;
                    const float d = ref_distance(state, a, b);
                    if (desc->particle_label && desc->particle_label[a] == desc->particle_label[b]) continue;
                    float limit = radius ? radius[a] + radius[b] + (float)desc->cutoff : (float)desc->cutoff;
                    if (desc->criterion == MD_CONTACT_CRITERION_TYPE_PAIR) {
                        limit = desc->type_cutoff[desc->particle_type[a] * desc->num_types + desc->particle_type[b]];
                    }
                    if (!(d < limit + REF_EPS)) continue;
                    if (ref_excluded(excl, sys, desc->exclude_bonds, a, b, alloc)) continue;
                    p.atom_pairs_hi += 1;
                    if (d < limit - REF_EPS) {
                        p.atom_pairs_lo += 1;
                        p.d_min = MIN(p.d_min, d);
                    }
                }
            }
            if (p.atom_pairs_hi) {
                md_array_push(result, p, alloc);
            }
        }
    }
    return result;
}

static inline bool same_contacts(const md_contact_set_t* set, const ref_pair_t* ref, size_t ref_count) {
    bool ok = true;
    size_t shown = 0;
    size_t a = 0, b = 0;
    for (size_t k = 1; k < set->count; ++k) {
        if (((uint64_t)set->i[k-1] << 32 | set->j[k-1]) >= ((uint64_t)set->i[k] << 32 | set->j[k])) {
            printf("Result not sorted at %zu\n", k);
            return false;
        }
    }
    while (a < set->count || b < ref_count) {
        const uint64_t ka = a < set->count ? ((uint64_t)set->i[a] << 32 | set->j[a]) : UINT64_MAX;
        const uint64_t kb = b < ref_count ? ((uint64_t)ref[b].i << 32 | ref[b].j) : UINT64_MAX;
        const char* err = NULL;
        if (ka == kb) {
            const ref_pair_t* r = &ref[b];
            if (set->atom_pairs[a] < r->atom_pairs_lo || set->atom_pairs[a] > r->atom_pairs_hi) err = "atom pair count";
            else if (r->atom_pairs_lo && fabsf(set->d_min[a] - r->d_min) > 1.0e-3f) err = "d_min";
            if (err && shown++ < 10) printf("  (%u, %u): %s, got %u pairs d_min %f, expected %u..%u pairs d_min %f\n", r->i, r->j, err, set->atom_pairs[a], set->d_min[a], r->atom_pairs_lo, r->atom_pairs_hi, r->d_min);
            ++a; ++b;
        } else if (ka < kb) {
            err = "only in result";
            if (shown++ < 10) printf("  (%u, %u): only in result, d_min %f\n", set->i[a], set->j[a], set->d_min[a]);
            ++a;
        } else {
            if (ref[b].atom_pairs_lo) {
                err = "only in reference";
                if (shown++ < 10) printf("  (%u, %u): only in reference, d_min %f\n", ref[b].i, ref[b].j, ref[b].d_min);
            }
            ++b;
        }
        if (err) ok = false;
    }
    return ok;
}

