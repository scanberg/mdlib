#include <md_contact.h>

#include <md_system.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_bitfield.h>
#include <core/md_coord_stream.h>
#include <core/md_hash.h>
#include <core/md_log.h>
#include <core/md_spatial_acc.h>

#include <float.h>
#include <math.h>
#include <stdlib.h>

// ### PARTICLE PAIRS ###

// The particles of a set, ascending
static bool extract_atoms(int32_t** out_atoms, size_t* out_count, const md_bitfield_t* set, uint32_t num_atoms, md_allocator_i* alloc) {
    if (md_bitfield_end_bit(set) > num_atoms) {
        MD_LOG_ERROR("Contact set refers to particles beyond the system (%u particles)", num_atoms);
        return false;
    }
    const size_t count = md_bitfield_popcount(set);
    int32_t* atoms = md_alloc(alloc, sizeof(int32_t) * MAX(count, 1));
    const size_t written = count ? md_bitfield_iter_extract_indices(atoms, count, md_bitfield_iter_create(set)) : 0;
    ASSERT(written == count);
    (void)written;
    *out_atoms = atoms;
    *out_count = count;
    return true;
}

// For every member particle, the member particles within max_bonds bonds of it (itself excluded), as sorted
// compressed rows over all particles. Particles which are not members never form a pair, so they are never listed.
static void build_exclusions(md_contact_pairs_t* p, const md_system_t* sys, uint32_t max_bonds, const uint8_t* is_member, md_allocator_i* alloc) {
    const uint32_t N = p->num_atoms;
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    uint8_t* depth = md_alloc(temp_alloc, N);    // 0 = not reached, otherwise 1 + number of bonds
    MEMSET(depth, 0, N);
    md_array(uint32_t) visited = 0;
    md_array(uint32_t) queue   = 0;
    md_array(uint32_t) partner = 0;

    uint32_t* off = md_alloc(alloc, sizeof(uint32_t) * (N + 1));
    off[0] = 0;
    const uint32_t depth_limit = MIN(max_bonds, 254) + 1;

    for (uint32_t k = 0; k < N; ++k) {
        if (is_member[k]) {
            // Breadth first from k up to max_bonds bonds
            md_array_shrink(queue, 0);
            md_array_shrink(visited, 0);
            depth[k] = 1;
            md_array_push(visited, k, temp_alloc);
            md_array_push(queue, k, temp_alloc);
            for (size_t head = 0; head < md_array_size(queue); ++head) {
                const uint32_t cur = queue[head];
                if (depth[cur] >= depth_limit) continue;
                md_bond_iter_t it = md_bond_iter(&sys->bond, cur);
                while (md_bond_iter_has_next(&it)) {
                    const uint32_t next = md_bond_iter_atom_index(&it);
                    md_bond_iter_next(&it);
                    if (next < N && depth[next] == 0) {
                        depth[next] = depth[cur] + 1;
                        md_array_push(visited, next, temp_alloc);
                        md_array_push(queue, next, temp_alloc);
                    }
                }
            }
            // Record the member partners, sorted
            const size_t row_beg = md_array_size(partner);
            for (size_t v = 0; v < md_array_size(visited); ++v) {
                const uint32_t a = visited[v];
                if (a != k && is_member[a]) {
                    md_array_push(partner, a, temp_alloc);
                }
                depth[a] = 0;
            }
            const size_t row_end = md_array_size(partner);
            for (size_t x = row_beg + 1; x < row_end; ++x) {
                // Insertion sort, the rows are short
                const uint32_t val = partner[x];
                size_t y = x;
                while (y > row_beg && partner[y - 1] > val) {
                    partner[y] = partner[y - 1];
                    --y;
                }
                partner[y] = val;
            }
        }
        off[k + 1] = (uint32_t)md_array_size(partner);
    }

    const size_t total = md_array_size(partner);
    uint32_t* atoms = md_alloc(alloc, sizeof(uint32_t) * MAX(total, 1));
    if (total) MEMCPY(atoms, partner, sizeof(uint32_t) * total);

    p->excl_off  = off;
    p->excl_atom = atoms;

    md_temp_end(temp);
}

static inline bool is_excluded(const md_contact_pairs_t* p, uint32_t a, uint32_t b) {
    if (!p->excl_off) return false;
    size_t lo = p->excl_off[a];
    size_t hi = p->excl_off[a + 1];
    while (lo < hi) {
        const size_t mid = (lo + hi) / 2;
        const uint32_t v = p->excl_atom[mid];
        if (v == b) return true;
        if (v < b) lo = mid + 1; else hi = mid;
    }
    return false;
}

void md_contact_pairs_free(md_contact_pairs_t* p) {
    if (!p || !p->alloc) return;
    md_allocator_i* alloc = p->alloc;
    const uint32_t N = p->num_atoms;
    if (p->a_atoms)   md_free(alloc, p->a_atoms, sizeof(int32_t) * MAX(p->num_a_atoms, 1));
    if (p->b_atoms)   md_free(alloc, p->b_atoms, sizeof(int32_t) * MAX(p->num_b_atoms, 1));
    if (p->excl_atom) md_free(alloc, p->excl_atom, sizeof(uint32_t) * MAX(p->excl_off[N], 1));
    if (p->excl_off)  md_free(alloc, p->excl_off, sizeof(uint32_t) * (N + 1));
    if (p->label)     md_free(alloc, p->label, sizeof(uint32_t) * MAX(N, 1));
    MEMSET(p, 0, sizeof(md_contact_pairs_t));
}

bool md_contact_pairs_init(md_contact_pairs_t* p, const md_contact_pairs_desc_t* desc, const md_system_t* sys, md_allocator_i* alloc) {
    ASSERT(p);
    MEMSET(p, 0, sizeof(md_contact_pairs_t));
    if (!desc || !desc->set_a || !sys || !alloc) {
        MD_LOG_ERROR("Contact pairs: missing set, system or allocator");
        return false;
    }
    if (sys->atom.count > UINT32_MAX - 1) {
        MD_LOG_ERROR("Contact pairs: too many particles");
        return false;
    }

    const bool self = desc->set_b == NULL;
    p->alloc     = alloc;
    p->num_atoms = (uint32_t)sys->atom.count;
    p->flags     = self ? MD_CONTACT_FLAG_SELF : MD_CONTACT_FLAG_NONE;
    p->radius    = (float)desc->radius;

    if (!extract_atoms(&p->a_atoms, &p->num_a_atoms, desc->set_a, p->num_atoms, alloc)) goto fail;
    if (!self && !extract_atoms(&p->b_atoms, &p->num_b_atoms, desc->set_b, p->num_atoms, alloc)) goto fail;

    if (desc->particle_label && p->num_atoms) {
        p->label = md_alloc(alloc, sizeof(uint32_t) * p->num_atoms);
        MEMCPY(p->label, desc->particle_label, sizeof(uint32_t) * p->num_atoms);
    }

    if (desc->exclude_bonds > 0 && sys->bond.conn.offset && p->num_atoms) {
        md_temp_scope_t temp = md_temp_begin_avoid(alloc);
        uint8_t* is_member = md_alloc(md_temp_allocator(temp), p->num_atoms);
        MEMSET(is_member, 0, p->num_atoms);
        for (size_t k = 0; k < p->num_a_atoms; ++k) is_member[p->a_atoms[k]] = 1;
        for (size_t k = 0; k < p->num_b_atoms; ++k) is_member[p->b_atoms[k]] = 1;
        build_exclusions(p, sys, desc->exclude_bonds, is_member, alloc);
        md_temp_end(temp);
    }
    return true;
fail:
    md_contact_pairs_free(p);
    return false;
}

typedef struct pair_stream_t {
    const md_contact_pairs_t* p;
    md_contact_pair_callback_t callback;
    void* user_param;
    md_array(uint32_t) a;
    md_array(uint32_t) b;
    md_array(float)    r;
    md_allocator_i* alloc;
} pair_stream_t;

static void pair_stream_callback(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    pair_stream_t* s = (pair_stream_t*)user_param;
    const md_contact_pairs_t* p = s->p;
    const bool self = p->flags & MD_CONTACT_FLAG_SELF;
    const float radius2 = p->radius * p->radius;

    md_array_ensure(s->a, num_pairs, s->alloc);
    md_array_ensure(s->b, num_pairs, s->alloc);
    md_array_ensure(s->r, num_pairs, s->alloc);
    size_t n = 0;
    for (size_t k = 0; k < num_pairs; ++k) {
        uint32_t x = i_idx[k];
        uint32_t y = j_idx[k];
        if (x == y) continue;
        if (!(ij_dist2[k] < radius2)) continue;
        if (p->label && p->label[x] == p->label[y]) continue;
        if (is_excluded(p, x, y)) continue;
        if (self && x > y) { const uint32_t t = x; x = y; y = t; }
        s->a[n] = x;
        s->b[n] = y;
        s->r[n] = sqrtf(ij_dist2[k]);
        ++n;
    }
    if (n) {
        s->callback(s->a, s->b, s->r, n, s->user_param);
    }
}

bool md_contact_pairs_for_each(const md_contact_pairs_t* p, const md_system_state_t* state, md_contact_pair_callback_t callback, void* user_param) {
    if (!p || !p->alloc || !state || !callback) {
        MD_LOG_ERROR("Contact pairs: missing pairs, state or callback");
        return false;
    }
    if (state->num_atoms != p->num_atoms || (p->num_atoms && !state->xyz)) {
        MD_LOG_ERROR("Contact pairs: the state does not match the system the pairs were prepared for");
        return false;
    }

    const bool self = p->flags & MD_CONTACT_FLAG_SELF;
    const size_t num_search = self ? p->num_a_atoms : p->num_b_atoms;
    if (p->radius <= 0.0f || p->num_a_atoms == 0 || num_search == 0) {
        return true;
    }

    // The scratch memory of the stream lives in an arena of its own. A temp scope could share its arena with
    // whatever the callback grows while the stream runs (a temp arena of the caller's), and ending the scope
    // would then rewind the arena beneath the caller's data.
    md_allocator_i* temp_alloc = md_vm_arena_create(GIGABYTES(4));

    pair_stream_t stream = {
        .p = p,
        .callback = callback,
        .user_param = user_param,
        .alloc = temp_alloc,
    };

    const md_spatial_acc_flags_t flags = MD_SPATIAL_ACC_FLAG_USE_COORD_STREAM_IDX;
    md_spatial_acc_t sa = { .alloc = temp_alloc };
    if (self) {
        md_coord_stream_t coords = md_coord_stream_from_aos((const float*)state->xyz, sizeof(vec3_t), p->a_atoms, p->num_a_atoms);
        md_spatial_acc_init(&sa, &coords, p->radius, &state->unitcell, flags);
        md_spatial_acc_for_each_internal_pair_within_cutoff(&sa, p->radius, pair_stream_callback, &stream);
    } else {
        md_coord_stream_t internal = md_coord_stream_from_aos((const float*)state->xyz, sizeof(vec3_t), p->b_atoms, p->num_b_atoms);
        md_coord_stream_t external = md_coord_stream_from_aos((const float*)state->xyz, sizeof(vec3_t), p->a_atoms, p->num_a_atoms);
        md_spatial_acc_init(&sa, &internal, p->radius, &state->unitcell, flags);
        md_spatial_acc_for_each_external_vs_internal_pair_within_cutoff(&sa, &external, p->radius, pair_stream_callback, &stream, flags);
    }

    md_vm_arena_destroy(temp_alloc);
    return true;
}

// ### GROUP CONTACTS ###

// Particle to group mapping in compressed rows over all particles: the groups of particle k are grp[off[k]] .. grp[off[k+1]-1].
// Also the union of the groups, as a bitfield allocated from temp_alloc.
static bool build_membership(uint32_t** out_off, uint32_t** out_grp, md_bitfield_t* out_union,
                             const md_bitfield_t* groups, size_t num_groups, uint32_t num_atoms, md_allocator_i* alloc, md_allocator_i* temp_alloc) {
    uint32_t* off = md_alloc(alloc, sizeof(uint32_t) * (num_atoms + 1));
    MEMSET(off, 0, sizeof(uint32_t) * (num_atoms + 1));
    *out_union = md_bitfield_create(temp_alloc);

    for (size_t g = 0; g < num_groups; ++g) {
        if (md_bitfield_end_bit(&groups[g]) > num_atoms) {
            MD_LOG_ERROR("Contact group %zu refers to particles beyond the system (%u particles)", g, num_atoms);
            md_free(alloc, off, sizeof(uint32_t) * (num_atoms + 1));
            return false;
        }
        md_bitfield_or_inplace(out_union, &groups[g]);
        md_bitfield_iter_t it = md_bitfield_iter_create(&groups[g]);
        while (md_bitfield_iter_next(&it)) {
            off[md_bitfield_iter_idx(&it) + 1] += 1;
        }
    }
    for (uint32_t k = 0; k < num_atoms; ++k) {
        off[k + 1] += off[k];
    }

    const size_t num_entries = off[num_atoms];
    uint32_t* grp = md_alloc(alloc, sizeof(uint32_t) * MAX(num_entries, 1));

    // Fill in group order, which leaves every row sorted
    uint32_t* cursor = md_alloc(temp_alloc, sizeof(uint32_t) * MAX(num_atoms, 1));
    if (num_atoms) MEMCPY(cursor, off, sizeof(uint32_t) * num_atoms);
    for (size_t g = 0; g < num_groups; ++g) {
        md_bitfield_iter_t it = md_bitfield_iter_create(&groups[g]);
        while (md_bitfield_iter_next(&it)) {
            grp[cursor[md_bitfield_iter_idx(&it)]++] = (uint32_t)g;
        }
    }

    *out_off = off;
    *out_grp = grp;
    return true;
}

static size_t membership_bytes(uint32_t num_atoms, const uint32_t* off) {
    return off ? sizeof(uint32_t) * MAX(off[num_atoms], 1) : 0;
}

void md_contact_query_free(md_contact_query_t* q) {
    if (!q || !q->alloc) return;
    md_allocator_i* alloc = q->alloc;
    const uint32_t N = q->num_atoms;
    md_contact_pairs_free(&q->pairs);
    if (q->a_grp)  md_free(alloc, q->a_grp, membership_bytes(N, q->a_off));
    if (q->a_off)  md_free(alloc, q->a_off, sizeof(uint32_t) * (N + 1));
    if (q->b_grp)  md_free(alloc, q->b_grp, membership_bytes(N, q->b_off));
    if (q->b_off)  md_free(alloc, q->b_off, sizeof(uint32_t) * (N + 1));
    if (q->radius) md_free(alloc, q->radius, sizeof(float) * MAX(N, 1));
    if (q->parent) md_free(alloc, q->parent, sizeof(uint32_t) * MAX(q->num_a, 1));
    if (q->type)   md_free(alloc, q->type, sizeof(uint32_t) * MAX(N, 1));
    if (q->type_cutoff) md_free(alloc, q->type_cutoff, sizeof(float) * MAX((size_t)q->num_types * q->num_types, 1));
    MEMSET(q, 0, sizeof(md_contact_query_t));
}

bool md_contact_query_init(md_contact_query_t* q, const md_contact_desc_t* desc, const md_system_t* sys, md_allocator_i* alloc) {
    ASSERT(q);
    MEMSET(q, 0, sizeof(md_contact_query_t));
    if (!desc || !sys || !alloc) {
        MD_LOG_ERROR("Contact query: missing description, system or allocator");
        return false;
    }
    if (desc->num_a && !desc->group_a) {
        MD_LOG_ERROR("Contact query: missing groups");
        return false;
    }
    if (desc->num_a > UINT32_MAX - 1 || desc->num_b > UINT32_MAX - 1 || sys->atom.count > UINT32_MAX - 1) {
        MD_LOG_ERROR("Contact query: too many groups or particles");
        return false;
    }
    if (desc->criterion != MD_CONTACT_CRITERION_DISTANCE && desc->criterion != MD_CONTACT_CRITERION_RADII && desc->criterion != MD_CONTACT_CRITERION_TYPE_PAIR) {
        MD_LOG_ERROR("Contact query: unknown criterion");
        return false;
    }
    if (desc->criterion == MD_CONTACT_CRITERION_TYPE_PAIR) {
        if (!desc->particle_type || !desc->type_cutoff || desc->num_types == 0 || desc->num_types > 65535) {
            MD_LOG_ERROR("Contact query: the type pair criterion needs particle types and a table of contact distances");
            return false;
        }
        const uint32_t T = desc->num_types;
        for (uint32_t ta = 0; ta < T; ++ta) {
            for (uint32_t tb = ta + 1; tb < T; ++tb) {
                if (desc->type_cutoff[ta * T + tb] != desc->type_cutoff[tb * T + ta]) {
                    MD_LOG_ERROR("Contact query: the table of contact distances is not symmetric (types %u and %u)", ta, tb);
                    return false;
                }
            }
        }
    }

    if (!desc->group_b && desc->num_b) {
        MD_LOG_ERROR("Contact query: a count of B groups without the groups");
        return false;
    }
    // An empty B is a set with no contacts, not the absence of B: that would be the contacts within A
    const bool self = desc->group_b == NULL;
    q->alloc     = alloc;
    q->num_atoms = (uint32_t)sys->atom.count;
    q->num_a     = (uint32_t)desc->num_a;
    q->num_b     = self ? (uint32_t)desc->num_a : (uint32_t)desc->num_b;
    q->flags     = self ? MD_CONTACT_FLAG_SELF : MD_CONTACT_FLAG_NONE;
    q->criterion = desc->criterion;
    q->cutoff    = (float)desc->cutoff;
    q->min_separation = self ? desc->min_separation : 0;

    const uint32_t N = q->num_atoms;
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    bool result = false;

    if (self && desc->group_parent && q->num_a) {
        q->parent = md_alloc(alloc, sizeof(uint32_t) * q->num_a);
        MEMCPY(q->parent, desc->group_parent, sizeof(uint32_t) * q->num_a);
    }

    md_bitfield_t set_a = {0};
    md_bitfield_t set_b = {0};
    if (!build_membership(&q->a_off, &q->a_grp, &set_a, desc->group_a, desc->num_a, N, alloc, temp_alloc)) goto done;
    if (!self && !build_membership(&q->b_off, &q->b_grp, &set_b, desc->group_b, desc->num_b, N, alloc, temp_alloc)) goto done;

    // The pairs are searched out to the furthest the criterion can reach
    float search_radius = q->cutoff;
    float max_radius = 0.0f;
    if (desc->criterion == MD_CONTACT_CRITERION_TYPE_PAIR) {
        const uint32_t T = desc->num_types;
        q->num_types = T;
        q->type = md_alloc(alloc, sizeof(uint32_t) * MAX(N, 1));
        if (N) MEMCPY(q->type, desc->particle_type, sizeof(uint32_t) * N);
        q->type_cutoff = md_alloc(alloc, sizeof(float) * T * T);
        MEMCPY(q->type_cutoff, desc->type_cutoff, sizeof(float) * T * T);

        // Only the types which are present bound the search
        uint8_t* present = md_alloc(temp_alloc, T);
        MEMSET(present, 0, T);
        for (uint32_t k = 0; k < N; ++k) {
            const bool member = (q->a_off[k + 1] > q->a_off[k]) || (q->b_off && q->b_off[k + 1] > q->b_off[k]);
            if (!member) continue;
            if (q->type[k] >= T) {
                MD_LOG_ERROR("Contact query: particle %u has type %u, beyond the %u types of the table", k, q->type[k], T);
                goto done;
            }
            present[q->type[k]] = 1;
        }
        search_radius = 0.0f;
        for (uint32_t ta = 0; ta < T; ++ta) {
            for (uint32_t tb = 0; tb < T; ++tb) {
                if (present[ta] && present[tb]) search_radius = MAX(search_radius, q->type_cutoff[ta * T + tb]);
            }
        }
    } else if (desc->criterion == MD_CONTACT_CRITERION_RADII) {
        q->radius = md_alloc(alloc, sizeof(float) * MAX(N, 1));
        if (desc->atom_radius) {
            MEMCPY(q->radius, desc->atom_radius, sizeof(float) * N);
        } else if (N) {
            md_atom_extract_radii(q->radius, 0, N, &sys->atom);
        }
        for (uint32_t k = 0; k < N; ++k) {
            const bool member = (q->a_off[k + 1] > q->a_off[k]) || (q->b_off && q->b_off[k + 1] > q->b_off[k]);
            if (member) max_radius = MAX(max_radius, q->radius[k]);
        }
        search_radius = q->cutoff + 2.0f * max_radius;
    }

    const md_contact_pairs_desc_t pairs_desc = {
        .set_a = &set_a,
        .set_b = self ? NULL : &set_b,
        .radius = search_radius,
        .exclude_bonds = desc->exclude_bonds,
        .particle_label = desc->particle_label,
    };
    result = md_contact_pairs_init(&q->pairs, &pairs_desc, sys, alloc);

done:
    md_temp_end(temp);
    if (!result) {
        md_contact_query_free(q);
    }
    return result;
}

typedef struct group_accum_t {
    const md_contact_query_t* q;
    md_hashmap32_t map;             // mix_pair_key(i << 32 | j) -> index into the arrays below
    md_array(uint64_t) key;
    md_array(uint32_t) atom_pairs;
    md_array(float)    d_min;
    md_allocator_i* alloc;
} group_accum_t;

// The hash map indexes by the low bits of (hi ^ lo) of the key. For a pair key that is i ^ j, and contacts are
// mostly between nearby groups, which would pile them into a few buckets. So the key is mixed first; the mix is a
// bijection, so distinct pairs keep distinct keys. The two values the map reserves are moved out of the way, which
// could only collide with another pair with a probability of 2^-63.
static inline uint64_t mix_pair_key(uint64_t x) {
    x ^= x >> 30;
    x *= 0xbf58476d1ce4e5b9ull;
    x ^= x >> 27;
    x *= 0x94d049bb133111ebull;
    x ^= x >> 31;
    if (x >= MD_HASH_TOMBSTONE) x ^= (1ull << 63);
    return x;
}

static inline void accumulate(group_accum_t* acc, uint32_t i, uint32_t j, float d) {
    const uint64_t key = ((uint64_t)i << 32) | j;
    const uint64_t hkey = mix_pair_key(key);
    uint32_t* idx = md_hashmap_get(&acc->map, hkey);
    if (idx) {
        acc->atom_pairs[*idx] += 1;
        acc->d_min[*idx] = MIN(acc->d_min[*idx], d);
    } else {
        const uint32_t new_idx = (uint32_t)md_array_size(acc->key);
        md_array_push(acc->key, key, acc->alloc);
        md_array_push(acc->atom_pairs, 1, acc->alloc);
        md_array_push(acc->d_min, d, acc->alloc);
        md_hashmap_add(&acc->map, hkey, new_idx);
    }
}

static void group_callback(const uint32_t* pa, const uint32_t* pb, const float* pr, size_t count, void* user_param) {
    group_accum_t* acc = (group_accum_t*)user_param;
    const md_contact_query_t* q = acc->q;
    const bool self = q->flags & MD_CONTACT_FLAG_SELF;
    const uint32_t* b_off = self ? q->a_off : q->b_off;
    const uint32_t* b_grp = self ? q->a_grp : q->b_grp;

    for (size_t k = 0; k < count; ++k) {
        const uint32_t x = pa[k];
        const uint32_t y = pb[k];
        const float d = pr[k];
        float limit = q->cutoff;
        if (q->criterion == MD_CONTACT_CRITERION_RADII) {
            limit = q->radius[x] + q->radius[y] + q->cutoff;
        } else if (q->criterion == MD_CONTACT_CRITERION_TYPE_PAIR) {
            limit = q->type_cutoff[q->type[x] * q->num_types + q->type[y]];
        }
        if (!(d < limit)) continue;

        for (uint32_t ga = q->a_off[x]; ga < q->a_off[x + 1]; ++ga) {
            for (uint32_t gb = b_off[y]; gb < b_off[y + 1]; ++gb) {
                uint32_t i = q->a_grp[ga];
                uint32_t j = b_grp[gb];
                if (self) {
                    if (i == j) continue;
                    if (i > j) { const uint32_t t = i; i = j; j = t; }
                    if (q->min_separation > 1 && j - i < q->min_separation && (!q->parent || q->parent[i] == q->parent[j])) continue;
                }
                accumulate(acc, i, j, d);
            }
        }
    }
}

typedef struct key_idx_t {
    uint64_t key;
    uint32_t idx;
} key_idx_t;

static int compare_key_idx(const void* a, const void* b) {
    const uint64_t ka = ((const key_idx_t*)a)->key;
    const uint64_t kb = ((const key_idx_t*)b)->key;
    return (ka > kb) - (ka < kb);
}

void md_contact_set_free(md_contact_set_t* set) {
    if (!set) return;
    if (set->alloc && set->count) {
        md_free(set->alloc, set->i, sizeof(uint32_t) * set->count);
        md_free(set->alloc, set->j, sizeof(uint32_t) * set->count);
        md_free(set->alloc, set->atom_pairs, sizeof(uint32_t) * set->count);
        md_free(set->alloc, set->d_min, sizeof(float) * set->count);
    }
    MEMSET(set, 0, sizeof(md_contact_set_t));
}

bool md_contact_query_eval(md_contact_set_t* out, const md_contact_query_t* q, const md_system_state_t* state, md_allocator_i* alloc) {
    ASSERT(out);
    MEMSET(out, 0, sizeof(md_contact_set_t));
    if (!q || !q->alloc || !state || !alloc) {
        MD_LOG_ERROR("Contact evaluation: missing query, state or allocator");
        return false;
    }

    out->num_a = q->num_a;
    out->num_b = q->num_b;
    out->flags = q->flags;
    out->alloc = alloc;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    group_accum_t acc = {
        .q = q,
        .map = { .allocator = temp_alloc },
        .alloc = temp_alloc,
    };

    bool result = md_contact_pairs_for_each(&q->pairs, state, group_callback, &acc);

    const size_t count = md_array_size(acc.key);
    if (result && count) {
        key_idx_t* order = md_alloc(temp_alloc, sizeof(key_idx_t) * count);
        for (size_t k = 0; k < count; ++k) {
            order[k] = (key_idx_t){ acc.key[k], (uint32_t)k };
        }
        qsort(order, count, sizeof(key_idx_t), compare_key_idx);

        out->count      = count;
        out->i          = md_alloc(alloc, sizeof(uint32_t) * count);
        out->j          = md_alloc(alloc, sizeof(uint32_t) * count);
        out->atom_pairs = md_alloc(alloc, sizeof(uint32_t) * count);
        out->d_min      = md_alloc(alloc, sizeof(float) * count);
        for (size_t k = 0; k < count; ++k) {
            const uint32_t src = order[k].idx;
            out->i[k]          = (uint32_t)(order[k].key >> 32);
            out->j[k]          = (uint32_t)(order[k].key & 0xFFFFFFFFu);
            out->atom_pairs[k] = acc.atom_pairs[src];
            out->d_min[k]      = acc.d_min[src];
        }
    }

    md_temp_end(temp);
    return result;
}

bool md_contact_compute(md_contact_set_t* out, const md_contact_desc_t* desc, const md_system_t* sys, const md_system_state_t* state, md_allocator_i* alloc) {
    ASSERT(out);
    MEMSET(out, 0, sizeof(md_contact_set_t));
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_contact_query_t q;
    bool result = md_contact_query_init(&q, desc, sys, md_temp_allocator(temp)) && md_contact_query_eval(out, &q, state, alloc);
    md_temp_end(temp);
    return result;
}

// ### CONTACT REGIONS ###

static int compare_u64(const void* a, const void* b) {
    const uint64_t x = *(const uint64_t*)a;
    const uint64_t y = *(const uint64_t*)b;
    return (x > y) - (x < y);
}

// Sorts and removes duplicates, returns the new count
static size_t sort_unique_u64(uint64_t* keys, size_t count) {
    if (count == 0) return 0;
    qsort(keys, count, sizeof(uint64_t), compare_u64);
    size_t n = 1;
    for (size_t k = 1; k < count; ++k) {
        if (keys[k] != keys[n - 1]) keys[n++] = keys[k];
    }
    return n;
}

bool md_contact_unit_adjacency_from_bonds(uint32_t** out_off, uint32_t** out_adj, const md_bitfield_t* units, size_t num_units, const md_system_t* sys, md_allocator_i* alloc) {
    if (!out_off || !out_adj || (num_units && !units) || !sys || !alloc || num_units > UINT32_MAX - 1) {
        MD_LOG_ERROR("Contact unit adjacency: invalid arguments");
        return false;
    }
    const size_t N = sys->atom.count;
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    // The unit of each particle, the first which holds it
    uint32_t* unit_of = md_alloc(temp_alloc, sizeof(uint32_t) * MAX(N, 1));
    for (size_t k = 0; k < N; ++k) unit_of[k] = UINT32_MAX;
    for (size_t u = 0; u < num_units; ++u) {
        md_bitfield_iter_t it = md_bitfield_iter_create(&units[u]);
        while (md_bitfield_iter_next(&it)) {
            const uint64_t k = md_bitfield_iter_idx(&it);
            if (k < N && unit_of[k] == UINT32_MAX) unit_of[k] = (uint32_t)u;
        }
    }

    // Both directions of every bond joining two units
    md_array(uint64_t) edge = 0;
    for (size_t b = 0; b < sys->bond.count && sys->bond.pairs; ++b) {
        const int64_t x = sys->bond.pairs[b].idx[0];
        const int64_t y = sys->bond.pairs[b].idx[1];
        if (x < 0 || y < 0 || (size_t)x >= N || (size_t)y >= N) continue;
        const uint32_t ux = unit_of[x];
        const uint32_t uy = unit_of[y];
        if (ux == UINT32_MAX || uy == UINT32_MAX || ux == uy) continue;
        md_array_push(edge, ((uint64_t)ux << 32) | uy, temp_alloc);
        md_array_push(edge, ((uint64_t)uy << 32) | ux, temp_alloc);
    }
    const size_t num_edges = sort_unique_u64(edge, md_array_size(edge));

    uint32_t* off = md_alloc(alloc, sizeof(uint32_t) * (num_units + 1));
    uint32_t* adj = md_alloc(alloc, sizeof(uint32_t) * MAX(num_edges, 1));
    MEMSET(off, 0, sizeof(uint32_t) * (num_units + 1));
    for (size_t e = 0; e < num_edges; ++e) {
        off[(edge[e] >> 32) + 1] += 1;
        adj[e] = (uint32_t)(edge[e] & 0xFFFFFFFFu);     // Sorted by (from, to): the rows come out in order
    }
    for (size_t u = 0; u < num_units; ++u) off[u + 1] += off[u];

    *out_off = off;
    *out_adj = adj;
    md_temp_end(temp);
    return true;
}

void md_contact_region_set_free(md_contact_region_set_t* set) {
    if (!set) return;
    if (set->alloc) {
        const size_t n = set->count;
        if (n) {
            md_free(set->alloc, set->body_i, sizeof(uint32_t) * n);
            md_free(set->alloc, set->body_j, sizeof(uint32_t) * n);
            md_free(set->alloc, set->unit_pairs, sizeof(uint32_t) * n);
            md_free(set->alloc, set->units_i, sizeof(uint32_t) * n);
            md_free(set->alloc, set->units_j, sizeof(uint32_t) * n);
            md_free(set->alloc, set->atom_pairs, sizeof(uint32_t) * n);
            md_free(set->alloc, set->d_min, sizeof(float) * n);
        }
        if (set->num_pairs) md_free(set->alloc, set->region, sizeof(uint32_t) * set->num_pairs);
    }
    MEMSET(set, 0, sizeof(md_contact_region_set_t));
}

// The units within reach steps of u in its body, u included, written to out (an md_array which is cleared first).
// stamp has one entry per unit, and marks the units visited in this search with the value stamp_value.
static void units_within_reach(md_array(uint32_t)* out, const md_contact_units_t* units, uint32_t u, uint32_t reach,
                               uint32_t* stamp, uint32_t stamp_value, md_array(uint32_t)* depth, md_allocator_i* alloc) {
    md_array_shrink(*out, 0);
    md_array_shrink(*depth, 0);
    md_array_push(*out, u, alloc);
    md_array_push(*depth, 0, alloc);
    stamp[u] = stamp_value;
    if (!units->adj_off || !units->adj) return;
    const uint32_t body = units->body[u];
    for (size_t head = 0; head < md_array_size(*out); ++head) {
        const uint32_t cur = (*out)[head];
        const uint32_t d = (*depth)[head];
        if (d >= reach) continue;
        for (uint32_t e = units->adj_off[cur]; e < units->adj_off[cur + 1]; ++e) {
            const uint32_t next = units->adj[e];
            if (next >= units->count || stamp[next] == stamp_value || units->body[next] != body) continue;
            stamp[next] = stamp_value;
            md_array_push(*out, next, alloc);
            md_array_push(*depth, d + 1, alloc);
        }
    }
}

static inline uint32_t uf_find(uint32_t* parent, uint32_t x) {
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

static inline void uf_union(uint32_t* parent, uint32_t a, uint32_t b) {
    a = uf_find(parent, a);
    b = uf_find(parent, b);
    if (a != b) {
        // The smaller index becomes the root, which keeps the result independent of the order of the unions
        if (a < b) parent[b] = a; else parent[a] = b;
    }
}

// The index of the unit pair (x, y) in the contact set, which is sorted by (i, j), or UINT32_MAX
static uint32_t find_unit_pair(const md_contact_set_t* set, uint32_t x, uint32_t y) {
    if (set->flags & MD_CONTACT_FLAG_SELF) {
        if (x == y) return UINT32_MAX;
        if (x > y) { const uint32_t t = x; x = y; y = t; }
    }
    const uint64_t key = ((uint64_t)x << 32) | y;
    size_t lo = 0, hi = set->count;
    while (lo < hi) {
        const size_t mid = (lo + hi) / 2;
        const uint64_t k = ((uint64_t)set->i[mid] << 32) | set->j[mid];
        if (k == key) return (uint32_t)mid;
        if (k < key) lo = mid + 1; else hi = mid;
    }
    return UINT32_MAX;
}

static bool units_valid(const md_contact_units_t* units, size_t expected) {
    if (!units || !units->body || units->count != expected) return false;
    if ((units->adj_off == NULL) != (units->adj == NULL)) return false;
    return true;
}

bool md_contact_regions(md_contact_region_set_t* out, const md_contact_set_t* contacts, const md_contact_units_t* units_a, const md_contact_units_t* units_b, uint32_t reach, md_allocator_i* alloc) {
    ASSERT(out);
    MEMSET(out, 0, sizeof(md_contact_region_set_t));
    if (!contacts || !alloc) {
        MD_LOG_ERROR("Contact regions: missing contact set or allocator");
        return false;
    }
    const bool self = contacts->flags & MD_CONTACT_FLAG_SELF;
    if (self) units_b = units_a;
    if (!units_valid(units_a, contacts->num_a) || !units_valid(units_b, contacts->num_b)) {
        MD_LOG_ERROR("Contact regions: the units do not match the groups of the contact set (a body per group, and both or neither adjacency array)");
        return false;
    }
    if (reach == 0) {
        MD_LOG_ERROR("Contact regions: reach has to be at least 1");
        return false;
    }

    const size_t n = contacts->count;
    out->alloc = alloc;
    out->num_pairs = n;
    if (n == 0) return true;
    if (n > UINT32_MAX - 1) {
        MD_LOG_ERROR("Contact regions: too many unit pairs");
        return false;
    }

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    // Join each unit pair with the unit pairs in contact among its neighbours on both sides
    uint32_t* parent = md_alloc(temp_alloc, sizeof(uint32_t) * n);
    for (uint32_t k = 0; k < n; ++k) parent[k] = k;

    uint32_t* stamp_a = md_alloc(temp_alloc, sizeof(uint32_t) * MAX(units_a->count, 1));
    uint32_t* stamp_b = self ? stamp_a : md_alloc(temp_alloc, sizeof(uint32_t) * MAX(units_b->count, 1));
    MEMSET(stamp_a, 0, sizeof(uint32_t) * MAX(units_a->count, 1));
    if (!self) MEMSET(stamp_b, 0, sizeof(uint32_t) * MAX(units_b->count, 1));
    uint32_t stamp_value = 0;

    md_array(uint32_t) near_x = 0;
    md_array(uint32_t) near_y = 0;
    md_array(uint32_t) depth  = 0;
    for (uint32_t k = 0; k < n; ++k) {
        const uint32_t x = contacts->i[k];
        const uint32_t y = contacts->j[k];
        units_within_reach(&near_x, units_a, x, reach, stamp_a, ++stamp_value, &depth, temp_alloc);
        units_within_reach(&near_y, units_b, y, reach, stamp_b, ++stamp_value, &depth, temp_alloc);
        for (size_t a = 0; a < md_array_size(near_x); ++a) {
            for (size_t b = 0; b < md_array_size(near_y); ++b) {
                const uint32_t m = find_unit_pair(contacts, near_x[a], near_y[b]);
                if (m != UINT32_MAX && m != k) uf_union(parent, k, m);
            }
        }
    }

    // Number the regions in order of their first unit pair
    uint32_t* region = md_alloc(alloc, sizeof(uint32_t) * n);
    uint32_t* root_region = md_alloc(temp_alloc, sizeof(uint32_t) * n);
    for (uint32_t k = 0; k < n; ++k) root_region[k] = UINT32_MAX;
    uint32_t num_regions = 0;
    for (uint32_t k = 0; k < n; ++k) {
        const uint32_t r = uf_find(parent, k);
        if (root_region[r] == UINT32_MAX) root_region[r] = num_regions++;
        region[k] = root_region[r];
    }

    out->count      = num_regions;
    out->region     = region;
    out->body_i     = md_alloc(alloc, sizeof(uint32_t) * num_regions);
    out->body_j     = md_alloc(alloc, sizeof(uint32_t) * num_regions);
    out->unit_pairs = md_alloc(alloc, sizeof(uint32_t) * num_regions);
    out->units_i    = md_alloc(alloc, sizeof(uint32_t) * num_regions);
    out->units_j    = md_alloc(alloc, sizeof(uint32_t) * num_regions);
    out->atom_pairs = md_alloc(alloc, sizeof(uint32_t) * num_regions);
    out->d_min      = md_alloc(alloc, sizeof(float) * num_regions);
    MEMSET(out->unit_pairs, 0, sizeof(uint32_t) * num_regions);
    MEMSET(out->atom_pairs, 0, sizeof(uint32_t) * num_regions);
    for (uint32_t r = 0; r < num_regions; ++r) out->d_min[r] = FLT_MAX;

    // The units on each side, as (region, unit) keys
    md_array(uint64_t) side_i = 0;
    md_array(uint64_t) side_j = 0;
    for (uint32_t k = 0; k < n; ++k) {
        const uint32_t r = region[k];
        uint32_t x = contacts->i[k];
        uint32_t y = contacts->j[k];
        uint32_t bx = units_a->body[x];
        uint32_t by = units_b->body[y];
        if (self && bx > by) {
            // Oriented so the i side is the lower body
            uint32_t t = x; x = y; y = t;
            t = bx; bx = by; by = t;
        }
        out->body_i[r] = bx;
        out->body_j[r] = by;
        out->unit_pairs[r] += 1;
        out->atom_pairs[r] += contacts->atom_pairs[k];
        out->d_min[r] = MIN(out->d_min[r], contacts->d_min[k]);
        if (self && bx == by) {
            // Within one body the two sides are the same units
            md_array_push(side_i, ((uint64_t)r << 32) | x, temp_alloc);
            md_array_push(side_i, ((uint64_t)r << 32) | y, temp_alloc);
        } else {
            md_array_push(side_i, ((uint64_t)r << 32) | x, temp_alloc);
            md_array_push(side_j, ((uint64_t)r << 32) | y, temp_alloc);
        }
    }
    MEMSET(out->units_i, 0, sizeof(uint32_t) * num_regions);
    MEMSET(out->units_j, 0, sizeof(uint32_t) * num_regions);
    const size_t ni = sort_unique_u64(side_i, md_array_size(side_i));
    const size_t nj = sort_unique_u64(side_j, md_array_size(side_j));
    for (size_t k = 0; k < ni; ++k) out->units_i[side_i[k] >> 32] += 1;
    for (size_t k = 0; k < nj; ++k) out->units_j[side_j[k] >> 32] += 1;
    for (uint32_t r = 0; r < num_regions; ++r) {
        if (self && out->body_i[r] == out->body_j[r]) out->units_j[r] = out->units_i[r];
    }

    md_temp_end(temp);
    return true;
}
