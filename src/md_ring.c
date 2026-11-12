#include "md_infer_internal.h"

// Ring inference: the rings of the bond graph and their flags (md_util.h).

static void ring_flags(md_ring_data_t* rings, const md_system_t* sys);

static inline bool compare_ring(const int* ring_a, int ring_a_size, const int* ring_b, int ring_b_size) {
    if (ring_a_size != ring_b_size) {
        return false;
    }
    for (int i = 0; i < ring_a_size; ++i) {
        if (ring_a[i] != ring_b[i]) {
            return false;
        }
    }
    return true;
}

static inline bool has_ring(const md_ring_data_t* ring_data, const int* ring, int ring_size) {
    const int64_t num_rings = md_ring_data_count(ring_data);
    for (int64_t i = 0; i < num_rings; ++i) {
        const int* ring_i     = md_ring_atoms(ring_data, i);
        const int ring_i_size = (int)md_ring_size(ring_data, i);
        if (compare_ring(ring_i, ring_i_size, ring, ring_size)) return true;
    }
    return false;
}

#define SORT_RING_PRINT 0

// This is a specialized deterministic sort function for reordering the indices of a ring in a particular way such that it can be subjected to md_memcmp
static void sort_ring(int* arr, int n) {
    MD_ASSERT(n >= 0);
    if (n == 0) return;

#if SORT_RING_PRINT
    printf("pre:  [");
    for (int i = 0; i < n; ++i) {
        char delim = i < (n-1) ? ' ' : ']';
        printf("%i%c", arr[i], delim);
    }
    printf("\n");
#endif

    int lowest_idx = 0;
    for (int i = 1; i < n; ++i) {
        if (arr[i] < arr[lowest_idx]) {
            lowest_idx = i;
        }
    }

    // Rotate ring such that the lowest index in the first slot
    if (lowest_idx != 0) {
        md_temp_scope_t temp = md_temp_begin();
        int* buf = md_temp_alloc_array(temp, int, n);
        int j = 0;
        for (int i = lowest_idx; i < n; ++i) {
            buf[j++] = arr[i];
        }
        for (int i = 0; i < lowest_idx; ++i) {
            buf[j++] = arr[i];
        }
        md_memcpy(arr, buf, sizeof(*arr) * n);
        md_temp_end(temp);
    }

#if SORT_RING_PRINT
    printf("mid:  [");
    for (int i = 0; i < n; ++i) {
        char delim = i < (n-1) ? ' ' : ']';
        printf("%i%c", arr[i], delim);
    }
    printf("\n");
#endif

    // And we want to store the lowest neighboring index in the second slot (right in this cyclic array)
    if (arr[n - 1] < arr[1]) {
        // Need to swap
        for (int i = 1, j = n - 1; i < j; ++i, --j) {
            int tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;
        }
    }
#if SORT_RING_PRINT
    printf("post: [");
    for (int i = 0; i < n; ++i) {
        char delim = i < (n-1) ? ' ' : ']';
        printf("%i%c", arr[i], delim);
    }
    printf("\n\n");
#endif
}

// Rings of up to 7 atoms: the 7 membered aromatic rings (tropolone, azulene, tropylium) count, larger ones are loops
#define MIN_RING_SIZE MD_RING_MIN_SIZE
#define MAX_RING_SIZE MD_RING_MAX_SIZE
#define MAX_DEPTH (MAX_RING_SIZE/2 + 1)

// Removes the rings which are sums of smaller ones: rings whose bonds are those of strictly smaller rings, each bond
// counted an odd number of times (the 7 membered ring around a cyclopropane fused to a cyclohexane, the 6 membered one
// around the two 5 membered rings of norbornane, the 4 membered one around the two 3 membered rings of
// bicyclo[1.1.0]butane, the 6 membered rings around two or three faces of cubane). What remains are the relevant rings
// of the graph up to the size limit: the smallest set of smallest rings, made symmetric (each of the three 6 membered
// rings of bicyclo[2.2.2]octane is the sum of the other two, which are not smaller, and all three are kept).
// The smaller rings considered are those sharing an atom with the ring, which every ring of such a sum does when the
// rings are this small.
#define ENVELOPE_MAX_CAND 64
#define ENVELOPE_MAX_BONDS ((ENVELOPE_MAX_CAND + 1) * MAX_RING_SIZE)
#define ENVELOPE_WORDS ((ENVELOPE_MAX_BONDS + 63) / 64)
static void remove_envelope_rings(md_ring_data_t* rings, const md_bond_data_t* bond, size_t num_atoms, md_allocator_i* temp_arena) {
    const size_t num_rings = md_ring_data_count(rings);
    if (num_rings < 2) return;

    // The rings of each atom
    uint32_t* off = md_alloc(temp_arena, (num_atoms + 1) * sizeof(uint32_t));
    md_memset(off, 0, (num_atoms + 1) * sizeof(uint32_t));
    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t* ring = md_ring_atoms(rings, r);
        const size_t n = md_ring_size(rings, r);
        for (size_t k = 0; k < n; ++k) off[ring[k] + 1] += 1;
    }
    for (size_t i = 0; i < num_atoms; ++i) off[i + 1] += off[i];
    uint32_t* idx = md_alloc(temp_arena, MD_MAX(off[num_atoms], 1) * sizeof(uint32_t));
    uint32_t* fill = md_alloc(temp_arena, num_atoms * sizeof(uint32_t));
    md_memcpy(fill, off, num_atoms * sizeof(uint32_t));
    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t* ring = md_ring_atoms(rings, r);
        const size_t n = md_ring_size(rings, r);
        for (size_t k = 0; k < n; ++k) idx[fill[ring[k]]++] = (uint32_t)r;
    }

    uint8_t* removed = md_alloc(temp_arena, num_rings);
    md_memset(removed, 0, num_rings);
    size_t removed_count = 0;

    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t* ring = md_ring_atoms(rings, r);
        const int n = (int)md_ring_size(rings, r);
        if (n < 4) continue;    // Two rings of 3 make one of 4 at the least

        // Smaller rings which share an atom with this one
        uint32_t cand[ENVELOPE_MAX_CAND];
        int num_cand = 0;
        for (int k = 0; k < n; ++k) {
            for (uint32_t t = off[ring[k]]; t < off[ring[k] + 1]; ++t) {
                const uint32_t s = idx[t];
                if (s == r || (int)md_ring_size(rings, s) >= n) continue;
                bool dup = false;
                for (int c = 0; c < num_cand; ++c) dup |= cand[c] == s;
                if (!dup && num_cand < ENVELOPE_MAX_CAND) cand[num_cand++] = s;
            }
        }
        if (num_cand < 2) continue;

        // The bonds, numbered locally, as bit vectors over GF(2): the ring is a sum of the candidates when it reduces to
        // nothing by them
        int32_t bonds[ENVELOPE_MAX_BONDS];
        int num_bonds = 0;
        uint64_t basis[ENVELOPE_MAX_BONDS][ENVELOPE_WORDS];   // basis[p]: the row whose lowest bit is p
        uint8_t  has_row[ENVELOPE_MAX_BONDS];
        md_memset(has_row, 0, sizeof(has_row));
        bool valid = true;
        bool envelope = false;
        for (int c = 0; c <= num_cand && valid; ++c) {
            const bool self = c == num_cand;
            const uint32_t s = self ? (uint32_t)r : cand[c];
            const int32_t* rs = md_ring_atoms(rings, s);
            const int ns = (int)md_ring_size(rings, s);
            uint64_t v[ENVELOPE_WORDS] = {0};
            for (int k = 0; k < ns; ++k) {
                const md_bond_idx_t b = md_bond_find(bond, rs[k], rs[(k + 1) % ns]);
                if (b < 0) { valid = false; break; }
                int j = 0;
                while (j < num_bonds && bonds[j] != b) ++j;
                if (j == num_bonds) bonds[num_bonds++] = b;
                v[j >> 6] ^= 1ull << (j & 63);
            }
            if (!valid) break;
            // Reduce by the rows so far
            int p = -1;
            for (;;) {
                p = -1;
                for (int w = 0; w < ENVELOPE_WORDS; ++w) {
                    if (v[w]) { p = w * 64 + (int)md_ctz64(v[w]); break; }
                }
                if (p < 0 || !has_row[p]) break;
                for (int w = 0; w < ENVELOPE_WORDS; ++w) v[w] ^= basis[p][w];
            }
            if (self) {
                envelope = p < 0;
            } else if (p >= 0) {
                md_memcpy(basis[p], v, sizeof(v));
                has_row[p] = 1;
            }
        }
        if (valid && envelope) {
            removed[r] = 1;
            removed_count += 1;
        }
    }
    if (!removed_count) return;

    size_t kept = 0;
    for (size_t r = 0; r < num_rings; ++r) {
        if (!removed[r]) rings->rings[kept++] = rings->rings[r];
    }
    md_array_shrink(rings->rings, kept);
}
#undef ENVELOPE_MAX_CAND
#undef ENVELOPE_MAX_BONDS
#undef ENVELOPE_WORDS

// Inspired by molstars implementation
// https://github.com/molstar/molstar/blob/master/src/mol-model/structure/structure/unit/rings/compute.ts#L249
// Copyright (c) 2018 mol* contributors, licensed under MIT, See LICENSE file for more info.
// @author David Sehnal <david.sehnal@gmail.com>

// The rings of the bond graph (bond, with its connectivity), pushed to out: see md_ring_infer. The atoms
// set in skip (optional) are neither roots nor passed through.
static void find_rings(md_ring_data_t* out, const md_bond_data_t* bond, size_t num_atoms, const uint64_t* skip, md_allocator_i* alloc, md_temp_scope_t temp) {
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    // Some typedefs to enable laboration of types
    typedef uint16_t color_t;
    typedef uint8_t  depth_t;
    typedef uint16_t mark_t;
    typedef  int32_t pred_t;
        
    color_t* color = md_temp_alloc_zero_array(temp, color_t, num_atoms);
    depth_t* depth = md_temp_alloc_zero_array(temp, depth_t, num_atoms);
    mark_t*  mark  = md_temp_alloc_zero_array(temp, mark_t, num_atoms);
    pred_t*  pred  = md_temp_alloc_zero_array(temp, pred_t, num_atoms);
    md_memset(pred, -1, num_atoms * sizeof(pred_t));    // We can do memset as the representation of -1 under two's complement is 0xFFFFFFFF

    // The capacity is arbitrary here, but will be resized if needed.
    md_fifo_t queue = md_fifo_create(64, temp_arena);

    md_hashset_t ring_set = {.allocator = temp_arena};
    
    color_t current_color = 1;
    mark_t  current_mark  = 1;

    for (int atom_idx = 0; atom_idx < (int)num_atoms; ++atom_idx) {
        if (skip && bitfield_test_bit(skip, atom_idx)) continue;

        // Skip any atom that has already been colored in the previous search
        if (color[atom_idx] == current_color) continue;

        // Set i as our root
        depth[atom_idx] = 1;
        pred [atom_idx] = -1;
        mark [atom_idx] = current_mark++;

        md_fifo_clear(&queue);
        md_fifo_push(&queue, atom_idx);
        while (!md_fifo_empty(&queue)) {
            int idx = md_fifo_pop(&queue);

            md_bond_iter_t it = md_bond_iter(bond, idx);
            while (md_bond_iter_has_next(&it)) {
                int next = md_bond_iter_atom_index(&it);
                md_bond_iter_next(&it);

                //if (next > idx) continue;
                if (next == pred[idx]) continue;  // avoid adding parent to search queue
                if (skip && bitfield_test_bit(skip, next)) continue;

                if (mark[next] == current_mark && depth[idx] >= MIN_RING_SIZE) {
                    // We found a junction point where the graph connects
                    // Now we need to traverse both branches of this graph back up
                    // In order to find the other junction where the graph branches
                    // We follow both branches up towards the root.

                    int l = idx;
                    int r = next;

                    // Only process one of the two branches/cases
                    if (r > l) continue;

                    int len = 0;
                    int ring[MAX_DEPTH * 2];

                    color_t col = current_color++;
                    int cur;

                    cur = l;
                    for (int d = 0; d < MAX_DEPTH; ++d) {
                        color[cur] = col;
                        cur = pred[cur];
                        if (cur < 0) break;
                    }

                    int target = -1;
                    cur = r;
                    for (int d = 0; d < MAX_DEPTH; ++d) {
                        if (color[cur] == col) {
                            target = cur;
                            break;
                        }
                        ring[len++] = cur;
                        cur = pred[cur];
                        if (cur == -1) break;
                    }

                    // No junction point found
                    if (target == -1) continue;

                    cur = l;
                    int left_beg = len;
                    for (int d = 0; d < MAX_DEPTH; ++d) {
                        ring[len++] = cur;
                        if (target == cur) break;
                        cur = pred[cur];
                        if (cur == -1) break;
                    }

                    // Otherwise we made a big whoopsie
                    MD_ASSERT(len < (int)MD_COUNTOF(ring));

                    if (MIN_RING_SIZE <= len && len <= MAX_RING_SIZE) {
                        // Ther order of the left branch is reversed so we need to swap that
                        // Such that the indices of the ring sequentially iterates the ring
                        for (int i = left_beg, j = len - 1; i < j; ++i, --j) {
                            int tmp = ring[i];
                            ring[i] = ring[j];
                            ring[j] = tmp;
                        }

                        // @TODO: Implement a custom sort for rings such that the indices within a ring
                        // Are predictable and can be easily md_memcmp or hashed
                        // But does not violate the sequential indexing such that it can be iterated
                        sort_ring(ring, len);
                        uint64_t key = md_hash64(ring, len * sizeof(*ring), 0);
                        if (!md_hashset_get(&ring_set, key)) {
                            md_hashset_add(&ring_set, key);
                            md_ring_data_push(out, ring, len, MD_RING_FLAG_NONE, alloc);
                        }
                    }
                } else {
                    // Avoid expanding too far from the root as we are only interested in rings up to a certain size
                    // Otherwise, we will may have alot of potential large cycles such as in the case of C60 or graphene
                    depth_t d = depth[idx] + 1;
                    if (d > MAX_DEPTH) continue;

                    depth[next] = d;
                    pred[next]  = idx;
                    mark[next]  = current_mark;
                    md_fifo_push(&queue, next);
                }
            }
        }
    }

    remove_envelope_rings(out, bond, num_atoms, temp_arena);
}

void md_ring_infer_from_bonds(md_ring_data_t* out_rings, const md_bond_data_t* bond, size_t num_atoms, md_allocator_i* alloc) {
    MD_ASSERT(out_rings && alloc);
    MD_ASSERT(bond);
    md_ring_data_clear(out_rings);
    if (num_atoms == 0 || bond->count == 0 || !bond->conn.offset) return;
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    find_rings(out_rings, bond, num_atoms, NULL, alloc, temp);
    md_temp_end(temp);
}

bool md_ring_infer(md_ring_data_t* out, const md_system_t* sys, md_allocator_i* alloc) {
    MD_ASSERT(out);
    MD_ASSERT(sys);
    MD_ASSERT(alloc);
    MD_ASSERT(out != &sys->ring);
    md_ring_data_clear(out);

    size_t num_atoms = md_system_atom_count(sys);
    size_t num_bonds = md_system_bond_count(sys);

    if (num_atoms == 0) {
        return false;
    }
    if (num_bonds == 0) {
        return false;
    }
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    // Water and ions have no rings, and neither do coarse grained beads or virtual sites: a ring is made of atoms.
    // Martini closes a few beads into a rigid triangle where the atoms had a ring (one bead per two or three ring
    // atoms, fused rings a triangulated patch), and the elastic network and secondary structure bonds of a coarse
    // grained protein close thousands of cycles which are no rings at all. Skipped atoms are neither roots nor
    // passed through, so a ring of atoms in a molecule which also holds beads is still found.
    uint64_t* skip = make_bitfield(num_atoms, temp_arena);
    for (size_t ci = 0; ci < sys->residue.count; ++ci) {
        const md_residue_kind_t kind = md_system_residue_kind(sys, ci);
        if (kind == MD_RESIDUE_KIND_WATER || kind == MD_RESIDUE_KIND_ION) {
            const md_urange_t range = md_system_residue_atom_range(sys, ci);
            for (uint32_t i = range.beg; i < range.end && i < num_atoms; ++i) bitfield_set_bit(skip, i);
        }
    }
    for (size_t i = 0; i < num_atoms; ++i) {
        if (md_atom_particle_kind(&sys->atom, i) != MD_PARTICLE_ATOM) bitfield_set_bit(skip, i);
    }

    // A ring is a cycle, and a cycle lies within what is left of the graph when its ends are cut off, one after
    // another: an atom with a single neighbour (or none) is on no cycle, and once it is gone its neighbour may have a
    // single one left. Most molecules have no ring at all (a solvent, an alkane chain) and nothing is left of them.
    {
        uint32_t* degree = md_temp_alloc_array(temp, uint32_t, num_atoms);
        uint32_t* stack  = md_temp_alloc_array(temp, uint32_t, num_atoms);
        size_t top = 0;
        for (size_t i = 0; i < num_atoms; ++i) {
            uint32_t d = 0;
            if (!bitfield_test_bit(skip, i)) {
                md_bond_iter_t it = md_bond_iter(&sys->bond, i);
                for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
                    d += !bitfield_test_bit(skip, md_bond_iter_atom_index(&it));
                }
                if (d <= 1) stack[top++] = (uint32_t)i;
            }
            degree[i] = d;
        }
        while (top) {
            const uint32_t i = stack[--top];
            if (bitfield_test_bit(skip, i)) continue;
            bitfield_set_bit(skip, i);
            md_bond_iter_t it = md_bond_iter(&sys->bond, i);
            for (; md_bond_iter_has_next(&it); md_bond_iter_next(&it)) {
                const uint32_t j = (uint32_t)md_bond_iter_atom_index(&it);
                if (bitfield_test_bit(skip, j)) continue;
                if (--degree[j] == 1) stack[top++] = j;
            }
        }
    }

    find_rings(out, &sys->bond, num_atoms, skip, alloc, temp);

    md_temp_end(temp);

    ring_flags(out, sys);

    return true;
}

// Is ring (n atoms in order) the ring of a monosaccharide (MD_RING_FLAG_SUGAR). Returns the position in the ring of
// its anomeric C, -1 if it is not one, and the position of the ring O in out_ring_o. The anomeric C is the neighbor of
// the ring O with an O or N outside of the ring (C1 of an aldose, C2 of a ketose); with both neighbors carrying one, the
// one with more, and the first of the two in the ring when they carry as many.
static int sugar_anomeric(const md_system_t* sys, const int32_t* ring, int n, int* out_ring_o) {
    if (n != 5 && n != 6) return -1;
    int ring_o = -1;
    for (int k = 0; k < n; ++k) {
        const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, ring[k]);
        if (z == MD_Z_O) {
            if (ring_o >= 0) return -1;
            ring_o = k;
        } else if (z != MD_Z_C) {
            return -1;
        }
    }
    if (ring_o < 0) return -1;

    int substituted = 0;
    int sub[MD_RING_MAX_SIZE] = {0};    // O and N outside of the ring, on each ring C
    for (int k = 0; k < n; ++k) {
        if (k == ring_o) continue;
        md_bond_iter_t it = md_bond_iter(&sys->bond, ring[k]);
        while (md_bond_iter_has_next(&it)) {
            const md_atom_idx_t x = md_bond_iter_atom_index(&it);
            const md_bond_flags_t f = (md_bond_flags_t)md_bond_iter_bond_flags(&it);
            md_bond_iter_next(&it);
            // No double or aromatic bond on the ring atoms, within the ring or out of it
            if (md_bond_order(f) >= 2 || (f & MD_BOND_FLAG_AROMATIC)) return -1;
            bool in_ring = false;
            for (int j = 0; j < n; ++j) in_ring |= ring[j] == x;
            if (in_ring) continue;
            const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, x);
            sub[k] += z == MD_Z_O || z == MD_Z_N;
        }
        substituted += sub[k] > 0;
    }
    const int next = (ring_o + 1) % n;
    const int prev = (ring_o + n - 1) % n;
    if (substituted < 2 || (!sub[next] && !sub[prev])) return -1;
    *out_ring_o = ring_o;
    if (sub[next] != sub[prev]) return sub[next] > sub[prev] ? next : prev;
    return MD_MIN(next, prev);
}

// Puts the atoms of a sugar ring in their canonical order: the ring O first, then the anomeric C and on around the
// ring away from the O (O5 C1 C2 C3 C4 C5 of a pyranose, O4' C1' C2' C3' C4' of the furanose of a nucleotide)
static void sugar_canonical_order(int32_t* ring, int n, int ring_o, int anomeric) {
    int32_t tmp[MD_RING_MAX_SIZE];
    const int step = anomeric == (ring_o + 1) % n ? 1 : n - 1;
    for (int k = 0; k < n; ++k) tmp[k] = ring[(ring_o + k * step) % n];
    md_memcpy(ring, tmp, n * sizeof(int32_t));
}

// The flags of the rings of rings, from the bonds and atoms of sys (whose rings they are or are to be)
static void ring_flags(md_ring_data_t* rings, const md_system_t* sys) {
    MD_ASSERT(sys);
    const size_t num_rings = md_ring_data_count(rings);
    if (num_rings == 0) return;
    for (size_t r = 0; r < num_rings; ++r) rings->rings[r].flags = 0;
    if (!sys->bond.conn.offset || !sys->bond.flags) return;

    // Bonds in more than one ring make those rings fused
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    uint8_t* bond_rings = md_temp_alloc(temp, MD_MAX(sys->bond.count, 1));
    md_memset(bond_rings, 0, MD_MAX(sys->bond.count, 1));
    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t* ring = md_ring_atoms(rings, r);
        const size_t n = md_ring_size(rings, r);
        for (size_t k = 0; k < n; ++k) {
            const md_bond_idx_t b = md_bond_find(&sys->bond, ring[k], ring[(k + 1) % n]);
            if (b >= 0 && bond_rings[b] < 255) bond_rings[b] += 1;
        }
    }

    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t* ring = md_ring_atoms(rings, r);
        const size_t n = md_ring_size(rings, r);
        md_ring_flags_t flags = MD_RING_FLAG_NONE;
        bool aromatic = n > 0;
        for (size_t k = 0; k < n; ++k) {
            if (md_atom_atomic_number(&sys->atom, ring[k]) != MD_Z_C) flags |= MD_RING_FLAG_HETERO;
            const md_bond_idx_t b = md_bond_find(&sys->bond, ring[k], ring[(k + 1) % n]);
            if (b < 0) {
                aromatic = false;
                continue;
            }
            aromatic &= (sys->bond.flags[b] & MD_BOND_FLAG_AROMATIC) != 0;
            if (bond_rings[b] > 1) flags |= MD_RING_FLAG_FUSED;
        }
        if (aromatic) flags |= MD_RING_FLAG_AROMATIC;
        int ring_o = -1;
        const int anomeric = sugar_anomeric(sys, ring, (int)n, &ring_o);
        if (anomeric >= 0) {
            flags |= MD_RING_FLAG_SUGAR;
            sugar_canonical_order(md_ring_atoms_mut(rings, r), (int)n, ring_o, anomeric);
        }
        rings->rings[r].flags = (uint16_t)flags;
    }

    // Bonds of aromatic rings fuse the other rings they are in to an aromatic ring: an aromatic ring which shares a
    // bond with another is of a fused aromatic system (naphthalene), one fused only to rings which are not aromatic
    // keeps a sextet of its own (tetralin)
    uint8_t* bond_aromatic = md_temp_alloc(temp, MD_MAX(sys->bond.count, 1));
    md_memset(bond_aromatic, 0, MD_MAX(sys->bond.count, 1));
    for (size_t r = 0; r < num_rings; ++r) {
        if (!(md_ring_flags(rings, r) & MD_RING_FLAG_AROMATIC)) continue;
        const int32_t* ring = md_ring_atoms(rings, r);
        const size_t n = md_ring_size(rings, r);
        for (size_t k = 0; k < n; ++k) {
            const md_bond_idx_t b = md_bond_find(&sys->bond, ring[k], ring[(k + 1) % n]);
            if (b >= 0 && bond_aromatic[b] < 255) bond_aromatic[b] += 1;
        }
    }
    for (size_t r = 0; r < num_rings; ++r) {
        const int32_t* ring = md_ring_atoms(rings, r);
        const size_t n = md_ring_size(rings, r);
        // Aromatic rings other than this one
        const uint8_t self = (md_ring_flags(rings, r) & MD_RING_FLAG_AROMATIC) ? 1 : 0;
        for (size_t k = 0; k < n; ++k) {
            const md_bond_idx_t b = md_bond_find(&sys->bond, ring[k], ring[(k + 1) % n]);
            if (b >= 0 && bond_aromatic[b] > self) {
                rings->rings[r].flags |= MD_RING_FLAG_FUSED_AROMATIC;
                break;
            }
        }
    }
    md_temp_end(temp);
}

void md_system_infer_ring_flags(md_system_t* sys) {
    MD_ASSERT(sys);
    ring_flags(&sys->ring, sys);
}

#undef MIN_RING_SIZE
#undef MAX_RING_SIZE
#undef MAX_DEPTH
