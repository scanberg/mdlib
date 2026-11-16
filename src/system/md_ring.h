#pragma once

#include <stdint.h>
#include <stddef.h>

#include <core/md_common.h>
#include <core/md_array.h>

struct md_allocator_i;

// Rings of a system (md_system_t::ring): 3 to 7 atoms in order around the ring, and what the ring is.

// The sizes of the rings of md_system_t::ring (md_ring_infer). Larger cycles are loops, not rings.
#define MD_RING_MIN_SIZE 3
#define MD_RING_MAX_SIZE 7

// Ring flags (md_ring_t::flags), what a ring of the system is. md_system_infer_ring_flags derives them
// from the rings, the bonds and their flags; md_ring_infer and md_chem_perceive call it.
//   AROMATIC  Every bond of the ring is aromatic (md_chem_perceive, or the bonds a file gives): an aromatic ring, alone
//             or fused with a neighbor.
//   HETERO    The ring has an atom other than carbon.
//   FUSED     The ring shares a bond with another ring.
//   SUGAR     The ring of a monosaccharide, a furanose or a pyranose: 5 or 6 atoms, one of them O and the rest C, no
//             bond of a known order above single or aromatic in or out of the ring, and at least two of the C carrying
//             an O or N outside of the ring (a hydroxyl, a glycosidic link, the N of a nucleoside), one of them next to
//             the ring O (the anomeric C). Deoxy sugars (the 2-deoxyribose of DNA) count; tetrahydrofuran does not.
//             The atoms of a sugar ring are in canonical order in md_system_t::ring: the ring O, the anomeric C, and on
//             around the ring away from the O (O5 C1 C2 C3 C4 C5 of a pyranose, O4' C1' C2' C3' C4' of a furanose).
//   FUSED_AROMATIC
//             The ring shares a bond with an aromatic ring. With AROMATIC: a ring of a fused aromatic system, whose
//             sextets are not assigned to particular rings (naphthalene, indole, a purine, azulene). An aromatic ring
//             without it holds a sextet of its own, alone or fused to rings which are not aromatic (benzene, the
//             benzene ring of tetralin or of estradiol, those of biphenylene).
typedef enum md_ring_flags_t {
    MD_RING_FLAG_NONE           = 0,
    MD_RING_FLAG_AROMATIC       = 1u << 0,
    MD_RING_FLAG_HETERO         = 1u << 1,
    MD_RING_FLAG_FUSED          = 1u << 2,
    MD_RING_FLAG_SUGAR          = 1u << 3,
    MD_RING_FLAG_FUSED_AROMATIC = 1u << 4,
} md_ring_flags_t;

MD_ENUM_FLAGS(md_ring_flags_t)

// A ring: MD_RING_MIN_SIZE to MD_RING_MAX_SIZE atoms in order around the ring, and its classification flags.
typedef struct md_ring_t {
    int32_t  atoms[MD_RING_MAX_SIZE];   // Only the first count are valid
    uint16_t count;
    uint16_t flags;                     // md_ring_flags_t bits
} md_ring_t;

MD_STATIC_ASSERT(sizeof(md_ring_t) == 32, "md_ring_t should be 32 bytes");

// The rings of a system.
typedef struct md_ring_data_t {
    md_array(md_ring_t) rings;
} md_ring_data_t;

static inline void md_ring_data_free(md_ring_data_t* data, struct md_allocator_i* alloc) {
    MD_ASSERT(data);
    if (data->rings) {
        MD_ASSERT(alloc);
        md_array_free(data->rings, alloc);
    }
#if MD_DEBUG
    md_memset(data, 0, sizeof(md_ring_data_t));
#endif
}

static inline void md_ring_data_clear(md_ring_data_t* data) {
    MD_ASSERT(data);
    md_array_shrink(data->rings, 0);
}

static inline size_t md_ring_data_count(const md_ring_data_t* data) {
    MD_ASSERT(data);
    return md_array_size(data->rings);
}

// Pushes a ring of count atoms (MD_RING_MIN_SIZE to MD_RING_MAX_SIZE) and returns its index
static inline size_t md_ring_data_push(md_ring_data_t* data, const int32_t* atoms, size_t count, md_ring_flags_t flags, struct md_allocator_i* alloc) {
    MD_ASSERT(data);
    MD_ASSERT(alloc);
    MD_ASSERT(atoms);
    MD_ASSERT(MD_RING_MIN_SIZE <= count && count <= MD_RING_MAX_SIZE);
    md_ring_t ring = {0};
    md_memcpy(ring.atoms, atoms, count * sizeof(int32_t));
    ring.count = (uint16_t)count;
    ring.flags = (uint16_t)flags;
    md_array_push(data->rings, ring, alloc);
    return md_array_size(data->rings) - 1;
}

// The atoms of a ring in order, NULL if the ring does not exist
static inline const int32_t* md_ring_atoms(const md_ring_data_t* data, size_t ring_idx) {
    MD_ASSERT(data);
    return ring_idx < md_array_size(data->rings) ? data->rings[ring_idx].atoms : NULL;
}

static inline int32_t* md_ring_atoms_mut(md_ring_data_t* data, size_t ring_idx) {
    MD_ASSERT(data);
    return ring_idx < md_array_size(data->rings) ? data->rings[ring_idx].atoms : NULL;
}

static inline size_t md_ring_size(const md_ring_data_t* data, size_t ring_idx) {
    MD_ASSERT(data);
    return ring_idx < md_array_size(data->rings) ? data->rings[ring_idx].count : 0;
}

static inline md_ring_flags_t md_ring_flags(const md_ring_data_t* data, size_t ring_idx) {
    MD_ASSERT(data);
    return ring_idx < md_array_size(data->rings) ? (md_ring_flags_t)data->rings[ring_idx].flags : MD_RING_FLAG_NONE;
}

struct md_system_t;
struct md_bond_data_t;

#ifdef __cplusplus
extern "C" {
#endif

// ### INFERENCE ###

// The rings of the system formed by its covalent bonds, into out (allocated with alloc, never &sys->ring), with their
// flags: the rings of 3 to 7 atoms which are not sums of smaller rings (the smallest set of smallest rings, made
// symmetric). Only atoms (MD_PARTICLE_ATOM) make rings: the cycles of coarse grained beads and virtual sites are not
// searched, nor water and ions. False, and out empty, when the system has no atoms or no bonds.
bool md_ring_infer(md_ring_data_t* out, const struct md_system_t* sys, struct md_allocator_i* alloc);

// The rings of a bond graph as md_ring_infer finds them, of every atom (nothing skipped), without flags: bond, with its
// connectivity (md_bond_build_connectivity), over num_atoms atoms. Into out, cleared first and allocated with alloc.
void md_ring_infer_from_bonds(md_ring_data_t* out, const struct md_bond_data_t* bond, size_t num_atoms, struct md_allocator_i* alloc);

// The flags of the rings of the system (md_ring_t::flags, md_ring_flags_t) anew, from the rings, the bonds and their
// flags, and the atoms of sugar rings in canonical order (MD_RING_FLAG_SUGAR). Again when the bond orders or aromaticity
// change (md_chem_perceive does).
void md_system_infer_ring_flags(struct md_system_t* sys);

#ifdef __cplusplus
}
#endif
