#include <md_system.h>

#include <core/md_array.h>
#include <core/md_allocator.h>

#ifdef __cplusplus
extern "C" {
#endif

#define SAFE_FREE(P)          \
    do {                         \
        if (P) {                 \
            md_array_free(P, alloc); \
            P = NULL;            \
        }                        \
    } while (0)

    void md_system_free(md_system_t* mol, struct md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(alloc);

    // Atom
    SAFE_FREE(mol->atom.x);
    SAFE_FREE(mol->atom.y);
    SAFE_FREE(mol->atom.z);
    SAFE_FREE(mol->atom.type_idx);
    SAFE_FREE(mol->atom.flags);

    // Atom Type
    SAFE_FREE(mol->atom.type.name);
    SAFE_FREE(mol->atom.type.z);
    SAFE_FREE(mol->atom.type.mass);
    SAFE_FREE(mol->atom.type.radius);
    SAFE_FREE(mol->atom.type.color);

    // Component
    SAFE_FREE(mol->comp.name);
    SAFE_FREE(mol->comp.seq_id);
    SAFE_FREE(mol->comp.atom_offset);
    SAFE_FREE(mol->comp.flags);

    // Instance
    SAFE_FREE(mol->inst.id);
    SAFE_FREE(mol->inst.auth_id);
    SAFE_FREE(mol->inst.entity_idx);
    SAFE_FREE(mol->inst.comp_offset);

    // Protein Backbone
    SAFE_FREE(mol->protein_backbone.range.offset);
    SAFE_FREE(mol->protein_backbone.range.inst_idx);
    SAFE_FREE(mol->protein_backbone.segment.atoms);
    SAFE_FREE(mol->protein_backbone.segment.angle);
    SAFE_FREE(mol->protein_backbone.segment.secondary_structure);
    SAFE_FREE(mol->protein_backbone.segment.rama_type);
    SAFE_FREE(mol->protein_backbone.segment.comp_idx);

    // Nucleic Backbone
    SAFE_FREE(mol->nucleic_backbone.range.offset);
    SAFE_FREE(mol->nucleic_backbone.range.inst_idx);
    SAFE_FREE(mol->nucleic_backbone.segment.atoms);
    SAFE_FREE(mol->nucleic_backbone.segment.comp_idx);

    // Bonds
    SAFE_FREE(mol->bond.pairs);
    SAFE_FREE(mol->bond.flags);
    
    SAFE_FREE(mol->bond.conn.atom_idx);
    SAFE_FREE(mol->bond.conn.bond_idx);
    SAFE_FREE(mol->bond.conn.offset);

    md_index_data_free(&mol->structure);
    md_index_data_free(&mol->ring);

    // Assembly
    SAFE_FREE(mol->assembly.atom_range);
    SAFE_FREE(mol->assembly.label);
    SAFE_FREE(mol->assembly.transform);

    MEMSET(mol, 0, sizeof(md_system_t));
}

// This is from Mathis work in Openspace modified slightly
// Some macros to help get the job done. The idea is to clone each field in src to dst,
// but some index fields need also to be incremented to point to the cloned data.

#define ARRAY_PUSH(A, B) \
        if (src->A.B) md_array_push_array(dst->A.B, src->A.B, src->A.count, alloc)

#define ARRAY_INCREMENT(A, B, C) \
        for (int64_t j = 0; j < src->A.count; j++) \
            if (src->A.B) dst->A.B[dst->A.count + j] += C

#define ARRAY_INCREMENT_FIELD(A, B, C, D) \
        for (int64_t j = 0; j < src->A.count; j++) \
            if (src->A.B) dst->A.B[dst->A.count + j].C += D

void md_system_copy(md_system_t* dst, const md_system_t* src, struct md_allocator_i* alloc) {
    // Make sure that the dst_molecule is clean!
    MEMSET(dst, 0, sizeof(md_system_t));

    ARRAY_PUSH(atom, x);
    ARRAY_PUSH(atom, y);
    ARRAY_PUSH(atom, z);
    ARRAY_PUSH(atom, type_idx);
    ARRAY_PUSH(atom, flags);

    ARRAY_PUSH(atom, type.name);
    ARRAY_PUSH(atom, type.z);
    ARRAY_PUSH(atom, type.mass);
    ARRAY_PUSH(atom, type.radius);

    ARRAY_PUSH(protein_backbone.segment, atoms);
    ARRAY_PUSH(protein_backbone.segment, angle);
    ARRAY_PUSH(protein_backbone.segment, secondary_structure);
    ARRAY_PUSH(protein_backbone.segment, rama_type);
    ARRAY_PUSH(protein_backbone.segment, comp_idx);

    md_array_push_array(dst->protein_backbone.range.offset, src->protein_backbone.range.offset, src->protein_backbone.range.count, alloc);
    md_array_push_array(dst->protein_backbone.range.inst_idx, src->protein_backbone.range.inst_idx, src->protein_backbone.range.count, alloc);

    ARRAY_PUSH(nucleic_backbone.segment, atoms);
    ARRAY_PUSH(nucleic_backbone.segment, comp_idx);

    md_array_push_array(dst->nucleic_backbone.range.offset, src->nucleic_backbone.range.offset, src->nucleic_backbone.range.count, alloc);
    md_array_push_array(dst->nucleic_backbone.range.inst_idx, src->nucleic_backbone.range.inst_idx, src->nucleic_backbone.range.count, alloc);

    ARRAY_PUSH(inst, id);
    ARRAY_PUSH(inst, auth_id);
    ARRAY_PUSH(inst, comp_offset);

    md_array_push_array(dst->bond.pairs, src->bond.pairs, src->bond.count, alloc);
    md_array_push_array(dst->bond.flags, src->bond.flags, src->bond.count, alloc);
    md_array_push_array(dst->bond.conn.atom_idx, src->bond.conn.atom_idx, src->bond.conn.count, alloc);
    md_array_push_array(dst->bond.conn.bond_idx, src->bond.conn.bond_idx, src->bond.conn.count, alloc);
    md_array_push_array(dst->bond.conn.offset, src->bond.conn.offset, src->bond.conn.offset_count, alloc);

    ARRAY_PUSH(comp, name);
    ARRAY_PUSH(comp, seq_id);
    ARRAY_PUSH(comp, atom_offset);
    ARRAY_PUSH(comp, flags);

    dst->atom.count           = src->atom.count;
    dst->atom.type.count      = src->atom.type.count;
    dst->protein_backbone.segment.count = src->protein_backbone.segment.count;
    dst->protein_backbone.range.count   = src->protein_backbone.range.count;
    dst->nucleic_backbone.segment.count = src->nucleic_backbone.segment.count;
    dst->nucleic_backbone.range.count   = src->nucleic_backbone.range.count;
    dst->inst.count           = src->inst.count;
    dst->comp.count           = src->comp.count;
    dst->bond.count           = src->bond.count;
    dst->bond.conn.count      = src->bond.conn.count;
    dst->bond.conn.offset_count = src->bond.conn.offset_count;
    dst->unitcell             = src->unitcell;
}

#undef ARRAY_PUSH
#undef ARRAY_INCREMENT
#undef ARRAY_INCREMENT_FIELD

static void build_connectivity(md_bond_conn_data_t* conn, const md_atom_pair_t* bond_pairs, size_t bond_pair_count, size_t atom_count, md_allocator_i* alloc) {
    ASSERT(conn);
    ASSERT(alloc);

    if (bond_pairs == NULL) return;
    if (bond_pair_count == 0) return;
    if (atom_count == 0) return;

    conn->offset_count = atom_count + 1;
    md_array_resize(conn->offset, conn->offset_count, alloc);
    MEMSET(conn->offset, 0, md_array_bytes(conn->offset));

    // This have length of 2 * bond_count (one for each direction of the bond)
    conn->count = 2 * bond_pair_count;
    md_array_resize(conn->atom_idx, conn->count, alloc);
    md_array_resize(conn->bond_idx, conn->count, alloc);

    typedef struct {
        uint16_t off[2];
    } offset_t;

    offset_t* local_offset = calloc(bond_pair_count, sizeof(offset_t));
    ASSERT(local_offset);

    // Two packed 16-bit local offsets for each of the bond idx
    // Use offsets as accumulators for length
    for (size_t i = 0; i < bond_pair_count; ++i) {
		ASSERT(bond_pairs[i].idx[0] < (md_atom_idx_t)atom_count);
		ASSERT(bond_pairs[i].idx[1] < (md_atom_idx_t)atom_count);
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
        const md_atom_pair_t p = bond_pairs[i];
        const int atom_a = p.idx[0];
        const int atom_b = p.idx[1];
        const int local_a = (int)local_offset[i].off[0];
        const int local_b = (int)local_offset[i].off[1];
        const int off_a = conn->offset[atom_a];
        const int off_b = conn->offset[atom_b];

        const int idx_a = off_a + local_a;
        const int idx_b = off_b + local_b;

        ASSERT(idx_a < (int)conn->count);
        ASSERT(idx_b < (int)conn->count);

        // Store the cross references to the 'other' atom index signified by the bond in the correct location
        conn->atom_idx[idx_a] = atom_b;
        conn->atom_idx[idx_b] = atom_a;

        conn->bond_idx[idx_a] = (md_bond_idx_t)i;
        conn->bond_idx[idx_b] = (md_bond_idx_t)i;
    }

    free(local_offset);
}

void md_bond_build_connectivity(md_bond_data_t* in_out_bond, size_t atom_count, md_allocator_i* alloc) {
    ASSERT(in_out_bond);
    ASSERT(alloc);
	build_connectivity(&in_out_bond->conn, in_out_bond->pairs, in_out_bond->count, atom_count, alloc);
}

void md_system_bond_build_connectivity(md_system_t* sys, md_allocator_i* alloc) {
    ASSERT(sys);
    ASSERT(alloc);
	md_bond_build_connectivity(&sys->bond, sys->atom.count, alloc);
}

#ifdef __cplusplus
}
#endif
