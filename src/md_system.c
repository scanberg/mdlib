#include <md_system.h>
#include <md_trajectory.h>

#include <core/md_log.h>
#include <core/md_array.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>

#ifdef __cplusplus
extern "C" {
#endif

void md_system_free(md_system_t* sys) {
    ASSERT(sys);
    ASSERT(sys->alloc);
    md_allocator_i* alloc = sys->alloc;
    md_trajectory_free(sys->trajectory);

    // ATOM
    md_array_free(sys->atom.x, alloc);
    md_array_free(sys->atom.y, alloc);
    md_array_free(sys->atom.z, alloc);
    md_array_free(sys->atom.type_idx, alloc);
    md_array_free(sys->atom.flags, alloc);

    // ATOM TYPE
    md_array_free(sys->atom.type.name, alloc);
    md_array_free(sys->atom.type.z, alloc);
    md_array_free(sys->atom.type.mass, alloc);
    md_array_free(sys->atom.type.radius, alloc);
    md_array_free(sys->atom.type.color, alloc);
    md_array_free(sys->atom.type.flags, alloc);

    // COMPONENT
    md_array_free(sys->component.name, alloc);
    md_array_free(sys->component.seq_id, alloc);
    md_array_free(sys->component.atom_offset, alloc);
    md_array_free(sys->component.flags, alloc);

    // INSTANCE
    md_array_free(sys->instance.id, alloc);
    md_array_free(sys->instance.auth_id, alloc);
    md_array_free(sys->instance.comp_offset, alloc);
    md_array_free(sys->instance.entity_idx, alloc);

    // ENTITY
    md_array_free(sys->entity.id, alloc);
    md_array_free(sys->entity.flags, alloc);
    for (size_t i = 0; i < sys->entity.count; ++i) {
        if (!str_empty(sys->entity.description[i])) {
            str_free(sys->entity.description[i], alloc);
        }
    }

    // PROTEIN BACKBONE
    md_array_free(sys->protein_backbone.range.offset, alloc);
    md_array_free(sys->protein_backbone.range.inst_idx, alloc);
    md_array_free(sys->protein_backbone.segment.atoms, alloc);
    md_array_free(sys->protein_backbone.segment.angle, alloc);
    md_array_free(sys->protein_backbone.segment.secondary_structure, alloc);
    md_array_free(sys->protein_backbone.segment.rama_type, alloc);
    md_array_free(sys->protein_backbone.segment.comp_idx, alloc);

    // NUCLEIC BACKBONE
    md_array_free(sys->nucleic_backbone.range.offset, alloc);
    md_array_free(sys->nucleic_backbone.range.inst_idx, alloc);
    md_array_free(sys->nucleic_backbone.segment.atoms, alloc);
    md_array_free(sys->nucleic_backbone.segment.comp_idx, alloc);

    // BONDS
    md_array_free(sys->bond.pairs, alloc);
    md_array_free(sys->bond.flags, alloc);
    md_array_free(sys->bond.conn.atom_idx, alloc);
    md_array_free(sys->bond.conn.bond_idx, alloc);
    md_array_free(sys->bond.conn.offset, alloc);

    // HYDROGEN BONDS
    md_array_free(sys->hydrogen_bond.candidate.acceptor.idx, alloc);
    md_array_free(sys->hydrogen_bond.candidate.acceptor.num_lone_pairs, alloc);
    md_array_free(sys->hydrogen_bond.candidate.donor.d_idx, alloc);
    md_array_free(sys->hydrogen_bond.candidate.donor.h_idx, alloc);

    md_index_data_free(&sys->ring);
    md_index_data_free(&sys->structure);

    // ASSEMBLY
    md_array_free(sys->assembly.atom_range, alloc);
    md_array_free(sys->assembly.label, alloc);
    md_array_free(sys->assembly.transform, alloc);

    if (!str_empty(sys->description)) {
        str_free(sys->description, alloc);
    }

    MEMSET(sys, 0, sizeof(md_system_t));
}

void md_system_reset(md_system_t* sys) {
    ASSERT(sys);
    md_allocator_i* alloc = sys->alloc;
    md_system_free(sys);
    sys->alloc = alloc;
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

bool md_system_copy(md_system_t* dst, const md_system_t* src) {
    ASSERT(dst);
    ASSERT(src);

    if (!dst->alloc) {
        MD_LOG_ERROR("System allocator not set");
        return false;
    }

    md_system_reset(dst);

    md_allocator_i* alloc = dst->alloc;

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

    ARRAY_PUSH(instance, id);
    ARRAY_PUSH(instance, auth_id);
    ARRAY_PUSH(instance, comp_offset);

    md_array_push_array(dst->bond.pairs, src->bond.pairs, src->bond.count, alloc);
    md_array_push_array(dst->bond.flags, src->bond.flags, src->bond.count, alloc);
    md_array_push_array(dst->bond.conn.atom_idx, src->bond.conn.atom_idx, src->bond.conn.count, alloc);
    md_array_push_array(dst->bond.conn.bond_idx, src->bond.conn.bond_idx, src->bond.conn.count, alloc);
    md_array_push_array(dst->bond.conn.offset, src->bond.conn.offset, src->bond.conn.offset_count, alloc);

    ARRAY_PUSH(component, name);
    ARRAY_PUSH(component, seq_id);
    ARRAY_PUSH(component, atom_offset);
    ARRAY_PUSH(component, flags);

    dst->atom.count           = src->atom.count;
    dst->atom.type.count      = src->atom.type.count;
    dst->protein_backbone.segment.count = src->protein_backbone.segment.count;
    dst->protein_backbone.range.count   = src->protein_backbone.range.count;
    dst->nucleic_backbone.segment.count = src->nucleic_backbone.segment.count;
    dst->nucleic_backbone.range.count   = src->nucleic_backbone.range.count;
    dst->instance.count       = src->instance.count;
    dst->component.count      = src->component.count;
    dst->bond.count           = src->bond.count;
    dst->bond.conn.count      = src->bond.conn.count;
    dst->bond.conn.offset_count = src->bond.conn.offset_count;
    dst->initial_unitcell     = src->initial_unitcell;
    dst->alloc                = alloc;
    dst->trajectory           = NULL;

    return true;
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

// Attach a trajectory to the system, freeing any existing attached trajectory.
void md_system_attach_trajectory(md_system_t* sys, struct md_trajectory_i* traj) {
    if (!sys) return;
    md_trajectory_free(sys->trajectory);
    sys->trajectory = traj;
}

#ifdef __cplusplus
}
#endif
