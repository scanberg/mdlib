#include <md_molecule.h>

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
    SAFE_FREE(mol->bond.order);
    
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
    md_array_push_array(dst->bond.order, src->bond.order, src->bond.count, alloc);
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

#ifdef __cplusplus
}
#endif
