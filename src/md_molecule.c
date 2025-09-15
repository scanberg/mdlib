#include <md_molecule.h>

#include <core/md_array.h>
#include <core/md_allocator.h>

#ifdef __cplusplus
extern "C" {
#endif

void md_atom_type_data_fill(md_atom_type_data_t* atom_types, const str_t atom_names[], const md_atomic_number_t atom_element[], const float atom_mass[], const float atom_radii[], size_t atom_count, md_allocator_i* alloc) {
    
}

    void md_molecule_free(md_molecule_t* mol, struct md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(alloc);

    // Atom
    if (mol->atom.x) md_array_free(mol->atom.x, alloc);
    if (mol->atom.y) md_array_free(mol->atom.y, alloc);
    if (mol->atom.z) md_array_free(mol->atom.z, alloc);
    if (mol->atom.type_idx) md_array_free(mol->atom.type_idx, alloc);
//    if (mol->atom.radius) md_array_free(mol->atom.radius, alloc);
//    if (mol->atom.mass) md_array_free(mol->atom.mass, alloc);
//    if (mol->atom.element) md_array_free(mol->atom.element, alloc);
//    if (mol->atom.type) md_array_free(mol->atom.type, alloc);
//    if (mol->atom.resid) md_array_free(mol->atom.resid, alloc);
//    if (mol->atom.resname) md_array_free(mol->atom.resname, alloc);
//    if (mol->atom.chainid) md_array_free(mol->atom.chainid, alloc);
//    if (mol->atom.res_idx) md_array_free(mol->atom.res_idx, alloc);
//    if (mol->atom.chain_idx) md_array_free(mol->atom.chain_idx, alloc);
    if (mol->atom.flags) md_array_free(mol->atom.flags, alloc);

    // Atom Type
    if (mol->atom.type_data.name)   md_array_free(mol->atom.type_data.name, alloc);
    if (mol->atom.type_data.z)      md_array_free(mol->atom.type_data.z, alloc);
    if (mol->atom.type_data.mass)   md_array_free(mol->atom.type_data.mass, alloc);
    if (mol->atom.type_data.radius) md_array_free(mol->atom.type_data.radius, alloc);

    // Residue
    if (mol->residue.name) md_array_free(mol->residue.name, alloc);
    if (mol->residue.id) md_array_free(mol->residue.id, alloc);
    if (mol->residue.atom_offset) md_array_free(mol->residue.atom_offset, alloc);

    // Chain
    if (mol->chain.id) md_array_free(mol->chain.id, alloc);
    if (mol->chain.res_range)  md_array_free(mol->chain.res_range, alloc);
    if (mol->chain.atom_range) md_array_free(mol->chain.atom_range, alloc);

    // Backbone
    if (mol->protein_backbone.range.offset) md_array_free(mol->protein_backbone.range.offset, alloc);
    if (mol->protein_backbone.range.chain_idx) md_array_free(mol->protein_backbone.range.chain_idx, alloc);
    if (mol->protein_backbone.atoms) md_array_free(mol->protein_backbone.atoms, alloc);
    if (mol->protein_backbone.angle) md_array_free(mol->protein_backbone.angle, alloc);
    if (mol->protein_backbone.secondary_structure) md_array_free(mol->protein_backbone.secondary_structure, alloc);
    if (mol->protein_backbone.ramachandran_type) md_array_free(mol->protein_backbone.ramachandran_type, alloc);
    if (mol->protein_backbone.residue_idx) md_array_free(mol->protein_backbone.residue_idx, alloc);

    // Nucleic Backbone
    if (mol->nucleic_backbone.range.offset) md_array_free(mol->nucleic_backbone.range.offset, alloc);
    if (mol->nucleic_backbone.range.chain_idx) md_array_free(mol->nucleic_backbone.range.chain_idx, alloc);
    if (mol->nucleic_backbone.atoms) md_array_free(mol->nucleic_backbone.atoms, alloc);
    if (mol->nucleic_backbone.residue_idx) md_array_free(mol->nucleic_backbone.residue_idx, alloc);

    // Bonds
    if (mol->bond.pairs) md_array_free(mol->bond.pairs, alloc);
    if (mol->bond.order) md_array_free(mol->bond.order, alloc);

    if (mol->bond.conn.atom_idx) md_array_free(mol->bond.conn.atom_idx, alloc);
    if (mol->bond.conn.bond_idx) md_array_free(mol->bond.conn.bond_idx, alloc);
    if (mol->bond.conn.offset) md_array_free(mol->bond.conn.offset, alloc);

    md_index_data_free(&mol->structure);
    md_index_data_free(&mol->ring);

    // Instance
    if (mol->instance.atom_range) md_array_free(mol->instance.atom_range, alloc);
    if (mol->instance.label) md_array_free(mol->instance.label, alloc);
    if (mol->instance.transform) md_array_free(mol->instance.transform, alloc);

    MEMSET(mol, 0, sizeof(md_molecule_t));
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

void md_molecule_copy(md_molecule_t* dst, const md_molecule_t* src, struct md_allocator_i* alloc) {
    // Make sure that the dst_molecule is clean!
    MEMSET(dst, 0, sizeof(md_molecule_t));

    ARRAY_PUSH(atom, x);
    ARRAY_PUSH(atom, y);
    ARRAY_PUSH(atom, z);
    ARRAY_PUSH(atom, type_idx);
    //ARRAY_PUSH(atom, radius);
    //ARRAY_PUSH(atom, mass);
    //ARRAY_PUSH(atom, element);
    //ARRAY_PUSH(atom, type);
    ARRAY_PUSH(atom, flags);
    //ARRAY_PUSH(atom, resid);
    //ARRAY_PUSH(atom, resname);
    //ARRAY_PUSH(atom, chainid);
    //ARRAY_PUSH(atom, res_idx);
    //ARRAY_PUSH(atom, chain_idx);

    ARRAY_PUSH(atom, type_data.name);
    ARRAY_PUSH(atom, type_data.z);
    ARRAY_PUSH(atom, type_data.mass);
    ARRAY_PUSH(atom, type_data.radius);

    ARRAY_PUSH(protein_backbone, atoms);
    ARRAY_PUSH(protein_backbone, angle);
    ARRAY_PUSH(protein_backbone, secondary_structure);
    ARRAY_PUSH(protein_backbone, ramachandran_type);
    ARRAY_PUSH(protein_backbone, residue_idx);

    md_array_push_array(dst->protein_backbone.range.offset, src->protein_backbone.range.offset, src->protein_backbone.range.count, alloc);
    md_array_push_array(dst->protein_backbone.range.chain_idx, src->protein_backbone.range.chain_idx, src->protein_backbone.range.count, alloc);

    ARRAY_PUSH(nucleic_backbone, atoms);
    ARRAY_PUSH(nucleic_backbone, residue_idx);

    md_array_push_array(dst->nucleic_backbone.range.offset, src->nucleic_backbone.range.offset, src->nucleic_backbone.range.count, alloc);
    md_array_push_array(dst->nucleic_backbone.range.chain_idx, src->nucleic_backbone.range.chain_idx, src->nucleic_backbone.range.count, alloc);

    ARRAY_PUSH(chain, id);
    ARRAY_PUSH(chain, res_range);
    ARRAY_PUSH(chain, atom_range);

    md_array_push_array(dst->bond.pairs, src->bond.pairs, src->bond.count, alloc);
    md_array_push_array(dst->bond.order, src->bond.order, src->bond.count, alloc);
    md_array_push_array(dst->bond.conn.atom_idx, src->bond.conn.atom_idx, src->bond.conn.count, alloc);
    md_array_push_array(dst->bond.conn.bond_idx, src->bond.conn.bond_idx, src->bond.conn.count, alloc);
    md_array_push_array(dst->bond.conn.offset, src->bond.conn.offset, src->bond.conn.offset_count, alloc);

    ARRAY_PUSH(residue, name);
    ARRAY_PUSH(residue, id);
    ARRAY_PUSH(residue, atom_offset);
    ARRAY_PUSH(residue, flags);

    dst->atom.count           = src->atom.count;
    dst->atom.type_data.count      = src->atom.type_data.count;
    dst->protein_backbone.count       = src->protein_backbone.count;
    dst->protein_backbone.range.count = src->protein_backbone.range.count;
    dst->nucleic_backbone.count       = src->nucleic_backbone.count;
    dst->nucleic_backbone.range.count = src->nucleic_backbone.range.count;
    dst->chain.count          = src->chain.count;
    dst->residue.count        = src->residue.count;
    dst->bond.count           = src->bond.count;
    dst->bond.conn.count      = src->bond.conn.count;
    dst->bond.conn.offset_count = src->bond.conn.offset_count;
    dst->unit_cell            = src->unit_cell;
}

#undef ARRAY_PUSH
#undef ARRAY_INCREMENT
#undef ARRAY_INCREMENT_FIELD

#ifdef __cplusplus
}
#endif
