#include <md_molecule.h>

#include <core/md_array.h>
#include <core/md_allocator.h>

#ifdef __cplusplus
extern "C" {
#endif

void md_molecule_free(md_molecule_t* mol, struct md_allocator_i* alloc) {
    ASSERT(mol);
    ASSERT(alloc);

    // Atom
    if (mol->atom.x) md_array_free(mol->atom.x, alloc);
    if (mol->atom.y) md_array_free(mol->atom.y, alloc);
    if (mol->atom.z) md_array_free(mol->atom.z, alloc);
    if (mol->atom.radius) md_array_free(mol->atom.radius, alloc);
    if (mol->atom.mass) md_array_free(mol->atom.mass, alloc);
    if (mol->atom.element) md_array_free(mol->atom.element, alloc);
    if (mol->atom.resid) md_array_free(mol->atom.resid, alloc);
    if (mol->atom.resname) md_array_free(mol->atom.resname, alloc);
    if (mol->atom.chainid) md_array_free(mol->atom.chainid, alloc);
    if (mol->atom.flags) md_array_free(mol->atom.flags, alloc);

    // Residue
    if (mol->residue.name) md_array_free(mol->residue.name, alloc);
    if (mol->residue.id) md_array_free(mol->residue.id, alloc);
    if (mol->residue.atom_offset) md_array_free(mol->residue.atom_offset, alloc);

    // Chain
    if (mol->chain.id) md_array_free(mol->chain.id, alloc);
    if (mol->chain.res_offset) md_array_free(mol->chain.res_offset, alloc);
    if (mol->chain.atom_offset) md_array_free(mol->chain.atom_offset, alloc);

    // Backbone
    if (mol->protein_backbone.range.offset) md_array_free(mol->protein_backbone.range.offset, alloc);
    if (mol->protein_backbone.atoms) md_array_free(mol->protein_backbone.atoms, alloc);
    if (mol->protein_backbone.angle) md_array_free(mol->protein_backbone.angle, alloc);
    if (mol->protein_backbone.secondary_structure) md_array_free(mol->protein_backbone.secondary_structure, alloc);
    if (mol->protein_backbone.ramachandran_type) md_array_free(mol->protein_backbone.ramachandran_type, alloc);
    if (mol->protein_backbone.residue_idx) md_array_free(mol->protein_backbone.residue_idx, alloc);

    // Bonds
    if (mol->bond.pairs) md_array_free(mol->bond.pairs, alloc);
    if (mol->bond.order) md_array_free(mol->bond.order, alloc);

    if (mol->bond.conn.atom_idx) md_array_free(mol->bond.conn.atom_idx, alloc);
    if (mol->bond.conn.bond_idx) md_array_free(mol->bond.conn.bond_idx, alloc);

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
    ARRAY_PUSH(atom, radius);
    ARRAY_PUSH(atom, mass);
    ARRAY_PUSH(atom, element);
    ARRAY_PUSH(atom, type);
    ARRAY_PUSH(atom, flags);
    ARRAY_PUSH(atom, resid);
    ARRAY_PUSH(atom, resname);
    ARRAY_PUSH(atom, chainid);

    ARRAY_PUSH(protein_backbone, atoms);
    ARRAY_PUSH(protein_backbone, angle);
    ARRAY_PUSH(protein_backbone, secondary_structure);
    ARRAY_PUSH(protein_backbone, ramachandran_type);
    ARRAY_PUSH(protein_backbone, residue_idx);

    md_array_push_array(dst->protein_backbone.range.offset, src->protein_backbone.range.offset, src->protein_backbone.range.count, alloc);

    ARRAY_PUSH(chain, id);
    ARRAY_PUSH(chain, res_offset);
    ARRAY_PUSH(chain, atom_offset);

    md_array_push_array(dst->bond.pairs, src->bond.pairs, src->bond.count, alloc);

    ARRAY_PUSH(residue, name);
    ARRAY_PUSH(residue, id);
    ARRAY_PUSH(residue, atom_offset);

    dst->atom.count           = src->atom.count;
    dst->protein_backbone.count       = src->protein_backbone.count;
    dst->protein_backbone.range.count = src->protein_backbone.range.count;
    dst->chain.count          = src->chain.count;
    dst->residue.count        = src->residue.count;
}

#undef ARRAY_PUSH
#undef ARRAY_INCREMENT
#undef ARRAY_INCREMENT_FIELD

#ifdef __cplusplus
}
#endif
