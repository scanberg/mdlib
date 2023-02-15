#include "md_molecule.h"
#include "core/md_array.h"
#include "core/md_allocator.h"

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
    if (mol->atom.vx) md_array_free(mol->atom.vx, alloc);
    if (mol->atom.vy) md_array_free(mol->atom.vy, alloc);
    if (mol->atom.vz) md_array_free(mol->atom.vz, alloc);
    if (mol->atom.radius) md_array_free(mol->atom.radius, alloc);
    if (mol->atom.mass) md_array_free(mol->atom.mass, alloc);
    if (mol->atom.valence) md_array_free(mol->atom.valence, alloc);
    if (mol->atom.element) md_array_free(mol->atom.element, alloc);
    if (mol->atom.flags) md_array_free(mol->atom.flags, alloc);
    if (mol->atom.residue_idx) md_array_free(mol->atom.residue_idx, alloc);
    if (mol->atom.chain_idx) md_array_free(mol->atom.chain_idx, alloc);

    // Residue
    if (mol->residue.name) md_array_free(mol->residue.name, alloc);
    if (mol->residue.id) md_array_free(mol->residue.id, alloc);
    if (mol->residue.atom_range) md_array_free(mol->residue.atom_range, alloc);

    // Chain
    if (mol->chain.id) md_array_free(mol->chain.id, alloc);
    if (mol->chain.residue_range) md_array_free(mol->chain.residue_range, alloc);
    if (mol->chain.atom_range) md_array_free(mol->chain.atom_range, alloc);

    // Backbone
    if (mol->backbone.range) md_array_free(mol->backbone.range, alloc);
    if (mol->backbone.atoms) md_array_free(mol->backbone.atoms, alloc);
    if (mol->backbone.angle) md_array_free(mol->backbone.angle, alloc);
    if (mol->backbone.secondary_structure) md_array_free(mol->backbone.secondary_structure, alloc);
    if (mol->backbone.ramachandran_type) md_array_free(mol->backbone.ramachandran_type, alloc);
    if (mol->backbone.residue_idx) md_array_free(mol->backbone.residue_idx, alloc);

    // Bonds
    if (mol->bonds) md_array_free(mol->bonds, alloc);
    //md_index_data_free(&mol->connectivity, alloc);
    md_index_data_free(&mol->structures, alloc);
    md_index_data_free(&mol->rings, alloc);
    
    if (mol->hydrogen_bonds) md_array_free(mol->hydrogen_bonds, alloc);

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
    ARRAY_PUSH(atom, vx);
    ARRAY_PUSH(atom, vy);
    ARRAY_PUSH(atom, vz);
    ARRAY_PUSH(atom, radius);
    ARRAY_PUSH(atom, mass);
    ARRAY_PUSH(atom, valence);
    ARRAY_PUSH(atom, element);
    ARRAY_PUSH(atom, name);
    ARRAY_PUSH(atom, flags);
    ARRAY_PUSH(atom, residue_idx);
    ARRAY_PUSH(atom, chain_idx);

    ARRAY_PUSH(backbone, atoms);
    ARRAY_PUSH(backbone, angle);
    ARRAY_PUSH(backbone, secondary_structure);
    ARRAY_PUSH(backbone, ramachandran_type);
    ARRAY_PUSH(backbone, residue_idx);

    md_array_push_array(dst->backbone.range, src->backbone.range, src->backbone.range_count, alloc);

    ARRAY_PUSH(chain, id);
    ARRAY_PUSH(chain, residue_range);
    ARRAY_PUSH(chain, atom_range);

    md_array_push_array(dst->hydrogen_bonds, src->hydrogen_bonds, md_array_size(src->hydrogen_bonds), alloc);
    md_array_push_array(dst->bonds, src->bonds, md_array_size(src->bonds), alloc);

    ARRAY_PUSH(residue, name);
    ARRAY_PUSH(residue, id);
    ARRAY_PUSH(residue, atom_range);

    dst->atom.count           = src->atom.count;
    dst->backbone.count       = src->backbone.count;
    dst->backbone.range_count = src->backbone.range_count;
    dst->chain.count          = src->chain.count;
    dst->residue.count        = src->residue.count;
}

#undef ARRAY_PUSH
#undef ARRAY_INCREMENT
#undef ARRAY_INCREMENT_FIELD

#ifdef __cplusplus
}
#endif
