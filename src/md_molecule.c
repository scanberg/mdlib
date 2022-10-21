#include "md_molecule.h"
#include "core/md_array.inl"

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
    if (mol->residue.internal_covalent_bond_range) md_array_free(mol->residue.internal_covalent_bond_range, alloc);
    if (mol->residue.complete_covalent_bond_range) md_array_free(mol->residue.complete_covalent_bond_range, alloc);

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
    if (mol->covalent_bond.bond) md_array_free(mol->covalent_bond.bond, alloc);
    if (mol->hydrogen_bond.bond) md_array_free(mol->hydrogen_bond.bond, alloc);

    // Instance
    if (mol->instance.atom_range) md_array_free(mol->instance.atom_range, alloc);
    if (mol->instance.label) md_array_free(mol->instance.label, alloc);
    if (mol->instance.transform) md_array_free(mol->instance.transform, alloc);

    memset(mol, 0, sizeof(md_molecule_t));
}

#ifdef __cplusplus
}
#endif