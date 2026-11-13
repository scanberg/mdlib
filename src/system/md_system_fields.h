#pragma once

// The arrays of md_system_t, as X macros: X(field), field a path from the system ("atom.type_idx"). Each array of the
// system is in exactly one of the lists, by what it is to inference (md_system_infer):
//
//   CORE       What the data gives: the atoms, their types, the residues, the instances and the assemblies. Inference
//              reads it and never writes it, so a system inferred on another thread can share it
//              (md_system_infer_prepare).
//   DERIVED    What inference writes: what it finds (the bonds, chains, parts, backbones, molecules, rings, lipids, fragments, the
//              chemistry, the kinds of the residues and the roles of their atoms), and what the data may give of it
//              which inference adds to or replaces (the bonds of a file, the part types of an mmCIF file).
//   DERIVED_STRINGS  Arrays of strings of DERIVED, each string allocated of its own
//
// A new array of the system goes in one of the lists, or it is neither copied (md_system_copy) nor published from a
// system inferred in the background (md_system_infer_publish). The unit test system_invariants.infer_write_set
// fails when inference writes an array of CORE.

#define MD_SYSTEM_CORE_ARRAYS(X) \
    X(atom.type_idx) \
    X(atom.type.name) X(atom.type.ff_type) X(atom.type.z) X(atom.type.mass) X(atom.type.radius) X(atom.type.color) X(atom.type.flags) \
    X(residue.name) X(residue.seq_id) X(residue.atom_offset) \
    X(instance.id) X(instance.auth_id) X(instance.res_offset) X(instance.part_type_idx) X(instance.mol_type_idx) \
    X(assembly.id) X(assembly.group_offset) X(assembly.group_inst_offset) X(assembly.inst) X(assembly.group_oper_offset) \
    X(assembly.oper) X(assembly.oper_id) X(assembly.oper_transform)

#define MD_SYSTEM_DERIVED_ARRAYS(X) \
    X(atom.flags) X(atom.formal_charge) X(atom.hydrogen_count) \
    X(residue.flags) X(residue.type_idx) \
    X(residue.type.name) X(residue.type.kind) X(residue.type.monomer) X(residue.type.color) \
    X(chain.id) X(chain.run_offset) X(chain.run) X(chain.res_chain) \
    X(part.res_offset) X(part.type_idx) X(part.type.id) X(part.type.flags) \
    X(protein_backbone.range.offset) X(protein_backbone.range.part_idx) X(protein_backbone.segment.atoms) X(protein_backbone.segment.res_idx) \
    X(nucleic_backbone.range.offset) X(nucleic_backbone.range.part_idx) X(nucleic_backbone.segment.atoms) X(nucleic_backbone.segment.res_idx) \
    X(lipid.residue_offset) X(lipid.residue_idx) X(lipid.head_offset) X(lipid.tail_offset) \
    X(lipid.head.offset) X(lipid.head.atom_idx) \
    X(lipid.tail.spine_offset) X(lipid.tail.spine_idx) X(lipid.tail.atom_offset) X(lipid.tail.atom_idx) \
    X(lipid.tail.anchor_beg) X(lipid.tail.anchor_end) X(lipid.tail.head_beg) X(lipid.tail.head_end) X(lipid.tail.kind) \
    X(bond.pairs) X(bond.flags) X(bond.conn.atom_idx) X(bond.conn.bond_idx) X(bond.conn.offset) \
    X(ring.rings) \
    X(molecule.offset) X(molecule.atom_idx) X(molecule.parent_idx) X(molecule.atom_slot) X(molecule.type_idx) X(molecule.type.flags) \
    X(fragment.mol_template) X(fragment.frag_offset) X(fragment.template_molecule) X(fragment.template_slot_offset) \
    X(fragment.slot_frag) X(fragment.template_frag_offset) X(fragment.frag_type) X(fragment.frag_slot_offset) \
    X(fragment.frag_slot) X(fragment.type.kind) X(fragment.type.color)

#define MD_SYSTEM_DERIVED_STRINGS(X) \
    X(part.type.description) X(molecule.type.name) X(fragment.type.name)

// The structs of DERIVED whole, with their counts: what md_system_infer_publish moves into the system. Of the atoms
// and the residues, which are CORE and DERIVED both, the DERIVED arrays one by one (MD_SYSTEM_DERIVED_ARRAYS) and these
// counts and structs.
#define MD_SYSTEM_DERIVED_STRUCTS(X) \
    X(residue.type) X(chain) X(part) X(protein_backbone) X(nucleic_backbone) X(lipid) X(bond) X(ring) X(molecule) X(fragment)
