#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <md_types.h>

// The levels a system is built of above its atoms: residues, instances, chains, parts and molecules, and their types.

// Residue type: the residues of the system of one name and kind. The ALA of every chain is one type, the HIE and
// the HID of a protein are two, which are both forms of the standard monomer HIS. Inferred from the residues
// (md_util_system_infer_residue_types), as the residues are classified late (a lipid is found by its structure
// after the residues are classified). Like the particle types they carry what the residues look like, their color.
typedef struct md_residue_type_data_t {
    size_t count;
    md_label_t*        name;
    md_residue_kind_t* kind;
    md_monomer_t*      monomer;     // The standard monomer it is a form of, MD_MONOMER_NONE if none (md_monomer_from_residue_name)
    uint32_t*          color;
} md_residue_type_data_t;

// Residue: a contiguous range of atoms
typedef struct md_residue_data_t {
    size_t count;
    md_label_t* name;
    md_sequence_id_t* seq_id;
    uint32_t* atom_offset;          // [count + 1]
    md_residue_flags_t* flags;      // Kind (amino acid, nucleotide, water, ion) and place in the chain, see md_types.h

    md_residue_type_idx_t* type_idx;    // [count] Once inferred (md_util_system_infer_residue_types), NULL before
    md_residue_type_data_t type;
} md_residue_data_t;

// THE HIERARCHY. A system is built of atoms, residues, parts and molecules, each nested in the next and each of a type:
//
//     atom       particle type      (md_atom_data_t)
//     residue    residue type       (md_residue_data_t)
//     part       part type          (md_part_data_t)       a strand, a glycan tree, a ligand, an ion, a water
//     molecule   molecule type      (md_molecule_data_t)   what the bonds hold together
//
// A chain is no level of it: it is a label, the author's chain id, which the data files atoms of any level under
// (md_chain_data_t). The instances below are no level either: they are what the data gives, the input the parts, the
// molecules and the chains are inferred from.

// Instance: a unit of the system AS THE DATA GIVES IT, a contiguous range of residues: an asym of an mmCIF file (a
// polymer chain, or one of the molecules of a non-polymer asym), a run of residues of one chain id of a PDB file
// (md_util_system_instances_from_chain_ids), a molecule of a topology (a GROMACS molecule, a molecule of a LAMMPS data
// file). Instances are never inferred: where the data gives no units (a GRO file, an XYZ file) the system has none. The
// instances are back to back, so a run of residues the data names nothing (the residues without a chain id between two
// chains of a PDB file) is an instance of no id.
//
// THE INSTANCES ARE INPUT. What they are differs with the data (an asym is a part of a molecule, a molecule of a topology
// a whole one), so nothing but what is inferred from them reads them: the parts (which never reach over an instance), the
// molecules (an instance of a molecule type is one molecule whatever its bonds, and md_match closes its molecules the
// same way), the chains (from their author ids), and the assemblies, which name the asyms of a file by index, as these stay put when the parts are inferred anew. Anything
// else uses the parts and the molecules.
//
// id is the instance's label (mmCIF label_asym_id) and auth_id the author's chain id (auth_asym_id, the chain id of a
// PDB file), "" when there is none: labels are the data's, never made up. The id is unique for polymer chains but
// SHARED by the small molecules of one asym: mmCIF puts all the waters of a chain in one asym. Look up instances by
// index, not by id. The id is not what a chain is called: that is the auth_id, which the chains are made from
// (md_chain_data_t).
typedef struct md_instance_data_t {
    size_t count;
    md_label_t* id;
    md_label_t* auth_id;
    uint32_t* res_offset;                   // [count + 1]
    md_part_type_idx_t* part_type_idx;      // [count] The part type the data gives it (an mmCIF entity), -1 for none
    md_molecule_type_idx_t* mol_type_idx;   // [count] The molecule type the data gives it (a molecule type of a topology),
                                            // -1 for none. Such an instance is one molecule, whatever its bonds.
} md_instance_data_t;

// Chain: the atoms which the data files under one chain identifier, the author's: the chain id of a PDB file, the
// auth_asym_id of an mmCIF file. Whatever carries the id is of the chain, a polymer and the ligands and waters filed
// with it alike, and an id is one chain however its atoms are spread: chain("A") of a script is one set.
//   - A chain is a label, never inferred: where the data gives no ids the system has no chains, and a strand of a
//     polymer is found as a part (polymer() of a script).
// A chain is runs of residues: a strand, and the ligands and waters filed after all of the polymers under
// its id. The chains are in the order of their first residue. They are derived from the instances
// (md_util_system_infer_chains), whose author ids hold the identifiers.
typedef struct md_chain_data_t {
    size_t count;
    md_label_t* id;                 // [count] The author's identifier
    uint32_t* run_offset;           // [count + 1] The runs of chain c: run[run_offset[c], run_offset[c + 1])
    md_urange_t* run;               // Runs of consecutive residues, ascending within each chain
    md_chain_idx_t* res_chain;      // [residue count] The chain of each residue, -1 for one in none
} md_chain_data_t;

// Part type: the parts of the system of one kind and one sequence of residues. Given by the data (an mmCIF entity:
// a polymer sequence, a chemical component, a branched glycan), else inferred with the parts
// (MD_PART_TYPE_FLAG_INFERRED), which gives a kind and a sequence of residue names one type. What a part is lives in the
// flags of its type (md_part_kind_t). The given types come first, the inferred after them.
typedef struct md_part_type_data_t {
    size_t count;
    md_label_t* id;                     // [count] The data's id (an mmCIF entity id), else its position
    md_part_type_flags_t* flags;        // [count] Kind (peptide, DNA, water, ...) and whether it was inferred, see md_types.h
    md_str_t* description;                 // [count]
} md_part_type_data_t;

// Part: one polymer strand, one glycan tree, one ligand, one lipid, one ion or one water: the parts the molecules are
// built of, each a run of consecutive residues, within one instance where the system has instances. A molecule may hold several (a
// glycoprotein, a strand and its glycans; the two strands of a duplex a topology gives as one molecule type), a
// chain likewise (a strand and the ligands and waters filed under its id). Every part is within one molecule.
//   - A polymer strand goes on over a gap (residues missing from a crystal structure) and ends where the chemistry says
//     the chain does (CHAIN BREAKS in md_util.h). A residue of no polymer bonded to it in the run (a cap, a modified
//     residue) is of it. A bond from a side chain (a disulfide, a glycosylation) joins no parts.
//   - Residues of no polymer are one part while each is bonded to one of it: a glycan tree, a lipid of several residues.
//     A coordination bond between residues joins none: a metal is a part of its own, as each ion and water is.
// They are derived from the instances, the bonds and the kinds of the residues (md_util_system_infer_parts). A part
// which is all of its instance is of the instance's part type; any other is of an inferred one.
typedef struct md_part_data_t {
    size_t count;
    uint32_t* res_offset;               // [count + 1]
    md_part_type_idx_t* type_idx;       // [count]
    md_part_type_data_t type;
} md_part_data_t;

// Molecule: one molecule of the system, molecule() of the script. A connected component of the graph of the
// persistent bonds, closed over what holds a molecule together where the bonds do not say so: the atoms of a residue
// (an atom no bond reaches, a virtual site) and the residues of a polymer chain (over a gap, residues missing from a
// crystal structure). So a residue and a polymer chain are always within one molecule; a molecule may hold several
// chains (joined by a disulfide), and a ligand bonded to a protein is of its molecule. A metal coordinated by other
// residues is not, nor are the waters of its shell: the coordination bonds between residues hold no molecule together,
// unless a topology gives them (a bond within one of its molecules, as GROMACS has it). A molecule of a topology is
// one molecule whatever its bonds (md_instance_data_t::mol_type_idx): the two strands of a DNA duplex which martinize-dna
// writes as one molecule type, held together by an elastic network. A molecule is whole parts. Coarse grained systems are held together by their
// hierarchy instead (md_util_system_infer_molecules).
// It also holds the iteration order of atoms and parent-child relationships for traversing the molecule as a tree,
// which is what unwrapping uses (md_system_molecule).
typedef struct md_molecule_t {
    const int32_t* atom_idx;
    const int32_t* parent_idx;
    size_t   count;
} md_molecule_t;

// The molecules of a system, stored contiguously: offset indexes the atom_idx and parent_idx arrays for each molecule,
// from which a single one (md_molecule_t) is read with md_system_molecule.
//
// The index into atom_idx / parent_idx is referred to as a SLOT. Both arrays are indexed by slot
// and both hold GLOBAL atom indices:
//   atom_idx[slot]   - the atom
//   parent_idx[slot] - the atom it was reached from during the traversal
// The root of each molecule is its own parent, so a root is parent_idx[slot] == atom_idx[slot] and
// consumers need no sentinel branch. There is exactly one root per molecule.
//
// atom_slot is the reverse map, indexed by GLOBAL atom index and covering every atom in the system:
//   atom_slot[atom] -> the slot which holds that atom
// It turns "what is the parent of this atom" into an O(1) lookup, which is what lets a caller walk
// the hierarchy upwards from an arbitrary subset of atoms rather than only downwards from a root
// (see md_util_unwrap_molecule).
//
// ROOT SELECTION:
//   The root is the most topologically central atom of the molecule - the atom whose greatest edge
//   distance to any other atom in the molecule is smallest, i.e. the graph center. This bounds the
//   depth of the traversal, which bounds both the cost of walking the hierarchy upwards and the
//   number of minimum image steps accumulated when unwrapping. It also makes the root a property of
//   the molecule rather than of the atom ordering in the source file, so two identical molecules
//   stored in different orders get the same traversal.
//
//   The center is located with the standard double sweep (farthest atom from an arbitrary seed,
//   then farthest atom from that, then the midpoint of the path between them). That is exact for
//   acyclic structures and a close approximation in the presence of rings, which are local and
//   small in practice.
//
// INVARIANT - load bearing, md_util_unwrap_molecule depends on it:
//   Within a molecule, atoms are stored in BFS order from the root, so every atom appears after
//   the one it was reached from: slot(parent) < slot(child). A single forward pass therefore always
//   sees a placed parent, and sorting an arbitrary set of slots ascending yields a valid
//   topological order.
//   Do not reorder atom_idx (for locality or otherwise) without reordering parent_idx and updating
//   atom_slot to match.
//
// MOLECULE TYPES. The molecule of a topology is of the molecule type the topology gives it (a GROMACS molecule type, by
// its name). Any other molecule is of an inferred type (MD_MOLECULE_TYPE_FLAG_INFERRED), shared by the molecules made
// alike: of the same part types in the same order, as many atoms of each particle type and as many bonds between each
// pair of particle types (counted, not placed: two molecules of these which are bonded differently share a type). The given
// types come first, the inferred after them, and are inferred anew with the molecules.
typedef struct md_molecule_type_data_t {
    size_t count;
    md_str_t* name;                        // [count] The topology's name of it, else what its parts are
    md_molecule_type_flags_t* flags;    // [count]
} md_molecule_type_data_t;

typedef struct md_molecule_data_t {
    size_t count;
    uint32_t* offset;    // [count + 1] The slots of molecule m: [offset[m], offset[m + 1])
    int32_t* atom_idx;   // [slot] -> global atom index
    int32_t* parent_idx; // [slot] -> global atom index of the atom it was reached from (a root is its own parent)
    int32_t* atom_slot;  // [global atom index] -> slot. Length is the system atom count

    md_molecule_type_idx_t* type_idx;   // [count]
    md_molecule_type_data_t type;
} md_molecule_data_t;

#ifdef __cplusplus
extern "C" {
#endif

// The index i of the range [offset[i], offset[i + 1]) which holds value, -1 when there is none.
// offset is non decreasing and count + 1 long. A binary search.
static inline int32_t md_offset_range_find(const uint32_t* offset, size_t count, size_t value) {
    if (!offset || count == 0 || value < offset[0] || value >= offset[count]) return -1;
    size_t lo = 0, hi = count;  // offset[lo] <= value < offset[hi]
    while (hi - lo > 1) {
        const size_t mid = lo + (hi - lo) / 2;
        if (offset[mid] <= value) lo = mid;
        else hi = mid;
    }
    return (int32_t)lo;
}

// RESIDUE TYPES (md_residue_type_data_t): a table, as the particle types are, read through the table
// (&sys->residue.type) like the particle types (&sys->atom.type)

static inline size_t md_residue_type_count(const md_residue_type_data_t* type) {
    MD_ASSERT(type);
    return type->count;
}

// The type of the given name and kind, -1 when there is none
static inline md_residue_type_idx_t md_residue_type_find(const md_residue_type_data_t* type, md_str_t name, md_residue_kind_t kind) {
    MD_ASSERT(type);
    for (size_t i = 0; i < type->count; ++i) {
        if (type->kind[i] == kind && md_str_eq(LBL_TO_STR(type->name[i]), name)) {
            return (md_residue_type_idx_t)i;
        }
    }
    return -1;
}

static inline md_residue_type_idx_t md_residue_type_add(md_residue_type_data_t* type, md_str_t name, md_residue_kind_t kind, md_monomer_t monomer, uint32_t color, struct md_allocator_i* alloc) {
    MD_ASSERT(type);
    MD_ASSERT(alloc);
    md_array_push(type->name, make_label(name), alloc);
    md_array_push(type->kind, kind, alloc);
    md_array_push(type->monomer, monomer, alloc);
    md_array_push(type->color, color, alloc);
    type->count++;
    return (md_residue_type_idx_t)(type->count - 1);
}

static inline md_str_t md_residue_type_name(const md_residue_type_data_t* type, size_t type_idx) {
    MD_ASSERT(type);
    if (type_idx < type->count) return LBL_TO_STR(type->name[type_idx]);
    return MD_STR_LIT("");
}

static inline md_residue_kind_t md_residue_type_kind(const md_residue_type_data_t* type, size_t type_idx) {
    MD_ASSERT(type);
    return type_idx < type->count ? type->kind[type_idx] : MD_RESIDUE_KIND_OTHER;
}

static inline md_monomer_t md_residue_type_monomer(const md_residue_type_data_t* type, size_t type_idx) {
    MD_ASSERT(type);
    return type_idx < type->count ? type->monomer[type_idx] : MD_MONOMER_NONE;
}

static inline uint32_t md_residue_type_color(const md_residue_type_data_t* type, size_t type_idx) {
    MD_ASSERT(type);
    return type_idx < type->count ? type->color[type_idx] : 0;
}

// The type tables of the parts and of the molecules (md_part_type_data_t, md_molecule_type_data_t), allocated from alloc,
// the system's. A push copies the strings and returns the index of the new type; an empty id of a part type is made its
// position ("1", "2", ...). A shrink drops the types from count on, which is how inference drops the inferred types it
// is about to infer anew (they come after the given ones).
void md_part_type_data_free(md_part_type_data_t* type, md_allocator_i* alloc);
void md_part_type_data_shrink(md_part_type_data_t* type, size_t count, md_allocator_i* alloc);
md_part_type_idx_t md_part_type_data_push(md_part_type_data_t* type, md_str_t id, md_part_type_flags_t flags, md_str_t description, md_allocator_i* alloc);

void md_molecule_type_data_free(md_molecule_type_data_t* type, md_allocator_i* alloc);
void md_molecule_type_data_shrink(md_molecule_type_data_t* type, size_t count, md_allocator_i* alloc);
md_molecule_type_idx_t md_molecule_type_data_push(md_molecule_type_data_t* type, md_str_t name, md_molecule_type_flags_t flags, md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
