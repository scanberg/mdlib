#ifndef _MD_UTIL_H_
#define _MD_UTIL_H_

#include "md_molecule.h"
#include <core/md_str.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator;
struct md_trajectory_box;

typedef struct md_util_backbone_args {
    struct {
        int64_t count;
        const char** name;
    } atom;

    struct {
        int64_t count;
        const md_range* atom_range;
    } residue;
} md_util_backbone_args_t;

typedef struct md_util_secondary_structure_args {
    struct {
        int64_t count;
        const float* x;
        const float* y;
        const float* z;
    } atom;

    struct {
        int64_t count;
        const md_backbone_atoms* atoms;
    } backbone;

    struct {
        int64_t count;
        const md_range* backbone_range;
    } chain;
} md_util_secondary_structure_args_t;

typedef struct md_util_backbone_angle_args {
    struct {
        int64_t count;
        const float* x;
        const float* y;
        const float* z;
    } atom;

    struct {
        int64_t count;
        const md_backbone_atoms* atoms;
    } backbone;

    struct {
        int64_t count;
        const md_range* backbone_range;
    } chain;
} md_util_backbone_angle_args_t;

typedef struct md_util_covalent_bond_args {
    struct {
        int64_t count;
        const float* x;
        const float* y;
        const float* z;
        const md_element* element;
    } atom;

    struct {
        int64_t count;
        const md_range* atom_range;
        // These are optional parameters and if it is not null, it will be filled with the corresponding ranges of bonds created
        md_range* internal_bond_range;
        md_range* complete_bond_range;
    } residue;
} md_util_covalent_bond_args_t;

// @TODO: Complete this with the relevant fields required for computing hydrogen bonds
typedef struct md_util_hydrogen_bond_args {
    struct {
        int64_t count;
        const float* x;
        const float* y;
        const float* z;
        const md_element* element;
    } atom;
} md_util_hydrogen_bond_args_t;

// This assumes the string exactly matches the value within the look up table
// The match is case sensitive and expects elements to be formatted with Big first letter and small second letter:
// E.g. H, He, Fe, Na, C
md_element md_util_lookup_element(str_t element_str);

// This is a more relaxed version which tries to extract the element from the atom type/name which usually contains alot more cruft.
// If Optionally can take a string to the residue name to resolve some ambiguities: Such as CA, is that Carbon Alpha or is it calcium?
// We can resolve that by looking at the residue name and in the case of Carbon Alpha, the residue name should be matched to an amino acid.
md_element md_util_decode_element(str_t atom_name, str_t res_name);

str_t md_util_element_symbol(md_element element);
str_t md_util_element_name(md_element element);

float md_util_element_vdw_radius(md_element element);
float md_util_element_covalent_radius(md_element element);
float md_util_element_atomic_mass(md_element element);
uint32_t md_util_element_cpk_color(md_element element);

bool md_util_is_resname_dna(str_t str);
bool md_util_is_resname_amino_acid(str_t str);


/*
TODO: implement these utility functions

// RESIDUE FUNCTIONS BASED ON RESIDUE NAMES
bool is_amino_acid(const char* str_ptr, int64_t str_len);
bool is_acidic(const char* str_ptr, int64_t str_len);
bool is_basic(const char* str_ptr, int64_t str_len);
bool is_neutral(const char* str_ptr, int64_t str_len);
bool is_water(const char* str_ptr, int64_t str_len);
bool is_hydrophobic(const char* str_ptr, int64_t str_len);
*/

static inline bool md_util_protein_backbone_atoms_valid(md_backbone_atoms prot) {
    return (prot.ca != prot.c) && (prot.ca != prot.o) && (prot.c != prot.o);
}

// Extracts the indices which make up the backbone from the supplied atom names within corresponding residues
bool md_util_extract_backbone_atoms(md_backbone_atoms backbone_atoms[], int64_t capacity, const md_util_backbone_args_t* args);

// Computes secondary structures from backbone atoms
// Does not allocate any data, it assumes that secondary_structures has the same length as args->backbone.count
bool md_util_compute_secondary_structure(md_secondary_structure secondary_structures[], int64_t capacity, const md_util_secondary_structure_args_t* args);

// Computes backbone angles from backbone atoms
// Does not allocate any data, assumes that backbone_angles has the same length as args->backbone.count
bool md_util_compute_backbone_angles(md_backbone_angles backbone_angles[], int64_t capacity, const md_util_backbone_angle_args_t* args);

// the result is an md_array allocated using the supplied allocator, since we cannot determine how many bonds will be formed.
md_bond* md_util_extract_covalent_bonds(const md_util_covalent_bond_args_t* args, struct md_allocator* alloc);
md_bond* md_util_extract_hydrogen_bonds(const md_util_hydrogen_bond_args_t* args, struct md_allocator* alloc);


typedef struct md_util_apply_pbc_args {
    struct {
        int64_t count;
        const float* x;
        const float* y;
        const float* z;
    } atom;

    struct {
        int64_t count;
        const md_range* atom_range;
    } residue;

    struct {
        int64_t count;
        const md_range* residue_range;
    } chain;

    struct {
        const struct md_trajectory_box* box;
    } pbc;
} md_util_apply_pbc_args_t;

bool md_util_apply_pbc(float* x, float* y, float* z, int64_t count, md_util_apply_pbc_args_t args);

#ifdef __cplusplus
}
#endif

#endif