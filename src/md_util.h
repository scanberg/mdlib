#pragma once

#include <md_molecule.h>

#include <core/md_str.h>
#include <core/md_vec_math.h>

struct md_allocator_i;
struct md_bitfield_t;

#ifdef __cplusplus
extern "C" {
#endif

enum {
    MD_UTIL_POSTPROCESS_NONE                = 0,
    MD_UTIL_POSTPROCESS_ELEMENT_BIT         = 0x0001,
    MD_UTIL_POSTPROCESS_RADIUS_BIT          = 0x0002,
    MD_UTIL_POSTPROCESS_MASS_BIT            = 0x0004,
    MD_UTIL_POSTPROCESS_BOND_BIT            = 0x0008,
    MD_UTIL_POSTPROCESS_CONNECTIVITY_BIT    = 0x0010,
    MD_UTIL_POSTPROCESS_CHAINS_BIT          = 0x0020,
    MD_UTIL_POSTPROCESS_BACKBONE_BIT        = 0x0040,
    MD_UTIL_POSTPROCESS_RESIDUE_BIT         = 0x0080,
    MD_UTIL_POSTPROCESS_STRUCTURE_BIT       = 0x0100,
    MD_UTIL_POSTPROCESS_ION_BIT             = 0x0200,

    MD_UTIL_POSTPROCESS_ALL                 = 0xFFFF,
    MD_UTIL_POSTPROCESS_COARSE_GRAINED      = MD_UTIL_POSTPROCESS_RADIUS_BIT | MD_UTIL_POSTPROCESS_MASS_BIT
};

typedef uint32_t md_util_postprocess_flags_t;

// This assumes the string exactly matches the value within the look up table
// The match is case sensitive and expects elements to be formatted with Big first letter and small second letter:
// E.g. H, He, Fe, Na, C
md_element_t md_util_element_lookup(str_t element_str);
md_element_t md_util_element_lookup_ignore_case(str_t element_str);

str_t md_util_element_symbol(md_element_t element);
str_t md_util_element_name(md_element_t element);

float md_util_element_vdw_radius(md_element_t element);
float md_util_element_covalent_radius(md_element_t element);
float md_util_element_atomic_mass(md_element_t element);
int   md_util_element_max_valence(md_element_t element);
uint32_t md_util_element_cpk_color(md_element_t element);

bool md_util_resname_dna(str_t str);
bool md_util_resname_acidic(str_t str);
bool md_util_resname_basic(str_t str);
bool md_util_resname_neutral(str_t str);
bool md_util_resname_water(str_t str);
bool md_util_resname_hydrophobic(str_t str);
bool md_util_resname_amino_acid(str_t str);

static inline bool md_util_backbone_atoms_valid(md_backbone_atoms_t prot) {
    return (prot.ca != prot.c) && (prot.ca != prot.o) && (prot.c != prot.o);
}

// This operation tries to deduce the element from the atom type/name which usually contains alot of cruft.
// It also tries resolve some ambiguities: Such as CA, is that Carbon Alpha or is it calcium?
// We can resolve that by looking at the residue name and in the case of Carbon Alpha, the residue name should be matched to an amino acid.
bool md_util_element_guess(md_element_t element[], int64_t capacity, const struct md_molecule_t* mol);

bool md_util_element_from_mass(md_element_t out_element[], const float in_mass[], int64_t count);

// Extracts the atom indices which are central for a segment within the backbone for a single residue
bool md_util_backbone_atoms_extract_from_residue_idx(md_backbone_atoms_t* backbone_atoms, md_residue_idx_t res_idx, const md_molecule_t* mol);

// Computes secondary structures from backbone atoms
// Does not allocate any data, it assumes that secondary_structures has the same length as args->backbone.count
bool md_util_backbone_secondary_structure_compute(md_secondary_structure_t secondary_structures[], int64_t capacity, const struct md_molecule_t* mol);

// Computes backbone angles from backbone atoms
// Does not allocate any data, assumes that backbone_angles has the same length as args->backbone.count
bool md_util_backbone_angles_compute(md_backbone_angles_t backbone_angles[], int64_t capacity, const struct md_molecule_t* mol);

// Classifies the ramachandran type (General / Glycine / Proline / Preproline) from the residue name
bool md_util_backbone_ramachandran_classify(md_ramachandran_type_t ramachandran_types[], int64_t capacity, const struct md_molecule_t* mol);

// Computes the covalent bonds based from a heuristic approach, uses the covalent radius (derived from element) to determine the appropriate bond
// length. atom_res_idx is an optional parameter and if supplied, it will limit the covalent bonds to only within the same or adjacent residues.
md_bond_data_t md_util_compute_covalent_bonds(const md_atom_data_t* atom_data, const md_residue_data_t* res_data, const md_unit_cell_t* cell, struct md_allocator_i* alloc);

// Grow a mask by bonds up to a certain extent (counted as number of bonds from the original mask)
// Viable mask is optional and if supplied, it will limit the growth to only within the viable mask
void md_util_grow_mask_by_bonds(struct md_bitfield_t* mask, const struct md_molecule_t* mol, int extent, const struct md_bitfield_t* viable_mask);

// Grow a mask by radius (in Angstrom)
// Viable mask is optional and if supplied, it will limit the growth to only within the viable mask
void md_util_grow_mask_by_radius(struct md_bitfield_t* mask, const struct md_molecule_t* mol, float radius, const struct md_bitfield_t* viable_mask);

//bool md_util_compute_hydrogen_bonds(md_bond_data_t* dst, const float* atom_x, const float* atom_y, const float* atom_z, const md_element_t* atom_elem, int64_t atom_count, vec3_t pbc_ext, struct md_allocator_i* alloc);

// Computes chains from connected residues
// The definition of a chain here is a linear sequence of residues which are connected by covalent bonds.
// This means that single residues which are not connected to any other residue will not classify as a chain.
//bool md_util_compute_chain_data(md_chain_data_t* chain_data, const md_residue_idx_t atom_residue_idx[], int64_t atom_count, const md_bond_t bonds[], int64_t covelent_bond_count, struct md_allocator_i* alloc);

// Compute the valence of atoms given the covalent bonds
//bool md_util_compute_atom_valence(md_valence_t atom_valence[], int64_t atom_count, const md_bond_t bonds[], int64_t bond_count);

// Compute rings formed by covalent bonds
//bool md_util_compute_rings(md_index_data_t* ring_data, int64_t atom_count, const md_bond_t bonds[], int64_t bond_count, struct md_allocator_i* alloc);

// Identify isolated structures by covalent bonds
//bool md_util_compute_structures(md_index_data_t* structures, int64_t atom_count, const md_bond_t bonds[], int64_t bond_count, struct md_allocator_i* alloc);

// Attempts to generate missing data such as covalent bonds, chains, secondary structures, backbone angles etc.
bool md_util_postprocess_molecule(struct md_molecule_t* mol, struct md_allocator_i* alloc, md_util_postprocess_flags_t flags);


// ### CELL ###

// Construct cells from orthographic extents
md_unit_cell_t md_util_unit_cell_from_extent(double x, double y, double z);

// Construct possibly a triclinic cell from extents a,b,c and axis angles alpha, beta, gamma (in degrees)
md_unit_cell_t md_util_unit_cell_from_extent_and_angles(double a, double b, double c, double alpha, double beta, double gamma);

// Construct cell from mat3 basis
md_unit_cell_t md_util_unit_cell_from_matrix(mat3_t M);

//Construct cell from triclinic basis
md_unit_cell_t md_util_unit_cell_from_triclinic(double x, double y, double z, double xy, double xz, double yz);

// Computes an array of distances between two sets of coordinates in a periodic domain (cell)
// out_dist:  Output array of distances, must have length of (num_a * num_b)
// coord_a:   Array of coordinates (a)
// num_a:     Length of coord_a
// coord_b:   Array of coordinates (b)
// num_b:     Length of coord_b
// cell:      Periodic boundary cell
void md_util_unit_cell_distance_array(float* out_dist_arr, const vec3_t* coord_a, int64_t num_a, const vec3_t* coord_b, int64_t num_b, const md_unit_cell_t* cell);

float md_util_unit_cell_min_distance(int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, int64_t num_a, const vec3_t* coord_b, int64_t num_b, const md_unit_cell_t* cell);
float md_util_unit_cell_max_distance(int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, int64_t num_a, const vec3_t* coord_b, int64_t num_b, const md_unit_cell_t* cell);

// Applies periodic boundary conditions to coordinates of atoms within molecule
// It ensures that residues and chains reside within the same period
//bool md_util_pbc_ortho(struct md_molecule_t* mol, vec3_t pbc_ext);


bool md_util_pbc_ortho(float* x, float* y, float* z, int64_t count, vec3_t box);

bool md_util_unwrap_ortho(float* x, float* y, float* z, md_index_data_t structures, vec3_t box);

// Deperiodizes the coordinates of an entire system and unwraps structures defined given by the covalent bonds across the periodic boundaries.
// If finally ensures that the center of mass of all structures (including individual atoms) reside within box.
bool md_util_deperiodize_system(float* x, float* y, float* z, const md_unit_cell_t* cell, const struct md_molecule_t* mol);

// Computes the minimum axis aligned bounding box for a set of points with a given radius
// Indices are optional and are used to select a subset of points, the count dictates the number of elements to process
void md_util_compute_aabb(vec3_t* aabb_min, vec3_t* aabb_max, const float* x, const float* y, const float* z, const float* r, const int32_t* indices, int64_t count);
void md_util_compute_aabb_vec4(vec3_t* aabb_min, vec3_t* aabb_max, const vec4_t* xyzr, const int32_t* indices, int64_t count);

// Computes the center of mass for a set of points with a given weight
// x,y,z / xyz: Arrays containing coordinates
// w:           Array of weights (optional): set as NULL to use equal weights
// indices:     Array of indices (optional): indices into the arrays (x,y,z,w)
// count:       Length of all arrays
vec3_t md_util_compute_com(const float *x, const float* y, const float* z, const float* w, const int32_t* indices, int64_t count);
vec3_t md_util_compute_com_vec4(const vec4_t* xyzw, const int32_t* indices, int64_t count);

// Computes the center of mass for a set of points with a given weight given in orthogonal periodic boundary conditions
// The indices used to access the arrays are given in the indices array
// x,y,z / xyz: Array of coordinates
// w:           Array of weights (optional): set as NULL to use equal weights
// indices:     Array of indices (optional): indices into the arrays (x,y,z,w)
// count:       Number of elements to process, either the length of the indices array or the length of the x,y,z,w arrays
// box:         Extent of periodic boundary box (optional per component): Set to zero if pbc does not apply in that dimension
vec3_t md_util_compute_com_ortho(const float *x, const float* y, const float* z, const float* w, const int32_t* indices, int64_t count, vec3_t box);
vec3_t md_util_compute_com_vec4_ortho(const vec4_t* xyzw, const int32_t* indices, int64_t count, vec3_t box);

// Computes the similarity between two sets of points with given weights.
// One of the sets is rotated and translated to match the other set in an optimal fashion before the similarity is computed.
// The rmsd is the root mean squared deviation between the two sets of aligned vectors.
// coords:  Coordinate arrays [2] (x0, y0, z0), (x1, y1, z1)
// com:     Center of mass [2] (xyz0), (xyz1)
// w:       Array of weights (optional): set as NULL to use equal weights
// count:   Length of all arrays (x0, y0, z0, x1, y1, z1, w)
double md_util_compute_rmsd(const md_vec3_soa_t coord[2], const vec3_t com[2], const float* w, int64_t count);

double md_util_compute_rmsd_vec4(const vec4_t* xyzw[2], const vec3_t com[2], int64_t count);

// Computes linear shape descriptor weights (linear, planar, isotropic) from a covariance matrix
vec3_t md_util_shape_weights(const mat3_t* covariance_matrix);

// Perform linear interpolation of supplied coordinates
// dst_coord:   Destination arrays (x,y,z)
// src_coord:   Source arrays [2] (x0, y0, z0), (x1, y1, z1)
// count:       Count of coordinates (this implies that all coordinate arrays must be equal in length)
// box:         Extent of periodic boundary box (optional) set to zero if should be ignored
// t: interpolation factor (0..1)
bool md_util_linear_interpolation(md_vec3_soa_t dst_coord, const md_vec3_soa_t src_coord[2], int64_t count, vec3_t box, float t);

// Perform cubic interpolation of supplied coordinates
// dst_coord:   Destination arrays (x,y,z)
// src_coord:   Source arrays [4] (x0, y0, z0), (x1, y1, z1), (x2, y2, z2), (x3, y3, z3)
// count:       Count of coordinates (this implies that all coordinate arrays must be equal in length)
// box:         Extent of periodic boundary (optional): Set to zero if pbc does not apply in that dimension
// t:           Interpolation factor (0..1)
// s:           Scaling factor (0..1), 0 is jerky, 0.5 is catmul rom, 1.0 is silky smooth
bool md_util_cubic_spline_interpolation(md_vec3_soa_t dst_coord, const md_vec3_soa_t src_coord[4], int64_t count, vec3_t box, float t, float s);

// Spatially sorts the input positions according to morton order. This makes it easy to create spatially coherent clusters, just select ranges within this space.
// There are some larger jumps within the morton order as well, so when creating clusters from consecutive ranges, this should be considered as well.
// The result (source_indices) is an array of remapping indices. It is assumed that the user has reserved space for this.
void md_util_spatial_sort_soa(uint32_t* source_indices, const float* x, const float* y, const float* z, int64_t count);

// Spatially sorts the input positions according to morton order. This makes it easy to create spatially coherent clusters, just select ranges within this space.
// There are some larger jumps within the morton order as well, so when creating clusters from consecutive ranges, this should be considered as well.
// The result (source_indices) is an array of remapping indices. It is assumed that the user has reserved space for this.
void md_util_spatial_sort(uint32_t* source_indices, const vec3_t* xyz, int64_t count);

// Structure matching operations
// In many of the cases, there will be multiple matches which contain the indices, only with slight permutations.
// This is due to the symmetry of extremities found in molecules.
// To filter the result into 1 permutation per matched set of indices, the user can pass a parameter to filter permutations by RMSD
// Which only selects the best matching permutation based on RMSD.

typedef enum {
    MD_UTIL_MATCH_LEVEL_STRUCTURE = 0,  // Match within complete structures
    MD_UTIL_MATCH_LEVEL_RESIDUE,        // Match within residues
    MD_UTIL_MATCH_LEVEL_CHAIN,          // Match within chains
} md_util_match_level_t;

typedef enum {
    MD_UTIL_MATCH_MODE_UNIQUE = 0,      // Store only unique matches
    MD_UTIL_MATCH_MODE_FIRST = 1,       // Store the first match
    MD_UTIL_MATCH_MODE_ALL = 2,		    // Store all matches
} md_util_match_mode_t;

// Performs complete structure matching within the given topology (mol) using a supplied reference structure.
md_index_data_t md_util_match_by_type(const int ref_indices[], int64_t ref_size, md_util_match_mode_t mode, md_util_match_level_t level, const md_molecule_t* mol, md_allocator_i* alloc);
md_index_data_t md_util_match_by_element(const int ref_indices[], int64_t ref_size, md_util_match_mode_t mode, md_util_match_level_t level, const md_molecule_t* mol, md_allocator_i* alloc);

// Performs complete structure matching within the given topology (mol) using a supplied reference structure given as a smiles string
md_index_data_t md_util_match_smiles(str_t smiles, md_util_match_mode_t mode, md_util_match_level_t level, const md_molecule_t* mol, md_allocator_i* alloc);

// Computes the maximum common subgraph between two structures
// The indices which maps from the source structure to the target structure is written to dst_idx_map
// The returned value is the number of common atoms
// It is assumed that the dst_idx_map has the same length as src_count
int64_t md_util_match_maximum_common_subgraph_by_type(int* dst_idx_map, const int* trg_indices, int64_t trg_count, const int* src_indices, int64_t src_count, const md_molecule_t* mol, md_allocator_i* alloc);
int64_t md_util_match_maximum_common_subgraph_by_element(int* dst_idx_map, const int* trg_indices, int64_t trg_count, const int* src_indices, int64_t src_count, const md_molecule_t* mol, md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
