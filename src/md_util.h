#pragma once

#include <md_molecule.h>
#include <core/md_atomic.h>

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
    MD_UTIL_POSTPROCESS_CHAINS_BIT          = 0x0010,
    MD_UTIL_POSTPROCESS_BACKBONE_BIT        = 0x0020,
    MD_UTIL_POSTPROCESS_RESIDUE_BIT         = 0x0040,
    MD_UTIL_POSTPROCESS_STRUCTURE_BIT       = 0x0080,
    MD_UTIL_POSTPROCESS_ION_BIT             = 0x0100,
    MD_UTIL_POSTPROCESS_ORDER_BIT           = 0x0200,

    MD_UTIL_POSTPROCESS_ALL                 = -1,
    MD_UTIL_POSTPROCESS_COARSE_GRAINED      = MD_UTIL_POSTPROCESS_RADIUS_BIT | MD_UTIL_POSTPROCESS_MASS_BIT
};

typedef uint32_t md_util_postprocess_flags_t;

// Access to the static arrays (preserved for direct access)
const str_t* md_util_element_symbols(void);
const str_t* md_util_element_names(void);
const float* md_util_element_vdw_radii(void);

// Element functions (now calling new atomic number API internally)
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

static inline bool md_util_backbone_atoms_valid(md_protein_backbone_atoms_t prot) {
    return (prot.ca != prot.c) && (prot.ca != prot.o) && (prot.c != prot.o);
}

// Element guess function - now delegates to the new inference system
bool md_util_element_guess(md_element_t element[], size_t capacity, const struct md_molecule_t* mol);

bool md_util_element_from_mass(md_element_t out_element[], const float in_mass[], size_t count);

// Conservative mass→element mapping for LAMMPS with CG/reduced-units detection
// Returns false if data appears to be CG or reduced-units (skips mapping)
bool md_util_lammps_element_from_mass(md_element_t out_element[], const float in_mass[], size_t count);

// Computes secondary structures from backbone atoms
// Does not allocate any data, it assumes that secondary_structures has the same length as mol.backbone.count
bool md_util_backbone_secondary_structure_compute(md_secondary_structure_t secondary_structures[], size_t capacity, const struct md_molecule_t* mol);

// Computes backbone angles from backbone atoms
// Does not allocate any data, assumes that backbone_angles has the same length as args->backbone.count
bool md_util_backbone_angles_compute(md_backbone_angles_t backbone_angles[], size_t capacity, const struct md_molecule_t* mol);

// Classifies the ramachandran type (General / Glycine / Proline / Preproline) from the residue name
bool md_util_backbone_ramachandran_classify(md_ramachandran_type_t ramachandran_types[], size_t capacity, const struct md_molecule_t* mol);

void md_util_covalent_bonds_compute_exp(md_bond_data_t* out_bonds, const float* in_x, const float* in_y, const float* in_z, const md_element_t* in_elem, size_t atom_count, const md_residue_data_t* in_res, const md_unit_cell_t* in_cell, struct md_allocator_i* alloc);

// Computes the covalent bonds based from a heuristic approach, uses the covalent radius (derived from element) to determine the appropriate bond
// length. atom_res_idx is an optional parameter and if supplied, it will limit the covalent bonds to only within the same or adjacent residues.
static inline void md_util_covalent_bonds_compute(md_bond_data_t* out_bonds, const md_molecule_t* mol, struct md_allocator_i* alloc) {
    md_util_covalent_bonds_compute_exp(out_bonds, mol->atom.x, mol->atom.y, mol->atom.z, mol->atom.element, mol->atom.count, &mol->residue, &mol->unit_cell, alloc);    
}

// Grow a mask by bonds up to a certain extent (counted as number of bonds from the original mask)
// Viable mask is optional and if supplied, it will limit the growth to only within the viable mask
void md_util_mask_grow_by_bonds(struct md_bitfield_t* mask, const struct md_molecule_t* mol, size_t extent, const struct md_bitfield_t* viable_mask);

// Grow a mask by radius (in Angstrom)
// Viable mask is optional and if supplied, it will limit the growth to only within the viable mask
void md_util_mask_grow_by_radius(struct md_bitfield_t* mask, const struct md_molecule_t* mol, double radius, const struct md_bitfield_t* viable_mask);

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
bool md_util_molecule_postprocess(struct md_molecule_t* mol, struct md_allocator_i* alloc, md_util_postprocess_flags_t flags);


// ### UNIT CELL ###

// Construct unit cell from orthographic extents
md_unit_cell_t md_util_unit_cell_from_extent(double x, double y, double z);

// Construct possibly a triclinic cell from extents a,b,c and axis angles alpha, beta, gamma (in degrees)
md_unit_cell_t md_util_unit_cell_from_extent_and_angles(double a, double b, double c, double alpha, double beta, double gamma);

// Construct cell from 3x3 matrix
md_unit_cell_t md_util_unit_cell_from_matrix(float M[3][3]);

//Construct cell from upper triangular matrix components
md_unit_cell_t md_util_unit_cell_from_triclinic(double x, double y, double z, double xy, double xz, double yz);

// Computes an array of distances between two sets of coordinates in a periodic domain (cell)
// out_dist:  Output array of distances, must have length of (num_a * num_b)
// coord_a:   Array of coordinates (a)
// num_a:     Length of coord_a
// coord_b:   Array of coordinates (b)
// num_b:     Length of coord_b
// cell:      Periodic boundary cell
void md_util_unit_cell_distance_array(float* out_dist_arr, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unit_cell_t* cell);

float md_util_unit_cell_min_distance(int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unit_cell_t* cell);
float md_util_unit_cell_max_distance(int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unit_cell_t* cell);

void md_util_min_image_vec3(vec3_t* in_out_dx, size_t count, const md_unit_cell_t* unit_cell);
void md_util_min_image_vec4(vec4_t* in_out_dx, size_t count, const md_unit_cell_t* unit_cell);

// Applies periodic boundary conditions to coordinates
bool md_util_pbc(float* in_out_x, float* in_out_y, float* in_out_z, const int32_t* in_idx, size_t count, const md_unit_cell_t* unit_cell);
bool md_util_pbc_vec4(vec4_t* in_out_xyzw, size_t count, const md_unit_cell_t* unit_cell);

// Unwraps a structure
// It implicitly uses the previous coordinate as a reference when deperiodizing
bool md_util_unwrap(float* in_out_x, float* in_out_y, float* in_out_z, const int32_t* in_idx, size_t count, const md_unit_cell_t* unit_cell);
bool md_util_unwrap_vec4(vec4_t* in_out_xyzw, size_t count, const md_unit_cell_t* unit_cell);

// Batch deperiodize a set of coordinates (vec4) with respect to a given reference
bool md_util_deperiodize_vec4(vec4_t* xyzw, size_t count, vec3_t ref_xyz, const md_unit_cell_t* cell);

// Computes the minimum axis aligned bounding box for a set of points with a given radius
// Indices are optional and are used to select a subset of points, the count dictates the number of elements to process
void md_util_aabb_compute     (float out_ext_min[3], float out_ext_max[3], const float* in_x, const float* in_y, const float* in_z, const float* in_r, const int32_t* in_idx, size_t count);
void md_util_aabb_compute_vec4(float out_ext_min[3], float out_ext_max[3], const vec4_t* in_xyzr, const int32_t* in_idx, size_t count);

// Computes an object oriented bounding box based on the PCA of the provided points (with optional radius)
void md_util_oobb_compute     (float out_rotation[3][3], float out_ext_min[3], float out_ext_max[3], const float* in_x, const float* in_y, const float* in_z, const float* in_r, const int32_t* in_idx, size_t count, const md_unit_cell_t* unit_cell);
void md_util_oobb_compute_vec4(float out_rotation[3][3], float out_ext_min[3], float out_ext_max[3], const vec4_t* in_xyzr, const int32_t* in_idx, size_t count, const md_unit_cell_t* unit_cell);

// Computes the center of mass for a set of points with a given weight
// x,y,z / xyz: Arrays containing coordinates
// w:           Array of weights (optional): set as NULL to use equal weights
// indices:     Array of indices (optional): indices into the arrays (x,y,z,w)
// count:       Length of all arrays
// unit_cell:   The unit_cell of the system [Optional]
vec3_t md_util_com_compute(const float* in_x, const float* in_y, const float* in_z, const float* in_w, const int32_t* in_idx, size_t count, const md_unit_cell_t* unit_cell);
vec3_t md_util_com_compute_vec4(const vec4_t* in_xyzw, const int32_t* in_idx, size_t count, const md_unit_cell_t* unit_cell);

// Computes the similarity between two sets of points with given weights.
// One of the sets is rotated and translated to match the other set in an optimal fashion before the similarity is computed.
// The rmsd is the root mean squared deviation between the two sets of aligned vectors.
// coords:  Coordinate arrays [2] (x0, y0, z0), (x1, y1, z1)
// com:     Center of mass [2] (xyz0), (xyz1)
// w:       Array of weights (optional): set as NULL to use equal weights
// count:   Length of all arrays (x0, y0, z0, x1, y1, z1, w)
double md_util_rmsd_compute(const float* const in_x[2], const float* const in_y[2], const float* const in_z[2], const float* const in_w[2], const int32_t* const in_idx[2], size_t count, const vec3_t in_com[2]);
double md_util_rmsd_compute_vec4(const vec4_t* const in_xyzw[2], const int32_t* const in_idx[2], size_t count, const vec3_t in_com[2]);

// Computes linear shape descriptor weights (linear, planar, isotropic) from a covariance matrix
vec3_t md_util_shape_weights(const mat3_t* covariance_matrix);

// Perform linear interpolation of supplied coordinates
// out_x/y/z:   Destination arrays (x,y,z)
// in_x/y/z:    Source arrays [2] (x0, x1), (y0, y1), (z0, z1)
// count:       Count of coordinates (this implies that all coordinate arrays must be equal in length)
// unit_cell:   The unit_cell of the system [Optional]
// t: interpolation factor (0..1)
// @NOTE: This will assume that the input and output are padded such that it can operate on full simd width
bool md_util_interpolate_linear(float* out_x, float* out_y, float* out_z, const float* const in_x[2], const float* const in_y[2], const float* const in_z[2], size_t count, const md_unit_cell_t* unit_cell, float t);

// Perform cubic interpolation of supplied coordinates
// out_x/y/z:   Destination arrays (x,y,z)
// in_x/y/z:    Source arrays [2] (x0, x1, x2, x3), (y0, y1, y2, y3), (z0, z1, z2, z3)
// count:       Count of coordinates (this implies that all coordinate arrays must be equal in length)
// unit_cell:   The unit_cell of the system [Optional]
// t:           Interpolation factor (0..1)
// s:           Scaling factor (0..1), 0 is jerky, 0.5 is catmul rom, 1.0 is silky smooth
bool md_util_interpolate_cubic_spline(float* out_x, float* out_y, float* out_z, const float* const in_x[4], const float* const in_y[4], const float* const in_z[4], size_t count, const md_unit_cell_t* unit_cell, float t, float s);

// Spatially sorts the input positions according to morton order. This makes it easy to create spatially coherent clusters, just select ranges within this space.
// There are some larger jumps within the morton order as well, so when creating clusters from consecutive ranges, this should be considered as well.
// The result (source_indices) is an array of remapping indices. It is assumed that the user has reserved space for this.
void md_util_sort_spatial(uint32_t* source_indices, const float* x, const float* y, const float* z, size_t count);

// Spatially sorts the input positions according to morton order. This makes it easy to create spatially coherent clusters, just select ranges within this space.
// There are some larger jumps within the morton order as well, so when creating clusters from consecutive ranges, this should be considered as well.
// The result (source_indices) is an array of remapping indices. It is assumed that the user has reserved space for this.
void md_util_sort_spatial_vec3(uint32_t* source_indices, const vec3_t* xyz, size_t count);

// Sort array of uint32_t in place using radix sort
void md_util_sort_radix_inplace_uint32(uint32_t* data, size_t count);

// Sort array of uint32_t by producing a remapping array of source indices
// The source_indices represents the indices of the sorted array, i.e. source_indices[0] is the index of the smallest element in data
//void md_util_sort_radix_uint32(uint32_t* source_indices, const uint32_t* data, size_t count);

// Structure matching operations
// In many of the cases, there will be multiple matches which contain the indices, only with slight permutations.
// This is due to the symmetry of extremities found in molecules.

typedef enum {
    MD_UTIL_MATCH_LEVEL_STRUCTURE = 0,  // Match within complete structures
    MD_UTIL_MATCH_LEVEL_RESIDUE,        // Match within residues
    MD_UTIL_MATCH_LEVEL_CHAIN,          // Match within chains
} md_util_match_level_t;

typedef enum {
    MD_UTIL_MATCH_MODE_UNIQUE = 0,      // Store only unique matches
    MD_UTIL_MATCH_MODE_FIRST,           // Store the first match
    MD_UTIL_MATCH_MODE_ALL,		        // Store all matches
} md_util_match_mode_t;

typedef enum {
    MD_UTIL_MATCH_FLAGS_NO_H  = 1,              // Disregard hydrogen
    MD_UTIL_MATCH_FLAGS_NO_CH = 2,              // Disregard hydrogen connected to carbon
    MD_UTIL_MATCH_FLAGS_STRICT_EDGE_COUNT = 4,  // Enforce a strict edge count for each matched atom pair
    MD_UTIL_MATCH_FLAGS_STRICT_EDGE_TYPE  = 8,  // Enforce a matching edge type between matches
} md_util_match_flags_t;

// Performs complete structure matching within the given topology (mol) using a supplied reference structure.
md_index_data_t md_util_match_by_type   (const int ref_indices[], size_t ref_size, md_util_match_mode_t mode, md_util_match_level_t level, const md_molecule_t* mol, md_allocator_i* alloc);
md_index_data_t md_util_match_by_element(const int ref_indices[], size_t ref_size, md_util_match_mode_t mode, md_util_match_level_t level, const md_molecule_t* mol, md_allocator_i* alloc);

// Performs complete structure matching within the given topology (mol) using a supplied reference structure given as a smiles string
// The matcing results are stored into supplied idx_data
// The returned value is the number of matches found
size_t md_util_match_smiles(md_index_data_t* idx_data, str_t smiles, md_util_match_mode_t mode, md_util_match_level_t level, md_util_match_flags_t flags, const md_molecule_t* mol, md_allocator_i* alloc);

// Computes the maximum common subgraph between two structures
// The indices which maps from the source structure to the target structure is written to dst_idx_map
// The returned value is the number of common atoms
// It is assumed that the dst_idx_map has the same length as src_count
size_t md_util_match_maximum_common_subgraph_by_type(int* dst_idx_map, const int* trg_indices, size_t trg_count, const int* src_indices, size_t src_count, const md_molecule_t* mol, md_allocator_i* alloc);
size_t md_util_match_maximum_common_subgraph_by_element(int* dst_idx_map, const int* trg_indices, size_t trg_count, const int* src_indices, size_t src_count, const md_molecule_t* mol, md_allocator_i* alloc);

#ifdef __cplusplus
}
#endif
