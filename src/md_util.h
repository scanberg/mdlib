#pragma once

#include <md_system.h>

#include <core/md_str.h>
#include <core/md_vec_math.h>

struct md_allocator_i;
struct md_bitfield_t;

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    MD_UTIL_INFER_NONE                = 0,
    MD_UTIL_INFER_COLOR_BIT           = 1u << 1,
    MD_UTIL_INFER_BOND_BIT            = 1u << 2,
    MD_UTIL_INFER_INSTANCE_BIT        = 1u << 3,
    MD_UTIL_INFER_BACKBONE_BIT        = 1u << 4,
    MD_UTIL_INFER_STRUCTURE_BIT       = 1u << 5,
	//MD_UTIL_INFER_SECONDARY_STRUCTURE_BIT = 0x0400,
	MD_UTIL_INFER_UNWRAP_STRUCTURE_BIT = 1u << 7,
    MD_UTIL_INFER_CHEMISTRY_BIT       = 1u << 8,     // Bond orders, aromaticity, formal charges, hydrogen counts and hybridization, see md_chem.h

    MD_UTIL_INFER_ALL                 = -1,
} md_infer_flags_t;

ENUM_FLAGS(md_infer_flags_t)

// Access to the static arrays (preserved for direct access)
const str_t* md_util_element_symbols(void);
const str_t* md_util_element_names(void);
const float* md_util_element_vdw_radii(void);

// Element functions (now calling new atomic number API internally)
md_element_t md_util_element_lookup(str_t element_str, bool ignore_case);
str_t md_util_element_symbol(md_element_t element);
str_t md_util_element_name(md_element_t element);
float md_util_element_vdw_radius(md_element_t element);
float md_util_element_covalent_radius(md_element_t element);
float md_util_element_atomic_mass(md_element_t element);
int   md_util_element_max_valence(md_element_t element);
uint32_t md_util_element_cpk_color(md_element_t element);

bool md_util_resname_dna(str_t str);
bool md_util_resname_rna(str_t str);
bool md_util_resname_acidic(str_t str);
bool md_util_resname_basic(str_t str);
bool md_util_resname_neutral(str_t str);
bool md_util_resname_water(str_t str);
bool md_util_resname_hydrophobic(str_t str);
bool md_util_resname_amino_acid(str_t str);
bool md_util_resname_nucleotide(str_t str);

void md_util_system_extract_xyzw_from_mask(vec4_t* out_xyzw, const struct md_bitfield_t* mask, const md_system_t* sys, const md_system_state_t* state);

// Classifies the components (md_component_kind_t: amino acid, nucleotide, water, ion) from their names, atoms and
// bonds. An amino acid or nucleotide whose backbone atoms are found by name and verified by their bonds, and which
// is either a standard residue by name or linked into a chain (a peptide or phosphodiester bond to another
// component), is MD_COMPONENT_FLAG_RESOLVED: its atoms are given their roles (MD_ATOM_FLAG_BACKBONE, _SIDE_CHAIN,
// ...) and its terminal groups are marked. The link keeps out the ligands which name their atoms like a residue
// (GTP, ATP, NAD) while modified residues in a chain (MSE, PSU) come in. A component which is only named like an
// amino acid or nucleotide gets the kind alone. Components a predefined coarse grained type already classified as
// water or ion are left alone. Requires the bonds.
bool md_util_system_infer_comp_flags(md_system_t* sys);

// Replaces the entities and instances of the system with inferred ones (MD_ENTITY_FLAG_INFERRED), from the
// classified components and the bonds (see md_instance_data_t for what an instance is). opt_comp_auth_asym_id is
// the author chain id of each component, when the format has one (a PDB file): it is kept as the instances'
// auth_id, a change of it always ends an instance, and a polymer chain continues over gaps while it stays the same.
bool md_util_system_infer_entity_and_instance(md_system_t* sys, const str_t opt_comp_auth_asym_id[]);

// The id of the instance with the given index when a loader generates them: A..Z, AA..ZZ, AAA.. (0 is A)
md_label_t md_util_instance_id_from_index(size_t idx);

// Classifies the entities of kind MD_ENTITY_KIND_UNKNOWN (which a topology names without saying what they are)
// from the components of their first instance. Requires the components to be classified.
void md_util_system_infer_entity_kinds(md_system_t* sys);

size_t md_util_element_from_mass(md_element_t out_element[], const float in_mass[], size_t count);

// The per segment quantities of a protein backbone (see md_protein_backbone_data_t), computed into the caller's arrays
// of one entry per segment (capacity at least backbone->segment.count). They do not allocate.

// Secondary structure (DSSP like) from the coordinates of a frame
bool md_util_backbone_secondary_structure_infer(md_secondary_structure_t secondary_structures[], size_t capacity, const vec3_t* xyz, const md_unitcell_t* cell, const md_protein_backbone_data_t* backbone);

// Backbone angles (phi, psi) from the coordinates of a frame
bool md_util_backbone_angles_compute(md_backbone_angles_t backbone_angles[], size_t capacity, const vec3_t* xyz, const md_unitcell_t* cell, const md_protein_backbone_data_t* backbone);

// Ramachandran type (General / Glycine / Proline / Preproline) from the residue names of sys->protein_backbone
bool md_util_backbone_ramachandran_classify(md_ramachandran_type_t ramachandran_types[], size_t capacity, const struct md_system_t* sys);

// THE BACKBONE OF A STATE. The angles and the secondary structure depend on the coordinates, so they belong to the
// state they were computed from and live in its attributes (md_system_state_t.attributes), one value per segment of
// sys->protein_backbone:
//   MD_BACKBONE_ANGLE_PATH                 F32 x 2 (phi, psi), radians
//   MD_BACKBONE_SECONDARY_STRUCTURE_PATH   I32, md_secondary_structure_t
// The same paths below a run (in the system's attributes) hold them for every frame of the run, with the frame as
// the first axis, and md_system_extract_frame slices a frame of them into a state when asked for them.
// A state carries them when its producer computed them, like the velocities of a TRR frame: a consumer handles their
// absence.
#define MD_BACKBONE_ANGLE_PATH                  "backbone/angle"
#define MD_BACKBONE_SECONDARY_STRUCTURE_PATH    "backbone/secondary_structure"

// The backbone of a state, NULL when it carries none or what it carries does not fit the backbone of sys
const md_backbone_angles_t*     md_util_state_backbone_angles(const md_system_state_t* state, const struct md_system_t* sys);
const md_secondary_structure_t* md_util_state_secondary_structure(const md_system_state_t* state, const struct md_system_t* sys);

// Storage in the state to write the backbone into, created when missing (or of another size) and reused otherwise,
// so that a producer which writes every frame does not allocate every frame. Each call marks the attribute as
// changed (md_attributes_touch): call it when about to write. NULL for a state without an allocator (a view) or a
// system without a protein backbone.
md_backbone_angles_t*     md_util_state_backbone_angles_write(md_system_state_t* state, const struct md_system_t* sys);
md_secondary_structure_t* md_util_state_secondary_structure_write(md_system_state_t* state, const struct md_system_t* sys);

// Computes the backbone angles and secondary structure from the state's own coordinates and stores them in it.
// False when the state has no coordinates, or the system no protein backbone.
bool md_util_state_backbone_compute(md_system_state_t* state, const struct md_system_t* sys);

void md_util_infer_covalent_bonds(md_bond_data_t* out_bond, const md_system_state_t* state, const md_system_t* sys, struct md_allocator_i* alloc);

// Computes the covalent bonds based from a heuristic approach, uses the covalent radius (derived from element) to determine the appropriate bond
// length. atom_res_idx is an optional parameter and if supplied, it will limit the covalent bonds to only within the same or adjacent residues.
void md_util_system_infer_covalent_bonds(md_system_t* sys, const md_system_state_t* state);

// Marks every bond between a metal and a non-metal MD_BOND_FLAG_COORDINATE, whatever its origin, and returns how many
// it marked. md_util_system_infer does it for each system; a bond added afterwards (by hand) wants it too.
size_t md_util_system_infer_coordination(md_system_t* sys);

// Grow a mask by bonds up to a certain extent (counted as number of bonds from the original mask)
// Viable mask is optional and if supplied, it will limit the growth to only within the viable mask
void md_util_mask_grow_by_bonds(struct md_bitfield_t* mask, const struct md_system_t* sys, size_t extent, const struct md_bitfield_t* viable_mask);

// Grow a mask by radius (in Angstrom)
// Viable mask is optional and if supplied, it will limit the growth to only within the viable mask
void md_util_mask_grow_by_radius(struct md_bitfield_t* mask, const md_system_state_t* state, double radius, const struct md_bitfield_t* viable_mask);

// Infer rings formed by covalent bonds
bool md_util_system_infer_rings(md_system_t* sys);

// Identify isolated structures by covalent bonds.
// For coarse grained systems (any particle a bead, md_system_is_coarse_grained), beads the bonds leave disconnected are
// joined through the component hierarchy instead: each bead to its component's backbone (or first) bead, and those
// anchors along consecutive polymer components. No coordinates are consulted and sys->bond is not modified.
bool md_util_system_infer_structures(md_system_t* sys);

// Identify atom types within the system
void md_util_system_infer_atom_types(md_system_t* sys, const str_t atom_labels[]);

// Applies the predefined atom type tables (coarse grained beads) to atoms that already have a type,
// for loaders that know more about their atoms than a name does (a tpr carries the mass and the LJ
// parameters of every bead). Only types without an element (z == 0) are considered. The table gives
// the atoms their role (backbone, side chain), the components their kind (amino acid, water) where
// they have none, and the types their particle kind where the loader did not say. Mass and radius
// are left alone: they are the loader's.
void md_util_system_augment_atom_types(md_system_t* sys);

// Attempts to generate missing data such as covalent bonds, chains, secondary structures, backbone angles etc.
// Infers the derivable parts of a system (covalent bonds, rings, structures, backbones, hydrogen
// bond roles) from the supplied state, and records that state as sys->reference.
//
// The state is an explicit parameter rather than read off the system so that the recorded
// reference is by construction the input which produced the topology. Requires coordinates for
// anything geometric; a state without them (a topology only format such as PSF) still yields the
// purely graph derived parts, provided the bonds were supplied by the format.
bool md_util_system_infer(struct md_system_t* sys, const md_system_state_t* state, md_infer_flags_t flags);


// Computes an array of distances between two sets of coordinates in a periodic domain (cell)
// out_dist:  Output array of distances, must have length of (num_a * num_b)
// coord_a:   Array of coordinates (a)
// num_a:     Length of coord_a
// coord_b:   Array of coordinates (b)
// num_b:     Length of coord_b
// cell:      Periodic boundary cell
void md_util_distance_array(float* out_dist_arr, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell);

// The minimum (maximum) distance between each group of points of a and the points of b, under the minimum image
// convention of cell (NULL, or a cell without periodic axes: plain euclidean). The minimum image is exact in triclinic
// cells too, also where distances reach past half the cell (see MINIMUM AND MAXIMUM DISTANCE BETWEEN SETS in md_util.c).
// Group g is coord_a[a_offsets[g] .. a_offsets[g + 1]), a_offsets holds num_groups + 1 entries.
// out_dist:   The distance per group, 0 for an empty group or an empty b
// out_idx_a:  (optional) Per group, the index into coord_a of the nearest (farthest) pair, -1 where there is none
// out_idx_b:  (optional) Per group, the index into coord_b of the pair, -1 where there is none
// Each group starts from the pair of the one before, so groups which follow each other in space are cheaper.
void md_util_min_distance_groups(float* out_dist, int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, const size_t* a_offsets, size_t num_groups, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell);
void md_util_max_distance_groups(float* out_dist, int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, const size_t* a_offsets, size_t num_groups, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell);

// The same for a single group. The minimum is FLT_MAX and the maximum 0 when a or b is empty, and the indices are then
// left unwritten (for the maximum also when it is 0).
float md_util_min_distance(int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell);
float md_util_max_distance(int64_t* out_idx_a, int64_t* out_idx_b, const vec3_t* coord_a, size_t num_a, const vec3_t* coord_b, size_t num_b, const md_unitcell_t* cell);

void md_util_min_image_vec3(vec3_t* in_out_dx, size_t count, const md_unitcell_t* cell);
void md_util_min_image_vec4(vec4_t* in_out_dx, size_t count, const md_unitcell_t* cell);

// Applies periodic boundary conditions to coordinates
bool md_util_pbc(vec3_t* in_out_xyz, const int32_t* in_idx, size_t count, const md_unitcell_t* cell);
bool md_util_pbc_vec4(vec4_t* in_out_xyzw, size_t count, const md_unitcell_t* cell);

// Applies periodic boundary conditions to all coordinates in a systems state (convenience function)
bool md_util_system_pbc(md_system_state_t* state);

// Unwraps a single structure by walking the parent hierarchy md_util_system_infer_structures
// already computed. Linear, allocation free, and does not consult bond connectivity.
void md_util_unwrap_structure(md_system_state_t* state, const md_structure_t* structure);

// Unwraps all structures in a system
void md_util_unwrap_system(md_system_state_t* state, const md_system_t* sys);

// Batch deperiodize a set of coordinates (vec4) with respect to a given reference
bool md_util_deperiodize_vec4(vec4_t* xyzw, size_t count, vec3_t ref_xyz, const md_unitcell_t* cell);

// Converts absolute coordinates (in unitcell) to relative coordinates (with respect to given reference point)
void md_util_convert_to_relative_coordinates_vec4(vec4_t* in_out_rel_xyzw, vec3_t ref_xyz, size_t count, const md_unitcell_t* cell);

// PERIODIC IMAGE SELECTION
//
// These let the estimator choose each point's periodic image, using its own objective as the
// criterion. Prefer them over unwrapping followed by estimation: unwrap has to commit to an image
// assignment before it knows what is being measured, and it commits using topology, which is
// unrelated to the measurement. These do not need topology at all, which is also why they work for
// sparse selections and for points which are not atoms.
//
// @NOTE: md_util_unwrap_* remains the right tool for making a whole connected structure whole - for
// rendering or export - where there is no estimator and no reference to align against, and where the
// structure may be more extended than half a cell. That is the only case it is for.

// Places a set of points in mutually consistent images and reports its centre.
//
// Seeded from the circular mean (which requires no image assignment, so it cannot be skewed by the
// images the input arrived in) and then alternated to a fixed point: place every point in the image
// nearest the centre, recompute the centre, repeat. Settles in one or two passes.
//
// The set is left in the image it arrived in: the points are made mutually consistent and the
// cluster as a whole is NOT moved into the reference cell, so the reported centre can be combined
// with coordinates that were not passed in. Call md_util_pbc afterwards if the cell is wanted.
// A set straddling a boundary has no image to preserve and may come back on either side.
//
// Correct while the set spans less than half a cell. Beyond that no per point image choice can
// represent it and the structure has to be unwrapped topologically instead.
bool md_util_deperiodize_self_vec4(vec4_t* in_out_xyzw, size_t count, const md_unitcell_t* cell, vec3_t* out_com);

// Optimal rigid rotation between a reference set and a target set under periodic boundaries.
//
// Solves for the rotation, the target centre AND each target point's periodic image together:
//     min over R, c, n_k  of  sum_k w_k | R (q_k + A n_k - c) - (p_k - ref_com) |^2
// Given R and c the best n_k is the image nearest the predicted position, per point and in closed
// form; given the n_k the best R and c is ordinary Kabsch. Alternating the two decreases the
// objective monotonically and terminates when no point changes image.
//
// out_rot maps the TARGET frame onto the REFERENCE frame: R * (q - out_com) ~= p - ref_com.
// out_trg_xyzw optionally receives the placed target points (weights preserved); pass NULL if only
// the transform is wanted, and it may alias trg_xyzw.
//
// out_com and the placed points are expressed in the image the TARGET arrived in, not folded into
// the reference cell - see md_util_deperiodize_self_vec4. ref_com is used only as the reference
// set's own centre, so its image does not matter.
//
// Returns the largest per point residual, which is the margin to an ambiguous image choice. Small
// against half a cell means the assignment is nowhere near flipping; approaching it means this
// frame's alignment should not be trusted.
float md_util_optimal_rotation_pbc_vec4_iter(mat3_t* out_rot, vec3_t* out_com, vec4_t* out_trg_xyzw,
                                             const vec4_t* ref_xyzw, vec3_t ref_com,
                                             const vec4_t* trg_xyzw, size_t count,
                                             const md_unitcell_t* cell, int max_iter, float tol);

mat3_t md_util_optimal_rotation_rel_vec4(const vec4_t* ref_rel_xyzw, const vec4_t* trg_rel_xyzw, size_t count);

// Computes the minimum axis aligned bounding box for a set of points with a given radius
// Indices are optional and are used to select a subset of points, the count dictates the number of elements to process
void md_util_aabb_compute     (float out_ext_min[3], float out_ext_max[3], const vec3_t* in_xyz, const float* in_r, const int32_t* in_idx, size_t count);
void md_util_aabb_compute_vec4(float out_ext_min[3], float out_ext_max[3], const vec4_t* in_xyzr, const int32_t* in_idx, size_t count);

// Computes an object oriented bounding box based on the PCA of the provided points (with optional radius)
void md_util_oobb_compute     (float out_rotation[3][3], float out_ext_min[3], float out_ext_max[3], const vec3_t* in_xyz, const float* in_r, const int32_t* in_idx, size_t count, const md_unitcell_t* cell);
void md_util_oobb_compute_vec4(float out_rotation[3][3], float out_ext_min[3], float out_ext_max[3], const vec4_t* in_xyzr, const int32_t* in_idx, size_t count, const md_unitcell_t* cell);

// Computes the center of mass for a set of points with a given weight
// xyz:         Packed coordinates
// w:           Array of weights (optional): set as NULL to use equal weights
// indices:     Array of indices (optional): indices into the arrays (xyz,w)
// count:       Length of all arrays
// unit_cell:   The unit_cell of the system [Optional]
vec3_t md_util_com_compute(const vec3_t* in_xyz, const float* in_w, const int32_t* in_idx, size_t count, const md_unitcell_t* cell);
vec3_t md_util_com_compute_vec4(const vec4_t* in_xyzw, const int32_t* in_idx, size_t count, const md_unitcell_t* cell);

// Computes the similarity between two sets of points with given weights.
// One of the sets is rotated and translated to match the other set in an optimal fashion before the similarity is computed.
// The rmsd is the root mean squared deviation between the two sets of aligned vectors.
// xyz:     Packed coordinate arrays [2]
// com:     Center of mass [2] (xyz0), (xyz1)
// w:       Array of weights (optional): set as NULL to use equal weights
// count:   Length of all arrays (xyz0, xyz1, w)
double md_util_rmsd_compute(const vec3_t* const in_xyz[2], const float* const in_w[2], const int32_t* const in_idx[2], size_t count, const vec3_t in_com[2]);
double md_util_rmsd_compute_vec4(const vec4_t* const in_xyzw[2], const int32_t* const in_idx[2], size_t count, const vec3_t in_com[2]);

// Computes linear shape descriptor weights (linear, planar, isotropic) from a covariance matrix
vec3_t md_util_shape_weights(const mat3_t* covariance_matrix);

// Perform linear interpolation of supplied coordinates
// out_xyz:     Destination, packed
// in_xyz:      Sources [2], packed
// count:       Count of coordinates (this implies that all coordinate arrays must be equal in length)
// unit_cell:   The unit_cell of the system [Optional]
// t: interpolation factor (0..1)
// @NOTE: Input and output must be padded to a multiple of eight atoms: it works on full simd width
bool md_util_interpolate_linear(vec3_t* out_xyz, const vec3_t* const in_xyz[2], size_t count, const md_unitcell_t* cell, float t);

// Perform cubic interpolation of supplied coordinates
// out_xyz:     Destination, packed
// in_xyz:      Sources [4], packed
// count:       Count of coordinates (this implies that all coordinate arrays must be equal in length)
// unit_cell:   The unit_cell of the system [Optional]
// t:           Interpolation factor (0..1)
// s:           Scaling factor (0..1), 0 is jerky, 0.5 is catmul rom, 1.0 is silky smooth
// @NOTE: Input and output must be padded to a multiple of eight atoms: it works on full simd width
bool md_util_interpolate_cubic_spline(vec3_t* out_xyz, const vec3_t* const in_xyz[4], size_t count, const md_unitcell_t* cell, float t, float s);

// Spatially sorts the input positions according to morton order. This makes it easy to create spatially coherent clusters, just select ranges within this space.
// There are some larger jumps within the morton order as well, so when creating clusters from consecutive ranges, this should be considered as well.
// The result (source_indices) is an array of remapping indices. It is assumed that the user has reserved space for this.
void md_util_sort_spatial(uint32_t* source_indices, const vec3_t* xyz, size_t count);

// Spatially sorts the input positions according to morton order. This makes it easy to create spatially coherent clusters, just select ranges within this space.
// There are some larger jumps within the morton order as well, so when creating clusters from consecutive ranges, this should be considered as well.
// The result (source_indices) is an array of remapping indices. It is assumed that the user has reserved space for this.
void md_util_sort_spatial_xyz(uint32_t* source_indices, const float* xyz, size_t stride_in_bytes, size_t count);

// Sort array of uint32_t in place using radix sort
void md_util_sort_radix_inplace_uint32(uint32_t* data, size_t count);

// Sort array of uint32_t by producing a remapping array of source indices
// The source_indices represents the indices of the sorted array, i.e. source_indices[0] is the index of the smallest element in data
void md_util_sort_radix_uint32(uint32_t* out_indices, const uint32_t* key, size_t count);

#ifdef __cplusplus
}
#endif
