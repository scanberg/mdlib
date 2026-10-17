#pragma once

//  MD GL - Molecule Dynamics GL
//  Copyright 2020 Robin Skånberg
// 
//  LICENSE: see license file

/*
@NOTE(Robin): This is a clusterf*ck API and I'm too lazy to change it by now since it will at some point be superseded by something that supports next gen rendering APIs.
*/

#include <stdint.h>
#include <stdbool.h>

#include <core/md_str.h>
#include <core/md_vec_math.h>
#include <md_types.h>

// Forward declarations
struct md_system_t;
struct md_unitcell_t;

#ifdef __cplusplus
extern "C" {
#endif

// ### MD GL ###
// Initialize the GL module and its resources

void md_gl_initialize(void);
void md_gl_shutdown(void);

/*

In order to interface with your application, the shaders have a customizable output.
A function called 'write_fragment' is called in the fragment shader to finalize and write the fragments to the respective outputs.
This part is customizable by providing a snipped which controls what parameters to write and to which channels of the framebuffer.

If no custom_shader_snippet is supplied to md_gl_context_init, then a default_shader_snippet is used which is listed below:

layout(location = 0) out vec4 out_color;

void write_fragment(vec3 view_coord, vec3 view_vel, vec3 view_normal, vec4 color, uint element_index) {
   out_color  = color;
}
*/

// The shaders are customizable in their output as mentioned in the comment above.
// You could have multiple versions of shaders each compiled for a specific output.
typedef struct {
    uint32_t id;
} md_gl_shaders_t;

md_gl_shaders_t md_gl_shaders_create(str_t str);
void md_gl_shaders_destroy(md_gl_shaders_t handle);


// ### GL MOLECULE ###
//
// The gl molecule represents the state of a single molecule which can be drawn.
// The required fields depends on the visual representation applied to the molecule, but to cover all cases the following fields are used:
//  - atom
//      - x      : f32
//      - y      : f32
//      - z      : f32
//      - radius : f32
//      - flags  : u8
//  - residue
//      - atom_range : u32[2]
//  - bond
//      - bond : md_bond_t = u32[2]
//  - backbone
//      - backbone_atoms      : md_backbone_atoms_t = u32[4] (N, CA, C, O)
//      - secondary_structure : md_secondary_structure_t = u8[4]

typedef struct {
    uint32_t id;
} md_gl_mol_t;

// Secondary Structure Type
// This is different from md_secondary_structure_t to allow interpolated values
// It essentially controls the blend (barycentric) between the geometrical representations of coil, helix, and sheet
// coil is implied when both helix and sheet are zero
typedef struct {
    float helix;  // 1.0 = full helix
    float sheet;  // 1.0 = full sheet
} md_gl_secondary_structure_t;

static inline md_gl_secondary_structure_t md_gl_secondary_structure_convert(md_secondary_structure_t ss) {
    md_gl_secondary_structure_t gl_ss = {0};
    switch (ss) {
        case MD_SECONDARY_STRUCTURE_HELIX_310:
        case MD_SECONDARY_STRUCTURE_HELIX_ALPHA:
        case MD_SECONDARY_STRUCTURE_HELIX_PI:
            gl_ss.helix = 1.0f;
            gl_ss.sheet = 0.0f;
            break;
        case MD_SECONDARY_STRUCTURE_BETA_SHEET:
        case MD_SECONDARY_STRUCTURE_BETA_BRIDGE:
            gl_ss.helix = 0.0f;
            gl_ss.sheet = 1.0f;
            break;
        case MD_SECONDARY_STRUCTURE_COIL:
        case MD_SECONDARY_STRUCTURE_TURN:
        case MD_SECONDARY_STRUCTURE_BEND:
        case MD_SECONDARY_STRUCTURE_UNKNOWN:
        default:
            gl_ss.helix = 0.0f;
            gl_ss.sheet = 0.0f;
            break;
    }
    return gl_ss;
}


md_gl_mol_t md_gl_mol_create(const struct md_system_t* mol);
// A molecule of bare particles without bonds or a backbone, for geometry which stands in for parts of a system (the
// cartoon of lipids and nucleic acids): set their positions, or derive them from the atoms of the system
// (md_gl_mol_set_derived), and their radii, and draw them as spheres (MD_GL_REP_SPACE_FILL, a radius of 0 is not
// drawn), as capsules between them (the weak bonds of a representation, MD_GL_REP_WEAK_BONDS) and as tubes through them
// (md_gl_mol_set_tubes, MD_GL_REP_TUBES). Positions and radii start at 0.
md_gl_mol_t md_gl_mol_create_particles(uint32_t count);
void md_gl_mol_destroy(md_gl_mol_t handle);

// ### DERIVED PARTICLES ###
// The particles of a molecule (md_gl_mol_create_particles) derived on the GPU from the atoms of another molecule, their
// source: particle i is the weighted sum of the atoms atom[offset[i] .. offset[i + 1]) of the source, with the weights
// weight[offset[i] .. offset[i + 1]) (their mean when the weights sum to 1), and its velocity the same sum of theirs.
// Each atom is taken in its periodic image nearest to the atom ref[i] of the source, so the particles which share a
// reference stay whole and together across the periodic boundary of the cell (the points of one molecule).
//  offset: [particle count + 1], atom and weight: [offset[particle count]], ref: [particle count]
// Atoms outside of the source are left out; a particle of no atoms is at the origin. Replaces a previous derivation.
// The source must outlive it.
void md_gl_mol_set_derived(md_gl_mol_t mol, md_gl_mol_t source, const uint32_t* offset, const uint32_t* atom, const float* weight, const uint32_t* ref);

// Computes the positions and velocities of the derived particles from the positions and velocities of the atoms of their
// source on the GPU (md_gl_mol_set_atom_position, md_gl_mol_compute_velocity, md_gl_mol_zero_velocity): call it after
// those change.
//  cell:       (optional) the unit cell of the positions, for the periodic images. Without it there are none.
//  source_xyz: (optional) the positions of the atoms of the source on the CPU, the same as on the GPU. They place the
//              tubes (md_gl_mol_set_tubes) for their culling and level of detail; without them the tubes are drawn in
//              full detail and never culled.
void md_gl_mol_compute_derived(md_gl_mol_t mol, const struct md_unitcell_t* cell, const vec3_t* source_xyz);

// ### TUBES ###
// Splines through the particles of a molecule (md_gl_mol_create_particles), drawn by MD_GL_REP_TUBES: the backbone of a
// molecule of particles, as the protein backbone is that of a system (MD_GL_REP_CARTOON), and computed the same way
// (uniform cubic B-splines, which pass near rather than through their control points, and smooth them: or through
// them, within about a tenth of what the B-spline is off by, for a tube with through set). Tube t has the control points
// [offset[t], offset[t + 1]), each a particle, in order.
// The profile is a stadium (a 2D capsule, as the cartoon's): a segment of half length flat widened by radius, two half
// circles joined by flat sides, a circle for a flat of 0. Both follow the control points smoothly along the tube
// (Ångström, at most 8 each). The flat sides face the particle face of a control point, where it has one: the backbone
// of a nucleic acid faces its bases. Otherwise the profile turns with the curve as the cartoon of a protein does.
// The ends are round: the Minkowski sum of the segment and a half sphere of the radius there, which also closes a tube
// next to a part of it that is hidden. A tube of fewer than two control points is not drawn.
// A part of a tube belongs to its control point: it picks as it (md_gl_picking_offset_t), takes its colour from the
// particle and is shown when the particle is (the alpha of its colour, the atom mask), so the particles of the control
// points can have colours and flags while a radius of 0 keeps them from being drawn as spheres.
typedef struct md_gl_tubes_t {
    uint32_t        count;      // Tubes
    const uint32_t* offset;     // [count + 1] The control points of each tube in the arrays below
    const uint32_t* particle;   // The particle of each control point
    const float*    radius;     // The radius of the profile, half its thickness
    const float*    flat;       // (optional) Half the length of the flat sides of the profile, NULL for 0: round tubes
    const uint32_t* face;       // (optional) The particle the flat sides face, UINT32_MAX for none
    const uint8_t*  through;    // (optional) [count] Per tube: it passes (nearly) through its control points rather than near them
} md_gl_tubes_t;

// Replaces the previous tubes; NULL or a count of 0 clears them
void md_gl_mol_set_tubes(md_gl_mol_t mol, const md_gl_tubes_t* tubes);

// ### SLABS ###
// Flat rounded polygons through the particles of a molecule (md_gl_mol_create_particles), drawn by MD_GL_REP_SLABS: the
// rings of a cartoon, the sugars and bases of nucleic acids. Slab s has the corners particle[offset[s] .. offset[s + 1]),
// 2 to 7 particles in order around it. It lies in the plane that fits its corners best, through their centroid, and is
// the polygon of the corners projected on that plane, made solid by a sphere: a rim of capsules along its edges and a
// plate between them (md_gl_draw_op_t args.slabs). The polygon is taken to be convex, as a ring is. A slab of two
// corners is the capsule between them, of the radius of the rims: a stick which joins two rings as smoothly as their
// own edges do (the sugar and the base of a nucleotide).
// A slab is raycast, exactly: it picks as md_gl_picking_offset_t ring_base + its index, takes the colour of its first
// corner and is shown when all its corners are. Replaces the previous slabs; a count of 0 clears them.
void md_gl_mol_set_slabs(md_gl_mol_t mol, uint32_t count, const uint32_t* offset, const uint32_t* particle);

// ### Set molecule data fields ###
void md_gl_mol_set_index_base   (md_gl_mol_t handle, uint32_t atom_index_base, uint32_t bond_index_base);
// Packed xyz, the buffers' own layout. Setting positions keeps the current ones as the previous.
void md_gl_mol_set_atom_position(md_gl_mol_t handle, uint32_t atom_offset, uint32_t atom_count, const vec3_t* xyz);
void md_gl_mol_set_atom_velocity(md_gl_mol_t handle, uint32_t atom_offset, uint32_t atom_count, const vec3_t* xyz);
void md_gl_mol_set_atom_radius  (md_gl_mol_t handle, uint32_t atom_offset, uint32_t atom_count, const float* radius,  uint32_t byte_stride);
void md_gl_mol_set_atom_flags   (md_gl_mol_t handle, uint32_t atom_offset, uint32_t atom_count, const uint8_t* flags, uint32_t byte_stride);

// Call this function after setting new atomic positions to update velocities
// It will compute a new velocity as the difference between new and old atomic positions
// pbc_ext is an optional parameter which corresponds to the extent of the periodic boundry conditions.
// If supplied, it will make sure to deperiodize the atomic positions before computing the velocity,
// which is crucial in order to resolve an accurate velocity in a periodic domain.
void md_gl_mol_compute_velocity(md_gl_mol_t handle, const float pbc_ext[3]);

// Clear the velocity to zero.
void md_gl_mol_zero_velocity(md_gl_mol_t handle);

// Changing the bonds clears the bond orders: set them again after (md_gl_mol_set_bond_orders)
void md_gl_mol_set_bonds(md_gl_mol_t mol, uint32_t offset, uint32_t count, const struct md_atom_pair_t* bond_pairs, uint32_t byte_stride);

// The orders of the bonds, for the licorice and ball and stick representations with bond_orders set, from the topology
// of the system: the order and the aromatic and delocalized flags of its bonds (md_chem_perceive), its rings and their
// flags and the hybridization of its atoms. The bonds of the system must be those of the molecule
// (md_gl_mol_set_bonds), in the same order. md_gl_mol_create sets them.
//   - double bond: two thinner sticks in the plane of its sp2 centers (the plane of their substituents)
//   - triple bond: three thinner sticks, perpendicular to the view
//   - aromatic ring (MD_RING_FLAG_AROMATIC): the sticks of its bonds and a thin solid circle inside it, in its plane,
//     one in each of two fused rings. The circle picks as the ring (md_gl_picking_offset_t).
//   - delocalized bond (carboxylate, nitro, guanidinium, phosphate, ...): the stick and a thin dashed line towards
//     the next bond of the group. So is an aromatic bond in no circle (in a ring larger than MD_RING_MAX_SIZE).
// A bond whose order is unknown is drawn as a single bond.
void md_gl_mol_set_bond_orders(md_gl_mol_t mol, const struct md_system_t* sys);

// The gaps of the chains, from the chain breaks of the system (MD_COMPONENT_FLAG_BREAK_AFTER, CHAIN BREAKS in md_util.h):
// for the representations with chain_gaps set, a straight dashed line over each break, from the residue before it to
// the one after it. Where the missing residues went is not known, so the line claims no more than that the chain goes
// on. The cartoon and ribbons draw it between the ends of the cartoon of the chain on either side (near the CA of the
// residue, where the spline ends) for the breaks of protein chains; licorice and ball and stick draw it between the
// atoms a link would join (md_util_chain_break_atoms: C and N, O3' and P, the backbone beads). Its dashes are of a
// length of their own (Ångström), as many as fit, so a long gap reads as a gap and not as a bond. It takes the colors
// of its two ends (the atoms it is drawn between), is shown when both of them are, and is drawn at any length: its
// ends are where the atoms are, not their minimum image. It picks as the component before the break
// (md_gl_picking_offset_t gap_base). The system must be the one of the molecule. md_gl_mol_create sets them.
void md_gl_mol_set_chain_gaps(md_gl_mol_t mol, const struct md_system_t* sys);

void md_gl_mol_set_backbone_secondary_structure(md_gl_mol_t mol, uint32_t offset, uint32_t count, const md_gl_secondary_structure_t* secondary_structure, uint32_t byte_stride);

// The backbone (ribbons, cartoon) orients its cross sections with temporal coherence: each draw continues from
// the orientation computed by the previous one, which keeps segments from flipping while the structure moves.
// Call this when the positions jump discontinuously (e.g. seeking to a distant frame) so the next draw starts
// from the structure alone.
void md_gl_mol_reset_backbone_history(md_gl_mol_t mol);

/*
 *  REPRESENTATIONS
 *  Interface for creating visual representations for molecules
 *  A representations is just a set of colors for each atom in the molecule
 */

typedef struct {
    uint32_t id;
} md_gl_rep_t;

md_gl_rep_t md_gl_rep_create(md_gl_mol_t mol);
void        md_gl_rep_destroy(md_gl_rep_t rep);

void        md_gl_rep_set_atom_colors(md_gl_rep_t rep, uint32_t offset, uint32_t count, const uint32_t* color, uint32_t byte_stride);

// Weak bonds (e.g. hydrogen bonds) of the representation, drawn by MD_GL_REP_WEAK_BONDS as dashed cylinders between the
// atoms. Unlike the covalent bonds of the molecule they belong to the representation, and are typically set every frame.
// Replaces the previous set; count 0 clears it. The arrays match md_hbond_set_t (hydrogen, acceptor, strength).
//  atom_a, atom_b: [count] atom indices of the bonds. Bonds with an atom outside of the molecule are left out.
//  weight:         [count] (optional, NULL for 1) scales the radius of each bond, so a bond can fade out continuously. A
//                  weight of 0 or less hides it.
// As for the other representations, a bond is shown only when both of its atoms are (alpha of the colors, atom mask).
void        md_gl_rep_set_weak_bonds(md_gl_rep_t rep, uint32_t count, const uint32_t* atom_a, const uint32_t* atom_b, const float* weight);

/*
 *  DRAW
 *  Interface for defining representations of molecules
 *
 */

typedef enum {
    MD_GL_REP_SPACE_FILL,
    MD_GL_REP_LICORICE,
    MD_GL_REP_BALL_AND_STICK,
    MD_GL_REP_RIBBONS,
    MD_GL_REP_CARTOON,
    MD_GL_REP_WEAK_BONDS,   // The weak bonds of the representation (md_gl_rep_set_weak_bonds), dashed
    MD_GL_REP_TUBES,        // The tubes of a molecule of particles (md_gl_mol_set_tubes)
    MD_GL_REP_SLABS,        // The slabs of a molecule of particles (md_gl_mol_set_slabs)
} md_gl_rep_type_t;

// Controls how the color of the bonds are calculated in the licorice and ball-and-stick representations
typedef enum {
	MD_GL_BOND_MODE_NEAREST,
	MD_GL_BOND_MODE_SMOOTH,
    MD_GL_BOND_MODE_UNIFORM,
} md_gl_bond_mode_t;

typedef struct md_gl_draw_op_t {
    md_gl_rep_type_t type;

    union {
        struct {
            float radius_scale;
        } space_fill;

        struct {
            float radius;
			md_gl_bond_mode_t color_mode;
            float sharpness; // 0.0 = fully smooth, 1.0 = fully sharp
			uint32_t uniform_color; // Used if color_mode is MD_GL_BOND_MODE_UNIFORM
            bool bond_orders;       // Show the orders of the bonds (md_gl_mol_set_bond_orders)
            bool chain_gaps;        // Show the gaps of the chains (md_gl_mol_set_chain_gaps)
        } licorice;

        struct {
            float ball_scale;
            float stick_radius;
            md_gl_bond_mode_t color_mode;
            float sharpness;
            uint32_t uniform_color; // Used if color_mode is MD_GL_BOND_MODE_UNIFORM
            bool bond_orders;       // Show the orders of the bonds (md_gl_mol_set_bond_orders)
            bool chain_gaps;        // Show the gaps of the chains (md_gl_mol_set_chain_gaps)
        } ball_and_stick;

        struct {
            float width_scale;
            float thickness_scale;
            bool  chain_gaps;       // Show the gaps of the chains (md_gl_mol_set_chain_gaps)
        } ribbons;

        struct {
            float coil_scale;
            float sheet_scale;
            float helix_scale;
            bool  raycast;          // Ray cast the profile (in closed form, a quad per segment of the spline) rather than draw a mesh of it
            bool  chain_gaps;       // Show the gaps of the chains (md_gl_mol_set_chain_gaps)
        } cartoon;

        // Dashed cylinders between the atoms: dash_count dashes with flat ends, the first starting at the first atom
        // and the last ending at the second, each covering the fraction dash_fill of the period of the pattern (1 is a
        // solid cylinder). The dashes are as many on every bond, so they scale with its length.
        struct {
            float radius;           // Ångström, scaled by the weight of each bond
            uint32_t dash_count;    // At least 1
            float dash_fill;        // (0, 1]
            md_gl_bond_mode_t color_mode;
            float sharpness;
            uint32_t uniform_color; // Used if color_mode is MD_GL_BOND_MODE_UNIFORM
        } weak_bonds;

        struct {
            float radius_scale;     // Of the profiles at the control points
            float band_scale;       // Of the radius of the profiles with flat sides (bands), on top of radius_scale: their
                                    // thickness, the flat sides keep their length. 0 or less is taken as 1
            bool  raycast;          // As for the cartoon
        } tubes;

        struct {
            float radius;           // Ångström, of the rim: half the thickness of the slab, and how far out of its corners it reaches
            float plate;            // The half thickness of the plate within the rim, as a fraction of the radius: 1 is flat
        } slabs;
    } args;

    md_gl_rep_t rep;
    const float* model_matrix;              // Column major float[4][4]
} md_gl_draw_op_t;

// Picking offsets: each renderable element writes base + its own index to the picking target, so the caller
// reserves one disjoint range per kind. Atoms (spacefill, licorice ends) write their atom index, bonds (licorice
// middle) their bond index, and the backbone representations (cartoon, ribbons) their protein backbone segment
// index (md_system_t::protein_backbone.segment). A cartoon is not made of atoms: it writes the segment it stems
// from, never an atom of it. A tube writes the index of its control point among those of all tubes (md_gl_mol_set_tubes)
// as a backbone, a slab its index (md_gl_mol_set_slabs) as a ring. A weak bond is not a bond of the molecule: it writes the atom index of its closer end.
// The circle of an aromatic ring writes the index of the ring in md_system_t::ring (all of its rings, so the range
// is md_index_data_num_ranges(&sys->ring) long, of which only the aromatic ones are drawn). The gap of a chain
// (md_gl_mol_set_chain_gaps) writes the index of the component before its break (a range of the system's component
// count, of which only those with MD_COMPONENT_FLAG_BREAK_AFTER are drawn).
typedef struct md_gl_picking_offset_t {
    uint32_t atom_base;
    uint32_t bond_base;
    uint32_t backbone_base;
    uint32_t ring_base;
    uint32_t gap_base;
} md_gl_picking_offset_t;

typedef struct md_gl_draw_args_t {
    md_gl_shaders_t shaders;

    struct {
        size_t count;
        md_gl_draw_op_t* ops;
    } draw_operations;

    struct {
        // @NOTE: Matrices are column major float[4][4]
        const float* view_matrix;
        const float* proj_matrix;

        // [Optional] These fields are only used to compute screen space velocity
        const float* prev_view_matrix; 
        const float* prev_proj_matrix;
    } view_transform;

    md_gl_picking_offset_t picking_offset;

    uint32_t atom_mask;
	float max_bond_length;	// Maximum bond length for rendering bonds, covalent and weak (in Ångström)
} md_gl_draw_args_t;

bool md_gl_draw(const md_gl_draw_args_t* args);

#ifdef __cplusplus
}
#endif
