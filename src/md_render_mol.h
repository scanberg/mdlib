/*
md_render_mol.h

The molecular layer over md_render: turns representations of a system into
md_render batches. A compiler, not a renderer -- it owns nothing with a
lifetime:

  * The caller's per-atom arrays (positions, velocities, radii, colours,
    visibility, selection) are read in place, shared with every other GPU
    consumer of the system (GTO evaluation, interpolation, analysis).

  * What is derived from the topology for drawing (bonds grouped by kind,
    bond-order sticks, ring polygons, chain gaps, backbone control points,
    bead tables) lives in ONE caller-owned block of GPU memory. It is built in
    two steps so that the system is analysed once and the caller decides where
    the memory comes from and how long it lives:

        md_render_mol_plan_t plan;
        md_render_mol_plan(&plan, sys, temp);                  // host, in a temp scope
        md_gpu_mem_t mem = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, plan.size);
        md_render_mol_topo_t topo;
        md_render_mol_write(s, &topo, &plan, mem.gpu, content); // the upload
        ...
        md_gpu_free(s, topo.base);                               // the only "destroy"

    The topology is a plain value: copy it, keep one per system. Rebuild when
    sys->topology_version differs from topo.topology_version.

  * A representation is a plain struct passed every frame. md_render_mol_emit
    is pure and host-only: it records no GPU work and appends batches, which
    the caller can inspect, change or add to before md_render_view_draw.

  * Clusters (md_render's first level of culling) come with the topology:
    the atoms and the bonds are stored in cluster order, and emit hands that
    order and the cluster offsets to md_render. Bonded structures are cut
    into whole residues, packed up to 64 atoms and never across molecules,
    once. Small molecules (solvent, ions: at most MD_RENDER_MOL_SMALL atoms)
    diffuse, so their clusters are cells of a grid, rebuilt on the GPU by
    md_render_mol_update when the caller asks:

        md_render_mol_t mol = md_render_mol_create(device, NULL);   // kernels, once
        md_gpu_addr_t   dyn = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, topo.dynamic_size).gpu;
        ...
        md_render_mol_upload_state(s, positions, &state, content);
        md_render_mol_update(mol, s, &topo, dyn, &atoms, MD_RENDER_MOL_UPDATE_REGROUP, content);

    When to regroup is the caller's call: on a new frame of a trajectory,
    not for the steps interpolated between two frames (the grouping of the
    first is still a good one for those). Skipping it is never wrong, only
    slower: the bounds follow the positions regardless, the clusters just
    grow as the molecules spread. Before the first regroup the small
    molecules are grouped in molecule order.

  * Radii come from the system's particle types: the topology holds each
    atom's type and radius (from sys->atom.type.radius when planned). A
    changed radius per type is a table upload and an update
    (MD_RENDER_MOL_UPDATE_RADII); a radius per atom (size by a property)
    is the caller's array, md_render_mol_atoms_t.radius.

Implemented: SPACEFILL, LICORICE, BALL_AND_STICK, BONDS over bonds grouped by
kind, solid (dashes per kind wait for capsule dashes in md_render). Declared
for the next phases: CARTOON, RIBBONS, lipid/nucleic glyphs, bond orders,
chain gaps, ring circles.
*/

#ifndef MD_RENDER_MOL_H
#define MD_RENDER_MOL_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#include <md_render.h>
#include <core/md_allocator.h>

struct md_system_t;
struct md_system_state_t;

#ifdef __cplusplus
extern "C" {
#endif

/* Bonds are drawn by class. The classes follow md_bond_kind_t, so a
   representation chooses which kinds it draws with a mask of classes. */
typedef enum md_render_mol_bond_class_t {
    MD_RENDER_MOL_BOND_COVALENT = 0,
    MD_RENDER_MOL_BOND_COORDINATE,
    MD_RENDER_MOL_BOND_ELASTIC,
    MD_RENDER_MOL_BOND_METALLIC,
    MD_RENDER_MOL_BOND_INTERACTION,
    MD_RENDER_MOL_BOND_HYDROGEN,
    MD_RENDER_MOL_BOND_CLASS_COUNT,
} md_render_mol_bond_class_t;

#define MD_RENDER_MOL_BONDS(c) (1u << (c))
#define MD_RENDER_MOL_BONDS_ALL ((1u << MD_RENDER_MOL_BOND_CLASS_COUNT) - 1u)

/* Molecules of at most this many atoms are "small": regrouped by position on
   every update instead of by the topology once */
#define MD_RENDER_MOL_SMALL 16

/* =========================================================================
   Kernels
   ========================================================================= */

/* The GPU work of the layer (md_render_mol_update) needs compiled kernels: a
   context holding those and nothing else. One per device, shared by every
   system and view. */
typedef struct md_render_mol* md_render_mol_t;

md_render_mol_t md_render_mol_create(md_gpu_device_t device, md_allocator_i* alloc);
void            md_render_mol_destroy(md_render_mol_t mol);

/* =========================================================================
   Topology: one caller-owned block per system
   ========================================================================= */

typedef struct md_render_mol_topo_t {
    md_gpu_addr_t base;            /* the block; md_gpu_free(stream, base) frees it    */
    uint64_t      size;
    uint64_t      topology_version;/* of the system it was planned from               */
    uint64_t      dynamic_size;    /* bytes of the per-frame block md_render_mol_update
                                      writes; 0: there is nothing to update          */
    uint32_t      atom_count;

    /* Particle types and radii, per atom */
    uint32_t      type_count;              /* entries a radius table has: sys->atom.type.count */
    md_gpu_addr_t atom_type;               /* uint16 [atom_count]                      */
    md_gpu_addr_t atom_radius;             /* float [atom_count]: the type's radius    */

    /* The atoms in cluster order: residue clusters first, then one cluster
       per grid cell (grid_cells of them, some empty) for the small molecules */
    md_gpu_addr_t atom_order;              /* uint32 [atom_count]                      */
    uint32_t      atom_cluster_count;
    md_gpu_addr_t atom_cluster_offset;     /* uint32 [atom_cluster_count + 1]          */

    /* Bonds, grouped by class: class c is bond_pairs[class_offset[c] ..
       class_offset[c + 1]), in cluster order: with the cluster of their
       first atom's (64 at most per cluster); the covalent bonds within small
       molecules last, by grid cell. Usable outside rendering too, but the
       order within a class changes with every update. */
    uint32_t      bond_count;
    uint32_t      bond_class_offset[MD_RENDER_MOL_BOND_CLASS_COUNT + 1];
    md_gpu_addr_t bond_pairs;              /* uint32 [2 * bond_count]                  */
    md_gpu_addr_t bond_ids;                /* uint32 [bond_count]: index in sys->bond  */
    uint32_t      bond_cluster_count[MD_RENDER_MOL_BOND_CLASS_COUNT];
    md_gpu_addr_t bond_cluster_offset[MD_RENDER_MOL_BOND_CLASS_COUNT];  /* uint32 [count + 1],
                                                       relative to the class's first bond */

    /* The small molecules, regrouped by md_render_mol_update */
    uint32_t      small_molecule_count;
    uint32_t      small_atom_count;        /* the tail of atom_order                    */
    uint32_t      small_bond_count;        /* the tail of the covalent class            */
    uint32_t      grid_cells;              /* the last clusters of the atoms (and of the
                                              covalent bonds, when small_bond_count)    */
    md_gpu_addr_t small_tables;            /* internal: their atoms and bonds by molecule */
} md_render_mol_topo_t;

/* Host side of a topology to be written: the block's image in the caller's
   temp scope, and the byte size to allocate. Valid until that scope ends. */
typedef struct md_render_mol_plan_t {
    uint64_t             size;
    md_render_mol_topo_t layout;   /* addresses are offsets from 0                     */
    const void*          data;     /* [size] bytes                                     */
} md_render_mol_plan_t;

bool md_render_mol_plan(md_render_mol_plan_t* out, const struct md_system_t* sys, md_temp_scope_t temp);

/* Uploads the plan into [base, base + plan->size) and writes the topology
   value. Bumps content[MD_RENDER_CONTENT_GEOMETRY] when `content` is given. */
bool md_render_mol_write(md_gpu_stream_t stream, md_render_mol_topo_t* out, const md_render_mol_plan_t* plan,
                         md_gpu_addr_t base, uint64_t content[MD_RENDER_CONTENT_COUNT]);

/* =========================================================================
   Atoms: the caller's arrays
   ========================================================================= */

typedef struct md_render_mol_atoms_t {
    md_gpu_addr_t position;        /* float3, `stride` bytes apart                    */
    uint32_t      stride;          /* 0: 12                                           */
    uint32_t      velocity_stride; /* 0: 8                                            */
    md_gpu_addr_t velocity;        /* packed half3 per atom, optional                 */
    md_gpu_addr_t radius;          /* float per atom, optional: instead of the types' */
    md_gpu_addr_t radius_table;    /* float per particle type [topo.type_count]: for
                                      MD_RENDER_MOL_UPDATE_RADII                      */
    md_gpu_addr_t visible;         /* bitset over atoms: globally shown, optional     */
    md_gpu_addr_t selection;       /* bitset over atoms, optional                     */
    md_gpu_addr_t highlight;       /* bitset over atoms, optional                     */
} md_render_mol_atoms_t;

/* The coordinates of a state into `dst` (num_atoms * 12 bytes): md_system_state_t
   stores them as packed float3 already, which is what md_render reads, so this
   is one upload and the GEOMETRY counter bump that goes with it. */
bool md_render_mol_upload_state(md_gpu_stream_t stream, md_gpu_addr_t dst, const struct md_system_state_t* state,
                                uint64_t content[MD_RENDER_CONTENT_COUNT]);

typedef uint32_t md_render_mol_update_flags_t;
enum {
    /* The small molecules regrouped into grid cells by the current positions:
       the tails of topo's atom_order and covalent bonds, and their cluster
       offsets. Nothing to do without small molecules. */
    MD_RENDER_MOL_UPDATE_REGROUP = 1u << 0,
    /* topo's atom radii from atoms->radius_table, by particle type */
    MD_RENDER_MOL_UPDATE_RADII   = 1u << 1,
};

/* GPU work only, no readback, ordered in `stream` like any other. `dynamic`
   ([topo->dynamic_size] bytes, caller-owned) is scratch: its contents need
   not survive between calls. Bumps content[MD_RENDER_CONTENT_GEOMETRY] when
   it did anything. */
bool md_render_mol_update(md_render_mol_t mol, md_gpu_stream_t stream, const md_render_mol_topo_t* topo,
                          md_gpu_addr_t dynamic, const md_render_mol_atoms_t* atoms,
                          md_render_mol_update_flags_t flags, uint64_t content[MD_RENDER_CONTENT_COUNT]);

/* =========================================================================
   Representations
   ========================================================================= */

typedef enum md_render_mol_rep_type_t {
    MD_RENDER_MOL_REP_INVALID = 0,
    MD_RENDER_MOL_REP_SPACEFILL,
    MD_RENDER_MOL_REP_LICORICE,
    MD_RENDER_MOL_REP_BALL_AND_STICK,
    MD_RENDER_MOL_REP_BONDS,       /* sticks only, no atoms                          */
    MD_RENDER_MOL_REP_RIBBONS,     /* not yet implemented                            */
    MD_RENDER_MOL_REP_CARTOON,     /* not yet implemented                            */
} md_render_mol_rep_type_t;

typedef struct md_render_mol_rep_t {
    /* Stable identity, < 2^56: the batches it emits get keys (key << 8) | part */
    uint64_t                 key;
    md_render_mol_rep_type_t type;
    union {
        struct { float radius_scale; } spacefill;
        struct { float radius; } licorice;
        struct { float ball_scale; float stick_radius; } ball_and_stick;
        struct { float radius; } bonds;
    } params;

    /* Which bond classes the sticks draw: MD_RENDER_MOL_BONDS(...) mask.
       0 takes the type's default (covalent and coordinate; BONDS: all). */
    uint32_t                 bond_classes;

    /* Colour of atoms and sticks. COLOR_POINTS reads `color` (uint32 RGBA8
       per atom), COLOR_PALETTE reads `palette_index` (uint16 per atom). */
    md_render_color_t        color;
    md_gpu_addr_t            atom_color;
    md_gpu_addr_t            atom_palette_index;

    md_gpu_addr_t            visible;          /* the representation's filter, bitset over atoms */

    const float*             transform;        /* as md_render_batch_t: copied by md_render_view_draw */
    uint32_t                 transform_count;
} md_render_mol_rep_t;

/* What emit writes, as the low byte of batch keys */
enum {
    MD_RENDER_MOL_PART_ATOMS = 0,
    MD_RENDER_MOL_PART_BONDS = 1,   /* + bond class */
};

/* Appends the batches that draw `rep` to out[0 .. cap). Returns how many it
   needs, which may exceed cap (then nothing past cap is written). Host only,
   no GPU work. 0 for an invalid representation (logged). */
uint32_t md_render_mol_emit(const md_render_mol_topo_t* topo, const md_render_mol_atoms_t* atoms,
                            const md_render_mol_rep_t* rep, md_render_batch_t* out, uint32_t cap);

/* =========================================================================
   Picking, in molecular terms
   ========================================================================= */

typedef enum md_render_mol_hit_kind_t {
    MD_RENDER_MOL_HIT_NONE = 0,
    MD_RENDER_MOL_HIT_ATOM,
    MD_RENDER_MOL_HIT_BOND,
} md_render_mol_hit_kind_t;

typedef struct md_render_mol_hit_t {
    md_render_mol_hit_kind_t kind;
    uint64_t rep_key;              /* the representation's key                        */
    uint32_t instance;
    uint32_t atom;                 /* the atom, or the bond's atom nearer the hit     */
    uint32_t bond;                 /* index in sys->bond; UINT32_MAX for an atom      */
    float    position[3];
} md_render_mol_hit_t;

/* A pick from a frame whose batches came from md_render_mol_emit. False for a
   pick of anything else (keys whose part it did not emit). */
bool md_render_mol_decode_pick(const md_render_pick_t* pick, md_render_mol_hit_t* out_hit);

#ifdef __cplusplus
}
#endif

#endif /* MD_RENDER_MOL_H */
