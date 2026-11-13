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
        n = md_render_mol_emit(&topo, dyn, &atoms, &rep, batches, cap);

    When to regroup is the caller's call: on a new frame of a trajectory,
    not for the steps interpolated between two frames (the grouping of the
    first is still a good one for those). Skipping it is never wrong, only
    slower: the bounds follow the positions regardless, the clusters just
    grow as the molecules spread. Before the first regroup the small
    molecules are grouped in molecule order.

  * What md_render draws is a copy of the atoms in cluster order, in the
    per-frame block: float4 (position, radius) and the velocity per point,
    gathered by every md_render_mol_update. Read through the cluster order
    instead, the caller's arrays in file order cost a scattered fetch in
    every pass over the elements; the copy costs one per position change.
    Everything per atom that is not geometry (colours, palette indices,
    visibility, selection) is read in place, through the point's atom
    (attr_index = atom_order), and picks report atoms. The update must run
    after every change of positions, velocities or radii.

  * The positions can be the caller's, or frames of a trajectory resident on
    the GPU (see Frames): each uploaded once, quantized, and interpolated by
    the update with the velocities of the curve, so nothing is uploaded or
    interpolated on the host per step.

  * Only the resident atoms have points: what the representations may show,
    given to md_render_mol_plan. The rest (water no representation shows)
    costs nothing past the plan: no point, no cluster, no bond, no frame
    data, no regroup, no copy. Per-atom arrays (colours, visibility,
    selection, the caller's positions) stay indexed by atom.

  * Radii come from the system's particle types: the topology holds each
    atom's type and a radius per type (sys->atom.type.radius when planned).
    A changed radius per type is a table and an update
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
#include <core/md_bitset.h>

struct md_system_t;
struct md_system_state_t;
struct md_unitcell_t;

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
    uint64_t      dynamic_size;    /* bytes of the per-frame block: the copy that is
                                      drawn, and the update's scratch                */
    uint64_t      dynamic_size_velocity; /* ... when md_render_mol_atoms_t.velocity is
                                      given: the velocities are copied too           */
    uint32_t      atom_count;      /* the system's: per-atom arrays are this long     */
    uint32_t      point_count;     /* the resident atoms: what has a point, what
                                      frames hold (atom_count when all are)          */
    md_gpu_addr_t resident;        /* internal: the resident atoms as bits and ranks;
                                      0 when all are                                 */
    md_gpu_addr_t regroups;        /* internal: a word the regroups count in (a copy
                                      gathered before the last one is in the old order) */

    /* Particle types and radii, per atom */
    uint32_t      type_count;              /* entries a radius table has: sys->atom.type.count */
    md_gpu_addr_t atom_type;               /* uint16 [atom_count]                      */
    md_gpu_addr_t type_radius;             /* float [type_count]: each type's radius   */

    /* The resident atoms in cluster order: residue clusters first, then one
       cluster per grid cell (grid_cells of them, some empty) for the small
       molecules. Point p (of the per-frame copy) is atom atom_order[p]. */
    md_gpu_addr_t atom_order;              /* uint32 [point_count]: point -> atom       */
    md_gpu_addr_t atom_point;              /* uint32 [atom_count]: atom -> point,
                                              UINT32_MAX for an atom not resident       */
    uint32_t      atom_cluster_count;
    md_gpu_addr_t atom_cluster_offset;     /* uint32 [atom_cluster_count + 1]          */

    /* Bonds, grouped by class: class c is bond_pairs[class_offset[c] ..
       class_offset[c + 1]), in cluster order: with the cluster of their
       first atom's (64 at most per cluster); the covalent bonds within small
       molecules last, by grid cell. Between points, not atoms (atom_order
       maps back); the order within a class changes with every regroup. */
    uint32_t      bond_count;
    uint32_t      bond_class_offset[MD_RENDER_MOL_BOND_CLASS_COUNT + 1];
    md_gpu_addr_t bond_pairs;              /* uint32 [2 * bond_count]: two points      */
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
    uint32_t      mixed_bond_count;        /* other bonds touching small molecules:
                                              their points follow every regroup         */
    md_gpu_addr_t small_tables;            /* internal: their atoms and bonds by molecule */
} md_render_mol_topo_t;

/* Host side of a topology to be written: the block's image in the caller's
   temp scope, and the byte size to allocate. Valid until that scope ends. */
typedef struct md_render_mol_plan_t {
    uint64_t             size;
    md_render_mol_topo_t layout;   /* addresses are offsets from 0                     */
    const void*          data;     /* [size] bytes                                     */
} md_render_mol_plan_t;

/* The resident atoms: the members of `resident` (NULL: every atom). The
   rest have no point: nothing of them is read, uploaded, interpolated or
   drawn, and bonds with an end among them are left out. Make it the union of
   what the representations may show (a selection that changes with the
   frame contributes its scope, not its current members): water no
   representation shows costs nothing. A small molecule keeps the members
   among its atoms. A change of the set is a new plan. */
bool md_render_mol_plan(md_render_mol_plan_t* out, const struct md_system_t* sys, const md_bitset_t* resident, md_temp_scope_t temp);

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
    md_gpu_addr_t velocity;        /* packed half3 per atom, optional (the same for
                                      update and emit: copied, or not); needs a
                                      per-frame block of topo.dynamic_size_velocity  */
    md_gpu_addr_t radius;          /* float per atom, optional: instead of the types' */
    md_gpu_addr_t radius_table;    /* float per particle type [topo.type_count]: for
                                      MD_RENDER_MOL_UPDATE_RADII                      */
    md_gpu_addr_t visible;         /* bitset over atoms: globally shown, optional     */
    md_gpu_addr_t selection;       /* bitset over atoms, optional                     */
    md_gpu_addr_t highlight;       /* bitset over atoms, optional                     */

    /* Instead of position and velocity: frames resident on the GPU (see
       Frames below), combined by the update. In a periodic cell frames [0]
       and [2] of four are taken at the minimum image of [1], and [3] at that
       of [2]; of two, [1] at that of [0] (the cell of [1] of four, of [0] of
       two). md_render_mol_frame_weights writes the weights. */
    const struct md_render_mol_frames_t* frames;
    uint32_t      frame_count;     /* 1: frame_slot[0] as it is; 2 or 4: weighted     */
    uint32_t      frame_slot[4];
    float         frame_w[4];      /* position = sum of w[k] frame[k]                 */
    float         frame_dw[4];     /* velocity = sum of dw[k] frame[k]. All 0: at rest,
                                      no velocities (as without `velocity`); else
                                      needs a block of topo.dynamic_size_velocity    */
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
    /* The radii changed: topo's radius per type copied from
       atoms->radius_table when given, and the copy's radii from those (or
       from atoms->radius). Also after changing atoms->radius (its contents
       or which array). */
    MD_RENDER_MOL_UPDATE_RADII   = 1u << 1,
};

/* A regroup reorders the topology: every per-frame block gathered from it
   before is stale (its copy in the old order) until it is updated again.
   One per-frame block per topology is the usual case. With atoms->frames
   the regroup belongs to md_render_mol_frame_upload, and REGROUP here
   fails.

   After positions, velocities or radii changed (and before the first draw):
   the copy that is drawn gathered into `dynamic` ([topo->dynamic_size]
   bytes, caller-owned, kept for the batches emitted from it), after what
   the flags ask. GPU work only, no readback, ordered in `stream` like any
   other. Bumps content[MD_RENDER_CONTENT_GEOMETRY]. */
bool md_render_mol_update(md_render_mol_t mol, md_gpu_stream_t stream, const md_render_mol_topo_t* topo,
                          md_gpu_addr_t dynamic, const md_render_mol_atoms_t* atoms,
                          md_render_mol_update_flags_t flags, uint64_t content[MD_RENDER_CONTENT_COUNT]);

/* =========================================================================
   Frames: a trajectory resident on the GPU, interpolated there
   ========================================================================= */

/* Instead of the caller's positions, the update can combine frames held in
   GPU memory: up to MD_RENDER_MOL_FRAME_SLOTS of them in one caller-owned
   block, each uploaded once and interpolated on the GPU for every step
   between them, with the velocities of the curve. What crosses the bus is
   one small frame per frame of the trajectory, nothing per step.

   A frame is quantized on the host: each axis to 16 bits over the bounds of
   the resident atoms, p = origin + q * scale, x, y and z per resident atom
   in increasing atom order, 6 bytes per atom. The step is the extent over
   65535: 0.01 A up to a 655 A system (the precision of an xtc at its
   default, so such a frame loses nothing the file had), 0.03 A at 2000 A;
   resident atoms that span less get a finer step. Atoms not resident are
   neither quantized nor uploaded.

   An upload reads the frame from the host (copied through the stream) or,
   better, from GPU memory: quantized by the host straight into
   MD_GPU_MEM_HOST_WRITE memory (on unified memory and with resizable BAR
   that is the transfer: nothing is copied), or written there by a GPU
   decoder. The upload reorders it into its slot; one gather of 6 bytes per
   atom.

   Resident frames are in the point order, so the update reads them in
   order: one coalesced pass, 6 bytes per frame and atom in, the copy out.
   A regroup reorders the points, so it happens at an upload, by the frame
   uploaded, and the frames already resident are reordered along (one pass
   over each). The usual pattern is a regroup with every new frame of the
   trajectory:

       md_render_mol_frames_t fr;
       md_gpu_addr_t fb = md_gpu_malloc(s, MD_GPU_MEM_DEVICE, md_render_mol_frames_size(&topo, 4)).gpu;
       md_render_mol_frames_init(&fr, &topo, fb, 4);
       md_gpu_mem_t up = md_gpu_malloc(s, MD_GPU_MEM_HOST_WRITE, md_render_mol_frame_size(&topo));
       ...
       md_render_mol_frame_quantize(&f, up.cpu, xyz, 0, n, cell, resident);   // host, any thread
       f.q_gpu = up.gpu;
       md_render_mol_frame_upload(mol, s, &topo, dyn, &fr, slot, &f, MD_RENDER_MOL_UPDATE_REGROUP);
       ...
       atoms.frames = &fr;                                            // per step
       atoms.frame_count = 4;  atoms.frame_slot = { ... };            // the slots of frames i-1 .. i+2
       md_render_mol_frame_weights(atoms.frame_w, atoms.frame_dw, 4, t, dt, 0.5f);
       md_render_mol_update(mol, s, &topo, dyn, &atoms, 0, content);

   Memory written for an upload must not be written again before the stream
   has passed the upload (two such buffers in turn, with a sync point each,
   let the host quantize the next frame meanwhile). Which slot holds which
   frame is the caller's business (a ring, a cache).

   The block: the frames, and for a regroup the old place of every small
   atom: slots * 6 bytes per resident atom, and 4 per resident atom of the
   small molecules. */

#define MD_RENDER_MOL_FRAME_SLOTS 4

/* A frame, as uploaded */
typedef struct md_render_mol_frame_t {
    float           origin[3];
    float           scale[3];      /* p = origin + q * scale, per axis               */
    float           cell[9];       /* columns a, b, c; a zero column: not periodic
                                      along it (all zero: not periodic)              */
    const uint16_t* q;             /* x, y, z per resident atom, [3 * point_count],
                                      on the host                                    */
    md_gpu_addr_t   q_gpu;         /* or (when nonzero) in GPU memory: 4-byte aligned,
                                      md_render_mol_frame_size bytes                 */
} md_render_mol_frame_t;

/* The positions of the resident atoms among `count` (float3, `stride`
   bytes apart; 0: 12) into q (x, y, z per resident atom, in increasing atom
   order: [3 * topo.point_count]; written in order and never read, so it may
   be write-combined memory) over their bounds, and the frame of them (out->q
   = q). `resident`: the plan's set (NULL: all `count`). NaN goes to the
   origin. `cell` may be NULL (not periodic). */
void md_render_mol_frame_quantize(md_render_mol_frame_t* out, uint16_t* q, const float* xyz, uint32_t stride,
                                  uint32_t count, const float cell[9], const md_bitset_t* resident);

/* An md_unitcell_t as such columns: the periodic axes' (PBC flags) */
void md_render_mol_frame_cell(float out[9], const struct md_unitcell_t* cell);

/* Bytes of a frame in GPU memory (q_gpu): 6 per resident atom, in whole words */
uint64_t md_render_mol_frame_size(const md_render_mol_topo_t* topo);

/* The block's frames, in the point order of one topology. A plain value
   that uploads change: keep one per topology. */
typedef struct md_render_mol_frames_t {
    md_gpu_addr_t base;            /* the block; md_gpu_free(stream, base) frees it  */
    uint64_t      size;
    md_gpu_addr_t topology_base;   /* whose point order: topo.base, topo.topology_version */
    uint64_t      topology_version;
    uint32_t      slots;
    uint32_t      point_count;
    uint32_t      valid;           /* bit k: slot k holds a frame                    */
    uint32_t      buffer[MD_RENDER_MOL_FRAME_SLOTS];   /* internal: each slot's buffer */
} md_render_mol_frames_t;

/* Bytes of a block of `slots` frames (1 .. MD_RENDER_MOL_FRAME_SLOTS) */
uint64_t md_render_mol_frames_size(const md_render_mol_topo_t* topo, uint32_t slots);

/* An empty block at `base` (md_render_mol_frames_size bytes, caller-owned) */
bool md_render_mol_frames_init(md_render_mol_frames_t* out, const md_render_mol_topo_t* topo, md_gpu_addr_t base, uint32_t slots);

/* `frame` into slot `slot`, replacing what it held. With
   MD_RENDER_MOL_UPDATE_REGROUP (other flags are ignored): the small
   molecules regrouped by this frame's positions, which needs `dynamic` (the
   regroup's scratch) and reorders the topology and the other resident
   frames along; then update before drawing. GPU work only, ordered in
   `stream`. A frame on the host is copied before this returns; one in GPU
   memory is read when the stream gets to it. */
bool md_render_mol_frame_upload(md_render_mol_t mol, md_gpu_stream_t stream, const md_render_mol_topo_t* topo,
                                md_gpu_addr_t dynamic, md_render_mol_frames_t* frames, uint32_t slot,
                                const md_render_mol_frame_t* frame, md_render_mol_update_flags_t flags);

/* Weights of `frame_count` frames (1, 2 or 4) at t in [0, 1] between the two
   interpolated ([0] and [1] of two, [1] and [2] of four), and of their
   derivative over `dt`, the time from one frame to the next in md_render's
   frame.time units (0: no velocities). Four: a cardinal spline of tangent
   scale `tension` (0.5: Catmull-Rom), as md_util_cubic_spline_weights. */
void md_render_mol_frame_weights(float w[4], float dw[4], uint32_t frame_count, float t, float dt, float tension);

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
uint32_t md_render_mol_emit(const md_render_mol_topo_t* topo, md_gpu_addr_t dynamic, const md_render_mol_atoms_t* atoms,
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
