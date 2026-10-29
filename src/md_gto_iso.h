#pragma once

// md_gto_iso - isosurface meshes of molecular orbitals and densities, built straight from the basis
//
// No volume is ever evaluated. The surface f = +iso / f = -iso of an orbital (f = psi) or a
// density (f = rho: SCF, spin, transition, attachment / detachment, any AO density matrix) is found
// with a lazily built octree whose cells are kept only where a rigorous bound on f admits the iso
// value, stopped at bricks of 4^3 lattice cells. The bricks' lattice corners are evaluated once;
// Surface Nets on them gives the topology, Newton steps along the analytic gradient put the
// vertices ON the true surface, tangential relaxation evens out the triangles and optional edge
// splitting refines where the surface curves.
//
// THE GUARANTEE
// A lattice cell whose corners show a crossing is never discarded: the bound is conservative, so
// the mesh is exactly the Surface Nets mesh of the full lattice at the same spacing, at a fraction
// of the evaluations (typically 10-100x fewer). What no method at spacing h sees - a surface piece
// that crosses no lattice edge - is not seen here either.
//
// TWO TIMESCALES
//   md_gto_iso_cache_t   per field and spacing, for an interval of iso values: the octree, the
//                        bricks and their corner values (~0.5 KB per brick). Building it is the
//                        expensive part; for a density mostly the corner values, which can be
//                        spread over threads (md_gto_iso_cache_desc_t.defer_corners).
//   md_gto_iso_mesh_t    per iso value inside the cache's interval: Surface Nets on the cached
//                        corners plus the vertex refinement. Dragging an iso slider inside the
//                        interval only pays for this.
//
// The field and the cache are immutable once created, so meshes can be extracted from several
// threads at once, and different fields can be processed in parallel.
//
// PRECISION AND UNITS
// Double precision throughout, positions in bohr, the same atom coordinates as the basis is
// evaluated at. Corner values are stored as float, which only decides on which side of the iso
// value a corner lies; vertex positions come from the double precision projection.

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#include <md_gto.h>

struct md_allocator_i;

#ifdef __cplusplus
extern "C" {
#endif

// ---------------------------------------------------------------------------
// FIELD: an orbital or a density, prepared for bounding and evaluation
// ---------------------------------------------------------------------------
// Everything below (cache, mesh) works the same for both kinds; f denotes the field's value. The
// positive sheet is f = +iso, the negative one f = -iso (an orbital's two phases; the negative
// regions of a difference or transition density).
//
// Common to the constructors: atom_xyz in bohr, stride in bytes (0 = tightly packed float3); AO
// vectors and matrices in the Cartesian AO order of the basis (md_gto_basis_num_ao entries per
// vector), as for md_gto_expand_with_ao_coeffs and the orbital/coefficient attributes. Each returns
// NULL on failure.
typedef struct md_gto_iso_field_t md_gto_iso_field_t;

// An orbital, f = psi = sum_i c_i phi_i.
md_gto_iso_field_t* md_gto_iso_field_create(const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride,
                                            const double* ao_coeffs, struct md_allocator_i* alloc);

// A density as a sum of orbital products, f = sum_k w_k a_k b_k, with a_k = left[k] . phi and
// b_k = right[k] . phi (count rows of num_ao each). right == NULL means squares, b_k = a_k.
// This is the form a density usually arrives in, so nothing is decomposed:
//   SCF density:            left = the occupied MOs, weight = the occupations
//   spin density:           alpha and beta MOs, weights +n and -n
//   transition density:     left = NTO holes, right = NTO particles, weight = the NTO amplitudes
//   attachment/detachment:  NTO particles / holes squared, weight = the NTO weights (lambda)
// Cost grows with count (each function is bounded and evaluated on its own), so pass only the
// significant products.
md_gto_iso_field_t* md_gto_iso_field_create_products(const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride,
                                                     const double* left, const double* right, const double* weight, size_t count,
                                                     struct md_allocator_i* alloc);

// A density from its AO density matrix, f = phi^T D phi (num_ao x num_ao, row major). Only the
// symmetric part of D contributes, so a non-symmetric (transition) matrix is fine. D is split into
// signed squares: pivoted LDL^T (pivoted Cholesky for a semidefinite D, whose factors are localized
// and give tight bounds), or an eigendecomposition where 1x1 pivots are unstable. The splitting's
// residual is bounded everywhere and carried in every bound (md_gto_iso_field_info: residual).
// The number of factors is the rank of D: the occupied count for an SCF density. A correlated
// density is close to full rank and costs accordingly; prefer _create_products where the products
// are at hand.
md_gto_iso_field_t* md_gto_iso_field_create_density(const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride,
                                                    const double* density_matrix, struct md_allocator_i* alloc);

void md_gto_iso_field_destroy(md_gto_iso_field_t* field);

// A hint: every local maximum of the field is at an atom, as for an SCF electron density (total, or of
// one spin). The GPU extraction's flood (md_gto_iso_gpu_switches_t.flood = 1, the default) relies on it:
// it evaluates a density's lattice only around the surface, finding each piece from a seed, and a piece
// with an atom in it always has one. Attachment, detachment and other densities whose maxima lie
// between the atoms leave it off (their small pieces could be missed), and are extracted exactly.
void md_gto_iso_field_set_maxima_at_atoms(md_gto_iso_field_t* field, bool maxima_at_atoms);

enum {
    MD_GTO_ISO_FIELD_ORBITAL = 0,
    MD_GTO_ISO_FIELD_DENSITY = 1,
};

enum {
    MD_GTO_ISO_SPLIT_NONE = 0,      // an orbital, or products given by the caller
    MD_GTO_ISO_SPLIT_LDL  = 1,      // a density matrix, pivoted LDL^T
    MD_GTO_ISO_SPLIT_EIGEN = 2,     // a density matrix, eigendecomposition
};

typedef struct md_gto_iso_field_info_t {
    uint32_t kind;                  // MD_GTO_ISO_FIELD_*
    uint32_t split;                 // MD_GTO_ISO_SPLIT_*
    uint32_t num_functions;         // orbital-like functions (1 for an orbital, the rank for a matrix)
    uint32_t num_terms;             // products (0 for an orbital)
    uint32_t num_groups;            // what a bound sums over: (atom, l) groups of an orbital, shells of a density
    double   residual;              // bound of |phi^T D phi - the split| over all space (matrix input)
    size_t   memory_bytes;
} md_gto_iso_field_info_t;

void md_gto_iso_field_info(md_gto_iso_field_info_t* out, const md_gto_iso_field_t* field);

// f (and optionally its gradient, 3 per point) at points in bohr.
void md_gto_iso_field_evaluate(const md_gto_iso_field_t* field, const float* xyz, size_t count, size_t stride,
                               double* out_psi, double* out_grad);

// ---------------------------------------------------------------------------
// CACHE: bricks and corner values for every iso value in [iso_min, iso_max]
// ---------------------------------------------------------------------------
typedef struct md_gto_iso_cache_t md_gto_iso_cache_t;

typedef struct md_gto_iso_cache_desc_t {
    double cell_size;   // lattice spacing in bohr; 0.1 - 0.15 is a good range for display
    double iso_min;     // > 0. The cache serves every iso in [iso_min, iso_max], both signs.
    double iso_max;     // >= iso_min; 0 means iso_min (a single value)
    bool   defer_corners;   // return after the octree; the caller fills the corner values with
                            // md_gto_iso_cache_evaluate, e.g. spread over threads
} md_gto_iso_cache_desc_t;

// The field must outlive the cache. Returns NULL on failure.
md_gto_iso_cache_t* md_gto_iso_cache_create(const md_gto_iso_field_t* field, const md_gto_iso_cache_desc_t* desc, struct md_allocator_i* alloc);
void md_gto_iso_cache_destroy(md_gto_iso_cache_t* cache);

// Evaluates the corner values of bricks [first_brick, first_brick + num_bricks) (num_bricks from
// md_gto_iso_cache_info). Building the cache does all of them unless desc.defer_corners was set;
// then every brick must be evaluated before a mesh is extracted. Disjoint ranges may be evaluated
// from several threads at once. This is the expensive part of a density cache: one AO evaluation
// and a num_ao x num_functions product per corner, 125 corners per brick.
void md_gto_iso_cache_evaluate(md_gto_iso_cache_t* cache, uint32_t first_brick, uint32_t num_bricks);

// Whether a mesh at this iso value can be extracted from the cache.
bool md_gto_iso_cache_covers(const md_gto_iso_cache_t* cache, double iso);

typedef struct md_gto_iso_cache_info_t {
    double   cell_size;
    double   iso_min, iso_max;
    double   origin[3];         // lattice origin, bohr
    uint32_t num_bricks;        // 4^3 cells each
    uint32_t cells_per_axis;    // of the (virtual) full lattice
    size_t   memory_bytes;
    double   seconds_octree;    // building it: the bounds
    double   seconds_corners;   // building it: the corner values of the bricks
} md_gto_iso_cache_info_t;

void md_gto_iso_cache_info(md_gto_iso_cache_info_t* out, const md_gto_iso_cache_t* cache);

// ---------------------------------------------------------------------------
// MESH
// ---------------------------------------------------------------------------
enum {
    MD_GTO_ISO_SHEET_POSITIVE = 1,   // f = +iso
    MD_GTO_ISO_SHEET_NEGATIVE = 2,   // f = -iso
    MD_GTO_ISO_SHEET_BOTH     = 3,
};

typedef struct md_gto_iso_mesh_desc_t {
    double   iso;                  // > 0 and covered by the cache
    uint32_t sheets;               // MD_GTO_ISO_SHEET_*, 0 = both
    // Zero initialized: the raw Surface Nets mesh with analytic normals - the cheap choice while an
    // iso value is being dragged. For display quality: project 8, relax 3.
    uint32_t project_iterations;   // Newton steps onto the true surface (stops early when converged)
    uint32_t relax_iterations;     // tangential relaxation passes, each followed by a projection
    uint32_t refine_passes;        // edge splits where the mesh strays from the surface (0-2)
} md_gto_iso_mesh_desc_t;

// Triangles are counter-clockwise seen from outside, outside being away from |f| > iso; normals
// point outward and are the analytic surface normals at the vertices. The positive sheet comes
// first, in both the vertex and the index arrays.
typedef struct md_gto_iso_mesh_t {
    uint32_t  vertex_count;
    uint32_t  index_count;
    float*    position;              // 3 per vertex, bohr
    float*    normal;                // 3 per vertex
    uint32_t* index;                 // 3 per triangle
    uint32_t  sheet_vertex_offset[2], sheet_vertex_count[2];   // [0] positive, [1] negative
    uint32_t  sheet_index_offset[2],  sheet_index_count[2];
    struct md_allocator_i* alloc;
} md_gto_iso_mesh_t;

bool md_gto_iso_mesh_extract(md_gto_iso_mesh_t* out, const md_gto_iso_cache_t* cache, const md_gto_iso_mesh_desc_t* desc, struct md_allocator_i* alloc);
void md_gto_iso_mesh_free(md_gto_iso_mesh_t* mesh);

#if MD_ENABLE_GPU
// ---------------------------------------------------------------------------
// GPU: the same meshes, built on md_gpu
// ---------------------------------------------------------------------------
// The CPU path above is the double precision reference; this one is the fast path. The pipeline is
// the same - octree from rigorous bounds, bricks of 4^3 cells, Surface Nets, Newton projection,
// relaxation - in float, with every rounding step of the bounds accounted for, so no lattice cell
// whose (float) corner values cross the iso value is missed.
//
// A SET is up to MD_GTO_ISO_GPU_MAX_FIELDS fields over the same basis and atoms - the orbitals of a
// grid of views, or one density - extracted together: the octree and the corner values are shared,
// so the AOs at a corner are evaluated once whatever the number of fields. Orbitals and densities
// mix freely.

#define MD_GTO_ISO_GPU_MAX_FIELDS 16

// Kernels are created with the first set (md_gto_iso_gpu_set_create, on its thread), so this costs nothing until
// then. A set may be made on any thread, on a stream of that thread's (md_gpu's streams are each one thread's):
// the GPU work of a set made on another stream is waited for (md_gpu_stream_wait) before it is extracted.
void md_gto_iso_gpu_initialize(md_gpu_device_t device);
void md_gto_iso_gpu_shutdown(void);

typedef struct md_gto_iso_gpu_set_t md_gto_iso_gpu_set_t;

// A field may be given more than once (one density at several iso values): its copies share their
// functions, so they cost nothing beyond the extra sheets.
// The fields must have been created from this basis and these atom positions, and must outlive the
// set (each extraction derives the root box from their tails). Uploads on 'stream'. Returns NULL on
// failure (too many fields, more than 1024 functions in all, and not all densities).
// A set of densities keeps each as its matrix too (phi^T D phi: the octree's bound, and with more than 1024
// functions in all the corners and vertices). With md_gto_gpu_initialize called on the device, the matrices
// are in tiles of the AO pairs that meet (md_gto_gpu_dmat), built on the device from the densities'
// functions, and the call waits for that; otherwise dense, num_ao^2 floats each, made on the host.
md_gto_iso_gpu_set_t* md_gto_iso_gpu_set_create(md_gpu_stream_t stream, const md_gto_basis_t* basis, const float* atom_xyz, size_t atom_xyz_stride,
                                                const md_gto_iso_field_t* const* fields, uint32_t num_fields, struct md_allocator_i* alloc);
void md_gto_iso_gpu_set_destroy(md_gpu_stream_t stream, md_gto_iso_gpu_set_t* set);

typedef struct md_gto_iso_gpu_desc_t {
    double   cell_size;                               // lattice spacing in bohr
    double   iso[MD_GTO_ISO_GPU_MAX_FIELDS];          // per field, > 0
    uint32_t sheets[MD_GTO_ISO_GPU_MAX_FIELDS];       // per field, MD_GTO_ISO_SHEET_*, 0 = both
    uint32_t project_iterations;                      // as md_gto_iso_mesh_desc_t
    uint32_t relax_iterations;
    bool     time_stages;                             // wait for and time every stage (slower)
    bool     read_back;                               // jobs: also copy the mesh to the host (mesh.host_*)
} md_gto_iso_gpu_desc_t;

enum {
    MD_GTO_ISO_GPU_STAGE_OCTREE,
    MD_GTO_ISO_GPU_STAGE_CORNERS,
    MD_GTO_ISO_GPU_STAGE_NETS,          // cells and quads
    MD_GTO_ISO_GPU_STAGE_PROJECT,
    MD_GTO_ISO_GPU_STAGE_RELAX,
    MD_GTO_ISO_GPU_STAGE_OUTPUT,        // regrouping, triangles, the final readback of the counts
    MD_GTO_ISO_GPU_STAGE_COUNT
};

// The meshes of a set, in device memory: one vertex and one index buffer for all fields, each field
// and sheet a contiguous range ([field][0] positive, [field][1] negative). Indices are into the
// whole vertex buffer. Triangles are counter-clockwise from outside, as on the CPU.
typedef struct md_gto_iso_gpu_mesh_t {
    md_gpu_addr_t position;            // float4 per vertex: xyz in bohr, w = 1
    md_gpu_addr_t normal;              // float4 per vertex: outward unit normal, w: |f| / |grad f|
    md_gpu_addr_t index;               // uint32 per triangle corner
    uint32_t vertex_count;
    uint32_t index_count;
    uint32_t vertex_offset[MD_GTO_ISO_GPU_MAX_FIELDS][2], vertex_count_of[MD_GTO_ISO_GPU_MAX_FIELDS][2];
    uint32_t index_offset[MD_GTO_ISO_GPU_MAX_FIELDS][2],  index_count_of[MD_GTO_ISO_GPU_MAX_FIELDS][2];
    uint32_t num_fields;
    double   origin[3];                // the lattice: origin (bohr), cell size, cells per axis of the root
    double   cell_size;
    uint32_t cells_per_axis;
    // what it took
    uint32_t num_bricks;               // union over the fields
    uint32_t num_corner_slots;         // (brick, field) pairs whose corners were stored
    uint32_t open_edges;               // crossing edges whose four cells were not all present (0)
    double   seconds[MD_GTO_ISO_GPU_STAGE_COUNT];   // with desc.time_stages
    double   seconds_total;            // wall clock from the start to the mesh
    // With desc.read_back (jobs): the mesh on the host, as on the device (position and normal float4 per
    // vertex, index uint32 per triangle corner); valid until md_gto_iso_gpu_mesh_free
    const float*    host_position;
    const float*    host_normal;
    const uint32_t* host_index;
    md_gpu_addr_t   host_mem;
} md_gto_iso_gpu_mesh_t;

// Blocking: a job (below) polled until it is done. Returns false on failure (lattice too large for
// the cell size, out of device memory).
bool md_gto_iso_gpu_extract(md_gpu_stream_t stream, const md_gto_iso_gpu_set_t* set, const md_gto_iso_gpu_desc_t* desc, md_gto_iso_gpu_mesh_t* out);
void md_gto_iso_gpu_mesh_free(md_gpu_stream_t stream, md_gto_iso_gpu_mesh_t* mesh);

// The extraction as a job that never blocks the caller: begin issues the GPU work up to the first
// point where the host needs a count from the device (the bricks, the vertices, the triangles: four
// or more such points, one per flood round besides); poll takes each count once the GPU has made it
// and issues the work after it. A job polled once a frame runs alongside the rendering; what the GPU
// does meanwhile is as without it. Its device memory is its own (md_gpu_malloc), so the stream may be
// used for other work between polls, temp scopes included.
//   md_gto_iso_gpu_job_t* job = md_gto_iso_gpu_job_begin(stream, set, &desc);
//   ... once a frame: if (md_gto_iso_gpu_job_poll(job, 2.0) != MD_GTO_ISO_GPU_JOB_RUNNING) {
//           md_gto_iso_gpu_job_end(job, &mesh); ... }
// The set must outlive the job. One thread at a time, the stream's (md_gpu.h).
typedef struct md_gto_iso_gpu_job_t md_gto_iso_gpu_job_t;

typedef enum md_gto_iso_gpu_job_status_t {
    MD_GTO_ISO_GPU_JOB_RUNNING,
    MD_GTO_ISO_GPU_JOB_DONE,
    MD_GTO_ISO_GPU_JOB_FAILED,
} md_gto_iso_gpu_job_status_t;

// NULL on invalid arguments or when the kernels are missing; the job may still fail later
md_gto_iso_gpu_job_t* md_gto_iso_gpu_job_begin(md_gpu_stream_t stream, const md_gto_iso_gpu_set_t* set, const md_gto_iso_gpu_desc_t* desc);

// Advances the job as far as the GPU allows, waiting for it at most budget_ms in all (0: not at all,
// < 0: as long as it takes, which blocks until the job is done)
md_gto_iso_gpu_job_status_t md_gto_iso_gpu_job_poll(md_gto_iso_gpu_job_t* job, double budget_ms);

// How far the job is, 0 to 1, as of the last poll: by the counts the host has taken (a few steps), and in
// a long job (its octree took 4 ms or more) by the chunks of its corner and vertex launches that have
// passed: those go in chunks of about 10-20 ms of work, each followed by a marker
float md_gto_iso_gpu_job_progress(const md_gto_iso_gpu_job_t* job);

// Ends the job. Done: the mesh goes to 'out' (free it with md_gto_iso_gpu_mesh_free) and true is
// returned. Running or failed: what it made is dropped (the GPU work in flight completes unseen), 'out'
// is zeroed and false is returned. Never waits for the GPU.
bool md_gto_iso_gpu_job_end(md_gto_iso_gpu_job_t* job, md_gto_iso_gpu_mesh_t* out);

// One field's mesh copied to the host in the CPU layout (blocking). Free with md_gto_iso_mesh_free.
bool md_gto_iso_gpu_mesh_read(md_gpu_stream_t stream, const md_gto_iso_gpu_mesh_t* mesh, uint32_t field, md_gto_iso_mesh_t* out, struct md_allocator_i* alloc);

// All of it as it is on the device, with one wait: position and normal float4 per vertex (vertex_count),
// index uint32 per triangle corner (index_count) into the whole vertex buffer. Any may be NULL.
bool md_gto_iso_gpu_mesh_read_all(md_gpu_stream_t stream, const md_gto_iso_gpu_mesh_t* mesh, float* position, float* normal, uint32_t* index);

// ---------------------------------------------------------------------------
// GPU: the fields of a set at arbitrary points
// ---------------------------------------------------------------------------
// The value and gradient of any field of a set wherever they are wanted - the vertices of a mesh (of
// another field, or of another surface altogether), probe points, a path - by the kernels that evaluate
// the extraction's vertices while projecting them onto their surface: the same evaluation, in float,
// whatever the field (orbitals, products of orbitals, densities as phi^T D phi). The set needs no
// extraction first and no iso value.
//
// What a field at a point costs depends on the AOs that reach it, found per group of consecutive points
// (16 to 32 of them) against the group's bounding box: points in a spatially coherent order - a mesh's
// vertices are - cost a fraction of the same points shuffled. The order never changes the values.
//
// Precision: float, as the extraction's vertices, with the extraction's evaluation cutoff (GTO_ISO_EVAL_EPS,
// 1e-6 by default): an AO is left out where it is below it, so the error has an absolute part of that
// order besides float rounding. Against md_gto_iso_field_evaluate in double, at 4000 random points around
// amide and a 26 atom molecule (DEF2-SVP): orbitals within 3.2e-6 absolute (|psi| up to 0.5), SCF
// densities within 4.1e-7 of their largest value, gradients alike.
//
// A value that cannot be made is NaN and counted (out_failed): a density in the D form (more than 1024
// functions in the set) at a point that too many shells reach for one group's AO list. Nothing else fails
// once the launch has been issued.
typedef struct md_gto_iso_gpu_eval_desc_t {
    md_gpu_addr_t xyz;              // positions in bohr, float: point i at xyz + i * xyz_stride
    uint32_t      xyz_stride;       // bytes from one point to the next, a multiple of 4; 0 = packed float3 (12).
                                    // A mesh's positions (md_gto_iso_gpu_mesh_t.position): 16
    uint32_t      num_xyz;
    uint32_t      num_fields;       // 1 .. MD_GTO_ISO_GPU_MAX_FIELDS
    uint32_t      fields[MD_GTO_ISO_GPU_MAX_FIELDS];    // the set's fields to evaluate, by index; repeats allowed
    // At least one of the two. Field k (the k:th of fields[]) of point i at [k * num_xyz + i]:
    md_gpu_addr_t out_value;        // float: f
    md_gpu_addr_t out_gradient;     // float4: df/dx, df/dy, df/dz (per bohr), f
    md_gpu_addr_t out_failed;       // optional uint32: how many values could not be made (NaN), over all fields
} md_gto_iso_gpu_eval_desc_t;

// Issues the evaluation into 'stream': asynchronous, no readbacks and no waits. The set, the points and the
// outputs must stay until the stream has passed it. Scratch is ~52 bytes per point, up to 2^20 points at a
// time. One field's evaluation costs about as much as another's whatever the set: a set of several fields
// evaluates the functions of all of them at every point, for every field asked for. Returns false, logged,
// when nothing was issued (invalid input, a kernel this device could not make, out of device memory).
bool md_gto_iso_gpu_evaluate_launch(md_gpu_stream_t stream, const md_gto_iso_gpu_set_t* set, const md_gto_iso_gpu_eval_desc_t* desc);

// The same from and to the host, waiting for the result. xyz: num_xyz points, xyz_stride bytes apart (0 =
// packed float3). out_value [num_fields][num_xyz] and out_gradient [num_fields][num_xyz][3] (df/dx, df/dy,
// df/dz), either may be NULL. The points are evaluated in a spatially coherent order (a Morton curve
// through their box) and the results put back in theirs, so any order costs the same. Returns false, logged,
// when nothing was evaluated, or when some value could not be made (NaN there).
bool md_gto_iso_gpu_evaluate(md_gpu_stream_t stream, const md_gto_iso_gpu_set_t* set, const uint32_t* fields, uint32_t num_fields,
                             const float* xyz, size_t xyz_stride, size_t num_xyz, float* out_value, float* out_gradient);

// ---------------------------------------------------------------------------
// BENCHMARK (temporary, branch gto-iso-gpu-opt): a step by step profile of the extractions
// ---------------------------------------------------------------------------
// With a profile set, md_gto_iso_gpu_extract waits for the GPU after every step (so it runs slower than
// without) and appends one entry per step: the wall clock since the previous entry and what the step
// worked on. Before every vertex evaluation it also reads back how many vertices are still moving;
// that readback is not timed.
#define MD_GTO_ISO_GPU_PROFILE_CAP 512

typedef struct md_gto_iso_gpu_profile_step_t {
    const char* name;       // static
    double   ms;            // since the previous entry: the work and one wait for it
    uint32_t items;         // what it worked on: nodes in, bricks, vertices still moving, ...
    uint32_t aux;           // nodes kept, active 64-vertex groups, ...
} md_gto_iso_gpu_profile_step_t;

typedef struct md_gto_iso_gpu_profile_t {
    uint32_t num_steps;
    md_gto_iso_gpu_profile_step_t step[MD_GTO_ISO_GPU_PROFILE_CAP];
    double   wait_ms;       // one wait with nothing to wait for: the floor under every step's ms
    // work counts
    uint64_t brick_ao_sum;  // over the bricks: the AOs whose shells reach the brick
    uint32_t brick_ao_max;
    uint64_t brick_ao_sq;   // over the bricks: the square of that (dform's work)
    uint32_t num_ao, num_fn, fn_stride;
    uint32_t corner_cols;   // function columns the corner kernel multiplied
    uint32_t corner_points; // lattice points the register-form corner kernels evaluated (0: not counted)
    uint64_t kernels_used;  // bit per kernel launched (the print lists what the driver reports of them)
    bool     orbitals_only;
    // the groups the list kernels take (consecutive points, consecutive vertices): the AOs whose shells
    // reach a group's box, what its kernel lists (larger than a brick's when its points are far apart)
    uint64_t cgroup_ao_sum; // corner groups of cgroup_size points
    uint32_t cgroup_count, cgroup_size;
    uint64_t vgroup_ao_sum; // vertex groups of vgroup_size vertices
    uint32_t vgroup_count, vgroup_size;
    // the bricks given corners, by what the corners show (with ordered vertices; 0 otherwise)
    uint32_t brick_surf;    // with surface cells
    uint32_t brick_in;      // without: every corner inside a sheet (the bound's low side kept it)
    uint32_t brick_out;     // without: every corner outside (the high side kept it)
} md_gto_iso_gpu_profile_t;

// NULL turns it off. One profile for every extraction (not thread safe); each extraction restarts it.
void md_gto_iso_gpu_profile_set(md_gto_iso_gpu_profile_t* profile);

// The steps summed by name, in order of first appearance, as a text table (nul terminated, truncated
// to cap). Returns the length written.
size_t md_gto_iso_gpu_profile_print(const md_gto_iso_gpu_profile_t* profile, char* buf, size_t cap);

// The revision and the switches the extraction runs with, as one line
const char* md_gto_iso_gpu_profile_variant(void);

// The switches (also read from GTO_ISO_* in the environment at initialization)
typedef struct md_gto_iso_gpu_switches_t {
    int  eval_eps_exp;      // the AO cutoff 1e-N of the evaluations and the octree's culling (4..12, default 6)
    bool dform;             // sets of densities: corners and vertices as phi^T D phi (default off; on regardless past 1024 functions)
    int  flood;             // APPROXIMATE corners from seeds and a flood: 0 off, 1 densities whose maxima are at the atoms (default), 2 every field
} md_gto_iso_gpu_switches_t;

void md_gto_iso_gpu_switches(md_gto_iso_gpu_switches_t* out);
void md_gto_iso_gpu_switches_set(const md_gto_iso_gpu_switches_t* in);
#endif

#ifdef __cplusplus
}
#endif
