#pragma once

#include <stdint.h>
#include <stdbool.h>

#include <md_types.h>
#include <core/md_unit.h>

struct md_allocator_i;

// ATTRIBUTES
//
// Auxiliary data attached to a system or a state, keyed by an HDF5 style path. The table is a flat
// array of LEAVES sorted by path; a group exists only as a prefix shared by the leaves below it and
// is never declared. Prefixes match at segment boundaries: "atom" covers "atom/charge", never
// "atomic_number/z". The path carries the meaning, the format only the layout.
//
//     atom/charge/mulliken        rank 1 {N}     components 1   N scalars
//     atom/velocity               rank 1 {N}     components 3   N 3-vectors
//     qm/atom/normal_mode         rank 2 {M,N}   components 3   N 3-vectors per mode
//     dipole/ground_state/vector  rank 1 {1}     components 3   one 3-vector
//     scf/energy                  rank 0         components 1   a single value
//
// CONVENTIONS, none of which the table enforces:
//   - Related data is a SIBLING in the same group, never an extra axis: an origin beside its vector
//     (dipole/x/{vector,origin}), a coordinate beside the values it indexes (script/rdf/bin), three
//     charge schemes as atom/charge/{mulliken,hirshfeld,lowdin}. Siblings share an INDEX SPACE, not
//     a shape; a member with fewer axes is constant over the ones it lacks (broadcasting).
//   - "run/<name>" holds one trajectory and everything sampled along it; "run/<name>/time" is its
//     frame axis. See FRAME AXES.
//   - In "atom/..." the atom axis is the LAST index axis. Nothing checks the extent against the
//     system: a consumer indexing by atom checks it.
//   - Distinct paths do not imply distinct data (see md_attributes_alias). Compare datums with
//     md_attribute_same_data, never by id.
//
// The design rationale lives in docs/attributes.md.

#define MD_ATTRIBUTE_INVALID  ((md_attribute_id_t)0)
#define MD_ATTRIBUTE_MAX_RANK 5

// Hash of the full path; zero is invalid. Derived from the path alone, so it is stable across a
// reload and can be resolved before the table exists (md_attributes_id_from_path). It identifies
// an ADDRESS, not a datum.
typedef uint64_t md_attribute_id_t;

typedef enum md_attribute_type_t {
    MD_ATTRIBUTE_TYPE_NONE = 0,
    MD_ATTRIBUTE_TYPE_F32,
    MD_ATTRIBUTE_TYPE_F64,
    MD_ATTRIBUTE_TYPE_I8,
    MD_ATTRIBUTE_TYPE_U8,
    MD_ATTRIBUTE_TYPE_I16,
    MD_ATTRIBUTE_TYPE_U16,
    MD_ATTRIBUTE_TYPE_I32,
    MD_ATTRIBUTE_TYPE_U32,
    MD_ATTRIBUTE_TYPE_I64,
    MD_ATTRIBUTE_TYPE_U64,

    // UTF-8 text. The stored element is a 4 byte handle into the table's interning string pool, so
    // every layout rule holds exactly as for a number. Published as str_t[], read back through
    // md_attribute_str / md_attribute_extract_str. Never virtual, never extracted as a number.
    MD_ATTRIBUTE_TYPE_STR,

    MD_ATTRIBUTE_TYPE_COUNT,
} md_attribute_type_t;

typedef enum md_attribute_flags_t {
    MD_ATTRIBUTE_FLAG_NONE     = 0,

    // The OUTERMOST index axis is a frame axis: shape[0] equals the extent of the attribute's axis,
    // the nearest "time" at or above its group (see FRAME AXES). Verified on create, so publish the
    // axis first. A virtual temporal attribute can only be extracted a frame at a time.
    MD_ATTRIBUTE_FLAG_TEMPORAL = 1,

    // The INNERMOST index axis is the upper triangle of a symmetric N x N matrix, packed row major -
    // row i holds columns i..N-1 - so its extent is N(N+1)/2. Nothing about indexing, slicing or
    // extracting changes: an outer axis still indexes whole matrices, and the values are read like
    // any other. The flag is what says they ARE a matrix, and md_attribute_packed_symmetric_dim gives
    // N back. For matrices large enough that holding the lower half as well, or as double, matters -
    // AO density matrices are N^2 in the basis. Verified on create: rank >= 1, one component, and a
    // triangular innermost extent.
    MD_ATTRIBUTE_FLAG_PACKED_SYMMETRIC = 2,
} md_attribute_flags_t;

// LAYOUT. Two independent things:
//   type + components   WHAT ONE VALUE IS. components (>= 1, never defaulted) type sized
//                       components form one atomic value, which is never indexed or sub-set.
//   rank + shape        WHERE THE VALUES LIVE. Every axis is an index axis. Row major: the last
//                       axis varies fastest, and the components of a value are contiguous below it.
// {S,N} with 1 component (N scalars per state) and {N} with 3 (N vectors) are different formats;
// folding components into shape would make them the same spelling.
typedef struct md_attribute_format_t {
    md_attribute_type_t type;                          // what one component is
    uint32_t            components;                    // components of ONE value, >= 1
    uint32_t            rank;                          // number of index axes, 0 means a single value
    uint32_t            shape[MD_ATTRIBUTE_MAX_RANK];  // extent of each index axis, zero past rank
} md_attribute_format_t;

// A contiguous window of an attribute: the first num_idx axes are fixed, every later axis is taken
// whole. num_idx 0 selects the whole attribute. It carries nothing about any attribute, so one slice
// applies across a group (clamp num_idx to each member's rank). Passed by value.
typedef struct md_attribute_slice_t {
    uint32_t idx[MD_ATTRIBUTE_MAX_RANK];  // index along each fixed leading axis
    uint32_t num_idx;                     // how many leading axes are fixed
} md_attribute_slice_t;

static inline md_attribute_slice_t md_attribute_slice_all(void) {
    md_attribute_slice_t s = {0};
    return s;
}

static inline md_attribute_slice_t md_attribute_slice_1(uint32_t i) {
    md_attribute_slice_t s = {0};
    s.idx[0]  = i;
    s.num_idx = 1;
    return s;
}

static inline md_attribute_slice_t md_attribute_slice_2(uint32_t i, uint32_t j) {
    md_attribute_slice_t s = {0};
    s.idx[0]  = i;
    s.idx[1]  = j;
    s.num_idx = 2;
    return s;
}

typedef struct md_attribute_t md_attribute_t;

// What a provider that streams from disk reads through (md_attribute_io_read_at). Owned by an
// extraction context (md_system_extract_begin), which keeps files open across frames; NULL otherwise.
typedef struct md_attribute_io_t md_attribute_io_t;

// VIRTUAL ATTRIBUTES are computed on every extract by a provider; the table never caches the result.
//
// The provider is handed the slice the caller asked for (NULL for the whole attribute, otherwise
// num_idx >= 1) and writes cap elements of the STORED type in the STORED unit into dst; conversion
// happens centrally afterwards. It returns cap on success and 0 on failure.
//
// It may read other attributes, but the dependency graph must be acyclic: nothing detects a cycle.
// Scratch comes from the calling thread's temp allocator. io is NULL outside an extraction context.
typedef size_t (*md_attribute_provider_fn)(
    void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, md_attribute_io_t* io
);

typedef struct md_attribute_virtual_t {
    md_attribute_provider_fn provider;

    // Either BORROWED (NULL, the system, a trajectory owned elsewhere; user_data_size 0), or OWNED:
    // allocated with md_attributes_alloc_user_data and freed by size with the owning attribute.
    void*  user_data;
    size_t user_data_size;
} md_attribute_virtual_t;

// Reading: path, label, description, format, unit, flags and version are for anyone. data and virt
// are what md_attribute_view and the extract functions read; prefer those.
struct md_attribute_t {
    md_attribute_id_t       id;
    uint64_t                version;        // see md_attributes_version
    md_attribute_flags_t    flags;
    str_t                   path;           // full path, owned by the table
    str_t                   label;          // presentation only, may be empty and need not be unique
    str_t                   description;    // presentation only, may be empty
    md_attribute_format_t   format;
    md_unit_t               unit;           // (md_unit_t){0} is dimensionless, a valid value
    void*                   data;           // resident storage, md_attribute_byte_size() bytes; NULL when virtual
    md_attribute_virtual_t  virt;           // how to compute it; zeroed when resident

    // The attribute that owns the storage. Equal to id except for an alias, which is flattened on
    // creation and so always names a real owner.
    md_attribute_id_t       root;
};

// What a producer fills in to publish one attribute. Designated initialisers; zero is the right
// default for everything except format.components, where it is rejected.
//
// At most one of {data, virt}:
//   data      copied in; byte_size must equal md_attribute_byte_size(&format) (for STR, the count
//             times sizeof(str_t)). It is the guard against a buffer whose type or count disagrees
//             with the format, so compute it from the SOURCE buffer, not from the format.
//   neither   reserve zeroed storage, filled later through md_attributes_data; byte_size 0.
//   virt      a virtual attribute; byte_size 0. *virt is copied.
typedef struct md_attribute_desc_t {
    str_t                         path;         // "atom/charge/mulliken"
    md_attribute_format_t         format;
    md_attribute_flags_t          flags;
    md_unit_t                     unit;         // md_unit_none() when there is nothing to say
    str_t                         label;        // optional
    str_t                         description;  // optional
    const void*                   data;
    size_t                        byte_size;
    const md_attribute_virtual_t* virt;
} md_attribute_desc_t;

// The table. Set alloc before the first create; md_attributes_free releases everything and clears
// it. Scoped to one dataset: two tables holding the same path is the normal case. Fields other than
// alloc are private.
typedef struct md_attributes_t {
    struct md_allocator_i*   alloc;
    md_array(md_attribute_t) attr;          // sorted by path
    md_array(uint32_t)       id_index;      // open addressing over attr, holds index+1; 0 is empty

    // String pool behind MD_ATTRIBUTE_TYPE_STR: append only, interning, entry 0 is "".
    md_array(char)           str_data;
    md_array(uint32_t)       str_offset;
    md_array(uint32_t)       str_index;

    uint64_t                 version_counter;
} md_attributes_t;

// Walks the attributes at or below a prefix in path order. See md_attributes_iter.
typedef struct md_attribute_iter_t {
    const md_attribute_t*  attr;    // the current attribute
    str_t                  child;       // its first path segment below the prefix; empty for the prefix itself
    str_t                  child_path;  // prefix/child, the full path of that segment

    const md_attributes_t* _table;
    str_t                  _prefix;
    size_t                 _next;
    bool                   _children;
} md_attribute_iter_t;

#ifdef __cplusplus
extern "C" {
#endif

// FORMAT
// Functions of the format alone, so a producer can size a buffer before creating anything.

size_t md_attribute_type_size(md_attribute_type_t type);          // bytes of one component, 0 for NONE
bool   md_attribute_type_is_numeric(md_attribute_type_t type);    // any type but NONE and STR
size_t md_attribute_value_count(const md_attribute_format_t* format);    // product of shape, 1 at rank 0
size_t md_attribute_element_count(const md_attribute_format_t* format);  // value_count * components
size_t md_attribute_byte_size(const md_attribute_format_t* format);      // element_count * type_size

// ATTRIBUTE

// "atom/charge/mulliken" -> "atom/charge" and "mulliken". Views into attr->path.
str_t md_attribute_group(const md_attribute_t* attr);
str_t md_attribute_leaf(const md_attribute_t* attr);

// A second name for another attribute's datum.
static inline bool md_attribute_is_alias(const md_attribute_t* attr) {
    return attr->root != attr->id;
}

// Computed by a provider on every read, including an alias of a computed attribute.
static inline bool md_attribute_is_virtual(const md_attribute_t* attr) {
    return attr->virt.provider != NULL;
}

// Do two attributes read the same datum? One compare, no table.
bool md_attribute_same_data(const md_attribute_t* a, const md_attribute_t* b);

// ZERO COPY ACCESS. The stored values in place, when attr is resident and its format is exactly
// {type, components, rank}; NULL otherwise, including for a NULL attr, so it chains with find:
//
//     const uint32_t* idx = (const uint32_t*)md_attribute_view(md_attributes_find(t, path), MD_ATTRIBUTE_TYPE_U32, 1, 1);
//
// Nothing is converted, so the pointer is only as good as the type asked for. STR is never viewed.
// The pointer follows the invalidation rule of md_attributes_get.
const void* md_attribute_view(const md_attribute_t* attr, md_attribute_type_t type, uint32_t components, uint32_t rank);

// SLICES
// Sizes are functions of the format alone, so they can be asked of a virtual attribute.

// Elements (components folded in) the slice selects: what to allocate for. Zero when it does not
// apply - more indices than axes, or an index out of range.
size_t md_attribute_slice_count(const md_attribute_t* attr, md_attribute_slice_t slice);

// The format of what the slice yields: the fixed leading axes removed.
bool md_attribute_slice_format(md_attribute_format_t* out, const md_attribute_t* attr, md_attribute_slice_t slice);

// N for a packed symmetric extent of N(N+1)/2 (MD_ATTRIBUTE_FLAG_PACKED_SYMMETRIC); 0 when the extent
// is not a triangular number.
size_t md_attribute_packed_symmetric_dim(size_t extent);

// EXTRACTION
// Copies a slice into dst, converting the stored type and the stored unit to dst_unit. Returns the
// number of elements written, 0 on failure; cap must be at least md_attribute_slice_count.
//
//   - dst_unit md_unit_none() means "as stored". Otherwise the dimensions must match or the call
//     FAILS: a silently rescaled quantity is undetectable downstream.
//   - Conversion is done in double and narrowed once.
//   - Integers convert numerically. STR is refused (md_attribute_extract_str reads text).
//   - A virtual temporal attribute refuses slice_all: fix the frame axis.
//   - An index out of range is an error, never a clamp.
//
// f64 is not a convenience: quantum chemistry data is double at the boundary on purpose. Use f32 for
// what heads for a colour ramp or a plot, f64 for what heads back into a computation.
size_t md_attribute_extract_f32(float  dst[], size_t cap, const md_attribute_t* attr, md_attribute_slice_t slice, md_unit_t dst_unit);
size_t md_attribute_extract_f64(double dst[], size_t cap, const md_attribute_t* attr, md_attribute_slice_t slice, md_unit_t dst_unit);

// The same, into a buffer allocated from alloc (count * sizeof(float|double) bytes). NULL on
// failure, with nothing left allocated. out_count, when given, receives the element count.
float*  md_attribute_extract_alloc_f32(size_t* out_count, const md_attribute_t* attr, md_attribute_slice_t slice, md_unit_t dst_unit, struct md_allocator_i* alloc);
double* md_attribute_extract_alloc_f64(size_t* out_count, const md_attribute_t* attr, md_attribute_slice_t slice, md_unit_t dst_unit, struct md_allocator_i* alloc);

// Text. The table is needed because the pool lives there. The returned strings are NUL terminated
// views into the pool, and follow the invalidation rule of md_attributes_get.
str_t  md_attribute_str(const md_attributes_t* attributes, const md_attribute_t* attr, size_t index);
size_t md_attribute_extract_str(str_t dst[], size_t cap, const md_attributes_t* attributes, const md_attribute_t* attr);

// TABLE

void   md_attributes_free (md_attributes_t* attributes);
size_t md_attributes_count(const md_attributes_t* attributes);

// The id a path has, whether or not anything is published there.
md_attribute_id_t md_attributes_id_from_path(str_t path);

// Private state for a provider (md_attribute_virtual_t.user_data) that the table frees by size with
// the owning attribute. NULL for size 0 or an unset allocator.
void* md_attributes_alloc_user_data(md_attributes_t* attributes, size_t size);

// Publishes an attribute. Returns its id, MD_ATTRIBUTE_INVALID on failure. Everything in desc is
// copied. Rejects: an unset allocator; an empty path, a leading or trailing '/', an empty segment;
// a duplicate path; type NONE; components 0; rank above MD_ATTRIBUTE_MAX_RANK; a zero extent or one
// set beyond rank; a byte_size disagreeing with the format; virt without a provider, or with data;
// a virtual STR; TEMPORAL without a matching frame axis.
md_attribute_id_t md_attributes_create(md_attributes_t* attributes, const md_attribute_desc_t* desc);

// create, replacing whatever is at desc->path (and its aliases). What a producer that can run twice
// wants. The id is unchanged by the replacement.
md_attribute_id_t md_attributes_replace(md_attributes_t* attributes, const md_attribute_desc_t* desc);

// A second path for an existing datum: shares the target's storage or provider, format, unit and
// version; path, label and description are its own. Aliasing an alias points at the owner. An alias
// cannot outlive its target: removing the target removes its aliases. MD_ATTRIBUTE_INVALID when the
// target is unknown or the path is bad or taken.
md_attribute_id_t md_attributes_alias(md_attributes_t* attributes, md_attribute_id_t target, str_t path, str_t label, str_t description);

// Removes one attribute and its aliases.
bool md_attributes_remove(md_attributes_t* attributes, md_attribute_id_t id);

// Removes everything at or below prefix, and aliases of it elsewhere. An empty prefix is refused.
// Returns the number removed, aliases included.
size_t md_attributes_remove_prefix(md_attributes_t* attributes, str_t prefix);

// Publishes a per atom float column, rank 1 {count} with 'components' floats per atom, interleaved.
// SKIPPED (returns MD_ATTRIBUTE_INVALID) when every value is identical or every value is NAN, so a
// loader may call it for every column its format defines. NAN marks "no value for this atom" and is
// carried through extraction bit for bit. mdlib builds with fast math, so a consumer compiled the
// same way cannot test for NAN with v != v or isnan().
md_attribute_id_t md_attributes_publish_atom_column(md_attributes_t* attributes, str_t path, md_unit_t unit, uint32_t components, const float values[], size_t count);

// LOOKUP
// Pointers are INVALIDATED by the next create, replace, alias or remove. Hold the id, not the pointer.

const md_attribute_t* md_attributes_get (const md_attributes_t* attributes, md_attribute_id_t id);
const md_attribute_t* md_attributes_find(const md_attributes_t* attributes, str_t path);

// find("<group>/<leaf>"); an empty group is find(leaf). For run relative paths and siblings.
const md_attribute_t* md_attributes_find_in(const md_attributes_t* attributes, str_t group, str_t leaf);

// The attribute named leaf in attr's own group: the origin beside a vector, the bin beside a series.
const md_attribute_t* md_attributes_sibling(const md_attributes_t* attributes, const md_attribute_t* attr, str_t leaf);

// Every attribute at or below prefix, in path order. An optional trailing '/' is ignored and an
// empty prefix matches everything. Attributes must not be added or removed while iterating.
//
//     for (md_attribute_iter_t it = md_attributes_iter(t, STR_LIT("atom")); md_attributes_next(&it);) {
//         use(it.attr);
//     }
md_attribute_iter_t md_attributes_iter(const md_attributes_t* attributes, str_t prefix);

// The immediate children of prefix, once each: "dipole" yields "ground_state", "transition", ...
// it.child is the child's name, it.child_path its full path ("dipole/ground_state"), and it.attr
// the first attribute below it.
md_attribute_iter_t md_attributes_iter_children(const md_attributes_t* attributes, str_t prefix);

bool md_attributes_next(md_attribute_iter_t* it);

// Ids of every attribute at or below prefix, in path order: a snapshot, for when the table is about
// to change. Returns the total number of matches and writes at most cap of them.
size_t md_attributes_query(md_attribute_id_t out_ids[], size_t cap, const md_attributes_t* attributes, str_t prefix);

// WRITING AND INVALIDATION
//
// A version belongs to the DATUM: every alias reports its target's version. 0 means no such
// attribute. create and replace bump it; a write through md_attributes_data does NOT, because the
// producer may take a long time to fill the buffer - it calls md_attributes_touch when done.
uint64_t md_attributes_version(const md_attributes_t* attributes, md_attribute_id_t id);

// "The contents changed": bumps the attribute and all its aliases to one new version and returns it
// (0 if unknown). Not thread safe; a parallel fill touches once, from one thread, when complete.
uint64_t md_attributes_touch(md_attributes_t* attributes, md_attribute_id_t id);

// Writable storage of a resident, non alias, non STR attribute whose type is expected_type; NULL
// otherwise. Same invalidation rule as md_attributes_get.
void* md_attributes_data(md_attributes_t* attributes, md_attribute_id_t id, md_attribute_type_t expected_type);

// FRAME AXES
//
// A temporal attribute's frame axis is the NEAREST attribute named "time" at or above its OWNER's
// group that is itself an axis: temporal, rank 1, one component, numeric. For "run/a/atom/position"
// that is "run/a/atom/time" if it exists, else "run/a/time", and so on. An axis is its own axis.
// Several axes coexist: an energy file sampled at its own rate keeps its own "time".

// NULL when attr is not temporal or its axis has been removed. Same invalidation rule as get.
const md_attribute_t* md_attributes_axis(const md_attributes_t* attributes, const md_attribute_t* attr);

// Where frame src_index of src_axis lands on dst_axis, matched by coordinate VALUE (converted to
// dst_axis' unit) within float precision or a thousandth of the local spacing. Axes must be non
// decreasing; an axis without a unit holds ordinals and only matches another such axis. False, and
// out_index untouched, when nothing matches.
bool md_attribute_axis_map(size_t* out_index, const md_attribute_t* src_axis, size_t src_index, const md_attribute_t* dst_axis);

// For a provider: read bytes at offset from the file at path, through io's cache of open files when
// io is non NULL. Returns the number of bytes read.
size_t md_attribute_io_read_at(md_attribute_io_t* io, str_t path, int64_t offset, void* dst, size_t bytes);

#ifdef __cplusplus
}
#endif
