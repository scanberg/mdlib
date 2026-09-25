#include <md_system.h>
#include <md_nonbonded.h>

#include <inttypes.h>
#include <stdio.h>

#include <core/md_log.h>
#include <core/md_simd.h>
#include <core/md_array.h>
#include <core/md_hash.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>

#ifdef __cplusplus
extern "C" {
#endif

void md_system_free(md_system_t* sys) {
    ASSERT(sys);
    ASSERT(sys->alloc);
    md_allocator_i* alloc = sys->alloc;

    if (sys->nonbonded) {
        md_nb_forcefield_free(sys->nonbonded);
        md_free(alloc, sys->nonbonded, sizeof(md_nb_forcefield_t));
        sys->nonbonded = NULL;
    }

    // ATOM
    md_array_free(sys->atom.type_idx, alloc);
    md_array_free(sys->atom.flags, alloc);

    // ATOM TYPE
    md_array_free(sys->atom.type.name, alloc);
    for (size_t i = 0; i < md_array_size(sys->atom.type.ff_type); ++i) {
        if (sys->atom.type.ff_type[i].ptr) str_free(sys->atom.type.ff_type[i], alloc);
    }
    md_array_free(sys->atom.type.ff_type, alloc);
    md_array_free(sys->atom.type.z, alloc);
    md_array_free(sys->atom.type.mass, alloc);
    md_array_free(sys->atom.type.radius, alloc);
    md_array_free(sys->atom.type.color, alloc);
    md_array_free(sys->atom.type.flags, alloc);

    // COMPONENT
    md_array_free(sys->component.name, alloc);
    md_array_free(sys->component.seq_id, alloc);
    md_array_free(sys->component.atom_offset, alloc);
    md_array_free(sys->component.flags, alloc);

    // INSTANCE
    md_array_free(sys->instance.id, alloc);
    md_array_free(sys->instance.auth_id, alloc);
    md_array_free(sys->instance.comp_offset, alloc);
    md_array_free(sys->instance.entity_idx, alloc);

    // ENTITY
    md_array_free(sys->entity.id, alloc);
    md_array_free(sys->entity.flags, alloc);
    for (size_t i = 0; i < sys->entity.count; ++i) {
        if (!str_empty(sys->entity.description[i])) {
            str_free(sys->entity.description[i], alloc);
        }
    }
    md_array_free(sys->entity.description, alloc);

    // PROTEIN BACKBONE
    md_array_free(sys->protein_backbone.range.offset, alloc);
    md_array_free(sys->protein_backbone.range.inst_idx, alloc);
    md_array_free(sys->protein_backbone.segment.atoms, alloc);
    md_array_free(sys->protein_backbone.segment.angle, alloc);
    md_array_free(sys->protein_backbone.segment.secondary_structure, alloc);
    md_array_free(sys->protein_backbone.segment.rama_type, alloc);
    md_array_free(sys->protein_backbone.segment.comp_idx, alloc);

    // NUCLEIC BACKBONE
    md_array_free(sys->nucleic_backbone.range.offset, alloc);
    md_array_free(sys->nucleic_backbone.range.inst_idx, alloc);
    md_array_free(sys->nucleic_backbone.segment.atoms, alloc);
    md_array_free(sys->nucleic_backbone.segment.comp_idx, alloc);

    // BONDS
    md_array_free(sys->bond.pairs, alloc);
    md_array_free(sys->bond.flags, alloc);
    md_array_free(sys->bond.conn.atom_idx, alloc);
    md_array_free(sys->bond.conn.bond_idx, alloc);
    md_array_free(sys->bond.conn.offset, alloc);

    // HYDROGEN BONDS
    md_array_free(sys->hydrogen_bond.candidate.acceptor.idx, alloc);
    md_array_free(sys->hydrogen_bond.candidate.acceptor.num_lone_pairs, alloc);
    md_array_free(sys->hydrogen_bond.candidate.donor.d_idx, alloc);
    md_array_free(sys->hydrogen_bond.candidate.donor.h_idx, alloc);

    md_index_data_free(&sys->ring);

    md_array_free(sys->structure.offset, alloc);
    md_array_free(sys->structure.atom_idx, alloc);
    md_array_free(sys->structure.parent_idx, alloc);
    md_array_free(sys->structure.atom_slot, alloc);

    // ASSEMBLY
    md_array_free(sys->assembly.atom_range, alloc);
    md_array_free(sys->assembly.label, alloc);
    md_array_free(sys->assembly.transform, alloc);

    if (!str_empty(sys->description)) {
        str_free(sys->description, alloc);
    }

    // REFERENCE STATE
    md_system_state_free(&sys->reference);

    // ATTRIBUTES
    md_attributes_free(&sys->attributes);

    MEMSET(sys, 0, sizeof(md_system_t));
}

bool md_system_state_init(md_system_state_t* state, size_t num_atoms) {
    ASSERT(state);

    if (!state->alloc) {
        MD_LOG_ERROR("State allocator not set");
        return false;
    }

    md_allocator_i* alloc = state->alloc;
    md_system_state_free(state);
    state->alloc = alloc;

    // The table allocates through its own handle, exactly as a system's does, so a producer can
    // publish into a freshly initialised state without a second setup step.
    state->attributes.alloc = alloc;

    // A freshly initialised state did not come from a run. Only md_system_extract_frame and the
    // interpolation which produces a state write a non negative frame. See md_system_state_t.
    state->frame = -1.0;

    if (num_atoms == 0) {
        return true;
    }

    const size_t capacity = ALIGN_TO(num_atoms, 16);

    md_array_resize(state->xyz, capacity, alloc);

    // Zero the padding past num_atoms. The capacity is rounded up so vectorised code may load whole
    // groups of atoms; an uninitialised tail puts garbage floats into those lanes.
    const size_t tail_bytes = (capacity - num_atoms) * sizeof(vec3_t);
    if (tail_bytes > 0) {
        MEMSET(state->xyz + num_atoms, 0, tail_bytes);
    }

    state->num_atoms = num_atoms;

    return true;
}

void md_system_state_free(md_system_state_t* state) {
    ASSERT(state);

    // A view owns nothing; zeroing it is the whole job.
    if (state->alloc) {
        md_array_free(state->xyz, state->alloc);
        md_attributes_free(&state->attributes);
    }
    md_allocator_i* alloc = state->alloc;
    MEMSET(state, 0, sizeof(md_system_state_t));
    state->alloc = alloc;
}

bool md_system_state_copy(md_system_state_t* dst, const md_system_state_t* src) {
    ASSERT(dst);
    ASSERT(src);

    if (dst == src) {
        return true;
    }
    if (!md_system_state_init(dst, src->num_atoms)) {
        return false;
    }
    if (src->num_atoms > 0 && src->xyz) {
        MEMCPY(dst->xyz, src->xyz, src->num_atoms * sizeof(vec3_t));
    }
    dst->unitcell = src->unitcell;
    dst->frame    = src->frame;

    // The attribute table is deliberately NOT carried over. The one caller is md_util_system_infer
    // writing sys->reference, and topology inference reads coordinates and the cell - nothing else.
    // Duplicating a frame's velocities into a reference snapshot would put a second copy of them in
    // the system with nothing keeping the two in step, which is the whole class of problem the
    // reference field's own rule exists to prevent.
    return true;
}

void md_system_reset(md_system_t* sys) {
    ASSERT(sys);
    md_allocator_i* alloc = sys->alloc;
    md_system_free(sys);
    sys->alloc = alloc;
    // The table allocates through its own handle, so every loader gets a usable table without
    // having to remember to wire this up itself.
    sys->attributes.alloc = alloc;
}

static void build_connectivity(md_bond_conn_data_t* conn, const md_atom_pair_t* bond_pairs, size_t bond_pair_count, size_t atom_count, md_allocator_i* alloc) {
    ASSERT(conn);
    ASSERT(alloc);

    if (bond_pairs == NULL) return;
    if (bond_pair_count == 0) return;
    if (atom_count == 0) return;

    conn->offset_count = atom_count + 1;
    md_array_resize(conn->offset, conn->offset_count, alloc);
    MEMSET(conn->offset, 0, md_array_bytes(conn->offset));

    // This have length of 2 * bond_count (one for each direction of the bond)
    conn->count = 2 * bond_pair_count;
    md_array_resize(conn->atom_idx, conn->count, alloc);
    md_array_resize(conn->bond_idx, conn->count, alloc);

    typedef struct {
        uint16_t off[2];
    } offset_t;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);

    offset_t* local_offset = md_temp_alloc_zero_array(temp, offset_t, bond_pair_count);
    ASSERT(local_offset);

    // Two packed 16-bit local offsets for each of the bond idx
    // Use offsets as accumulators for length
    for (size_t i = 0; i < bond_pair_count; ++i) {
		ASSERT(bond_pairs[i].idx[0] < (md_atom_idx_t)atom_count);
		ASSERT(bond_pairs[i].idx[1] < (md_atom_idx_t)atom_count);
        local_offset[i].off[0] = (uint16_t)conn->offset[bond_pairs[i].idx[0]]++;
        local_offset[i].off[1] = (uint16_t)conn->offset[bond_pairs[i].idx[1]]++;
    }

    // Compute complete edge offsets (exclusive scan)
    uint32_t off = 0;
    for (size_t i = 0; i < conn->offset_count; ++i) {
        const uint32_t len = conn->offset[i];
        conn->offset[i] = off;
        off += len;
    }

    // Write edge indices to correct location
    for (size_t i = 0; i < bond_pair_count; ++i) {
        const md_atom_pair_t p = bond_pairs[i];
        const int atom_a = p.idx[0];
        const int atom_b = p.idx[1];
        const int local_a = (int)local_offset[i].off[0];
        const int local_b = (int)local_offset[i].off[1];
        const int off_a = conn->offset[atom_a];
        const int off_b = conn->offset[atom_b];

        const int idx_a = off_a + local_a;
        const int idx_b = off_b + local_b;

        ASSERT(idx_a < (int)conn->count);
        ASSERT(idx_b < (int)conn->count);

        // Store the cross references to the 'other' atom index signified by the bond in the correct location
        conn->atom_idx[idx_a] = atom_b;
        conn->atom_idx[idx_b] = atom_a;

        conn->bond_idx[idx_a] = (md_bond_idx_t)i;
        conn->bond_idx[idx_b] = (md_bond_idx_t)i;
    }

    md_temp_end(temp);
}

void md_bond_build_connectivity(md_bond_data_t* in_out_bond, size_t atom_count, md_allocator_i* alloc) {
    ASSERT(in_out_bond);
    ASSERT(alloc);
	build_connectivity(&in_out_bond->conn, in_out_bond->pairs, in_out_bond->count, atom_count, alloc);
}

void md_system_bond_build_connectivity(md_system_t* sys) {
    ASSERT(sys);
	md_bond_build_connectivity(&sys->bond, sys->atom.count, sys->alloc);
}



// ATTRIBUTES

static const size_t attr_type_size[MD_ATTRIBUTE_TYPE_COUNT] = {
    [MD_ATTRIBUTE_TYPE_NONE] = 0,
    [MD_ATTRIBUTE_TYPE_F32]  = 4,
    [MD_ATTRIBUTE_TYPE_F64]  = 8,
    [MD_ATTRIBUTE_TYPE_I8]   = 1,
    [MD_ATTRIBUTE_TYPE_U8]   = 1,
    [MD_ATTRIBUTE_TYPE_I16]  = 2,
    [MD_ATTRIBUTE_TYPE_U16]  = 2,
    [MD_ATTRIBUTE_TYPE_I32]  = 4,
    [MD_ATTRIBUTE_TYPE_U32]  = 4,
    [MD_ATTRIBUTE_TYPE_I64]  = 8,
    [MD_ATTRIBUTE_TYPE_U64]  = 8,
    // The handle, not the text. That is the whole trick: the element stays fixed width, so every
    // layout rule keeps working, and the variable length part lives in the pool.
    [MD_ATTRIBUTE_TYPE_STR]  = 4,
};

// Absence is decided on the BIT PATTERN and never on a float comparison - see the note in
// md_attributes_publish_atom_column for why every float spelling of this test is unusable in a
// build with /fp:fast or -ffast-math. Internal, and staying internal: what a NAN MEANS is the
// producer's and the consumer's business, and the table only has to avoid destroying it.
static inline uint32_t attr_f32_bits(float v) {
    uint32_t u;
    MEMCPY(&u, &v, sizeof(u));
    return u;
}

static inline bool attr_f32_bits_absent(uint32_t bits) {
    return (bits & 0x7fffffffu) > 0x7f800000u;   // exponent all ones and a non zero mantissa
}

size_t md_attribute_type_size(md_attribute_type_t type) {
    if (type <= MD_ATTRIBUTE_TYPE_NONE || type >= MD_ATTRIBUTE_TYPE_COUNT) {
        return 0;
    }
    return attr_type_size[type];
}

// components lives in the format rather than in a trailing axis, so this is a read and not an
// interpretation. A format which reached the table has been validated to hold at least 1.
size_t md_attribute_components(const md_attribute_format_t* format) {
    ASSERT(format);
    return (size_t)format->components;
}

// Every axis in shape is an index axis, so this is the whole product. Rank 0 is the empty
// product, which is 1 - the single value case falls out rather than being branched on.
size_t md_attribute_value_count(const md_attribute_format_t* format) {
    ASSERT(format);
    size_t count = 1;
    for (uint32_t i = 0; i < format->rank; ++i) {
        count *= (size_t)format->shape[i];
    }
    return count;
}

size_t md_attribute_element_count(const md_attribute_format_t* format) {
    ASSERT(format);
    return md_attribute_value_count(format) * (size_t)format->components;
}

size_t md_attribute_byte_size(const md_attribute_format_t* format) {
    ASSERT(format);
    return md_attribute_element_count(format) * md_attribute_type_size(format->type);
}

// Splits at the LAST separator. loc is the separator offset, false if there is none.
static bool attr_path_split(size_t* loc, str_t path) {
    return str_rfind_char(loc, path, '/');
}

str_t md_attribute_group(const md_attribute_t* attr) {
    ASSERT(attr);
    size_t loc;
    if (!attr_path_split(&loc, attr->path)) {
        return (str_t){0};
    }
    return str_substr(attr->path, 0, loc);
}

str_t md_attribute_leaf(const md_attribute_t* attr) {
    ASSERT(attr);
    size_t loc;
    if (!attr_path_split(&loc, attr->path)) {
        return attr->path;
    }
    return str_substr(attr->path, loc + 1, SIZE_MAX);
}

// The shared core, generated once per destination type. first and count are in ELEMENTS
// (components already folded in), so every entry point differs only in how it computes the window -
// which is the point: there is exactly one place where a stored type becomes a number and a unit
// becomes a factor.
//
// f32 and f64 destinations exist because both are load bearing. A colour ramp or a plot wants
// floats; AO coefficients and total energies are double at the boundary on purpose, and narrowing
// them to ask a question and widening them again would throw away the precision that boundary
// exists to keep.
// first/count/slice all describe the SAME window: first and count are its resolved offset and
// length in elements, slice is the (possibly NULL) selection that produced them, forwarded on
// unchanged so a virtual attribute's provider sees exactly what the caller asked for rather than
// an offset it would have to reverse back into indices.
#define MD_ATTR_DEFINE_EXTRACT_RANGE(SUFFIX, DST_T)                                                 \
static size_t attr_extract_range_##SUFFIX(DST_T dst[], size_t cap, const md_attribute_t* attr,      \
                                          size_t first, size_t count, const md_attribute_slice_t* slice, md_unit_t dst_unit, md_attribute_io_t* io) { \
    ASSERT(attr);                                                                                   \
                                                                                                    \
    if (!dst) {                                                                                     \
        return 0;                                                                                   \
    }                                                                                               \
    if (count > cap) {                                                                              \
        MD_LOG_ERROR("Attribute '" STR_FMT "' needs %zu values, %zu supplied", STR_ARG(attr->path), count, cap); \
        return 0;                                                                                   \
    }                                                                                               \
    if (count == 0) {                                                                               \
        return 0;                                                                                   \
    }                                                                                               \
    /* Refused, not converted: a pool handle reinterpreted as a number is the failure where a  */    \
    /* wrong answer still looks like data. md_attribute_extract_str is the way to read these.  */    \
    if (attr->format.type == MD_ATTRIBUTE_TYPE_STR) {                                               \
        MD_LOG_ERROR("Attribute '" STR_FMT "' is a string; read it with md_attribute_extract_str", STR_ARG(attr->path)); \
        return 0;                                                                                   \
    }                                                                                               \
                                                                                                    \
    /* md_unit_none() as the target means "as stored". Anything else has to be convertible, and a */ \
    /* refusal here is the point: a silently rescaled quantity is not detectable downstream.      */ \
    double factor = 1.0;                                                                            \
    if (!md_unit_is_none(dst_unit) && !md_unit_conversion_factor(&factor, attr->unit, dst_unit)) {   \
        char from[64], to[64];                                                                      \
        size_t from_len = md_unit_print(from, sizeof(from), attr->unit);                            \
        size_t to_len   = md_unit_print(to,   sizeof(to),   dst_unit);                              \
        MD_LOG_ERROR("Attribute '" STR_FMT "' is '%.*s' and cannot be expressed as '%.*s'",         \
            STR_ARG(attr->path), (int)from_len, from, (int)to_len, to);                             \
        return 0;                                                                                   \
    }                                                                                               \
                                                                                                    \
    /* A computed attribute has no attr->data: it is read through its provider into a scratch     */ \
    /* buffer of the STORED type, and from there on is indistinguishable from a resident one - the*/ \
    /* conversion below is shared by both.                                                        */ \
    /*                                                                                            */ \
    /* The test is the PROVIDER and not the storage tag, and that distinction is load bearing: an */ \
    /* ALIAS is a second NAME for a datum, not a different kind of storage, so how it is read is  */ \
    /* its target's business - which is exactly why it inherits both 'data' and 'virt'. Switching */ \
    /* on storage sent an alias of a computed attribute down the resident path, where it found no */ \
    /* data and returned nothing at all. */                                                          \
    /* When the stored type IS the destination type and nothing is rescaled, the scratch buffer  */ \
    /* would be written once and copied once for nothing: the provider writes straight into dst.  */ \
    /* For a coordinate array that copy was 5-8% of decoding the frame it came from.              */ \
    if (attr->virt.provider && attr->format.type == MD_ATTRIBUTE_TYPE_##SUFFIX && factor == 1.0) {   \
        size_t written = attr->virt.provider(dst, count, attr, slice, attr->virt.user_data, io);    \
        if (written != count) {                                                                     \
            MD_LOG_ERROR("Attribute '" STR_FMT "' provider wrote %zu of %zu requested values", STR_ARG(attr->path), written, count); \
            return 0;                                                                               \
        }                                                                                           \
        return count;                                                                               \
    }                                                                                               \
    const void* src = NULL;                                                                         \
    md_temp_scope_t temp = {0};                                                                     \
    bool own_temp = false;                                                                          \
    if (attr->virt.provider) {                                                                      \
        temp = md_temp_begin();                                                                     \
        own_temp = true;                                                                            \
        void* buf = md_temp_alloc(temp, count * md_attribute_type_size(attr->format.type));         \
        if (!buf) {                                                                                 \
            /* count is bounded by the caller's own cap, so this is a real allocation failure and  */ \
            /* not a runaway slice - but the provider contract says dst is cap elements, and it is */ \
            /* entitled to write without checking. */                                               \
            MD_LOG_ERROR("Failed to allocate scratch for %zu values of virtual attribute '" STR_FMT "'", count, STR_ARG(attr->path)); \
            md_temp_end(temp);                                                                      \
            return 0;                                                                               \
        }                                                                                           \
        size_t written = attr->virt.provider(buf, count, attr, slice, attr->virt.user_data, io);    \
        if (written != count) {                                                                     \
            MD_LOG_ERROR("Attribute '" STR_FMT "' provider wrote %zu of %zu requested values", STR_ARG(attr->path), written, count); \
            md_temp_end(temp);                                                                      \
            return 0;                                                                               \
        }                                                                                           \
        src = buf;                                                                                  \
    } else {                                                                                        \
        if (attr->storage == MD_ATTRIBUTE_STORAGE_VIRTUAL) {                                        \
            MD_LOG_ERROR("Attribute '" STR_FMT "' is virtual but has no provider", STR_ARG(attr->path)); \
            return 0;                                                                               \
        }                                                                                           \
        if (!attr->data) {                                                                          \
            /* Reserved but never filled in. Silence here is what hid the alias bug above, so say */ \
            /* it: a caller getting 0 back has no other way to tell this from an empty slice. */      \
            MD_LOG_ERROR("Attribute '" STR_FMT "' has no data to read", STR_ARG(attr->path));        \
            return 0;                                                                               \
        }                                                                                           \
        src = (const uint8_t*)attr->data + first * md_attribute_type_size(attr->format.type);       \
    }                                                                                               \
                                                                                                    \
    size_t result = count;                                                                          \
    /* The stored type already matching the destination, with nothing to rescale, is a memcpy. */   \
    if (attr->format.type == MD_ATTRIBUTE_TYPE_##SUFFIX && factor == 1.0) {                         \
        MEMCPY(dst, src, count * sizeof(DST_T));                                                    \
    } else {                                                                                        \
        /* The scale is applied in double and narrowed once, at the end. */                         \
        switch (attr->format.type) {                                                                \
        case MD_ATTRIBUTE_TYPE_F32: MD_ATTR_CONVERT(float,    DST_T); break;                        \
        case MD_ATTRIBUTE_TYPE_F64: MD_ATTR_CONVERT(double,   DST_T); break;                        \
        case MD_ATTRIBUTE_TYPE_I8:  MD_ATTR_CONVERT(int8_t,   DST_T); break;                        \
        case MD_ATTRIBUTE_TYPE_U8:  MD_ATTR_CONVERT(uint8_t,  DST_T); break;                        \
        case MD_ATTRIBUTE_TYPE_I16: MD_ATTR_CONVERT(int16_t,  DST_T); break;                        \
        case MD_ATTRIBUTE_TYPE_U16: MD_ATTR_CONVERT(uint16_t, DST_T); break;                        \
        case MD_ATTRIBUTE_TYPE_I32: MD_ATTR_CONVERT(int32_t,  DST_T); break;                        \
        case MD_ATTRIBUTE_TYPE_U32: MD_ATTR_CONVERT(uint32_t, DST_T); break;                        \
        case MD_ATTRIBUTE_TYPE_I64: MD_ATTR_CONVERT(int64_t,  DST_T); break;                        \
        case MD_ATTRIBUTE_TYPE_U64: MD_ATTR_CONVERT(uint64_t, DST_T); break;                        \
        default:                                                                                    \
            MD_LOG_ERROR("Attribute '" STR_FMT "' has no readable type", STR_ARG(attr->path));      \
            result = 0;                                                                             \
            break;                                                                                  \
        }                                                                                           \
    }                                                                                               \
                                                                                                    \
    if (own_temp) {                                                                                 \
        md_temp_end(temp);                                                                          \
    }                                                                                               \
                                                                                                    \
    return result;                                                                                  \
}

#define MD_ATTR_CONVERT(SRC_T, DST_T)                               \
    do {                                                            \
        const SRC_T* s = (const SRC_T*)src;                         \
        for (size_t i = 0; i < count; ++i) {                        \
            dst[i] = (DST_T)((double)s[i] * factor);                \
        }                                                           \
    } while (0)

// Only the two the header exposes. The integer destinations were instantiated as well, and
// nothing ever called them: a consumer wanting an index column takes it through f64, which is
// exact to 2^53 and therefore lossless for every integer type an attribute can hold.
MD_ATTR_DEFINE_EXTRACT_RANGE(F32, float)
MD_ATTR_DEFINE_EXTRACT_RANGE(F64, double)

#undef MD_ATTR_CONVERT
#undef MD_ATTR_DEFINE_EXTRACT_RANGE

// Turns a slice into the contiguous window it selects. Row major, so fixing the FIRST num_idx axes
// selects one block whose size is the product of the axes left free times the components of a
// value; a slice fixing nothing leaves every axis free and the window is the whole attribute.
//
// This is the one piece of layout arithmetic in the library, and it is here rather than in each
// caller because it is where a wrong answer still looks like data. It reads only the FORMAT, never
// the storage, which is why a slice can be sized against an attribute whose data is not resident.
static bool attr_slice_window(size_t* out_first, size_t* out_count, const md_attribute_t* attr, const md_attribute_slice_t* slice) {
    ASSERT(attr);

    const md_attribute_format_t* fmt = &attr->format;
    const uint32_t num_idx = slice ? slice->num_idx : 0;

    if (num_idx > fmt->rank) {
        MD_LOG_ERROR("Attribute '" STR_FMT "' has rank %u, %u indices supplied", STR_ARG(attr->path), fmt->rank, num_idx);
        return false;
    }

    size_t block = (size_t)fmt->components;
    for (uint32_t i = num_idx; i < fmt->rank; ++i) {
        block *= (size_t)fmt->shape[i];
    }

    // Horner over the fixed axes gives the block ordinal; the block size turns it into elements.
    size_t ordinal = 0;
    for (uint32_t i = 0; i < num_idx; ++i) {
        if (slice->idx[i] >= fmt->shape[i]) {
            MD_LOG_ERROR("Attribute '" STR_FMT "': index %u out of range on axis %u of extent %u",
                STR_ARG(attr->path), slice->idx[i], i, fmt->shape[i]);
            return false;
        }
        ordinal = ordinal * (size_t)fmt->shape[i] + (size_t)slice->idx[i];
    }

    *out_first = ordinal * block;
    *out_count = block;
    return true;
}

size_t md_attribute_slice_count(const md_attribute_t* attr, const md_attribute_slice_t* slice) {
    ASSERT(attr);
    size_t first, count;
    return attr_slice_window(&first, &count, attr, slice) ? count : 0;
}

bool md_attribute_slice_format(md_attribute_format_t* out, const md_attribute_t* attr, const md_attribute_slice_t* slice) {
    ASSERT(out);
    ASSERT(attr);

    const md_attribute_format_t* fmt = &attr->format;
    const uint32_t num_idx = slice ? slice->num_idx : 0;

    size_t first, count;
    if (!attr_slice_window(&first, &count, attr, slice)) {
        return false;
    }

    // The value is untouched by slicing - fixing an index picks values, it never splits one - so
    // only the index axes change.
    MEMSET(out, 0, sizeof(*out));
    out->type       = fmt->type;
    out->components = fmt->components;
    out->rank       = fmt->rank - num_idx;
    for (uint32_t i = 0; i < out->rank; ++i) {
        out->shape[i] = fmt->shape[num_idx + i];
    }
    return true;
}

// "The whole attribute" of a temporal one is every frame of it, and whether that is reasonable
// depends on what producing it COSTS, not on the frame axis itself.
//
// Resident: the bytes already exist and the extract is a copy of an array the caller could have
// memcpy'd. 'frame/time' is exactly this - a plot legitimately wants every frame time, and it is a
// few hundred kilobytes.
//
// Virtual: every frame has to be produced, which for a per atom quantity means decoding the whole
// trajectory into one buffer. That is never what a caller meant to type, and refusing it is the
// point of the rule. An ALIAS is read through whatever it inherited, so it is judged by its
// provider, not by its storage tag - same reason the extract branches on the provider throughout.
static bool attr_reject_whole_temporal(const md_attribute_t* attr) {
    if ((attr->flags & MD_ATTRIBUTE_FLAG_TEMPORAL) && attr->virt.provider) {
        MD_LOG_ERROR("Attribute '" STR_FMT "' is temporal and computed on demand: fix the frame axis with a slice rather than asking for every frame",
            STR_ARG(attr->path));
        return true;
    }
    return false;
}

size_t md_attribute_extract_f32(float dst[], size_t cap, const md_attribute_t* attr, md_unit_t dst_unit) {
    ASSERT(attr);
    if (attr_reject_whole_temporal(attr)) return 0;
    return attr_extract_range_F32(dst, cap, attr, 0, md_attribute_element_count(&attr->format), NULL, dst_unit, NULL);
}

size_t md_attribute_extract_f64(double dst[], size_t cap, const md_attribute_t* attr, md_unit_t dst_unit) {
    ASSERT(attr);
    if (attr_reject_whole_temporal(attr)) return 0;
    return attr_extract_range_F64(dst, cap, attr, 0, md_attribute_element_count(&attr->format), NULL, dst_unit, NULL);
}

// The slice extracts with an io to hand to a provider. Only an extraction context has one to give,
// so this stays internal and the public functions pass NULL.
static size_t attr_extract_slice_io_F32(float dst[], size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, md_unit_t dst_unit, md_attribute_io_t* io) {
    ASSERT(attr);
    if ((!slice || slice->num_idx == 0) && attr_reject_whole_temporal(attr)) return 0;
    size_t first, count;
    if (!attr_slice_window(&first, &count, attr, slice)) {
        return 0;
    }
    return attr_extract_range_F32(dst, cap, attr, first, count, slice, dst_unit, io);
}

static size_t attr_extract_slice_io_F64(double dst[], size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, md_unit_t dst_unit, md_attribute_io_t* io) {
    ASSERT(attr);
    if ((!slice || slice->num_idx == 0) && attr_reject_whole_temporal(attr)) return 0;
    size_t first, count;
    if (!attr_slice_window(&first, &count, attr, slice)) {
        return 0;
    }
    return attr_extract_range_F64(dst, cap, attr, first, count, slice, dst_unit, io);
}

size_t md_attribute_extract_slice_f32(float dst[], size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, md_unit_t dst_unit) {
    return attr_extract_slice_io_F32(dst, cap, attr, slice, dst_unit, NULL);
}

size_t md_attribute_extract_slice_f64(double dst[], size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, md_unit_t dst_unit) {
    return attr_extract_slice_io_F64(dst, cap, attr, slice, dst_unit, NULL);
}

static bool attr_path_valid(str_t path) {
    if (str_empty(path)) {
        return false;
    }
    if (path.ptr[0] == '/' || path.ptr[path.len - 1] == '/') {
        return false;
    }
    for (size_t i = 1; i < path.len; ++i) {
        if (path.ptr[i] == '/' && path.ptr[i - 1] == '/') {
            return false;
        }
    }
    return true;
}

// A trailing separator on the prefix is ignored, so "atom" and "atom/" behave identically.
static str_t attr_prefix_trim(str_t prefix) {
    while (prefix.len > 0 && prefix.ptr[prefix.len - 1] == '/') {
        prefix.len -= 1;
    }
    return prefix;
}

// Matches at a segment boundary: "atom" covers "atom" and "atom/charge", never "atomic/z".
// prefix is expected to have been trimmed.
static bool attr_path_covered_by(str_t path, str_t prefix) {
    if (str_empty(prefix)) {
        return true;
    }
    if (!str_begins_with(path, prefix)) {
        return false;
    }
    return path.len == prefix.len || path.ptr[prefix.len] == '/';
}

// First index whose name is not ordered before name. The array is sorted, so this is where
// name belongs and, for a prefix, where its run of matches begins.
static size_t attr_lower_bound(const md_attributes_t* attributes, str_t name) {
    size_t lo = 0;
    size_t hi = md_array_size(attributes->attr);
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (str_cmp_lex(attributes->attr[mid].path, name) < 0) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo;
}

static size_t attr_index_from_id(const md_attributes_t* attributes, md_attribute_id_t id) {
    if (id == MD_ATTRIBUTE_INVALID) {
        return SIZE_MAX;
    }
    for (size_t i = 0; i < md_array_size(attributes->attr); ++i) {
        if (attributes->attr[i].id == id) {
            return i;
        }
    }
    return SIZE_MAX;
}

// FRAME AXES. What qualifies as an axis is a property of the attribute alone: temporal, one index
// axis, one component, a number, and called "time". Anything else called "time" is walked past.
static bool attr_leaf_is_time(str_t path) {
    size_t loc;
    str_t leaf = attr_path_split(&loc, path) ? str_substr(path, loc + 1, SIZE_MAX) : path;
    return str_eq(leaf, STR_LIT("time"));
}

static bool attr_format_is_axis(const md_attribute_format_t* format, md_attribute_flags_t flags) {
    return (flags & MD_ATTRIBUTE_FLAG_TEMPORAL) &&
        format->rank == 1 && format->components == 1 &&
        format->type != MD_ATTRIBUTE_TYPE_NONE && format->type != MD_ATTRIBUTE_TYPE_STR;
}

// The nearest valid axis at or above the GROUP of path, never path itself: create asks this for a
// path that is not in the table yet, and an axis answers for itself before getting here. Walks
// group by group towards the root, "run/a/atom/time", "run/a/time", "run/time", "time".
static const md_attribute_t* attr_axis_above(const md_attributes_t* attributes, str_t path) {
    char buf[512];
    size_t loc;
    str_t group = attr_path_split(&loc, path) ? str_substr(path, 0, loc) : (str_t){0};

    for (;;) {
        size_t len = 0;
        if (!str_empty(group)) {
            if (group.len + 1 + 4 >= sizeof(buf)) {
                MD_LOG_ERROR("Attribute path '" STR_FMT "' is too long to search for a frame axis", STR_ARG(path));
                return NULL;
            }
            MEMCPY(buf, group.ptr, group.len);
            len = group.len;
            buf[len++] = '/';
        }
        MEMCPY(buf + len, "time", 4);
        len += 4;

        const md_attribute_t* axis = md_attributes_find(attributes, (str_t){buf, len});
        if (axis && attr_format_is_axis(&axis->format, axis->flags)) {
            return axis;
        }
        if (str_empty(group)) {
            return NULL;
        }
        group = attr_path_split(&loc, group) ? str_substr(group, 0, loc) : (str_t){0};
    }
}

// The storage tag is STATED by the producer rather than derived from the other fields, and stating
// it is what makes this check possible: every alternative reading of an attribute - who owns the
// bytes, who computes them, whose name this is - is already recorded elsewhere in the struct, so
// the tag and those fields can be held against each other. A derived tag could never disagree and
// could never catch anything either.
//
// What each tag claims, and what would contradict it:
//   RESIDENT  owns its bytes, under its own name, nothing computes it.
//   VIRTUAL   computed under its own name, so there are no resident bytes to own.
//   ALIAS     a second NAME for somebody else's datum, so root names that owner and never itself,
//             and it owns no provider state (user_data_size is the ownership marker, and an alias
//             zeroes it when it inherits virt).
// Note an ALIAS is deliberately unconstrained in data/provider: it inherits whichever its target
// had, which is exactly why read paths branch on the provider and not on this tag.
static bool attr_storage_consistent(const md_attribute_t* attr) {
    ASSERT(attr);
    switch (attr->storage) {
    case MD_ATTRIBUTE_STORAGE_RESIDENT:
        return attr->root == attr->id && attr->virt.provider == NULL;
    case MD_ATTRIBUTE_STORAGE_VIRTUAL:
        return attr->root == attr->id && attr->virt.provider != NULL && attr->data == NULL;
    case MD_ATTRIBUTE_STORAGE_ALIAS:
        return attr->root != attr->id && attr->virt.user_data_size == 0;
    default:
        return false;
    }
}

// Everything a new path has to clear before anything is allocated for it, shared by the two
// producers so they cannot drift on what "this path is available" means. Hands back where the
// entry belongs in the sorted array and the id it will have.
static bool attr_reserve_slot(md_attributes_t* attributes, str_t path, size_t* out_idx, md_attribute_id_t* out_id) {
    ASSERT(attributes);

    if (!attributes->alloc) {
        MD_LOG_ERROR("Attribute table allocator not set");
        return false;
    }
    if (!attr_path_valid(path)) {
        MD_LOG_ERROR("Invalid attribute path '" STR_FMT "': expected non empty segments separated by '/'", STR_ARG(path));
        return false;
    }

    size_t idx = attr_lower_bound(attributes, path);
    if (idx < md_array_size(attributes->attr) && str_eq(attributes->attr[idx].path, path)) {
        MD_LOG_ERROR("Attribute '" STR_FMT "' already exists", STR_ARG(path));
        return false;
    }

    md_attribute_id_t id = md_attributes_id_from_path(path);
    if (attr_index_from_id(attributes, id) != SIZE_MAX) {
        MD_LOG_ERROR("Hash collision for attribute '" STR_FMT "'", STR_ARG(path));
        return false;
    }

    *out_idx = idx;
    *out_id  = id;
    return true;
}

// The one place an attribute enters the table, which is why the tag is checked here: both producers
// pass through it and neither can install a struct whose storage disagrees with the rest of it.
// idx comes from attr_reserve_slot, so opening the hole keeps the array sorted by path.
static void attr_insert_at(md_attributes_t* attributes, size_t idx, const md_attribute_t* attr) {
    ASSERT(attr_storage_consistent(attr));

    md_attribute_t empty = {0};
    md_array_push(attributes->attr, empty, attributes->alloc);

    size_t count = md_array_size(attributes->attr);
    if (idx + 1 < count) {
        MEMMOVE(attributes->attr + idx + 1, attributes->attr + idx, (count - 1 - idx) * sizeof(md_attribute_t));
    }
    attributes->attr[idx] = *attr;
}

// Everything one attribute owns. Removing one and tearing the whole table down have to agree about
// this, and when they were two copies of the list they were one edit away from disagreeing - which
// is a leak in one path or a double free in the other, neither visible until it is.
static void attr_release(md_attribute_t* attr, md_allocator_i* alloc) {
    ASSERT(attr);
    ASSERT(alloc);
    // A tag disagreeing with the rest of the struct is a leak on one branch and a double free on
    // the other, and neither shows up where it was caused.
    ASSERT(attr_storage_consistent(attr));
    str_free(attr->path, alloc);
    // Both are optional and an absent one is a zeroed str_t, which is not something to hand to an
    // allocator.
    if (attr->label.ptr)       str_free(attr->label, alloc);
    if (attr->description.ptr) str_free(attr->description, alloc);
    // Only an owner releases storage; an alias borrows both the buffer and the provider state.
    if (attr->storage == MD_ATTRIBUTE_STORAGE_RESIDENT && attr->data) {
        md_free(alloc, attr->data, md_attribute_byte_size(&attr->format));
    }
    // user_data_size is 0 for a borrowed pointer, so this only ever frees memory this table itself
    // handed out through md_attributes_alloc_user_data.
    if (attr->virt.user_data && attr->virt.user_data_size) {
        md_free(alloc, attr->virt.user_data, attr->virt.user_data_size);
    }
}

void md_attributes_free(md_attributes_t* attributes) {
    ASSERT(attributes);
    md_allocator_i* alloc = attributes->alloc;
    if (alloc) {
        for (size_t i = 0; i < md_array_size(attributes->attr); ++i) {
            attr_release(attributes->attr + i, alloc);
        }
        md_array_free(attributes->attr, alloc);
        // The pool is append only for the table's whole life, so this is the one place it goes.
        md_array_free(attributes->str_data,   alloc);
        md_array_free(attributes->str_offset, alloc);
        md_array_free(attributes->str_index,  alloc);
    }
    MEMSET(attributes, 0, sizeof(md_attributes_t));
}

size_t md_attributes_count(const md_attributes_t* attributes) {
    ASSERT(attributes);
    return md_array_size(attributes->attr);
}

void* md_attributes_alloc_user_data(md_attributes_t* attributes, size_t size) {
    ASSERT(attributes);
    if (size == 0) {
        return NULL;
    }
    if (!attributes->alloc) {
        MD_LOG_ERROR("Attribute table allocator not set");
        return NULL;
    }
    return md_alloc(attributes->alloc, size);
}

md_attribute_id_t md_attributes_id_from_path(str_t path) {
    md_attribute_id_t id = (md_attribute_id_t)md_hash64_str(path, 0);
    if (id == MD_ATTRIBUTE_INVALID) {
        // Zero is the invalid value, so the one path that hashes to it borrows the next id.
        id = 1;
    }
    return id;
}

// ---------------------------------------------------------------------------
// The string pool
// ---------------------------------------------------------------------------
// Append only, interning, and freed with the table. See the STRINGS note in md_system.h for why the
// element is a handle rather than the text.

// Entry 0 is the empty string, so a zeroed handle reads as "" instead of as garbage. Called before
// every intern rather than at table creation, because a table is zero initialised by its owner and
// there is no init hook to hang this on.
static bool attr_str_pool_init(md_attributes_t* attributes) {
    if (md_array_size(attributes->str_offset) > 0) {
        return true;
    }
    md_array_push(attributes->str_data,   '\0', attributes->alloc);
    md_array_push(attributes->str_offset, 0u,   attributes->alloc);
    md_array_push(attributes->str_offset, 1u,   attributes->alloc);
    return md_array_size(attributes->str_offset) == 2;
}

static str_t attr_str_pool_get(const md_attributes_t* attributes, uint32_t handle) {
    const size_t count = md_array_size(attributes->str_offset);
    if (count < 2 || (size_t)handle + 1 >= count) {
        return (str_t){0};
    }
    const uint32_t beg = attributes->str_offset[handle];
    const uint32_t end = attributes->str_offset[handle + 1];
    // The stored NUL is not part of the string, but it IS there, so str_ptr can be handed to a C
    // API without a copy.
    return (str_t){ attributes->str_data + beg, (size_t)(end - beg - 1) };
}

static void attr_str_index_insert(md_attributes_t* attributes, uint64_t hash, uint32_t handle) {
    const size_t mask = md_array_size(attributes->str_index) - 1;
    size_t slot = (size_t)hash & mask;
    while (attributes->str_index[slot] != 0) {
        slot = (slot + 1) & mask;
    }
    attributes->str_index[slot] = handle + 1;
}

// Grown at 2/3 load. Rehashing walks the pool rather than storing the hashes: the entries are right
// there and this happens O(log n) times over a table's life.
static bool attr_str_index_reserve(md_attributes_t* attributes, size_t needed) {
    const size_t cap = md_array_size(attributes->str_index);
    if (cap != 0 && needed * 3 <= cap * 2) {
        return true;
    }
    size_t new_cap = cap ? cap * 2 : 64;
    while (needed * 3 > new_cap * 2) {
        new_cap *= 2;
    }

    md_array(uint32_t) old_index = attributes->str_index;
    attributes->str_index = 0;
    md_array_resize(attributes->str_index, new_cap, attributes->alloc);
    if (md_array_size(attributes->str_index) != new_cap) {
        attributes->str_index = old_index;
        return false;
    }
    MEMSET(attributes->str_index, 0, sizeof(uint32_t) * new_cap);

    const size_t num_entries = md_array_size(attributes->str_offset) - 1;
    for (uint32_t h = 0; h < (uint32_t)num_entries; ++h) {
        const str_t s = attr_str_pool_get(attributes, h);
        attr_str_index_insert(attributes, md_hash64_str(s, 0), h);
    }
    md_array_free(old_index, attributes->alloc);
    return true;
}

// The one way text enters the table. Returns the handle; 0 (the empty string) is a valid answer and
// also what a failure degrades to, which keeps a caller from having to branch on an error it cannot
// do anything about.
static uint32_t attr_str_intern(md_attributes_t* attributes, str_t str) {
    if (!attr_str_pool_init(attributes) || str_empty(str)) {
        return 0;
    }

    const uint64_t hash = md_hash64_str(str, 0);
    if (md_array_size(attributes->str_index) > 0) {
        const size_t mask = md_array_size(attributes->str_index) - 1;
        size_t slot = (size_t)hash & mask;
        while (attributes->str_index[slot] != 0) {
            const uint32_t handle = attributes->str_index[slot] - 1;
            if (str_eq(attr_str_pool_get(attributes, handle), str)) {
                return handle;
            }
            slot = (slot + 1) & mask;
        }
    }

    const uint32_t handle = (uint32_t)(md_array_size(attributes->str_offset) - 1);
    md_array_push_array(attributes->str_data, str.ptr, str.len, attributes->alloc);
    md_array_push(attributes->str_data, '\0', attributes->alloc);
    md_array_push(attributes->str_offset, (uint32_t)md_array_size(attributes->str_data), attributes->alloc);

    if (!attr_str_index_reserve(attributes, (size_t)handle + 1)) {
        return handle;  // interned and readable; only the dedup lookup is degraded
    }
    attr_str_index_insert(attributes, hash, handle);
    return handle;
}

str_t md_attribute_str(const md_attributes_t* attributes, const md_attribute_t* attr, size_t index) {
    ASSERT(attributes);
    ASSERT(attr);
    if (attr->format.type != MD_ATTRIBUTE_TYPE_STR || !attr->data) {
        MD_LOG_ERROR("Attribute '" STR_FMT "' is not a resident string attribute", STR_ARG(attr->path));
        return (str_t){0};
    }
    if (index >= md_attribute_element_count(&attr->format)) {
        return (str_t){0};
    }
    return attr_str_pool_get(attributes, ((const uint32_t*)attr->data)[index]);
}

size_t md_attribute_extract_str(str_t dst[], size_t cap, const md_attributes_t* attributes, const md_attribute_t* attr) {
    ASSERT(attributes);
    ASSERT(attr);
    if (!dst) {
        return 0;
    }
    if (attr->format.type != MD_ATTRIBUTE_TYPE_STR || !attr->data) {
        MD_LOG_ERROR("Attribute '" STR_FMT "' is not a resident string attribute", STR_ARG(attr->path));
        return 0;
    }
    const size_t count = md_attribute_element_count(&attr->format);
    if (count > cap) {
        MD_LOG_ERROR("Attribute '" STR_FMT "' needs %zu values, %zu supplied", STR_ARG(attr->path), count, cap);
        return 0;
    }
    const uint32_t* handles = (const uint32_t*)attr->data;
    for (size_t i = 0; i < count; ++i) {
        dst[i] = attr_str_pool_get(attributes, handles[i]);
    }
    return count;
}

md_attribute_id_t md_attributes_create(md_attributes_t* attributes, const md_attribute_desc_t* desc) {
    ASSERT(attributes);

    if (!desc) {
        MD_LOG_ERROR("Attribute descriptor is NULL");
        return MD_ATTRIBUTE_INVALID;
    }

    md_attribute_format_t format = desc->format;
    const str_t path = desc->path;

    // The path is settled before the format is looked at: it is the cheaper rejection, and "that
    // path is taken" is a more actionable message than a format complaint about a create which was
    // never going to land anyway.
    size_t idx;
    md_attribute_id_t id;
    if (!attr_reserve_slot(attributes, path, &idx, &id)) {
        return MD_ATTRIBUTE_INVALID;
    }

    if (md_attribute_type_size(format.type) == 0) {
        MD_LOG_ERROR("Invalid type for attribute '" STR_FMT "'", STR_ARG(path));
        return MD_ATTRIBUTE_INVALID;
    }
    // Zero components is a producer who declared the extents and forgot to say how wide a value
    // is. It is not defaulted to 1: that would make the mistake legal and silent, and a wrong
    // components is indistinguishable from a right one everywhere downstream.
    if (format.components == 0) {
        MD_LOG_ERROR("Attribute '" STR_FMT "' declares 0 components; one value is at least 1 component wide", STR_ARG(path));
        return MD_ATTRIBUTE_INVALID;
    }
    if (format.rank > MD_ATTRIBUTE_MAX_RANK) {
        MD_LOG_ERROR("Rank %u exceeds MD_ATTRIBUTE_MAX_RANK for attribute '" STR_FMT "'", format.rank, STR_ARG(path));
        return MD_ATTRIBUTE_INVALID;
    }
    for (uint32_t i = 0; i < format.rank; ++i) {
        if (format.shape[i] == 0) {
            MD_LOG_ERROR("Zero extent in axis %u of attribute '" STR_FMT "'", i, STR_ARG(path));
            return MD_ATTRIBUTE_INVALID;
        }
    }
    // An extent past the declared rank is not harmless padding, it is a producer who wrote the
    // old spelling: a trailing component axis left in shape with rank not narrowed to match.
    // Nothing reads those slots, so accepting it would store a format which reads one way in a
    // debugger and behaves another.
    for (uint32_t i = format.rank; i < MD_ATTRIBUTE_MAX_RANK; ++i) {
        if (format.shape[i] != 0) {
            MD_LOG_ERROR("Attribute '" STR_FMT "' sets extent %u on axis %u, beyond its rank of %u",
                STR_ARG(path), format.shape[i], i, format.rank);
            return MD_ATTRIBUTE_INVALID;
        }
    }

    size_t required = md_attribute_byte_size(&format);

    // A string attribute is published as the TEXT and stored as handles, so what the caller hands
    // over and what the table keeps are different sizes. The guard is checked against the input,
    // where the caller's mistake would be.
    const bool is_str = (format.type == MD_ATTRIBUTE_TYPE_STR);
    const size_t input_size = is_str ? md_attribute_element_count(&format) * sizeof(str_t) : required;

    if (is_str && desc->virt) {
        // A provider writes the STORED type, and the stored type here is a handle into a pool the
        // provider has no way to intern into.
        MD_LOG_ERROR("Attribute '" STR_FMT "' is a string and cannot be virtual", STR_ARG(path));
        return MD_ATTRIBUTE_INVALID;
    }

    if (desc->virt) {
        if (!desc->virt->provider) {
            MD_LOG_ERROR("Attribute '" STR_FMT "' declares virt with no provider", STR_ARG(path));
            return MD_ATTRIBUTE_INVALID;
        }
        if (desc->data || desc->byte_size != 0) {
            MD_LOG_ERROR("Attribute '" STR_FMT "' declares virt and resident data at the same time", STR_ARG(path));
            return MD_ATTRIBUTE_INVALID;
        }
    } else if (desc->data) {
        if (desc->byte_size != input_size) {
            MD_LOG_ERROR("Attribute '" STR_FMT "' declares %zu bytes but %zu were supplied", STR_ARG(path), input_size, desc->byte_size);
            return MD_ATTRIBUTE_INVALID;
        }
    } else if (desc->byte_size != 0) {
        // A size without a pointer is a caller who meant to pass one.
        MD_LOG_ERROR("Attribute '" STR_FMT "' supplied %zu bytes with no data pointer", STR_ARG(path), desc->byte_size);
        return MD_ATTRIBUTE_INVALID;
    }

    // TEMPORAL is a claim about the outermost axis, so it is checked rather than believed. This is
    // the whole reason for tagging instead of inferring: a shape that merely looks frame sized is a
    // coincidence, while a tag that disagrees with its axis is a bug, and catching it here beats
    // finding it partway through an extract. An axis is its own axis and has nothing to agree with.
    if (desc->flags & MD_ATTRIBUTE_FLAG_TEMPORAL) {
        if (desc->format.rank == 0) {
            MD_LOG_ERROR("Attribute '" STR_FMT "' is temporal but has no index axes", STR_ARG(path));
            return MD_ATTRIBUTE_INVALID;
        }
        if (!(attr_leaf_is_time(path) && attr_format_is_axis(&format, desc->flags))) {
            const md_attribute_t* axis = attr_axis_above(attributes, path);
            if (!axis) {
                MD_LOG_ERROR("Attribute '" STR_FMT "' is temporal but no frame axis ('time') exists at or above it", STR_ARG(path));
                return MD_ATTRIBUTE_INVALID;
            }
            if (desc->format.shape[0] != axis->format.shape[0]) {
                MD_LOG_ERROR("Attribute '" STR_FMT "' is temporal with an outermost extent of %u, but its axis '" STR_FMT "' has %u frames",
                    STR_ARG(path), desc->format.shape[0], STR_ARG(axis->path), axis->format.shape[0]);
                return MD_ATTRIBUTE_INVALID;
            }
        }
    }

    md_allocator_i* alloc = attributes->alloc;

    void* storage = NULL;
    if (!desc->virt) {
        storage = md_alloc(alloc, required);
        if (!storage) {
            MD_LOG_ERROR("Failed to allocate %zu bytes for attribute '" STR_FMT "'", required, STR_ARG(path));
            return MD_ATTRIBUTE_INVALID;
        }
        if (desc->data && is_str) {
            // Interned one at a time. A zeroed handle is the empty string, so a failure to intern
            // degrades to "" rather than to a dangling index.
            const str_t* values = (const str_t*)desc->data;
            uint32_t* handles = (uint32_t*)storage;
            const size_t count = md_attribute_element_count(&format);
            for (size_t i = 0; i < count; ++i) {
                handles[i] = attr_str_intern(attributes, values[i]);
            }
        } else if (desc->data) {
            MEMCPY(storage, desc->data, required);
        } else {
            MEMSET(storage, 0, required);
        }
    }

    str_t stored_path = str_copy(path, alloc);
    if (str_empty(stored_path)) {
        if (storage) {
            md_free(alloc, storage, required);
        }
        MD_LOG_ERROR("Failed to copy attribute path '" STR_FMT "'", STR_ARG(path));
        return MD_ATTRIBUTE_INVALID;
    }

    // Presentation only, so an empty one is a valid state and a failed copy is not worth failing
    // the create over: the consumer's fallback for "no label" is the same either way.
    str_t stored_label = str_empty(desc->label) ? (str_t){0} : str_copy(desc->label, alloc);
    str_t stored_desc  = str_empty(desc->description) ? (str_t){0} : str_copy(desc->description, alloc);

    const md_attribute_t entry = {
        .id          = id,
        .path        = stored_path,
        .label       = stored_label,
        .description = stored_desc,
        .format      = format,
        .unit        = desc->unit,
        .flags       = desc->flags,
        .version     = ++attributes->version_counter,
        .storage     = desc->virt ? MD_ATTRIBUTE_STORAGE_VIRTUAL : MD_ATTRIBUTE_STORAGE_RESIDENT,
        .data        = storage,
        .virt        = desc->virt ? *desc->virt : (md_attribute_virtual_t){0},
        .root        = id,   // it owns its own storage; an alias is what points elsewhere
    };
    attr_insert_at(attributes, idx, &entry);

    return id;
}

md_attribute_id_t md_attributes_alias(md_attributes_t* attributes, md_attribute_id_t target, str_t path, str_t label, str_t description) {
    ASSERT(attributes);

    size_t idx;
    md_attribute_id_t id;
    if (!attr_reserve_slot(attributes, path, &idx, &id)) {
        return MD_ATTRIBUTE_INVALID;
    }

    size_t target_idx = attr_index_from_id(attributes, target);
    if (target_idx == SIZE_MAX) {
        MD_LOG_ERROR("Cannot alias '" STR_FMT "': no such target attribute", STR_ARG(path));
        return MD_ATTRIBUTE_INVALID;
    }

    // Everything inherited is copied out BEFORE the array grows: md_array_push may reallocate, and
    // a pointer into the old block is exactly the stale pointer the header warns callers about.
    // Aliasing an alias flattens here - root already names the owner - so a chain is never walked.
    const md_attribute_t*        tgt      = attributes->attr + target_idx;
    const md_attribute_format_t  format   = tgt->format;
    const md_unit_t              unit     = tgt->unit;
    void* const                  data     = tgt->data;
    const md_attribute_flags_t   flags    = tgt->flags;
    const md_attribute_id_t      root     = tgt->root;
    const uint64_t               version  = tgt->version;
    md_attribute_virtual_t       virt     = tgt->virt;

    // The target keeps ownership of its provider's private state. Zero the size so teardown of the
    // alias never reaches it - the size IS the ownership marker, as md_attribute_virtual_t says.
    virt.user_data_size = 0;

    md_allocator_i* alloc = attributes->alloc;

    str_t stored_path = str_copy(path, alloc);
    if (str_empty(stored_path)) {
        MD_LOG_ERROR("Failed to copy attribute path '" STR_FMT "'", STR_ARG(path));
        return MD_ATTRIBUTE_INVALID;
    }
    str_t stored_label = str_empty(label)       ? (str_t){0} : str_copy(label, alloc);
    str_t stored_desc  = str_empty(description) ? (str_t){0} : str_copy(description, alloc);

    const md_attribute_t entry = {
        .id          = id,
        .path        = stored_path,
        .label       = stored_label,
        .description = stored_desc,
        .format      = format,
        .unit        = unit,
        .flags       = flags,
        // Inherited, not stamped: naming a datum a second time does not change its contents, and a
        // fresh counter value here would read to a consumer as "this changed" on every reload that
        // re-established the alias.
        .version     = version,
        .storage     = MD_ATTRIBUTE_STORAGE_ALIAS,
        .data        = data,
        .virt        = virt,
        .root        = root,
    };
    attr_insert_at(attributes, idx, &entry);

    return id;
}

bool md_attribute_same_data(const md_attribute_t* a, const md_attribute_t* b) {
    ASSERT(a);
    ASSERT(b);
    return a->root == b->root;
}

static void attr_remove_at(md_attributes_t* attributes, size_t idx) {
    md_allocator_i* alloc = attributes->alloc;
    ASSERT(alloc);

    attr_release(attributes->attr + idx, alloc);

    size_t count = md_array_size(attributes->attr);
    if (idx + 1 < count) {
        MEMMOVE(attributes->attr + idx, attributes->attr + idx + 1, (count - 1 - idx) * sizeof(md_attribute_t));
    }
    md_array_pop(attributes->attr);
}

md_attribute_id_t md_attributes_publish_atom_column(md_attributes_t* attributes, str_t path, md_unit_t unit, uint32_t components, const float values[], size_t count) {
    ASSERT(attributes);

    if (count == 0 || components == 0 || !values) {
        return MD_ATTRIBUTE_INVALID;
    }

    // Uniformity is a property of the whole VALUE, not of each component: a velocity column where
    // every atom moves the same way is as uninformative as a constant occupancy, and one where only
    // the x components happen to agree is not.
    const size_t num_elements = count * (size_t)components;

    // NAN is how a loader marks "this atom has no value", for a column its format defines but this
    // particular file leaves blank per atom. A column that is entirely absent is exactly as
    // uninformative as a constant one, so it is skipped for the same reason.
    //
    // BOTH TESTS BELOW ARE ON BITS, NOT ON FLOATS, and that is not a micro optimisation. This
    // library is built with /fp:fast on MSVC and -ffast-math on Clang, which license the compiler
    // to assume no operand is ever NaN - and it takes the licence: 'x != x' folds to false,
    // isnan() folds to 0, and even a comparison AGAINST a NaN reports equal. Written as float
    // comparisons this loop dropped any column whose FIRST atom was blank, which is a partly
    // filled column silently becoming no column at all.
    //
    // Bit equality also decides uniformity, which is stricter than '==' in exactly one place:
    // +0.0 and -0.0 read as different values. A column mixing the two therefore publishes. That
    // is the harmless direction, and no format in the tree distinguishes them anyway.
    bool uniform = true;
    bool all_absent = true;
    for (size_t i = 0; i < num_elements; ++i) {
        const size_t c = i % components;
        const uint32_t bits_i = attr_f32_bits(values[i]);
        const uint32_t bits_0 = attr_f32_bits(values[c]);
        if (!attr_f32_bits_absent(bits_i)) {
            all_absent = false;
        }
        // Covers both questions at once: a present value against an absent one differs in bits,
        // and so do two different present values.
        if (bits_i != bits_0) {
            uniform = false;
        }
    }
    if (uniform || all_absent) {
        return MD_ATTRIBUTE_INVALID;
    }

    // One scalar per atom: the atom axis is the only index axis, and a value is one component wide.
    // components is stated rather than left to a default; see the ATTRIBUTES note above.
    const md_attribute_desc_t desc = {
        .path   = path,
        .format = {
            .type       = MD_ATTRIBUTE_TYPE_F32,
            .components = components,
            .rank       = 1,
            .shape      = {(uint32_t)count},
        },
        .unit      = unit,
        .data      = values,
        .byte_size = num_elements * sizeof(float),
    };
    return md_attributes_create(attributes, &desc);
}

uint64_t md_attributes_version(const md_attributes_t* attributes, md_attribute_id_t id) {
    ASSERT(attributes);
    size_t idx = attr_index_from_id(attributes, id);
    return idx == SIZE_MAX ? 0 : attributes->attr[idx].version;
}

uint64_t md_attributes_touch(md_attributes_t* attributes, md_attribute_id_t id) {
    ASSERT(attributes);
    size_t idx = attr_index_from_id(attributes, id);
    if (idx == SIZE_MAX) {
        return 0;
    }

    // A version describes a DATUM and an alias is a second name for one, so a touch has to reach
    // every name of it. Bumping only the name the producer happened to call is the invalidation
    // hole aliases were always going to open: the whole point of aliasing is that a consumer reads
    // through the neutral path, and that consumer would have cached against a version which never
    // moved again. One counter value for all of them, so they compare equal as well as fresh.
    const md_attribute_id_t root = attributes->attr[idx].root;
    const uint64_t version = ++attributes->version_counter;
    for (size_t i = 0; i < md_array_size(attributes->attr); ++i) {
        if (attributes->attr[i].root == root) {
            attributes->attr[i].version = version;
        }
    }
    return version;
}

md_attribute_id_t md_attributes_replace(md_attributes_t* attributes, const md_attribute_desc_t* desc) {
    ASSERT(attributes);
    ASSERT(desc);

    const md_attribute_t* existing = md_attributes_find(attributes, desc->path);
    if (existing) {
        md_attributes_remove(attributes, existing->id);
    }
    return md_attributes_create(attributes, desc);
}

bool md_attributes_remove(md_attributes_t* attributes, md_attribute_id_t id) {
    ASSERT(attributes);

    if (attr_index_from_id(attributes, id) == SIZE_MAX) {
        return false;
    }

    // Aliases read this attribute's storage, so they cannot outlive it. They go first, one at a
    // time because every removal shifts the array; aliases are flattened at creation, so nothing
    // aliases an alias and one pass over the survivors always terminates.
    for (;;) {
        size_t alias_idx = SIZE_MAX;
        for (size_t i = 0; i < md_array_size(attributes->attr); ++i) {
            const md_attribute_t* a = attributes->attr + i;
            if (a->storage == MD_ATTRIBUTE_STORAGE_ALIAS && a->root == id) {
                alias_idx = i;
                break;
            }
        }
        if (alias_idx == SIZE_MAX) {
            break;
        }
        attr_remove_at(attributes, alias_idx);
    }

    // Re-find: removing the aliases shifted everything after them.
    size_t idx = attr_index_from_id(attributes, id);
    if (idx == SIZE_MAX) {
        return false;
    }
    attr_remove_at(attributes, idx);

    return true;
}

size_t md_attributes_remove_prefix(md_attributes_t* attributes, str_t prefix) {
    ASSERT(attributes);

    if (str_empty(attr_prefix_trim(prefix))) {
        MD_LOG_ERROR("Refusing to remove attributes under an empty prefix");
        return 0;
    }

    // In batches, because every removal shifts the array and may take aliases elsewhere with it -
    // which is also why the count is taken from the table rather than from the removals.
    const size_t before = md_array_size(attributes->attr);
    md_attribute_id_t ids[64];
    for (;;) {
        const size_t n = md_attributes_query(ids, ARRAY_SIZE(ids), attributes, prefix);
        if (n == 0) {
            break;
        }
        for (size_t i = 0; i < MIN(n, ARRAY_SIZE(ids)); ++i) {
            md_attributes_remove(attributes, ids[i]);
        }
    }
    return before - md_array_size(attributes->attr);
}

const md_attribute_t* md_attributes_get(const md_attributes_t* attributes, md_attribute_id_t id) {
    ASSERT(attributes);
    size_t idx = attr_index_from_id(attributes, id);
    return idx == SIZE_MAX ? NULL : attributes->attr + idx;
}

const md_attribute_t* md_attributes_find(const md_attributes_t* attributes, str_t path) {
    ASSERT(attributes);
    size_t idx = attr_lower_bound(attributes, path);
    if (idx < md_array_size(attributes->attr) && str_eq(attributes->attr[idx].path, path)) {
        return attributes->attr + idx;
    }
    return NULL;
}

void* md_attributes_data(md_attributes_t* attributes, md_attribute_id_t id, md_attribute_type_t expected_type) {
    ASSERT(attributes);
    size_t idx = attr_index_from_id(attributes, id);
    if (idx == SIZE_MAX) {
        return NULL;
    }
    md_attribute_t* attr = attributes->attr + idx;
    if (attr->storage == MD_ATTRIBUTE_STORAGE_ALIAS) {
        // Writing through a second name would be writing to somebody else's attribute behind its
        // back. Fill the owner in and every name sees it.
        MD_LOG_ERROR("Attribute '" STR_FMT "' is an alias; fill in the attribute which owns the storage", STR_ARG(attr->path));
        return NULL;
    }
    if (attr->storage != MD_ATTRIBUTE_STORAGE_RESIDENT) {
        MD_LOG_ERROR("Attribute '" STR_FMT "' is virtual and has no resident storage", STR_ARG(attr->path));
        return NULL;
    }
    if (attr->format.type == MD_ATTRIBUTE_TYPE_STR) {
        // The storage is handles into the pool. Handing it out writable invites a producer to
        // fabricate one, which is how a table ends up pointing at text it does not own. Publish the
        // strings through md_attributes_create instead - that is the only way text gets in.
        MD_LOG_ERROR("Attribute '" STR_FMT "' is a string; publish its values rather than writing handles", STR_ARG(attr->path));
        return NULL;
    }
    if (attr->format.type != expected_type) {
        MD_LOG_ERROR("Type mismatch for attribute '" STR_FMT "'", STR_ARG(attr->path));
        return NULL;
    }
    return attr->data;
}

size_t md_attributes_query(md_attribute_id_t out_ids[], size_t cap, const md_attributes_t* attributes, str_t prefix) {
    return md_attributes_query_flags(out_ids, cap, attributes, prefix, MD_ATTRIBUTE_FLAG_NONE, MD_ATTRIBUTE_FLAG_NONE);
}

size_t md_attributes_query_flags(md_attribute_id_t out_ids[], size_t cap, const md_attributes_t* attributes, str_t prefix,
                                 md_attribute_flags_t mask, md_attribute_flags_t value) {
    ASSERT(attributes);

    str_t base = attr_prefix_trim(prefix);
    size_t count = 0;

    // Sorted by name, so the run starts at the lower bound of the prefix. A path may begin
    // with the prefix and continue with a character below '/' (say "atom-x" under "atom"),
    // which sorts inside the run without being covered by it, hence the two conditions.
    for (size_t i = str_empty(base) ? 0 : attr_lower_bound(attributes, base); i < md_array_size(attributes->attr); ++i) {
        str_t name = attributes->attr[i].path;
        if (!str_empty(base) && !str_begins_with(name, base)) {
            break;
        }
        if (!attr_path_covered_by(name, base)) {
            continue;
        }
        if ((attributes->attr[i].flags & mask) != value) {
            continue;
        }
        if (out_ids && count < cap) {
            out_ids[count] = attributes->attr[i].id;
        }
        count += 1;
    }

    return count;
}

size_t md_attributes_query_children(str_t out_names[], size_t cap, const md_attributes_t* attributes, str_t prefix) {
    ASSERT(attributes);

    str_t base = attr_prefix_trim(prefix);
    size_t count = 0;
    str_t prev = {0};

    for (size_t i = str_empty(base) ? 0 : attr_lower_bound(attributes, base); i < md_array_size(attributes->attr); ++i) {
        str_t name = attributes->attr[i].path;
        if (!str_empty(base) && !str_begins_with(name, base)) {
            break;
        }
        if (!attr_path_covered_by(name, base)) {
            continue;
        }
        // The prefix itself is an attribute, not a child of one.
        if (name.len == base.len) {
            continue;
        }

        size_t offset = str_empty(base) ? 0 : base.len + 1;
        str_t rest = str_substr(name, offset, SIZE_MAX);
        size_t loc;
        str_t child = str_find_char(&loc, rest, '/') ? str_substr(rest, 0, loc) : rest;

        // Names are sorted, so every path under one child is contiguous and a duplicate can
        // only ever be the one just emitted.
        if (!str_empty(prev) && str_eq(child, prev)) {
            continue;
        }
        prev = child;

        if (out_names && count < cap) {
            out_names[count] = child;
        }
        count += 1;
    }

    return count;
}

const md_attribute_t* md_attributes_axis(const md_attributes_t* attributes, const md_attribute_t* attr) {
    ASSERT(attributes);
    ASSERT(attr);

    // An alias is a second name for a datum, and the datum's axis is decided by where its owner
    // lives. Searching from the alias' own path would pair it with whatever "time" happens to sit
    // above the new name.
    const md_attribute_t* owner = attr;
    if (attr->root != attr->id) {
        owner = md_attributes_get(attributes, attr->root);
        if (!owner) {
            return NULL;
        }
    }
    if (!(owner->flags & MD_ATTRIBUTE_FLAG_TEMPORAL)) {
        return NULL;
    }
    if (attr_leaf_is_time(owner->path) && attr_format_is_axis(&owner->format, owner->flags)) {
        return owner;
    }
    return attr_axis_above(attributes, owner->path);
}

static inline double attr_abs(double x) {
    return x < 0.0 ? -x : x;
}

// One coordinate, read through the extract so a computed axis works exactly like a resident one.
// unit is md_unit_none() for "as stored".
static bool attr_axis_value(double* out, const md_attribute_t* axis, size_t i, md_unit_t unit) {
    const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)i);
    return md_attribute_extract_slice_f64(out, 1, axis, &slice, unit) == 1;
}

bool md_attribute_axis_map(size_t* out_index, const md_attribute_t* src_axis, size_t src_index, const md_attribute_t* dst_axis) {
    ASSERT(out_index);
    ASSERT(src_axis);
    ASSERT(dst_axis);

    if (!attr_format_is_axis(&src_axis->format, src_axis->flags) || !attr_format_is_axis(&dst_axis->format, dst_axis->flags)) {
        MD_LOG_ERROR("Mapping between '" STR_FMT "' and '" STR_FMT "': both have to be frame axes", STR_ARG(src_axis->path), STR_ARG(dst_axis->path));
        return false;
    }

    const size_t src_count = src_axis->format.shape[0];
    const size_t dst_count = dst_axis->format.shape[0];
    if (src_index >= src_count) {
        return false;
    }

    if (md_attribute_same_data(src_axis, dst_axis)) {
        *out_index = src_index;
        return true;
    }

    // Ordinals only meet ordinals. With a unit on both sides the extract converts, and refuses
    // outright when the dimensions differ, which is the same answer for the same reason.
    if (md_unit_is_none(src_axis->unit) != md_unit_is_none(dst_axis->unit)) {
        return false;
    }

    double t;
    if (!attr_axis_value(&t, src_axis, src_index, dst_axis->unit)) {
        return false;
    }

    // First coordinate not below t; the match is that one or the one before it.
    size_t lo = 0;
    size_t hi = dst_count;
    while (lo < hi) {
        const size_t mid = lo + (hi - lo) / 2;
        double v;
        if (!attr_axis_value(&v, dst_axis, mid, md_unit_none())) {
            return false;
        }
        if (v < t) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }

    size_t best = SIZE_MAX;
    double best_dist = 0.0;
    double best_val  = 0.0;
    for (size_t c = (lo > 0 ? lo - 1 : 0); c <= lo && c < dst_count; ++c) {
        double v;
        if (!attr_axis_value(&v, dst_axis, c, md_unit_none())) {
            return false;
        }
        const double d = attr_abs(v - t);
        if (best == SIZE_MAX || d < best_dist) {
            best = c;
            best_dist = d;
            best_val = v;
        }
    }
    if (best == SIZE_MAX) {
        return false;
    }

    // Float precision of the value itself: an XTC stores its times as float, so 12345.6 ps comes
    // back a few ulps off, and a fixed absolute margin either rejects that or accepts neighbours
    // on a finely sampled axis. A thousandth of the local spacing covers the other side.
    double tol = 1.0e-6 * attr_abs(t);
    double spacing = 0.0;
    if (best + 1 < dst_count) {
        double v;
        if (attr_axis_value(&v, dst_axis, best + 1, md_unit_none())) {
            spacing = attr_abs(v - best_val);
        }
    }
    if (best > 0) {
        double v;
        if (attr_axis_value(&v, dst_axis, best - 1, md_unit_none())) {
            const double s = attr_abs(best_val - v);
            if (s > 0.0 && (spacing == 0.0 || s < spacing)) {
                spacing = s;
            }
        }
    }
    tol = MAX(tol, 1.0e-3 * spacing);

    if (best_dist > tol) {
        return false;
    }
    *out_index = best;
    return true;
}

// ### IO ###

#define ATTR_IO_MAX_FILES 16

// A small cache of open files, keyed by path. A handful is enough: an extraction context reads one
// run, and a run is one file per source. The limit is what keeps a long ensemble, one context per
// thread, from walking into the process limit on open descriptors.
struct md_attribute_io_t {
    struct {
        uint64_t  hash;       // of the path; 0 marks an empty slot
        str_t     path;       // owned by alloc
        md_file_t file;
        uint64_t  last_use;
    } slot[ATTR_IO_MAX_FILES];
    uint64_t        tick;
    md_allocator_i* alloc;
};

static void attr_io_close_slot(md_attribute_io_t* io, size_t i) {
    if (io->slot[i].hash) {
        md_file_close(&io->slot[i].file);
        str_free(io->slot[i].path, io->alloc);
        MEMSET(&io->slot[i], 0, sizeof(io->slot[i]));
    }
}

static void attr_io_close_all(md_attribute_io_t* io) {
    for (size_t i = 0; i < ATTR_IO_MAX_FILES; ++i) {
        attr_io_close_slot(io, i);
    }
}

size_t md_attribute_io_read_at(md_attribute_io_t* io, str_t path, int64_t offset, void* dst, size_t bytes) {
    if (!io) {
        md_file_t file = {0};
        if (!md_file_open(&file, path, MD_FILE_READ)) {
            MD_LOG_ERROR("Failed to open '" STR_FMT "'", STR_ARG(path));
            return 0;
        }
        const size_t read = md_file_read_at(file, offset, dst, bytes);
        md_file_close(&file);
        return read;
    }

    uint64_t hash = md_hash64(path.ptr, path.len, 0);
    if (hash == 0) hash = 1;

    size_t idx = SIZE_MAX;
    size_t lru = 0;
    for (size_t i = 0; i < ATTR_IO_MAX_FILES; ++i) {
        if (io->slot[i].hash == hash && str_eq(io->slot[i].path, path)) {
            idx = i;
            break;
        }
        if (io->slot[i].last_use < io->slot[lru].last_use) {
            lru = i;
        }
    }

    if (idx == SIZE_MAX) {
        // An empty slot has last_use 0, so the least recently used is also the first empty one.
        idx = lru;
        attr_io_close_slot(io, idx);
        md_file_t file = {0};
        if (!md_file_open(&file, path, MD_FILE_READ)) {
            MD_LOG_ERROR("Failed to open '" STR_FMT "'", STR_ARG(path));
            return 0;
        }
        io->slot[idx].hash = hash;
        io->slot[idx].path = str_copy(path, io->alloc);
        io->slot[idx].file = file;
    }

    io->slot[idx].last_use = ++io->tick;
    return md_file_read_at(io->slot[idx].file, offset, dst, bytes);
}

// ### EXTRACTION ###

typedef enum extract_kind_t {
    EXTRACT_POSITION,   // atom/position into x, y, z
    EXTRACT_CELL,       // unitcell into unitcell
    EXTRACT_OTHER,      // anything else into out->attributes
} extract_kind_t;

typedef struct extract_entry_t {
    str_t             path;      // relative to the run, owned
    extract_kind_t    kind;
    md_attribute_id_t id;
    md_attribute_id_t axis_id;
} extract_entry_t;

struct md_system_extract_t {
    const md_system_t*     sys;
    md_allocator_i*        arena;       // everything the context owns
    md_attribute_id_t      run_axis_id;
    size_t                 num_frames;

    extract_entry_t*       entries;
    size_t                 num_entries;

    md_attribute_io_t      io;
    md_thread_id_t         owner;       // the thread that extracts; 0 until the first frame
};

static str_t extract_run_path(char* buf, size_t cap, str_t run, str_t leaf) {
    const int len = snprintf(buf, cap, STR_FMT "/" STR_FMT, STR_ARG(run), STR_ARG(leaf));
    return (len > 0 && (size_t)len < cap) ? (str_t){buf, (size_t)len} : (str_t){0};
}

md_system_extract_t* md_system_extract_begin(const md_system_t* sys, str_t run, const str_t paths[], size_t num_paths, md_allocator_i* alloc) {
    ASSERT(sys);
    ASSERT(alloc);
    ASSERT(paths || num_paths == 0);

    md_allocator_i* arena = md_arena_allocator_create(alloc, KILOBYTES(16));
    md_system_extract_t* ex = md_alloc(arena, sizeof(md_system_extract_t));
    MEMSET(ex, 0, sizeof(md_system_extract_t));
    ex->sys      = sys;
    ex->arena    = arena;
    ex->io.alloc = arena;
    ex->entries  = md_alloc(arena, MAX(num_paths, 1) * sizeof(extract_entry_t));

    const md_attributes_t* attributes = &sys->attributes;
    char buf[512];

    {
        const md_attribute_t* axis = str_empty(run) ? NULL : md_attributes_find(attributes, extract_run_path(buf, sizeof(buf), run, STR_LIT("time")));
        if (!axis || md_attributes_axis(attributes, axis) != axis) {
            MD_LOG_ERROR("No run '" STR_FMT "' to extract from", STR_ARG(run));
            goto fail;
        }
        ex->run_axis_id = axis->id;
        ex->num_frames  = axis->format.shape[0];
    }

    for (size_t i = 0; i < num_paths; ++i) {
        extract_entry_t* e = &ex->entries[ex->num_entries];
        MEMSET(e, 0, sizeof(*e));
        e->path = str_copy(paths[i], arena);
        e->kind = str_eq(paths[i], STR_LIT("atom/position")) ? EXTRACT_POSITION :
                  str_eq(paths[i], STR_LIT("unitcell"))      ? EXTRACT_CELL : EXTRACT_OTHER;

        const md_attribute_t* attr = md_attributes_find(attributes, extract_run_path(buf, sizeof(buf), run, paths[i]));
        if (!attr) {
            MD_LOG_ERROR("Nothing at '" STR_FMT "' in '" STR_FMT "' to extract", STR_ARG(paths[i]), STR_ARG(run));
            goto fail;
        }

        if (!(attr->flags & MD_ATTRIBUTE_FLAG_TEMPORAL)) {
            MD_LOG_ERROR("'" STR_FMT "' does not vary over the run; read it from the system directly", STR_ARG(attr->path));
            goto fail;
        }
        if (attr->format.type == MD_ATTRIBUTE_TYPE_STR) {
            MD_LOG_ERROR("'" STR_FMT "' is text and cannot be carried by a state", STR_ARG(attr->path));
            goto fail;
        }
        if (e->kind == EXTRACT_POSITION && !(attr->format.rank == 2 && attr->format.components == 3)) {
            MD_LOG_ERROR("'" STR_FMT "' is not {F,N} with three components", STR_ARG(attr->path));
            goto fail;
        }
        if (e->kind == EXTRACT_CELL && !(attr->format.rank == 3 && attr->format.shape[1] == 3 && attr->format.shape[2] == 3 && attr->format.components == 1)) {
            MD_LOG_ERROR("'" STR_FMT "' is not {F,3,3}", STR_ARG(attr->path));
            goto fail;
        }
        const md_attribute_t* axis = md_attributes_axis(attributes, attr);
        if (!axis) {
            MD_LOG_ERROR("'" STR_FMT "' has no frame axis", STR_ARG(attr->path));
            goto fail;
        }
        e->id      = attr->id;
        e->axis_id = axis->id;
        ex->num_entries += 1;
    }

    return ex;

fail:
    md_arena_allocator_destroy(arena);
    return NULL;
}

void md_system_extract_end(md_system_extract_t* ex) {
    if (!ex) return;
    attr_io_close_all(&ex->io);
    md_arena_allocator_destroy(ex->arena);
}

// The value of attr at one row in its STORED type, which a state keeps it in: secondary structure
// labels stay integers. Resident storage is copied; a provider is asked, with the context's io.
static bool extract_raw_row(void* dst, const md_attribute_t* attr, const md_attribute_slice_t* slice, md_attribute_io_t* io) {
    size_t first, count;
    if (!attr_slice_window(&first, &count, attr, slice)) {
        return false;
    }
    if (attr->virt.provider) {
        return attr->virt.provider(dst, count, attr, slice, attr->virt.user_data, io) == count;
    }
    if (!attr->data) {
        return false;
    }
    const size_t type_size = md_attribute_type_size(attr->format.type);
    MEMCPY(dst, (const uint8_t*)attr->data + first * type_size, count * type_size);
    return true;
}

bool md_system_extract_frame(md_system_extract_t* ex, int64_t frame, md_system_state_t* out) {
    ASSERT(ex);
    ASSERT(out);

    // One thread at a time. Recorded on first use rather than at begin, so a context may be made on
    // one thread and handed to the one that uses it.
    const md_thread_id_t tid = md_thread_id();
    if (ex->owner == 0) {
        ex->owner = tid;
    }
    ASSERT(ex->owner == tid && "an extraction context is used by one thread at a time");

    if (frame < 0 || (size_t)frame >= ex->num_frames) {
        MD_LOG_ERROR("Frame %" PRId64 " is outside the %zu frames being extracted from", frame, ex->num_frames);
        return false;
    }

    const md_attributes_t* attributes = &ex->sys->attributes;
    const bool want_coords = out->xyz != NULL;

    const md_attribute_t* run_axis = md_attributes_get(attributes, ex->run_axis_id);
    if (!run_axis) {
        MD_LOG_ERROR("The run being extracted from is gone; end the context before removing it");
        return false;
    }

    md_temp_scope_t temp = md_temp_begin();
    bool result = true;

    for (size_t i = 0; i < ex->num_entries && result; ++i) {
        const extract_entry_t* e = &ex->entries[i];
        const md_attribute_t* attr = md_attributes_get(attributes, e->id);
        const md_attribute_t* axis = md_attributes_get(attributes, e->axis_id);
        if (!attr || !axis) {
            MD_LOG_ERROR("'" STR_FMT "' is gone from the run being extracted from", STR_ARG(e->path));
            result = false;
            break;
        }

        size_t row = (size_t)frame;
        if (e->axis_id != ex->run_axis_id && !md_attribute_axis_map(&row, run_axis, (size_t)frame, axis)) {
            if (e->kind != EXTRACT_OTHER) {
                MD_LOG_ERROR("'" STR_FMT "' has no value at frame %" PRId64, STR_ARG(e->path), frame);
                result = false;
                break;
            }
            // Sampled at other times than the frames (velocities written every fifth frame): this
            // frame has none, and the state says so by not carrying it - never the value of another
            // frame, which a state reused from the last one would otherwise still hold.
            if (out->attributes.alloc) {
                const md_attribute_t* prev = md_attributes_find(&out->attributes, e->path);
                if (prev) md_attributes_remove(&out->attributes, prev->id);
            }
            continue;
        }
        const md_attribute_slice_t slice = md_attribute_slice_1((uint32_t)row);

        switch (e->kind) {
        case EXTRACT_POSITION: {
            if (!want_coords) break;
            const size_t N = attr->format.shape[1];
            if (out->num_atoms != 0 && out->num_atoms != N) {
                MD_LOG_ERROR("The state holds %zu atoms, '" STR_FMT "' %zu", out->num_atoms, STR_ARG(attr->path), N);
                result = false;
                break;
            }
            if (attr_extract_slice_io_F32((float*)out->xyz, N * 3, attr, &slice, md_unit_angstrom(), &ex->io) != N * 3) {
                result = false;
                break;
            }
            out->num_atoms = N;
            break;
        }
        case EXTRACT_CELL: {
            float box[3][3];
            if (attr_extract_slice_io_F32(&box[0][0], 9, attr, &slice, md_unit_angstrom(), &ex->io) != 9) {
                result = false;
                break;
            }
            const bool empty = box[0][0] == 0.0f && box[1][1] == 0.0f && box[2][2] == 0.0f;
            out->unitcell = empty ? (md_unitcell_t){0} : md_unitcell_from_matrix_float(box);
            break;
        }
        case EXTRACT_OTHER: {
            if (!out->attributes.alloc) {
                MD_LOG_ERROR("'" STR_FMT "' goes into the state's attributes, and the state has none (no allocator)", STR_ARG(e->path));
                result = false;
                break;
            }
            md_attribute_format_t fmt;
            md_attribute_slice_format(&fmt, attr, &slice);
            const size_t bytes = md_attribute_byte_size(&fmt);
            void* data = md_temp_alloc(temp, MAX(bytes, 1));
            if (!data || !extract_raw_row(data, attr, &slice, &ex->io)) {
                result = false;
                break;
            }
            const md_attribute_desc_t desc = {
                .path   = e->path,
                .format = fmt,
                .unit   = attr->unit,
                .label  = attr->label,
                .data   = data,
                .byte_size = bytes,
            };
            if (!md_attributes_replace(&out->attributes, &desc)) {
                result = false;
            }
            break;
        }
        }
    }

    md_temp_end(temp);
    if (result) {
        out->frame = (double)frame;
    }
    return result;
}

// ### RUNS ###

static str_t run_cache_path(char* buf, size_t cap, str_t source_path) {
    const int len = snprintf(buf, cap, STR_FMT ".cache", STR_ARG(source_path));
    return (len > 0 && (size_t)len < cap) ? (str_t){buf, (size_t)len} : (str_t){0};
}

bool md_run_cache_open(md_file_t* out_file, md_run_cache_header_t* out_header, str_t source_path, uint64_t magic, uint64_t version) {
    ASSERT(out_file);
    ASSERT(out_header);
    MEMSET(out_file, 0, sizeof(*out_file));

    md_file_info_t info = {0};
    if (!md_file_info_extract_from_path(source_path, &info)) {
        return false;
    }
    char buf[4096];
    const str_t cache_path = run_cache_path(buf, sizeof(buf), source_path);
    md_file_t file = {0};
    if (str_empty(cache_path) || !md_file_open(&file, cache_path, MD_FILE_READ)) {
        return false;
    }

    const char* stale = NULL;
    if (md_file_read(file, out_header, sizeof(*out_header)) != sizeof(*out_header)) {
        stale = "incomplete header";
    } else if (out_header->magic != magic) {
        stale = "not this format's";
    } else if (out_header->version != version) {
        stale = "an older version";
    } else if (out_header->source_size != (uint64_t)info.size) {
        stale = "made from a file of another size";
    } else if (out_header->source_modified != info.modified_time) {
        stale = "made from a file modified at another time";
    } else if (out_header->num_frames == 0 || out_header->num_atoms == 0) {
        stale = "empty";
    }
    if (stale) {
        MD_LOG_INFO("'" STR_FMT "' is %s; it will be made again", STR_ARG(cache_path), stale);
        md_file_close(&file);
        return false;
    }
    *out_file = file;
    return true;
}

bool md_run_cache_create(md_file_t* out_file, str_t source_path, const md_file_info_t* scanned, uint64_t magic, uint64_t version, size_t num_atoms, size_t num_frames) {
    ASSERT(out_file);
    ASSERT(scanned);
    MEMSET(out_file, 0, sizeof(*out_file));
    const md_file_info_t info = *scanned;
    char buf[4096];
    const str_t cache_path = run_cache_path(buf, sizeof(buf), source_path);
    md_file_t file = {0};
    if (str_empty(cache_path) || !md_file_open(&file, cache_path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
        MD_LOG_INFO("Could not create '" STR_FMT "'", STR_ARG(cache_path));
        return false;
    }
    const md_run_cache_header_t header = {
        .magic           = magic,
        .version         = version,
        .source_size     = (uint64_t)info.size,
        .source_modified = info.modified_time,
        .num_atoms       = num_atoms,
        .num_frames      = num_frames,
    };
    if (md_file_write(file, &header, sizeof(header)) != sizeof(header)) {
        md_file_close(&file);
        return false;
    }
    *out_file = file;
    return true;
}

str_t md_run_path(char* buf, size_t cap, str_t run, str_t leaf) {
    return extract_run_path(buf, cap, run, leaf);
}

bool md_run_publish(md_system_t* sys, str_t run, const md_run_desc_t* d) {
    ASSERT(sys);
    ASSERT(d);
    if (str_empty(run) || !d->time || !d->source_offset || !d->source_size || !d->position_virt || str_empty(d->source_path)) {
        MD_LOG_ERROR("An incomplete description of the run '" STR_FMT "'", STR_ARG(run));
        return false;
    }
    const size_t F = d->num_frames;
    const size_t N = d->num_atoms;
    if (F == 0 || N == 0 || F > UINT32_MAX || N > UINT32_MAX) {
        MD_LOG_ERROR("The run '" STR_FMT "' has %zu frames of %zu atoms, which cannot be published", STR_ARG(run), F, N);
        return false;
    }
    if (sys->atom.count != 0 && sys->atom.count != N) {
        MD_LOG_ERROR("The run '" STR_FMT "' holds %zu atoms, the system %zu", STR_ARG(run), N, sys->atom.count);
        return false;
    }

    md_attributes_t* attributes = &sys->attributes;
    if (!attributes->alloc) {
        attributes->alloc = sys->alloc;
    }

    const md_attribute_format_t series_f64 = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)F } };
    const md_attribute_format_t series_i64 = { .type = MD_ATTRIBUTE_TYPE_I64, .components = 1, .rank = 1, .shape = { (uint32_t)F } };
    const md_attribute_format_t cell_format = { .type = MD_ATTRIBUTE_TYPE_F32, .components = 1, .rank = 3, .shape = { (uint32_t)F, 3, 3 } };
    char buf[512];

    // The axis first: everything temporal below is checked against it.
    bool ok = md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), run, STR_LIT("time")), .format = series_f64, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
        .unit = d->time_unit, .label = STR_LIT("Time"),
        .data = d->time, .byte_size = F * sizeof(double)});

    if (ok && d->step) {
        ok = md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), run, STR_LIT("step")), .format = series_i64, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
            .unit = md_unit_none(), .label = STR_LIT("Step"),
            .description = STR_LIT("The file's own notion of where a frame sits in the run, not the frame ordinal"),
            .data = d->step, .byte_size = F * sizeof(int64_t)});
    }

    if (ok && (d->unitcell || d->unitcell_virt)) {
        ok = md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), run, STR_LIT("unitcell")), .format = cell_format,
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = md_unit_angstrom(), .label = STR_LIT("Unit Cell"),
            .data = d->unitcell_virt ? NULL : d->unitcell, .byte_size = d->unitcell_virt ? 0 : F * 9 * sizeof(float),
            .virt = d->unitcell_virt});
    }

    ok = ok && md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), run, STR_LIT("source/path")),
        .format = { .type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 0 },
        .unit = md_unit_none(), .data = &d->source_path, .byte_size = sizeof(str_t)});

    ok = ok && md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), run, STR_LIT("source/offset")), .format = series_i64, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
        .unit = md_unit_none(), .data = d->source_offset, .byte_size = F * sizeof(int64_t)});

    ok = ok && md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), run, STR_LIT("source/size")), .format = series_i64, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
        .unit = md_unit_none(), .data = d->source_size, .byte_size = F * sizeof(int64_t)});

    ok = ok && md_attributes_replace(attributes, &(md_attribute_desc_t){
        .path = md_run_path(buf, sizeof(buf), run, STR_LIT("atom/position")),
        .format = { .type = MD_ATTRIBUTE_TYPE_F32, .components = 3, .rank = 2, .shape = { (uint32_t)F, (uint32_t)N } },
        .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = md_unit_angstrom(), .label = STR_LIT("Position"),
        .virt = d->position_virt});

    if (!ok) {
        MD_LOG_ERROR("Failed to publish the run '" STR_FMT "' from '" STR_FMT "'", STR_ARG(run), STR_ARG(d->source_path));
        md_attributes_remove_prefix(attributes, run);
    }
    return ok;
}

bool md_run_source(md_run_source_t* out, const md_attributes_t* attributes, const md_attribute_t* attr, str_t leaf) {
    ASSERT(out);
    ASSERT(attributes);
    ASSERT(attr);
    MEMSET(out, 0, sizeof(*out));

    const str_t p = attr->path;
    if (p.len <= leaf.len + 1 || !str_ends_with(p, leaf) || p.ptr[p.len - leaf.len - 1] != '/') {
        return false;
    }
    out->run = str_substr(p, 0, p.len - leaf.len - 1);

    char buf[512];
    const md_attribute_t* path   = md_attributes_find(attributes, md_run_path(buf, sizeof(buf), out->run, STR_LIT("source/path")));
    const md_attribute_t* offset = md_attributes_find(attributes, md_run_path(buf, sizeof(buf), out->run, STR_LIT("source/offset")));
    const md_attribute_t* size   = md_attributes_find(attributes, md_run_path(buf, sizeof(buf), out->run, STR_LIT("source/size")));
    if (!path || path->format.type != MD_ATTRIBUTE_TYPE_STR ||
        !offset || offset->format.type != MD_ATTRIBUTE_TYPE_I64 || !offset->data ||
        !size   || size->format.type   != MD_ATTRIBUTE_TYPE_I64 || !size->data ||
        offset->format.shape[0] != size->format.shape[0]) {
        MD_LOG_ERROR("The run '" STR_FMT "' has lost its source attributes", STR_ARG(out->run));
        return false;
    }
    out->path       = md_attribute_str(attributes, path, 0);
    out->offset     = (const int64_t*)offset->data;
    out->size       = (const int64_t*)size->data;
    out->num_frames = offset->format.shape[0];
    return true;
}

// Lower case letters and digits, anything else one '_' between them
static size_t run_series_slug(char* buf, size_t cap, str_t name) {
    size_t len = 0;
    bool pending_sep = false;
    for (size_t i = 0; i < name.len && len + 2 < cap; ++i) {
        char c = name.ptr[i];
        if (c >= 'A' && c <= 'Z') c = (char)(c - 'A' + 'a');
        const bool keep = (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9');
        if (!keep) {
            pending_sep = (len > 0);
            continue;
        }
        if (pending_sep) {
            buf[len++] = '_';
            pending_sep = false;
        }
        buf[len++] = c;
    }
    buf[len] = '\0';
    return len;
}

bool md_run_publish_series(md_system_t* sys, str_t run, const md_run_series_desc_t* d) {
    ASSERT(sys);
    ASSERT(d);
    md_attributes_t* attributes = &sys->attributes;
    if (!attributes->alloc) {
        attributes->alloc = sys->alloc;
    }

    char buf[512];
    const md_attribute_t* run_axis = str_empty(run) ? NULL : md_attributes_find(attributes, md_run_path(buf, sizeof(buf), run, STR_LIT("time")));
    if (!run_axis || md_attributes_axis(attributes, run_axis) != run_axis) {
        MD_LOG_ERROR("No run '" STR_FMT "' to publish '" STR_FMT "' along", STR_ARG(run), STR_ARG(d->group));
        return false;
    }
    const size_t F = run_axis->format.shape[0];
    const size_t R = d->num_rows;
    if (R == 0 || d->num_columns == 0 || R > UINT32_MAX || !d->names || !d->columns) {
        MD_LOG_ERROR("'" STR_FMT "' holds no values", STR_ARG(d->group));
        return false;
    }

    char group_buf[512];
    const str_t group = md_run_path(group_buf, sizeof(group_buf), run, d->group);
    if (str_empty(group) || str_empty(d->group)) {
        MD_LOG_ERROR("No group to publish the series under");
        return false;
    }

    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    bool result = false;

    md_unit_t time_unit = d->time_unit;
    if (d->time) {
        if (md_unit_is_none(time_unit)) {
            time_unit = run_axis->unit;
        }
        for (size_t i = 1; i < R; ++i) {
            if (d->time[i] < d->time[i - 1]) {
                MD_LOG_ERROR("'" STR_FMT "': time decreases at row %zu", STR_ARG(d->group), i);
                goto done;
            }
        }
        // Whether the rows cover the run is a question about the values, answered against a scratch
        // table before the real one is touched.
        md_attributes_t probe = { .alloc = temp_alloc };
        const md_attribute_id_t probe_id = md_attributes_create(&probe, &(md_attribute_desc_t){
            .path = STR_LIT("time"), .format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)R } },
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = time_unit, .data = d->time, .byte_size = R * sizeof(double)});
        const md_attribute_t* axis = probe_id ? md_attributes_get(&probe, probe_id) : NULL;
        if (!axis) goto done;
        for (size_t f = 0; f < F; ++f) {
            size_t row;
            if (!md_attribute_axis_map(&row, run_axis, f, axis)) {
                MD_LOG_ERROR("'" STR_FMT "' has no row at the time of frame %zu of '" STR_FMT "'; it belongs to another run", STR_ARG(d->group), f, STR_ARG(run));
                goto done;
            }
        }
    } else if (R != F) {
        MD_LOG_ERROR("'" STR_FMT "' has %zu rows and no time, and '" STR_FMT "' %zu frames: the rows cannot be the frames", STR_ARG(d->group), R, STR_ARG(run), F);
        goto done;
    }

    md_attributes_remove_prefix(attributes, group);

    const md_attribute_format_t series = { .type = MD_ATTRIBUTE_TYPE_F32, .components = 1, .rank = 1, .shape = { (uint32_t)R } };
    bool ok = true;
    if (d->time) {
        ok = md_attributes_create(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), group, STR_LIT("time")),
            .format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)R } },
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = time_unit, .label = STR_LIT("Time"),
            .data = d->time, .byte_size = R * sizeof(double)}) != MD_ATTRIBUTE_INVALID;
    }

    md_array(str_t) taken = 0;
    md_array_push(taken, STR_LIT("time"),   temp_alloc);
    md_array_push(taken, STR_LIT("source"), temp_alloc);
    for (size_t k = 0; ok && k < d->num_columns; ++k) {
        char slug[128];
        size_t len = run_series_slug(slug, sizeof(slug) - 8, d->names[k]);
        if (len == 0) {
            len = (size_t)snprintf(slug, sizeof(slug), "column%zu", k + 1);
        }
        const size_t base = len;
        for (int n = 2;; ++n) {
            bool clash = false;
            for (size_t t = 0; t < md_array_size(taken); ++t) {
                if (str_eq(taken[t], (str_t){slug, len})) { clash = true; break; }
            }
            if (!clash) break;
            len = base + (size_t)snprintf(slug + base, sizeof(slug) - base, "_%d", n);
        }
        md_array_push(taken, str_copy((str_t){slug, len}, temp_alloc), temp_alloc);

        ok = md_attributes_create(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), group, (str_t){slug, len}),
            .format = series, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
            .unit = d->units ? d->units[k] : md_unit_none(), .label = d->names[k],
            .data = d->columns[k], .byte_size = R * sizeof(float)}) != MD_ATTRIBUTE_INVALID;
    }
    if (ok && !str_empty(d->source_path)) {
        ok = md_attributes_create(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), group, STR_LIT("source")),
            .format = { .type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 0 },
            .unit = md_unit_none(), .data = &d->source_path, .byte_size = sizeof(str_t)}) != MD_ATTRIBUTE_INVALID;
    }
    if (!ok) {
        MD_LOG_ERROR("Failed to publish '" STR_FMT "'", STR_ARG(group));
        md_attributes_remove_prefix(attributes, group);
        goto done;
    }
    result = true;

done:
    md_temp_end(temp);
    return result;
}

#ifdef __cplusplus
}
#endif
