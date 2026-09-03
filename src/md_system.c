#include <md_system.h>
#include <md_trajectory.h>

#include <core/md_log.h>
#include <core/md_array.h>
#include <core/md_hash.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>

#ifdef __cplusplus
extern "C" {
#endif

void md_system_free(md_system_t* sys) {
    ASSERT(sys);
    ASSERT(sys->alloc);
    md_allocator_i* alloc = sys->alloc;
    md_trajectory_free(sys->trajectory);

    // ATOM
    md_array_free(sys->atom.type_idx, alloc);
    md_array_free(sys->atom.flags, alloc);

    // ATOM TYPE
    md_array_free(sys->atom.type.name, alloc);
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

    // A freshly initialised state did not come from a trajectory. Only md_trajectory_reader_load_frame
    // and the interpolation which produces a state write a non negative frame. See md_system_state_t.
    state->frame = -1.0;

    if (num_atoms == 0) {
        return true;
    }

    const size_t capacity = ALIGN_TO(num_atoms, 16);

    md_array_resize(state->x, capacity, alloc);
    md_array_resize(state->y, capacity, alloc);
    md_array_resize(state->z, capacity, alloc);

    // Zero the padding past num_atoms. The capacity is rounded up so vectorised code may load whole
    // 16 wide chunks; an uninitialised tail puts garbage floats into those lanes.
    const size_t tail_bytes = (capacity - num_atoms) * sizeof(float);
    if (tail_bytes > 0) {
        MEMSET(state->x + num_atoms, 0, tail_bytes);
        MEMSET(state->y + num_atoms, 0, tail_bytes);
        MEMSET(state->z + num_atoms, 0, tail_bytes);
    }

    state->num_atoms = num_atoms;

    return true;
}

void md_system_state_free(md_system_state_t* state) {
    ASSERT(state);

    // A view owns nothing; zeroing it is the whole job.
    if (state->alloc) {
        md_array_free(state->x, state->alloc);
        md_array_free(state->y, state->alloc);
        md_array_free(state->z, state->alloc);
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
    if (src->num_atoms > 0 && src->x && src->y && src->z) {
        MEMCPY(dst->x, src->x, src->num_atoms * sizeof(float));
        MEMCPY(dst->y, src->y, src->num_atoms * sizeof(float));
        MEMCPY(dst->z, src->z, src->num_atoms * sizeof(float));
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

// Attach a trajectory to the system, freeing any existing attached trajectory.
void md_system_attach_trajectory(md_system_t* sys, struct md_trajectory_i* traj) {
    if (!sys) return;
    md_trajectory_free(sys->trajectory);
    sys->trajectory = traj;
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
                                          size_t first, size_t count, const md_attribute_slice_t* slice, md_unit_t dst_unit) { \
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
        size_t written = attr->virt.provider(buf, count, attr, slice, attr->virt.user_data);        \
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

MD_ATTR_DEFINE_EXTRACT_RANGE(I32, int32_t)
MD_ATTR_DEFINE_EXTRACT_RANGE(I64, int64_t)
MD_ATTR_DEFINE_EXTRACT_RANGE(U32, uint32_t)
MD_ATTR_DEFINE_EXTRACT_RANGE(U64, uint64_t)
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
    return attr_extract_range_F32(dst, cap, attr, 0, md_attribute_element_count(&attr->format), NULL, dst_unit);
}

size_t md_attribute_extract_f64(double dst[], size_t cap, const md_attribute_t* attr, md_unit_t dst_unit) {
    ASSERT(attr);
    if (attr_reject_whole_temporal(attr)) return 0;
    return attr_extract_range_F64(dst, cap, attr, 0, md_attribute_element_count(&attr->format), NULL, dst_unit);
}

size_t md_attribute_extract_slice_f32(float dst[], size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, md_unit_t dst_unit) {
    ASSERT(attr);
    if ((!slice || slice->num_idx == 0) && attr_reject_whole_temporal(attr)) return 0;
    size_t first, count;
    if (!attr_slice_window(&first, &count, attr, slice)) {
        return 0;
    }
    return attr_extract_range_F32(dst, cap, attr, first, count, slice, dst_unit);
}

size_t md_attribute_extract_slice_f64(double dst[], size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, md_unit_t dst_unit) {
    ASSERT(attr);
    if ((!slice || slice->num_idx == 0) && attr_reject_whole_temporal(attr)) return 0;
    size_t first, count;
    if (!attr_slice_window(&first, &count, attr, slice)) {
        return 0;
    }
    return attr_extract_range_F64(dst, cap, attr, first, count, slice, dst_unit);
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

void md_attributes_free(md_attributes_t* attributes) {
    ASSERT(attributes);
    md_allocator_i* alloc = attributes->alloc;
    if (alloc) {
        for (size_t i = 0; i < md_array_size(attributes->attr); ++i) {
            md_attribute_t* attr = attributes->attr + i;
            str_free(attr->path, alloc);
            // Both are optional and an absent one is a zeroed str_t, which is not something
            // to hand to an allocator.
            if (attr->label.ptr)       str_free(attr->label, alloc);
            if (attr->description.ptr) str_free(attr->description, alloc);
            // An alias shares the owner's buffer and frees none of it. Its user_data_size is zero
            // for the same reason, so the guard below already leaves the provider state alone.
            if (attr->storage == MD_ATTRIBUTE_STORAGE_RESIDENT && attr->data) {
                md_free(alloc, attr->data, md_attribute_byte_size(&attr->format));
            }
            // user_data_size is 0 for a borrowed pointer, so this only ever frees memory this
            // table itself handed out through md_attributes_alloc_user_data.
            if (attr->virt.user_data && attr->virt.user_data_size) {
                md_free(alloc, attr->virt.user_data, attr->virt.user_data_size);
            }
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

    if (!attributes->alloc) {
        MD_LOG_ERROR("Attribute table allocator not set");
        return MD_ATTRIBUTE_INVALID;
    }
    if (!attr_path_valid(path)) {
        MD_LOG_ERROR("Invalid attribute path '" STR_FMT "': expected non empty segments separated by '/'", STR_ARG(path));
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
    // coincidence, while a tag that disagrees with num_frames is a bug, and catching it here beats
    // finding it partway through an extract.
    if (desc->flags & MD_ATTRIBUTE_FLAG_TEMPORAL) {
        if (desc->format.rank == 0) {
            MD_LOG_ERROR("Attribute '" STR_FMT "' is temporal but has no index axes", STR_ARG(path));
            return MD_ATTRIBUTE_INVALID;
        }
        if (attributes->num_frames != 0 && desc->format.shape[0] != attributes->num_frames) {
            MD_LOG_ERROR("Attribute '" STR_FMT "' is temporal with an outermost extent of %u, but the table is indexed against %u frames",
                STR_ARG(path), desc->format.shape[0], attributes->num_frames);
            return MD_ATTRIBUTE_INVALID;
        }
    }

    size_t idx = attr_lower_bound(attributes, path);
    if (idx < md_array_size(attributes->attr) && str_eq(attributes->attr[idx].path, path)) {
        MD_LOG_ERROR("Attribute '" STR_FMT "' already exists", STR_ARG(path));
        return MD_ATTRIBUTE_INVALID;
    }

    md_attribute_id_t id = md_attributes_id_from_path(path);
    if (attr_index_from_id(attributes, id) != SIZE_MAX) {
        MD_LOG_ERROR("Hash collision for attribute '" STR_FMT "'", STR_ARG(path));
        return MD_ATTRIBUTE_INVALID;
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

    // Grow by one, then open a hole at idx so the array stays sorted by path.
    md_attribute_t empty = {0};
    md_array_push(attributes->attr, empty, alloc);

    size_t count = md_array_size(attributes->attr);
    if (idx + 1 < count) {
        MEMMOVE(attributes->attr + idx + 1, attributes->attr + idx, (count - 1 - idx) * sizeof(md_attribute_t));
    }

    attributes->attr[idx] = (md_attribute_t){
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

    return id;
}

md_attribute_id_t md_attributes_alias(md_attributes_t* attributes, md_attribute_id_t target, str_t path, str_t label, str_t description) {
    ASSERT(attributes);

    if (!attributes->alloc) {
        MD_LOG_ERROR("Attribute table allocator not set");
        return MD_ATTRIBUTE_INVALID;
    }
    if (!attr_path_valid(path)) {
        MD_LOG_ERROR("Invalid attribute path '" STR_FMT "': expected non empty segments separated by '/'", STR_ARG(path));
        return MD_ATTRIBUTE_INVALID;
    }

    size_t target_idx = attr_index_from_id(attributes, target);
    if (target_idx == SIZE_MAX) {
        MD_LOG_ERROR("Cannot alias '" STR_FMT "': no such target attribute", STR_ARG(path));
        return MD_ATTRIBUTE_INVALID;
    }

    size_t idx = attr_lower_bound(attributes, path);
    if (idx < md_array_size(attributes->attr) && str_eq(attributes->attr[idx].path, path)) {
        MD_LOG_ERROR("Attribute '" STR_FMT "' already exists", STR_ARG(path));
        return MD_ATTRIBUTE_INVALID;
    }

    md_attribute_id_t id = md_attributes_id_from_path(path);
    if (attr_index_from_id(attributes, id) != SIZE_MAX) {
        MD_LOG_ERROR("Hash collision for attribute '" STR_FMT "'", STR_ARG(path));
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

    md_attribute_t empty = {0};
    md_array_push(attributes->attr, empty, alloc);

    size_t count = md_array_size(attributes->attr);
    if (idx + 1 < count) {
        MEMMOVE(attributes->attr + idx + 1, attributes->attr + idx, (count - 1 - idx) * sizeof(md_attribute_t));
    }

    attributes->attr[idx] = (md_attribute_t){
        .id          = id,
        .path        = stored_path,
        .label       = stored_label,
        .description = stored_desc,
        .format      = format,
        .unit        = unit,
        .flags       = flags,
        .version     = ++attributes->version_counter,
        .storage     = MD_ATTRIBUTE_STORAGE_ALIAS,
        .data        = data,
        .virt        = virt,
        .root        = root,
    };

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

    md_attribute_t* attr = attributes->attr + idx;
    str_free(attr->path, alloc);
    if (attr->label.ptr)       str_free(attr->label, alloc);
    if (attr->description.ptr) str_free(attr->description, alloc);
    // Only an owner releases storage; an alias borrows both the buffer and the provider state.
    if (attr->storage == MD_ATTRIBUTE_STORAGE_RESIDENT && attr->data) {
        md_free(alloc, attr->data, md_attribute_byte_size(&attr->format));
    }
    if (attr->virt.user_data && attr->virt.user_data_size) {
        md_free(alloc, attr->virt.user_data, attr->virt.user_data_size);
    }

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
    // particular file leaves blank per atom. Two NANs are the same absence and IEEE does not say so,
    // hence the explicit test - and a column that is entirely absent is exactly as uninformative as
    // a constant one, so it is skipped for the same reason.
    bool uniform = true;
    bool all_absent = true;
    for (size_t i = 0; i < num_elements; ++i) {
        const size_t c = i % components;
        const bool absent_i = (values[i] != values[i]);
        const bool absent_0 = (values[c] != values[c]);
        if (!absent_i) {
            all_absent = false;
        }
        if (absent_i != absent_0 || (!absent_i && values[i] != values[c])) {
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
    attributes->attr[idx].version = ++attributes->version_counter;
    return attributes->attr[idx].version;
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

#ifdef __cplusplus
}
#endif
