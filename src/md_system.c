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
    return true;
}

void md_system_reset(md_system_t* sys) {
    ASSERT(sys);
    md_allocator_i* alloc = sys->alloc;
    md_system_free(sys);
    sys->alloc = alloc;
    // The table allocates through its own handle, so every loader gets a usable table
    // without having to remember to wire this up itself.
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
};

size_t md_attribute_type_size(md_attribute_type_t type) {
    if (type <= MD_ATTRIBUTE_TYPE_NONE || type >= MD_ATTRIBUTE_TYPE_COUNT) {
        return 0;
    }
    return attr_type_size[type];
}

size_t md_attribute_components(const md_attribute_format_t* format) {
    ASSERT(format);
    return format->rank < 2 ? 1 : (size_t)format->shape[format->rank - 1];
}

size_t md_attribute_value_count(const md_attribute_format_t* format) {
    ASSERT(format);
    if (format->rank < 2) {
        return format->rank == 0 ? 1 : (size_t)format->shape[0];
    }
    size_t count = 1;
    for (uint32_t i = 0; i < format->rank - 1; ++i) {
        count *= (size_t)format->shape[i];
    }
    return count;
}

size_t md_attribute_element_count(const md_attribute_format_t* format) {
    ASSERT(format);
    size_t count = 1;
    for (uint32_t i = 0; i < format->rank; ++i) {
        count *= (size_t)format->shape[i];
    }
    return count;
}

size_t md_attribute_byte_size(const md_attribute_format_t* format) {
    ASSERT(format);
    return md_attribute_element_count(format) * md_attribute_type_size(format->type);
}

// Splits at the LAST separator. loc is the separator offset, false if there is none.
static bool attr_path_split(size_t* loc, str_t path) {
    return str_rfind_char(loc, path, '/');
}

str_t md_attribute_category(const md_attribute_t* attr) {
    ASSERT(attr);
    size_t loc;
    if (!attr_path_split(&loc, attr->name)) {
        return (str_t){0};
    }
    return str_substr(attr->name, 0, loc);
}

str_t md_attribute_field(const md_attribute_t* attr) {
    ASSERT(attr);
    size_t loc;
    if (!attr_path_split(&loc, attr->name)) {
        return attr->name;
    }
    return str_substr(attr->name, loc + 1, SIZE_MAX);
}

size_t md_attribute_extract_f32(float dst[], size_t cap, const md_attribute_t* attr, md_unit_t dst_unit) {
    ASSERT(attr);

    if (!dst) {
        return 0;
    }

    size_t count = md_attribute_element_count(&attr->format);
    if (count > cap) {
        MD_LOG_ERROR("Attribute '" STR_FMT "' needs %zu floats, %zu supplied", STR_ARG(attr->name), count, cap);
        return 0;
    }
    if (!attr->data || count == 0) {
        return 0;
    }

    // md_unit_none() as the target means "as stored". Anything else has to be convertible, and a
    // refusal here is the point: a silently rescaled quantity is not detectable downstream.
    double factor = 1.0;
    if (!md_unit_is_none(dst_unit) && !md_unit_conversion_factor(&factor, attr->unit, dst_unit)) {
        char from[64], to[64];
        size_t from_len = md_unit_print(from, sizeof(from), attr->unit);
        size_t to_len   = md_unit_print(to,   sizeof(to),   dst_unit);
        MD_LOG_ERROR("Attribute '" STR_FMT "' is '%.*s' and cannot be expressed as '%.*s'",
            STR_ARG(attr->name), (int)from_len, from, (int)to_len, to);
        return 0;
    }

    if (attr->format.type == MD_ATTRIBUTE_TYPE_F32 && factor == 1.0) {
        MEMCPY(dst, attr->data, count * sizeof(float));
        return count;
    }

    // The scale is applied in double and narrowed once, at the end.
#define MD_ATTR_CONVERT(T)                                          \
    do {                                                            \
        const T* src = (const T*)attr->data;                        \
        for (size_t i = 0; i < count; ++i) {                        \
            dst[i] = (float)((double)src[i] * factor);              \
        }                                                           \
    } while (0)

    switch (attr->format.type) {
    case MD_ATTRIBUTE_TYPE_F32: MD_ATTR_CONVERT(float);    break;
    case MD_ATTRIBUTE_TYPE_F64: MD_ATTR_CONVERT(double);   break;
    case MD_ATTRIBUTE_TYPE_I8:  MD_ATTR_CONVERT(int8_t);   break;
    case MD_ATTRIBUTE_TYPE_U8:  MD_ATTR_CONVERT(uint8_t);  break;
    case MD_ATTRIBUTE_TYPE_I16: MD_ATTR_CONVERT(int16_t);  break;
    case MD_ATTRIBUTE_TYPE_U16: MD_ATTR_CONVERT(uint16_t); break;
    case MD_ATTRIBUTE_TYPE_I32: MD_ATTR_CONVERT(int32_t);  break;
    case MD_ATTRIBUTE_TYPE_U32: MD_ATTR_CONVERT(uint32_t); break;
    case MD_ATTRIBUTE_TYPE_I64: MD_ATTR_CONVERT(int64_t);  break;
    case MD_ATTRIBUTE_TYPE_U64: MD_ATTR_CONVERT(uint64_t); break;
    default:
        MD_LOG_ERROR("Attribute '" STR_FMT "' has no readable type", STR_ARG(attr->name));
        return 0;
    }

#undef MD_ATTR_CONVERT

    return count;
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
        if (str_cmp_lex(attributes->attr[mid].name, name) < 0) {
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
            str_free(attr->name, alloc);
            if (attr->data) {
                md_free(alloc, attr->data, md_attribute_byte_size(&attr->format));
            }
        }
        md_array_free(attributes->attr, alloc);
    }
    MEMSET(attributes, 0, sizeof(md_attributes_t));
}

size_t md_attributes_count(const md_attributes_t* attributes) {
    ASSERT(attributes);
    return md_array_size(attributes->attr);
}

md_attribute_id_t md_attributes_create(md_attributes_t* attributes, str_t name, md_attribute_format_t format, md_unit_t unit, const void* data, size_t byte_size) {
    ASSERT(attributes);

    if (!attributes->alloc) {
        MD_LOG_ERROR("Attribute table allocator not set");
        return MD_ATTRIBUTE_INVALID;
    }
    if (!attr_path_valid(name)) {
        MD_LOG_ERROR("Invalid attribute path '" STR_FMT "': expected non empty segments separated by '/'", STR_ARG(name));
        return MD_ATTRIBUTE_INVALID;
    }
    if (md_attribute_type_size(format.type) == 0) {
        MD_LOG_ERROR("Invalid type for attribute '" STR_FMT "'", STR_ARG(name));
        return MD_ATTRIBUTE_INVALID;
    }
    if (format.rank > MD_ATTRIBUTE_MAX_RANK) {
        MD_LOG_ERROR("Rank %u exceeds MD_ATTRIBUTE_MAX_RANK for attribute '" STR_FMT "'", format.rank, STR_ARG(name));
        return MD_ATTRIBUTE_INVALID;
    }
    for (uint32_t i = 0; i < format.rank; ++i) {
        if (format.shape[i] == 0) {
            MD_LOG_ERROR("Zero extent in axis %u of attribute '" STR_FMT "'", i, STR_ARG(name));
            return MD_ATTRIBUTE_INVALID;
        }
    }
    // Axes beyond the declared rank are not read anywhere, so make sure they cannot be
    // mistaken for data if someone inspects the format in a debugger.
    size_t required = md_attribute_byte_size(&format);
    if (data) {
        if (byte_size != required) {
            MD_LOG_ERROR("Attribute '" STR_FMT "' declares %zu bytes but %zu were supplied", STR_ARG(name), required, byte_size);
            return MD_ATTRIBUTE_INVALID;
        }
    } else if (byte_size != 0) {
        // A size without a pointer is a caller who meant to pass one.
        MD_LOG_ERROR("Attribute '" STR_FMT "' supplied %zu bytes with no data pointer", STR_ARG(name), byte_size);
        return MD_ATTRIBUTE_INVALID;
    }

    for (uint32_t i = format.rank; i < MD_ATTRIBUTE_MAX_RANK; ++i) {
        format.shape[i] = 0;
    }

    size_t idx = attr_lower_bound(attributes, name);
    if (idx < md_array_size(attributes->attr) && str_eq(attributes->attr[idx].name, name)) {
        MD_LOG_ERROR("Attribute '" STR_FMT "' already exists", STR_ARG(name));
        return MD_ATTRIBUTE_INVALID;
    }

    md_attribute_id_t id = (md_attribute_id_t)md_hash64_str(name, 0);
    if (id == MD_ATTRIBUTE_INVALID) {
        // Zero is the invalid value, so the one path that hashes to it borrows the next id.
        id = 1;
    }
    if (attr_index_from_id(attributes, id) != SIZE_MAX) {
        MD_LOG_ERROR("Hash collision for attribute '" STR_FMT "'", STR_ARG(name));
        return MD_ATTRIBUTE_INVALID;
    }

    md_allocator_i* alloc = attributes->alloc;

    void* storage = md_alloc(alloc, required);
    if (!storage) {
        MD_LOG_ERROR("Failed to allocate %zu bytes for attribute '" STR_FMT "'", required, STR_ARG(name));
        return MD_ATTRIBUTE_INVALID;
    }
    if (data) {
        MEMCPY(storage, data, required);
    } else {
        MEMSET(storage, 0, required);
    }

    str_t stored_name = str_copy(name, alloc);
    if (str_empty(stored_name)) {
        md_free(alloc, storage, required);
        MD_LOG_ERROR("Failed to copy attribute path '" STR_FMT "'", STR_ARG(name));
        return MD_ATTRIBUTE_INVALID;
    }

    // Grow by one, then open a hole at idx so the array stays sorted by name.
    md_attribute_t empty = {0};
    md_array_push(attributes->attr, empty, alloc);

    size_t count = md_array_size(attributes->attr);
    if (idx + 1 < count) {
        MEMMOVE(attributes->attr + idx + 1, attributes->attr + idx, (count - 1 - idx) * sizeof(md_attribute_t));
    }

    attributes->attr[idx] = (md_attribute_t){
        .id     = id,
        .name   = stored_name,
        .format = format,
        .unit   = unit,
        .data   = storage,
    };

    return id;
}

bool md_attributes_remove(md_attributes_t* attributes, md_attribute_id_t id) {
    ASSERT(attributes);

    size_t idx = attr_index_from_id(attributes, id);
    if (idx == SIZE_MAX) {
        return false;
    }

    md_allocator_i* alloc = attributes->alloc;
    ASSERT(alloc);

    md_attribute_t* attr = attributes->attr + idx;
    str_free(attr->name, alloc);
    if (attr->data) {
        md_free(alloc, attr->data, md_attribute_byte_size(&attr->format));
    }

    size_t count = md_array_size(attributes->attr);
    if (idx + 1 < count) {
        MEMMOVE(attributes->attr + idx, attributes->attr + idx + 1, (count - 1 - idx) * sizeof(md_attribute_t));
    }
    md_array_pop(attributes->attr);

    return true;
}

const md_attribute_t* md_attributes_get(const md_attributes_t* attributes, md_attribute_id_t id) {
    ASSERT(attributes);
    size_t idx = attr_index_from_id(attributes, id);
    return idx == SIZE_MAX ? NULL : attributes->attr + idx;
}

const md_attribute_t* md_attributes_find(const md_attributes_t* attributes, str_t name) {
    ASSERT(attributes);
    size_t idx = attr_lower_bound(attributes, name);
    if (idx < md_array_size(attributes->attr) && str_eq(attributes->attr[idx].name, name)) {
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
    if (attr->format.type != expected_type) {
        MD_LOG_ERROR("Type mismatch for attribute '" STR_FMT "'", STR_ARG(attr->name));
        return NULL;
    }
    return attr->data;
}

size_t md_attributes_query(md_attribute_id_t out_ids[], size_t cap, const md_attributes_t* attributes, str_t prefix) {
    ASSERT(attributes);

    str_t base = attr_prefix_trim(prefix);
    size_t count = 0;

    // Sorted by name, so the run starts at the lower bound of the prefix. A path may begin
    // with the prefix and continue with a character below '/' (say "atom-x" under "atom"),
    // which sorts inside the run without being covered by it, hence the two conditions.
    for (size_t i = str_empty(base) ? 0 : attr_lower_bound(attributes, base); i < md_array_size(attributes->attr); ++i) {
        str_t name = attributes->attr[i].name;
        if (!str_empty(base) && !str_begins_with(name, base)) {
            break;
        }
        if (!attr_path_covered_by(name, base)) {
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
        str_t name = attributes->attr[i].name;
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
