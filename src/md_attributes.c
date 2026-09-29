#include <md_attributes.h>
#include "md_attributes_internal.h"

#include <stdio.h>

#include <core/md_log.h>
#include <core/md_array.h>
#include <core/md_hash.h>
#include <core/md_allocator.h>
#include <core/md_os.h>

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

bool md_attribute_type_is_numeric(md_attribute_type_t type) {
    return type > MD_ATTRIBUTE_TYPE_NONE && type < MD_ATTRIBUTE_TYPE_COUNT && type != MD_ATTRIBUTE_TYPE_STR;
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

// The shared core, generated once per destination type: the one place where a stored type becomes a
// number and a unit becomes a factor. first and count are the window in ELEMENTS (components folded
// in); slice is the selection that produced them, forwarded to a provider as asked (NULL for the
// whole attribute) rather than as an offset it would have to reverse back into indices.
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
    /* Refused, not converted: a pool handle read as a number still looks like data. */            \
    if (attr->format.type == MD_ATTRIBUTE_TYPE_STR) {                                               \
        MD_LOG_ERROR("Attribute '" STR_FMT "' is a string; read it with md_attribute_extract_str", STR_ARG(attr->path)); \
        return 0;                                                                                   \
    }                                                                                               \
                                                                                                    \
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
    /* A virtual attribute (an alias of one included) is read through its provider into scratch  */ \
    /* of the STORED type, and from there on is converted exactly like a resident one. When no    */ \
    /* conversion is needed the provider writes straight into dst: for a coordinate array the     */ \
    /* extra copy was 5-8% of decoding the frame.                                                 */ \
    const bool as_stored = (attr->format.type == MD_ATTRIBUTE_TYPE_##SUFFIX && factor == 1.0);      \
    if (md_attribute_is_virtual(attr) && as_stored) {                                               \
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
    if (md_attribute_is_virtual(attr)) {                                                            \
        temp = md_temp_begin();                                                                     \
        own_temp = true;                                                                            \
        void* buf = md_temp_alloc(temp, count * md_attribute_type_size(attr->format.type));         \
        if (!buf) {                                                                                 \
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
        if (!attr->data) {                                                                          \
            MD_LOG_ERROR("Attribute '" STR_FMT "' has no data to read", STR_ARG(attr->path));        \
            return 0;                                                                               \
        }                                                                                           \
        src = (const uint8_t*)attr->data + first * md_attribute_type_size(attr->format.type);       \
    }                                                                                               \
                                                                                                    \
    size_t result = count;                                                                          \
    if (as_stored) {                                                                                \
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

// Only f32 and f64: an integer column is read through f64, which is exact to 2^53.
MD_ATTR_DEFINE_EXTRACT_RANGE(F32, float)
MD_ATTR_DEFINE_EXTRACT_RANGE(F64, double)

#undef MD_ATTR_CONVERT
#undef MD_ATTR_DEFINE_EXTRACT_RANGE

// The one piece of layout arithmetic in the library. Row major, so fixing the first num_idx axes
// selects one contiguous block. Reads only the FORMAT, never the storage.
bool md_attribute_slice_window(size_t* out_first, size_t* out_count, const md_attribute_t* attr, md_attribute_slice_t slice) {
    ASSERT(attr);

    const md_attribute_format_t* fmt = &attr->format;

    if (slice.num_idx > fmt->rank) {
        MD_LOG_ERROR("Attribute '" STR_FMT "' has rank %u, %u indices supplied", STR_ARG(attr->path), fmt->rank, slice.num_idx);
        return false;
    }

    size_t block = (size_t)fmt->components;
    for (uint32_t i = slice.num_idx; i < fmt->rank; ++i) {
        block *= (size_t)fmt->shape[i];
    }

    // Horner over the fixed axes gives the block ordinal; the block size turns it into elements.
    size_t ordinal = 0;
    for (uint32_t i = 0; i < slice.num_idx; ++i) {
        if (slice.idx[i] >= fmt->shape[i]) {
            MD_LOG_ERROR("Attribute '" STR_FMT "': index %u out of range on axis %u of extent %u",
                STR_ARG(attr->path), slice.idx[i], i, fmt->shape[i]);
            return false;
        }
        ordinal = ordinal * (size_t)fmt->shape[i] + (size_t)slice.idx[i];
    }

    *out_first = ordinal * block;
    *out_count = block;
    return true;
}

size_t md_attribute_slice_count(const md_attribute_t* attr, md_attribute_slice_t slice) {
    ASSERT(attr);
    size_t first, count;
    return md_attribute_slice_window(&first, &count, attr, slice) ? count : 0;
}

bool md_attribute_slice_format(md_attribute_format_t* out, const md_attribute_t* attr, md_attribute_slice_t slice) {
    ASSERT(out);
    ASSERT(attr);

    size_t first, count;
    if (!md_attribute_slice_window(&first, &count, attr, slice)) {
        return false;
    }

    // Slicing picks values and never splits one, so only the index axes change.
    const md_attribute_format_t* fmt = &attr->format;
    MEMSET(out, 0, sizeof(*out));
    out->type       = fmt->type;
    out->components = fmt->components;
    out->rank       = fmt->rank - slice.num_idx;
    for (uint32_t i = 0; i < out->rank; ++i) {
        out->shape[i] = fmt->shape[slice.num_idx + i];
    }
    return true;
}

// The whole of a VIRTUAL temporal attribute is every frame of it decoded into one buffer, which is
// never what a caller meant. The whole of a resident one is a copy of bytes that already exist
// ('run/x/time' for a plot), which is fine.
static bool attr_reject_whole_temporal(const md_attribute_t* attr, md_attribute_slice_t slice) {
    if (slice.num_idx == 0 && (attr->flags & MD_ATTRIBUTE_FLAG_TEMPORAL) && md_attribute_is_virtual(attr)) {
        MD_LOG_ERROR("Attribute '" STR_FMT "' is temporal and computed on demand: fix the frame axis with a slice rather than asking for every frame",
            STR_ARG(attr->path));
        return true;
    }
    return false;
}

#define MD_ATTR_DEFINE_EXTRACT(SUFFIX, suffix, DST_T)                                                     \
static size_t attr_extract_##SUFFIX(DST_T dst[], size_t cap, const md_attribute_t* attr, md_attribute_slice_t slice, md_unit_t dst_unit, md_attribute_io_t* io) { \
    ASSERT(attr);                                                                                   \
    if (attr_reject_whole_temporal(attr, slice)) return 0;                                          \
    size_t first, count;                                                                            \
    if (!md_attribute_slice_window(&first, &count, attr, slice)) return 0;                          \
    return attr_extract_range_##SUFFIX(dst, cap, attr, first, count, slice.num_idx ? &slice : NULL, dst_unit, io); \
}                                                                                                   \
                                                                                                    \
DST_T* md_attribute_extract_alloc_##suffix(size_t* out_count, const md_attribute_t* attr, md_attribute_slice_t slice, md_unit_t dst_unit, md_allocator_i* alloc) { \
    ASSERT(attr);                                                                                   \
    ASSERT(alloc);                                                                                  \
    if (out_count) *out_count = 0;                                                                  \
    const size_t count = md_attribute_slice_count(attr, slice);                                     \
    if (count == 0) return NULL;                                                                    \
    DST_T* dst = (DST_T*)md_alloc(alloc, count * sizeof(DST_T));                                    \
    if (!dst) return NULL;                                                                          \
    if (attr_extract_##SUFFIX(dst, count, attr, slice, dst_unit, NULL) != count) {                  \
        md_free(alloc, dst, count * sizeof(DST_T));                                                 \
        return NULL;                                                                                \
    }                                                                                               \
    if (out_count) *out_count = count;                                                              \
    return dst;                                                                                     \
}

MD_ATTR_DEFINE_EXTRACT(F32, f32, float)
MD_ATTR_DEFINE_EXTRACT(F64, f64, double)

#undef MD_ATTR_DEFINE_EXTRACT

size_t md_attribute_extract_f32(float dst[], size_t cap, const md_attribute_t* attr, md_attribute_slice_t slice, md_unit_t dst_unit) {
    return attr_extract_F32(dst, cap, attr, slice, dst_unit, NULL);
}

size_t md_attribute_extract_f64(double dst[], size_t cap, const md_attribute_t* attr, md_attribute_slice_t slice, md_unit_t dst_unit) {
    return attr_extract_F64(dst, cap, attr, slice, dst_unit, NULL);
}

size_t md_attribute_extract_io_f32(float dst[], size_t cap, const md_attribute_t* attr, md_attribute_slice_t slice, md_unit_t dst_unit, md_attribute_io_t* io) {
    return attr_extract_F32(dst, cap, attr, slice, dst_unit, io);
}

bool md_attribute_read_stored(void* dst, const md_attribute_t* attr, md_attribute_slice_t slice, md_attribute_io_t* io) {
    ASSERT(dst);
    ASSERT(attr);
    size_t first, count;
    if (!md_attribute_slice_window(&first, &count, attr, slice)) {
        return false;
    }
    if (md_attribute_is_virtual(attr)) {
        return attr->virt.provider(dst, count, attr, slice.num_idx ? &slice : NULL, attr->virt.user_data, io) == count;
    }
    if (!attr->data) {
        return false;
    }
    const size_t type_size = md_attribute_type_size(attr->format.type);
    MEMCPY(dst, (const uint8_t*)attr->data + first * type_size, count * type_size);
    return true;
}

const void* md_attribute_view(const md_attribute_t* attr, md_attribute_type_t type, uint32_t components, uint32_t rank) {
    if (!attr || md_attribute_is_virtual(attr) || !attr->data || type == MD_ATTRIBUTE_TYPE_STR) {
        return NULL;
    }
    const md_attribute_format_t* fmt = &attr->format;
    if (fmt->type != type || fmt->components != components || fmt->rank != rank) {
        return NULL;
    }
    return attr->data;
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

// Ids are what consumers are told to hold, so resolving one must not cost a scan of the table. The
// index is open addressing over the attribute array, keyed by the id (already a hash), and rebuilt
// after every insert and remove - both of which shift the array anyway. Should a rebuild fail to
// allocate, the index is left empty and lookup falls back to scanning.
static void attr_id_index_rebuild(md_attributes_t* attributes) {
    const size_t count = md_array_size(attributes->attr);
    size_t cap = 16;
    while (cap < count * 2) {
        cap *= 2;
    }
    if ((size_t)md_array_size(attributes->id_index) != cap) {
        md_array_resize(attributes->id_index, cap, attributes->alloc);
        if ((size_t)md_array_size(attributes->id_index) != cap) {
            md_array_free(attributes->id_index, attributes->alloc);
            attributes->id_index = 0;
            return;
        }
    }
    MEMSET(attributes->id_index, 0, cap * sizeof(uint32_t));
    const size_t mask = cap - 1;
    for (size_t i = 0; i < count; ++i) {
        size_t slot = (size_t)attributes->attr[i].id & mask;
        while (attributes->id_index[slot] != 0) {
            slot = (slot + 1) & mask;
        }
        attributes->id_index[slot] = (uint32_t)(i + 1);
    }
}

static size_t attr_index_from_id(const md_attributes_t* attributes, md_attribute_id_t id) {
    if (id == MD_ATTRIBUTE_INVALID) {
        return SIZE_MAX;
    }
    const size_t cap = md_array_size(attributes->id_index);
    if (cap == 0) {
        for (size_t i = 0; i < md_array_size(attributes->attr); ++i) {
            if (attributes->attr[i].id == id) {
                return i;
            }
        }
        return SIZE_MAX;
    }
    const size_t mask = cap - 1;
    for (size_t slot = (size_t)id & mask; attributes->id_index[slot] != 0; slot = (slot + 1) & mask) {
        const size_t i = attributes->id_index[slot] - 1;
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

// Ownership is derived rather than stated: an alias is root != id, and a virtual attribute is one
// with a provider. The one combination nothing can produce is storage AND a provider.
static bool attr_consistent(const md_attribute_t* attr) {
    return !(attr->data && attr->virt.provider);
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

// The one place an attribute enters the table. idx comes from attr_reserve_slot, so opening the
// hole keeps the array sorted by path.
static void attr_insert_at(md_attributes_t* attributes, size_t idx, const md_attribute_t* attr) {
    ASSERT(attr_consistent(attr));

    md_attribute_t empty = {0};
    md_array_push(attributes->attr, empty, attributes->alloc);

    size_t count = md_array_size(attributes->attr);
    if (idx + 1 < count) {
        MEMMOVE(attributes->attr + idx + 1, attributes->attr + idx, (count - 1 - idx) * sizeof(md_attribute_t));
    }
    attributes->attr[idx] = *attr;
    attr_id_index_rebuild(attributes);
}

// Everything one attribute owns. Removing one and tearing the whole table down have to agree about
// this, and when they were two copies of the list they were one edit away from disagreeing - which
// is a leak in one path or a double free in the other, neither visible until it is.
static void attr_release(md_attribute_t* attr, md_allocator_i* alloc) {
    ASSERT(attr);
    ASSERT(alloc);
    ASSERT(attr_consistent(attr));
    str_free(attr->path, alloc);
    // Both are optional and an absent one is a zeroed str_t, which is not something to hand to an
    // allocator.
    if (attr->label.ptr)       str_free(attr->label, alloc);
    if (attr->description.ptr) str_free(attr->description, alloc);
    // Only an owner releases storage; an alias borrows both the buffer and the provider state.
    if (md_attribute_is_alias(attr)) {
        return;
    }
    if (attr->data) {
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
        md_array_free(attributes->id_index, alloc);
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
    const md_attribute_virtual_t virt     = tgt->virt;   // still owned by the target: see attr_release

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
    attr_id_index_rebuild(attributes);
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
            if (md_attribute_is_alias(a) && a->root == id) {
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
    if (md_attribute_is_alias(attr)) {
        // Writing through a second name would be writing to somebody else's attribute behind its
        // back. Fill the owner in and every name sees it.
        MD_LOG_ERROR("Attribute '" STR_FMT "' is an alias; fill in the attribute which owns the storage", STR_ARG(attr->path));
        return NULL;
    }
    if (md_attribute_is_virtual(attr)) {
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

// Sorted by path, so the matches of a prefix start at its lower bound and run contiguously - except
// that a path may continue the prefix with a character below '/' ("atom-x" under "atom"), which sorts
// inside the run without being covered by it and is skipped.
static md_attribute_iter_t attr_iter_begin(const md_attributes_t* attributes, str_t prefix, bool children) {
    ASSERT(attributes);
    const str_t base = attr_prefix_trim(prefix);
    md_attribute_iter_t it = {
        ._table    = attributes,
        ._prefix   = base,
        ._next     = str_empty(base) ? 0 : attr_lower_bound(attributes, base),
        ._children = children,
    };
    return it;
}

md_attribute_iter_t md_attributes_iter(const md_attributes_t* attributes, str_t prefix) {
    return attr_iter_begin(attributes, prefix, false);
}

md_attribute_iter_t md_attributes_iter_children(const md_attributes_t* attributes, str_t prefix) {
    return attr_iter_begin(attributes, prefix, true);
}

bool md_attributes_next(md_attribute_iter_t* it) {
    ASSERT(it);
    const md_attributes_t* attributes = it->_table;
    if (!attributes) {
        return false;
    }
    const str_t  base  = it->_prefix;
    const size_t count = md_array_size(attributes->attr);

    while (it->_next < count) {
        const md_attribute_t* attr = attributes->attr + it->_next++;
        const str_t name = attr->path;
        if (!str_empty(base) && !str_begins_with(name, base)) {
            break;
        }
        if (!attr_path_covered_by(name, base)) {
            continue;
        }

        str_t child = {0};
        str_t child_path = {0};
        if (name.len > base.len) {
            const size_t offset = str_empty(base) ? 0 : base.len + 1;
            const str_t rest = str_substr(name, offset, SIZE_MAX);
            size_t loc;
            child = str_find_char(&loc, rest, '/') ? str_substr(rest, 0, loc) : rest;
            child_path = str_substr(name, 0, offset + child.len);
        }

        if (it->_children) {
            // The prefix itself is not a child of itself. Paths are sorted, so a child's paths are
            // contiguous and a repeat is the child just yielded - with one exception: a leaf "C" and
            // the paths below "C/" can be split by a sibling "C-x", since '-' sorts before '/'. The
            // leaf always comes first and has been yielded already, so its presence settles it.
            if (str_empty(child) || (it->attr && str_eq(child, it->child))) {
                continue;
            }
            if (name.len > child_path.len && md_attributes_find(attributes, child_path)) {
                continue;
            }
        }

        it->attr       = attr;
        it->child      = child;
        it->child_path = child_path;
        return true;
    }

    it->_table     = NULL;
    it->attr       = NULL;
    it->child      = (str_t){0};
    it->child_path = (str_t){0};
    return false;
}

size_t md_attributes_query(md_attribute_id_t out_ids[], size_t cap, const md_attributes_t* attributes, str_t prefix) {
    size_t count = 0;
    for (md_attribute_iter_t it = md_attributes_iter(attributes, prefix); md_attributes_next(&it);) {
        if (out_ids && count < cap) {
            out_ids[count] = it.attr->id;
        }
        count += 1;
    }
    return count;
}

const md_attribute_t* md_attributes_find_in(const md_attributes_t* attributes, str_t group, str_t leaf) {
    ASSERT(attributes);
    group = attr_prefix_trim(group);
    if (str_empty(group)) {
        return md_attributes_find(attributes, leaf);
    }
    char buf[512];
    if (group.len + 1 + leaf.len > sizeof(buf)) {
        MD_LOG_ERROR("Attribute path '" STR_FMT "/" STR_FMT "' is too long", STR_ARG(group), STR_ARG(leaf));
        return NULL;
    }
    MEMCPY(buf, group.ptr, group.len);
    buf[group.len] = '/';
    MEMCPY(buf + group.len + 1, leaf.ptr, leaf.len);
    return md_attributes_find(attributes, (str_t){buf, group.len + 1 + leaf.len});
}

const md_attribute_t* md_attributes_sibling(const md_attributes_t* attributes, const md_attribute_t* attr, str_t leaf) {
    ASSERT(attributes);
    if (!attr) {
        return NULL;
    }
    return md_attributes_find_in(attributes, md_attribute_group(attr), leaf);
}

const md_attribute_t* md_attributes_axis(const md_attributes_t* attributes, const md_attribute_t* attr) {
    ASSERT(attributes);
    ASSERT(attr);

    // An alias is a second name for a datum, and the datum's axis is decided by where its owner
    // lives. Searching from the alias' own path would pair it with whatever "time" happens to sit
    // above the new name.
    const md_attribute_t* owner = attr;
    if (md_attribute_is_alias(attr)) {
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
    return md_attribute_extract_f64(out, 1, axis, md_attribute_slice_1((uint32_t)i), unit) == 1;
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

static void attr_io_close_slot(md_attribute_io_t* io, size_t i) {
    if (io->slot[i].hash) {
        md_file_close(&io->slot[i].file);
        str_free(io->slot[i].path, io->alloc);
        MEMSET(&io->slot[i], 0, sizeof(io->slot[i]));
    }
}

void md_attribute_io_close_all(md_attribute_io_t* io) {
    for (size_t i = 0; i < MD_ATTRIBUTE_IO_MAX_FILES; ++i) {
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
    for (size_t i = 0; i < MD_ATTRIBUTE_IO_MAX_FILES; ++i) {
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

