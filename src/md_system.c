#include <md_system.h>
#include "md_attributes_internal.h"
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
    md_array_free(sys->atom.formal_charge, alloc);
    md_array_free(sys->atom.hydrogen_count, alloc);

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

    {
        const md_attribute_t* axis = str_empty(run) ? NULL : md_attributes_find_in(attributes, run, STR_LIT("time"));
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

        const md_attribute_t* attr = md_attributes_find_in(attributes, run, paths[i]);
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
    md_attribute_io_close_all(&ex->io);
    md_arena_allocator_destroy(ex->arena);
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
            if (md_attribute_extract_io_f32((float*)out->xyz, N * 3, attr, slice, md_unit_angstrom(), &ex->io) != N * 3) {
                result = false;
                break;
            }
            out->num_atoms = N;
            break;
        }
        case EXTRACT_CELL: {
            float box[3][3];
            if (md_attribute_extract_io_f32(&box[0][0], 9, attr, slice, md_unit_angstrom(), &ex->io) != 9) {
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
            md_attribute_slice_format(&fmt, attr, slice);
            const size_t bytes = md_attribute_byte_size(&fmt);
            void* data = md_temp_alloc(temp, MAX(bytes, 1));
            if (!data || !md_attribute_read_stored(data, attr, slice, &ex->io)) {
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
        .unit = d->time_unit, .label = STR_INIT("Time"),
        .data = d->time, .byte_size = F * sizeof(double)});

    if (ok && d->step) {
        ok = md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), run, STR_LIT("step")), .format = series_i64, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
            .unit = md_unit_none(), .label = STR_INIT("Step"),
            .description = STR_INIT("The file's own notion of where a frame sits in the run, not the frame ordinal"),
            .data = d->step, .byte_size = F * sizeof(int64_t)});
    }

    if (ok && (d->unitcell || d->unitcell_virt)) {
        ok = md_attributes_replace(attributes, &(md_attribute_desc_t){
            .path = md_run_path(buf, sizeof(buf), run, STR_LIT("unitcell")), .format = cell_format,
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = md_unit_angstrom(), .label = STR_INIT("Unit Cell"),
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
        .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = md_unit_angstrom(), .label = STR_INIT("Position"),
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

    const md_attribute_t* path   = md_attributes_find_in(attributes, out->run, STR_LIT("source/path"));
    const md_attribute_t* offset = md_attributes_find_in(attributes, out->run, STR_LIT("source/offset"));
    const md_attribute_t* size   = md_attributes_find_in(attributes, out->run, STR_LIT("source/size"));
    const int64_t* offsets = (const int64_t*)md_attribute_view(offset, MD_ATTRIBUTE_TYPE_I64, 1, 1);
    const int64_t* sizes   = (const int64_t*)md_attribute_view(size,   MD_ATTRIBUTE_TYPE_I64, 1, 1);
    if (!path || path->format.type != MD_ATTRIBUTE_TYPE_STR || !offsets || !sizes ||
        offset->format.shape[0] != size->format.shape[0]) {
        MD_LOG_ERROR("The run '" STR_FMT "' has lost its source attributes", STR_ARG(out->run));
        return false;
    }
    out->path       = md_attribute_str(attributes, path, 0);
    out->offset     = offsets;
    out->size       = sizes;
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
    const md_attribute_t* run_axis = str_empty(run) ? NULL : md_attributes_find_in(attributes, run, STR_LIT("time"));
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
            .path = STR_INIT("time"), .format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)R } },
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
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = time_unit, .label = STR_INIT("Time"),
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
