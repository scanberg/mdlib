// The public evaluation API (md_vial.h): a program evaluated over one run with the interpreter, its published values
// made attributes.

#include "interp.h"
#include "check.h"

#include <md_system.h>
#include <md_attributes.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_intrinsics.h>

#include <stdio.h>
#include <string.h>

struct md_vial_eval_t {
    md_allocator_i*             arena;
    const md_vial_program_t*    prog;
    const vial_ir_t*            ir;
    const md_system_t*          sys;
    str_t                       run;
    vial_sys_t                  vs;
    vial_interp_t               in;
    md_system_extract_t*        ex;
    md_system_state_t           state;
    double*                     times;      // the run's time axis, in ps where it has a unit
    int64_t                     num_frames;
    volatile uint64_t           stop;
    md_attributes_t             attributes;
    md_vial_diagnostic_t*       diags;      // md_array
};

static bool load_frame(void* user, int64_t index, vial_frame_t* out) {
    md_vial_eval_t* e = user;
    if (index < 0 || index >= e->num_frames) return false;
    if (!e->ex) {
        const str_t paths[] = {STR_LIT("atom/position"), STR_LIT("unitcell")};
        e->ex = md_system_extract_begin(e->sys, e->run, paths, 2, e->arena);
        if (!e->ex) return false;
    }
    if (!md_system_extract_frame(e->ex, index, &e->state)) return false;
    *out = (vial_frame_t){.xyz = (const float*)e->state.xyz, .num_atoms = e->state.num_atoms, .cell = e->state.unitcell, .index = index,
                          .time = e->times ? e->times[index] : (double)index};
    return true;
}

md_vial_eval_t* md_vial_eval_create(const md_vial_program_t* prog, const md_vial_eval_desc_t* desc, md_allocator_i* alloc) {
    if (!prog || !desc || !desc->sys) return NULL;
    const vial_ir_t* ir = md_vial__run_ir(prog);
    if (!ir) return NULL;
    const md_vial_diagnostic_t* d;
    const size_t nd = md_vial_diagnostics(prog, &d);
    for (size_t i = 0; i < nd; ++i) if (d[i].severity == MD_VIAL_SEVERITY_ERROR) return NULL;
    md_allocator_i* arena = md_arena_allocator_create(alloc, MEGABYTES(1));
    if (!arena) return NULL;
    md_vial_eval_t* e = md_alloc(arena, sizeof(md_vial_eval_t));
    memset(e, 0, sizeof(*e));
    e->arena = arena;
    e->prog = prog;
    e->ir = ir;
    e->sys = desc->sys;
    e->run = str_copy(desc->run, arena);
    e->attributes.alloc = arena;
    // The run's frames, and their times
    char path[256];
    snprintf(path, sizeof(path), "%.*s/time", (int)e->run.len, e->run.ptr);
    const md_attribute_t* time = md_attributes_find(&desc->sys->attributes, str_from_cstr(path));
    if (time && time->format.rank >= 1) {
        e->num_frames = time->format.shape[0];
        e->times = md_alloc(arena, (e->num_frames ? e->num_frames : 1) * sizeof(double));
        double f;
        const md_unit_t unit = !md_unit_is_none(time->unit) && md_unit_conversion_factor(&f, time->unit, md_unit_picosecond()) ? md_unit_picosecond() : time->unit;
        if (md_attribute_extract_f64(e->times, (size_t)e->num_frames, time, md_attribute_slice_all(), unit) != (size_t)e->num_frames) {
            for (int64_t i = 0; i < e->num_frames; ++i) e->times[i] = (double)i;
        }
    }
    e->state = (md_system_state_t){.alloc = arena};
    md_system_state_init(&e->state, desc->sys->atom.count);
    vial_sys_init(&e->vs, desc->sys, arena);
    if (!vial_interp_init(&e->in, ir, &e->vs, desc->sys, alloc)) { md_arena_allocator_destroy(arena); return NULL; }
    e->in.max_elements = desc->max_elements;
    e->in.load_frame = load_frame;
    e->in.load_user = e;
    e->in.stop = &e->stop;
    return e;
}

void md_vial_eval_free(md_vial_eval_t* e) {
    if (!e) return;
    if (e->ex) md_system_extract_end(e->ex);
    vial_interp_free(&e->in);
    md_arena_allocator_destroy(e->arena);
}

double md_vial_eval_progress(const md_vial_eval_t* e) {
    if (!e) return 0.0;
    const uint64_t total = md_atomic_load64(&e->in.frames_total), done = md_atomic_load64(&e->in.frames_done);
    return total ? (double)done / (double)total : 0.0;
}

void md_vial_eval_interrupt(md_vial_eval_t* e) {
    if (e) md_atomic_fetch_add64(&e->stop, 1);
}

size_t md_vial_eval_diagnostics(const md_vial_eval_t* e, const md_vial_diagnostic_t** out) {
    if (!e) { if (out) *out = NULL; return 0; }
    if (out) *out = e->diags;
    return md_array_size(e->diags);
}

const struct md_attributes_t* md_vial_eval_attributes(const md_vial_eval_t* e) { return e ? &e->attributes : NULL; }

// ### PUBLISHING ###

static md_attribute_type_t attr_type(uint8_t dtype) {
    switch (dtype) {
    case VIAL_DT_BOOL:    return MD_ATTRIBUTE_TYPE_U8;
    case VIAL_DT_INT:     return MD_ATTRIBUTE_TYPE_I32;
    case VIAL_DT_FLOAT:   return MD_ATTRIBUTE_TYPE_F32;
    case VIAL_DT_FLOAT64: return MD_ATTRIBUTE_TYPE_F64;
    default:              return MD_ATTRIBUTE_TYPE_NONE;
    }
}

static void publish(md_vial_eval_t* e) {
    md_attributes_t* attrs = &e->attributes;
    const int64_t nf = e->num_frames;
    if (nf > 0) {
        const md_attribute_desc_t axis = {
            .path = STR_LIT("time"), .format = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = {(uint32_t)nf}},
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = md_unit_none(), .label = STR_LIT("Frame"),
        };
        const md_attribute_id_t id = md_attributes_create(attrs, &axis);
        double* ord = id ? md_attributes_data(attrs, id, MD_ATTRIBUTE_TYPE_F64) : NULL;
        for (int64_t i = 0; ord && i < nf; ++i) ord[i] = (double)i;
    }
    const vial_ir_t* ir = e->ir;
    for (size_t k = 0; k < md_array_size(ir->outputs); ++k) {
        const vial_ir_output_t* o = &ir->outputs[k];
        if (!o->name.len || o->name.ptr[0] == '_') continue;
        const vial_type_t* t = &ir->values[o->value].type;
        const md_attribute_type_t at = attr_type(t->dtype);
        if (at == MD_ATTRIBUTE_TYPE_NONE || t->ndim > 2) continue;
        const uint32_t op = ir->values[o->value].op;
        const bool per_frame = op != VIAL_IR_NO_OP && e->in.op_class && e->in.op_class[op] == VIAL_OP_PER_FRAME;
        const vial_val_t* first = vial_interp_at(&e->in, o->value, 0);
        if (!first) continue;
        vial_val_t one = *first;
        const bool series_whole = (first->flags & VIAL_VAL_FRAMES) != 0;    // smooth: one value holds every frame
        if (series_whole) { one.ndim -= 1; memmove(one.dims, one.dims + 1, one.ndim * sizeof(uint32_t)); }
        char path[160];
        snprintf(path, sizeof(path), "script/%.*s", (int)o->name.len, o->name.ptr);
        md_attribute_desc_t desc = {.path = str_from_cstr(path), .unit = t->unit_kind == VIAL_UNIT_KNOWN ? t->unit : md_unit_none(),
                                    .label = o->name, .format = {.type = at, .components = 1}};
        const bool temporal = (per_frame || series_whole) && nf > 0;
        if (temporal) { desc.flags = MD_ATTRIBUTE_FLAG_TEMPORAL; desc.format.shape[desc.format.rank++] = (uint32_t)nf; }
        for (uint8_t d = 0; d < one.ndim; ++d) desc.format.shape[desc.format.rank++] = one.dims[d];
        const md_attribute_id_t id = md_attributes_create(attrs, &desc);
        uint8_t* dst = id ? md_attributes_data(attrs, id, at) : NULL;
        if (!dst) continue;
        const size_t inner = vial_val_count(&one), es = vial_dtype_size(t->dtype);
        // A missing value is published as NaN for floats and 0 otherwise; the series shows a gap there
        for (int64_t f = 0; f < (temporal ? nf : 1); ++f) {
            const vial_val_t* v = series_whole ? first : (temporal ? vial_interp_at(&e->in, o->value, f) : first);
            uint8_t* row = dst + (size_t)f * inner * es;
            if (!v || !v->data || vial_val_count(v) < (series_whole ? (size_t)(f + 1) * inner : inner)) {
                if (t->dtype == VIAL_DT_FLOAT) for (size_t i = 0; i < inner; ++i) ((float*)row)[i] = NAN;
                continue;
            }
            const size_t base = series_whole ? (size_t)f * inner : 0;
            memcpy(row, (const uint8_t*)v->data + base * es, inner * es);
            for (size_t i = 0; i < inner; ++i) {
                if (!vial_val_is_missing(v, base + i)) continue;
                if (t->dtype == VIAL_DT_FLOAT) ((float*)row)[i] = NAN;
                else if (t->dtype == VIAL_DT_FLOAT64) ((double*)row)[i] = NAN;
            }
        }
    }
}

bool md_vial_eval_run(md_vial_eval_t* e) {
    if (!e) return false;
    const bool ok = vial_interp_run(&e->in, e->num_frames > 0 ? e->num_frames : 0);
    md_array_shrink(e->diags, 0);
    for (size_t i = 0; i < md_array_size(e->in.errors); ++i) {
        const vial_rt_error_t* r = &e->in.errors[i];
        md_vial_diagnostic_t d = {.severity = MD_VIAL_SEVERITY_RUNTIME, .span = r->span, .message = str_copy(str_from_cstr(r->message), e->arena)};
        memcpy(d.code, r->code, 6);
        md_array_push(e->diags, d, e->arena);
    }
    if (!e->in.interrupted) publish(e);
    return ok && !e->in.interrupted;
}

// ### HOVER ###
// The hovered expression's value at a frame: the script as lowered, before inlining, so every expression written is
// there, evaluated for that value alone in an interpreter of its own, with its own reader of the run and its own view
// of the system. It shares nothing that changes with md_vial_eval_run.

typedef struct hover_frames_t {
    const md_system_t*      sys;
    str_t                   run;
    const double*           times;
    int64_t                 num_frames;
    md_system_extract_t*    ex;
    md_system_state_t       state;
    md_allocator_i*         alloc;
} hover_frames_t;

static bool hover_load(void* user, int64_t index, vial_frame_t* out) {
    hover_frames_t* h = user;
    if (index < 0 || index >= h->num_frames) return false;
    if (!h->ex) {
        const str_t paths[] = {STR_LIT("atom/position"), STR_LIT("unitcell")};
        h->ex = md_system_extract_begin(h->sys, h->run, paths, 2, h->alloc);
        if (!h->ex) return false;
    }
    if (!md_system_extract_frame(h->ex, index, &h->state)) return false;
    *out = (vial_frame_t){.xyz = (const float*)h->state.xyz, .num_atoms = h->state.num_atoms, .cell = h->state.unitcell, .index = index,
                          .time = h->times ? h->times[index] : (double)index};
    return true;
}

bool md_vial_eval_hover(const md_vial_eval_t* e, uint32_t offset, int64_t frame, char* buf, size_t cap) {
    if (buf && cap) buf[0] = 0;
    if (!e || !buf || !cap) return false;
    const vial_ir_t* ir = md_vial__ir(e->prog);
    const vial_ssa_t v = md_vial__hover_value(e->prog, offset, NULL);
    if (!ir || !v || v >= md_array_size(ir->values)) return false;
    // A value over the trajectory needs every frame: not one frame's work
    if (ir->values[v].type.stage == VIAL_STAGE_TRAJECTORY) return false;
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    hover_frames_t h = {.sys = e->sys, .run = e->run, .times = e->times, .num_frames = e->num_frames, .alloc = arena,
                        .state = {.alloc = arena}};
    md_system_state_init(&h.state, e->sys->atom.count);
    vial_sys_t vs;
    vial_sys_init(&vs, e->sys, arena);
    vial_interp_t in;
    bool ok = vial_interp_init(&in, ir, &vs, e->sys, md_get_heap_allocator());
    if (ok) {
        vial_interp_only(&in, &v, 1);
        in.max_elements = (size_t)1 << 24;     // hover computes at most 16 million entries
        in.load_frame = hover_load;
        in.load_user = &h;
        // Static values first: at() reads its own frame, through the same reader
        ok = vial_interp_static(&in);
        if (ok && ir->values[v].type.stage != VIAL_STAGE_STATIC) {
            vial_frame_t fr;
            ok = hover_load(&h, frame, &fr) && vial_interp_frame(&in, &fr);
        }
        const vial_val_t* val = ok ? vial_interp_value(&in, v) : NULL;
        ok = val != NULL;
        if (ok) vial_val_print(buf, cap, val, 6);
        vial_interp_free(&in);
    }
    if (h.ex) md_system_extract_end(h.ex);
    md_arena_allocator_destroy(arena);
    return ok;
}
