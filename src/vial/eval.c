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
#include <core/md_os.h>
#include <core/md_common.h>

#include <math.h>
#include <stdio.h>
#include <string.h>

// A worker's own reader of the run, and view of the system
typedef struct eval_worker_t {
    struct md_vial_eval_t*      e;
    md_allocator_i*             arena;      // its own: workers allocate at once
    vial_sys_t                  vs;
    md_system_extract_t*        ex;
    md_system_state_t           state;
} eval_worker_t;

// An output taken over from an earlier evaluation: its attribute as that one published it, renamed
typedef struct eval_held_t {
    md_attribute_desc_t desc;
    char                path[160];
} eval_held_t;

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
    volatile uint64_t           done;       // md_vial_eval_run finished every frame: its values may be read
    eval_held_t*                held;       // md_array: outputs taken over, published as they were
    md_attribute_id_t*          streamed;   // md_array: outputs written as frames finish (decided 31 October)
    // Frames on several threads: one interpreter and one reader each
    vial_interp_t*              workers;
    eval_worker_t*              readers;
    size_t                      num_readers;
    size_t                      num_workers;
    md_vial_parallel_for_t      parallel_for;
    void*                       parallel_user;
};

static bool read_frame(md_vial_eval_t* e, md_system_extract_t** ex, md_system_state_t* state, md_allocator_i* alloc, int64_t index, vial_frame_t* out) {
    if (index < 0 || index >= e->num_frames) return false;
    if (!*ex) {
        const str_t paths[] = {STR_LIT("atom/position"), STR_LIT("unitcell")};
        *ex = md_system_extract_begin(e->sys, e->run, paths, 2, alloc);
        if (!*ex) return false;
    }
    if (!md_system_extract_frame(*ex, index, state)) return false;
    *out = (vial_frame_t){.xyz = (const float*)state->xyz, .num_atoms = state->num_atoms, .cell = state->unitcell, .index = index,
                          .time = e->times ? e->times[index] : (double)index};
    return true;
}

static bool load_frame(void* user, int64_t index, vial_frame_t* out) {
    md_vial_eval_t* e = user;
    return read_frame(e, &e->ex, &e->state, e->arena, index, out);
}

static bool worker_load_frame(void* user, int64_t index, vial_frame_t* out) {
    eval_worker_t* w = user;
    return read_frame(w->e, &w->ex, &w->state, w->arena, index, out);
}

// A reader belongs to one thread at a time, and a worker moves between threads: it ends its reader when a thread is
// done with it, and the next thread begins one
static void worker_release(void* user) {
    eval_worker_t* w = user;
    if (w->ex) { md_system_extract_end(w->ex); w->ex = NULL; }
}

static void main_release(void* user) {
    md_vial_eval_t* e = user;
    if (e->ex) { md_system_extract_end(e->ex); e->ex = NULL; }
}

// Without the caller's pool: num_workers - 1 threads of our own and the calling one, taking indices in turn
typedef struct own_pool_t {
    size_t              count;
    void              (*fn)(void* arg, size_t i);
    void*               arg;
    volatile uint64_t   next;
} own_pool_t;

static void own_pool_work(void* data) {
    own_pool_t* p = data;
    for (;;) {
        const uint64_t i = md_atomic_fetch_add64(&p->next, 1);
        if (i >= p->count) break;
        p->fn(p->arg, (size_t)i);
    }
}

static void own_parallel_for(void* user, size_t count, void (*fn)(void* arg, size_t i), void* arg) {
    const md_vial_eval_t* e = user;
    own_pool_t pool = {.count = count, .fn = fn, .arg = arg};
    md_thread_t* threads[64];
    size_t nt = 0;
    for (size_t t = 1; t < e->num_workers && t < 64 && t < count; ++t) {
        md_thread_t* th = md_thread_create(own_pool_work, &pool);
        if (th) threads[nt++] = th;
    }
    own_pool_work(&pool);
    for (size_t t = 0; t < nt; ++t) md_thread_join(threads[t]);
}

// The output of prev that out is taken over from: one computed alike, the one of the same name first. NULL for none.
static const md_vial_output_t* held_from(const md_vial_eval_t* prev, const md_vial_output_t* out) {
    if (!prev || !md_atomic_load64(&prev->done) || !out->fingerprint || out->selection) return NULL;
    const md_vial_output_t* po;
    const size_t np = md_vial_outputs(prev->prog, &po);
    const md_vial_output_t* any = NULL;
    for (size_t i = 0; i < np; ++i) {
        if (po[i].fingerprint != out->fingerprint || po[i].selection) continue;
        if (str_eq(po[i].name, out->name)) return &po[i];
        if (!any) any = &po[i];
    }
    return any;
}

bool md_vial_eval_plan(md_vial_plan_t* out, size_t cap, const md_vial_eval_t* prev, const md_vial_program_t* prog, size_t num_frames) {
    const md_vial_output_t* o;
    const size_t n = md_vial_outputs(prog, &o);
    if (!prog || cap < n) return false;
    const md_vial_output_t* po = NULL;
    const size_t np = prev ? md_vial_outputs(prev->prog, &po) : 0;
    const bool same_run = prev && (size_t)prev->num_frames == num_frames;
    for (size_t i = 0; i < n; ++i) {
        const md_vial_output_t* h = same_run ? held_from(prev, &o[i]) : NULL;
        md_vial_plan_t p = {.status = MD_VIAL_PLAN_NEW};
        if (h) {
            p.status = MD_VIAL_PLAN_HELD;
            p.from = h->name;
            p.num_frames = (uint32_t)num_frames;
        } else {
            for (size_t k = 0; k < np; ++k) if (str_eq(po[k].name, o[i].name)) p.status = MD_VIAL_PLAN_CHANGED;
            p.num_pending = o[i].stage == MD_VIAL_STAGE_STATIC ? 0 : (uint32_t)num_frames;
        }
        out[i] = p;
    }
    return true;
}

md_vial_eval_t* md_vial_eval_create(const md_vial_program_t* prog, const md_vial_eval_desc_t* desc, md_allocator_i* alloc) {
    return md_vial_eval_create_from(NULL, prog, desc, alloc);
}

md_vial_eval_t* md_vial_eval_create_from(const md_vial_eval_t* prev, const md_vial_program_t* prog, const md_vial_eval_desc_t* desc, md_allocator_i* alloc) {
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
    // What an earlier evaluation of the same run holds: those outputs are published as it published them, and only
    // the others are computed
    if (prev && prev->num_frames == e->num_frames && str_eq(prev->run, e->run)) {
        const md_vial_output_t* o;
        const size_t n = md_vial_outputs(prog, &o);
        vial_ssa_t* roots = NULL;
        for (size_t i = 0; i < n; ++i) {
            const md_vial_output_t* h = held_from(prev, &o[i]);
            const md_attribute_t* at = NULL;
            if (h) {
                char from[160];
                snprintf(from, sizeof(from), "script/%.*s", (int)h->name.len, h->name.ptr);
                at = md_attributes_find(&prev->attributes, str_from_cstr(from));
            }
            if (at && at->data) {
                eval_held_t hd = {0};
                snprintf(hd.path, sizeof(hd.path), "script/%.*s", (int)o[i].name.len, o[i].name.ptr);
                const size_t bytes = md_attribute_byte_size(&at->format);
                void* data = md_alloc(arena, bytes ? bytes : 1);
                memcpy(data, at->data, bytes);
                hd.desc = (md_attribute_desc_t){.format = at->format, .flags = at->flags, .unit = at->unit,
                                                .label = str_copy(o[i].name, arena), .data = data, .byte_size = bytes};
                md_array_push(e->held, hd, arena);
                continue;
            }
        }
        // What is live: every output but the hidden ones and those taken over
        for (size_t k = 0; md_array_size(e->held) && k < md_array_size(ir->outputs); ++k) {
            const str_t nm = ir->outputs[k].name;
            if (!nm.len || nm.ptr[0] == '_') continue;
            bool held = false;
            for (size_t h = 0; h < md_array_size(e->held) && !held; ++h) held = str_eq(e->held[h].desc.label, nm);
            if (!held) md_array_push(roots, ir->outputs[k].value, arena);
        }
        if (md_array_size(e->held)) vial_interp_only(&e->in, roots, md_array_size(roots));
    }
    e->in.max_elements = desc->max_elements;
    e->in.load_frame = load_frame;
    e->in.load_user = e;
    e->in.load_release = main_release;
    e->in.stop = &e->stop;
    // Frames on several threads: an interpreter, a reader and a view of the system per thread
    if (desc->num_threads > 1 && e->num_frames > 1) {
        const size_t nw = MIN(desc->num_threads, (size_t)64);
        e->workers = md_alloc(arena, nw * sizeof(vial_interp_t));
        e->readers = md_alloc(arena, nw * sizeof(eval_worker_t));
        memset(e->workers, 0, nw * sizeof(vial_interp_t));
        memset(e->readers, 0, nw * sizeof(eval_worker_t));
        for (size_t w = 0; w < nw; ++w) {
            eval_worker_t* r = &e->readers[w];
            r->e = e;
            // From the heap, which threads may share, not from alloc, which may be an arena
            r->arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
            if (!r->arena) { md_vial_eval_free(e); return NULL; }
            e->num_readers = w + 1;
            r->state = (md_system_state_t){.alloc = r->arena};
            md_system_state_init(&r->state, desc->sys->atom.count);
            vial_sys_init(&r->vs, desc->sys, r->arena);
            if (!vial_interp_init(&e->workers[w], ir, &r->vs, desc->sys, md_get_heap_allocator())) { md_vial_eval_free(e); return NULL; }
            e->workers[w].load_frame = worker_load_frame;
            e->workers[w].load_user = r;
            e->workers[w].load_release = worker_release;
            e->num_workers = w + 1;
        }
        e->parallel_for = desc->parallel_for ? desc->parallel_for : own_parallel_for;
        e->parallel_user = desc->parallel_for ? desc->parallel_user : e;
        e->in.workers = e->workers;
        e->in.num_workers = e->num_workers;
        e->in.parallel_for = e->parallel_for;
        e->in.parallel_user = e->parallel_user;
    }
    return e;
}

void md_vial_eval_free(md_vial_eval_t* e) {
    if (!e) return;
    for (size_t w = 0; w < e->num_workers; ++w) {
        if (e->readers[w].ex) md_system_extract_end(e->readers[w].ex);
        vial_interp_free(&e->workers[w]);
    }
    for (size_t w = 0; w < e->num_readers; ++w) md_arena_allocator_destroy(e->readers[w].arena);
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

// Before the frames: the frame axis, what was taken over, and the per-frame outputs whose shape the type knows, each
// with its attribute made now, its frames not computed yet NaN (0 for ints), to be written as frames finish
static void prepare(md_vial_eval_t* e) {
    md_attributes_t* attrs = &e->attributes;
    const int64_t nf = e->num_frames;
    if (nf > 0 && !md_attributes_find(attrs, STR_LIT("time"))) {
        const md_attribute_desc_t axis = {
            .path = STR_LIT("time"), .format = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = {(uint32_t)nf}},
            .flags = MD_ATTRIBUTE_FLAG_TEMPORAL, .unit = md_unit_none(), .label = STR_LIT("Frame"),
        };
        const md_attribute_id_t id = md_attributes_create(attrs, &axis);
        double* ord = id ? md_attributes_data(attrs, id, MD_ATTRIBUTE_TYPE_F64) : NULL;
        for (int64_t i = 0; ord && i < nf; ++i) ord[i] = (double)i;
    }
    for (size_t i = 0; i < md_array_size(e->held); ++i) {
        md_attribute_desc_t d = e->held[i].desc;
        d.path = str_from_cstr(e->held[i].path);
        if (!md_attributes_find(attrs, d.path)) md_attributes_create(attrs, &d);
    }
    vial_interp_t* in = &e->in;
    vial_interp_classify(in);
    for (int s = 0; s < 8; ++s) {
        if (!in->frame_done[s]) in->frame_done[s] = md_alloc(e->arena, nf > 0 ? (size_t)nf : 1);
        memset(in->frame_done[s], 0, nf > 0 ? (size_t)nf : 1);
    }
    const vial_ir_t* ir = e->ir;
    const size_t nv = md_array_size(ir->values);
    if (!in->stream) {
        in->stream = md_alloc(e->arena, (nv ? nv : 1) * sizeof(uint8_t*));
        in->stream_count = md_alloc(e->arena, (nv ? nv : 1) * sizeof(uint32_t));
    }
    memset(in->stream, 0, (nv ? nv : 1) * sizeof(uint8_t*));
    md_array_shrink(e->streamed, 0);
    if (nf <= 0) return;
    for (size_t k = 0; k < md_array_size(ir->outputs); ++k) {
        const vial_ir_output_t* o = &ir->outputs[k];
        if (!o->name.len || o->name.ptr[0] == '_') continue;
        const uint32_t op = ir->values[o->value].op;
        if (op == VIAL_IR_NO_OP || !in->live[op] || in->op_class[op] != VIAL_OP_PER_FRAME) continue;
        const vial_type_t* t = &ir->values[o->value].type;
        const md_attribute_type_t at = attr_type(t->dtype);
        if (at == MD_ATTRIBUTE_TYPE_NONE || t->ndim > 2) continue;
        bool known = true;
        uint32_t inner = 1;
        for (uint8_t d = 0; d < t->ndim; ++d) { known = known && t->dims[d].kind == VIAL_EXT_STATIC; inner *= t->dims[d].n; }
        if (!known) continue;
        char path[160];
        snprintf(path, sizeof(path), "script/%.*s", (int)o->name.len, o->name.ptr);
        md_attribute_desc_t desc = {.path = str_from_cstr(path), .unit = t->unit_kind == VIAL_UNIT_KNOWN ? t->unit : md_unit_none(),
                                    .label = o->name, .format = {.type = at, .components = 1}, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL};
        desc.format.shape[desc.format.rank++] = (uint32_t)nf;
        for (uint8_t d = 0; d < t->ndim; ++d) desc.format.shape[desc.format.rank++] = t->dims[d].n;
        const md_attribute_id_t id = md_attributes_replace(attrs, &desc);
        uint8_t* dst = id ? md_attributes_data(attrs, id, at) : NULL;
        if (!dst) continue;
        const size_t n = (size_t)nf * inner;
        if (t->dtype == VIAL_DT_FLOAT) for (size_t i = 0; i < n; ++i) ((float*)dst)[i] = NAN;
        else if (t->dtype == VIAL_DT_FLOAT64) for (size_t i = 0; i < n; ++i) ((double*)dst)[i] = NAN;
        else memset(dst, 0, n * vial_dtype_size(t->dtype));
        in->stream[o->value] = dst;
        in->stream_count[o->value] = inner;
        md_array_push(e->streamed, id, e->arena);
    }
}

static void publish(md_vial_eval_t* e) {
    md_attributes_t* attrs = &e->attributes;
    const int64_t nf = e->num_frames;
    // What was written as frames finished is complete: a new version, read once
    for (size_t i = 0; i < md_array_size(e->streamed); ++i) md_attributes_touch(attrs, e->streamed[i]);
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
        if (e->in.stream && e->in.stream[o->value]) continue;            // written as frames finished
        bool held = false;
        for (size_t h = 0; h < md_array_size(e->held) && !held; ++h) held = str_eq(str_from_cstr(e->held[h].path), str_from_cstr(path));
        if (held) continue;
        md_attribute_desc_t desc = {.path = str_from_cstr(path), .unit = t->unit_kind == VIAL_UNIT_KNOWN ? t->unit : md_unit_none(),
                                    .label = o->name, .format = {.type = at, .components = 1}};
        const bool temporal = (per_frame || series_whole) && nf > 0;
        if (temporal) { desc.flags = MD_ATTRIBUTE_FLAG_TEMPORAL; desc.format.shape[desc.format.rank++] = (uint32_t)nf; }
        for (uint8_t d = 0; d < one.ndim; ++d) desc.format.shape[desc.format.rank++] = one.dims[d];
        const md_attribute_id_t id = md_attributes_replace(attrs, &desc);
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
    e->done = 0;    // a run again: its values are not read until it has finished
    prepare(e);
    const bool ok = vial_interp_run(&e->in, e->num_frames > 0 ? e->num_frames : 0);
    md_array_shrink(e->diags, 0);
    for (size_t i = 0; i < md_array_size(e->in.errors); ++i) {
        const vial_rt_error_t* r = &e->in.errors[i];
        md_vial_diagnostic_t d = {.severity = MD_VIAL_SEVERITY_RUNTIME, .span = r->span, .message = str_copy(str_from_cstr(r->message), e->arena)};
        memcpy(d.code, r->code, 6);
        md_array_push(e->diags, d, e->arena);
    }
    if (!e->in.interrupted) {
        publish(e);
        md_atomic_fetch_add64(&e->done, 1);
    }
    return ok && !e->in.interrupted;
}

size_t md_vial_eval_output_frames(const md_vial_eval_t* e, str_t name, uint8_t* out, size_t cap) {
    if (!e) return 0;
    const size_t nf = e->num_frames > 0 ? (size_t)e->num_frames : 0;
    if (!out || cap < nf) return nf;
    memset(out, 0, nf);
    char path[160];
    snprintf(path, sizeof(path), "script/%.*s", (int)name.len, name.ptr);
    for (size_t h = 0; h < md_array_size(e->held); ++h) {
        if (str_eq(str_from_cstr(e->held[h].path), str_from_cstr(path))) { memset(out, 1, nf); return nf; }
    }
    const vial_ir_t* ir = e->ir;
    for (size_t k = 0; k < md_array_size(ir->outputs); ++k) {
        if (!str_eq(ir->outputs[k].name, name)) continue;
        const uint32_t op = ir->values[ir->outputs[k].value].op;
        const vial_interp_t* in = &e->in;
        if (md_atomic_load64(&e->done)) { memset(out, 1, nf); return nf; }
        if (op == VIAL_IR_NO_OP || !in->op_class || in->op_class[op] != VIAL_OP_PER_FRAME) return nf;
        const uint8_t s = in->op_sweep[op];
        if (s < 8 && in->frame_done[s]) memcpy(out, in->frame_done[s], nf);
        return nf;
    }
    return nf;
}

size_t md_vial_eval_aggregate(const md_vial_eval_t* e, str_t name, const uint8_t* frames, float* out, size_t cap) {
    if (!e || !frames || !md_atomic_load64(&e->done)) return 0;
    const vial_ir_t* ir = e->ir;
    vial_ssa_t v = 0;
    for (size_t k = 0; k < md_array_size(ir->outputs) && !v; ++k) if (str_eq(ir->outputs[k].name, name)) v = ir->outputs[k].value;
    if (!v) return 0;
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    vial_val_t r;
    size_t n = 0;
    if (vial_interp_over_frames(&e->in, v, frames, &r, arena) && (r.dtype == VIAL_DT_FLOAT || r.dtype == VIAL_DT_FLOAT64 || r.dtype == VIAL_DT_INT) && r.data) {
        n = vial_val_count(&r);
        for (size_t i = 0; out && n <= cap && i < n; ++i) {
            const double x = r.dtype == VIAL_DT_FLOAT ? ((const float*)r.data)[i] : r.dtype == VIAL_DT_FLOAT64 ? ((const double*)r.data)[i] : ((const int32_t*)r.data)[i];
            out[i] = vial_val_is_missing(&r, i) ? NAN : (float)x;
        }
    }
    md_arena_allocator_destroy(arena);
    return n;
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
    // A value over the trajectory needs every frame: the evaluation's, once it has finished and if it computed it
    if (ir->values[v].type.stage == VIAL_STAGE_TRAJECTORY) {
        if (!md_atomic_load64(&e->done)) return false;
        const vial_ssa_t rv = md_vial__run_value(e->prog, v);
        const vial_val_t* val = rv ? vial_interp_at(&e->in, rv, frame) : NULL;
        if (!val) return false;
        // A series kept whole, as smooth gives: its frame
        if (val->flags & VIAL_VAL_FRAMES) {
            if (frame < 0 || frame >= (int64_t)val->dims[0]) return false;
            vial_val_t one = *val;
            one.flags &= ~VIAL_VAL_FRAMES;
            one.ndim -= 1;
            memmove(one.dims, one.dims + 1, one.ndim * sizeof(uint32_t));
            size_t inner = 1;
            for (uint8_t q = 0; q < one.ndim; ++q) inner *= one.dims[q];
            if (val->missing) return false;
            one.data = (uint8_t*)val->data + (size_t)frame * inner * vial_dtype_size(val->dtype);
            vial_val_print(buf, cap, &one, 6);
            return true;
        }
        vial_val_print(buf, cap, val, 6);
        return true;
    }
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
