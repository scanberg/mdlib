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
#include <stdarg.h>
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

// The run's time axis, <run>/time, in ps where it has a unit, and its number of frames; NULL and 0 when it has none
static double* run_times(const md_system_t* sys, str_t run, md_allocator_i* alloc, int64_t* num_frames) {
    *num_frames = 0;
    char path[256];
    snprintf(path, sizeof(path), "%.*s/time", (int)run.len, run.ptr);
    const md_attribute_t* time = md_attributes_find(&sys->attributes, str_from_cstr(path));
    if (!time || time->format.rank < 1) return NULL;
    const int64_t n = time->format.shape[0];
    double* times = md_alloc(alloc, (n ? n : 1) * sizeof(double));
    double f;
    const md_unit_t unit = !md_unit_is_none(time->unit) && md_unit_conversion_factor(&f, time->unit, md_unit_picosecond()) ? md_unit_picosecond() : time->unit;
    if (md_attribute_extract_f64(times, (size_t)n, time, md_attribute_slice_all(), unit) != (size_t)n) {
        for (int64_t i = 0; i < n; ++i) times[i] = (double)i;
    }
    *num_frames = n;
    return times;
}

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
    e->times = run_times(desc->sys, e->run, arena, &e->num_frames);
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

void md_vial_eval_reset_interrupt(md_vial_eval_t* e) {
    if (!e) return;
    e->stop = 0;
    e->in.interrupted = false;
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

// One value of the IR as lowered, alone: its own interpreter, reader and view of the system
typedef struct point_t {
    md_allocator_i* arena;
    hover_frames_t  h;
    vial_sys_t      vs;
    vial_interp_t   in;
    bool            have_in;
} point_t;

static void point_free(point_t* pt) {
    if (pt->have_in) vial_interp_free(&pt->in);
    if (pt->h.ex) md_system_extract_end(pt->h.ex);
    if (pt->arena) md_arena_allocator_destroy(pt->arena);
    memset(pt, 0, sizeof(*pt));
}

// v at frame of the run, or at the coordinates of state when given, which are then that frame's; at most 16 million
// entries. While drawing, what it draws is attached (vial_val_t.drawn).
static const vial_val_t* point_eval(point_t* pt, const md_vial_program_t* prog, const md_system_t* sys, str_t run, const double* times,
                                    int64_t num_frames, int64_t frame, const md_system_state_t* state, vial_ssa_t v, bool draw) {
    memset(pt, 0, sizeof(*pt));
    const vial_ir_t* ir = md_vial__ir(prog);
    pt->arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    pt->h = (hover_frames_t){.sys = sys, .run = run, .times = times, .num_frames = num_frames, .alloc = pt->arena, .state = {.alloc = pt->arena}};
    md_system_state_init(&pt->h.state, sys->atom.count);
    vial_sys_init(&pt->vs, sys, pt->arena);
    if (!vial_interp_init(&pt->in, ir, &pt->vs, sys, md_get_heap_allocator())) return NULL;
    pt->have_in = true;
    vial_interp_t* in = &pt->in;
    vial_interp_only(in, &v, 1);
    in->max_elements = (size_t)1 << 24;     // hover computes at most 16 million entries
    in->load_frame = hover_load;
    in->load_user = &pt->h;
    in->draw = draw;
    // Static values first: at() reads its own frame, through the same reader
    bool ok = vial_interp_static(in);
    if (ok && ir->values[v].type.stage != VIAL_STAGE_STATIC) {
        vial_frame_t fr;
        if (state) {
            if (state->num_atoms != sys->atom.count || !state->xyz) return NULL;
            fr = (vial_frame_t){.xyz = (const float*)state->xyz, .num_atoms = state->num_atoms, .cell = state->unitcell, .index = frame,
                                .time = times && frame >= 0 && frame < num_frames ? times[frame] : (double)frame};
        } else if (!hover_load(&pt->h, frame, &fr)) {
            return NULL;
        }
        ok = vial_interp_frame(in, &fr);
    }
    return ok ? vial_interp_value(in, v) : NULL;
}

// A value over the trajectory, from a finished evaluation, at frame when it is a series kept whole: into one
static bool trajectory_value(const md_vial_eval_t* e, vial_ssa_t v, int64_t frame, vial_val_t* one) {
    if (!e || !md_atomic_load64(&e->done)) return false;
    const vial_ssa_t rv = md_vial__run_value(e->prog, v);
    const vial_val_t* val = rv ? vial_interp_at(&e->in, rv, frame) : NULL;
    if (!val) return false;
    *one = *val;
    if (!(val->flags & VIAL_VAL_FRAMES)) return true;
    // A series kept whole, as smooth gives: its frame
    if (frame < 0 || frame >= (int64_t)val->dims[0] || val->missing) return false;
    one->flags &= ~VIAL_VAL_FRAMES;
    one->ndim -= 1;
    memmove(one->dims, one->dims + 1, one->ndim * sizeof(uint32_t));
    size_t inner = 1;
    for (uint8_t q = 0; q < one->ndim; ++q) inner *= one->dims[q];
    one->data = (uint8_t*)val->data + (size_t)frame * inner * vial_dtype_size(val->dtype);
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
        vial_val_t one;
        if (!trajectory_value(e, v, frame, &one)) return false;
        vial_val_print(buf, cap, &one, 6);
        return true;
    }
    point_t pt;
    const vial_val_t* val = point_eval(&pt, e->prog, e->sys, e->run, e->times, e->num_frames, frame, NULL, v, false);
    if (val) vial_val_print(buf, cap, val, 6);
    point_free(&pt);
    return val != NULL;
}

// ### DRAWING ###

#define DRAW_MAX_PRIMS 10000

typedef struct draw_out_t {
    md_vial_drawing_t*  d;
    md_allocator_i*     alloc;
    int64_t             only;   // the element drawn alone, or -1
    size_t              count;  // primitives and labels so far
} draw_out_t;

static bool draw_room(draw_out_t* o) {
    if (o->count >= DRAW_MAX_PRIMS) { o->d->num_omitted += 1; return false; }
    o->count += 1;
    return true;
}

static void draw_label(draw_out_t* o, const float pos[3], uint32_t element, const vial_prim_t* p) {
    if (!draw_room(o)) return;
    md_vial_label_t l = {.pos = {pos[0], pos[1], pos[2]}, .element = element, .text = str_copy(p->text, o->alloc), .value = p->value, .unit = p->unit};
    md_array_push(o->d->labels, l, o->alloc);
}

static void draw_group(draw_out_t* o, uint32_t element, const uint32_t* atoms, uint32_t n) {
    md_vial_group_t g = {.element = element, .beg = (uint32_t)md_array_size(o->d->atoms)};
    if (n) md_array_push_array(o->d->atoms, atoms, n, o->alloc);
    g.end = (uint32_t)md_array_size(o->d->atoms);
    md_array_push(o->d->groups, g, o->alloc);
}

static void draw_prims(draw_out_t* o, const vial_drawing_t* dr, uint32_t element) {
    if (!dr) return;
    for (uint32_t i = 0; i < dr->count; ++i) {
        const vial_prim_t* p = &dr->prims[i];
        switch (p->kind) {
        case VIAL_PRIM_HIGHLIGHT: draw_group(o, element, p->atoms, p->num_atoms); continue;
        case VIAL_PRIM_LABEL:     draw_label(o, p->p[0], element, p); continue;
        default: break;
        }
        if (!draw_room(o)) continue;
        static const md_vial_prim_kind_t kinds[] = {
            [VIAL_PRIM_MARKER] = MD_VIAL_PRIM_MARKER, [VIAL_PRIM_LINE] = MD_VIAL_PRIM_LINE, [VIAL_PRIM_ARROW] = MD_VIAL_PRIM_ARROW,
            [VIAL_PRIM_SPHERE] = MD_VIAL_PRIM_SPHERE, [VIAL_PRIM_ARC] = MD_VIAL_PRIM_ARC, [VIAL_PRIM_DISK] = MD_VIAL_PRIM_DISK,
        };
        md_vial_prim_t q = {.kind = kinds[p->kind], .element = element, .r = p->r};
        memcpy(q.p, p->p, sizeof(q.p));
        md_array_push(o->d->prims, q, o->alloc);
        if (!p->labelled) continue;
        // A line's label at its middle, an arc's inside the angle
        float at[3];
        for (int k = 0; k < 3; ++k) {
            at[k] = p->kind == VIAL_PRIM_ARC ? p->p[0][k] + 0.25f * ((p->p[1][k] - p->p[0][k]) + (p->p[2][k] - p->p[0][k]))
                                             : 0.5f * (p->p[0][k] + p->p[1][k]);
        }
        draw_label(o, at, element, p);
    }
}

// A value by its type: sets as their particles, points (Å, a last dimension of 3) as markers, anything else nothing
static void draw_by_type(draw_out_t* o, const vial_val_t* v, const vial_type_t* t) {
    const size_t n = v->dtype == VIAL_DT_SET ? (v->ndim ? vial_val_count(v) : 1) : vial_val_count(v);
    if (v->dtype == VIAL_DT_SET && v->sets) {
        o->d->num_elements = n;
        for (size_t i = 0; i < n; ++i) {
            if (o->only >= 0 && (int64_t)i != o->only) continue;
            draw_group(o, (uint32_t)i, vial_sets_items(v->sets, (uint32_t)i), vial_sets_size(v->sets, (uint32_t)i));
        }
        return;
    }
    md_unit_t length;
    vial_unit_name(STR_LIT("angstrom"), &length);
    const bool pts = (v->dtype == VIAL_DT_FLOAT || v->dtype == VIAL_DT_FLOAT64) && v->ndim >= 1 && v->dims[v->ndim - 1] == 3 && v->data &&
                     t->unit_kind == VIAL_UNIT_KNOWN && md_unit_base_equal(t->unit, length);
    o->d->num_elements = pts ? n / 3 : n;
    if (!pts) return;
    for (size_t i = 0; i < n / 3; ++i) {
        if (o->only >= 0 && (int64_t)i != o->only) continue;
        md_vial_prim_t q = {.kind = MD_VIAL_PRIM_MARKER, .element = (uint32_t)i};
        bool ok = true;
        for (int k = 0; k < 3; ++k) {
            const size_t at = 3 * i + (size_t)k;
            q.p[0][k] = v->dtype == VIAL_DT_FLOAT ? ((const float*)v->data)[at] : (float)((const double*)v->data)[at];
            ok = ok && !vial_val_is_missing(v, at) && isfinite(q.p[0][k]);
        }
        if (ok && draw_room(o)) md_array_push(o->d->prims, q, o->alloc);
    }
}

static void draw_value(draw_out_t* o, const vial_val_t* v, const vial_type_t* t) {
    // The drawing attached, or the value itself when it is a drawing
    const vial_val_t* d = v->drawn ? v->drawn : v->dtype == VIAL_DT_DRAWING ? v : NULL;
    if (!d) { draw_by_type(o, v, t); return; }
    const size_t n = vial_val_count(d);
    o->d->num_elements = n;
    if (!d->data) return;
    for (size_t i = 0; i < n; ++i) {
        if (o->only >= 0 && (int64_t)i != o->only) continue;
        if (vial_val_is_missing(d, i)) continue;
        draw_prims(o, ((const vial_drawing_t* const*)d->data)[i], (uint32_t)i);
    }
}

static bool draw_done(md_vial_drawing_t* d, bool ok) {
    d->num_prims = md_array_size(d->prims);
    d->num_labels = md_array_size(d->labels);
    d->num_groups = md_array_size(d->groups);
    return ok;
}

bool md_vial_draw(md_vial_drawing_t* out, const md_vial_draw_desc_t* desc, uint32_t expr, md_allocator_i* alloc) {
    if (!out) return false;
    memset(out, 0, sizeof(*out));
    if (!desc || !desc->prog || !desc->sys || !alloc) return false;
    const vial_ir_t* ir = md_vial__ir(desc->prog);
    const vial_ssa_t v = md_vial__expr_value(desc->prog, expr);
    if (!ir || !v || v >= md_array_size(ir->values)) return false;
    const vial_type_t* t = &ir->values[v].type;
    draw_out_t o = {.d = out, .alloc = alloc, .only = desc->element};
    if (t->stage == VIAL_STAGE_TRAJECTORY) {
        // Over the trajectory: the evaluation's value, drawn by its type
        vial_val_t one;
        if (!trajectory_value(desc->eval, v, desc->frame, &one)) return false;
        one.drawn = NULL;
        draw_value(&o, &one, t);
        return draw_done(out, true);
    }
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    int64_t num_frames = 0;
    const double* times = desc->eval ? desc->eval->times : run_times(desc->sys, desc->run, arena, &num_frames);
    if (desc->eval) num_frames = desc->eval->num_frames;
    point_t pt;
    const vial_val_t* val = point_eval(&pt, desc->prog, desc->sys, desc->run, times, num_frames, desc->frame, desc->state, v, true);
    if (val) draw_value(&o, val, t);
    point_free(&pt);
    md_arena_allocator_destroy(arena);
    return draw_done(out, val != NULL);
}

size_t md_vial_value_text(char* buf, size_t cap, const md_vial_draw_desc_t* desc, uint32_t expr) {
    if (buf && cap) buf[0] = 0;
    if (!buf || !cap || !desc || !desc->prog || !desc->sys) return 0;
    const vial_ir_t* ir = md_vial__ir(desc->prog);
    const vial_ssa_t v = md_vial__expr_value(desc->prog, expr);
    if (!ir || !v || v >= md_array_size(ir->values)) return 0;
    if (ir->values[v].type.stage == VIAL_STAGE_TRAJECTORY) {
        vial_val_t one;
        if (!trajectory_value(desc->eval, v, desc->frame, &one)) return 0;
        return vial_val_print(buf, cap, &one, 6);
    }
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    int64_t num_frames = 0;
    const double* times = desc->eval ? desc->eval->times : run_times(desc->sys, desc->run, arena, &num_frames);
    if (desc->eval) num_frames = desc->eval->num_frames;
    point_t pt;
    const vial_val_t* val = point_eval(&pt, desc->prog, desc->sys, desc->run, times, num_frames, desc->frame, desc->state, v, false);
    const size_t len = val ? vial_val_print(buf, cap, val, 6) : 0;
    point_free(&pt);
    md_arena_allocator_destroy(arena);
    return len;
}

// ### FILTERS ###

#define FILTER_HEAD "_filter = (\n"

// Whether expr reads name: as a word of its own
static bool expr_reads(str_t expr, str_t name) {
    for (size_t i = 0; name.len && i + name.len <= expr.len; ++i) {
        if (strncmp(expr.ptr + i, name.ptr, name.len)) continue;
        const char a = i ? expr.ptr[i - 1] : ' ', z = i + name.len < expr.len ? expr.ptr[i + name.len] : ' ';
        const bool wa = (a >= 'a' && a <= 'z') || (a >= 'A' && a <= 'Z') || (a >= '0' && a <= '9') || a == '_' || a == '.';
        const bool wz = (z >= 'a' && z <= 'z') || (z >= 'A' && z <= 'Z') || (z >= '0' && z <= '9') || z == '_';
        if (!wa && !wz) return true;
    }
    return false;
}

// The expression as a statement of its own, with the stored selections and the program's selections as names: its
// static ones, and the recipes it reads, evaluated at the coordinates given, which make it dynamic
static md_vial_program_t* filter_compile(str_t expr, const md_vial_filter_desc_t* desc, md_allocator_i* arena, bool* reads_recipe) {
    if (reads_recipe) *reads_recipe = false;
    md_vial_named_set_t* sets = NULL;
    for (size_t i = 0; i < desc->num_sets && desc->sets; ++i) md_array_push(sets, desc->sets[i], arena);
    if (desc->prog) {
        const md_vial_output_t* o;
        const size_t n = md_vial_outputs(desc->prog, &o);
        const vial_ir_t* ir = md_vial__ir(desc->prog);
        for (size_t i = 0; ir && i < n; ++i) {
            if (!o[i].selection || (o[i].stage != MD_VIAL_STAGE_STATIC && o[i].stage != MD_VIAL_STAGE_FRAME)) continue;
            bool taken = false;
            for (size_t k = 0; k < md_array_size(sets); ++k) taken |= str_eq(sets[k].name, o[i].name);
            if (taken) continue;
            const vial_ssa_t v = md_vial__hover_value(desc->prog, o[i].span.beg, NULL);
            if (!v || v >= md_array_size(ir->values) || ir->values[v].op == VIAL_IR_NO_OP) continue;
            if (o[i].stage == MD_VIAL_STAGE_FRAME) {
                // A recipe: evaluated at the coordinates, only when read
                if (!desc->state || !expr_reads(expr, o[i].name)) continue;
                point_t pt;
                const vial_val_t* rv = point_eval(&pt, desc->prog, desc->sys, desc->run, NULL, 0, 0, desc->state, v, false);
                if (rv && rv->dtype == VIAL_DT_SET && rv->sets && !rv->ndim && rv->sets->count == 1) {
                    const uint32_t n = vial_sets_size(rv->sets, 0);
                    uint32_t* at = md_alloc(arena, (n ? n : 1) * sizeof(uint32_t));
                    if (n) memcpy(at, vial_sets_items(rv->sets, 0), n * sizeof(uint32_t));
                    md_vial_named_set_t ns = {.name = o[i].name, .atoms = at, .num_atoms = n};
                    md_array_push(sets, ns, arena);
                    if (reads_recipe) *reads_recipe = true;
                }
                point_free(&pt);
                continue;
            }
            const vial_ir_op_t* op = &ir->ops[ir->values[v].op];
            if (op->kind != VIAL_IR_CONST) continue;
            const vial_val_t* c = &ir->consts[op->aux];
            if (c->dtype != VIAL_DT_SET || !c->sets || c->ndim || c->sets->count != 1) continue;
            md_vial_named_set_t ns = {.name = o[i].name, .atoms = vial_sets_items(c->sets, 0), .num_atoms = vial_sets_size(c->sets, 0)};
            md_array_push(sets, ns, arena);
        }
    }
    const size_t head = sizeof(FILTER_HEAD) - 1;
    char* src = md_alloc(arena, head + expr.len + 3);
    memcpy(src, FILTER_HEAD, head);
    if (expr.len) memcpy(src + head, expr.ptr, expr.len);
    memcpy(src + head + expr.len, "\n)\n", 3);
    const md_vial_compile_desc_t cd = {.sys = desc->sys, .run = desc->run, .sets = sets, .num_sets = md_array_size(sets)};
    return md_vial_compile((str_t){src, head + expr.len + 3}, &cd, arena);
}

static bool filter_fail(md_vial_filter_t* out, md_allocator_i* alloc, str_t expr, const char code[6], md_vial_span_t span, uint32_t base, const char* fmt, ...) {
    out->ok = false;
    memcpy(out->error.code, code, 5);
    out->error.severity = code[0] == 'R' ? MD_VIAL_SEVERITY_RUNTIME : MD_VIAL_SEVERITY_ERROR;
    // The span in expr; empty where it is about the statement around it
    if (span.beg >= base && span.end <= base + expr.len && span.end >= span.beg) out->error.span = (md_vial_span_t){span.beg - base, span.end - base};
    char buf[256];
    va_list args;
    va_start(args, fmt);
    vsnprintf(buf, sizeof(buf), fmt, args);
    va_end(args);
    out->error.message = str_copy(str_from_cstr(buf), alloc);
    return false;
}

bool md_vial_filter(md_vial_filter_t* out, str_t expr, const md_vial_filter_desc_t* desc, md_allocator_i* alloc) {
    if (!out) return false;
    memset(out, 0, sizeof(*out));
    const uint32_t base = (uint32_t)(sizeof(FILTER_HEAD) - 1);
    const md_vial_span_t whole = {0, 0};
    if (!desc || !desc->sys || !alloc) return filter_fail(out, alloc ? alloc : md_get_heap_allocator(), expr, "E0101", whole, base, "no system");
    if (str_empty(str_trim(expr))) return filter_fail(out, alloc, expr, "E0101", whole, base, "the filter is empty");
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    bool reads_recipe = false;
    md_vial_program_t* fp = filter_compile(expr, desc, arena, &reads_recipe);
    bool ok = false;
    point_t pt = {0};
    if (!fp) { filter_fail(out, alloc, expr, "E0101", whole, base, "out of memory"); goto done; }
    {
        const md_vial_diagnostic_t* dg;
        const size_t nd = md_vial_diagnostics(fp, &dg);
        for (size_t i = 0; i < nd; ++i) {
            if (dg[i].severity != MD_VIAL_SEVERITY_ERROR) continue;
            filter_fail(out, alloc, expr, dg[i].code, dg[i].span, base, "%.*s", (int)dg[i].message.len, dg[i].message.ptr);
            goto done;
        }
    }
    const vial_ir_t* ir = md_vial__ir(fp);
    const vial_ssa_t v = md_vial__hover_value(fp, 0, NULL);
    if (!ir || !v || v >= md_array_size(ir->values)) { filter_fail(out, alloc, expr, "E0101", whole, base, "the filter has no value"); goto done; }
    const vial_type_t* t = &ir->values[v].type;
    if (t->dtype != VIAL_DT_SET) {
        char ty[128];
        vial_type_print(ty, sizeof(ty), t);
        filter_fail(out, alloc, expr, "E0313", (md_vial_span_t){base, base + (uint32_t)expr.len}, base, "a filter picks particles; this is %s", ty);
        goto done;
    }
    if (t->ndim > 2) { filter_fail(out, alloc, expr, "E0313", whole, base, "a filter gives a set, or an array of sets of one or two dimensions"); goto done; }
    if (t->stage == VIAL_STAGE_TRAJECTORY) { filter_fail(out, alloc, expr, "E0401", whole, base, "a filter picks at one frame; this needs the trajectory"); goto done; }
    out->dynamic = t->stage != VIAL_STAGE_STATIC || reads_recipe;
    if (out->dynamic && !desc->state) { filter_fail(out, alloc, expr, "E0401", whole, base, "it depends on the coordinates, and there are none"); goto done; }
    const vial_val_t* val = point_eval(&pt, fp, desc->sys, desc->run, NULL, 0, 0, desc->state, v, false);
    if (!val || val->dtype != VIAL_DT_SET || !val->sets) {
        if (pt.have_in && pt.in.failed) {
            filter_fail(out, alloc, expr, pt.in.error.code, pt.in.error.span, base, "%s", pt.in.error.message);
        } else {
            filter_fail(out, alloc, expr, "R0103", whole, base, "the filter could not be evaluated");
        }
        goto done;
    }
    {
        // Its sets, cut to the scope
        const uint32_t na = (uint32_t)desc->sys->atom.count;
        uint8_t* in_scope = NULL;
        if (desc->scope) {
            in_scope = md_alloc(arena, na ? na : 1);
            memset(in_scope, 0, na);
            for (size_t i = 0; i < desc->num_scope; ++i) if (desc->scope[i] < na) in_scope[desc->scope[i]] = 1;
        }
        const size_t n = val->ndim ? vial_val_count(val) : 1;
        const size_t inner = val->ndim == 2 ? val->dims[1] : 1;
        uint32_t* offsets = NULL;
        uint32_t* atoms = NULL;
        md_array_push(offsets, 0, alloc);
        size_t rows = 0;
        for (size_t r = 0; r < n / inner; ++r) {
            const size_t mark_a = md_array_size(atoms), mark_o = md_array_size(offsets);
            bool keep = true;
            for (size_t q = 0; q < inner; ++q) {
                const uint32_t si = (uint32_t)(r * inner + q);
                const uint32_t* it = vial_sets_items(val->sets, si);
                const uint32_t sz = vial_sets_size(val->sets, si);
                size_t got = 0;
                for (uint32_t k = 0; k < sz; ++k) {
                    if (in_scope && (it[k] >= na || !in_scope[it[k]])) continue;
                    md_array_push(atoms, it[k], alloc);
                    ++got;
                }
                md_array_push(offsets, (uint32_t)md_array_size(atoms), alloc);
                // An array's set left empty by the scope is dropped, a pair with it
                if (in_scope && val->ndim && !got) keep = false;
            }
            if (!keep) {
                md_array_shrink(atoms, mark_a);
                md_array_shrink(offsets, mark_o);
                continue;
            }
            ++rows;
        }
        out->ndim = val->ndim;
        if (val->ndim >= 1) out->dims[0] = (uint32_t)rows;
        if (val->ndim == 2) out->dims[1] = (uint32_t)inner;
        out->num_sets = md_array_size(offsets) - 1;
        out->offsets = offsets;
        out->atoms = atoms ? atoms : md_alloc(alloc, sizeof(uint32_t));
        out->ok = ok = true;
    }
done:
    point_free(&pt);
    if (fp) md_vial_program_free(fp);
    md_arena_allocator_destroy(arena);
    return ok;
}

bool md_vial_filter_complete(md_vial_completions_t* out, str_t expr, uint32_t cursor, const md_vial_filter_desc_t* desc, md_allocator_i* alloc) {
    if (!out) return false;
    memset(out, 0, sizeof(*out));
    if (!desc || !desc->sys || !alloc || cursor > expr.len) return false;
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(1));
    md_vial_program_t* fp = filter_compile(expr, desc, arena, NULL);
    const uint32_t base = (uint32_t)(sizeof(FILTER_HEAD) - 1);
    const bool ok = fp && md_vial_complete(out, fp, base + cursor, alloc);
    if (ok) {
        out->span.beg = out->span.beg >= base ? out->span.beg - base : 0;
        out->span.end = out->span.end >= base ? MIN(out->span.end - base, (uint32_t)expr.len) : 0;
        // The statement's own name is no name to write
        size_t k = 0;
        for (size_t i = 0; i < out->count; ++i) if (!str_eq(out->items[i].label, STR_LIT("_filter"))) out->items[k++] = out->items[i];
        out->count = k;
    }
    if (fp) md_vial_program_free(fp);
    md_arena_allocator_destroy(arena);
    return ok;
}
