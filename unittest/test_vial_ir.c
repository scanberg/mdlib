#include "utest.h"

#include <core/md_common.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

#include <md_vial.h>
#include "../src/vial/vial_internal.h"
#include "../src/vial/ir.h"
#include "../src/vial/interp.h"
#include <core/md_os.h>
#include <core/md_intrinsics.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "vial_test_systems.h"

// Phase 3: the array IR. Lowering, the passes and the interpreter are tested here; the reference's examples are in
// test_vial_spec.c.

static md_vial_program_t* vial_ir__compile(const char* src, char sys, md_allocator_i* alloc) {
    const md_vial_compile_desc_t desc = {.sys = vial_test_system(sys)};
    return md_vial_compile(str_from_cstr(src), &desc, alloc);
}

// VIAL_IR_SCRIPT=file [VIAL_IR_SYS=T|P]: prints the IR of a script and what was not lowered. A tool, not a test.
UTEST(vial_ir, dump_script) {
    const char* path = getenv("VIAL_IR_SCRIPT");
    if (!path) return;
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const str_t text = load_textfile(str_from_cstr(path), alloc);
    const char* sys = getenv("VIAL_IR_SYS");
    const md_vial_compile_desc_t desc = {.sys = vial_test_system(sys ? sys[0] : 0)};
    md_vial_program_t* prog = md_vial_compile(text, &desc, alloc);
    const vial_ir_t* ir = md_vial__ir(prog);
    const str_t d = vial_ir_dump(ir, text, alloc);
    printf("%.*s", (int)d.len, d.ptr);
    char err[256];
    if (!vial_ir_verify(ir, err, sizeof(err))) printf("does not verify: %s\n", err);
    for (size_t k = 0; k < md_array_size(ir->failures); ++k) {
        const md_vial_span_t sp = ir->failures[k].span;
        printf("not lowered: %s: '%.*s'\n", ir->failures[k].what, (int)(sp.end - sp.beg), text.ptr + sp.beg);
    }
    if (sys && md_vial__run_ir(prog) && vial_test_state(sys[0])) {
        // Evaluated over the structure's own frame, as the fuzz check does
        const vial_ir_t* rir = md_vial__run_ir(prog);
        vial_sys_t vs;
        vial_sys_init(&vs, desc.sys, alloc);
        vial_interp_t in;
        vial_interp_init(&in, rir, &vs, desc.sys, md_get_heap_allocator());
        const md_system_state_t* st = vial_test_state(sys[0]);
        const vial_frame_t fr = {.xyz = (const float*)st->xyz, .num_atoms = st->num_atoms, .cell = st->unitcell};
        vial_interp_frame(&in, &fr);
        for (size_t k = 0; k < md_array_size(rir->outputs); ++k) {
            char tb[128], vb[200] = "(not computed)";
            vial_type_print(tb, sizeof(tb), &rir->values[rir->outputs[k].value].type);
            const vial_val_t* v = vial_interp_value(&in, rir->outputs[k].value);
            if (v) vial_val_print(vb, sizeof(vb), v, 4);
            printf("= %.*s : %s = %s (dtype %d, rank %u)\n", (int)rir->outputs[k].name.len, rir->outputs[k].name.ptr, tb, vb, v ? v->dtype : -1, v ? v->ndim : 0);
        }
        for (size_t k = 0; k < md_array_size(in.errors); ++k) printf("%s %s\n", in.errors[k].code, in.errors[k].message);
        vial_interp_free(&in);
    }
    const md_vial_diagnostic_t* diags;
    const size_t nd = md_vial_diagnostics(prog, &diags);
    for (size_t k = 0; k < nd; ++k) printf("%.5s %.*s\n", diags[k].code, (int)diags[k].message.len, diags[k].message.ptr);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

UTEST(vial_ir, lowers_and_verifies) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    md_vial_program_t* prog = vial_ir__compile("a = 1 + 2\nb = residues(chain='A')\nc = distance(b.com, b[0].com)\n", 'T', alloc);
    const vial_ir_t* ir = md_vial__ir(prog);
    char err[256];
    EXPECT_TRUE(vial_ir_verify(ir, err, sizeof(err)));
    EXPECT_EQ(0, (int)md_array_size(ir->failures));
    EXPECT_EQ(3, (int)md_array_size(ir->outputs));
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

// ### EVALUATION ###

#include "../src/vial/interp.h"
#include <md_pdb.h>
#include <md_gro.h>
#include <md_util.h>
#include <core/md_vec_math.h>

typedef struct vial_ir__run_t {
    md_allocator_i*     alloc;
    md_vial_program_t*  prog;
    vial_sys_t          vs;
    vial_interp_t       in;
    bool                ok;
} vial_ir__run_t;

static bool vial_ir__start(vial_ir__run_t* r, const char* src, const md_system_t* sys) {
    memset(r, 0, sizeof(*r));
    r->alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const md_vial_compile_desc_t desc = {.sys = sys};
    r->prog = md_vial_compile(str_from_cstr(src), &desc, r->alloc);
    const vial_ir_t* ir = md_vial__run_ir(r->prog);
    if (!ir) {
        const vial_ir_t* lowered = md_vial__ir(r->prog);
        printf("not lowered: %s\n", md_array_size(lowered->failures) ? lowered->failures[0].what : "does not inline");
        return false;
    }
    if (sys) vial_sys_init(&r->vs, sys, r->alloc);
    if (!vial_interp_init(&r->in, ir, sys ? &r->vs : NULL, sys, r->alloc)) return false;
    return true;
}

static void vial_ir__end(vial_ir__run_t* r) {
    vial_interp_free(&r->in);
    md_vial_program_free(r->prog);
    md_arena_allocator_destroy(r->alloc);
}

// The printed value of a published name, or what failed
static const char* vial_ir__out(vial_ir__run_t* r, const char* name, char* buf, size_t cap) {
    const vial_ir_t* ir = md_vial__run_ir(r->prog);
    for (size_t i = 0; i < md_array_size(ir->outputs); ++i) {
        if (!str_eq(ir->outputs[i].name, str_from_cstr(name))) continue;
        const vial_val_t* v = vial_interp_value(&r->in, ir->outputs[i].value);
        if (!v) { snprintf(buf, cap, "(not computed: %s)", r->in.failed ? r->in.error.message : "-"); return buf; }
        vial_val_print(buf, cap, v, 6);
        return buf;
    }
    snprintf(buf, cap, "(no output %s)", name);
    return buf;
}

UTEST(vial_ir, evaluates_static) {
    vial_ir__run_t r;
    const char* src =
        "a = 1 + 2 * 3\n"
        "b = 7 // 2\n"
        "c = -7 % 3\n"
        "d = [1.5, 2.5, 3.5]\n"
        "e = sum(d)\n"
        "f = d[1:]\n"
        "g = mean(d)\n"
        "h = max([3, 9, 2])\n"
        "i = argmin([3, 9, 2])\n"
        "j = [x * 2 for x in range(4)]\n"
        "k = [x for x in range(6) if x % 2 == 0]\n"
        "l = zeros(3)\n"
        "m = sqrt(16.0)\n"
        "n = len(residues(chain='A'))\n"
        "o = count(residues(chain='A'))[0:3]\n"
        "p = where([True, False], [1, 2], [3, 4])\n"
        "q = concatenate(([1, 2], [3]))\n";
    ASSERT_TRUE(vial_ir__start(&r, src, vial_test_system('T')));
    const bool ok = vial_interp_static(&r.in);
    if (!ok) printf("failed: %s %s\n", r.in.error.code, r.in.error.message);
    EXPECT_TRUE(ok);
    char buf[256];
    EXPECT_STREQ("7", vial_ir__out(&r, "a", buf, sizeof(buf)));
    EXPECT_STREQ("3", vial_ir__out(&r, "b", buf, sizeof(buf)));
    EXPECT_STREQ("2", vial_ir__out(&r, "c", buf, sizeof(buf)));
    EXPECT_STREQ("[1.5, 2.5, 3.5]", vial_ir__out(&r, "d", buf, sizeof(buf)));
    EXPECT_STREQ("7.5", vial_ir__out(&r, "e", buf, sizeof(buf)));
    EXPECT_STREQ("[2.5, 3.5]", vial_ir__out(&r, "f", buf, sizeof(buf)));
    EXPECT_STREQ("2.5", vial_ir__out(&r, "g", buf, sizeof(buf)));
    EXPECT_STREQ("9", vial_ir__out(&r, "h", buf, sizeof(buf)));
    EXPECT_STREQ("2", vial_ir__out(&r, "i", buf, sizeof(buf)));
    EXPECT_STREQ("[0, 2, 4, 6]", vial_ir__out(&r, "j", buf, sizeof(buf)));
    EXPECT_STREQ("[0, 2, 4]", vial_ir__out(&r, "k", buf, sizeof(buf)));
    EXPECT_STREQ("[0.0, 0.0, 0.0]", vial_ir__out(&r, "l", buf, sizeof(buf)));
    EXPECT_STREQ("4.0", vial_ir__out(&r, "m", buf, sizeof(buf)));
    EXPECT_STREQ("442", vial_ir__out(&r, "n", buf, sizeof(buf)));
    EXPECT_STREQ("[8, 11, 9]", vial_ir__out(&r, "o", buf, sizeof(buf)));
    EXPECT_STREQ("[1, 4]", vial_ir__out(&r, "p", buf, sizeof(buf)));
    EXPECT_STREQ("[1, 2, 3]", vial_ir__out(&r, "q", buf, sizeof(buf)));
    vial_ir__end(&r);
}

// The run of P: 1ALA-560ns.pdb's models as frames
static const md_system_t* vial_ir__run_p(int64_t* num_frames) {
    static md_system_t sys;
    static bool loaded;
    static int64_t frames;
    if (!loaded) {
        loaded = true;
        md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(16));
        sys = (md_system_t){.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        const str_t path = STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");
        if (md_pdb_system_init_from_file(&sys, &state, path, MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE)) {
            md_util_system_infer(&sys, &state, MD_UTIL_INFER_ALL);
            if (md_pdb_system_publish_run(&sys, path, STR_LIT("run"), MD_RUN_FLAG_DISABLE_CACHE_WRITE)) {
                const md_attribute_t* t = md_attributes_find(&sys.attributes, STR_LIT("run/time"));
                frames = t ? (int64_t)t->format.shape[0] : 0;
            }
        }
    }
    *num_frames = frames;
    return frames ? &sys : NULL;
}

// Every frame of the run through the interpreter; per frame, out[f] gets the value of name
static bool vial_ir__frames(vial_ir__run_t* r, const md_system_t* sys, int64_t num_frames, const char* name, double* out) {
    md_system_state_t state = {.alloc = r->alloc};
    if (!md_system_state_init(&state, sys->atom.count)) return false;
    const str_t paths[] = {STR_LIT("atom/position"), STR_LIT("unitcell")};
    md_system_extract_t* ex = md_system_extract_begin(sys, STR_LIT("run"), paths, 2, r->alloc);
    if (!ex) return false;
    const vial_ir_t* ir = md_vial__run_ir(r->prog);
    vial_ssa_t v = 0;
    for (size_t i = 0; i < md_array_size(ir->outputs); ++i) if (str_eq(ir->outputs[i].name, str_from_cstr(name))) v = ir->outputs[i].value;
    bool ok = v != 0;
    for (int64_t f = 0; ok && f < num_frames; ++f) {
        ok = md_system_extract_frame(ex, f, &state);
        const vial_frame_t fr = {.xyz = (const float*)state.xyz, .num_atoms = state.num_atoms, .cell = state.unitcell, .index = f, .time = (double)f};
        if (ok) ok = vial_interp_frame(&r->in, &fr);
        if (!ok && r->in.failed) printf("frame %lld: %s %s\n", (long long)f, r->in.error.code, r->in.error.message);
        const vial_val_t* val = ok ? vial_interp_value(&r->in, v) : NULL;
        if (!val || !val->data) { ok = false; break; }
        out[f] = val->dtype == VIAL_DT_FLOAT ? *(const float*)val->data : val->dtype == VIAL_DT_INT ? *(const int32_t*)val->data : 0.0;
    }
    md_system_extract_end(ex);
    return ok;
}

UTEST(vial_ir, evaluates_frames) {
    int64_t nf = 0;
    const md_system_t* sys = vial_ir__run_p(&nf);
    ASSERT_TRUE(sys != NULL);
    ASSERT_GT(nf, 1);
    vial_ir__run_t r;
    const char* src =
        "ca = particles(name='CA')\n"
        "d  = distance(ca[0], ca[3])\n"
        "a  = angle(ca[0], ca[1], ca[2])\n"
        "t  = dihedral(ca[0], ca[1], ca[2], ca[3])\n"
        "c  = distance(ca.com, ca[0])\n"
        "dd = distance(ca.xyz, ca[0])\n"
        "n5 = sum(dd > 5.0)\n"
        "far = argmax(dd)\n";
    ASSERT_TRUE(vial_ir__start(&r, src, sys));
    double* d = malloc(nf * sizeof(double));
    double* n5 = malloc(nf * sizeof(double));
    ASSERT_TRUE(vial_ir__frames(&r, sys, nf, "d", d));
    // By hand, from the frames' positions
    md_system_state_t state = {.alloc = r.alloc};
    md_system_state_init(&state, sys->atom.count);
    const str_t paths[] = {STR_LIT("atom/position"), STR_LIT("unitcell")};
    md_system_extract_t* ex = md_system_extract_begin(sys, STR_LIT("run"), paths, 2, r.alloc);
    uint32_t ca[4], k = 0;
    for (uint32_t i = 0; i < sys->atom.count && k < 4; ++i) if (str_eq(md_atom_name(&sys->atom, i), STR_LIT("CA"))) ca[k++] = i;
    ASSERT_EQ(4u, k);
    for (int64_t f = 0; f < nf; ++f) {
        md_system_extract_frame(ex, f, &state);
        vec3_t dx = vec3_sub(state.xyz[ca[3]], state.xyz[ca[0]]);
        md_util_min_image_vec3(&dx, 1, &state.unitcell);
        const double want = sqrtf(dx.x * dx.x + dx.y * dx.y + dx.z * dx.z);
        EXPECT_NEAR(want, d[f], 1e-5);
    }
    md_system_extract_end(ex);
    ASSERT_TRUE(vial_ir__frames(&r, sys, nf, "n5", n5));
    char buf[256];
    printf("frame %lld: d = %s, a = %s, t = %s, c = %s, n5 = %s, far = %s\n", (long long)nf - 1, vial_ir__out(&r, "d", buf, sizeof(buf)),
        vial_ir__out(&r, "a", buf + 40, sizeof(buf) - 40), vial_ir__out(&r, "t", buf + 80, sizeof(buf) - 80),
        vial_ir__out(&r, "c", buf + 120, sizeof(buf) - 120), vial_ir__out(&r, "n5", buf + 160, sizeof(buf) - 160),
        vial_ir__out(&r, "far", buf + 200, sizeof(buf) - 200));
    free(d);
    free(n5);
    vial_ir__end(&r);
}

UTEST(vial_ir, test_runs_load) {
    EXPECT_GT(vial_test_run('P'), 1);
    EXPECT_EQ(0, vial_test_run('T'));
    ASSERT_TRUE(vial_test_state('T') != NULL);
    EXPECT_EQ(vial_test_system('T')->atom.count, vial_test_state('T')->num_atoms);
}

typedef struct vial_ir__loader_t {
    const md_system_t*  sys;
    md_system_state_t   state;
} vial_ir__loader_t;

static bool vial_ir__load(void* user, int64_t index, vial_frame_t* out) {
    vial_ir__loader_t* l = user;
    const str_t paths[] = {STR_LIT("atom/position"), STR_LIT("unitcell")};
    md_system_extract_t* ex = md_system_extract_begin(l->sys, STR_LIT("run"), paths, 2, md_get_heap_allocator());
    const bool ok = ex && md_system_extract_frame(ex, index, &l->state);
    if (ex) md_system_extract_end(ex);
    *out = (vial_frame_t){.xyz = (const float*)l->state.xyz, .num_atoms = l->state.num_atoms, .cell = l->state.unitcell, .index = index, .time = (double)index};
    return ok;
}

static double vial_ir__num(const vial_val_t* v) {
    if (!v || !v->data || vial_val_is_missing(v, 0)) return NAN;
    return v->dtype == VIAL_DT_FLOAT ? *(const float*)v->data : v->dtype == VIAL_DT_INT ? *(const int32_t*)v->data : *(const double*)v->data;
}

static vial_ssa_t vial_ir__output(const vial_ir_t* ir, const char* name) {
    for (size_t i = 0; i < md_array_size(ir->outputs); ++i) if (str_eq(ir->outputs[i].name, str_from_cstr(name))) return ir->outputs[i].value;
    return 0;
}

UTEST(vial_ir, evaluates_trajectory) {
    const md_system_t* sys = vial_test_system('P');
    const int64_t nf = vial_test_run('P');
    ASSERT_TRUE(sys && nf > 2);
    vial_ir__run_t r;
    const char* src =
        "ca   = particles(name='CA')\n"
        "ends = distance(ca[0], ca[14])\n"
        "m    = mean(ends, axis='frame')\n"
        "mx   = max(ends, axis='frame')\n"
        "sm   = smooth(ends, 3)\n"
        "prev = previous(ends)\n"
        "def step(s, x):\n"
        "    return s + x\n"
        "cum  = fold(step, ends, init=zeros_like(ends))\n"
        "dev  = ends - m\n";
    ASSERT_TRUE(vial_ir__start(&r, src, sys));
    vial_ir__loader_t l = {.sys = sys, .state = {.alloc = r.alloc}};
    md_system_state_init(&l.state, sys->atom.count);
    r.in.load_frame = vial_ir__load;
    r.in.load_user = &l;
    const bool ok = vial_interp_run(&r.in, nf);
    for (size_t i = 0; i < md_array_size(r.in.errors); ++i) printf("failed: %s %s\n", r.in.errors[i].code, r.in.errors[i].message);
    ASSERT_TRUE(ok);
    const vial_ir_t* ir = md_vial__run_ir(r.prog);
    if (getenv("VIAL_IR_DUMP")) {
        const str_t d = vial_ir_dump(ir, str_from_cstr(src), r.alloc);
        printf("%.*s", (int)d.len, d.ptr);
        for (size_t v = 1; v < md_array_size(ir->values); ++v) {
            const uint32_t op = ir->values[v].op;
            if (op == VIAL_IR_NO_OP || ir->ops[op].parent != ir->funcs[0].body) continue;
            printf("%%%zu: class %u sweep %u live %u state %u series %p\n", v, r.in.op_class[op], r.in.op_sweep[op], r.in.live[op], r.in.state[v], (void*)r.in.series[v]);
        }
    }
    const vial_ssa_t ends = vial_ir__output(ir, "ends"), m = vial_ir__output(ir, "m"), mx = vial_ir__output(ir, "mx"), sm = vial_ir__output(ir, "sm");
    const vial_ssa_t prev = vial_ir__output(ir, "prev"), cum = vial_ir__output(ir, "cum"), dev = vial_ir__output(ir, "dev");
    double sum = 0.0, best = -1e30;
    float run = 0.0f;
    for (int64_t f = 0; f < nf; ++f) {
        const double e = vial_ir__num(vial_interp_at(&r.in, ends, f));
        sum += e;
        if (e > best) best = e;
        run += (float)e;
        EXPECT_NEAR(run, vial_ir__num(vial_interp_at(&r.in, cum, f)), 1e-3);
        if (f == 0) EXPECT_TRUE(isnan(vial_ir__num(vial_interp_at(&r.in, prev, f))));
        else EXPECT_EQ((float)vial_ir__num(vial_interp_at(&r.in, ends, f - 1)), (float)vial_ir__num(vial_interp_at(&r.in, prev, f)));
        if (f > 0 && f + 1 < nf) {
            const double want = (vial_ir__num(vial_interp_at(&r.in, ends, f - 1)) + e + vial_ir__num(vial_interp_at(&r.in, ends, f + 1))) / 3.0;
            const vial_val_t* sv = vial_interp_at(&r.in, sm, f);
            ASSERT_TRUE(sv != NULL);
            EXPECT_NEAR(want, ((const float*)sv->data)[f], 1e-4);     // smooth keeps the whole series
        }
    }
    EXPECT_NEAR(sum / nf, vial_ir__num(vial_interp_at(&r.in, m, 0)), 1e-4);
    EXPECT_EQ((float)best, (float)vial_ir__num(vial_interp_at(&r.in, mx, 0)));
    const double mean = vial_ir__num(vial_interp_at(&r.in, m, 0));
    EXPECT_NEAR(vial_ir__num(vial_interp_at(&r.in, ends, 5)) - mean, vial_ir__num(vial_interp_at(&r.in, dev, 5)), 1e-4);
    md_system_state_free(&l.state);
    vial_ir__end(&r);
}

UTEST(vial_ir, eval_api) {
    const md_system_t* sys = vial_test_system('P');
    const int64_t nf = vial_test_run('P');
    ASSERT_TRUE(sys && nf > 2);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const char* src =
        "ca   = particles(name='CA')\n"
        "ends = distance(ca[0], ca[14])\n"
        "m    = mean(ends, axis='frame')\n"
        "dd   = distance(ca.xyz, ca[0])\n"
        "_hid = ends * 2.0\n"
        "big  = zeros(100) + ends\n";
    const md_vial_compile_desc_t cdesc = {.sys = sys, .run = STR_LIT("run")};
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &cdesc, alloc);
    const md_vial_eval_desc_t edesc = {.sys = sys, .run = STR_LIT("run"), .max_elements = 50};
    md_vial_eval_t* eval = md_vial_eval_create(prog, &edesc, alloc);
    ASSERT_TRUE(eval != NULL);
    {
        // Before the evaluation has run: no value over the trajectory yet
        char hv[64];
        EXPECT_FALSE(md_vial_eval_hover(eval, (uint32_t)(strstr(src, "mean(") - src) + 4, 0, hv, sizeof(hv)));
    }
    // A value of 100 elements is over this evaluation's limit: R0105, and the rest is computed
    EXPECT_FALSE(md_vial_eval_run(eval));
    EXPECT_EQ(1.0, md_vial_eval_progress(eval));
    const md_vial_diagnostic_t* d;
    ASSERT_EQ(1u, (unsigned)md_vial_eval_diagnostics(eval, &d));
    EXPECT_EQ(0, strncmp(d[0].code, "R0105", 5));
    EXPECT_EQ(MD_VIAL_SEVERITY_RUNTIME, d[0].severity);
    const md_attributes_t* attrs = md_vial_eval_attributes(eval);
    const md_attribute_t* ends = md_attributes_find(attrs, STR_LIT("script/ends"));
    const md_attribute_t* m = md_attributes_find(attrs, STR_LIT("script/m"));
    const md_attribute_t* dd = md_attributes_find(attrs, STR_LIT("script/dd"));
    ASSERT_TRUE(ends && m && dd);
    EXPECT_TRUE(ends->flags & MD_ATTRIBUTE_FLAG_TEMPORAL);
    EXPECT_FALSE(m->flags & MD_ATTRIBUTE_FLAG_TEMPORAL);
    EXPECT_EQ(2u, dd->format.rank);
    EXPECT_EQ((uint32_t)nf, dd->format.shape[0]);
    EXPECT_EQ(15u, dd->format.shape[1]);
    EXPECT_TRUE(md_attributes_find(attrs, STR_LIT("script/_hid")) == NULL);
    EXPECT_TRUE(md_attributes_find(attrs, STR_LIT("script/ca")) == NULL);      // a selection is a recipe
    double sum = 0.0, mean = 0.0;
    for (int64_t f = 0; f < nf; ++f) {
        double v;
        md_attribute_extract_f64(&v, 1, ends, md_attribute_slice_1((uint32_t)f), ends->unit);
        sum += v;
    }
    md_attribute_extract_f64(&mean, 1, m, md_attribute_slice_all(), m->unit);
    EXPECT_NEAR(sum / nf, mean, 1e-4);

    // Hover: an expression's value at a frame, the bits the evaluation published
    {
        char hv[160], want[64];
        const uint32_t call = (uint32_t)(strstr(src, "distance(ca[0]") - src) + 8;     // its '('
        for (int64_t f = 0; f < nf; f += nf / 3) {
            float v;
            md_attribute_extract_f32(&v, 1, ends, md_attribute_slice_1((uint32_t)f), ends->unit);
            vial_print_float(want, sizeof(want), v, true);
            EXPECT_TRUE(md_vial_eval_hover(eval, call, f, hv, sizeof(hv)));
            EXPECT_STREQ(want, hv);
        }
        // A static value, and one past the evaluation's budget but within hover's
        EXPECT_TRUE(md_vial_eval_hover(eval, (uint32_t)(strstr(src, "ca[0]") - src) + 2, 0, hv, sizeof(hv)));
        EXPECT_TRUE(hv[0] != 0);
        EXPECT_TRUE(md_vial_eval_hover(eval, (uint32_t)(strstr(src, "zeros(100)") - src) + 5, 1, hv, sizeof(hv)));
        EXPECT_EQ('[', hv[0]);
        // Over the trajectory: the evaluation's value, once it has run every frame. One value failed (R0105); the
        // rest, the mean included, was computed and published
        {
            double mv;
            md_attribute_extract_f64(&mv, 1, m, md_attribute_slice_all(), m->unit);
            char wm[64];
            vial_print_float(wm, sizeof(wm), (float)mv, true);
            EXPECT_TRUE(md_vial_eval_hover(eval, (uint32_t)(strstr(src, "mean(") - src) + 4, 0, hv, sizeof(hv)));
            EXPECT_STREQ(wm, hv);
        }
        EXPECT_FALSE(md_vial_eval_hover(eval, 0xFFFFFFu, 0, hv, sizeof(hv)));
    }
    md_vial_eval_free(eval);

    // Interrupted before it starts: a value per frame is there, as written so far, which is nothing; the rest is not
    eval = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    md_vial_eval_interrupt(eval);
    EXPECT_FALSE(md_vial_eval_run(eval));
    {
        const md_attribute_t* e2 = md_attributes_find(md_vial_eval_attributes(eval), STR_LIT("script/ends"));
        ASSERT_TRUE(e2 != NULL);
        for (int64_t f = 0; f < nf; ++f) EXPECT_TRUE(isnan(((const float*)e2->data)[f]));
        EXPECT_TRUE(md_attributes_find(md_vial_eval_attributes(eval), STR_LIT("script/m")) == NULL);
        uint8_t held[256];
        ASSERT_EQ((size_t)nf, md_vial_eval_output_frames(eval, STR_LIT("ends"), held, sizeof(held)));
        for (int64_t f = 0; f < nf; ++f) EXPECT_EQ(0, held[f]);
    }
    md_vial_eval_free(eval);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

static uint64_t vial_ir__fp(const md_vial_output_t* o, size_t n, const char* name) {
    for (size_t i = 0; i < n; ++i) if (str_eq(o[i].name, str_from_cstr(name))) return o[i].fingerprint;
    return 0;
}

// A caller's pool for md_vial_eval: threads of its own per call, taking indices in turn
typedef struct vial_ir__pool_t {
    size_t              count;
    void              (*fn)(void* arg, size_t i);
    void*               arg;
    volatile uint64_t   next;
} vial_ir__pool_t;

static void vial_ir__pool_work(void* data) {
    vial_ir__pool_t* p = data;
    for (uint64_t i; (i = md_atomic_fetch_add64(&p->next, 1)) < p->count;) p->fn(p->arg, (size_t)i);
}

static void vial_ir__parallel_for(void* user, size_t count, void (*fn)(void* arg, size_t i), void* arg) {
    (void)user;
    vial_ir__pool_t pool = {.count = count, .fn = fn, .arg = arg};
    md_thread_t* th[2];
    for (int t = 0; t < 2; ++t) th[t] = md_thread_create(vial_ir__pool_work, &pool);
    vial_ir__pool_work(&pool);
    for (int t = 0; t < 2; ++t) if (th[t]) md_thread_join(th[t]);
}

// The evaluation API on several threads, its own and a caller's: the same attributes, byte for byte, and the same
// diagnostics in the same order as on one
UTEST(vial_ir, eval_threads) {
    const md_system_t* sys = vial_test_system('P');
    const int64_t nf = vial_test_run('P');
    ASSERT_TRUE(sys && nf > 2);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const char* src =
        "ca   = particles(name='CA')\n"
        "ends = distance(ca[0], ca[14])\n"
        "m    = mean(ends, axis='frame')\n"
        "dev  = ends - m\n"
        "near = len(within(6.0, of=ca[0:3]))\n"
        "cp   = distance(ca[0:5], ca[8:15], between='closest')\n"
        "big  = zeros(100) + ends\n"
        "acc  = cumulative(ends, op='sum')\n";
    const md_vial_compile_desc_t cdesc = {.sys = sys, .run = STR_LIT("run")};
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &cdesc, alloc);
    ASSERT_TRUE(prog != NULL);
    md_vial_eval_t* ev[3];
    ev[0] = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run"), .max_elements = 50}, alloc);
    ev[1] = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run"), .max_elements = 50, .num_threads = 4}, alloc);
    ev[2] = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run"), .max_elements = 50, .num_threads = 3,
                                                             .parallel_for = vial_ir__parallel_for}, alloc);
    for (int i = 0; i < 3; ++i) {
        ASSERT_TRUE(ev[i] != NULL);
        EXPECT_FALSE(md_vial_eval_run(ev[i]));     // big is past the limit: R0105 at every frame, once
        EXPECT_EQ(1.0, md_vial_eval_progress(ev[i]));
    }
    const md_vial_diagnostic_t* d0;
    const size_t n0 = md_vial_eval_diagnostics(ev[0], &d0);
    EXPECT_GE(n0, (size_t)1);
    const md_attributes_t* a0 = md_vial_eval_attributes(ev[0]);
    for (int i = 1; i < 3; ++i) {
        const md_vial_diagnostic_t* di;
        ASSERT_EQ(n0, md_vial_eval_diagnostics(ev[i], &di));
        for (size_t k = 0; k < n0; ++k) {
            EXPECT_EQ(0, strncmp(d0[k].code, di[k].code, 5));
            EXPECT_EQ(d0[k].span.beg, di[k].span.beg);
        }
        const md_attributes_t* ai = md_vial_eval_attributes(ev[i]);
        EXPECT_EQ(md_attributes_count(a0), md_attributes_count(ai));
        int compared = 0;
        for (md_attribute_iter_t it = md_attributes_iter(a0, STR_LIT("script")); md_attributes_next(&it);) {
            const md_attribute_t* x = it.attr;
            const md_attribute_t* y = md_attributes_find(ai, x->path);
            ASSERT_TRUE(y != NULL);
            const size_t bytes = md_attribute_byte_size(&x->format);
            ASSERT_EQ(bytes, md_attribute_byte_size(&y->format));
            EXPECT_TRUE(x->data && y->data && memcmp(x->data, y->data, bytes) == 0);
            compared += 1;
        }
        EXPECT_GE(compared, 6);
    }
    for (int i = 0; i < 3; ++i) md_vial_eval_free(ev[i]);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

// Fingerprints: what a value is computed from, not what it is called or where it is written
UTEST(vial_ir, output_fingerprints) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    const int64_t nf = vial_test_run('P');     // publishes the run, before anything is compiled against it
    ASSERT_GT(nf, 2);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const char* a =
        "ca = particles(name='CA')\n"
        "d  = distance(ca[0], ca[1])\n"
        "e  = d * 2.0\n"
        "f  = distance(ca[0], ca[2])\n"
        "h  = histogram(d, axis='frame', bins=10, range=(0.0, 20.0))\n"
        "m  = mean(d, axis='frame')\n";
    const char* b =
        "# the same, renamed and moved\n"
        "ca = particles(name='CA')\n\n"
        "dist = distance(ca[0],   ca[1])\n"
        "m  = mean(dist, axis='frame')\n"
        "e  = dist * 3.0\n"
        "f  = distance(ca[0], ca[3])\n"
        "h  = histogram(dist, axis='frame', bins=10, range=(0.0, 20.0))\n";
    const md_vial_compile_desc_t desc = {.sys = sys, .run = STR_LIT("run")};
    md_vial_program_t* pa = md_vial_compile(str_from_cstr(a), &desc, alloc);
    md_vial_program_t* pb = md_vial_compile(str_from_cstr(b), &desc, alloc);
    const md_vial_output_t* oa;
    const md_vial_output_t* ob;
    const size_t na = md_vial_outputs(pa, &oa), nb = md_vial_outputs(pb, &ob);
    #define FP(o, n, nm) vial_ir__fp(o, n, nm)
    EXPECT_NE(0u, (unsigned)(FP(oa, na, "d") != 0));
    EXPECT_TRUE(FP(oa, na, "d") == FP(ob, nb, "dist"));     // renamed, spaced: the same
    EXPECT_TRUE(FP(oa, na, "m") == FP(ob, nb, "m"));        // reads the same, bound elsewhere
    EXPECT_TRUE(FP(oa, na, "h") == FP(ob, nb, "h"));
    EXPECT_TRUE(FP(oa, na, "e") != FP(ob, nb, "e"));        // another constant
    EXPECT_TRUE(FP(oa, na, "f") != FP(ob, nb, "f"));        // another particle
    EXPECT_TRUE(FP(oa, na, "d") != FP(oa, na, "f"));
    for (size_t i = 0; i < na; ++i) {
        if (str_eq(oa[i].name, STR_LIT("d"))) EXPECT_TRUE(oa[i].flags & MD_VIAL_OUTPUT_TEMPORAL);
        if (str_eq(oa[i].name, STR_LIT("h"))) { EXPECT_TRUE(oa[i].flags & MD_VIAL_OUTPUT_DISTRIBUTION); EXPECT_FALSE(oa[i].flags & MD_VIAL_OUTPUT_TEMPORAL); }
        if (str_eq(oa[i].name, STR_LIT("m"))) EXPECT_EQ(0u, oa[i].flags);
    }
    #undef FP

    // Taken over: an evaluation of b made from one of a publishes what a computed alike as a published it, and
    // computes the rest
    {
        md_vial_eval_t* ea = md_vial_eval_create(pa, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
        ASSERT_TRUE(ea != NULL);
        md_vial_plan_t plan[16];
        // Before a has run every frame, nothing is held
        ASSERT_TRUE(md_vial_eval_plan(plan, 16, ea, pb, (size_t)nf));
        for (size_t i = 0; i < nb; ++i) EXPECT_NE((int)MD_VIAL_PLAN_HELD, (int)plan[i].status);
        EXPECT_TRUE(md_vial_eval_run(ea));
        ASSERT_TRUE(md_vial_eval_plan(plan, 16, ea, pb, (size_t)nf));
        for (size_t i = 0; i < nb; ++i) {
            const char* want = str_eq(ob[i].name, STR_LIT("e")) || str_eq(ob[i].name, STR_LIT("f")) ? "changed" : ob[i].selection ? "new" : "held";
            const char* got = plan[i].status == MD_VIAL_PLAN_HELD ? "held" : plan[i].status == MD_VIAL_PLAN_CHANGED ? "changed" : "new";
            if (str_eq(ob[i].name, STR_LIT("ca"))) continue;   // a selection: never taken over, evaluated as it is
            EXPECT_STREQ(want, got);
            if (str_eq(ob[i].name, STR_LIT("dist"))) EXPECT_TRUE(str_eq(plan[i].from, STR_LIT("d")));
        }
        md_vial_eval_t* eb = md_vial_eval_create_from(ea, pb, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
        ASSERT_TRUE(eb != NULL);
        // Taken over: published from when it is made, before any run, with the histogram's range
        EXPECT_TRUE(md_attributes_find(md_vial_eval_attributes(eb), STR_LIT("script/dist")) != NULL);
        EXPECT_TRUE(md_attributes_find(md_vial_eval_attributes(eb), STR_LIT("script/h/range")) != NULL);
        EXPECT_TRUE(md_vial_eval_run(eb));
        const md_attributes_t* ta = md_vial_eval_attributes(ea);
        const md_attributes_t* tb = md_vial_eval_attributes(eb);
        const char* pairs[][2] = {{"script/d", "script/dist"}, {"script/m", "script/m"}, {"script/h", "script/h"}};
        for (int k = 0; k < 3; ++k) {
            const md_attribute_t* x = md_attributes_find(ta, str_from_cstr(pairs[k][0]));
            const md_attribute_t* y = md_attributes_find(tb, str_from_cstr(pairs[k][1]));
            ASSERT_TRUE(x && y);
            const size_t bytes = md_attribute_byte_size(&x->format);
            ASSERT_EQ(bytes, md_attribute_byte_size(&y->format));
            EXPECT_EQ(0, memcmp(x->data, y->data, bytes));
            EXPECT_EQ(x->flags, y->flags);
        }
        // Changed: computed again, e = dist * 3 where a had d * 2
        const md_attribute_t* ae = md_attributes_find(ta, STR_LIT("script/e"));
        const md_attribute_t* be = md_attributes_find(tb, STR_LIT("script/e"));
        ASSERT_TRUE(ae && be);
        EXPECT_NEAR(((const float*)ae->data)[3] * 1.5f, ((const float*)be->data)[3], 1e-4f);
        md_vial_eval_free(eb);
        md_vial_eval_free(ea);
    }
    md_vial_program_free(pa);
    md_vial_program_free(pb);
    md_arena_allocator_destroy(alloc);
}

// A plot period: outputs over the frames computed again over chosen frames from the kept series
UTEST(vial_ir, eval_aggregate) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    const int64_t nf = vial_test_run('P');
    ASSERT_GT(nf, 3);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const char* src =
        "ca = particles(name='CA')\n"
        "d  = distance(ca[0], ca[1])\n"
        "h  = histogram(d, axis='frame', bins=10, range=(0.0, 20.0))\n"
        "m  = mean(d, axis='frame')\n"
        "x  = max(d * 2.0, axis='frame')\n";
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    ASSERT_TRUE(prog != NULL);
    for (int threads = 1; threads <= 3; threads += 2) {
        md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run"), .num_threads = threads}, alloc);
        ASSERT_TRUE(e != NULL);
        uint8_t* all = md_alloc(alloc, (size_t)nf);
        uint8_t* odd = md_alloc(alloc, (size_t)nf);
        for (int64_t f = 0; f < nf; ++f) { all[f] = 1; odd[f] = (uint8_t)(f & 1); }
        float h[10];
        EXPECT_EQ(0u, md_vial_eval_aggregate(e, STR_LIT("h"), all, h, 10));   // not run yet
        // Stopped before it starts, it stays stopped until reset
        md_vial_eval_interrupt(e);
        EXPECT_FALSE(md_vial_eval_run(e));
        EXPECT_FALSE(md_vial_eval_run(e));
        md_vial_eval_reset_interrupt(e);
        EXPECT_TRUE(md_vial_eval_run(e));
        const md_attributes_t* attr = md_vial_eval_attributes(e);
        const md_attribute_t* ad = md_attributes_find(attr, STR_LIT("script/d"));
        const md_attribute_t* ah = md_attributes_find(attr, STR_LIT("script/h"));
        const md_attribute_t* am = md_attributes_find(attr, STR_LIT("script/m"));
        ASSERT_TRUE(ad && ah && am);
        const float* d = (const float*)ad->data;
        // Beside the histogram, its range: the bin edges, in the unit of what it counts
        const md_attribute_t* ar = md_attributes_find(attr, STR_LIT("script/h/range"));
        ASSERT_TRUE(ar != NULL);
        double rr[2];
        ASSERT_EQ(2u, md_attribute_extract_f64(rr, 2, ar, md_attribute_slice_all(), ar->unit));
        EXPECT_EQ(0.0, rr[0]);
        EXPECT_EQ(20.0, rr[1]);
        EXPECT_TRUE(md_unit_base_equal(ar->unit, ad->unit));

        // Over every frame: what was published
        ASSERT_EQ(10u, md_vial_eval_aggregate(e, STR_LIT("h"), all, h, 10));
        for (int i = 0; i < 10; ++i) EXPECT_EQ(((const float*)ah->data)[i], h[i]);
        float m = 0;
        ASSERT_EQ(1u, md_vial_eval_aggregate(e, STR_LIT("m"), all, &m, 1));
        EXPECT_EQ(((const float*)am->data)[0], m);

        // Over the odd frames: counted by hand from the per-frame distances
        float want[10] = {0};
        double sum = 0, mx = -INFINITY;
        int n = 0;
        for (int64_t f = 1; f < nf; f += 2) {
            const int bin = (int)(d[f] / 2.0f);
            if (bin >= 0 && bin < 10) want[bin] += 1.0f;
            sum += d[f]; mx = MAX(mx, 2.0 * d[f]); ++n;
        }
        ASSERT_EQ(10u, md_vial_eval_aggregate(e, STR_LIT("h"), odd, h, 10));
        for (int i = 0; i < 10; ++i) EXPECT_EQ(want[i], h[i]);
        ASSERT_EQ(1u, md_vial_eval_aggregate(e, STR_LIT("m"), odd, &m, 1));
        EXPECT_NEAR((float)(sum / n), m, 1e-4f);
        float x = 0;
        ASSERT_EQ(1u, md_vial_eval_aggregate(e, STR_LIT("x"), odd, &x, 1));
        EXPECT_NEAR((float)mx, x, 1e-4f);

        // Too small an out: the count, nothing written; not over the frames, or unknown: 0
        h[0] = -1.0f;
        EXPECT_EQ(10u, md_vial_eval_aggregate(e, STR_LIT("h"), odd, h, 4));
        EXPECT_EQ(-1.0f, h[0]);
        EXPECT_EQ(0u, md_vial_eval_aggregate(e, STR_LIT("d"), odd, h, 10));
        EXPECT_EQ(0u, md_vial_eval_aggregate(e, STR_LIT("nope"), odd, h, 10));
        md_vial_eval_free(e);
    }
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

static uint32_t vial_test_count_ca(const md_system_t* sys) {
    uint32_t n = 0;
    for (size_t i = 0; i < sys->atom.count; ++i) n += str_eq(md_atom_name(&sys->atom, i), STR_LIT("CA"));
    return n;
}

// What hover draws: the drawing shown() attached, kept through names and calls; built-ins that draw themselves; and
// values by their type
UTEST(vial_ir, draw) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    const int64_t nf = vial_test_run('P');
    ASSERT_GT(nf, 3);
    const md_system_state_t* st = vial_test_state('P');
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const char* src =
        "ca = particles(name='CA')\n"
        "d  = distance(ca[0], ca[1])\n"
        "ang = angle(ca[0], ca[1], ca[2])\n"
        "tor = dihedral(ca[0], ca[1], ca[2], ca[3])\n"
        "w  = within(4.0, of=ca[0])\n"
        "p  = ca.xyz\n"
        "n  = count(ca)\n"
        "both = line(ca[0], ca[1], label='x') + sphere(ca[2], 2.0)\n"
        "rs = residues()\n"
        "dc = distance(rs[0:3], ca[4], between='closest')\n"
        "dm = distance(rs[:, None], rs, between='closest')\n"
        "def ends(s: set):\n"
        "    u = s[0]\n"
        "    v = s[-1]\n"
        "    q = distance(u, v)\n"
        "    return shown(q, marker(u) + marker(v))\n"
        "e = ends(ca)\n"
        "f = e * 2.0\n";
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    ASSERT_TRUE(prog != NULL);
    {
        const md_vial_diagnostic_t* dg;
        const size_t nd = md_vial_diagnostics(prog, &dg);
        for (size_t i = 0; i < nd; ++i) if (dg[i].severity == MD_VIAL_SEVERITY_ERROR) printf("%s %.*s\n", dg[i].code, (int)dg[i].message.len, dg[i].message.ptr);
    }
    #define AT(text, plus) ((uint32_t)(strstr(src, text) - src) + (plus))
    md_vial_draw_desc_t desc = {.prog = prog, .sys = sys, .run = STR_LIT("run"), .frame = 0, .state = st, .element = -1};
    md_vial_drawing_t dr;

    // A set: its particles, and they are the CA atoms
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("particles(", 9)), alloc));
    ASSERT_EQ(1u, dr.num_groups);
    const uint32_t* ca = dr.atoms + dr.groups[0].beg;
    const uint32_t nca = dr.groups[0].end - dr.groups[0].beg;
    ASSERT_GT(nca, 4u);
    const float* x = (const float*)st->xyz;
    #define XYZ(atom, k) x[3 * (atom) + (k)]

    // distance: the prelude's shown line, labelled with the distance in Å
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("distance(ca[0]", 8)), alloc));
    ASSERT_EQ(1u, dr.num_prims);
    EXPECT_EQ((int)MD_VIAL_PRIM_LINE, (int)dr.prims[0].kind);
    for (int k = 0; k < 3; ++k) { EXPECT_EQ(XYZ(ca[0], k), dr.prims[0].p[0][k]); EXPECT_NEAR(XYZ(ca[1], k), dr.prims[0].p[1][k], 1e-4f); }
    ASSERT_EQ(1u, dr.num_labels);
    float len = 0;
    for (int k = 0; k < 3; ++k) len += (dr.prims[0].p[1][k] - dr.prims[0].p[0][k]) * (dr.prims[0].p[1][k] - dr.prims[0].p[0][k]);
    EXPECT_NEAR(sqrtf(len), dr.labels[0].value, 1e-4);
    EXPECT_FALSE(md_unit_is_none(dr.labels[0].unit));
    EXPECT_GT(dr.labels[0].text.len, 4u);
    // The name it is bound to draws the same
    md_vial_drawing_t by_name;
    ASSERT_TRUE(md_vial_draw(&by_name, &desc, md_vial_expr_at(prog, AT("d  =", 0)), alloc));
    EXPECT_EQ(1u, by_name.num_prims);

    // Marks: every expression with a type, outer before inner. ca[0] is marked apart from ca, which starts where it
    // does; drawing the call's mark draws the line, and the value beside it is the distance
    {
        const md_vial_mark_t* mk;
        const size_t nm = md_vial_marks(prog, &mk);
        ASSERT_GT(nm, 10u);
        for (size_t i = 1; i < nm; ++i) EXPECT_TRUE(mk[i - 1].span.beg < mk[i].span.beg || (mk[i - 1].span.beg == mk[i].span.beg && mk[i - 1].span.end >= mk[i].span.end));
        const uint32_t at = AT("ca[0], ca[1])", 0);
        size_t sub = SIZE_MAX, name = SIZE_MAX, call = SIZE_MAX;
        for (size_t i = 0; i < nm; ++i) {
            if (mk[i].span.beg == AT("distance(ca[0]", 0)) call = i;
            if (mk[i].span.beg != at) continue;
            if (sub == SIZE_MAX) sub = i; else if (name == SIZE_MAX) name = i;
        }
        ASSERT_TRUE(sub != SIZE_MAX && name != SIZE_MAX && call != SIZE_MAX);
        EXPECT_EQ(5u, mk[sub].span.end - mk[sub].span.beg);
        EXPECT_EQ(2u, mk[name].span.end - mk[name].span.beg);
        EXPECT_EQ(mk[sub].depth + 1, mk[name].depth);
        EXPECT_EQ(mk[call].depth + 1, mk[sub].depth);
        EXPECT_EQ((uint32_t)strlen("distance(ca[0], ca[1])"), mk[call].span.end - mk[call].span.beg);
        md_vial_drawing_t dm;
        ASSERT_TRUE(md_vial_draw(&dm, &desc, mk[call].expr, alloc));
        EXPECT_EQ(1u, dm.num_prims);
        md_vial_hover_t info;
        ASSERT_TRUE(md_vial_expr_info(prog, mk[call].expr, &info));
        EXPECT_TRUE(strstr(info.type, "float") != NULL);
        char vt[64];
        ASSERT_GT(md_vial_value_text(vt, sizeof(vt), &desc, mk[call].expr), 0u);
        EXPECT_NEAR(dm.labels[0].value, atof(vt), 1e-3);
        EXPECT_EQ(0u, md_vial_value_text(vt, sizeof(vt), &desc, 0));
    }

    // Moved coordinates move the drawing: values per frame come from the state given
    {
        md_system_state_t moved = {.alloc = alloc};
        md_system_state_copy(&moved, st);
        for (size_t i = 0; i < moved.num_atoms; ++i) moved.xyz[i].x += 10.0f;
        md_vial_draw_desc_t d2 = desc;
        d2.state = &moved;
        md_vial_drawing_t m;
        ASSERT_TRUE(md_vial_draw(&m, &d2, md_vial_expr_at(prog, AT("distance(ca[0]", 8)), alloc));
        ASSERT_EQ(1u, m.num_prims);
        EXPECT_NEAR(XYZ(ca[0], 0) + 10.0f, m.prims[0].p[0][0], 1e-4f);
    }

    // angle: an arc at the vertex; dihedral: the two planes and a label
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("angle(", 5)), alloc));
    ASSERT_EQ(1u, dr.num_prims);
    EXPECT_EQ((int)MD_VIAL_PRIM_ARC, (int)dr.prims[0].kind);
    for (int k = 0; k < 3; ++k) EXPECT_EQ(XYZ(ca[1], k), dr.prims[0].p[0][k]);
    EXPECT_EQ(1u, dr.num_labels);
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("dihedral(", 8)), alloc));
    ASSERT_EQ(3u, dr.num_prims);
    EXPECT_EQ((int)MD_VIAL_PRIM_DISK, (int)dr.prims[0].kind);
    EXPECT_EQ((int)MD_VIAL_PRIM_DISK, (int)dr.prims[1].kind);
    EXPECT_EQ((int)MD_VIAL_PRIM_ARC, (int)dr.prims[2].kind);
    EXPECT_EQ(1u, dr.num_labels);
    {
        // The arc's arms are at right angles to the axis, and the angle between them is the dihedral's
        const float* c = dr.prims[2].p[0];
        float ax[3], u[3], v[3];
        for (int k = 0; k < 3; ++k) { ax[k] = XYZ(ca[2], k) - XYZ(ca[1], k); u[k] = dr.prims[2].p[1][k] - c[k]; v[k] = dr.prims[2].p[2][k] - c[k]; }
        const float uu = sqrtf(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]), vv = sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        EXPECT_NEAR(0.0f, (u[0] * ax[0] + u[1] * ax[1] + u[2] * ax[2]) / uu, 1e-3f);
        EXPECT_NEAR(0.0f, (v[0] * ax[0] + v[1] * ax[1] + v[2] * ax[2]) / vv, 1e-3f);
        const double deg = acos((u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / (uu * vv)) * 180.0 / 3.14159265358979;
        EXPECT_NEAR(fabs(dr.labels[0].value), deg, 0.05);
    }

    // within: a sphere around the point, the particles in it highlighted
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("within(", 6)), alloc));
    ASSERT_EQ(1u, dr.num_prims);
    EXPECT_EQ((int)MD_VIAL_PRIM_SPHERE, (int)dr.prims[0].kind);
    EXPECT_EQ(4.0f, dr.prims[0].r);
    ASSERT_EQ(1u, dr.num_groups);
    bool has = false;
    for (uint32_t i = dr.groups[0].beg; i < dr.groups[0].end; ++i) has |= dr.atoms[i] == ca[0];
    EXPECT_TRUE(has);

    // Points as markers; a number draws nothing
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT(".xyz", 1)), alloc));
    EXPECT_EQ((size_t)nca, dr.num_prims);
    EXPECT_EQ((size_t)nca, dr.num_elements);
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("count(", 5)), alloc));
    EXPECT_EQ(0u, dr.num_prims + dr.num_labels + dr.num_groups);

    // Primitives combined, a text label
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("both", 0)), alloc));
    ASSERT_EQ(2u, dr.num_prims);
    EXPECT_EQ((int)MD_VIAL_PRIM_LINE, (int)dr.prims[0].kind);
    EXPECT_EQ((int)MD_VIAL_PRIM_SPHERE, (int)dr.prims[1].kind);
    ASSERT_EQ(1u, dr.num_labels);
    EXPECT_TRUE(str_eq(dr.labels[0].text, STR_LIT("x")));
    EXPECT_TRUE(isnan(dr.labels[0].value));

    // A def mapped over groups: one line per element, and one element alone
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("distance(rs", 8)), alloc));
    EXPECT_EQ(3u, dr.num_elements);
    EXPECT_EQ(3u, dr.num_prims);
    for (size_t i = 0; i < dr.num_prims; ++i) EXPECT_EQ((uint32_t)i, dr.prims[i].element);
    desc.element = 1;
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("distance(rs", 8)), alloc));
    ASSERT_EQ(1u, dr.num_prims);
    EXPECT_EQ(1u, dr.prims[0].element);
    desc.element = -1;

    // A map over the pairs of an array with itself: one line per pair. An element of a drawn value is drawn by its
    // type, so no call carries the whole array's drawing on
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("distance(rs[:", 8)), alloc));
    {
        const size_t nr = md_system_residue_count(sys);
        EXPECT_EQ(nr * nr, dr.num_elements);
        EXPECT_EQ(MIN(nr * nr, (size_t)5000), dr.num_prims);
    }

    // A script's def: what it attached, through the call and the name; an operation on it is drawn by its type
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("ends(ca)", 4)), alloc));
    ASSERT_EQ(2u, dr.num_prims);
    for (int k = 0; k < 3; ++k) { EXPECT_EQ(XYZ(ca[0], k), dr.prims[0].p[0][k]); EXPECT_EQ(XYZ(ca[nca - 1], k), dr.prims[1].p[0][k]); }
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("e = ", 0)), alloc));
    EXPECT_EQ(2u, dr.num_prims);
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("e * 2", 2)), alloc));
    EXPECT_EQ(0u, dr.num_prims);

    // A frame of the run, read: the line is as long as the distance published for that frame
    {
        md_vial_eval_t* ev = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
        ASSERT_TRUE(ev != NULL);
        EXPECT_TRUE(md_vial_eval_run(ev));
        const md_attribute_t* ad = md_attributes_find(md_vial_eval_attributes(ev), STR_LIT("script/d"));
        ASSERT_TRUE(ad != NULL);
        md_vial_draw_desc_t d3 = {.prog = prog, .sys = sys, .run = STR_LIT("run"), .frame = 3, .eval = ev, .element = -1};
        ASSERT_TRUE(md_vial_draw(&dr, &d3, md_vial_expr_at(prog, AT("distance(ca[0]", 8)), alloc));
        ASSERT_EQ(1u, dr.num_labels);
        EXPECT_NEAR(((const float*)ad->data)[3], dr.labels[0].value, 1e-5);
        md_vial_eval_free(ev);
    }
    #undef XYZ
    #undef AT
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

// Completion at the cursor, written | in the source
static md_vial_completions_t vial_ir__complete(const char* marked, md_allocator_i* alloc) {
    const char* bar = strchr(marked, '|');
    char* src = md_alloc(alloc, strlen(marked));
    const size_t at = (size_t)(bar - marked);
    memcpy(src, marked, at);
    strcpy(src + at, bar + 1);
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &(md_vial_compile_desc_t){.sys = vial_test_system('P'), .run = STR_LIT("run")}, alloc);
    md_vial_completions_t c;
    md_vial_complete(&c, prog, (uint32_t)at, alloc);
    // Strings are alloc's: the program may go
    md_vial_program_free(prog);
    return c;
}

static const md_vial_completion_t* vial_ir__find(const md_vial_completions_t* c, const char* label) {
    for (size_t i = 0; i < c->count; ++i) if (str_eq(c->items[i].label, str_from_cstr(label))) return &c->items[i];
    return NULL;
}

UTEST(vial_ir, editor) {
    ASSERT_TRUE(vial_test_system('P') != NULL);
    ASSERT_GT(vial_test_run('P'), 0);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));

    // Lists for highlighting, and names
    EXPECT_EQ(35u, md_vial_keywords(NULL, 0));
    str_t names[512];
    const size_t nb = md_vial_builtin_names(names, 512);
    ASSERT_LT(nb, 512u);
    bool has_particles = false, has_distance = false, has_hidden = false;
    for (size_t i = 0; i < nb; ++i) {
        has_particles |= str_eq(names[i], STR_LIT("particles"));
        has_distance |= str_eq(names[i], STR_LIT("distance"));
        has_hidden |= names[i].len && names[i].ptr[0] == '_';
    }
    EXPECT_TRUE(has_particles && has_distance && !has_hidden);
    EXPECT_TRUE(md_vial_name_valid(STR_LIT("my_sel2")));
    EXPECT_TRUE(md_vial_name_valid(STR_LIT("_x")));
    EXPECT_FALSE(md_vial_name_valid(STR_LIT("def")));
    EXPECT_FALSE(md_vial_name_valid(STR_LIT("particles")));
    EXPECT_FALSE(md_vial_name_valid(STR_LIT("distance")));
    EXPECT_FALSE(md_vial_name_valid(STR_LIT("2a")));
    EXPECT_FALSE(md_vial_name_valid(STR_LIT("a-b")));
    EXPECT_TRUE(str_eq(STR_LIT("energy/total"), md_vial_attr_path(STR_LIT("run"), STR_LIT("run/energy/total"))));
    EXPECT_EQ(0u, md_vial_attr_path(STR_LIT("run"), STR_LIT("runs/x")).len);

    md_vial_completions_t c;
    const md_vial_completion_t* it;
    // A name being typed: what starts so, names before keywords
    c = vial_ir__complete("ca = particles(name='CA')\nu = ce|", alloc);
    EXPECT_TRUE(str_eq(c.prefix, STR_LIT("ce")));
    EXPECT_TRUE(vial_ir__find(&c, "centroid") != NULL);
    for (size_t i = 0; i < c.count; ++i) EXPECT_TRUE(c.items[i].label.len >= 2 && !strncmp(c.items[i].label.ptr, "ce", 2));
    c = vial_ir__complete("a1 = 1\nb = a|\na2 = 2\n", alloc);
    ASSERT_TRUE((it = vial_ir__find(&c, "a1")) != NULL);
    EXPECT_EQ((int)MD_VIAL_COMPLETION_VARIABLE, (int)it->kind);
    EXPECT_TRUE(vial_ir__find(&c, "a2") == NULL);     // bound after
    EXPECT_TRUE(vial_ir__find(&c, "b") == NULL);      // being bound
    ASSERT_TRUE((it = vial_ir__find(&c, "and")) != NULL);
    EXPECT_EQ((int)MD_VIAL_COMPLETION_KEYWORD, (int)it->kind);
    EXPECT_EQ((int)MD_VIAL_COMPLETION_VARIABLE, (int)c.items[0].kind);
    // In a def: its parameters and locals; not outside it
    c = vial_ir__complete("def f(s: set, scale=1.0):\n    inner = count(s)\n    return in|\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "inner") != NULL);
    EXPECT_TRUE(vial_ir__find(&c, "in") != NULL);
    c = vial_ir__complete("def f(s: set, scale=1.0):\n    inner = count(s)\n    return inner\nq = in|\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "inner") == NULL);
    c = vial_ir__complete("def f(s: set, scale=1.0):\n    return s|\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "scale") != NULL);

    // The parameters of the call around the cursor: a prelude def's, a script def's, a built-in's; not one given
    c = vial_ir__complete("ca = particles(name='CA')\nd = distance(ca[0], ca[1], |)\n", alloc);
    ASSERT_TRUE((it = vial_ir__find(&c, "between")) != NULL);
    EXPECT_TRUE(str_eq(it->text, STR_LIT("between=")));
    EXPECT_EQ((int)MD_VIAL_COMPLETION_PARAMETER, (int)c.items[0].kind);
    c = vial_ir__complete("def f(s: set, scale=1.0):\n    return count(s) * scale\nx = f(particles(), sc|)\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "scale") != NULL && vial_ir__find(&c, "scale")->kind == MD_VIAL_COMPLETION_PARAMETER);
    c = vial_ir__complete("x = histogram(particles().x, axis='frame', |)\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "bins") != NULL);
    EXPECT_TRUE(vial_ir__find(&c, "axis") == NULL);
    c = vial_ir__complete("ca = particles(name='CA')\nd = angle(ca[0], ca[1], |)\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "c") != NULL);

    // A selector's keywords, those of its level
    c = vial_ir__complete("x = particles(res|)\n", alloc);
    ASSERT_TRUE((it = vial_ir__find(&c, "resname")) != NULL);
    EXPECT_EQ((int)MD_VIAL_COMPLETION_PARAMETER, (int)it->kind);
    c = vial_ir__complete("x = residues(|)\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "resname") == NULL);
    EXPECT_TRUE(vial_ir__find(&c, "id") != NULL);

    // After a dot: properties
    c = vial_ir__complete("ca = particles(name='CA')\nd = ca.c|\n", alloc);
    ASSERT_TRUE((it = vial_ir__find(&c, "com")) != NULL);
    EXPECT_EQ((int)MD_VIAL_COMPLETION_PROPERTY, (int)it->kind);
    EXPECT_TRUE(vial_ir__find(&c, "count") == NULL || vial_ir__find(&c, "count")->kind == MD_VIAL_COMPLETION_PROPERTY);

    // In a string: the system's values for the keyword, and the run's attributes for attr()
    c = vial_ir__complete("x = particles(resname='|')\n", alloc);
    ASSERT_TRUE((it = vial_ir__find(&c, "ALA")) != NULL);
    EXPECT_EQ((int)MD_VIAL_COMPLETION_VALUE, (int)it->kind);
    c = vial_ir__complete("x = particles(name='C|')\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "CA") != NULL && vial_ir__find(&c, "CB") != NULL);
    EXPECT_TRUE(vial_ir__find(&c, "N") == NULL);
    EXPECT_TRUE(str_eq(c.prefix, STR_LIT("C")));
    c = vial_ir__complete("x = particles('C|')\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "CA") != NULL);
    c = vial_ir__complete("x = residues(kind='a|')\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "amino_acid") != NULL);
    c = vial_ir__complete("w = attr('|')\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "time") != NULL);
    md_arena_allocator_destroy(alloc);
}

// Stored selections read by name, and filters: an expression compiled alone, with names
UTEST(vial_ir, filters) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    const md_system_state_t* st = vial_test_state('P');
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const uint32_t mine[] = {2, 0, 1, 1};     // out of order, once twice: the set is 0, 1, 2
    const md_vial_named_set_t stored[] = {{STR_LIT("mysel"), mine, 4}, {STR_LIT("bad name"), mine, 1}};

    // In a script: a static set, counted when compiling; its atoms are in its fingerprint
    {
        const char* src = "n = count(mysel)\nm = mysel.com\n";
        md_vial_program_t* p = md_vial_compile(str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .sets = stored, .num_sets = 2}, alloc);
        md_vial_hover_t h;
        ASSERT_TRUE(md_vial_hover(p, 0, &h));
        EXPECT_STREQ("3", h.value);
        const uint32_t other[] = {0, 1, 3};
        const md_vial_named_set_t moved[] = {{STR_LIT("mysel"), other, 3}};
        md_vial_program_t* q = md_vial_compile(str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .sets = moved, .num_sets = 1}, alloc);
        const md_vial_output_t* op;
        const md_vial_output_t* oq;
        md_vial_outputs(p, &op);
        md_vial_outputs(q, &oq);
        EXPECT_TRUE(op[1].fingerprint != oq[1].fingerprint);
        // A name the script binds hides it
        md_vial_program_t* r = md_vial_compile(STR_LIT("mysel = particles(name='CA')\nn = count(mysel)\n"), &(md_vial_compile_desc_t){.sys = sys, .sets = stored, .num_sets = 1}, alloc);
        ASSERT_TRUE(md_vial_hover(r, (uint32_t)strlen("mysel = particles(name='CA')\n"), &h));
        EXPECT_STRNE("3", h.value);
        md_vial_program_free(p); md_vial_program_free(q); md_vial_program_free(r);
    }

    md_vial_program_t* prog = md_vial_compile(STR_LIT("ca_s = particles(name='CA')\nnear = within(3.0, of=ca_s)\n"), &(md_vial_compile_desc_t){.sys = sys}, alloc);
    md_vial_filter_desc_t d = {.sys = sys, .state = st, .sets = stored, .num_sets = 1, .prog = prog};
    md_vial_filter_t f;
    // A predefined set, a selector's array, a stored selection, the program's static selection
    ASSERT_TRUE(md_vial_filter(&f, STR_LIT("protein"), &d, alloc));
    EXPECT_FALSE(f.dynamic);
    EXPECT_EQ(0, f.ndim);
    EXPECT_EQ(1u, f.num_sets);
    EXPECT_GT(f.offsets[1], 0u);
    ASSERT_TRUE(md_vial_filter(&f, STR_LIT("residues()"), &d, alloc));
    EXPECT_EQ(1, f.ndim);
    EXPECT_EQ(md_system_residue_count(sys), (size_t)f.dims[0]);
    EXPECT_EQ((size_t)f.dims[0], f.num_sets);
    ASSERT_TRUE(md_vial_filter(&f, STR_LIT("mysel"), &d, alloc));
    ASSERT_EQ(3u, f.offsets[1]);
    EXPECT_EQ(0u, f.atoms[0]); EXPECT_EQ(1u, f.atoms[1]); EXPECT_EQ(2u, f.atoms[2]);
    ASSERT_TRUE(md_vial_filter(&f, STR_LIT("ca_s and not mysel"), &d, alloc));
    const uint32_t nca = f.offsets[1];
    EXPECT_GT(nca, 0u);
    // Bonds: pairs of one-particle sets
    ASSERT_TRUE(md_vial_filter(&f, STR_LIT("bonds(residues()[0])"), &d, alloc));
    EXPECT_EQ(2, f.ndim);
    EXPECT_EQ(2u, f.dims[1]);
    EXPECT_GT(f.dims[0], 0u);
    for (size_t i = 0; i < f.num_sets; ++i) EXPECT_EQ(1u, f.offsets[i + 1] - f.offsets[i]);
    // Dynamic: from the coordinates, and none without them
    ASSERT_TRUE(md_vial_filter(&f, STR_LIT("within(3.0, of=mysel)"), &d, alloc));
    EXPECT_TRUE(f.dynamic);
    EXPECT_GE(f.offsets[1], 3u);
    md_vial_filter_desc_t nd = d;
    nd.state = NULL;
    EXPECT_FALSE(md_vial_filter(&f, STR_LIT("within(3.0, of=mysel)"), &nd, alloc));
    // A recipe of the program's: evaluated at the coordinates, so dynamic; without them it is not there
    ASSERT_TRUE(md_vial_filter(&f, STR_LIT("near and not ca_s"), &d, alloc));
    EXPECT_TRUE(f.dynamic);
    EXPECT_FALSE(md_vial_filter(&f, STR_LIT("near"), &nd, alloc));
    EXPECT_FALSE(md_vial_filter(&f, STR_LIT("nearby"), &d, alloc));
    // Arrays of sets the program publishes: the name gives the array, and a subscript one of it
    {
        md_vial_program_t* pa = md_vial_compile(STR_LIT("res = residues()\nca_r = residues() and particles(name='CA')\nbb = bonds(residues()[0])\n"),
            &(md_vial_compile_desc_t){.sys = sys}, alloc);
        md_vial_filter_desc_t ad = {.sys = sys, .state = st, .prog = pa};
        ASSERT_TRUE(md_vial_filter(&f, STR_LIT("res"), &ad, alloc));
        EXPECT_EQ(1, f.ndim);
        EXPECT_EQ(md_system_residue_count(sys), (size_t)f.dims[0]);
        md_vial_filter_t r3;
        ASSERT_TRUE(md_vial_filter(&r3, STR_LIT("residues()[3]"), &ad, alloc));
        ASSERT_TRUE(md_vial_filter(&f, STR_LIT("res[3]"), &ad, alloc));
        EXPECT_EQ(0, f.ndim);
        ASSERT_EQ(r3.offsets[1], f.offsets[1]);
        for (uint32_t i = 0; i < f.offsets[1]; ++i) EXPECT_EQ(r3.atoms[i], f.atoms[i]);
        ASSERT_TRUE(md_vial_filter(&f, STR_LIT("ca_r[2]"), &ad, alloc));
        EXPECT_EQ(1u, f.offsets[1]);
        ASSERT_TRUE(md_vial_filter(&f, STR_LIT("bb"), &ad, alloc));
        EXPECT_EQ(2, f.ndim);
        EXPECT_EQ(2u, f.dims[1]);
        ASSERT_TRUE(md_vial_filter(&f, STR_LIT("any(res[0:2])"), &ad, alloc));
        EXPECT_EQ(0, f.ndim);

        // As a stored selection given to a script: an array, known when compiling
        const uint32_t at[] = {0, 1, 2, 5, 4};
        const uint32_t off[] = {0, 3, 5};
        const md_vial_named_set_t arr[] = {{.name = STR_LIT("pair"), .atoms = at, .num_atoms = 5, .ndim = 1, .dims = {2}, .offsets = off}};
        md_vial_program_t* q = md_vial_compile(STR_LIT("n = len(pair)\nm = count(pair[1])\n"), &(md_vial_compile_desc_t){.sys = sys, .sets = arr, .num_sets = 1}, alloc);
        md_vial_hover_t h;
        ASSERT_TRUE(md_vial_hover(q, 0, &h));
        EXPECT_STREQ("2", h.value);
        ASSERT_TRUE(md_vial_hover(q, (uint32_t)strlen("n = len(pair)\n"), &h));
        EXPECT_STREQ("2", h.value);
        md_vial_program_free(q);
        md_vial_program_free(pa);
    }
    // Within a scope: cut to it, emptied sets dropped
    {
        md_vial_filter_t r0;
        ASSERT_TRUE(md_vial_filter(&r0, STR_LIT("residues()[0]"), &d, alloc));
        md_vial_filter_desc_t sd = d;
        sd.scope = r0.atoms;
        sd.num_scope = r0.offsets[1];
        ASSERT_TRUE(md_vial_filter(&f, STR_LIT("residues()"), &sd, alloc));
        EXPECT_EQ(1u, f.dims[0]);
        EXPECT_EQ(r0.offsets[1], f.offsets[1]);
        ASSERT_TRUE(md_vial_filter(&f, STR_LIT("protein"), &sd, alloc));
        EXPECT_EQ(r0.offsets[1], f.offsets[1]);
    }
    // Errors, with spans in the expression
    EXPECT_FALSE(md_vial_filter(&f, STR_LIT("1 + 2"), &d, alloc));
    EXPECT_STREQ("E0313", f.error.code);
    EXPECT_EQ(0u, f.error.span.beg);
    EXPECT_EQ(5u, f.error.span.end);
    EXPECT_FALSE(md_vial_filter(&f, STR_LIT("particles(nme='CA')"), &d, alloc));
    EXPECT_GE(f.error.span.beg, 0u);
    EXPECT_LE(f.error.span.end, 19u);
    EXPECT_GT(f.error.message.len, 0u);
    EXPECT_FALSE(md_vial_filter(&f, STR_LIT("   "), &d, alloc));
    EXPECT_FALSE(md_vial_filter(&f, STR_LIT("protein)"), &d, alloc));

    // Completion, in the expression's own offsets
    md_vial_completions_t c;
    ASSERT_TRUE(md_vial_filter_complete(&c, STR_LIT("particles(resn"), 14, &d, alloc));
    EXPECT_TRUE(vial_ir__find(&c, "resname") != NULL);
    EXPECT_EQ(10u, c.span.beg);
    EXPECT_EQ(14u, c.span.end);
    ASSERT_TRUE(md_vial_filter_complete(&c, STR_LIT("mys"), 3, &d, alloc));
    EXPECT_TRUE(vial_ir__find(&c, "mysel") != NULL);
    ASSERT_TRUE(md_vial_filter_complete(&c, STR_LIT("ca"), 2, &d, alloc));
    EXPECT_TRUE(vial_ir__find(&c, "ca_s") != NULL);
    ASSERT_TRUE(md_vial_filter_complete(&c, STR_LIT("_"), 1, &d, alloc));
    EXPECT_TRUE(vial_ir__find(&c, "_filter") == NULL);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

// A periodic output publishes its period and seam beside it, from when the evaluation is made, as v1 does
// The predefined sets of the structure's anatomy: a lipid's head and tails, a nucleotide's base and nucleoside
static bool vial_ir__load_file(md_system_t* sys, md_system_state_t* st, const char* file, md_allocator_i* alloc) {
    char path[512];
    snprintf(path, sizeof(path), "%s/%s", MD_UNITTEST_DATA_DIR, file);
    *sys = (md_system_t){.alloc = alloc};
    *st = (md_system_state_t){.alloc = alloc};
    const str_t p = str_from_cstr(path);
    const bool ok = strstr(file, ".gro") ? md_gro_system_init_from_file(sys, st, p) : md_pdb_system_init_from_file(sys, st, p, MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE);
    return ok && md_util_system_infer(sys, st, MD_UTIL_INFER_ALL);
}

static uint32_t vial_ir__count(const char* expr, const md_system_t* sys, const md_system_state_t* st, md_allocator_i* alloc) {
    md_vial_filter_t f;
    if (!md_vial_filter(&f, str_from_cstr(expr), &(md_vial_filter_desc_t){.sys = sys, .state = st}, alloc) || f.num_sets != 1) return UINT32_MAX;
    return f.offsets[1];
}

UTEST(vial_ir, anatomy_sets) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(16));
    md_system_t sys;
    md_system_state_t st;
    ASSERT_TRUE(vial_ir__load_file(&sys, &st, "dppc64.pdb", alloc));
    const uint32_t lipid = vial_ir__count("lipid", &sys, &st, alloc);
    const uint32_t head = vial_ir__count("particles(heads)", &sys, &st, alloc), tail = vial_ir__count("particles(tails)", &sys, &st, alloc);
    EXPECT_GT(head, 0u);
    EXPECT_GT(tail, 0u);
    EXPECT_EQ(lipid, head + tail);
    EXPECT_EQ(0u, vial_ir__count("particles(heads) and particles(tails)", &sys, &st, alloc));
    EXPECT_EQ(lipid, vial_ir__count("particles(heads) or particles(tails)", &sys, &st, alloc));
    EXPECT_EQ(UINT32_MAX, vial_ir__count("head", &sys, &st, alloc));     // units, not kinds: no head set
    // One head per DPPC, with its hydrogens
    {
        md_vial_filter_t h;
        ASSERT_TRUE(md_vial_filter(&h, STR_LIT("heads"), &(md_vial_filter_desc_t){.sys = &sys, .state = &st}, alloc));
        EXPECT_EQ(1, h.ndim);
        EXPECT_EQ(sys.lipid.head.count, (size_t)h.num_sets);
        EXPECT_EQ(sys.lipid.count, (size_t)h.num_sets);
        EXPECT_EQ(head, h.offsets[h.num_sets]);
    }
    // One set per tail; all of them are the tail particles
    md_vial_filter_t f;
    ASSERT_TRUE(md_vial_filter(&f, STR_LIT("tails"), &(md_vial_filter_desc_t){.sys = &sys, .state = &st}, alloc));
    EXPECT_EQ(1, f.ndim);
    EXPECT_EQ(sys.lipid.tail.count, (size_t)f.num_sets);
    EXPECT_EQ(tail, f.offsets[f.num_sets]);
    EXPECT_EQ(tail, vial_ir__count("particles(tails)", &sys, &st, alloc));
    ASSERT_TRUE(md_vial_filter(&f, STR_LIT("tails(residues[0])"), &(md_vial_filter_desc_t){.sys = &sys, .state = &st}, alloc));
    EXPECT_EQ(2u, (unsigned)f.num_sets);     // DPPC: two acyl chains
    {
        md_vial_program_t* p = md_vial_compile(STR_LIT("l = labels(tails)[1]\n"), &(md_vial_compile_desc_t){.sys = &sys}, alloc);
        md_vial_hover_t h;
        ASSERT_TRUE(md_vial_hover(p, 0, &h));
        EXPECT_TRUE(strstr(h.value, "tail2/") != NULL);
        md_vial_program_free(p);
    }

    ASSERT_TRUE(vial_ir__load_file(&sys, &st, "nucl-dna.gro", alloc));
    const uint32_t base = vial_ir__count("nucleobase", &sys, &st, alloc), side = vial_ir__count("nucleoside", &sys, &st, alloc);
    EXPECT_GT(base, 0u);
    EXPECT_GT(side, base);
    EXPECT_EQ(base, vial_ir__count("nucleobase and nucleoside", &sys, &st, alloc));
    EXPECT_EQ(0u, vial_ir__count("nucleoside and not nucleic", &sys, &st, alloc));
    md_arena_allocator_destroy(alloc);
}

// A mask over a set's particles picks those where it is true, per frame: as many as the mask counts
UTEST(vial_ir, set_mask) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    ASSERT_GT(vial_test_run('P'), 3);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(8));
    const char* src =
        "n = count(particles[particles.z > particles[0].z])\n"
        "m = sum(particles.z > particles[0].z)\n"
        "k = count(protein[(protein.x in -5.0..5.0) and (protein.y > 0.0)])\n"
        "j = sum((protein.x in -5.0..5.0) and (protein.y > 0.0))\n";
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    ASSERT_EQ(0u, (unsigned)md_vial_diagnostics(prog, NULL));
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    ASSERT_TRUE(e && md_vial_eval_run(e));
    const md_attributes_t* at = md_vial_eval_attributes(e);
    const md_attribute_t* n = md_attributes_find(at, STR_LIT("script/n"));
    const md_attribute_t* m = md_attributes_find(at, STR_LIT("script/m"));
    const md_attribute_t* k = md_attributes_find(at, STR_LIT("script/k"));
    const md_attribute_t* j = md_attributes_find(at, STR_LIT("script/j"));
    ASSERT_TRUE(n && m && k && j);
    const size_t nf = n->format.shape[0];
    int varies = 0;
    for (size_t f = 0; f < nf; ++f) {
        EXPECT_EQ(((const int32_t*)m->data)[f], ((const int32_t*)n->data)[f]);
        EXPECT_EQ(((const int32_t*)j->data)[f], ((const int32_t*)k->data)[f]);
        varies |= ((const int32_t*)n->data)[f] != ((const int32_t*)n->data)[0];
    }
    EXPECT_TRUE(varies);
    // As a filter: dynamic
    md_vial_filter_t fl;
    ASSERT_TRUE(md_vial_filter(&fl, STR_LIT("particles[particles.z > 0.0]"), &(md_vial_filter_desc_t){.sys = sys, .state = vial_test_state('P')}, alloc));
    EXPECT_TRUE(fl.dynamic);
    md_vial_eval_free(e);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

// A per-frame population has its summary over the members published beside it, as v1's: /mean, /variance, /extent;
// a periodic one around its circle. Its members are labelled by the units its dimension runs over.
UTEST(vial_ir, population_summaries) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    ASSERT_GT(vial_test_run('P'), 3);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(8));
    const char* src =
        "ca  = particles(name='CA')\n"
        "d   = distance(ca[0], ca.xyz)\n"
        "psi = dihedral(residues() and particles(name='N'), residues() and particles(name='CA'), residues() and particles(name='C'), residues() and particles(name='O'))\n"
        "one = distance(ca[0], ca[1])\n";
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    ASSERT_EQ(0u, (unsigned)md_vial_diagnostics(prog, NULL));
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    ASSERT_TRUE(e != NULL);
    ASSERT_TRUE(md_vial_eval_run(e));
    const md_attributes_t* at = md_vial_eval_attributes(e);
    const md_attribute_t* d = md_attributes_find(at, STR_LIT("script/d"));
    const md_attribute_t* dm = md_attributes_find(at, STR_LIT("script/d/mean"));
    const md_attribute_t* dv = md_attributes_find(at, STR_LIT("script/d/variance"));
    const md_attribute_t* dx = md_attributes_find(at, STR_LIT("script/d/extent"));
    ASSERT_TRUE(d && dm && dv && dx);
    EXPECT_TRUE(md_attributes_find(at, STR_LIT("script/one/mean")) == NULL);    // one value has no population
    const size_t nf = d->format.shape[0], n = d->format.shape[1];
    ASSERT_EQ(15u, (unsigned)n);
    EXPECT_EQ(nf, (size_t)dm->format.shape[0]);
    EXPECT_EQ(2u, dx->format.components);
    for (size_t f = 0; f < nf; ++f) {
        const float* row = (const float*)d->data + f * n;
        double sum = 0, lo = row[0], hi = row[0];
        for (size_t i = 0; i < n; ++i) { sum += row[i]; lo = row[i] < lo ? row[i] : lo; hi = row[i] > hi ? row[i] : hi; }
        const double mean = sum / n;
        double var = 0;
        for (size_t i = 0; i < n; ++i) var += (row[i] - mean) * (row[i] - mean);
        EXPECT_NEAR(mean, ((const float*)dm->data)[f], 1e-3);
        EXPECT_NEAR(var / n, ((const float*)dv->data)[f], 1e-2);
        EXPECT_NEAR(lo, ((const float*)dx->data)[2 * f], 1e-5);
        EXPECT_NEAR(hi, ((const float*)dx->data)[2 * f + 1], 1e-5);
    }
    // Periodic: the mean lies in the period, the arc starts in it
    const md_attribute_t* pm = md_attributes_find(at, STR_LIT("script/psi/mean"));
    const md_attribute_t* pp = md_attributes_find(at, STR_LIT("script/psi/period"));
    const md_attribute_t* px = md_attributes_find(at, STR_LIT("script/psi/extent"));
    ASSERT_TRUE(pm && pp && px);
    const float beg = ((const float*)pp->data)[0], end = ((const float*)pp->data)[1];
    for (size_t f = 0; f < nf; ++f) {
        const float m = ((const float*)pm->data)[f], lo = ((const float*)px->data)[2 * f], hi = ((const float*)px->data)[2 * f + 1];
        EXPECT_TRUE(m >= beg && m < end);
        EXPECT_TRUE(lo >= beg && lo < end);
        EXPECT_TRUE(hi >= lo && hi - lo <= end - beg);
    }
    // The labels of the CAs, as the data names them
    const md_attribute_t* dl = md_attributes_find(at, STR_LIT("script/d/labels"));
    ASSERT_TRUE(dl != NULL);
    ASSERT_EQ(15u, dl->format.shape[0]);
    EXPECT_TRUE(str_eq(md_attribute_str(at, dl, 0), STR_LIT("CA/ALA1")));

    // An evaluation that takes d over keeps its summary
    md_vial_program_t* p2 = md_vial_compile(str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    md_vial_eval_t* e2 = md_vial_eval_create_from(e, p2, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    ASSERT_TRUE(e2 != NULL);
    const md_attribute_t* dm2 = md_attributes_find(md_vial_eval_attributes(e2), STR_LIT("script/d/mean"));
    ASSERT_TRUE(dm2 != NULL);
    dm = md_attributes_find(at, STR_LIT("script/d/mean"));
    EXPECT_EQ(0, memcmp(dm->data, dm2->data, nf * sizeof(float)));
    md_vial_eval_free(e2);
    md_vial_eval_free(e);
    md_vial_program_free(p2);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

UTEST(vial_ir, periodic_outputs) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    ASSERT_GT(vial_test_run('P'), 3);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const char* src =
        "ca = particles(name='CA')\n"
        "t  = dihedral(ca[0], ca[1], ca[2], ca[3])\n"
        "d  = distance(ca[0], ca[1])\n";
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    const md_vial_output_t* o;
    const size_t no = md_vial_outputs(prog, &o);
    for (size_t i = 0; i < no; ++i) {
        if (str_eq(o[i].name, STR_LIT("t"))) EXPECT_TRUE(o[i].flags & MD_VIAL_OUTPUT_PERIODIC);
        if (str_eq(o[i].name, STR_LIT("d"))) EXPECT_FALSE(o[i].flags & MD_VIAL_OUTPUT_PERIODIC);
    }
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    ASSERT_TRUE(e != NULL);
    const md_attributes_t* at = md_vial_eval_attributes(e);
    const md_attribute_t* per = md_attributes_find(at, STR_LIT("script/t/period"));
    const md_attribute_t* seam = md_attributes_find(at, STR_LIT("script/t/seam"));
    ASSERT_TRUE(per && seam);
    EXPECT_TRUE(md_attributes_find(at, STR_LIT("script/d/period")) == NULL);
    double p[2], sm;
    ASSERT_EQ(2u, md_attribute_extract_f64(p, 2, per, md_attribute_slice_all(), per->unit));
    ASSERT_EQ(1u, md_attribute_extract_f64(&sm, 1, seam, md_attribute_slice_all(), seam->unit));
    EXPECT_LT(p[0], p[1]);
    EXPECT_GE(sm, p[0]);
    EXPECT_LT(sm, p[1]);
    // In the unit the values are published in, and the values lie in the period
    EXPECT_TRUE(md_vial_eval_run(e));
    const md_attribute_t* vt = md_attributes_find(at, STR_LIT("script/t"));
    ASSERT_TRUE(vt != NULL);
    per = md_attributes_find(at, STR_LIT("script/t/period"));     // the table has grown: pointers into it move
    ASSERT_TRUE(per != NULL);
    EXPECT_TRUE(md_unit_equal(per->unit, vt->unit));
    const size_t nf = (size_t)vt->format.shape[0];
    for (size_t f = 0; f < nf; ++f) {
        const float v = ((const float*)vt->data)[f];
        EXPECT_TRUE(v >= p[0] - 1e-3 && v <= p[1] + 1e-3);
    }
    md_vial_eval_free(e);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

// R0101: a position outside a dynamic extent makes the value missing in that frame and is said once; nothing fails, and
// threads say the same as one
UTEST(vial_ir, missing_is_said) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    const int64_t nf = vial_test_run('P');
    ASSERT_GT(nf, 20);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    const char* src =
        "ca = particles(name='CA')\n"
        "x = distance(ca[0], ca[frame()])\n"
        "d = distance(ca[0], ca[1])\n";
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    for (int threads = 1; threads <= 4; threads += 3) {
        md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run"), .num_threads = threads}, alloc);
        ASSERT_TRUE(e != NULL);
        EXPECT_TRUE(md_vial_eval_run(e));
        const md_vial_diagnostic_t* dg;
        ASSERT_EQ(1u, md_vial_eval_diagnostics(e, &dg));
        EXPECT_STREQ("R0101", dg[0].code);
        EXPECT_EQ(MD_VIAL_SEVERITY_RUNTIME, dg[0].severity);
        const md_attribute_t* ax = md_attributes_find(md_vial_eval_attributes(e), STR_LIT("script/x"));
        ASSERT_TRUE(ax != NULL);
        // The first frame past the CA particles names it, and from there the value is missing
        char want[32];
        snprintf(want, sizeof(want), "frame %u:", (unsigned)(vial_test_count_ca(sys)));
        EXPECT_TRUE(strncmp(dg[0].message.ptr, want, strlen(want)) == 0);
        EXPECT_TRUE(isnan(((const float*)ax->data)[nf - 1]));
        EXPECT_FALSE(isnan(((const float*)ax->data)[1]));
        md_vial_eval_free(e);
    }
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

UTEST(vial_ir, contraction) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    // Fused: a product used once, floats throughout. Not fused: ints, and a product that is also published.
    md_vial_program_t* prog = vial_ir__compile("a = [1.5, 2.5]\nb = [3.0, 4.0]\nc = a * b + 1.0\nd = 2 * 3 + 1\ne = a * b\nf = e + 1.0\n", 0, alloc);
    const vial_ir_t* ir = md_vial__run_ir(prog);
    ASSERT_TRUE(ir != NULL);
    int fmas = 0;
    for (size_t i = 0; i < md_array_size(ir->ops); ++i) fmas += ir->ops[i].kind == VIAL_IR_FMA;
    EXPECT_EQ(1, fmas);
    char err[256];
    EXPECT_TRUE(vial_ir_verify(ir, err, sizeof(err)));
    vial_interp_t in;
    vial_interp_init(&in, ir, NULL, NULL, alloc);
    EXPECT_TRUE(vial_interp_static(&in));
    const vial_val_t* c = vial_interp_value(&in, vial_ir__output(ir, "c"));
    ASSERT_TRUE(c != NULL);
    EXPECT_EQ(fmaf(1.5f, 3.0f, 1.0f), ((const float*)c->data)[0]);
    EXPECT_EQ(fmaf(2.5f, 4.0f, 1.0f), ((const float*)c->data)[1]);
    vial_interp_free(&in);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

UTEST(vial_ir, if_per_frame) {
    const md_system_t* sys = vial_test_system('P');
    const int64_t nf = vial_test_run('P');
    ASSERT_TRUE(sys && nf > 2);
    vial_ir__run_t r;
    const char* src =
        "ca = particles(name='CA')\n"
        "def clamp(d):\n"
        "    if d > 30.0:\n"
        "        return 30.0\n"
        "    else:\n"
        "        return d\n"
        "def bump(d):\n"
        "    e = d\n"
        "    if d > 30.0:\n"
        "        e = d - 10.0\n"
        "    return e\n"
        // One branch returns, the other goes on: the rest of the body is that branch's
        "def early(d):\n"
        "    if d > 30.0:\n"
        "        return 30.0\n"
        "    e = d * 2.0\n"
        "    return e\n"
        "def late(d):\n"
        "    if d <= 30.0:\n"
        "        e = d + 1.0\n"
        "    else:\n"
        "        return 0.0\n"
        "    return e\n"
        "def chain(d):\n"
        "    if d > 40.0:\n"
        "        return 40.0\n"
        "    if d > 30.0:\n"
        "        return 30.0\n"
        "    return d\n"
        "ends = distance(ca[0], ca[14])\n"
        "c = clamp(ends)\n"
        "b = bump(ends)\n"
        "t = 30.0 if ends > 30.0 else ends\n"
        "m1 = early(ends)\n"
        "m2 = late(ends)\n"
        "m3 = chain(ends)\n";
    ASSERT_TRUE(vial_ir__start(&r, src, sys));
    vial_ir__loader_t l = {.sys = sys, .state = {.alloc = r.alloc}};
    md_system_state_init(&l.state, sys->atom.count);
    r.in.load_frame = vial_ir__load;
    r.in.load_user = &l;
    const bool ok = vial_interp_run(&r.in, nf);
    for (size_t i = 0; i < md_array_size(r.in.errors); ++i) printf("failed: %s %s\n", r.in.errors[i].code, r.in.errors[i].message);
    ASSERT_TRUE(ok);
    const vial_ir_t* ir = md_vial__run_ir(r.prog);
    int above = 0;
    for (int64_t f = 0; f < nf; ++f) {
        const double e = vial_ir__num(vial_interp_at(&r.in, vial_ir__output(ir, "ends"), f));
        above += e > 30.0;
        EXPECT_EQ((float)(e > 30.0 ? 30.0 : e), (float)vial_ir__num(vial_interp_at(&r.in, vial_ir__output(ir, "c"), f)));
        EXPECT_EQ(e > 30.0 ? (float)e - 10.0f : (float)e, (float)vial_ir__num(vial_interp_at(&r.in, vial_ir__output(ir, "b"), f)));
        EXPECT_EQ((float)(e > 30.0 ? 30.0 : e), (float)vial_ir__num(vial_interp_at(&r.in, vial_ir__output(ir, "t"), f)));
        EXPECT_EQ(e > 30.0 ? 30.0f : (float)e * 2.0f, (float)vial_ir__num(vial_interp_at(&r.in, vial_ir__output(ir, "m1"), f)));
        EXPECT_EQ(e <= 30.0 ? (float)e + 1.0f : 0.0f, (float)vial_ir__num(vial_interp_at(&r.in, vial_ir__output(ir, "m2"), f)));
        EXPECT_EQ(e > 40.0 ? 40.0f : e > 30.0 ? 30.0f : (float)e, (float)vial_ir__num(vial_interp_at(&r.in, vial_ir__output(ir, "m3"), f)));
    }
    EXPECT_GT(above, 0);
    EXPECT_LT(above, (int)nf);
    md_system_state_free(&l.state);
    vial_ir__end(&r);
}

// What the check folds is what evaluation gives, bit for bit: each a is known when compiling and its b, the same
// arithmetic on a value of the frame, is evaluated. Float is 32 bits in both; a bare number beside a float64 is one.
UTEST(vial_ir, folds_as_evaluated) {
    vial_ir__run_t r;
    const char* src =
        "one = unitless(particles[0].x * 0.0) + 1.0\n"
        "a1 = 0.1 + 0.2 - 0.3\n"
        "b1 = one * 0.1 + 0.2 - 0.3\n"
        "a2 = 1.1 ** 3\n"
        "b2 = (one * 1.1) ** 3\n"
        "a3 = 1.0 // 0.1\n"
        "b3 = one // 0.1\n"
        "a4 = 0.3 == 0.1 + 0.2\n"
        "b4 = one * 0.3 == one * 0.1 + 0.2\n"
        "a5 = 0.0 Å + 1.3 nm\n"
        "b5 = one * 0.0 Å + 1.3 nm\n"
        "a6 = int(2.9999999999)\n"
        "b6 = int(one * 2.9999999999)\n"
        "a7 = float(16777217)\n"
        "b7 = float(int(one) * 16777217)\n"
        "a8 = 7.0 % 0.7\n"
        "b8 = (one * 7.0) % 0.7\n";
    const md_system_state_t* st = vial_test_state('T');
    ASSERT_TRUE(st != NULL);
    ASSERT_TRUE(vial_ir__start(&r, src, vial_test_system('T')));
    const vial_frame_t fr = {.xyz = (const float*)st->xyz, .num_atoms = st->num_atoms, .cell = st->unitcell};
    EXPECT_TRUE(vial_interp_frame(&r.in, &fr));
    const vial_ir_t* ir = md_vial__run_ir(r.prog);
    for (int i = 1; i <= 8; ++i) {
        char an[8], bn[8], av[128], bv[128], wrapped[140];
        snprintf(an, sizeof(an), "a%d", i);
        snprintf(bn, sizeof(bn), "b%d", i);
        // a is a constant: folded, not computed
        const vial_ssa_t as = vial_ir__output(ir, an);
        ASSERT_TRUE(as != 0);
        EXPECT_EQ(VIAL_IR_CONST, ir->ops[ir->values[as].op].kind);
        vial_ir__out(&r, an, av, sizeof(av));
        vial_ir__out(&r, bn, bv, sizeof(bv));
        snprintf(wrapped, sizeof(wrapped), "[%s]", av);
        if (strcmp(wrapped, bv)) printf("%s = %s, but %s = %s\n", an, av, bn, bv);
        EXPECT_STREQ(wrapped, bv);
    }
    vial_ir__end(&r);
}

// A map stops between calls when asked, and a map's calls do not keep their intermediates: one over every pair of
// residues runs in the memory of its result and of one call
UTEST(vial_ir, map_interrupt_and_memory) {
    vial_ir__run_t r;
    const char* src =
        "res = residues(chain='A', kind='amino_acid')\n"
        "d = distance(res[:, None], res, between='closest')\n";
    const md_system_state_t* st = vial_test_state('T');
    ASSERT_TRUE(st != NULL);
    ASSERT_TRUE(vial_ir__start(&r, src, vial_test_system('T')));
    const vial_frame_t fr = {.xyz = (const float*)st->xyz, .num_atoms = st->num_atoms, .cell = st->unitcell};
    volatile uint64_t stop = 1;
    r.in.stop = &stop;
    // Batched for the built-in, and call by call
    for (int b = 0; b < 2; ++b) {
        r.in.plain = b == 1;
        r.in.interrupted = false;
        EXPECT_FALSE(vial_interp_frame(&r.in, &fr));
        EXPECT_TRUE(r.in.interrupted);
        EXPECT_EQ(0, (int)md_array_size(r.in.errors));
    }
    stop = 0;
    r.in.interrupted = false;
    EXPECT_TRUE(vial_interp_frame(&r.in, &fr));
    const vial_ir_t* ir = md_vial__run_ir(r.prog);
    const vial_val_t* d = vial_interp_value(&r.in, vial_ir__output(ir, "d"));
    ASSERT_TRUE(d != NULL);
    EXPECT_EQ(2, (int)d->ndim);
    EXPECT_EQ(d->dims[0], d->dims[1]);
    // The scratch arena of the map is back where it began
    ASSERT_TRUE(r.in.map_scratch[0] != NULL);
    const size_t pos = md_vm_arena_get_pos(r.in.map_scratch[0]);
    EXPECT_TRUE(vial_interp_frame(&r.in, &fr));
    EXPECT_EQ(pos, md_vm_arena_get_pos(r.in.map_scratch[0]));
    vial_ir__end(&r);
}

// A map whose calls are only a closest or farthest pair gives the built-in every call at once: the same bits as the
// calls one by one, on T and over P's frames
UTEST(vial_ir, pair_maps_batch) {
    const char* src =
        "res = residues(chain='A', kind='amino_acid')\n"
        "c = distance(res[:, None], res, between='closest')\n"
        "f = distance(res[0:40, None], res[10:60], between='farthest')\n"
        "e = distance(res[0:30], res[5:35], between='closest')\n";
    const char* src_p =
        "res = residues\n"
        "c = distance(res[:, None], res, between='closest')\n"
        "f = distance(res[:, None], res, between='farthest')\n";
    for (int pass = 0; pass < 2; ++pass) {
        const char structure = pass ? 'P' : 'T';
        int64_t run_frames = 0;
        const md_system_t* sys = pass ? vial_ir__run_p(&run_frames) : vial_test_system(structure);
        ASSERT_TRUE(sys != NULL);
        vial_ir__run_t r[2];
        for (int b = 0; b < 2; ++b) {
            ASSERT_TRUE(vial_ir__start(&r[b], pass ? src_p : src, sys));
            r[b].in.plain = b == 1;
        }
        vial_frame_t fr;
        vial_ir__loader_t l = {.sys = sys, .state = {.alloc = r[0].alloc}};
        const md_system_state_t* st = vial_test_state(structure);
        const int64_t nf = pass ? 3 : 1;
        if (pass) md_system_state_init(&l.state, sys->atom.count);
        for (int64_t f = 0; f < nf; ++f) {
            if (pass) ASSERT_TRUE(vial_ir__load(&l, f, &fr));
            else fr = (vial_frame_t){.xyz = (const float*)st->xyz, .num_atoms = st->num_atoms, .cell = st->unitcell};
            for (int b = 0; b < 2; ++b) EXPECT_TRUE(vial_interp_frame(&r[b].in, &fr));
            const vial_ir_t* ir = md_vial__run_ir(r[0].prog);
            for (size_t o = 0; o < md_array_size(ir->outputs); ++o) {
                const vial_val_t* x = vial_interp_value(&r[0].in, ir->outputs[o].value);
                const vial_val_t* y = vial_interp_value(&r[1].in, ir->outputs[o].value);
                ASSERT_TRUE(x && y);
                if (!vial_val_equal(x, y)) printf("%c frame %d: %.*s differs batched\n", structure, (int)f, (int)ir->outputs[o].name.len, ir->outputs[o].name.ptr);
                EXPECT_TRUE(vial_val_equal(x, y));
            }
        }
        // The batched run never ran a call of the map: it made no scratch arena for one
        EXPECT_TRUE(r[0].in.map_scratch[0] == NULL);
        EXPECT_TRUE(r[1].in.map_scratch[0] != NULL);
        if (pass) md_system_state_free(&l.state);
        for (int b = 0; b < 2; ++b) vial_ir__end(&r[b]);
    }
}

// Reductions over every pair (30 October) are lowered to VIAL_IR_PAIRS and never build the matrix; each gives the
// bits the matrix gives, reduced here in C from the matrix the script also publishes. On T, CA against a ligand, and
// on P, CA against CA over three frames.
static float vial_ir__f(const vial_val_t* v, size_t i) { return ((const float*)v->data)[i]; }

UTEST(vial_ir, pairs_reductions) {
    const char* body =
        "dm  = distance(p[:, None], q)\n"
        "a1  = min(dm)\n"
        "a2  = max(dm)\n"
        "a3  = min(dm, axis=1)\n"
        "a4  = max(dm, axis=0)\n"
        "a5  = argmin(dm, axis=1)\n"
        "a6  = argmax(dm, axis=0)\n"
        "a7  = sum(dm < 8.0)\n"
        "a8  = sum(dm <= 8.0, axis=1)\n"
        "a9  = any(dm < 5.0, axis=0)\n"
        "a10 = min(dm, axis=-1)\n"
        "_h  = distance(p[:, None], q)\n"
        "a11 = sum(_h < 6.0, axis=-1)\n"
        "a12 = sum(dm < 1000.0, axis=1)\n";     // past half a period: every line computed, no search
    for (int pass = 0; pass < 2; ++pass) {
        char src[1024];
        snprintf(src, sizeof(src), "%s%s", pass ? "p = particles(name='CA').xyz\nq = p\n" :
                 "p = particles(name='CA').xyz\nq = particles(resname='VLB').xyz\n", body);
        int64_t run_frames = 0;
        const md_system_t* sys = pass ? vial_ir__run_p(&run_frames) : vial_test_system('T');
        ASSERT_TRUE(sys != NULL);
        vial_ir__run_t r;
        ASSERT_TRUE(vial_ir__start(&r, src, sys));
        const vial_ir_t* ir = md_vial__run_ir(r.prog);
        int npairs = 0;
        for (size_t i = 0; i < md_array_size(ir->ops); ++i) npairs += ir->ops[i].kind == VIAL_IR_PAIRS;
        EXPECT_EQ(12, npairs);
        vial_ir__loader_t l = {.sys = sys, .state = {.alloc = r.alloc}};
        if (pass) md_system_state_init(&l.state, sys->atom.count);
        const md_system_state_t* st = vial_test_state('T');
        for (int64_t f = 0; f < (pass ? 3 : 1); ++f) {
            vial_frame_t fr;
            if (pass) ASSERT_TRUE(vial_ir__load(&l, f * 7, &fr));
            else fr = (vial_frame_t){.xyz = (const float*)st->xyz, .num_atoms = st->num_atoms, .cell = st->unitcell};
            EXPECT_TRUE(vial_interp_frame(&r.in, &fr));
            const vial_val_t* dm = vial_interp_value(&r.in, vial_ir__output(ir, "dm"));
            ASSERT_TRUE(dm && dm->ndim == 2);
            // _h is hidden and only reduced: never built
            EXPECT_TRUE(vial_interp_value(&r.in, vial_ir__output(ir, "_h")) == NULL);
            const uint32_t n = dm->dims[0], m = dm->dims[1];
            ASSERT_TRUE(n > 1 && m > 1);
            const vial_val_t* a[13];
            for (int k = 1; k <= 12; ++k) {
                char nm[8];
                snprintf(nm, sizeof(nm), "a%d", k);
                a[k] = vial_interp_value(&r.in, vial_ir__output(ir, nm));
                ASSERT_TRUE(a[k] != NULL);
            }
            float mn = vial_ir__f(dm, 0), mx = mn;
            int32_t c7 = 0;
            for (size_t e = 0; e < (size_t)n * m; ++e) {
                const float d = vial_ir__f(dm, e);
                if (d < mn) mn = d;
                if (d > mx) mx = d;
                c7 += (double)d < 8.0;
            }
            EXPECT_EQ(mn, vial_ir__f(a[1], 0));
            EXPECT_EQ(mx, vial_ir__f(a[2], 0));
            EXPECT_EQ(c7, *(const int32_t*)a[7]->data);
            for (uint32_t i = 0; i < n; ++i) {
                float rmin = vial_ir__f(dm, (size_t)i * m);
                int32_t at = 0, c8 = 0, c11 = 0, c12 = 0;
                for (uint32_t j = 0; j < m; ++j) {
                    const float d = vial_ir__f(dm, (size_t)i * m + j);
                    if (d < rmin) { rmin = d; at = (int32_t)j; }
                    c8 += (double)d <= 8.0;
                    c11 += (double)d < 6.0;
                    c12 += (double)d < 1000.0;
                }
                EXPECT_EQ(c12, ((const int32_t*)a[12]->data)[i]);
                EXPECT_EQ(rmin, vial_ir__f(a[3], i));
                EXPECT_EQ(at, ((const int32_t*)a[5]->data)[i]);
                EXPECT_EQ(c8, ((const int32_t*)a[8]->data)[i]);
                EXPECT_EQ(rmin, vial_ir__f(a[10], i));
                EXPECT_EQ(c11, ((const int32_t*)a[11]->data)[i]);
            }
            for (uint32_t j = 0; j < m; ++j) {
                float cmax = vial_ir__f(dm, j);
                int32_t at = 0;
                bool near = false;
                for (uint32_t i = 0; i < n; ++i) {
                    const float d = vial_ir__f(dm, (size_t)i * m + j);
                    if (d > cmax) { cmax = d; at = (int32_t)i; }
                    near = near || (double)d < 5.0;
                }
                EXPECT_EQ(cmax, vial_ir__f(a[4], j));
                EXPECT_EQ(at, ((const int32_t*)a[6]->data)[j]);
                EXPECT_EQ(near, ((const uint8_t*)a[9]->data)[j] != 0);
            }
        }
        if (pass) md_system_state_free(&l.state);
        vial_ir__end(&r);
    }
}

// within() by every pair and by a spatial search, against the matrix's entries: the same particles, also at a tie,
// where r is the distance of a pair exactly
UTEST(vial_ir, within_matches_matrix) {
    int64_t nf = 0;
    const md_system_t* sys = vial_ir__run_p(&nf);
    ASSERT_TRUE(sys != NULL);
    const char* src =
        "_e = distance(particles.xyz[:, None], particles[0:3].xyz)\n"
        "_f = distance(particles.xyz[:, None], particles(name='CA').xyz)\n"
        "bad = abs(len(within(4.5, of=particles[0:3])) - sum(any(_e <= 4.5, axis=1)))"
        " + abs(len(within(6.0, of=particles(name='CA'))) - sum(any(_f <= 6.0, axis=1)))"
        " + abs(len(within(0.5, of=particles[0:3])) - sum(any(_e <= 0.5, axis=1)))\n";
    for (int plain = 0; plain < 2; ++plain) {
        vial_ir__run_t r;
        ASSERT_TRUE(vial_ir__start(&r, src, sys));
        r.in.plain = plain == 1;
        double bad[64] = {0};
        const int64_t n = nf < 64 ? nf : 64;
        EXPECT_TRUE(vial_ir__frames(&r, sys, n, "bad", bad));
        for (int64_t f = 0; f < n; ++f) EXPECT_EQ(0.0, bad[f]);
        vial_ir__end(&r);
    }
}

// The closest and farthest pairs, which mdlib's kernels find, against the matrix: the same entry, and the first pair in
// row order at it. Coordinates on a grid give many exact ties, and many within an ulp.
UTEST(vial_ir, pair_extremes_match_matrix) {
    const char* body =
        "dm = distance(p[:, None], q)\n"
        "c, ci, cj = closest_pair(p, q)\n"
        "f, fi, fj = farthest_pair(p, q)\n";
    const char* sets[3] = {
        "p = particles(name='CA').xyz\nq = particles.xyz[100:400]\n",
        "p = round(particles(name='CA').xyz / 2.0)\nq = round(particles(name='CA').xyz / 2.0 + 0.25)\n",
        "p = round(particles(name='CA').xyz * 3.0) / 3.0 + 0.1\nq = round(particles.xyz[0:300] * 3.0) / 3.0\n",
    };
    // The trajectory's system is small, so its sets are every particle; its cell is periodic
    const char* traj_sets[3] = {
        "p = particles.xyz\nq = particles.xyz + 0.37\n",
        "p = round(particles.xyz / 2.0)\nq = round(particles.xyz / 2.0 + 0.25)\n",
        "p = round(particles.xyz * 2.0) / 2.0\nq = round(particles.xyz * 2.0) / 2.0 + 0.5\n",
    };
    for (int pass = 0; pass < 6; ++pass) {
        const bool traj = pass >= 3;
        char src[1024];
        snprintf(src, sizeof(src), "%s%s", traj ? traj_sets[pass % 3] : sets[pass % 3], body);
        int64_t run_frames = 0;
        const md_system_t* sys = traj ? vial_ir__run_p(&run_frames) : vial_test_system('T');
        ASSERT_TRUE(sys != NULL);
        vial_ir__run_t r;
        ASSERT_TRUE(vial_ir__start(&r, src, sys));
        const vial_ir_t* ir = md_vial__run_ir(r.prog);
        vial_ir__loader_t l = {.sys = sys, .state = {.alloc = r.alloc}};
        if (traj) md_system_state_init(&l.state, sys->atom.count);
        const md_system_state_t* st = vial_test_state('T');
        for (int64_t f = 0; f < (traj ? 3 : 1); ++f) {
            vial_frame_t fr;
            if (traj) ASSERT_TRUE(vial_ir__load(&l, f * 7, &fr));
            else fr = (vial_frame_t){.xyz = (const float*)st->xyz, .num_atoms = st->num_atoms, .cell = st->unitcell};
            EXPECT_TRUE(vial_interp_frame(&r.in, &fr));
            const vial_val_t* dm = vial_interp_value(&r.in, vial_ir__output(ir, "dm"));
            ASSERT_TRUE(dm && dm->ndim == 2);
            const uint32_t n = dm->dims[0], m = dm->dims[1];
            EXPECT_TRUE((uint64_t)n * m >= 64);     // as many as the kernel is used for
            float mn = vial_ir__f(dm, 0), mx = mn;
            int32_t mi = 0, mj = 0, xi = 0, xj = 0;
            for (uint32_t i = 0; i < n; ++i) {
                for (uint32_t j = 0; j < m; ++j) {
                    const float d = vial_ir__f(dm, (size_t)i * m + j);
                    if (d < mn) { mn = d; mi = (int32_t)i; mj = (int32_t)j; }
                    if (d > mx) { mx = d; xi = (int32_t)i; xj = (int32_t)j; }
                }
            }
            const char* names[6] = {"c", "ci", "cj", "f", "fi", "fj"};
            const vial_val_t* v[6];
            for (int k = 0; k < 6; ++k) {
                v[k] = vial_interp_value(&r.in, vial_ir__output(ir, names[k]));
                ASSERT_TRUE(v[k] != NULL);
            }
            EXPECT_EQ(mn, vial_ir__f(v[0], 0));
            EXPECT_EQ(mi, *(const int32_t*)v[1]->data);
            EXPECT_EQ(mj, *(const int32_t*)v[2]->data);
            EXPECT_EQ(mx, vial_ir__f(v[3], 0));
            EXPECT_EQ(xi, *(const int32_t*)v[4]->data);
            EXPECT_EQ(xj, *(const int32_t*)v[5]->data);
        }
        if (traj) md_system_state_free(&l.state);
        vial_ir__end(&r);
    }
}

// A map whose calls give arrays of sets of different lengths is R0103, not a ragged array (found by fuzzing)
UTEST(vial_ir, ragged_sets_fail) {
    vial_ir__run_t r;
    ASSERT_TRUE(vial_ir__start(&r, "n = [[a for a in r] for r in residues[0:3]]\n", vial_test_system('T')));
    const md_system_state_t* st = vial_test_state('T');
    const vial_frame_t fr = {.xyz = (const float*)st->xyz, .num_atoms = st->num_atoms, .cell = st->unitcell};
    vial_interp_frame(&r.in, &fr);
    const vial_ir_t* ir = md_vial__run_ir(r.prog);
    const vial_val_t* n = vial_interp_value(&r.in, vial_ir__output(ir, "n"));
    // T's first residues have 8, 11 and 9 particles
    EXPECT_TRUE(n == NULL);
    ASSERT_EQ(1, (int)md_array_size(r.in.errors));
    EXPECT_EQ(0, strncmp(r.in.errors[0].code, "R0103", 5));
    vial_ir__end(&r);
}

// eigh: m v = v w and v's columns orthonormal, to float precision, on a symmetric 120 x 120 matrix; eigenvalues
// ascending
UTEST(vial_ir, eigh_decomposes) {
    vial_ir__run_t r;
    const char* src =
        "c = com(residues(chain='A', kind='amino_acid')[0:120])\n"
        "m = unitless(distance(c[:, None], c))\n"
        "w, v = eigh(m)\n"
        "res = max(abs(m @ v - v * w)) / max(abs(w))\n"
        "orth = max(abs(v.T @ v - eye(120)))\n"
        "up = int(all(w[1:] >= w[:-1]))\n";
    ASSERT_TRUE(vial_ir__start(&r, src, vial_test_system('T')));
    const md_system_state_t* st = vial_test_state('T');
    const vial_frame_t fr = {.xyz = (const float*)st->xyz, .num_atoms = st->num_atoms, .cell = st->unitcell};
    const md_tick_t t0 = md_tick_now();
    EXPECT_TRUE(vial_interp_frame(&r.in, &fr));
    const double ms = md_tick_to_nanoseconds(md_tick_now() - t0) * 1e-6;
    const vial_ir_t* ir = md_vial__run_ir(r.prog);
    for (size_t i = 0; i < md_array_size(r.in.errors); ++i) printf("failed: %s %s\n", r.in.errors[i].code, r.in.errors[i].message);
    EXPECT_LT(vial_ir__num(vial_interp_value(&r.in, vial_ir__output(ir, "res"))), 1e-5);
    EXPECT_LT(vial_ir__num(vial_interp_value(&r.in, vial_ir__output(ir, "orth"))), 1e-5);
    EXPECT_EQ(1.0, vial_ir__num(vial_interp_value(&r.in, vial_ir__output(ir, "up"))));
    EXPECT_LT(ms, 2000.0);
    vial_ir__end(&r);
}


// at() inside a def: computed from the def's own values at that frame, also when the def maps over groups, whose sizes
// differ; static where the def's value is
UTEST(vial_ir, at_in_a_def) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    ASSERT_GT(vial_test_run('P'), 3);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(8));
    const char* src =
        "def first_x(r: group):\n    return at(r.xyz, frame=0)[0, 0]\n"
        "def drift(r: group):\n    return rmsd(r, at(r.xyz, frame=0))\n"
        "x0 = first_x(residues)\n"
        "x  = residues[1].xyz[0, 0]\n"
        "per_res = drift(residues)\n"
        "r1 = rmsd(residues[1], at(residues[1].xyz, frame=0))\n";
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    ASSERT_EQ(0u, (unsigned)md_vial_diagnostics(prog, NULL));
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run"), .num_threads = 3}, alloc);
    ASSERT_TRUE(e && md_vial_eval_run(e));
    const md_attributes_t* at = md_vial_eval_attributes(e);
    const md_attribute_t* x0 = md_attributes_find(at, STR_LIT("script/x0"));
    const md_attribute_t* x = md_attributes_find(at, STR_LIT("script/x"));
    const md_attribute_t* pr = md_attributes_find(at, STR_LIT("script/per_res"));
    const md_attribute_t* r1 = md_attributes_find(at, STR_LIT("script/r1"));
    ASSERT_TRUE(x0 && x && pr && r1);
    EXPECT_EQ(1, x0->format.rank);                  // static: one value per residue
    EXPECT_NEAR(((const float*)x->data)[0], ((const float*)x0->data)[1], 1e-5);
    const size_t nf = r1->format.shape[0], nr = pr->format.shape[1];
    ASSERT_EQ(15u, (unsigned)nr);
    for (size_t f = 0; f < nf; ++f) EXPECT_NEAR(((const float*)r1->data)[f], ((const float*)pr->data)[f * nr + 1], 1e-5);
    md_vial_eval_free(e);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

// ### VOLUMES ###

static const md_vial_output_t* vial_ir__pub(const md_vial_program_t* prog, const char* name) {
    const md_vial_output_t* o;
    const size_t n = md_vial_outputs(prog, &o);
    for (size_t i = 0; i < n; ++i) if (str_eq(o[i].name, str_from_cstr(name))) return &o[i];
    return NULL;
}

static double vial_ir__attr1(const md_attributes_t* at, const char* path, size_t i) {
    const md_attribute_t* a = md_attributes_find(at, str_from_cstr(path));
    if (!a || !a->data) return NAN;
    return ((const double*)a->data)[i];
}

// sdf and density (decided 1 November): grids that follow what they cover, a density per Å³ that keeps the count it
// was made from, and what places them published beside them
UTEST(vial_ir, volumes) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    const int64_t nf = vial_test_run('P');
    ASSERT_GT(nf, 3);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(16));
    const char* src =
        "s  = sdf(residues[1:5], particles(name='O'), 6.0)\n"
        "sw = sdf(residues[1:5], particles(name='O'), (6.0, 4.0, 3.0), samples=2, weights=mass(particles(name='O')))\n"
        "d  = density(particles(), box=((-100.0, -100.0, -100.0), (100.0, 100.0, 100.0)), samples=0.25)\n"
        "ds = density(particles(), box=((-100.0, -100.0, -100.0), (100.0, 100.0, 100.0)), samples=0.25, sigma=6.0)\n"
        "dc = density(particles(name='CA'), samples=1)\n"
        "n  = max(s)\n";
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    ASSERT_TRUE(prog != NULL);
    const md_vial_diagnostic_t* dg;
    const size_t ndg = md_vial_diagnostics(prog, &dg);
    for (size_t i = 0; i < ndg; ++i) printf("%.5s %.*s\n", dg[i].code, (int)dg[i].message.len, dg[i].message.ptr);
    ASSERT_EQ(0u, (unsigned)ndg);
    // Types: 12 Å at 4 per Å is 48 a side; 12, 8 and 6 Å at 2 per Å are 24, 16 and 12, rounded up to 24, 16 and 16
    const md_vial_output_t* os = vial_ir__pub(prog, "s");
    const md_vial_output_t* ow = vial_ir__pub(prog, "sw");
    const md_vial_output_t* od = vial_ir__pub(prog, "d");
    const md_vial_output_t* oc = vial_ir__pub(prog, "dc");
    ASSERT_TRUE(os && ow && od && oc);
    EXPECT_STREQ("float[48][48][48] count/Å³", os->type);
    EXPECT_STREQ("float[24][16][16] u/Å³", ow->type);
    EXPECT_STREQ("float[56][56][56] count/Å³", od->type);
    EXPECT_EQ(MD_VIAL_STAGE_TRAJECTORY, os->stage);
    EXPECT_EQ((uint32_t)(MD_VIAL_OUTPUT_VOLUME | MD_VIAL_OUTPUT_SDF), os->flags);
    EXPECT_EQ((uint32_t)MD_VIAL_OUTPUT_VOLUME, od->flags);
    EXPECT_EQ((uint32_t)MD_VIAL_OUTPUT_VOLUME, oc->flags);

    for (int threads = 1; threads <= 4; threads += 3) {
        md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run"), .num_threads = threads}, alloc);
        ASSERT_TRUE(e != NULL);
        ASSERT_TRUE(md_vial_eval_run(e));
        const md_attributes_t* at = md_vial_eval_attributes(e);
        // The box holds every particle, so a density over it, by voxel volume, holds them all: cloud in cell keeps
        // what it spreads, and so does smoothing away from the edges
        const md_attribute_t* ad = md_attributes_find(at, STR_LIT("script/d"));
        const md_attribute_t* as = md_attributes_find(at, STR_LIT("script/ds"));
        ASSERT_TRUE(ad && as);
        ASSERT_EQ(3u, ad->format.rank);
        EXPECT_EQ(56u, ad->format.shape[0]);
        const double step = vial_ir__attr1(at, "script/d/basis", 0);
        EXPECT_NEAR(200.0 / 56.0, step, 1e-9);
        EXPECT_EQ(0.0, vial_ir__attr1(at, "script/d/basis", 1));
        EXPECT_EQ(-100.0, vial_ir__attr1(at, "script/d/origin", 0));
        EXPECT_EQ((double)nf, vial_ir__attr1(at, "script/d/samples", 0));
        double sum = 0, ssum = 0;
        for (size_t i = 0; i < 56 * 56 * 56; ++i) { sum += ((const float*)ad->data)[i]; ssum += ((const float*)as->data)[i]; }
        EXPECT_NEAR((double)sys->atom.count, sum * step * step * step, 1e-2);
        EXPECT_NEAR((double)sys->atom.count, ssum * step * step * step, 1e-2);
        // The bulk density: the particles over the cell
        const double cellv = 46.645 * 96.666 * 48.362;
        EXPECT_NEAR(sys->atom.count / cellv, vial_ir__attr1(at, "script/d/bulk", 0), 1e-6);
        // Over the cell: every CA in it, wrapped
        const md_attribute_t* ac = md_attributes_find(at, STR_LIT("script/dc"));
        ASSERT_TRUE(ac != NULL);
        const size_t nc = (size_t)ac->format.shape[0] * ac->format.shape[1] * ac->format.shape[2];
        EXPECT_EQ(48u, ac->format.shape[0]);
        EXPECT_EQ(104u, ac->format.shape[1]);
        EXPECT_EQ(56u, ac->format.shape[2]);
        double csum = 0;
        for (size_t i = 0; i < nc; ++i) csum += ((const float*)ac->data)[i];
        EXPECT_NEAR(15.0, csum * cellv / (double)nc, 1e-2);
        EXPECT_NEAR(46.645 / 48, vial_ir__attr1(at, "script/dc/basis", 0), 1e-4);
        EXPECT_NEAR(96.666 / 104, vial_ir__attr1(at, "script/dc/basis", 4), 1e-4);
        // An sdf: four structures each frame, in their own frame, centred
        const md_attribute_t* asd = md_attributes_find(at, STR_LIT("script/s"));
        ASSERT_TRUE(asd != NULL);
        EXPECT_EQ(-6.0, vial_ir__attr1(at, "script/s/origin", 2));
        EXPECT_EQ(0.25, vial_ir__attr1(at, "script/s/basis", 8));
        EXPECT_EQ(4.0 * nf, vial_ir__attr1(at, "script/s/samples", 0));
        double smax = 0, stot = 0;
        for (size_t i = 0; i < 48 * 48 * 48; ++i) { smax = MAX(smax, ((const float*)asd->data)[i]); stot += ((const float*)asd->data)[i]; }
        EXPECT_GT(stot, 0.0);
        const md_attribute_t* an = md_attributes_find(at, STR_LIT("script/n"));
        ASSERT_TRUE(an != NULL);
        EXPECT_EQ((float)smax, ((const float*)an->data)[0]);
        // Weighted by mass: an O weighs 16, so the volume is 16 times the count, where the grids are alike
        const md_attribute_t* aw = md_attributes_find(at, STR_LIT("script/sw"));
        ASSERT_TRUE(aw != NULL);
        double wtot = 0;
        for (size_t i = 0; i < 24 * 16 * 16; ++i) wtot += ((const float*)aw->data)[i];
        EXPECT_GT(wtot, 0.0);

        // Over chosen frames: every frame is what was published; its size is known without reading a frame
        uint8_t* all = md_alloc(alloc, (size_t)nf);
        memset(all, 1, (size_t)nf);
        const size_t nv = 48 * 48 * 48;
        EXPECT_EQ(nv, md_vial_eval_aggregate(e, STR_LIT("s"), all, NULL, 0));
        float* over = md_alloc(alloc, nv * sizeof(float));
        ASSERT_EQ(nv, md_vial_eval_aggregate(e, STR_LIT("s"), all, over, nv));
        double maxdiff = 0;
        for (size_t i = 0; i < nv; ++i) maxdiff = MAX(maxdiff, fabs((double)over[i] - ((const float*)asd->data)[i]));
        EXPECT_LT(maxdiff, 1e-5 * (smax + 1));
        // The first frame alone is not the whole
        uint8_t* one = md_alloc(alloc, (size_t)nf);
        memset(one, 0, (size_t)nf);
        one[0] = 1;
        ASSERT_EQ(nv, md_vial_eval_aggregate(e, STR_LIT("s"), one, over, nv));
        double osum = 0;
        for (size_t i = 0; i < nv; ++i) osum += over[i];
        EXPECT_GT(osum, 0.0);
        EXPECT_NE(osum, stot);
        md_vial_eval_free(e);
    }

    // Drawing an sdf: its four structures where the frame has them, each with the pose that places the volume
    const uint32_t expr = md_vial_expr_at(prog, (uint32_t)(strstr(src, "sdf(") - src));
    md_vial_drawing_t dr;
    ASSERT_TRUE(md_vial_draw(&dr, &(md_vial_draw_desc_t){.prog = prog, .sys = sys, .run = STR_LIT("run"), .frame = 2, .element = -1}, expr, alloc));
    ASSERT_EQ(4u, (unsigned)dr.num_poses);
    for (size_t i = 0; i < dr.num_poses; ++i) {
        EXPECT_EQ(10u, dr.poses[i].end - dr.poses[i].beg);
        // A rigid motion: its rotation's columns have length 1
        const float* c0 = dr.poses[i].m[0];
        EXPECT_NEAR(1.0, sqrt(c0[0] * c0[0] + c0[1] * c0[1] + c0[2] * c0[2]), 1e-4);
    }
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

UTEST(vial_ir, volume_diagnostics) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    vial_test_run('P');
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(8));
    struct { const char* src; const char* code; } cases[] = {
        {"v = sdf(residues[0:5], particles(name='O'), 6.0)\n", "E0320"},           // the first residue has three Hs on N
        {"v = sdf(residues[1:5], particles(name='O'), 300.0)\n", "W0407"},         // 2,400 a side at 4 per Å
        {"v = [sdf(r, particles(name='O'), 6.0) for r in residues[1:3]]\n", "E0406"},
        {"v = sdf(residues[1:5], residues, 6.0)\n", "E0313"},                      // an array of sets is no target
        {"v = sdf(residues[1:5], particles(name='O'), 6.0, samples=0)\n", "E0313"},
        {"v = density(particles(), box=((0, 0, 0), (0, 10, 10)))\n", "E0313"},
        {"x = 4.0\nv = sdf(residues[1:5], particles(name='O'), 6.0, ref=particles(name='O').xyz)\n", "E0313"},
    };
    for (size_t i = 0; i < ARRAY_SIZE(cases); ++i) {
        md_vial_program_t* prog = md_vial_compile(str_from_cstr(cases[i].src), &(md_vial_compile_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
        const md_vial_diagnostic_t* d;
        const size_t n = md_vial_diagnostics(prog, &d);
        bool found = false;
        for (size_t k = 0; k < n; ++k) found = found || !strncmp(d[k].code, cases[i].code, 5);
        EXPECT_TRUE(found);
        if (!found) {
            printf("%s", cases[i].src);
            for (size_t k = 0; k < n; ++k) printf("  %.5s %.*s\n", d[k].code, (int)d[k].message.len, d[k].message.ptr);
        }
        md_vial_program_free(prog);
    }
    md_arena_allocator_destroy(alloc);
}

// mcs(a, b): which particles of b play the part of which of a, typed when compiling
UTEST(vial_ir, mcs) {
    const md_system_t* sys = vial_test_system('T');
    ASSERT_TRUE(sys != NULL);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(8));
    const char* src =
        "gtp = particles(resname='GTP')\n"
        "gdp = particles(resname='GDP')\n"
        "a, b = mcs(gtp, gdp)\n"
        "r = rmsd(b.xyz, a.xyz)\n"
        "phe = residues('PHE')\n"
        "tyr = residues('TYR')\n"
        "c = mcs(phe[0], tyr[0], by='element')\n";
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys}, alloc);
    const md_vial_diagnostic_t* d;
    const size_t nd = md_vial_diagnostics(prog, &d);
    for (size_t i = 0; i < nd; ++i) printf("%.5s %.*s\n", d[i].code, (int)d[i].message.len, d[i].message.ptr);
    EXPECT_EQ(0u, (unsigned)nd);
    md_vial__result_t r;
    ASSERT_TRUE(md_vial__line_result(prog, 3, &r));
    EXPECT_STREQ("set[28]", r.types[0]);
    EXPECT_STREQ("set[28]", r.types[1]);
    ASSERT_TRUE(md_vial__line_result(prog, 4, &r));
    EXPECT_STREQ("float Å", r.types[0]);
    // Of the particles, as filters read them: GDP's in the order of GTP's they play
    md_vial_filter_t f;
    ASSERT_TRUE(md_vial_filter(&f, STR_LIT("mcs(particles(resname='GTP'), particles(resname='GDP'))[1]"), &(md_vial_filter_desc_t){.sys = sys}, alloc));
    EXPECT_EQ(28u, (unsigned)f.num_sets);
    // Too little in common: E0321, unless least= asks for less
    struct { const char* src; const char* code; } bad[] = {
        {"m = mcs(particles(resname='GTP'), particles(resname='VLB'))\n", "E0321"},
        {"m = mcs(residues('PHE'), residues('TYR'))\n", "E0313"},
        {"m = mcs(particles(resname='GTP'), particles(resname='GDP'), by='colour')\n", "E0313"},
    };
    for (size_t i = 0; i < ARRAY_SIZE(bad); ++i) {
        md_vial_program_t* p = md_vial_compile(str_from_cstr(bad[i].src), &(md_vial_compile_desc_t){.sys = sys}, alloc);
        const size_t n = md_vial_diagnostics(p, &d);
        bool found = false;
        for (size_t k = 0; k < n; ++k) found = found || !strncmp(d[k].code, bad[i].code, 5);
        EXPECT_TRUE(found);
        md_vial_program_free(p);
    }
    md_vial_program_t* ok = md_vial_compile(STR_LIT("m = mcs(particles(resname='GTP'), particles(resname='VLB'), least=0.2)\n"), &(md_vial_compile_desc_t){.sys = sys}, alloc);
    EXPECT_EQ(0u, (unsigned)md_vial_diagnostics(ok, NULL));
    md_vial_program_free(ok);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

static double vial_ir__max_diff(const md_attributes_t* at, const char* p, const char* q, size_t* n_out) {
    const md_attribute_t* a = md_attributes_find(at, str_from_cstr(p));
    const md_attribute_t* b = md_attributes_find(at, str_from_cstr(q));
    if (!a || !b || !a->data || !b->data) return INFINITY;
    const size_t n = md_attribute_element_count(&a->format);
    if (n != md_attribute_element_count(&b->format)) return INFINITY;
    double d = 0;
    for (size_t i = 0; i < n; ++i) d = MAX(d, fabs((double)((const float*)a->data)[i] - ((const float*)b->data)[i]));
    if (n_out) *n_out = n;
    return d;
}

// An sdf over structures paired column by column, and over kinds of structure placed on the first kind's reference:
// what is alike gives what one kind gives
UTEST(vial_ir, sdf_kinds) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    ASSERT_GT(vial_test_run('P'), 3);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(16));
    const char* src =
        "ala = residues[1:3] and not particles(element='H')\n"
        "cb  = particles(name='CB')\n"
        "s   = sdf(ala, cb, 5.0)\n"
        "a, b = mcs(residues[1], residues[2])\n"
        "so  = sdf([a, b], cb, 5.0)\n"
        "rv  = sdf([a[::-1], b[::-1]], cb, 5.0)\n"
        "s2  = sdf((ala, ala), cb, 5.0)\n"
        "sp  = sdf((ala, ala), cb, 5.0, poses=(None, (eye(3), zeros(3, unit='Å'))))\n"
        "ct  = sdf((residues[1:5], residues[14:15]), cb, 5.0)\n";
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    const md_vial_diagnostic_t* dg;
    const size_t ndg = md_vial_diagnostics(prog, &dg);
    for (size_t i = 0; i < ndg; ++i) printf("%.5s %.*s\n", dg[i].code, (int)dg[i].message.len, dg[i].message.ptr);
    ASSERT_EQ(0u, (unsigned)ndg);
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run"), .num_threads = 3}, alloc);
    ASSERT_TRUE(e != NULL);
    ASSERT_TRUE(md_vial_eval_run(e));
    const md_vial_diagnostic_t* rd;
    const size_t nr = md_vial_eval_diagnostics(e, &rd);
    for (size_t i = 0; i < nr; ++i) printf("%.5s %.*s\n", rd[i].code, (int)rd[i].message.len, rd[i].message.ptr);
    EXPECT_EQ(0u, (unsigned)nr);
    const md_attributes_t* at = md_vial_eval_attributes(e);
    size_t n = 0;
    // Paired column by column as mcs gives them, which is the residues' own order: the same volume; and with the
    // columns the other way round, which pairs the same particles
    EXPECT_LT(vial_ir__max_diff(at, "script/s", "script/so", &n), 1e-6);
    EXPECT_EQ(40u * 40u * 40u, (unsigned)n);
    EXPECT_LT(vial_ir__max_diff(at, "script/s", "script/rv", NULL), 2e-4);    // the fits sum in another order
    // The same kind twice: placed on itself, the same density per structure, twice as much added
    EXPECT_LT(vial_ir__max_diff(at, "script/s", "script/s2", NULL), 1e-4);
    EXPECT_LT(vial_ir__max_diff(at, "script/s", "script/sp", NULL), 1e-4);
    EXPECT_EQ(2 * vial_ir__attr1(at, "script/s/samples", 0), vial_ir__attr1(at, "script/s2/samples", 0));
    // A kind not alike, the last residue with its OXT, placed by the five heavy atoms it has in common
    const md_attribute_t* act = md_attributes_find(at, STR_LIT("script/ct"));
    ASSERT_TRUE(act != NULL);
    double tot = 0;
    for (size_t i = 0; i < 40 * 40 * 40; ++i) tot += ((const float*)act->data)[i];
    EXPECT_GT(tot, 0.0);
    // Drawn: every structure of every kind with its pose
    const uint32_t expr = md_vial_expr_at(prog, (uint32_t)(strstr(src, "sdf((residues") - src));
    md_vial_drawing_t dr;
    ASSERT_TRUE(md_vial_draw(&dr, &(md_vial_draw_desc_t){.prog = prog, .sys = sys, .run = STR_LIT("run"), .frame = 3, .element = -1}, expr, alloc));
    EXPECT_EQ(5u, (unsigned)dr.num_poses);
    md_vial_eval_free(e);
    md_vial_program_free(prog);

    // What the check says
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T != NULL);
    struct { const md_system_t* sys; const char* src; const char* code; } cases[] = {
        {T, "v = sdf((particles(resname='GTP'), particles(resname='VLB')), particles(name='O'), 8.0)\n", "E0321"},
        {T, "v = sdf((residues('PHE')[0:2], residues('TYR')[0:2]), particles(name='O'), 8.0)\n", NULL},
        {T, "v = sdf((residues('PHE'), residues('TYR')), particles(name='O'), 8.0, poses=(None, eye(3)))\n", "E0313"},
        {T, "v = sdf((residues('PHE'), residues('TYR')), particles(name='O'), 8.0, poses=((eye(3), zeros(3)), None))\n", "E0313"},
        {T, "v = sdf((residues('PHE'), residues('TYR')), particles(name='O'), 8.0, ref=(None,))\n", "E0313"},
    };
    for (size_t i = 0; i < ARRAY_SIZE(cases); ++i) {
        md_vial_program_t* p = md_vial_compile(str_from_cstr(cases[i].src), &(md_vial_compile_desc_t){.sys = cases[i].sys}, alloc);
        const size_t nd = md_vial_diagnostics(p, &dg);
        bool found = cases[i].code == NULL ? nd == 0 : false;
        for (size_t k = 0; cases[i].code && k < nd; ++k) found = found || !strncmp(dg[k].code, cases[i].code, 5);
        EXPECT_TRUE(found);
        if (!found) {
            printf("%s", cases[i].src);
            for (size_t k = 0; k < nd; ++k) printf("  %.5s %.*s\n", dg[k].code, (int)dg[k].message.len, dg[k].message.ptr);
        }
        md_vial_program_free(p);
    }
    md_arena_allocator_destroy(alloc);
}
