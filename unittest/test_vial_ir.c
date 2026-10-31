#include "utest.h"

#include <core/md_common.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

#include <md_vial.h>
#include "../src/vial/vial_internal.h"
#include "../src/vial/ir.h"
#include "../src/vial/interp.h"

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
    md_vial_eval_free(eval);

    // Interrupted before it starts: nothing is published
    eval = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    md_vial_eval_interrupt(eval);
    EXPECT_FALSE(md_vial_eval_run(eval));
    EXPECT_TRUE(md_attributes_find(md_vial_eval_attributes(eval), STR_LIT("script/ends")) == NULL);
    md_vial_eval_free(eval);
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
        "ends = distance(ca[0], ca[14])\n"
        "c = clamp(ends)\n"
        "b = bump(ends)\n"
        "t = 30.0 if ends > 30.0 else ends\n";
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
    }
    EXPECT_GT(above, 0);
    EXPECT_LT(above, (int)nf);
    md_system_state_free(&l.state);
    vial_ir__end(&r);
}
