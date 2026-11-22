#include "utest.h"
#include <md_io.h>
#include "test_executor.h"

#include <core/md_common.h>
#include <md_run.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

#include <md_vial.h>
#include "../src/vial/vial_internal.h"
#include "../src/vial/ir.h"
#include "../src/vial/interp.h"
#include <core/md_os.h>
#include <core/md_file.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "vial_test_systems.h"

// Phase 3: the array IR. Lowering, the passes and the interpreter are tested here; the reference's examples are in
// test_vial_spec.c.

static md_vial_program_t* vial_ir__compile(const char* src, char sys, md_allocator_i* alloc) {
    const md_vial_compile_desc_t desc = {.sys = vial_test_system(sys)};
    return md_vial_compile(md_str_from_cstr(src), &desc, alloc);
}

// VIAL_IR_SCRIPT=file [VIAL_IR_SYS=T|P]: prints the IR of a script and what was not lowered. A tool, not a test.
UTEST(vial_ir, dump_script) {
    const char* path = getenv("VIAL_IR_SCRIPT");
    if (!path) return;
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const md_str_t text = md_load_textfile(md_str_from_cstr(path), alloc);
    const char* sys = getenv("VIAL_IR_SYS");
    const md_vial_compile_desc_t desc = {.sys = vial_test_system(sys ? sys[0] : 0)};
    md_vial_program_t* prog = md_vial_compile(text, &desc, alloc);
    const vial_ir_t* ir = md_vial__ir(prog);
    const md_str_t d = vial_ir_dump(ir, text, alloc);
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
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    md_vial_program_t* prog = vial_ir__compile("a = 1 + 2\nb = residue(chain='A')\nc = distance(b.com, b[0].com)\n", 'T', alloc);
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
#include <md_tpr.h>
#include <md_trr.h>
#include <md_hbond.h>
#include <md_mmcif.h>
#include <md_attributes.h>
#include <md_infer.h>
#include <core/md_geom.h>
#include <core/md_vec_math.h>
#include <core/md_linalg.h>

typedef struct vial_ir__run_t {
    md_allocator_i*     alloc;
    md_vial_program_t*  prog;
    vial_sys_t          vs;
    vial_interp_t       in;
    bool                ok;
} vial_ir__run_t;

static bool vial_ir__start(vial_ir__run_t* r, const char* src, const md_system_t* sys) {
    memset(r, 0, sizeof(*r));
    r->alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const md_vial_compile_desc_t desc = {.sys = sys};
    r->prog = md_vial_compile(md_str_from_cstr(src), &desc, r->alloc);
    const vial_ir_t* ir = md_vial__run_ir(r->prog);
    if (!ir) {
        const md_vial_diagnostic_t* dg;
        const size_t ndg = md_vial_diagnostics(r->prog, &dg);
        for (size_t i = 0; i < ndg; ++i) printf("%.5s %.*s\n", dg[i].code, (int)dg[i].message.len, dg[i].message.ptr);
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
        if (!md_str_eq(ir->outputs[i].name, md_str_from_cstr(name))) continue;
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
        "n = len(residue(chain='A'))\n"
        "o = count(residue(chain='A'))[0:3]\n"
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
        md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
        sys = (md_system_t){.alloc = alloc};
        md_system_state_t state = {.alloc = alloc};
        const md_str_t path = MD_STR_LIT(MD_UNITTEST_DATA_DIR "/1ALA-560ns.pdb");
        if (md_pdb_system_init_from_file(&sys, &state, path, NULL)) {
            md_system_infer(&sys, &state, NULL);
            if (md_pdb_system_publish_run(&sys, path, MD_STR_LIT("run"), MD_RUN_FLAG_DISABLE_CACHE_WRITE)) {
                const md_attribute_t* t = md_attributes_find(&sys.attributes, MD_STR_LIT("run/time"));
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
    const md_str_t paths[] = {MD_STR_LIT("atom/position"), MD_STR_LIT("unitcell")};
    md_system_extract_t* ex = md_system_extract_begin(sys, MD_STR_LIT("run"), paths, 2, r->alloc);
    if (!ex) return false;
    const vial_ir_t* ir = md_vial__run_ir(r->prog);
    vial_ssa_t v = 0;
    for (size_t i = 0; i < md_array_size(ir->outputs); ++i) if (md_str_eq(ir->outputs[i].name, md_str_from_cstr(name))) v = ir->outputs[i].value;
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
        "ca = particle(name='CA')\n"
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
    const md_str_t paths[] = {MD_STR_LIT("atom/position"), MD_STR_LIT("unitcell")};
    md_system_extract_t* ex = md_system_extract_begin(sys, MD_STR_LIT("run"), paths, 2, r.alloc);
    uint32_t ca[4], k = 0;
    for (uint32_t i = 0; i < sys->atom.count && k < 4; ++i) if (md_str_eq(md_atom_name(&sys->atom, i), MD_STR_LIT("CA"))) ca[k++] = i;
    ASSERT_EQ(4u, k);
    for (int64_t f = 0; f < nf; ++f) {
        md_system_extract_frame(ex, f, &state);
        md_float3_t dx = md_float3_sub(state.xyz[ca[3]], state.xyz[ca[0]]);
        md_geom_min_image(dx.elem, 0, 1, &state.unitcell);
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
    const md_str_t paths[] = {MD_STR_LIT("atom/position"), MD_STR_LIT("unitcell")};
    md_system_extract_t* ex = md_system_extract_begin(l->sys, MD_STR_LIT("run"), paths, 2, md_get_heap_allocator());
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
    for (size_t i = 0; i < md_array_size(ir->outputs); ++i) if (md_str_eq(ir->outputs[i].name, md_str_from_cstr(name))) return ir->outputs[i].value;
    return 0;
}

// __all__: exactly the names it lists are published; the others stay for hover, as _names do
UTEST(vial_ir, all_publishes) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    const char* src = "ca = particle(name='CA')\nd = distance(ca[0], ca[14])\n_h = d * 2\nshown_too = _h + 1\n__all__ = ['d', '_h']\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys}, alloc);
    ASSERT_EQ(0u, (unsigned)md_vial_diagnostics(prog, NULL));
    const md_vial_output_t* o;
    const size_t n = md_vial_outputs(prog, &o);
    bool d = false, h = false, s = false, ca = false, all = false;
    for (size_t i = 0; i < n; ++i) {
        d = d || md_str_eq(o[i].name, MD_STR_LIT("d"));
        h = h || md_str_eq(o[i].name, MD_STR_LIT("_h"));
        s = s || md_str_eq(o[i].name, MD_STR_LIT("shown_too"));
        ca = ca || md_str_eq(o[i].name, MD_STR_LIT("ca"));
        all = all || md_str_eq(o[i].name, MD_STR_LIT("__all__"));
    }
    EXPECT_TRUE(d);
    EXPECT_TRUE(h);         // listed, so published despite its _
    EXPECT_FALSE(s);
    EXPECT_FALSE(ca);
    EXPECT_FALSE(all);
    const vial_ir_t* ir = md_vial__run_ir(prog);
    ASSERT_TRUE(ir != NULL);
    for (size_t i = 0; i < md_array_size(ir->outputs); ++i) {
        const bool listed = md_str_eq(ir->outputs[i].name, MD_STR_LIT("d")) || md_str_eq(ir->outputs[i].name, MD_STR_LIT("_h"));
        EXPECT_EQ(!listed, ir->outputs[i].hidden);
    }
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

// pi, prod, median, percentile, cumsum, sort and argsort: NumPy's answers on small arrays
UTEST(vial_ir, more_builtins) {
    vial_ir__run_t r;
    const char* src =
        "a   = [3.0, 1.0, 2.0]\n"
        "s   = sort(a)\n"
        "i   = argsort(a)\n"
        "c   = cumsum([1, 2, 3])\n"
        "cb  = cumsum([True, False, True])\n"
        "m   = median(a)\n"
        "m2  = median([4.0, 1.0, 3.0, 2.0])\n"
        "p   = percentile(a, 50)\n"
        "p0  = percentile([1.0, 2.0, 3.0, 4.0], 25)\n"
        "p1  = percentile([1.0, 2.0, 3.0, 4.0], 100)\n"
        "pr  = prod([2, 3, 4])\n"
        "tp  = 2 * pi\n"
        "sm  = sort([[3, 1], [2, 0]])\n"
        "sa  = sort([[3, 1], [2, 0]], axis=0)\n"
        "ia  = argsort([[3, 1], [2, 0]], axis=0)\n"
        "c2  = cumsum([[1, 2], [3, 4]], axis=1)\n"
        "mr  = median([[1.0, 5.0], [3.0, 2.0]], axis=1)\n";
    ASSERT_TRUE(vial_ir__start(&r, src, NULL));
    const bool ok = vial_interp_static(&r.in);
    if (!ok) printf("failed: %s %s\n", r.in.error.code, r.in.error.message);
    EXPECT_TRUE(ok);
    char buf[256];
    EXPECT_STREQ("[1.0, 2.0, 3.0]", vial_ir__out(&r, "s", buf, sizeof(buf)));
    EXPECT_STREQ("[1, 2, 0]", vial_ir__out(&r, "i", buf, sizeof(buf)));
    EXPECT_STREQ("[1, 3, 6]", vial_ir__out(&r, "c", buf, sizeof(buf)));
    EXPECT_STREQ("[1, 1, 2]", vial_ir__out(&r, "cb", buf, sizeof(buf)));
    EXPECT_STREQ("2.0", vial_ir__out(&r, "m", buf, sizeof(buf)));
    EXPECT_STREQ("2.5", vial_ir__out(&r, "m2", buf, sizeof(buf)));
    EXPECT_STREQ("2.0", vial_ir__out(&r, "p", buf, sizeof(buf)));
    EXPECT_STREQ("1.75", vial_ir__out(&r, "p0", buf, sizeof(buf)));
    EXPECT_STREQ("4.0", vial_ir__out(&r, "p1", buf, sizeof(buf)));
    EXPECT_STREQ("24", vial_ir__out(&r, "pr", buf, sizeof(buf)));
    EXPECT_STREQ("6.2831855", vial_ir__out(&r, "tp", buf, sizeof(buf)));
    EXPECT_STREQ("[[1, 3], [0, 2]]", vial_ir__out(&r, "sm", buf, sizeof(buf)));
    EXPECT_STREQ("[[2, 0], [3, 1]]", vial_ir__out(&r, "sa", buf, sizeof(buf)));
    EXPECT_STREQ("[[1, 1], [0, 0]]", vial_ir__out(&r, "ia", buf, sizeof(buf)));
    EXPECT_STREQ("[[1, 3], [3, 7]]", vial_ir__out(&r, "c2", buf, sizeof(buf)));
    EXPECT_STREQ("[3.0, 2.5]", vial_ir__out(&r, "mr", buf, sizeof(buf)));
    vial_ir__end(&r);
}

// weights= means the same everywhere: 'mass', None, or one number per point, so the array forms give what the words give
UTEST(vial_ir, weights_uniform) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    ASSERT_GT(vial_test_run('P'), 2);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(8));
    const char* src =
        "s   = residue[3]\n"
        "ref = at(s.xyz, frame=0)\n"
        "c1  = com(s)[0]\n"
        "c2  = com(s, weights=mass(s))[0]\n"
        "c3  = com(s, weights=None)[0]\n"
        "c4  = com(s, weights=ones(count(s)))[0]\n"
        "g1  = gyration(s)\n"
        "g2  = gyration(s, weights=mass(s))\n"
        "g3  = gyration(s, weights=None)\n"
        "g4  = gyration(s, weights=ones(count(s)))\n"
        "r1  = rmsd(s, ref)\n"
        "r2  = rmsd(s, ref, weights=ones(count(s)))\n"
        "r3  = rmsd(s, ref, weights='mass')\n"
        "r4  = rmsd(s, ref, weights=mass(s))\n"
        "f1  = fit(s.xyz, ref, weights='mass')[0, 0]\n"
        "f2  = fit(s.xyz, ref, weights=mass(s))[0, 0]\n"
        "w1  = weight(s)[0]\n"
        "w2  = mass(s)[0]\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    const md_vial_diagnostic_t* dg;
    const size_t ndg = md_vial_diagnostics(prog, &dg);
    for (size_t i = 0; i < ndg; ++i) printf("%.5s %.*s\n", dg[i].code, (int)dg[i].message.len, dg[i].message.ptr);
    ASSERT_EQ(0u, (unsigned)ndg);
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run"), .executor = test_executor_of_threads(2)}, alloc);
    ASSERT_TRUE(e != NULL);
    const bool ran = md_vial_eval_run(e);
    const md_vial_diagnostic_t* ed;
    const size_t ned = md_vial_eval_diagnostics(e, &ed);
    for (size_t i = 0; i < ned; ++i) printf("%.5s %.*s\n", ed[i].code, (int)ed[i].message.len, ed[i].message.ptr);
    ASSERT_TRUE(ran);
    const md_attributes_t* at = md_vial_eval_attributes(e);
    static const char* pairs[][2] = {{"c1", "c2"}, {"c3", "c4"}, {"g1", "g2"}, {"g3", "g4"}, {"r1", "r2"}, {"r3", "r4"}, {"f1", "f2"}, {"w1", "w2"}};
    for (size_t q = 0; q < MD_COUNTOF(pairs); ++q) {
        char pa[32], pb[32];
        snprintf(pa, sizeof(pa), "script/%s", pairs[q][0]);
        snprintf(pb, sizeof(pb), "script/%s", pairs[q][1]);
        const md_attribute_t* a = md_attributes_find(at, md_str_from_cstr(pa));
        const md_attribute_t* b = md_attributes_find(at, md_str_from_cstr(pb));
        ASSERT_TRUE(a && b && a->data && b->data);
        const size_t n = a->format.rank ? a->format.shape[0] : 1;
        ASSERT_EQ(n, b->format.rank ? b->format.shape[0] : 1);
        for (size_t f = 0; f < n; ++f) {
            const float va = ((const float*)a->data)[f], vb = ((const float*)b->data)[f];
            if (fabsf(va - vb) > 1e-4f) printf("%s and %s differ at %zu: %g and %g\n", pairs[q][0], pairs[q][1], f, va, vb);
            EXPECT_NEAR(va, vb, 1e-4);
        }
    }
    // Different weights give a different number
    const md_attribute_t* g1 = md_attributes_find(at, MD_STR_LIT("script/g1"));
    const md_attribute_t* g3 = md_attributes_find(at, MD_STR_LIT("script/g3"));
    ASSERT_TRUE(g1 && g3);
    EXPECT_GT(fabsf(((const float*)g1->data)[0] - ((const float*)g3->data)[0]), 1e-4f);
    md_vial_eval_free(e);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

// A for is unrolled: its body runs once per element, what it binds carries from turn to turn and out of the loop,
// and break and continue on static conditions cut it short
UTEST(vial_ir, loops_unrolled) {
    vial_ir__run_t r;
    const char* src =
        "def f(n):\n    t = 0\n    for i in range(n):\n        t = t + i\n    return t\n"
        "x = f(4)\n"
        "t = 0\nfor q in range(4):\n    t = t + q\n"
        "u = t\n"
        "s = 0\nfor q in range(10):\n    if q == 3:\n        break\n    s += q\n"
        "c = 0\nfor q in range(6):\n    if q % 2 == 0:\n        continue\n    c += q\n"
        "n = 0\nfor res in residue[0:3]:\n    n += count(res)\n"
        "m = 0.0\nfor p in residue[0]:\n    m += mass(p)[0]\n"
        "z = 0\nfor a, b in zip([1, 2, 3], [10, 20, 30]):\n    z += a * b\n"
        "w = 0\nfor v in [5, 6, 7]:\n    w = w * 10 + v\n"
        "e = 0\nfor q in range(0):\n    e = 99\n"
        "n2 = 0\nfor p in residue[0]:\n    n2 += 1\n"
        "n3 = 0\nfor p in residue[0]:\n    n3 = n3 + count(p)\n"
        "mr = sum(mass(residue[0]))\n"
        "g = 0\nfor i in range(3):\n    for j in range(3):\n        if j > i:\n            break\n        g += 1\n"
        "def h(k):\n    acc = 0\n    for i in range(k):\n        acc += i * i\n    return acc\n"
        "hh = h(4)\n";
    ASSERT_TRUE(vial_ir__start(&r, src, vial_test_system('T')));
    const bool ok = vial_interp_static(&r.in);
    if (!ok) printf("failed: %s %s\n", r.in.error.code, r.in.error.message);
    EXPECT_TRUE(ok);
    char buf[256];
    EXPECT_STREQ("6", vial_ir__out(&r, "x", buf, sizeof(buf)));
    EXPECT_STREQ("6", vial_ir__out(&r, "u", buf, sizeof(buf)));
    EXPECT_STREQ("6", vial_ir__out(&r, "t", buf, sizeof(buf)));
    EXPECT_STREQ("3", vial_ir__out(&r, "s", buf, sizeof(buf)));
    EXPECT_STREQ("9", vial_ir__out(&r, "c", buf, sizeof(buf)));
    EXPECT_STREQ("28", vial_ir__out(&r, "n", buf, sizeof(buf)));       // MET, ARG, GLU of chain A: 8 + 11 + 9 heavy atoms
    EXPECT_STREQ("8", vial_ir__out(&r, "n2", buf, sizeof(buf)));
    EXPECT_STREQ("8", vial_ir__out(&r, "n3", buf, sizeof(buf)));
    EXPECT_STREQ("140", vial_ir__out(&r, "z", buf, sizeof(buf)));
    EXPECT_STREQ("567", vial_ir__out(&r, "w", buf, sizeof(buf)));
    EXPECT_STREQ("0", vial_ir__out(&r, "e", buf, sizeof(buf)));
    char ref[256];
    EXPECT_STREQ(vial_ir__out(&r, "mr", ref, sizeof(ref)), vial_ir__out(&r, "m", buf, sizeof(buf)));
    EXPECT_STREQ("6", vial_ir__out(&r, "g", buf, sizeof(buf)));         // (i, j) with j <= i: 1 + 2 + 3
    EXPECT_STREQ("14", vial_ir__out(&r, "hh", buf, sizeof(buf)));
    vial_ir__end(&r);
}

UTEST(vial_ir, evaluates_trajectory) {
    const md_system_t* sys = vial_test_system('P');
    const int64_t nf = vial_test_run('P');
    ASSERT_TRUE(sys && nf > 2);
    vial_ir__run_t r;
    const char* src =
        "ca   = particle(name='CA')\n"
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
        const md_str_t d = vial_ir_dump(ir, md_str_from_cstr(src), r.alloc);
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
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const char* src =
        "ca   = particle(name='CA')\n"
        "ends = distance(ca[0], ca[14])\n"
        "m    = mean(ends, axis='frame')\n"
        "dd   = distance(ca.xyz, ca[0])\n"
        "_hid = ends * 2.0\n"
        "big  = zeros(100) + ends\n";
    const md_vial_compile_desc_t cdesc = {.sys = sys, .run = MD_STR_LIT("run")};
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &cdesc, alloc);
    const md_vial_eval_desc_t edesc = {.sys = sys, .run = MD_STR_LIT("run"), .max_elements = 50};
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
    const md_attribute_t* ends = md_attributes_find(attrs, MD_STR_LIT("script/ends"));
    const md_attribute_t* m = md_attributes_find(attrs, MD_STR_LIT("script/m"));
    const md_attribute_t* dd = md_attributes_find(attrs, MD_STR_LIT("script/dd"));
    ASSERT_TRUE(ends && m && dd);
    EXPECT_TRUE(ends->flags & MD_ATTRIBUTE_FLAG_TEMPORAL);
    EXPECT_FALSE(m->flags & MD_ATTRIBUTE_FLAG_TEMPORAL);
    EXPECT_EQ(2u, dd->format.rank);
    EXPECT_EQ((uint32_t)nf, dd->format.shape[0]);
    EXPECT_EQ(15u, dd->format.shape[1]);
    EXPECT_TRUE(md_attributes_find(attrs, MD_STR_LIT("script/_hid")) == NULL);
    EXPECT_TRUE(md_attributes_find(attrs, MD_STR_LIT("script/ca")) == NULL);      // a selection is a recipe
    double sum = 0.0, mean = 0.0;
    for (int64_t f = 0; f < nf; ++f) {
        double v;
        md_attribute_extract_f64(&v, 1, attrs, ends, md_attribute_slice_1((uint32_t)f), ends->unit);
        sum += v;
    }
    md_attribute_extract_f64(&mean, 1, attrs, m, md_attribute_slice_all(), m->unit);
    EXPECT_NEAR(sum / nf, mean, 1e-4);

    // Hover: an expression's value at a frame, the bits the evaluation published
    {
        char hv[160], want[64];
        const uint32_t call = (uint32_t)(strstr(src, "distance(ca[0]") - src) + 8;     // its '('
        for (int64_t f = 0; f < nf; f += nf / 3) {
            float v;
            md_attribute_extract_f32(&v, 1, attrs, ends, md_attribute_slice_1((uint32_t)f), ends->unit);
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
            md_attribute_extract_f64(&mv, 1, attrs, m, md_attribute_slice_all(), m->unit);
            char wm[64];
            vial_print_float(wm, sizeof(wm), (float)mv, true);
            EXPECT_TRUE(md_vial_eval_hover(eval, (uint32_t)(strstr(src, "mean(") - src) + 4, 0, hv, sizeof(hv)));
            EXPECT_STREQ(wm, hv);
        }
        EXPECT_FALSE(md_vial_eval_hover(eval, 0xFFFFFFu, 0, hv, sizeof(hv)));
    }
    md_vial_eval_free(eval);

    // Interrupted before it starts: a value per frame is there, as written so far, which is nothing; the rest is not
    eval = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    md_vial_eval_interrupt(eval);
    EXPECT_FALSE(md_vial_eval_run(eval));
    {
        const md_attribute_t* e2 = md_attributes_find(md_vial_eval_attributes(eval), MD_STR_LIT("script/ends"));
        ASSERT_TRUE(e2 != NULL);
        for (int64_t f = 0; f < nf; ++f) EXPECT_TRUE(isnan(((const float*)e2->data)[f]));
        EXPECT_TRUE(md_attributes_find(md_vial_eval_attributes(eval), MD_STR_LIT("script/m")) == NULL);
        uint8_t held[256];
        ASSERT_EQ((size_t)nf, md_vial_eval_output_frames(eval, MD_STR_LIT("ends"), held, sizeof(held)));
        for (int64_t f = 0; f < nf; ++f) EXPECT_EQ(0, held[f]);
    }
    md_vial_eval_free(eval);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

static uint64_t vial_ir__fp(const md_vial_output_t* o, size_t n, const char* name) {
    for (size_t i = 0; i < n; ++i) if (md_str_eq(o[i].name, md_str_from_cstr(name))) return o[i].fingerprint;
    return 0;
}

// The evaluation API on several threads, its own and a caller's: the same attributes, byte for byte, and the same
// diagnostics in the same order as on one
UTEST(vial_ir, eval_threads) {
    const md_system_t* sys = vial_test_system('P');
    const int64_t nf = vial_test_run('P');
    ASSERT_TRUE(sys && nf > 2);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const char* src =
        "ca   = particle(name='CA')\n"
        "ends = distance(ca[0], ca[14])\n"
        "m    = mean(ends, axis='frame')\n"
        "dev  = ends - m\n"
        "near = len(within(6.0, of=ca[0:3]))\n"
        "cp   = distance(ca[0:5], ca[8:15], between='closest')\n"
        "big  = zeros(100) + ends\n"
        "acc  = cumulative(ends, op='sum')\n";
    const md_vial_compile_desc_t cdesc = {.sys = sys, .run = MD_STR_LIT("run")};
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &cdesc, alloc);
    ASSERT_TRUE(prog != NULL);
    md_vial_eval_t* ev[3];
    test_executor_t shuffled = test_executor_shuffled(3);   // out of order, thread_idx in turn
    ev[0] = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run"), .max_elements = 50}, alloc);
    ev[1] = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run"), .max_elements = 50, .executor = test_executor_of_threads(4)}, alloc);
    ev[2] = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run"), .max_elements = 50,
                                                             .executor = test_executor(&shuffled)}, alloc);
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
        for (md_attribute_iter_t it = md_attributes_iter(a0, MD_STR_LIT("script")); md_attributes_next(&it);) {
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
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const char* a =
        "ca = particle(name='CA')\n"
        "d  = distance(ca[0], ca[1])\n"
        "e  = d * 2.0\n"
        "f  = distance(ca[0], ca[2])\n"
        "h  = histogram(d, axis='frame', bins=10, range=(0.0, 20.0))\n"
        "m  = mean(d, axis='frame')\n";
    const char* b =
        "# the same, renamed and moved\n"
        "ca = particle(name='CA')\n\n"
        "dist = distance(ca[0],   ca[1])\n"
        "m  = mean(dist, axis='frame')\n"
        "e  = dist * 3.0\n"
        "f  = distance(ca[0], ca[3])\n"
        "h  = histogram(dist, axis='frame', bins=10, range=(0.0, 20.0))\n";
    const md_vial_compile_desc_t desc = {.sys = sys, .run = MD_STR_LIT("run")};
    md_vial_program_t* pa = md_vial_compile(md_str_from_cstr(a), &desc, alloc);
    md_vial_program_t* pb = md_vial_compile(md_str_from_cstr(b), &desc, alloc);
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
        if (md_str_eq(oa[i].name, MD_STR_LIT("d"))) EXPECT_TRUE(oa[i].flags & MD_VIAL_OUTPUT_TEMPORAL);
        if (md_str_eq(oa[i].name, MD_STR_LIT("h"))) { EXPECT_TRUE(oa[i].flags & MD_VIAL_OUTPUT_DISTRIBUTION); EXPECT_FALSE(oa[i].flags & MD_VIAL_OUTPUT_TEMPORAL); }
        if (md_str_eq(oa[i].name, MD_STR_LIT("m"))) EXPECT_EQ(0u, oa[i].flags);
    }
    #undef FP

    // Taken over: an evaluation of b made from one of a publishes what a computed alike as a published it, and
    // computes the rest
    {
        md_vial_eval_t* ea = md_vial_eval_create(pa, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
        ASSERT_TRUE(ea != NULL);
        md_vial_plan_t plan[16];
        // Before a has run every frame, nothing is held
        ASSERT_TRUE(md_vial_eval_plan(plan, 16, ea, pb, (size_t)nf));
        for (size_t i = 0; i < nb; ++i) EXPECT_NE((int)MD_VIAL_PLAN_HELD, (int)plan[i].status);
        EXPECT_TRUE(md_vial_eval_run(ea));
        ASSERT_TRUE(md_vial_eval_plan(plan, 16, ea, pb, (size_t)nf));
        for (size_t i = 0; i < nb; ++i) {
            const char* want = md_str_eq(ob[i].name, MD_STR_LIT("e")) || md_str_eq(ob[i].name, MD_STR_LIT("f")) ? "changed" : ob[i].selection ? "new" : "held";
            const char* got = plan[i].status == MD_VIAL_PLAN_HELD ? "held" : plan[i].status == MD_VIAL_PLAN_CHANGED ? "changed" : "new";
            if (md_str_eq(ob[i].name, MD_STR_LIT("ca"))) continue;   // a selection: never taken over, evaluated as it is
            EXPECT_STREQ(want, got);
            if (md_str_eq(ob[i].name, MD_STR_LIT("dist"))) EXPECT_TRUE(md_str_eq(plan[i].from, MD_STR_LIT("d")));
        }
        md_vial_eval_t* eb = md_vial_eval_create_from(ea, pb, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
        ASSERT_TRUE(eb != NULL);
        // Taken over: published from when it is made, before any run, with the histogram's range
        EXPECT_TRUE(md_attributes_find(md_vial_eval_attributes(eb), MD_STR_LIT("script/dist")) != NULL);
        EXPECT_TRUE(md_attributes_find(md_vial_eval_attributes(eb), MD_STR_LIT("script/h/range")) != NULL);
        EXPECT_TRUE(md_vial_eval_run(eb));
        const md_attributes_t* ta = md_vial_eval_attributes(ea);
        const md_attributes_t* tb = md_vial_eval_attributes(eb);
        const char* pairs[][2] = {{"script/d", "script/dist"}, {"script/m", "script/m"}, {"script/h", "script/h"}};
        for (int k = 0; k < 3; ++k) {
            const md_attribute_t* x = md_attributes_find(ta, md_str_from_cstr(pairs[k][0]));
            const md_attribute_t* y = md_attributes_find(tb, md_str_from_cstr(pairs[k][1]));
            ASSERT_TRUE(x && y);
            const size_t bytes = md_attribute_byte_size(&x->format);
            ASSERT_EQ(bytes, md_attribute_byte_size(&y->format));
            EXPECT_EQ(0, memcmp(x->data, y->data, bytes));
            EXPECT_EQ(x->flags, y->flags);
        }
        // Changed: computed again, e = dist * 3 where a had d * 2
        const md_attribute_t* ae = md_attributes_find(ta, MD_STR_LIT("script/e"));
        const md_attribute_t* be = md_attributes_find(tb, MD_STR_LIT("script/e"));
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
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const char* src =
        "ca = particle(name='CA')\n"
        "d  = distance(ca[0], ca[1])\n"
        "h  = histogram(d, axis='frame', bins=10, range=(0.0, 20.0))\n"
        "m  = mean(d, axis='frame')\n"
        "x  = max(d * 2.0, axis='frame')\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(prog != NULL);
    for (int threads = 1; threads <= 3; threads += 2) {
        md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run"), .executor = test_executor_of_threads((uint32_t)(threads))}, alloc);
        ASSERT_TRUE(e != NULL);
        uint8_t* all = md_alloc(alloc, (size_t)nf);
        uint8_t* odd = md_alloc(alloc, (size_t)nf);
        for (int64_t f = 0; f < nf; ++f) { all[f] = 1; odd[f] = (uint8_t)(f & 1); }
        float h[10];
        EXPECT_EQ(0u, md_vial_eval_aggregate(e, MD_STR_LIT("h"), all, h, 10));   // not run yet
        // Stopped before it starts, it stays stopped until reset
        md_vial_eval_interrupt(e);
        EXPECT_FALSE(md_vial_eval_run(e));
        EXPECT_FALSE(md_vial_eval_run(e));
        md_vial_eval_reset_interrupt(e);
        EXPECT_TRUE(md_vial_eval_run(e));
        const md_attributes_t* attr = md_vial_eval_attributes(e);
        const md_attribute_t* ad = md_attributes_find(attr, MD_STR_LIT("script/d"));
        const md_attribute_t* ah = md_attributes_find(attr, MD_STR_LIT("script/h"));
        const md_attribute_t* am = md_attributes_find(attr, MD_STR_LIT("script/m"));
        ASSERT_TRUE(ad && ah && am);
        const float* d = (const float*)ad->data;
        // Beside the histogram, its range: the bin edges, in the unit of what it counts
        const md_attribute_t* ar = md_attributes_find(attr, MD_STR_LIT("script/h/range"));
        ASSERT_TRUE(ar != NULL);
        double rr[2];
        ASSERT_EQ(2u, md_attribute_extract_f64(rr, 2, attr, ar, md_attribute_slice_all(), ar->unit));
        EXPECT_EQ(0.0, rr[0]);
        EXPECT_EQ(20.0, rr[1]);
        EXPECT_TRUE(md_unit_base_equal(ar->unit, ad->unit));

        // Over every frame: what was published
        ASSERT_EQ(10u, md_vial_eval_aggregate(e, MD_STR_LIT("h"), all, h, 10));
        for (int i = 0; i < 10; ++i) EXPECT_EQ(((const float*)ah->data)[i], h[i]);
        float m = 0;
        ASSERT_EQ(1u, md_vial_eval_aggregate(e, MD_STR_LIT("m"), all, &m, 1));
        EXPECT_EQ(((const float*)am->data)[0], m);

        // Over the odd frames: counted by hand from the per-frame distances
        float want[10] = {0};
        double sum = 0, mx = -INFINITY;
        int n = 0;
        for (int64_t f = 1; f < nf; f += 2) {
            const int bin = (int)(d[f] / 2.0f);
            if (bin >= 0 && bin < 10) want[bin] += 1.0f;
            sum += d[f]; mx = MD_MAX(mx, 2.0 * d[f]); ++n;
        }
        ASSERT_EQ(10u, md_vial_eval_aggregate(e, MD_STR_LIT("h"), odd, h, 10));
        for (int i = 0; i < 10; ++i) EXPECT_EQ(want[i], h[i]);
        ASSERT_EQ(1u, md_vial_eval_aggregate(e, MD_STR_LIT("m"), odd, &m, 1));
        EXPECT_NEAR((float)(sum / n), m, 1e-4f);
        float x = 0;
        ASSERT_EQ(1u, md_vial_eval_aggregate(e, MD_STR_LIT("x"), odd, &x, 1));
        EXPECT_NEAR((float)mx, x, 1e-4f);

        // Too small an out: the count, nothing written; not over the frames, or unknown: 0
        h[0] = -1.0f;
        EXPECT_EQ(10u, md_vial_eval_aggregate(e, MD_STR_LIT("h"), odd, h, 4));
        EXPECT_EQ(-1.0f, h[0]);
        EXPECT_EQ(0u, md_vial_eval_aggregate(e, MD_STR_LIT("d"), odd, h, 10));
        EXPECT_EQ(0u, md_vial_eval_aggregate(e, MD_STR_LIT("nope"), odd, h, 10));
        md_vial_eval_free(e);
    }
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

static uint32_t vial_test_count_ca(const md_system_t* sys) {
    uint32_t n = 0;
    for (size_t i = 0; i < sys->atom.count; ++i) n += md_str_eq(md_atom_name(&sys->atom, i), MD_STR_LIT("CA"));
    return n;
}

// Every form of angle draws its arms and arc: three points, the points stacked in one argument, over arrays, and in a
// comprehension, where hovering the body draws every element
static uint32_t vial_ir__prims(const char* src, const char* at_text, md_allocator_i* alloc) {
    const md_system_t* sys = vial_test_system('P');
    const md_system_state_t* st = vial_test_state('P');
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys}, alloc);
    if (md_vial_diagnostics(prog, NULL)) { md_vial_program_free(prog); return UINT32_MAX; }
    md_vial_draw_desc_t desc = {.prog = prog, .sys = sys, .frame = 0, .state = st, .element = -1};
    const uint32_t at = md_vial_expr_at(prog, (uint32_t)(strstr(src, at_text) - src) + 1);
    md_vial_drawing_t d = {0};
    const bool ok = md_vial_draw(&d, &desc, at, alloc);
    md_vial_program_free(prog);
    return ok ? d.num_prims : UINT32_MAX;
}

UTEST(vial_ir, angle_draws_every_form) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    EXPECT_EQ(3u, vial_ir__prims("r = residue[2]\na = angle(r[2], r[8], r[9])\n", "angle(", alloc));
    EXPECT_EQ(3u, vial_ir__prims("r = residue[2]\na = angle(r[[2, 8, 9]])\n", "angle(", alloc));
    EXPECT_EQ(3u, vial_ir__prims("a = angle(particle[[0, 1, 2]])\n", "angle(", alloc));
    EXPECT_EQ(9u, vial_ir__prims("rs = residue[2:5]\na = angle(rs.particle[[2, 8, 9]])\n", "angle(", alloc));
    EXPECT_EQ(9u, vial_ir__prims("rs = residue[2:5]\na = angle(rs['N'], rs['CA'], rs['C'])\n", "angle(", alloc));
    EXPECT_EQ(9u, vial_ir__prims("rs = residue[2:5]\nd = dihedral(rs.particle[[0, 2, 4, 6]])\n", "dihedral(", alloc));     // two disks and an arc each
    // Inside a comprehension's body: the comprehension, every element
    EXPECT_EQ(9u, vial_ir__prims("a = [angle(r[2], r[8], r[9]) for r in residue[2:5]]\n", "angle(", alloc));
    EXPECT_EQ(9u, vial_ir__prims("a = [angle(r[2], r[8], r[9]) for r in residue[2:5]]\n", "r[8]", alloc));
    md_arena_allocator_destroy(alloc);
}

// The atoms hover highlights at at_text: their count, and each one summed, to tell two selections apart
static bool vial_ir__atoms(const char* src, const char* at_text, md_allocator_i* alloc, size_t* groups, size_t* count, uint64_t* sum) {
    const md_system_t* sys = vial_test_system('P');
    const md_system_state_t* st = vial_test_state('P');
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys}, alloc);
    const md_vial_diagnostic_t* diags;
    const size_t nd = md_vial_diagnostics(prog, &diags);
    for (size_t k = 0; k < nd; ++k) printf("%.5s %.*s\n", diags[k].code, (int)diags[k].message.len, diags[k].message.ptr);
    if (nd) { md_vial_program_free(prog); return false; }
    md_vial_draw_desc_t desc = {.prog = prog, .sys = sys, .frame = 0, .state = st, .element = -1};
    const uint32_t at = md_vial_expr_at(prog, (uint32_t)(strstr(src, at_text) - src) + 1);
    md_vial_drawing_t d = {0};
    const bool ok = md_vial_draw(&d, &desc, at, alloc);
    md_vial_program_free(prog);
    *groups = d.num_groups;
    *count = 0;
    *sum = 0;
    for (size_t g = 0; g < d.num_groups; ++g) {
        for (uint32_t i = d.groups[g].beg; i < d.groups[g].end; ++i) { ++*count; *sum += d.atoms[i]; }
    }
    return ok;
}

// A literal is shown through what it selects: 'ALA' in residue('ALA') is those residues, 'CB' in rs['CB'] or
// rs.particle('CB') those particles, 5 in within(5, of=x) what is within, and 'CB' in rs[['O', 'CB']] its column.
// Elsewhere, a literal is itself: nothing is drawn
UTEST(vial_ir, literal_draws_its_context) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    size_t g0, n0, g, n;
    uint64_t s0, s;
    const char* src = "rs = residue('ALA')\ncz = rs['CB']\ncols = rs[['O', 'CB']]\nw = within(3.0, of=cz)\nq = percentile(distance(rs.com, cz), 90)\n";
    ASSERT_TRUE(vial_ir__atoms(src, "residue(", alloc, &g0, &n0, &s0));
    ASSERT_TRUE(g0 > 1 && n0 > g0);
    ASSERT_TRUE(vial_ir__atoms(src, "'ALA'", alloc, &g, &n, &s));
    EXPECT_EQ(g0, g); EXPECT_EQ(n0, n); EXPECT_EQ(s0, s);
    ASSERT_TRUE(vial_ir__atoms(src, "cz)", alloc, &g0, &n0, &s0));
    EXPECT_EQ(g, g0); EXPECT_EQ(g0, n0);
    ASSERT_TRUE(vial_ir__atoms(src, "'CB']\n", alloc, &g, &n, &s));
    EXPECT_EQ(g0, g); EXPECT_EQ(n0, n); EXPECT_EQ(s0, s);
    // One column of rs[['O', 'CB']]: the CB of every ALA, as rs['CB'] has them, and the O the other
    ASSERT_TRUE(vial_ir__atoms(src, "'CB']]", alloc, &g, &n, &s));
    EXPECT_EQ(g0, g); EXPECT_EQ(n0, n); EXPECT_EQ(s0, s);
    ASSERT_TRUE(vial_ir__atoms(src, "'O'", alloc, &g, &n, &s));
    EXPECT_EQ(g0, g); EXPECT_EQ(n0, n); EXPECT_NE(s0, s);
    ASSERT_TRUE(vial_ir__atoms(src, "[['O'", alloc, &g, &n, &s));
    EXPECT_EQ(2 * g0, g); EXPECT_EQ(2 * n0, n);
    ASSERT_TRUE(vial_ir__atoms(src, "within(", alloc, &g0, &n0, &s0));
    ASSERT_TRUE(vial_ir__atoms(src, "3.0", alloc, &g, &n, &s));
    EXPECT_EQ(g0, g); EXPECT_EQ(n0, n); EXPECT_EQ(s0, s);
    ASSERT_TRUE(vial_ir__atoms(src, "90", alloc, &g, &n, &s));
    EXPECT_EQ(0u, g);
    // One of a list argument: the selection with that item alone, so 'CB' in particle(['N', 'CB']) is the CB atoms, as
    // rs['CB'] has them (in one set, as particle() gives), by position or by keyword; the list itself is both
    const char* lsrc = "rs = residue('ALA')\ncb = rs['CB']\nps = particle(['N', 'CB'])\nqs = particle(name=['N', 'CB'])\n";
    ASSERT_TRUE(vial_ir__atoms(lsrc, "cb = rs", alloc, &g0, &n0, &s0));
    EXPECT_EQ(15u, g0);
    ASSERT_TRUE(vial_ir__atoms(lsrc, "'CB'])\nqs", alloc, &g, &n, &s));
    EXPECT_EQ(1u, g); EXPECT_EQ(n0, n); EXPECT_EQ(s0, s);
    ASSERT_TRUE(vial_ir__atoms(lsrc, "'CB'])\n\0", alloc, &g, &n, &s));
    EXPECT_EQ(1u, g); EXPECT_EQ(n0, n); EXPECT_EQ(s0, s);
    ASSERT_TRUE(vial_ir__atoms(lsrc, "'N', 'CB'])\nqs", alloc, &g, &n, &s));
    EXPECT_EQ(1u, g); EXPECT_EQ(n0, n); EXPECT_NE(s0, s);
    ASSERT_TRUE(vial_ir__atoms(lsrc, "particle([", alloc, &g, &n, &s));
    EXPECT_EQ(1u, g); EXPECT_EQ(2 * n0, n);
    // A list indexing: the subscript, every column, as hovering its brackets draws it; one of its items its column
    const char* isrc = "a = angle(particle[[0, 1, 2]])\n";
    ASSERT_TRUE(vial_ir__atoms(isrc, "[[0", alloc, &g, &n, &s));
    EXPECT_EQ(3u, g); EXPECT_EQ(3u, n); EXPECT_EQ(3u, s);
    ASSERT_TRUE(vial_ir__atoms(isrc, ", 2]", alloc, &g, &n, &s));
    EXPECT_EQ(3u, g); EXPECT_EQ(3u, n); EXPECT_EQ(3u, s);
    ASSERT_TRUE(vial_ir__atoms(isrc, " 1,", alloc, &g, &n, &s));
    EXPECT_EQ(1u, g); EXPECT_EQ(1u, n); EXPECT_EQ(1u, s);
    // Through what the script bound before it: the view call on rs, one CB per residue
    ASSERT_TRUE(vial_ir__atoms("rs = residue('ALA')\nvs = rs.particle(['N', 'CB'])\n", "'CB'", alloc, &g, &n, &s));
    EXPECT_EQ(g0, g); EXPECT_EQ(n0, n); EXPECT_EQ(s0, s);
    md_arena_allocator_destroy(alloc);
}

// What hover draws: the drawing shown() attached, kept through names and calls; built-ins that draw themselves; and
// values by their type
// Of more than 1,000 spheres, 1,000 are drawn: distinct, the same at the same frame, others at another, and every
// particle inside still highlighted
UTEST(vial_ir, draw_spheres_capped) {
    const md_system_t* sys = vial_test_system('T');
    ASSERT_TRUE(sys != NULL);
    const md_system_state_t* st = vial_test_state('T');
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    const char* src = "w = within(2.0, of=particle[0:3000])\nv = within(2.0, of=particle[0:10])\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys}, alloc);
    ASSERT_TRUE(prog != NULL);
    md_vial_draw_desc_t desc = {.prog = prog, .sys = sys, .frame = 0, .state = st, .element = -1};
    const uint32_t at = md_vial_expr_at(prog, (uint32_t)(strstr(src, "within(") - src) + 1);
    md_vial_drawing_t a, b, c;
    ASSERT_TRUE(md_vial_draw(&a, &desc, at, alloc));
    ASSERT_TRUE(md_vial_draw(&b, &desc, at, alloc));
    desc.frame = 1;
    ASSERT_TRUE(md_vial_draw(&c, &desc, at, alloc));
    ASSERT_EQ(1000u, a.num_prims);
    EXPECT_EQ(2000u, a.num_omitted);
    ASSERT_EQ(1u, a.num_groups);
    EXPECT_GE(a.groups[0].end - a.groups[0].beg, 3000u);
    // Distinct centres, each one of the 3,000 points
    const float* x = (const float*)st->xyz;
    uint8_t* seen = md_alloc(alloc, 3000);
    memset(seen, 0, 3000);
    uint32_t lo = 3000, hi = 0;
    for (size_t i = 0; i < a.num_prims; ++i) {
        ASSERT_EQ((int)MD_VIAL_PRIM_SPHERE, (int)a.prims[i].kind);
        uint32_t j = 0;
        while (j < 3000 && (x[3 * j] != a.prims[i].p[0][0] || x[3 * j + 1] != a.prims[i].p[0][1] || x[3 * j + 2] != a.prims[i].p[0][2])) ++j;
        ASSERT_LT(j, 3000u);
        EXPECT_EQ(0, seen[j]);
        seen[j] = 1;
        lo = MD_MIN(lo, j); hi = MD_MAX(hi, j);
    }
    // Spread over the points, not the first 1,000
    EXPECT_LT(lo, 50u);
    EXPECT_GT(hi, 2950u);
    // The same frame draws the same; another frame other spheres
    ASSERT_EQ(a.num_prims, b.num_prims);
    EXPECT_EQ(0, memcmp(a.prims, b.prims, a.num_prims * sizeof(md_vial_prim_t)));
    ASSERT_EQ(a.num_prims, c.num_prims);
    size_t same = 0;
    for (size_t i = 0; i < c.num_prims; ++i) same += memcmp(&a.prims[i], &c.prims[i], sizeof(md_vial_prim_t)) == 0;
    EXPECT_LT(same, 100u);
    // Few spheres are all drawn
    desc.frame = 0;
    ASSERT_TRUE(md_vial_draw(&a, &desc, md_vial_expr_at(prog, (uint32_t)(strstr(src, "v = within(") - src) + 5), alloc));
    EXPECT_EQ(10u, a.num_prims);
    EXPECT_EQ(0u, a.num_omitted);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

UTEST(vial_ir, draw) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    const int64_t nf = vial_test_run('P');
    ASSERT_GT(nf, 3);
    const md_system_state_t* st = vial_test_state('P');
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const char* src =
        "ca = particle(name='CA')\n"
        "d  = distance(ca[0], ca[1])\n"
        "ang = angle(ca[0], ca[1], ca[2])\n"
        "tor = dihedral(ca[0], ca[1], ca[2], ca[3])\n"
        "w  = within(4.0, of=ca[0])\n"
        "p  = ca.xyz\n"
        "n  = count(ca)\n"
        "both = line(ca[0], ca[1], text='x') + sphere(ca[2], 2.0)\n"
        "rs = residue()\n"
        "dc = distance(rs[0:3], ca[4], between='closest')\n"
        "dm = distance(rs[:, None], rs, between='closest')\n"
        "def ends(s: set):\n"
        "    u = s[0]\n"
        "    v = s[-1]\n"
        "    q = distance(u, v)\n"
        "    return shown(q, marker(u) + marker(v))\n"
        "e = ends(ca)\n"
        "f = e * 2.0\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(prog != NULL);
    {
        const md_vial_diagnostic_t* dg;
        const size_t nd = md_vial_diagnostics(prog, &dg);
        for (size_t i = 0; i < nd; ++i) if (dg[i].severity == MD_VIAL_SEVERITY_ERROR) printf("%s %.*s\n", dg[i].code, (int)dg[i].message.len, dg[i].message.ptr);
    }
    #define AT(text, plus) ((uint32_t)(strstr(src, text) - src) + (plus))
    md_vial_draw_desc_t desc = {.prog = prog, .sys = sys, .run = MD_STR_LIT("run"), .frame = 0, .state = st, .element = -1};
    md_vial_drawing_t dr;

    // A set: its particles, and they are the CA atoms
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("particle(", 9)), alloc));
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
        // In another unit's scale: Å shown as nm
        md_vial_draw_desc_t in_nm = desc;
        in_nm.value_scale = 0.1;
        ASSERT_GT(md_vial_value_text(vt, sizeof(vt), &in_nm, mk[call].expr), 0u);
        EXPECT_NEAR(dm.labels[0].value * 0.1, atof(vt), 1e-4);
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

    // angle: both arms and an arc at the vertex; dihedral: the two planes and a label
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("angle(", 5)), alloc));
    ASSERT_EQ(3u, dr.num_prims);
    EXPECT_EQ((int)MD_VIAL_PRIM_LINE, (int)dr.prims[0].kind);
    EXPECT_EQ((int)MD_VIAL_PRIM_LINE, (int)dr.prims[1].kind);
    EXPECT_EQ((int)MD_VIAL_PRIM_ARC, (int)dr.prims[2].kind);
    for (int k = 0; k < 3; ++k) {
        EXPECT_EQ(XYZ(ca[1], k), dr.prims[0].p[0][k]);
        EXPECT_NEAR(XYZ(ca[0], k), dr.prims[0].p[1][k], 1e-4f);
        EXPECT_EQ(XYZ(ca[1], k), dr.prims[1].p[0][k]);
        EXPECT_NEAR(XYZ(ca[2], k), dr.prims[1].p[1][k], 1e-4f);
        EXPECT_EQ(XYZ(ca[1], k), dr.prims[2].p[0][k]);
    }
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
    EXPECT_TRUE(md_str_eq(dr.labels[0].text, MD_STR_LIT("x")));
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
        EXPECT_EQ(MD_MIN(nr * nr, (size_t)5000), dr.num_prims);
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
        md_vial_eval_t* ev = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
        ASSERT_TRUE(ev != NULL);
        EXPECT_TRUE(md_vial_eval_run(ev));
        const md_attribute_t* ad = md_attributes_find(md_vial_eval_attributes(ev), MD_STR_LIT("script/d"));
        ASSERT_TRUE(ad != NULL);
        md_vial_draw_desc_t d3 = {.prog = prog, .sys = sys, .run = MD_STR_LIT("run"), .frame = 3, .eval = ev, .element = -1};
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
static md_vial_completions_t vial_ir__complete_on(const char* marked, const md_system_t* sys, md_allocator_i* alloc) {
    const char* bar = strchr(marked, '|');
    char* src = md_alloc(alloc, strlen(marked));
    const size_t at = (size_t)(bar - marked);
    memcpy(src, marked, at);
    strcpy(src + at, bar + 1);
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    md_vial_completions_t c;
    md_vial_complete(&c, prog, (uint32_t)at, alloc);
    // Strings are alloc's: the program may go
    md_vial_program_free(prog);
    return c;
}

static md_vial_completions_t vial_ir__complete(const char* marked, md_allocator_i* alloc) {
    return vial_ir__complete_on(marked, vial_test_system('P'), alloc);
}

static const md_vial_completion_t* vial_ir__find(const md_vial_completions_t* c, const char* label) {
    for (size_t i = 0; i < c->count; ++i) if (md_str_eq(c->items[i].label, md_str_from_cstr(label))) return &c->items[i];
    return NULL;
}

UTEST(vial_ir, editor) {
    ASSERT_TRUE(vial_test_system('P') != NULL);
    ASSERT_GT(vial_test_run('P'), 0);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));

    // Lists for highlighting, and names
    EXPECT_EQ(35u, md_vial_keywords(NULL, 0));
    md_str_t names[512];
    const size_t nb = md_vial_builtin_names(names, 512);
    ASSERT_LT(nb, 512u);
    bool has_particles = false, has_distance = false, has_hidden = false;
    for (size_t i = 0; i < nb; ++i) {
        has_particles |= md_str_eq(names[i], MD_STR_LIT("particle"));
        has_distance |= md_str_eq(names[i], MD_STR_LIT("distance"));
        has_hidden |= names[i].len && names[i].ptr[0] == '_';
    }
    EXPECT_TRUE(has_particles && has_distance && !has_hidden);
    EXPECT_TRUE(md_vial_name_valid(MD_STR_LIT("my_sel2")));
    EXPECT_TRUE(md_vial_name_valid(MD_STR_LIT("_x")));
    EXPECT_FALSE(md_vial_name_valid(MD_STR_LIT("def")));
    EXPECT_FALSE(md_vial_name_valid(MD_STR_LIT("particle")));
    EXPECT_FALSE(md_vial_name_valid(MD_STR_LIT("distance")));
    EXPECT_FALSE(md_vial_name_valid(MD_STR_LIT("2a")));
    EXPECT_FALSE(md_vial_name_valid(MD_STR_LIT("a-b")));
    EXPECT_TRUE(md_str_eq(MD_STR_LIT("energy/total"), md_vial_attr_path(MD_STR_LIT("run"), MD_STR_LIT("run/energy/total"))));
    EXPECT_EQ(0u, md_vial_attr_path(MD_STR_LIT("run"), MD_STR_LIT("runs/x")).len);

    md_vial_completions_t c;
    const md_vial_completion_t* it;
    // A name being typed: what starts so, names before keywords
    c = vial_ir__complete("ca = particle(name='CA')\nu = ce|", alloc);
    EXPECT_TRUE(md_str_eq(c.prefix, MD_STR_LIT("ce")));
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
    c = vial_ir__complete("ca = particle(name='CA')\nd = distance(ca[0], ca[1], |)\n", alloc);
    ASSERT_TRUE((it = vial_ir__find(&c, "between")) != NULL);
    EXPECT_TRUE(md_str_eq(it->text, MD_STR_LIT("between=")));
    EXPECT_EQ((int)MD_VIAL_COMPLETION_PARAMETER, (int)c.items[0].kind);
    c = vial_ir__complete("def f(s: set, scale=1.0):\n    return count(s) * scale\nx = f(particle(), sc|)\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "scale") != NULL && vial_ir__find(&c, "scale")->kind == MD_VIAL_COMPLETION_PARAMETER);
    c = vial_ir__complete("x = histogram(particle().x, axis='frame', |)\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "bins") != NULL);
    EXPECT_TRUE(vial_ir__find(&c, "axis") == NULL);
    c = vial_ir__complete("ca = particle(name='CA')\nd = angle(ca[0], ca[1], |)\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "c") != NULL);

    // A selector's keywords, those of its level
    c = vial_ir__complete("x = particle(res|)\n", alloc);
    ASSERT_TRUE((it = vial_ir__find(&c, "resname")) != NULL);
    EXPECT_EQ((int)MD_VIAL_COMPLETION_PARAMETER, (int)it->kind);
    c = vial_ir__complete("x = residue(|)\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "resname") == NULL);
    EXPECT_TRUE(vial_ir__find(&c, "id") != NULL);
    EXPECT_TRUE(vial_ir__find(&c, "gap") != NULL && vial_ir__find(&c, "resolved") != NULL && vial_ir__find(&c, "role") != NULL);
    c = vial_ir__complete("x = atom(|)\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "kind") == NULL || vial_ir__find(&c, "kind")->kind != MD_VIAL_COMPLETION_PARAMETER);    // kind() is a builtin
    EXPECT_TRUE(vial_ir__find(&c, "serial") == NULL);
    EXPECT_TRUE(vial_ir__find(&c, "hybridization") != NULL && vial_ir__find(&c, "formal_charge") != NULL);
    c = vial_ir__complete("x = fragment(|)\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "name") != NULL);

    // After a dot: properties
    c = vial_ir__complete("ca = particle(name='CA')\nd = ca.c|\n", alloc);
    ASSERT_TRUE((it = vial_ir__find(&c, "com")) != NULL);
    EXPECT_EQ((int)MD_VIAL_COMPLETION_PROPERTY, (int)it->kind);
    EXPECT_TRUE(vial_ir__find(&c, "count") == NULL || vial_ir__find(&c, "count")->kind == MD_VIAL_COMPLETION_PROPERTY);
    c = vial_ir__complete("rs = residue[1:6]\np = rs.p|\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "particle") != NULL);
    // A subscript's string: a particle's name
    c = vial_ir__complete("rs = residue[1:6]\nca = rs['C|']\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "CA") != NULL && vial_ir__find(&c, "CB") != NULL);
    c = vial_ir__complete("ca = residue[1:6][0]['|']\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "CA") != NULL);

    // In a string: the system's values for the keyword, and the run's attributes for attr()
    c = vial_ir__complete("x = particle(resname='|')\n", alloc);
    ASSERT_TRUE((it = vial_ir__find(&c, "ALA")) != NULL);
    EXPECT_EQ((int)MD_VIAL_COMPLETION_VALUE, (int)it->kind);
    c = vial_ir__complete("x = particle(name='C|')\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "CA") != NULL && vial_ir__find(&c, "CB") != NULL);
    EXPECT_TRUE(vial_ir__find(&c, "N") == NULL);
    EXPECT_TRUE(md_str_eq(c.prefix, MD_STR_LIT("C")));
    c = vial_ir__complete("x = particle('C|')\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "CA") != NULL);
    c = vial_ir__complete("x = residue(kind='a|')\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "amino_acid") != NULL);
    // The words of a closed set: a residue's roles and gaps, a ring's kinds, a particle's hybridizations
    c = vial_ir__complete("x = residue(role='n|')\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "n_terminus") != NULL && vial_ir__find(&c, "c_terminus") == NULL);
    c = vial_ir__complete("x = residue(gap='|')\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "before") != NULL && vial_ir__find(&c, "after") != NULL);
    c = vial_ir__complete("x = ring(kind='f|')\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "fused") != NULL && vial_ir__find(&c, "fused_aromatic") != NULL);
    c = vial_ir__complete("x = particle(hybridization='sp|')\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "sp2") != NULL);
    // a part's name: its description, else its kind; T's parts have none
    c = vial_ir__complete_on("x = part(name='p|')\n", vial_test_system('T'), alloc);
    EXPECT_TRUE(vial_ir__find(&c, "peptide") != NULL);
    c = vial_ir__complete("w = attr('|')\n", alloc);
    EXPECT_TRUE(vial_ir__find(&c, "time") != NULL);
    // and the names of the per-atom attributes, for attr(name, of=s)
    c = vial_ir__complete_on("o = attr('|', of=protein)\n", vial_test_system('T'), alloc);
    EXPECT_TRUE(vial_ir__find(&c, "occupancy") != NULL && vial_ir__find(&c, "b_factor") != NULL);
    md_arena_allocator_destroy(alloc);
}

// Stored selections read by name, and filters: an expression compiled alone, with names
UTEST(vial_ir, filters) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    const md_system_state_t* st = vial_test_state('P');
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const uint32_t mine[] = {2, 0, 1, 1};     // out of order, once twice: the set is 0, 1, 2
    const md_vial_named_set_t stored[] = {{MD_STR_LIT("mysel"), mine, 4}, {MD_STR_LIT("bad name"), mine, 1}};

    // In a script: a static set, counted when compiling; its atoms are in its fingerprint
    {
        const char* src = "n = count(mysel)\nm = mysel.com\n";
        md_vial_program_t* p = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .sets = stored, .num_sets = 2}, alloc);
        md_vial_hover_t h;
        ASSERT_TRUE(md_vial_hover(p, 0, &h));
        EXPECT_STREQ("3", h.value);
        const uint32_t other[] = {0, 1, 3};
        const md_vial_named_set_t moved[] = {{MD_STR_LIT("mysel"), other, 3}};
        md_vial_program_t* q = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .sets = moved, .num_sets = 1}, alloc);
        const md_vial_output_t* op;
        const md_vial_output_t* oq;
        md_vial_outputs(p, &op);
        md_vial_outputs(q, &oq);
        EXPECT_TRUE(op[1].fingerprint != oq[1].fingerprint);
        // A name the script binds hides it
        md_vial_program_t* r = md_vial_compile(MD_STR_LIT("mysel = particle(name='CA')\nn = count(mysel)\n"), &(md_vial_compile_desc_t){.sys = sys, .sets = stored, .num_sets = 1}, alloc);
        ASSERT_TRUE(md_vial_hover(r, (uint32_t)strlen("mysel = particle(name='CA')\n"), &h));
        EXPECT_STRNE("3", h.value);
        md_vial_program_free(p); md_vial_program_free(q); md_vial_program_free(r);
    }

    md_vial_program_t* prog = md_vial_compile(MD_STR_LIT("ca_s = particle(name='CA')\nnear = within(3.0, of=ca_s)\n"), &(md_vial_compile_desc_t){.sys = sys}, alloc);
    md_vial_filter_desc_t d = {.sys = sys, .state = st, .sets = stored, .num_sets = 1, .prog = prog};
    md_vial_filter_t f;
    // A predefined set, a selector's array, a stored selection, the program's static selection
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("protein"), &d, alloc));
    EXPECT_FALSE(f.dynamic);
    EXPECT_EQ(0, f.ndim);
    EXPECT_EQ(1u, f.num_sets);
    EXPECT_GT(f.offsets[1], 0u);
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("residue()"), &d, alloc));
    EXPECT_EQ(1, f.ndim);
    EXPECT_EQ(md_system_residue_count(sys), (size_t)f.dims[0]);
    EXPECT_EQ((size_t)f.dims[0], f.num_sets);
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("mysel"), &d, alloc));
    ASSERT_EQ(3u, f.offsets[1]);
    EXPECT_EQ(0u, f.atoms[0]); EXPECT_EQ(1u, f.atoms[1]); EXPECT_EQ(2u, f.atoms[2]);
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("ca_s and not mysel"), &d, alloc));
    const uint32_t nca = f.offsets[1];
    EXPECT_GT(nca, 0u);
    // Bonds: pairs of one-particle sets
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("bond(residue()[0])"), &d, alloc));
    EXPECT_EQ(2, f.ndim);
    EXPECT_EQ(2u, f.dims[1]);
    EXPECT_GT(f.dims[0], 0u);
    for (size_t i = 0; i < f.num_sets; ++i) EXPECT_EQ(1u, f.offsets[i + 1] - f.offsets[i]);
    // Dynamic: from the coordinates, and none without them
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("within(3.0, of=mysel)"), &d, alloc));
    EXPECT_TRUE(f.dynamic);
    EXPECT_GE(f.offsets[1], 3u);
    md_vial_filter_desc_t nd = d;
    nd.state = NULL;
    EXPECT_FALSE(md_vial_filter(&f, MD_STR_LIT("within(3.0, of=mysel)"), &nd, alloc));
    // A recipe of the program's: evaluated at the coordinates, so dynamic; without them it is not there
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("near and not ca_s"), &d, alloc));
    EXPECT_TRUE(f.dynamic);
    EXPECT_FALSE(md_vial_filter(&f, MD_STR_LIT("near"), &nd, alloc));
    EXPECT_FALSE(md_vial_filter(&f, MD_STR_LIT("nearby"), &d, alloc));
    // Arrays of sets the program publishes: the name gives the array, and a subscript one of it
    {
        md_vial_program_t* pa = md_vial_compile(MD_STR_LIT("res = residue()\nca_r = residue() and particle(name='CA')\nbb = bond(residue()[0])\n"),
            &(md_vial_compile_desc_t){.sys = sys}, alloc);
        md_vial_filter_desc_t ad = {.sys = sys, .state = st, .prog = pa};
        ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("res"), &ad, alloc));
        EXPECT_EQ(1, f.ndim);
        EXPECT_EQ(md_system_residue_count(sys), (size_t)f.dims[0]);
        md_vial_filter_t r3;
        ASSERT_TRUE(md_vial_filter(&r3, MD_STR_LIT("residue()[3]"), &ad, alloc));
        ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("res[3]"), &ad, alloc));
        EXPECT_EQ(0, f.ndim);
        ASSERT_EQ(r3.offsets[1], f.offsets[1]);
        for (uint32_t i = 0; i < f.offsets[1]; ++i) EXPECT_EQ(r3.atoms[i], f.atoms[i]);
        ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("ca_r[2]"), &ad, alloc));
        EXPECT_EQ(1u, f.offsets[1]);
        ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("bb"), &ad, alloc));
        EXPECT_EQ(2, f.ndim);
        EXPECT_EQ(2u, f.dims[1]);
        ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("any(res[0:2])"), &ad, alloc));
        EXPECT_EQ(0, f.ndim);

        // As a stored selection given to a script: an array, known when compiling
        const uint32_t at[] = {0, 1, 2, 5, 4};
        const uint32_t off[] = {0, 3, 5};
        const md_vial_named_set_t arr[] = {{.name = MD_STR_LIT("pair"), .atoms = at, .num_atoms = 5, .ndim = 1, .dims = {2}, .offsets = off}};
        md_vial_program_t* q = md_vial_compile(MD_STR_LIT("n = len(pair)\nm = count(pair[1])\n"), &(md_vial_compile_desc_t){.sys = sys, .sets = arr, .num_sets = 1}, alloc);
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
        ASSERT_TRUE(md_vial_filter(&r0, MD_STR_LIT("residue()[0]"), &d, alloc));
        md_vial_filter_desc_t sd = d;
        sd.scope = r0.atoms;
        sd.num_scope = r0.offsets[1];
        ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("residue()"), &sd, alloc));
        EXPECT_EQ(1u, f.dims[0]);
        EXPECT_EQ(r0.offsets[1], f.offsets[1]);
        ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("protein"), &sd, alloc));
        EXPECT_EQ(r0.offsets[1], f.offsets[1]);
        // Nothing written: the scope, whole
        ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT(""), &sd, alloc));
        EXPECT_EQ(0, f.ndim);
        ASSERT_EQ(1u, f.num_sets);
        ASSERT_EQ(r0.offsets[1], f.offsets[1]);
        for (uint32_t i = 0; i < f.offsets[1]; ++i) EXPECT_EQ(r0.atoms[i], f.atoms[i]);
    }
    // Nothing written picks every particle, as an empty md_filter does: a representation's empty filter
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT(""), &d, alloc));
    EXPECT_FALSE(f.dynamic);
    EXPECT_EQ(0, f.ndim);
    ASSERT_EQ(1u, f.num_sets);
    ASSERT_EQ((uint32_t)sys->atom.count, f.offsets[1]);
    for (uint32_t i = 0; i < f.offsets[1]; ++i) EXPECT_EQ(i, f.atoms[i]);
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("  \n "), &d, alloc));
    EXPECT_EQ((uint32_t)sys->atom.count, f.offsets[1]);
    // Errors, with spans in the expression
    EXPECT_FALSE(md_vial_filter(&f, MD_STR_LIT("1 + 2"), &d, alloc));
    EXPECT_STREQ("E0313", f.error.code);
    EXPECT_EQ(0u, f.error.span.beg);
    EXPECT_EQ(5u, f.error.span.end);
    EXPECT_FALSE(md_vial_filter(&f, MD_STR_LIT("particle(nme='CA')"), &d, alloc));
    EXPECT_GE(f.error.span.beg, 0u);
    EXPECT_LE(f.error.span.end, 19u);
    EXPECT_GT(f.error.message.len, 0u);
    EXPECT_FALSE(md_vial_filter(&f, MD_STR_LIT("protein)"), &d, alloc));

    // Completion, in the expression's own offsets
    md_vial_completions_t c;
    ASSERT_TRUE(md_vial_filter_complete(&c, MD_STR_LIT("particle(resn"), 13, &d, alloc));
    EXPECT_TRUE(vial_ir__find(&c, "resname") != NULL);
    EXPECT_EQ(9u, c.span.beg);
    EXPECT_EQ(13u, c.span.end);
    ASSERT_TRUE(md_vial_filter_complete(&c, MD_STR_LIT("mys"), 3, &d, alloc));
    EXPECT_TRUE(vial_ir__find(&c, "mysel") != NULL);
    ASSERT_TRUE(md_vial_filter_complete(&c, MD_STR_LIT("ca"), 2, &d, alloc));
    EXPECT_TRUE(vial_ir__find(&c, "ca_s") != NULL);
    ASSERT_TRUE(md_vial_filter_complete(&c, MD_STR_LIT("_"), 1, &d, alloc));
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
    const md_str_t p = md_str_from_cstr(path);
    const bool ok = strstr(file, ".gro") ? md_gro_system_init_from_file(sys, st, p, NULL) :
                    strstr(file, ".tpr") ? md_tpr_system_init_from_file(sys, st, p, NULL) : md_pdb_system_init_from_file(sys, st, p, NULL);
    return ok && md_system_infer(sys, st, NULL);
}

static uint32_t vial_ir__count(const char* expr, const md_system_t* sys, const md_system_state_t* st, md_allocator_i* alloc) {
    md_vial_filter_t f;
    if (!md_vial_filter(&f, md_str_from_cstr(expr), &(md_vial_filter_desc_t){.sys = sys, .state = st}, alloc) || f.num_sets != 1) return UINT32_MAX;
    return f.offsets[1];
}

static bool vial_ir__same_sets(const char* a, const char* b, const md_system_t* sys, const md_system_state_t* st, md_allocator_i* alloc, uint32_t* n) {
    md_vial_filter_t fa, fb;
    if (!md_vial_filter(&fa, md_str_from_cstr(a), &(md_vial_filter_desc_t){.sys = sys, .state = st}, alloc)) return false;
    if (!md_vial_filter(&fb, md_str_from_cstr(b), &(md_vial_filter_desc_t){.sys = sys, .state = st}, alloc)) return false;
    if (fa.num_sets != fb.num_sets || fa.offsets[fa.num_sets] != fb.offsets[fb.num_sets]) return false;
    for (uint32_t i = 0; i <= fa.num_sets; ++i) if (fa.offsets[i] != fb.offsets[i]) return false;
    for (uint32_t i = 0; i < fa.offsets[fa.num_sets]; ++i) if (fa.atoms[i] != fb.atoms[i]) return false;
    *n = fa.num_sets;
    return true;
}

// and/or/not in a keyword's value, known when compiling and evaluated per frame: the same sets, as the check's filter
// goes with a selection made per frame (VIAL_VAL_FILTER)
UTEST(vial_ir, selector_logic_per_frame) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(64));
    const md_system_t* sys = vial_test_system('T');
    const md_system_state_t* st = vial_test_state('T');
    ASSERT_TRUE(sys && st);
    // Per frame: every particle within 1000 Å of GTP, which is all of it but GTP
    static const struct { const char* sel; const char* args; uint32_t n; } cases[] = {
        {"residue", "role='terminus' and not 'n_terminus'", 2},                // the two C-terminal residues
        {"residue", "id=1..50 and not 10..20, chain='A'", 32},          // A has no 39 to 45
        {"residue", "gap='before' or 'after'", 2},
        {"ring", "kind='aromatic' and not 'fused'", 98},
        {"atom", "name='C*' and not 'CA', resname='PHE'", 1},                   // one set
        {"atom", "hybridization=not 'sp3', resname='TRP'", 1},
        {"residue", "'ALA' or 'GLY'", 123},
    };
    for (size_t i = 0; i < MD_COUNTOF(cases); ++i) {
        char a[256], b[256];
        snprintf(a, sizeof(a), "%s(%s)", cases[i].sel, cases[i].args);
        snprintf(b, sizeof(b), "%s(%s, scope=within(1000.0, of=particle(resname='GTP')))", cases[i].sel, cases[i].args);
        uint32_t n = 0;
        EXPECT_TRUE(vial_ir__same_sets(a, b, sys, st, alloc, &n));
        EXPECT_EQ(cases[i].n, n);
    }
    md_arena_allocator_destroy(alloc);
}

// bond()'s hydrogen bonds, found per frame, with and/or/not as with the list it stands for
UTEST(vial_ir, bond_logic_per_frame) {
    int64_t nf = 0;
    const md_system_t* sys = vial_ir__run_p(&nf);
    ASSERT_TRUE(sys != NULL);
    vial_ir__run_t r;
    const char* src =
        "a = len(bond(kind='interaction' and not 'contact'))\n"
        "b = len(bond(kind=['hydrogen', 'ionic', 'halogen']))\n"
        "c = len(bond(kind='hydrogen' or 'covalent'))\n"
        "d = len(bond(kind=['hydrogen', 'covalent']))\n";
    ASSERT_TRUE(vial_ir__start(&r, src, sys));
    double* v[4];
    for (int i = 0; i < 4; ++i) v[i] = malloc(nf * sizeof(double));
    ASSERT_TRUE(vial_ir__frames(&r, sys, nf, "a", v[0]));
    ASSERT_TRUE(vial_ir__frames(&r, sys, nf, "b", v[1]));
    ASSERT_TRUE(vial_ir__frames(&r, sys, nf, "c", v[2]));
    ASSERT_TRUE(vial_ir__frames(&r, sys, nf, "d", v[3]));
    double sum = 0.0;
    for (int64_t f = 0; f < nf; ++f) {
        EXPECT_EQ(v[1][f], v[0][f]);
        EXPECT_EQ(v[3][f], v[2][f]);
        EXPECT_GT(v[2][f], v[0][f]);
        sum += v[0][f];
    }
    EXPECT_GT(sum, 0.0);    // P's peptide has hydrogen bonds
    for (int i = 0; i < 4; ++i) free(v[i]);
    vial_ir__end(&r);
}

// s.residue[i] and the like, known when compiling and evaluated per frame: the same sets
UTEST(vial_ir, selector_views) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(64));
    const md_system_t* sys = vial_test_system('T');
    const md_system_state_t* st = vial_test_state('T');
    ASSERT_TRUE(sys && st);
    // Per frame: each chain and every particle within 1000 Å of GTP, which is all of it but GTP
    const char* near = "(chain and within(1000.0, of=particle(resname='GTP')))";
    char b[256];
    uint32_t n = 0;
    snprintf(b, sizeof(b), "%s.residue[0]", near);
    EXPECT_TRUE(vial_ir__same_sets("chain.residue[0]", b, sys, st, alloc, &n));
    EXPECT_EQ(3u, n);
    snprintf(b, sizeof(b), "%s.residue[[0, 1, 2]]", near);
    EXPECT_TRUE(vial_ir__same_sets("chain.residue[[0, 1, 2]]", b, sys, st, alloc, &n));
    EXPECT_EQ(9u, n);
    snprintf(b, sizeof(b), "%s.atom[3]", near);
    EXPECT_TRUE(vial_ir__same_sets("chain.atom[3]", b, sys, st, alloc, &n));
    EXPECT_TRUE(vial_ir__same_sets("residue('PHE').ring[0]", "(residue('PHE') and within(1000.0, of=particle(resname='GTP'))).ring[0]", sys, st, alloc, &n));
    EXPECT_EQ(42u, n);
    // A view alone on one set is the selector within it, when compiling and per frame
    EXPECT_TRUE(vial_ir__same_sets("residue(scope=chain[0])", "chain[0].residue", sys, st, alloc, &n));
    EXPECT_EQ(442u, n);
    snprintf(b, sizeof(b), "(chain[0] and within(1000.0, of=particle(resname='GTP'))).residue");
    EXPECT_TRUE(vial_ir__same_sets("residue(scope=(chain[0] and within(1000.0, of=particle(resname='GTP'))))", b, sys, st, alloc, &n));
    EXPECT_EQ(442u, n);
    EXPECT_TRUE(vial_ir__same_sets("residue(scope=chain[0])[0:2]", "chain[0].residue[0:2]", sys, st, alloc, &n));
    EXPECT_EQ(2u, n);
    // s.particle['CA'] is s['CA'] written in full, when compiling and per frame
    EXPECT_TRUE(vial_ir__same_sets("residue('PHE')['CZ']", "residue('PHE').particle['CZ']", sys, st, alloc, &n));
    EXPECT_EQ(42u, n);
    EXPECT_TRUE(vial_ir__same_sets("(residue('PHE') and within(1000.0, of=particle(resname='GTP')))[['CG', 'CZ']]",
                                   "(residue('PHE') and within(1000.0, of=particle(resname='GTP'))).atom[['CG', 'CZ']]", sys, st, alloc, &n));
    EXPECT_EQ(84u, n);
    // Positions along the first dimension of a stacked array take its rows whole (found by the fuzzer)
    const md_system_t* P = vial_test_system('P');
    const md_system_state_t* Ps = vial_test_state('P');
    ASSERT_TRUE(P && Ps);
    EXPECT_TRUE(vial_ir__same_sets("residue[2:6][['N', 'CA', 'C']]", "residue[1:6][['N', 'CA', 'C']][1:]", P, Ps, alloc, &n));
    EXPECT_EQ(12u, n);
    EXPECT_TRUE(vial_ir__same_sets("residue[[5, 3]][['N', 'CA', 'C']]", "residue[1:6][['N', 'CA', 'C']][[4, 2]]", P, Ps, alloc, &n));
    EXPECT_EQ(6u, n);
    md_arena_allocator_destroy(alloc);
}

// A particle's roles (role=) and the values of its keywords (kind(s), hybridization(s), ...): each value function picks
// what its keyword does, and the roles give what the predefined sets name
UTEST(vial_ir, particle_roles_and_values) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    const md_system_t* P = vial_test_system('P');
    const md_system_state_t* Ps = vial_test_state('P');
    ASSERT_TRUE(P && Ps);
    EXPECT_EQ(vial_ir__count("sidechain", P, Ps, alloc), vial_ir__count("particle(role='side')", P, Ps, alloc));
    EXPECT_EQ(vial_ir__count("backbone", P, Ps, alloc), vial_ir__count("particle(role='backbone')", P, Ps, alloc));
    const uint32_t n_end = vial_ir__count("particle(role='n_terminus')", P, Ps, alloc), c_end = vial_ir__count("particle(role='c_terminus')", P, Ps, alloc);
    EXPECT_EQ(4u, n_end);   // N and its three hydrogens
    EXPECT_EQ(3u, c_end);
    EXPECT_EQ(n_end + c_end, vial_ir__count("particle(role='terminus')", P, Ps, alloc));
    EXPECT_EQ(0u, vial_ir__count("particle(role='5_terminus' or '3_terminus' or 'qm')", P, Ps, alloc));
    // main: N, CA, C, O and OXT with H, HA and the N-terminal hydrogens; with side, every atom of a peptide without caps
    EXPECT_EQ(93u, vial_ir__count("particle(role='main')", P, Ps, alloc));   // 45 backbone, 16 oxygens, 17 H on N, 15 HA
    EXPECT_EQ(15u, vial_ir__count("particle(role='main', name='HA')", P, Ps, alloc));
    EXPECT_EQ(0u, vial_ir__count("particle(role='side', name='HA')", P, Ps, alloc));
    EXPECT_EQ(153u, vial_ir__count("particle(role='main' or 'side')", P, Ps, alloc));
    EXPECT_EQ(0u, vial_ir__count("particle(role='main' and 'side')", P, Ps, alloc));
    EXPECT_EQ(vial_ir__count("backbone", P, Ps, alloc), vial_ir__count("particle(role='main' and 'backbone')", P, Ps, alloc));
    EXPECT_EQ(153u, vial_ir__count("particle(role=not 'qm')", P, Ps, alloc));

    vial_ir__run_t r;
    const char* src =
        "q  = sum(formal_charge(particle))\n"
        "h  = sum(hydrogens(atom(element='C')))\n"
        "s1 = count(atom(hybridization='sp2'))\n"
        "s2 = count(particle[hybridization(particle) == 'sp2'])\n"
        "f1 = count(atom(formal_charge=not 0))\n"
        "f2 = count(particle[formal_charge(particle) != 0])\n"
        "a1 = count(particle[aromatic(particle)])\n"
        "m1 = count(particle[metallic(particle)])\n"
        "k1 = count(particle[kind(particle) == 'atom'])\n"
        "h3 = count(particle[hydrogens(particle) == 3])\n"
        "h4 = count(atom(hydrogens=3))\n";
    ASSERT_TRUE(vial_ir__start(&r, src, P));
    const bool ok = vial_interp_static(&r.in);
    if (!ok) printf("failed: %s %s\n", r.in.error.code, r.in.error.message);
    EXPECT_TRUE(ok);
    char buf[256], buf2[256];
    EXPECT_STREQ("0", vial_ir__out(&r, "q", buf, sizeof(buf)));      // the N-terminal N+ and a C-terminal O-
    EXPECT_STREQ("60", vial_ir__out(&r, "h", buf, sizeof(buf)));     // CA's one and CB's three, of fifteen alanines
    EXPECT_STREQ(vial_ir__out(&r, "s1", buf, sizeof(buf)), vial_ir__out(&r, "s2", buf2, sizeof(buf2)));
    EXPECT_STREQ("2", vial_ir__out(&r, "f1", buf, sizeof(buf)));
    EXPECT_STREQ("2", vial_ir__out(&r, "f2", buf, sizeof(buf)));
    EXPECT_STREQ("0", vial_ir__out(&r, "a1", buf, sizeof(buf)));
    EXPECT_STREQ("0", vial_ir__out(&r, "m1", buf, sizeof(buf)));
    EXPECT_STREQ("153", vial_ir__out(&r, "k1", buf, sizeof(buf)));
    EXPECT_STREQ(vial_ir__out(&r, "h3", buf, sizeof(buf)), vial_ir__out(&r, "h4", buf2, sizeof(buf2)));
    vial_ir__end(&r);

    // A nucleotide's side is its nucleoside, its sugar and base, and its ends are 5' and 3'
    md_system_t sys;
    md_system_state_t st;
    ASSERT_TRUE(vial_ir__load_file(&sys, &st, "nucl-dna.gro", alloc));
    const uint32_t base = vial_ir__count("nucleobase", &sys, &st, alloc), side = vial_ir__count("nucleoside", &sys, &st, alloc);
    EXPECT_GT(base, 0u);
    EXPECT_EQ(base, vial_ir__count("particle(role='base')", &sys, &st, alloc));
    EXPECT_EQ(side, vial_ir__count("particle(role='side') and nucleic", &sys, &st, alloc));
    EXPECT_EQ(side, vial_ir__count("particle(role='sugar' or 'base')", &sys, &st, alloc));
    EXPECT_GT(vial_ir__count("particle(role='side' and 'backbone')", &sys, &st, alloc), 0u);    // C3' and C4'
    // main: the sugar-phosphate chain with its phosphate oxygens; with side, every atom of the strands
    EXPECT_GT(vial_ir__count("particle(role='main' and not 'backbone') and nucleic", &sys, &st, alloc), 0u);
    EXPECT_EQ(0u, vial_ir__count("nucleic and not particle(role='main' or 'side')", &sys, &st, alloc));
    const uint32_t e5 = vial_ir__count("particle(role='5_terminus')", &sys, &st, alloc), e3 = vial_ir__count("particle(role='3_terminus')", &sys, &st, alloc);
    EXPECT_GT(e5, 0u);
    EXPECT_GT(e3, 0u);
    EXPECT_EQ(e5 + e3, vial_ir__count("particle(role='terminus') and nucleic", &sys, &st, alloc));
    EXPECT_EQ(0u, vial_ir__count("particle(role='n_terminus' or 'c_terminus') and nucleic", &sys, &st, alloc));
    md_system_free(&sys);
    md_arena_allocator_destroy(alloc);
}

static md_attribute_id_t vial_ir__publish_f64(md_attributes_t* t, const char* path, md_attribute_format_t fmt, md_attribute_flags_t flags,
                                               md_unit_t unit, const double* data) {
    return md_attributes_create(t, &(md_attribute_desc_t){.path = md_str_from_cstr(path), .format = fmt, .flags = flags, .unit = unit,
                                                          .data = data, .byte_size = md_attribute_element_count(&fmt) * sizeof(double)});
}

// attr() reads the run at the frame being evaluated: a series on its own axis at the row the frame's time lands on, a
// frame it has no sample at missing; a static attribute whole; the particles' per-frame values with of=, the run's
// before the system's
UTEST(vial_ir, attr_reads_the_run_per_frame) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(64));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t st = {.alloc = alloc};
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.gro"), NULL));
    ASSERT_TRUE(md_trr_system_publish_run(&sys, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.trr"), MD_STR_LIT("run/t"), MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    md_attributes_t* t = &sys.attributes;
    const md_attribute_t* frames = md_attributes_find(t, MD_STR_LIT("run/t/time"));
    ASSERT_TRUE(frames != NULL);
    const uint32_t F = frames->format.shape[0];
    ASSERT_GT(F, 3u);
    double* ft = md_alloc(alloc, F * sizeof(double));
    ASSERT_EQ((size_t)F, md_attribute_extract_f64(ft, F, t, frames, md_attribute_slice_all(), md_unit_picosecond()));
    const uint32_t N = (uint32_t)sys.atom.count;

    // A series at twice the rate of the frames, on an axis of its own: every frame lands on one of its rows
    const uint32_t R = 2 * F - 1;
    double* ot = md_alloc(alloc, R * sizeof(double));
    double* ov = md_alloc(alloc, R * sizeof(double));
    double* ow = md_alloc(alloc, R * 3 * sizeof(double));
    for (uint32_t r = 0; r < R; ++r) {
        ot[r] = r % 2 == 0 ? ft[r / 2] : 0.5 * (ft[r / 2] + ft[r / 2 + 1]);
        ov[r] = 2.0 * ot[r];
        for (int c = 0; c < 3; ++c) ow[r * 3 + c] = ot[r] + c;
    }
    const md_attribute_format_t f1 = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = {R}};
    const md_attribute_format_t f3 = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 3, .rank = 1, .shape = {R}};
    ASSERT_NE(MD_ATTRIBUTE_INVALID, vial_ir__publish_f64(t, "run/t/obs/time", f1, MD_ATTRIBUTE_FLAG_TEMPORAL, md_unit_picosecond(), ot));
    ASSERT_NE(MD_ATTRIBUTE_INVALID, vial_ir__publish_f64(t, "run/t/obs/value", f1, MD_ATTRIBUTE_FLAG_TEMPORAL, md_unit_kelvin(), ov));
    ASSERT_NE(MD_ATTRIBUTE_INVALID, vial_ir__publish_f64(t, "run/t/obs/vec", f3, MD_ATTRIBUTE_FLAG_TEMPORAL, md_unit_none(), ow));
    double* ov3 = md_alloc(alloc, R * sizeof(double));
    for (uint32_t r = 0; r < R; ++r) ov3[r] = 3.0 * ot[r];
    ASSERT_NE(MD_ATTRIBUTE_INVALID, vial_ir__publish_f64(t, "run/t/obs/value_b", f1, MD_ATTRIBUTE_FLAG_TEMPORAL, md_unit_kelvin(), ov3));
    // A series at every other frame: the frames between have no sample
    const uint32_t S = (F + 1) / 2;
    double* st_ = md_alloc(alloc, S * sizeof(double));
    for (uint32_t r = 0; r < S; ++r) st_[r] = ft[2 * r];
    const md_attribute_format_t fs = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = {S}};
    ASSERT_NE(MD_ATTRIBUTE_INVALID, vial_ir__publish_f64(t, "run/t/sparse/time", fs, MD_ATTRIBUTE_FLAG_TEMPORAL, md_unit_picosecond(), st_));
    ASSERT_NE(MD_ATTRIBUTE_INVALID, vial_ir__publish_f64(t, "run/t/sparse/v", fs, MD_ATTRIBUTE_FLAG_TEMPORAL, md_unit_picosecond(), st_));
    // A static attribute of the run
    const double kv[2] = {3.0, 4.0};
    ASSERT_NE(MD_ATTRIBUTE_INVALID, vial_ir__publish_f64(t, "run/t/k", (md_attribute_format_t){.type = MD_ATTRIBUTE_TYPE_F64, .components = 1,
        .rank = 1, .shape = {2}}, 0, md_unit_none(), kv));
    // The particles' values per frame, on the run's axis, and the system's of the same name, which the run's hides
    double* q = md_alloc(alloc, (size_t)F * N * sizeof(double));
    for (uint32_t f = 0; f < F; ++f) for (uint32_t a = 0; a < N; ++a) q[(size_t)f * N + a] = 100.0 * f + a;
    ASSERT_NE(MD_ATTRIBUTE_INVALID, vial_ir__publish_f64(t, "run/t/atom/q", (md_attribute_format_t){.type = MD_ATTRIBUTE_TYPE_F64,
        .components = 1, .rank = 2, .shape = {F, N}}, MD_ATTRIBUTE_FLAG_TEMPORAL, md_unit_elementary_charge(), q));
    float* q0 = md_alloc(alloc, N * sizeof(float));
    for (uint32_t a = 0; a < N; ++a) q0[a] = -1.0f - (float)a;
    ASSERT_NE(MD_ATTRIBUTE_INVALID, md_attributes_publish_atom_column(t, MD_STR_LIT("atom/q"), md_unit_elementary_charge(), 1, q0, N));

    const char* src = "v = attr('obs/value')\nw = attr('obs/vec')[2]\ns = attr('sparse/v')\nk = attr('k')[1]\n"
                      "q = attr('q', of=particle[0:3])[2]\nvs = attr('obs/value*')\nvb = vs[1]\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = &sys, .run = MD_STR_LIT("run/t")}, alloc);
    ASSERT_TRUE(prog != NULL);
    const md_vial_diagnostic_t* dg;
    const size_t ndg = md_vial_diagnostics(prog, &dg);
    for (size_t i = 0; i < ndg; ++i) printf("%.5s %.*s\n", dg[i].code, (int)dg[i].message.len, dg[i].message.ptr);
    // Only the sparse series may be missing
    ASSERT_EQ(1u, (unsigned)ndg);
    EXPECT_EQ(0, strncmp(dg[0].code, "W0405", 5));
    md_vial_hover_t h;
    ASSERT_TRUE(md_vial_hover(prog, (uint32_t)(strstr(src, "s = ") - src), &h));
    EXPECT_TRUE(strstr(h.type, "missing") != NULL);
    ASSERT_TRUE(md_vial_hover(prog, (uint32_t)(strstr(src, "v = ") - src), &h));
    EXPECT_TRUE(strstr(h.type, "missing") == NULL);
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = &sys, .run = MD_STR_LIT("run/t"), .executor = test_executor_of_threads(2)}, alloc);
    ASSERT_TRUE(e != NULL);
    const bool ran = md_vial_eval_run(e);
    const md_vial_diagnostic_t* ed;
    const size_t ned = md_vial_eval_diagnostics(e, &ed);
    for (size_t i = 0; i < ned; ++i) printf("%.5s %.*s\n", ed[i].code, (int)ed[i].message.len, ed[i].message.ptr);
    ASSERT_TRUE(ran);
    const md_attributes_t* out = md_vial_eval_attributes(e);
    const md_attribute_t* av = md_attributes_find(out, MD_STR_LIT("script/v"));
    const md_attribute_t* aw = md_attributes_find(out, MD_STR_LIT("script/w"));
    const md_attribute_t* as = md_attributes_find(out, MD_STR_LIT("script/s"));
    const md_attribute_t* aq = md_attributes_find(out, MD_STR_LIT("script/q"));
    const md_attribute_t* ak = md_attributes_find(out, MD_STR_LIT("script/k"));
    const md_attribute_t* ab = md_attributes_find(out, MD_STR_LIT("script/vb"));
    ASSERT_TRUE(av && aw && as && aq && ak && ab && av->data && aw->data && as->data && aq->data && ak->data && ab->data);
    ASSERT_EQ(F, av->format.shape[0]);
    for (uint32_t f = 0; f < F; ++f) {
        EXPECT_NEAR(2.0 * ft[f], ((const float*)av->data)[f], 1e-3 * (1.0 + ft[f]));
        EXPECT_NEAR(ft[f] + 2.0, ((const float*)aw->data)[f], 1e-3 * (1.0 + ft[f]));
        const float sv = ((const float*)as->data)[f];
        if (f % 2 == 0) EXPECT_NEAR(ft[f], sv, 1e-3 * (1.0 + ft[f]));
        else EXPECT_TRUE(sv != sv);     // missing, published as NaN
        EXPECT_NEAR(100.0 * f + 2.0, ((const float*)aq->data)[f], 1e-3);
        EXPECT_NEAR(3.0 * ft[f], ((const float*)ab->data)[f], 1e-3 * (1.0 + ft[f]));
    }
    EXPECT_NEAR(4.0, ((const float*)ak->data)[0], 1e-6);
    md_vial_eval_free(e);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

// at(x, time=t) reads the frame closest to t on the run's axis, as at(x, frame=k) reads frame k; hover says which frame
// and the time it has, for either
UTEST(vial_ir, at_time) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(64));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t st = {.alloc = alloc};
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.gro"), NULL));
    ASSERT_TRUE(md_trr_system_publish_run(&sys, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/tryptophan-md.trr"), MD_STR_LIT("run/t"), MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    const md_attribute_t* frames = md_attributes_find(&sys.attributes, MD_STR_LIT("run/t/time"));
    ASSERT_TRUE(frames != NULL);
    const uint32_t F = frames->format.shape[0];
    ASSERT_GT(F, 3u);
    double* ft = md_alloc(alloc, F * sizeof(double));
    ASSERT_EQ((size_t)F, md_attribute_extract_f64(ft, F, &sys.attributes, frames, md_attribute_slice_all(), md_unit_picosecond()));
    const double dt = ft[3] - ft[2];
    ASSERT_GT(dt, 0.0);

    // Frame 2 three ways: by its number, by a time a third of a frame past it in ps, and the same in ns
    char src[512], tf[64];
    snprintf(src, sizeof(src), "a = at(particle[0].xyz, frame=2)\nb = at(particle[0].xyz, time=%.6f)\nc = at(particle[0].xyz, time=%.9f ns)\n",
        ft[2] + dt / 3.0, (ft[2] - dt / 3.0) / 1000.0);
    snprintf(tf, sizeof(tf), "frame 2, t = %g ps", ft[2]);
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = &sys, .run = MD_STR_LIT("run/t"), .num_frames = F}, alloc);
    ASSERT_TRUE(prog != NULL);
    const md_vial_diagnostic_t* dg;
    const size_t ndg = md_vial_diagnostics(prog, &dg);
    for (size_t i = 0; i < ndg; ++i) printf("%.5s %.*s\n", dg[i].code, (int)dg[i].message.len, dg[i].message.ptr);
    ASSERT_EQ(0u, (unsigned)ndg);
    const char* calls[3] = {"at(particle[0].xyz, frame", "at(particle[0].xyz, time=", "at(particle[0].xyz, time=%"};
    for (int i = 0; i < 3; ++i) {
        const char* at = i < 2 ? strstr(src, calls[i]) : strstr(strstr(src, calls[1]) + 1, calls[1]);
        ASSERT_TRUE(at != NULL);
        md_vial_hover_t h;
        ASSERT_TRUE(md_vial_expr_info(prog, md_vial_expr_at(prog, (uint32_t)(at - src) + 1), &h));
        EXPECT_STREQ(tf, h.note);
    }
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = &sys, .run = MD_STR_LIT("run/t")}, alloc);
    ASSERT_TRUE(e != NULL);
    const bool ran = md_vial_eval_run(e);
    const md_vial_diagnostic_t* ed;
    const size_t ned = md_vial_eval_diagnostics(e, &ed);
    for (size_t i = 0; i < ned; ++i) printf("%.5s %.*s\n", ed[i].code, (int)ed[i].message.len, ed[i].message.ptr);
    ASSERT_TRUE(ran);
    const md_attributes_t* out = md_vial_eval_attributes(e);
    const md_attribute_t* oa = md_attributes_find(out, MD_STR_LIT("script/a"));
    const md_attribute_t* ob = md_attributes_find(out, MD_STR_LIT("script/b"));
    const md_attribute_t* oc = md_attributes_find(out, MD_STR_LIT("script/c"));
    ASSERT_TRUE(oa && ob && oc && oa->data && ob->data && oc->data);
    for (int k = 0; k < 3; ++k) {
        EXPECT_EQ(((const float*)oa->data)[k], ((const float*)ob->data)[k]);
        EXPECT_EQ(((const float*)oa->data)[k], ((const float*)oc->data)[k]);
    }
    md_vial_eval_free(e);
    md_vial_program_free(prog);

    // A time outside the run, and frame= with time= or neither
    const char* bad[3] = {"a = at(particle[0].xyz, time=1e9)\n", "a = at(particle[0].xyz, frame=1, time=2.0)\n", "a = at(particle[0].xyz)\n"};
    const char* code[3] = {"E0405", "E0204", "E0204"};
    for (int i = 0; i < 3; ++i) {
        md_vial_program_t* p = md_vial_compile(md_str_from_cstr(bad[i]), &(md_vial_compile_desc_t){.sys = &sys, .run = MD_STR_LIT("run/t"), .num_frames = F}, alloc);
        ASSERT_TRUE(p != NULL);
        const size_t n = md_vial_diagnostics(p, &dg);
        if (n) printf("%.5s %.*s\n", dg[0].code, (int)dg[0].message.len, dg[0].message.ptr);
        ASSERT_EQ(1u, (unsigned)n);
        EXPECT_EQ(0, strncmp(dg[0].code, code[i], 5));
        md_vial_program_free(p);
    }
    md_arena_allocator_destroy(alloc);
}

// bond() of a set that changes per frame is evaluated per frame, and gives what it would when compiling
UTEST(vial_ir, bonds_per_frame) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    const md_system_t* T = vial_test_system('T');
    const md_system_state_t* Ts = vial_test_state('T');
    ASSERT_TRUE(T && Ts);
    uint32_t n = 0;
    const char* near = "within(1000.0, of=particle(resname='GTP'))";
    char b[256];
    snprintf(b, sizeof(b), "bond(particle(resname='GTP') and %s)", near);
    EXPECT_TRUE(vial_ir__same_sets("bond(particle(resname='GTP'))", b, T, Ts, alloc, &n));
    EXPECT_GT(n, 0u);
    snprintf(b, sizeof(b), "bond(residue('PHE')[0] and %s, order=['aromatic'])", near);
    EXPECT_TRUE(vial_ir__same_sets("bond(residue('PHE')[0], order='aromatic')", b, T, Ts, alloc, &n));
    EXPECT_EQ(12u, n);      // 6 bonds of two particles
    snprintf(b, sizeof(b), "bond(%s, kind=['coordinate'])", near);
    EXPECT_TRUE(vial_ir__same_sets("bond(kind='coordinate')", b, T, Ts, alloc, &n));
    EXPECT_EQ(14u, n);
    md_arena_allocator_destroy(alloc);
}

// bond(kind='hydrogen') finds the hydrogen bonds of each frame from the geometry (md_hbond), of origin inferred
UTEST(vial_ir, hydrogen_bonds_per_frame) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    // Two waters, the hydrogen H1 of the first pointing at the oxygen of the second, 2.8 Å apart
    static const char pdb[] =
        "HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n"
        "HETATM    2  H1  HOH A   1       0.960   0.000   0.000  1.00  0.00           H\n"
        "HETATM    3  H2  HOH A   1      -0.240   0.930   0.000  1.00  0.00           H\n"
        "HETATM    4  O   HOH A   2       2.800   0.000   0.000  1.00  0.00           O\n"
        "HETATM    5  H1  HOH A   2       3.040   0.930   0.000  1.00  0.00           H\n"
        "HETATM    6  H2  HOH A   2       3.040  -0.930   0.000  1.00  0.00           H\n"
        "END\n";
    md_system_t sys = {.alloc = alloc};
    md_system_state_t st = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_str(&sys, &st, md_str_from_cstr(pdb), NULL));
    ASSERT_TRUE(md_system_infer(&sys, &st, NULL));
    const md_vial_filter_desc_t at0 = {.sys = &sys, .state = &st};
    md_vial_filter_t f;
    // One, H...A, the hydrogen first
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("bond(kind='hydrogen')"), &at0, alloc));
    ASSERT_EQ(2u, f.num_sets);
    EXPECT_EQ(1u, f.atoms[f.offsets[0]]);
    EXPECT_EQ(3u, f.atoms[f.offsets[1]]);
    // Its flags, as the filter found it: the system has no bond of the pair
    ASSERT_TRUE(f.pair_flags != NULL);
    EXPECT_EQ((uint32_t)md_bond_flags_set_origin((md_bond_flags_t)MD_BOND_KIND_HYDROGEN, MD_BOND_ORIGIN_INFERRED), f.pair_flags[0]);
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("particle[1]"), &at0, alloc));
    EXPECT_TRUE(f.pair_flags == NULL);
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("bond(kind='interaction', origin='inferred')"), &at0, alloc));
    EXPECT_EQ(2u, f.num_sets);
    // Both ends in the set
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("bond(residue[1], kind='hydrogen')"), &at0, alloc));
    EXPECT_EQ(0u, f.num_sets);
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("bond(particle[1] or residue[1], kind='hydrogen')"), &at0, alloc));
    EXPECT_EQ(2u, f.num_sets);
    // Of origin inferred and no order: origin='user' and order= keep them out
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("bond(kind='hydrogen', origin='user')"), &at0, alloc));
    EXPECT_EQ(0u, f.num_sets);
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("bond(kind=['hydrogen', 'covalent'], order=1)"), &at0, alloc));
    EXPECT_EQ(8u, f.num_sets);      // the four O-H
    ASSERT_TRUE(f.pair_flags != NULL);
    for (size_t i = 0; i < 4; ++i) EXPECT_EQ(MD_BOND_KIND_COVALENT, md_bond_kind((md_bond_flags_t)f.pair_flags[i]));
    md_system_free(&sys);

    // Per frame over a run: what md_hbond finds in each of P's 38 frames
    const md_system_t* P = vial_test_system('P');
    const int64_t nf = vial_test_run('P');
    ASSERT_TRUE(P && nf > 0);
    md_vial_program_t* prog = md_vial_compile(MD_STR_LIT("n = len(bond(kind='hydrogen'))\n"), &(md_vial_compile_desc_t){.sys = P, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(prog != NULL);
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = P, .run = MD_STR_LIT("run"), .executor = test_executor_of_threads(2)}, alloc);
    ASSERT_TRUE(e && md_vial_eval_run(e));
    const md_attribute_t* a = md_attributes_find(md_vial_eval_attributes(e), MD_STR_LIT("script/n"));
    ASSERT_TRUE(a && a->data);
    ASSERT_EQ((uint32_t)nf, a->format.shape[0]);
    const md_str_t paths[] = {MD_STR_LIT("atom/position"), MD_STR_LIT("unitcell")};
    md_system_extract_t* ex = md_system_extract_begin(P, MD_STR_LIT("run"), paths, 2, alloc);
    ASSERT_TRUE(ex != NULL);
    md_system_state_t fs = {.alloc = alloc};
    md_system_state_init(&fs, P->atom.count);
    int32_t total = 0;
    for (int64_t fi = 0; fi < nf; ++fi) {
        ASSERT_TRUE(md_system_extract_frame(ex, fi, &fs));
        md_hbond_set_t hs = {0};
        ASSERT_TRUE(md_hbond_compute(&hs, &(md_hbond_desc_t){0}, P, &fs, alloc));
        EXPECT_EQ((int32_t)hs.count, ((const int32_t*)a->data)[fi]);
        total += (int32_t)hs.count;
    }
    md_system_extract_end(ex);
    EXPECT_GT(total, 0);
    md_vial_eval_free(e);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

// formed() and broken() per frame: a proton moves from one water to the other, in the second of two models
UTEST(vial_ir, bonds_formed_and_broken) {
    static const char pdb[] =
        "MODEL        1\n"
        "HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n"
        "HETATM    2  H1  HOH A   1       0.960   0.000   0.000  1.00  0.00           H\n"
        "HETATM    3  H2  HOH A   1      -0.240   0.930   0.000  1.00  0.00           H\n"
        "HETATM    4  O   HOH A   2       2.800   0.000   0.000  1.00  0.00           O\n"
        "HETATM    5  H1  HOH A   2       3.040   0.930   0.000  1.00  0.00           H\n"
        "HETATM    6  H2  HOH A   2       3.040  -0.930   0.000  1.00  0.00           H\n"
        "ENDMDL\n"
        "MODEL        2\n"
        "HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n"
        "HETATM    2  H1  HOH A   1       1.800   0.000   0.000  1.00  0.00           H\n"
        "HETATM    3  H2  HOH A   1      -0.240   0.930   0.000  1.00  0.00           H\n"
        "HETATM    4  O   HOH A   2       2.800   0.000   0.000  1.00  0.00           O\n"
        "HETATM    5  H1  HOH A   2       3.040   0.930   0.000  1.00  0.00           H\n"
        "HETATM    6  H2  HOH A   2       3.040  -0.930   0.000  1.00  0.00           H\n"
        "ENDMDL\n"
        "END\n";
    const char* path = "vial_ir_proton_transfer.pdb";
    FILE* fp = fopen(path, "wb");
    ASSERT_TRUE(fp != NULL);
    fwrite(pdb, 1, sizeof(pdb) - 1, fp);
    fclose(fp);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t st = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, md_str_from_cstr(path), NULL));
    ASSERT_TRUE(md_system_infer(&sys, &st, NULL));
    ASSERT_TRUE(md_pdb_system_publish_run(&sys, md_str_from_cstr(path), MD_STR_LIT("run"), MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    const char* src =
        "f = len(formed())\nb = len(broken())\n"
        "f2 = len(formed(particle[3:6]))\nb2 = len(broken(particle[3:6]))\n"
        "b3 = len(broken(breakable='topology'))\n"
        "fh = formed()[0, 0]\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = &sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(prog != NULL);
    const md_vial_diagnostic_t* dg;
    const size_t ndg = md_vial_diagnostics(prog, &dg);
    for (size_t i = 0; i < ndg; ++i) printf("%.5s %.*s\n", dg[i].code, (int)dg[i].message.len, dg[i].message.ptr);
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = &sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(e != NULL);
    const bool ran = md_vial_eval_run(e);
    const md_vial_diagnostic_t* ed;
    const size_t ned = md_vial_eval_diagnostics(e, &ed);
    for (size_t i = 0; i < ned; ++i) printf("%.5s %.*s\n", ed[i].code, (int)ed[i].message.len, ed[i].message.ptr);
    ASSERT_TRUE(ran);
    const md_attributes_t* out = md_vial_eval_attributes(e);
    // Per frame: none in the first, the transfer in the second; within the second water only the bond it forms; and
    // nothing breaks where only a topology's bonds may
    static const struct { const char* name; int32_t f0, f1; } want[] = {{"script/f", 0, 1}, {"script/b", 0, 1}, {"script/f2", 0, 1},
                                                                       {"script/b2", 0, 0}, {"script/b3", 0, 0}};
    for (size_t i = 0; i < MD_COUNTOF(want); ++i) {
        const md_attribute_t* a = md_attributes_find(out, md_str_from_cstr(want[i].name));
        ASSERT_TRUE(a && a->data);
        ASSERT_EQ(2u, a->format.shape[0]);
        const int32_t* v = (const int32_t*)a->data;
        EXPECT_EQ(want[i].f0, v[0]);
        EXPECT_EQ(want[i].f1, v[1]);
    }
    md_vial_eval_free(e);
    md_vial_program_free(prog);
    // A filter at the second frame: the bond formed has its flags as formed() found it, the bond broken the system's
    {
        const md_str_t paths[] = {MD_STR_LIT("atom/position"), MD_STR_LIT("unitcell")};
        md_system_extract_t* ex = md_system_extract_begin(&sys, MD_STR_LIT("run"), paths, 2, alloc);
        ASSERT_TRUE(ex != NULL);
        md_system_state_t fs = {.alloc = alloc};
        md_system_state_init(&fs, sys.atom.count);
        ASSERT_TRUE(md_system_extract_frame(ex, 1, &fs));
        const md_vial_filter_desc_t at1 = {.sys = &sys, .state = &fs};
        md_vial_filter_t f;
        ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("formed()"), &at1, alloc));
        ASSERT_EQ(2u, f.num_sets);
        ASSERT_TRUE(f.pair_flags != NULL);
        EXPECT_EQ(MD_BOND_KIND_COVALENT, md_bond_kind((md_bond_flags_t)f.pair_flags[0]));
        EXPECT_EQ(MD_BOND_ORIGIN_INFERRED, md_bond_origin((md_bond_flags_t)f.pair_flags[0]));
        ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("broken()"), &at1, alloc));
        ASSERT_EQ(2u, f.num_sets);
        ASSERT_TRUE(f.pair_flags != NULL);
        EXPECT_EQ((uint32_t)md_system_bond_flags(&sys, md_system_bond_find(&sys, (md_atom_idx_t)f.atoms[0], (md_atom_idx_t)f.atoms[1])), f.pair_flags[0]);
        md_system_extract_end(ex);
    }
    // Values formed() and broken() refuse
    const char* bad[3] = {"f = formed(tolerance=0.7)\n", "f = broken(breakable='pdb')\n", "f = formed(particle, color=1)\n"};
    const char* code[3] = {"E0313", "E0206", "E0202"};
    for (int i = 0; i < 3; ++i) {
        md_vial_program_t* p = md_vial_compile(md_str_from_cstr(bad[i]), &(md_vial_compile_desc_t){.sys = &sys, .run = MD_STR_LIT("run")}, alloc);
        ASSERT_TRUE(p != NULL);
        const size_t n = md_vial_diagnostics(p, &dg);
        if (n) printf("%.5s %.*s\n", dg[0].code, (int)dg[0].message.len, dg[0].message.ptr);
        ASSERT_GE(n, 1u);
        EXPECT_EQ(0, strncmp(dg[0].code, code[i], 5));
        md_vial_program_free(p);
    }
    md_system_free(&sys);
    md_arena_allocator_destroy(alloc);
    remove(path);
}

// planarity(): the RMS distance from the plane of least spread. Four atoms of a square, alternately 0.1 Å above and
// below its plane in the first model, flat in the second.
UTEST(vial_ir, planarity) {
    static const char pdb[] =
        "MODEL        1\n"
        "HETATM    1  C1  LIG A   1       1.000   1.000   0.100  1.00  0.00           C\n"
        "HETATM    2  C2  LIG A   1      -1.000   1.000  -0.100  1.00  0.00           C\n"
        "HETATM    3  C3  LIG A   1      -1.000  -1.000   0.100  1.00  0.00           C\n"
        "HETATM    4  C4  LIG A   1       1.000  -1.000  -0.100  1.00  0.00           C\n"
        "ENDMDL\n"
        "MODEL        2\n"
        "HETATM    1  C1  LIG A   1       1.000   1.000   0.000  1.00  0.00           C\n"
        "HETATM    2  C2  LIG A   1      -1.000   1.000   0.000  1.00  0.00           C\n"
        "HETATM    3  C3  LIG A   1      -1.000  -1.000   0.000  1.00  0.00           C\n"
        "HETATM    4  C4  LIG A   1       1.000  -1.000   0.000  1.00  0.00           C\n"
        "ENDMDL\n"
        "END\n";
    const char* path = "vial_ir_planarity.pdb";
    FILE* fp = fopen(path, "wb");
    ASSERT_TRUE(fp != NULL);
    fwrite(pdb, 1, sizeof(pdb) - 1, fp);
    fclose(fp);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t st = {.alloc = alloc};
    ASSERT_TRUE(md_pdb_system_init_from_file(&sys, &st, md_str_from_cstr(path), NULL));
    ASSERT_TRUE(md_pdb_system_publish_run(&sys, md_str_from_cstr(path), MD_STR_LIT("run"), MD_RUN_FLAG_DISABLE_CACHE_WRITE));
    md_vial_program_t* prog = md_vial_compile(MD_STR_LIT("pl = planarity(particle)\n"), &(md_vial_compile_desc_t){.sys = &sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(prog != NULL);
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = &sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(e && md_vial_eval_run(e));
    const md_attributes_t* out = md_vial_eval_attributes(e);
    const md_attribute_t* a = md_attributes_find(out, MD_STR_LIT("script/pl"));
    ASSERT_TRUE(a != NULL);
    double v[2] = {0};
    ASSERT_EQ(2u, md_attribute_extract_f64(v, 2, out, a, md_attribute_slice_all(), md_unit_none()));
    EXPECT_NEAR(0.1, v[0], 1e-5);
    EXPECT_NEAR(0.0, v[1], 1e-5);
    md_vial_eval_free(e);
    md_vial_program_free(prog);
    md_system_free(&sys);
    md_arena_allocator_destroy(alloc);
    remove(path);
}

// Fragments, as when compiling so per frame: a scope that changes per frame takes the selector to evaluation
UTEST(vial_ir, fragments) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    const md_system_t* T = vial_test_system('T');
    const md_system_state_t* Ts = vial_test_state('T');
    ASSERT_TRUE(T && Ts);
    uint32_t n = 0;
    EXPECT_TRUE(vial_ir__same_sets("fragment('phenyl')", "fragment('phenyl', scope=within(1000.0, of=particle(resname='GTP')))", T, Ts, alloc, &n));
    EXPECT_EQ(75u, n);
    // A positional name per frame, as for every selector
    EXPECT_TRUE(vial_ir__same_sets("residue('TYR')", "residue('TYR', scope=within(1000.0, of=particle(resname='GTP')))", T, Ts, alloc, &n));
    EXPECT_TRUE(vial_ir__same_sets("chain('B')", "chain('B', scope=within(1000.0, of=particle(resname='GTP')))", T, Ts, alloc, &n));
    EXPECT_EQ(1u, n);
    EXPECT_TRUE(vial_ir__same_sets("residue('TYR')[0].fragment", "(residue('TYR')[0] and within(1000.0, of=particle(resname='GTP'))).fragment",
                                   T, Ts, alloc, &n));
    EXPECT_EQ(6u, n);
    // Every atom with an element is in one fragment
    EXPECT_EQ(vial_ir__count("atom", T, Ts, alloc), vial_ir__count("particle(fragment)", T, Ts, alloc));
    md_arena_allocator_destroy(alloc);
}

// A force field's type (type=, type(s)), of atoms, beads and virtual sites alike, as the topology gives it; and the
// B-factors (bfactor(s)) as the file gives them
UTEST(vial_ir, particle_types_and_bfactors) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    md_system_t sys;
    md_system_state_t st;
    // AMBER: a type is finer than an element and is no name: CA is the alpha carbon by name, an aromatic carbon by type
    ASSERT_TRUE(vial_ir__load_file(&sys, &st, "tpr/peptide_tip4p.tpr", alloc));
    EXPECT_EQ(2074u, vial_ir__count("particle(type='*')", &sys, &st, alloc));
    EXPECT_EQ(7u, vial_ir__count("particle(type='CT')", &sys, &st, alloc));
    EXPECT_EQ(4u, vial_ir__count("particle(type='CA')", &sys, &st, alloc));
    EXPECT_EQ(2u, vial_ir__count("particle(name='CA')", &sys, &st, alloc));
    EXPECT_EQ(0u, vial_ir__count("particle(name='CA', type='CA')", &sys, &st, alloc));
    EXPECT_EQ(506u, vial_ir__count("particle(type='OW*')", &sys, &st, alloc));
    EXPECT_EQ(4u, vial_ir__count("particle(type='Na' or 'Cl')", &sys, &st, alloc));
    // TIP4P's M is a virtual site with a type of its own
    EXPECT_EQ(506u, vial_ir__count("particle(type='MW')", &sys, &st, alloc));
    EXPECT_EQ(506u, vial_ir__count("particle(kind='virtual', type='MW')", &sys, &st, alloc));
    EXPECT_EQ(0u, vial_ir__count("atom(type='MW')", &sys, &st, alloc));
    EXPECT_EQ(2074u - 506u, vial_ir__count("particle(type=not 'MW')", &sys, &st, alloc));
    {
        vial_ir__run_t r;
        const char* src =
            "t1 = count(particle[type(particle) == 'CT'])\n"
            "t0 = count(particle[type(particle) == ''])\n"
            "tm = type(particle(name='MW'))[0]\n";
        ASSERT_TRUE(vial_ir__start(&r, src, &sys));
        const bool ok = vial_interp_static(&r.in);
        if (!ok) printf("failed: %s %s\n", r.in.error.code, r.in.error.message);
        EXPECT_TRUE(ok);
        char buf[256];
        EXPECT_STREQ("7", vial_ir__out(&r, "t1", buf, sizeof(buf)));
        EXPECT_STREQ("0", vial_ir__out(&r, "t0", buf, sizeof(buf)));
        EXPECT_STREQ("'MW'", vial_ir__out(&r, "tm", buf, sizeof(buf)));
        vial_ir__end(&r);
    }
    md_system_free(&sys);
    // Martini: a bead's name repeats across residues, and its type says which bead it is
    ASSERT_TRUE(vial_ir__load_file(&sys, &st, "martini/martini2_dna.tpr", alloc));
    EXPECT_EQ(24u, vial_ir__count("bead(name='SC2')", &sys, &st, alloc));
    EXPECT_EQ(8u, vial_ir__count("bead(name='SC2', type='TG2')", &sys, &st, alloc));
    EXPECT_EQ(44u, vial_ir__count("bead(type='Q*')", &sys, &st, alloc));
    EXPECT_EQ(48u, vial_ir__count("bead(type='T*' and not 'TN*')", &sys, &st, alloc));
    md_system_free(&sys);

    // T's B-factors, as the file gives them
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T);
    const float* b = md_attribute_view(md_attributes_find(&T->attributes, MD_STR_LIT("atom/b_factor")), MD_ATTRIBUTE_TYPE_F32, 1, 1);
    ASSERT_TRUE(b);
    double bmax = 0.0;
    uint32_t over = 0;
    for (size_t a = 0; a < T->atom.count; ++a) { bmax = b[a] > bmax ? b[a] : bmax; over += b[a] > 150.0f; }
    {
        vial_ir__run_t r;
        ASSERT_TRUE(vial_ir__start(&r, "m = max(bfactor(particle))\nn = count(particle[bfactor(particle) > 150.0])\n", T));
        const bool ok = vial_interp_static(&r.in);
        if (!ok) printf("failed: %s %s\n", r.in.error.code, r.in.error.message);
        EXPECT_TRUE(ok);
        char buf[256], want[64];
        snprintf(want, sizeof(want), "%u", over);
        EXPECT_GT(over, 0u);
        EXPECT_LT(over, (uint32_t)T->atom.count);
        EXPECT_STREQ(want, vial_ir__out(&r, "n", buf, sizeof(buf)));
        const double got = atof(vial_ir__out(&r, "m", buf, sizeof(buf)));
        EXPECT_NEAR(bmax, got, 1e-3 * bmax);
        vial_ir__end(&r);
    }
    // attr(name, of=s) reads the same column as bfactor(s)
    {
        vial_ir__run_t r;
        ASSERT_TRUE(vial_ir__start(&r, "s = particle(chain='A')\nn = count(s)\nm = count(s[attr('b_factor', of=s) == bfactor(s)])\n"
                                       "o = mean(attr('occupancy', of=s))\n", T));
        const bool ok = vial_interp_static(&r.in);
        if (!ok) printf("failed: %s %s\n", r.in.error.code, r.in.error.message);
        EXPECT_TRUE(ok);
        char buf[256], buf2[256];
        EXPECT_STREQ(vial_ir__out(&r, "n", buf, sizeof(buf)), vial_ir__out(&r, "m", buf2, sizeof(buf2)));
        const float* occ = md_attribute_view(md_attributes_find(&T->attributes, MD_STR_LIT("atom/occupancy")), MD_ATTRIBUTE_TYPE_F32, 1, 1);
        ASSERT_TRUE(occ);
        EXPECT_GT(atof(vial_ir__out(&r, "o", buf, sizeof(buf))), 0.0);
        vial_ir__end(&r);
    }
    // An atom the file gives no B-factor ('?') is missing, and reductions skip it
    static const char cif[] =
        "data_TEST\n"
        "loop_\n"
        "_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_alt_id\n"
        "_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n"
        "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.B_iso_or_equiv\n"
        "ATOM 1 N N . ALA A 1 1 0.000 0.000 0.000 10.0\n"
        "ATOM 2 C CA . ALA A 1 1 1.458 0.000 0.000 20.0\n"
        "ATOM 3 C C . ALA A 1 1 2.009 1.420 0.000 ?\n"
        "ATOM 4 O O . ALA A 1 1 1.251 2.390 0.000 30.0\n"
        "#\n";
    sys = (md_system_t){.alloc = alloc};
    st = (md_system_state_t){.alloc = alloc};
    ASSERT_TRUE(md_mmcif_system_init_from_str(&sys, &st, (md_str_t){cif, sizeof(cif) - 1}, NULL));
    ASSERT_TRUE(md_system_infer(&sys, &st, NULL));
    {
        vial_ir__run_t r;
        ASSERT_TRUE(vial_ir__start(&r, "b = bfactor(particle)\nn = sum(valid(b))\nm = mean(b)\n", &sys));
        const bool ok = vial_interp_static(&r.in);
        if (!ok) printf("failed: %s %s\n", r.in.error.code, r.in.error.message);
        EXPECT_TRUE(ok);
        char buf[256];
        EXPECT_STREQ("3", vial_ir__out(&r, "n", buf, sizeof(buf)));
        EXPECT_STREQ("20.0", vial_ir__out(&r, "m", buf, sizeof(buf)));
        vial_ir__end(&r);
    }
    // A per-atom attribute of components, and one of variants: dimensions after the particles', in that order
    {
        const float dip[4 * 3] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
        ASSERT_NE(MD_ATTRIBUTE_INVALID, md_attributes_publish_atom_column(&sys.attributes, MD_STR_LIT("atom/dipole"), md_unit_none(), 3, dip, 4));
        const double q[2 * 4] = {0.1, 0.2, 0.3, 0.4, -0.1, -0.2, -0.3, -0.4};    // [state][atom]
        const md_attribute_format_t fmt = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2, .shape = {2, 4}};
        ASSERT_NE(MD_ATTRIBUTE_INVALID, md_attributes_create(&sys.attributes, &(md_attribute_desc_t){.path = MD_STR_LIT("atom/state_charge"),
            .format = fmt, .unit = md_unit_elementary_charge(), .data = q, .byte_size = sizeof(q)}));
        vial_ir__run_t r;
        ASSERT_TRUE(vial_ir__start(&r, "d = attr('dipole', of=particle[1:3])\nq = attr('state_charge', of=particle[2:4])\n", &sys));
        const bool ok = vial_interp_static(&r.in);
        if (!ok) printf("failed: %s %s\n", r.in.error.code, r.in.error.message);
        EXPECT_TRUE(ok);
        char buf[256];
        EXPECT_STREQ("[[4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]", vial_ir__out(&r, "d", buf, sizeof(buf)));
        EXPECT_STREQ("[[0.3, -0.3], [0.4, -0.4]]", vial_ir__out(&r, "q", buf, sizeof(buf)));
        vial_ir__end(&r);
    }
    // A per-frame attribute of the particles in the run comes before the system's of that name: per frame, its variants
    // and components after the particles' as the system's
    {
        const double t[2] = {0.0, 1.0};
        ASSERT_NE(MD_ATTRIBUTE_INVALID, md_attributes_create(&sys.attributes, &(md_attribute_desc_t){.path = MD_STR_LIT("run/time"),
            .format = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = {2}}, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
            .unit = md_unit_picosecond(), .data = t, .byte_size = sizeof(t)}));
        const double v[2 * 4 * 3] = {0};
        ASSERT_NE(MD_ATTRIBUTE_INVALID, md_attributes_create(&sys.attributes, &(md_attribute_desc_t){.path = MD_STR_LIT("run/atom/dipole"),
            .format = {.type = MD_ATTRIBUTE_TYPE_F64, .components = 3, .rank = 2, .shape = {2, 4}}, .flags = MD_ATTRIBUTE_FLAG_TEMPORAL,
            .unit = md_unit_none(), .data = v, .byte_size = sizeof(v)}));
        md_vial_program_t* prog = md_vial_compile(MD_STR_LIT("d = attr('dipole', of=particle)\n"),
            &(md_vial_compile_desc_t){.sys = &sys, .run = MD_STR_LIT("run")}, alloc);
        ASSERT_TRUE(prog != NULL);
        const md_vial_diagnostic_t* dg;
        const size_t ndg = md_vial_diagnostics(prog, &dg);
        for (size_t i = 0; i < ndg; ++i) printf("%.5s %.*s\n", dg[i].code, (int)dg[i].message.len, dg[i].message.ptr);
        EXPECT_EQ(0u, ndg);
        md_vial_hover_t h;
        ASSERT_TRUE(md_vial_hover(prog, 0, &h));
        EXPECT_STREQ("float[4][3]", h.type);
        EXPECT_EQ((int)MD_VIAL_STAGE_FRAME, (int)h.stage);
        md_vial_program_free(prog);
    }
    md_system_free(&sys);
    md_arena_allocator_destroy(alloc);
}

UTEST(vial_ir, anatomy_sets) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    md_system_t sys;
    md_system_state_t st;
    ASSERT_TRUE(vial_ir__load_file(&sys, &st, "dppc64.pdb", alloc));
    const uint32_t lipid = vial_ir__count("lipid", &sys, &st, alloc);
    const uint32_t head = vial_ir__count("particle(head)", &sys, &st, alloc), tail = vial_ir__count("particle(tail)", &sys, &st, alloc);
    EXPECT_GT(head, 0u);
    EXPECT_GT(tail, 0u);
    EXPECT_EQ(lipid, head + tail);
    EXPECT_EQ(0u, vial_ir__count("particle(head) and particle(tail)", &sys, &st, alloc));
    EXPECT_EQ(lipid, vial_ir__count("particle(head) or particle(tail)", &sys, &st, alloc));
    EXPECT_EQ(UINT32_MAX, vial_ir__count("head", &sys, &st, alloc));     // units, not kinds: no head set
    // One head per DPPC, with its hydrogens
    {
        md_vial_filter_t h;
        ASSERT_TRUE(md_vial_filter(&h, MD_STR_LIT("head"), &(md_vial_filter_desc_t){.sys = &sys, .state = &st}, alloc));
        EXPECT_EQ(1, h.ndim);
        EXPECT_EQ(sys.lipid.head.count, (size_t)h.num_sets);
        EXPECT_EQ(sys.lipid.count, (size_t)h.num_sets);
        EXPECT_EQ(head, h.offsets[h.num_sets]);
    }
    // One set per tail; all of them are the tail particles
    md_vial_filter_t f;
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("tail"), &(md_vial_filter_desc_t){.sys = &sys, .state = &st}, alloc));
    EXPECT_EQ(1, f.ndim);
    EXPECT_EQ(sys.lipid.tail.count, (size_t)f.num_sets);
    EXPECT_EQ(tail, f.offsets[f.num_sets]);
    EXPECT_EQ(tail, vial_ir__count("particle(tail)", &sys, &st, alloc));
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("tail(residue[0])"), &(md_vial_filter_desc_t){.sys = &sys, .state = &st}, alloc));
    EXPECT_EQ(2u, (unsigned)f.num_sets);     // DPPC: two acyl chains
    {
        md_vial_program_t* p = md_vial_compile(MD_STR_LIT("l = label(tail)[1]\n"), &(md_vial_compile_desc_t){.sys = &sys}, alloc);
        md_vial_hover_t h;
        ASSERT_TRUE(md_vial_hover(p, 0, &h));
        EXPECT_TRUE(strstr(h.value, "tail2/") != NULL);
        md_vial_program_free(p);
    }

    md_system_free(&sys);
    ASSERT_TRUE(vial_ir__load_file(&sys, &st, "nucl-dna.gro", alloc));
    const uint32_t base = vial_ir__count("nucleobase", &sys, &st, alloc), side = vial_ir__count("nucleoside", &sys, &st, alloc);
    EXPECT_GT(base, 0u);
    EXPECT_GT(side, base);
    EXPECT_EQ(base, vial_ir__count("nucleobase and nucleoside", &sys, &st, alloc));
    EXPECT_EQ(0u, vial_ir__count("nucleoside and not nucleic", &sys, &st, alloc));
    md_system_free(&sys);
    md_arena_allocator_destroy(alloc);
}

// A mask over a set's particles picks those where it is true, per frame: as many as the mask counts
UTEST(vial_ir, set_mask) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    ASSERT_GT(vial_test_run('P'), 3);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(8));
    const char* src =
        "n = count(particle[particle.z > particle[0].z])\n"
        "m = sum(particle.z > particle[0].z)\n"
        "k = count(protein[(protein.x in -5.0..5.0) and (protein.y > 0.0)])\n"
        "j = sum((protein.x in -5.0..5.0) and (protein.y > 0.0))\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_EQ(0u, (unsigned)md_vial_diagnostics(prog, NULL));
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(e && md_vial_eval_run(e));
    const md_attributes_t* at = md_vial_eval_attributes(e);
    const md_attribute_t* n = md_attributes_find(at, MD_STR_LIT("script/n"));
    const md_attribute_t* m = md_attributes_find(at, MD_STR_LIT("script/m"));
    const md_attribute_t* k = md_attributes_find(at, MD_STR_LIT("script/k"));
    const md_attribute_t* j = md_attributes_find(at, MD_STR_LIT("script/j"));
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
    ASSERT_TRUE(md_vial_filter(&fl, MD_STR_LIT("particle[particle.z > 0.0]"), &(md_vial_filter_desc_t){.sys = sys, .state = vial_test_state('P')}, alloc));
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
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(8));
    const char* src =
        "ca  = particle(name='CA')\n"
        "d   = distance(ca[0], ca.xyz)\n"
        "psi = dihedral(residue() and particle(name='N'), residue() and particle(name='CA'), residue() and particle(name='C'), residue() and particle(name='O'))\n"
        "one = distance(ca[0], ca[1])\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_EQ(0u, (unsigned)md_vial_diagnostics(prog, NULL));
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(e != NULL);
    ASSERT_TRUE(md_vial_eval_run(e));
    const md_attributes_t* at = md_vial_eval_attributes(e);
    const md_attribute_t* d = md_attributes_find(at, MD_STR_LIT("script/d"));
    const md_attribute_t* dm = md_attributes_find(at, MD_STR_LIT("script/d/mean"));
    const md_attribute_t* dv = md_attributes_find(at, MD_STR_LIT("script/d/variance"));
    const md_attribute_t* dx = md_attributes_find(at, MD_STR_LIT("script/d/extent"));
    ASSERT_TRUE(d && dm && dv && dx);
    EXPECT_TRUE(md_attributes_find(at, MD_STR_LIT("script/one/mean")) == NULL);    // one value has no population
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
    const md_attribute_t* pm = md_attributes_find(at, MD_STR_LIT("script/psi/mean"));
    const md_attribute_t* pp = md_attributes_find(at, MD_STR_LIT("script/psi/period"));
    const md_attribute_t* px = md_attributes_find(at, MD_STR_LIT("script/psi/extent"));
    ASSERT_TRUE(pm && pp && px);
    const float beg = ((const float*)pp->data)[0], end = ((const float*)pp->data)[1];
    for (size_t f = 0; f < nf; ++f) {
        const float m = ((const float*)pm->data)[f], lo = ((const float*)px->data)[2 * f], hi = ((const float*)px->data)[2 * f + 1];
        EXPECT_TRUE(m >= beg && m < end);
        EXPECT_TRUE(lo >= beg && lo < end);
        EXPECT_TRUE(hi >= lo && hi - lo <= end - beg);
    }
    // The labels of the CAs, as the data names them
    const md_attribute_t* dl = md_attributes_find(at, MD_STR_LIT("script/d/labels"));
    ASSERT_TRUE(dl != NULL);
    ASSERT_EQ(15u, dl->format.shape[0]);
    EXPECT_TRUE(md_str_eq(md_attribute_str(at, dl, 0), MD_STR_LIT("CA/ALA1")));

    // An evaluation that takes d over keeps its summary
    md_vial_program_t* p2 = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    md_vial_eval_t* e2 = md_vial_eval_create_from(e, p2, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(e2 != NULL);
    const md_attribute_t* dm2 = md_attributes_find(md_vial_eval_attributes(e2), MD_STR_LIT("script/d/mean"));
    ASSERT_TRUE(dm2 != NULL);
    dm = md_attributes_find(at, MD_STR_LIT("script/d/mean"));
    EXPECT_EQ(0, memcmp(dm->data, dm2->data, nf * sizeof(float)));
    md_vial_eval_free(e2);
    md_vial_eval_free(e);
    md_vial_program_free(p2);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

// A script whose outputs an earlier evaluation holds, moved to other lines or renamed: taken over whole, it is done as
// it is made, and its plan says so; a run is not needed for it to be up to date
UTEST(vial_ir, taken_over_whole) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    const size_t nf = (size_t)vial_test_run('P');
    ASSERT_GT(nf, 3u);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const char* a = "_ca = particle(name='CA')\nd = distance(_ca[0], _ca[1])\nt = dihedral(_ca[0], _ca[1], _ca[2], _ca[3])\n";
    const char* b = "# moved down\n\n_ca = particle(name='CA')\nt = dihedral(_ca[0], _ca[1], _ca[2], _ca[3])\n\nd2 = distance(_ca[0], _ca[1])\n";
    const md_vial_compile_desc_t cd = {.sys = sys, .run = MD_STR_LIT("run")};
    const md_vial_eval_desc_t ed = {.sys = sys, .run = MD_STR_LIT("run")};
    md_vial_program_t* pa = md_vial_compile(md_str_from_cstr(a), &cd, alloc);
    md_vial_program_t* pb = md_vial_compile(md_str_from_cstr(b), &cd, alloc);
    md_vial_eval_t* ea = md_vial_eval_create(pa, &ed, alloc);
    ASSERT_TRUE(ea && pb);
    ASSERT_TRUE(md_vial_eval_run(ea));
    md_vial_eval_t* eb = md_vial_eval_create_from(ea, pb, &ed, alloc);
    ASSERT_TRUE(eb != NULL);
    const md_vial_output_t* o;
    const size_t n = md_vial_outputs(pb, &o);
    md_vial_plan_t plan[8];
    ASSERT_LE(n, 8u);
    ASSERT_TRUE(md_vial_eval_plan(plan, 8, eb, pb, nf));
    for (size_t i = 0; i < n; ++i) {
        EXPECT_EQ(MD_VIAL_PLAN_HELD, plan[i].status);
        EXPECT_EQ(0u, plan[i].num_pending);
        uint8_t held[64];
        ASSERT_LE(nf, sizeof(held));
        EXPECT_EQ(nf, md_vial_eval_output_frames(eb, o[i].name, held, sizeof(held)));
        for (size_t f = 0; f < nf; ++f) EXPECT_EQ(1, held[f]);
    }
    EXPECT_TRUE(md_attributes_find(md_vial_eval_attributes(eb), MD_STR_LIT("script/d2")) != NULL);
    // With a published selection as well, the same script compiled again: a selection is the program's, nothing for
    // the evaluation to compute, so it is done as it is made all the same
    {
        const char* s = "ca = particle(name='CA')\nd = distance(ca[0], ca[1])\n";
        md_vial_program_t* ps1 = md_vial_compile(md_str_from_cstr(s), &cd, alloc);
        md_vial_program_t* ps2 = md_vial_compile(md_str_from_cstr(s), &cd, alloc);
        md_vial_eval_t* es1 = md_vial_eval_create(ps1, &ed, alloc);
        ASSERT_TRUE(es1 && ps2);
        ASSERT_TRUE(md_vial_eval_run(es1));
        md_vial_eval_t* es2 = md_vial_eval_create_from(es1, ps2, &ed, alloc);
        ASSERT_TRUE(es2 != NULL);
        const size_t ns = md_vial_outputs(ps2, &o);
        ASSERT_TRUE(md_vial_eval_plan(plan, 8, es2, ps2, nf));
        for (size_t i = 0; i < ns; ++i) {
            if (o[i].selection) continue;
            EXPECT_EQ(MD_VIAL_PLAN_HELD, plan[i].status);
            EXPECT_EQ(0u, plan[i].num_pending);
        }
        // And taken over again from it, as an edit after an edit is
        md_vial_eval_t* es3 = md_vial_eval_create_from(es2, ps2, &ed, alloc);
        ASSERT_TRUE(es3 != NULL);
        ASSERT_TRUE(md_vial_eval_plan(plan, 8, es3, ps2, nf));
        for (size_t i = 0; i < ns; ++i) if (!o[i].selection) EXPECT_EQ(MD_VIAL_PLAN_HELD, plan[i].status);
        md_vial_eval_free(es3);
        md_vial_eval_free(es2);
        md_vial_eval_free(es1);
        md_vial_program_free(ps2);
        md_vial_program_free(ps1);
    }
    // And partly: what is new is computed by a run, and the plan says how much
    const char* c = "_ca = particle(name='CA')\nd = distance(_ca[0], _ca[1])\nnew = distance(_ca[0], _ca[2])\n";
    md_vial_program_t* pc = md_vial_compile(md_str_from_cstr(c), &cd, alloc);
    md_vial_eval_t* ec = md_vial_eval_create_from(ea, pc, &ed, alloc);
    ASSERT_TRUE(ec != NULL);
    ASSERT_TRUE(md_vial_eval_plan(plan, 8, ec, pc, nf));
    const size_t nc = md_vial_outputs(pc, &o);
    for (size_t i = 0; i < nc; ++i) if (md_str_eq(o[i].name, MD_STR_LIT("new"))) EXPECT_EQ((uint32_t)nf, plan[i].num_pending);
    md_vial_eval_free(ec);
    md_vial_eval_free(eb);
    md_vial_eval_free(ea);
    md_vial_program_free(pc);
    md_vial_program_free(pb);
    md_vial_program_free(pa);
    md_arena_allocator_destroy(alloc);
}

// A for over residues of a selection that changes per frame: residues in every frame, however many (one included)
UTEST(vial_ir, for_over_dynamic_residues) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    ASSERT_GT(vial_test_run('P'), 3);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const char* src =
        "_d = [distance(r[1], r[2]) for r in residue(within(5.0, of=residue[14]))]\n"
        "m = mean(_d)\n"
        "n = len(residue(within(5.0, of=residue[14])))\n"
        "_one = [distance(r[1], r[2]) for r in residue(within(5.0, of=residue[14]))[0:1]]\n"
        "n1 = len(_one)\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    const md_vial_diagnostic_t* dg;
    const size_t nd = md_vial_diagnostics(prog, &dg);
    for (size_t i = 0; i < nd; ++i) printf("%.5s %.*s\n", dg[i].code, (int)dg[i].message.len, dg[i].message.ptr);
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(e != NULL);
    EXPECT_TRUE(md_vial_eval_run(e));
    const size_t ne = md_vial_eval_diagnostics(e, &dg);
    for (size_t i = 0; i < ne; ++i) printf("%.5s %.*s\n", dg[i].code, (int)dg[i].message.len, dg[i].message.ptr);
    EXPECT_EQ(0u, ne);
    // As many values as residues near residue 14 in each frame, and one where one is taken
    const md_attributes_t* at = md_vial_eval_attributes(e);
    const md_attribute_t* an = md_attributes_find(at, MD_STR_LIT("script/n"));
    const md_attribute_t* a1 = md_attributes_find(at, MD_STR_LIT("script/n1"));
    ASSERT_TRUE(an && a1);
    for (uint32_t f = 0; f < an->format.shape[0]; ++f) {
        const int32_t k = ((const int32_t*)an->data)[f];
        EXPECT_GT(k, 0);
        EXPECT_EQ(1, ((const int32_t*)a1->data)[f]);
    }
    md_vial_eval_free(e);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

UTEST(vial_ir, periodic_outputs) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    ASSERT_GT(vial_test_run('P'), 3);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const char* src =
        "ca = particle(name='CA')\n"
        "t  = dihedral(ca[0], ca[1], ca[2], ca[3])\n"
        "d  = distance(ca[0], ca[1])\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    const md_vial_output_t* o;
    const size_t no = md_vial_outputs(prog, &o);
    for (size_t i = 0; i < no; ++i) {
        if (md_str_eq(o[i].name, MD_STR_LIT("t"))) EXPECT_TRUE(o[i].flags & MD_VIAL_OUTPUT_PERIODIC);
        if (md_str_eq(o[i].name, MD_STR_LIT("d"))) EXPECT_FALSE(o[i].flags & MD_VIAL_OUTPUT_PERIODIC);
    }
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(e != NULL);
    const md_attributes_t* at = md_vial_eval_attributes(e);
    const md_attribute_t* per = md_attributes_find(at, MD_STR_LIT("script/t/period"));
    const md_attribute_t* seam = md_attributes_find(at, MD_STR_LIT("script/t/seam"));
    ASSERT_TRUE(per && seam);
    EXPECT_TRUE(md_attributes_find(at, MD_STR_LIT("script/d/period")) == NULL);
    double p[2], sm;
    ASSERT_EQ(2u, md_attribute_extract_f64(p, 2, at, per, md_attribute_slice_all(), per->unit));
    ASSERT_EQ(1u, md_attribute_extract_f64(&sm, 1, at, seam, md_attribute_slice_all(), seam->unit));
    EXPECT_LT(p[0], p[1]);
    EXPECT_GE(sm, p[0]);
    EXPECT_LT(sm, p[1]);
    // In the unit the values are published in, and the values lie in the period
    EXPECT_TRUE(md_vial_eval_run(e));
    const md_attribute_t* vt = md_attributes_find(at, MD_STR_LIT("script/t"));
    ASSERT_TRUE(vt != NULL);
    per = md_attributes_find(at, MD_STR_LIT("script/t/period"));     // the table has grown: pointers into it move
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
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const char* src =
        "ca = particle(name='CA')\n"
        "x = distance(ca[0], ca[frame()])\n"
        "d = distance(ca[0], ca[1])\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    for (int threads = 1; threads <= 4; threads += 3) {
        md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run"), .executor = test_executor_of_threads((uint32_t)(threads))}, alloc);
        ASSERT_TRUE(e != NULL);
        EXPECT_TRUE(md_vial_eval_run(e));
        const md_vial_diagnostic_t* dg;
        ASSERT_EQ(1u, md_vial_eval_diagnostics(e, &dg));
        EXPECT_STREQ("R0101", dg[0].code);
        EXPECT_EQ(MD_VIAL_SEVERITY_RUNTIME, dg[0].severity);
        const md_attribute_t* ax = md_attributes_find(md_vial_eval_attributes(e), MD_STR_LIT("script/x"));
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
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
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
        "ca = particle(name='CA')\n"
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
        "one = unitless(particle[0].x * 0.0) + 1.0\n"
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
        "res = residue(chain='A', kind='amino_acid')\n"
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
        "res = residue(chain='A', kind='amino_acid')\n"
        "c = distance(res[:, None], res, between='closest')\n"
        "f = distance(res[0:40, None], res[10:60], between='farthest')\n"
        "e = distance(res[0:30], res[5:35], between='closest')\n";
    const char* src_p =
        "res = residue\n"
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
        snprintf(src, sizeof(src), "%s%s", pass ? "p = particle(name='CA').xyz\nq = p\n" :
                 "p = particle(name='CA').xyz\nq = particle(resname='VLB').xyz\n", body);
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
        "_e = distance(particle.xyz[:, None], particle[0:3].xyz)\n"
        "_f = distance(particle.xyz[:, None], particle(name='CA').xyz)\n"
        "bad = abs(len(within(4.5, of=particle[0:3])) - sum(any(_e <= 4.5, axis=1)))"
        " + abs(len(within(6.0, of=particle(name='CA'))) - sum(any(_f <= 6.0, axis=1)))"
        " + abs(len(within(0.5, of=particle[0:3])) - sum(any(_e <= 0.5, axis=1)))\n";
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

// The kernels' fast paths against their general ones, bit for bit: a flat elementwise loop (each operand the whole
// shape or one value) against the index per dimension an outer broadcast takes; a sum over trailing dimensions with
// nothing missing against one with a mask (nothing in it missing); sqrt and abs of floats; displacement of all the
// points from one against the same with a mask
UTEST(vial_ir, fast_paths_agree) {
    int64_t nf = 0;
    const md_system_t* sys = vial_ir__run_p(&nf);
    ASSERT_TRUE(sys != NULL);
    const char* src =
        "_d = distance(particle.xyz, particle[0].xyz)\n"
        "_m = missing_where(_d < -1.0, _d)\n"
        "_g = _d[:, None] * [0.37, 1.9, -2.5] - _d[:, None] / 7.0\n"
        "_x = particle.xyz\n"
        "_xm = missing_where(_x < -1.0e9, _x)\n"
        "bad = abs(sum(_d) - sum(_m)) + abs(mean(_d) - mean(_m)) + abs(std(_d) - std(_m))"
        " + sum(abs(_g[:, 1] - (_d * 1.9 - _d / 7.0)))"
        " + sum(abs(sqrt(_d) - sqrt(_m))) + sum(abs(abs(_d - 3.0) - abs(_m - 3.0)))"
        " + sum(abs(sum(_x, axis=1) - sum(_xm, axis=1)))"
        " + sum(abs(displacement(_x, particle[0].xyz) - displacement(_xm, particle[0].xyz)))"
        " + sum(abs(norm(_x) - norm(_xm)))"
        " + sum(sum(_x[:, 0:0], axis=1))\n";     // an empty reduction, which the trailing path must not walk
    vial_ir__run_t r;
    ASSERT_TRUE(vial_ir__start(&r, src, sys));
    double bad[64] = {0};
    const int64_t n = nf < 64 ? nf : 64;
    EXPECT_TRUE(vial_ir__frames(&r, sys, n, "bad", bad));
    for (int64_t f = 0; f < n; ++f) EXPECT_EQ(0.0, bad[f]);
    vial_ir__end(&r);
}

// s['CA'] every way it is evaluated: known when compiling, in a mapped def, in a comprehension, and as the set algebra it
// stands for; and a name some sets lack, missing there in each
UTEST(vial_ir, particle_names_agree) {
    int64_t nf = 0;
    const md_system_t* sys = vial_ir__run_p(&nf);
    ASSERT_TRUE(sys != NULL);
    const char* src =
        "rs = residue[1:6]\n"
        "_a = angle(rs['N'], rs['CA'], rs['C'])\n"
        "def f(r: set):\n"
        "    return angle(r['N'], r['CA'], r['C'])\n"
        "_b = f(rs)\n"
        "_c = [angle(r['N'], r['CA'], r['C']) for r in rs]\n"
        "_d = angle(rs and particle(name='N'), rs and particle(name='CA'), rs and particle(name='C'))\n"
        "def g(r: set):\n"
        "    return distance(r['H'], r['N'])\n"
        "_h = g(residue)\n"
        "_h2 = distance(residue['H'], residue['N'])\n"
        "_near = residue(within(6.0, of=residue[0]))\n"
        "_e = distance(_near['CA'], _near['N'])\n"
        "def cn(r: set):\n"
        "    return distance((r and particle(name='CA')).com, (r and particle(name='N')).com)\n"
        "_p = angle(rs.particle[1], rs.particle[0], rs.particle[2])\n"
        "_q = [angle(r[1], r[0], r[2]) for r in rs]\n"
        "def ph(r: set):\n"
        "    return angle(r.particle[1], r.particle[0], r.particle[2])\n"
        "_r = ph(rs)\n"
        "def nh(r: set):\n"
        "    return count(r.particle(element='H'))\n"
        "_s = nh(_near)\n"
        "_s2 = count(_near and particle(element='H'))\n"
        "_t1 = angle(rs[['N', 'CA', 'C']])\n"
        "_t2 = angle(rs.particle[[1, 0, 2]])\n"
        "def st(r: set):\n"
        "    return angle(r[['N', 'CA', 'C']])\n"
        "def sp(r: set):\n"
        "    return angle(r[[1, 0, 2]])\n"
        "_t3 = st(_near)\n"
        "_t4 = angle(_near['N'], _near['CA'], _near['C'])\n"
        "_t5 = sp(_near)\n"
        "_t6 = angle(_near.particle[[1, 0, 2]])\n"
        "_t7 = distance(rs[['N', 'CA']])\n"
        "_t8 = distance(rs['N'], rs['CA'])\n"
        "bad = unitless(sum(abs(_a - _b)) + sum(abs(_a - _c)) + sum(abs(_a - _d))) + unitless(abs(mean(_h) - mean(_h2)))"
        " + unitless(sum(abs(_p - _q)) + sum(abs(_p - _r))) + sum(abs(_s - _s2))"
        " + unitless(sum(abs(_t1 - _a)) + sum(abs(_t2 - _p)) + abs(sum(_t3) - sum(_t4)) + abs(sum(_t5) - sum(_t6)))"
        " + unitless(sum(abs(_t7 - _t8)))\n"
        "near = abs(sum(_e) - sum(cn(_near)))\n";
    vial_ir__run_t r;
    ASSERT_TRUE(vial_ir__start(&r, src, sys));
    double bad[64] = {0}, near[64] = {0};
    const int64_t n = nf < 64 ? nf : 64;
    EXPECT_TRUE(vial_ir__frames(&r, sys, n, "bad", bad));
    vial_ir__end(&r);
    ASSERT_TRUE(vial_ir__start(&r, src, sys));
    EXPECT_TRUE(vial_ir__frames(&r, sys, n, "near", near));
    vial_ir__end(&r);
    // The same particles: the same bits, but for a centre of one particle, which goes through a weighted mean
    for (int64_t fr = 0; fr < n; ++fr) { EXPECT_EQ(0.0, bad[fr]); EXPECT_NEAR(0.0, near[fr], 1e-4); }
}

// The closest and farthest pairs, which mdlib's kernels find, against the matrix: the same entry, and the first pair in
// row order at it. Coordinates on a grid give many exact ties, and many within an ulp.
UTEST(vial_ir, pair_extremes_match_matrix) {
    const char* body =
        "dm = distance(p[:, None], q)\n"
        "c, ci, cj = closest_pair(p, q)\n"
        "f, fi, fj = farthest_pair(p, q)\n";
    const char* sets[3] = {
        "p = particle(name='CA').xyz\nq = particle.xyz[100:400]\n",
        "p = round(particle(name='CA').xyz / 2.0)\nq = round(particle(name='CA').xyz / 2.0 + 0.25)\n",
        "p = round(particle(name='CA').xyz * 3.0) / 3.0 + 0.1\nq = round(particle.xyz[0:300] * 3.0) / 3.0\n",
    };
    // The trajectory's system is small, so its sets are every particle; its cell is periodic
    const char* traj_sets[3] = {
        "p = particle.xyz\nq = particle.xyz + 0.37\n",
        "p = round(particle.xyz / 2.0)\nq = round(particle.xyz / 2.0 + 0.25)\n",
        "p = round(particle.xyz * 2.0) / 2.0\nq = round(particle.xyz * 2.0) / 2.0 + 0.5\n",
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
    ASSERT_TRUE(vial_ir__start(&r, "n = [[a for a in r] for r in residue[0:3]]\n", vial_test_system('T')));
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
        "c = com(residue(chain='A', kind='amino_acid')[0:120])\n"
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
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(8));
    const char* src =
        "def first_x(r: body):\n    return at(r.xyz, frame=0)[0, 0]\n"
        "def drift(r: body):\n    return rmsd(r, at(r.xyz, frame=0))\n"
        "x0 = first_x(residue)\n"
        "x  = residue[1].xyz[0, 0]\n"
        "per_res = drift(residue)\n"
        "r1 = rmsd(residue[1], at(residue[1].xyz, frame=0))\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_EQ(0u, (unsigned)md_vial_diagnostics(prog, NULL));
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run"), .executor = test_executor_of_threads(3)}, alloc);
    ASSERT_TRUE(e && md_vial_eval_run(e));
    const md_attributes_t* at = md_vial_eval_attributes(e);
    const md_attribute_t* x0 = md_attributes_find(at, MD_STR_LIT("script/x0"));
    const md_attribute_t* x = md_attributes_find(at, MD_STR_LIT("script/x"));
    const md_attribute_t* pr = md_attributes_find(at, MD_STR_LIT("script/per_res"));
    const md_attribute_t* r1 = md_attributes_find(at, MD_STR_LIT("script/r1"));
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
    for (size_t i = 0; i < n; ++i) if (md_str_eq(o[i].name, md_str_from_cstr(name))) return &o[i];
    return NULL;
}

static double vial_ir__attr1(const md_attributes_t* at, const char* path, size_t i) {
    const md_attribute_t* a = md_attributes_find(at, md_str_from_cstr(path));
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
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    const char* src =
        "s  = sdf(residue[1:5], particle(name='O'), 6.0)\n"
        "sw = sdf(residue[1:5], particle(name='O'), (6.0, 4.0, 3.0), samples=2, weights=mass(particle(name='O')))\n"
        "d  = density(particle(), box=((-100.0, -100.0, -100.0), (100.0, 100.0, 100.0)), samples=0.25)\n"
        "ds = density(particle(), box=((-100.0, -100.0, -100.0), (100.0, 100.0, 100.0)), samples=0.25, sigma=6.0)\n"
        "dc = density(particle(name='CA'), samples=1)\n"
        "dm = density(particle(name='CA'), samples=1, weights='mass')\n"
        "n  = max(s)\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
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
    const md_vial_output_t* om = vial_ir__pub(prog, "dm");
    ASSERT_TRUE(om != NULL);
    EXPECT_STREQ(u8"float[n][n][n] u/Å³", om->type);     // weights='mass' as everywhere: the masses' unit per Å³
    EXPECT_EQ(MD_VIAL_STAGE_TRAJECTORY, os->stage);
    EXPECT_EQ((uint32_t)(MD_VIAL_OUTPUT_VOLUME | MD_VIAL_OUTPUT_SDF), os->flags);
    EXPECT_EQ((uint32_t)MD_VIAL_OUTPUT_VOLUME, od->flags);
    EXPECT_EQ((uint32_t)MD_VIAL_OUTPUT_VOLUME, oc->flags);

    for (int threads = 1; threads <= 4; threads += 3) {
        md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run"), .executor = test_executor_of_threads((uint32_t)(threads))}, alloc);
        ASSERT_TRUE(e != NULL);
        ASSERT_TRUE(md_vial_eval_run(e));
        const md_attributes_t* at = md_vial_eval_attributes(e);
        // The box holds every particle, so a density over it, by voxel volume, holds them all: cloud in cell keeps
        // what it spreads, and so does smoothing away from the edges
        const md_attribute_t* ad = md_attributes_find(at, MD_STR_LIT("script/d"));
        const md_attribute_t* as = md_attributes_find(at, MD_STR_LIT("script/ds"));
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
        const md_attribute_t* ac = md_attributes_find(at, MD_STR_LIT("script/dc"));
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
        const md_attribute_t* asd = md_attributes_find(at, MD_STR_LIT("script/s"));
        ASSERT_TRUE(asd != NULL);
        EXPECT_EQ(-6.0, vial_ir__attr1(at, "script/s/origin", 2));
        EXPECT_EQ(0.25, vial_ir__attr1(at, "script/s/basis", 8));
        EXPECT_EQ(4.0 * nf, vial_ir__attr1(at, "script/s/samples", 0));
        double smax = 0, stot = 0;
        for (size_t i = 0; i < 48 * 48 * 48; ++i) { smax = MD_MAX(smax, ((const float*)asd->data)[i]); stot += ((const float*)asd->data)[i]; }
        EXPECT_GT(stot, 0.0);
        const md_attribute_t* an = md_attributes_find(at, MD_STR_LIT("script/n"));
        ASSERT_TRUE(an != NULL);
        EXPECT_EQ((float)smax, ((const float*)an->data)[0]);
        // Weighted by mass: an O weighs 16, so the volume is 16 times the count, where the grids are alike
        const md_attribute_t* aw = md_attributes_find(at, MD_STR_LIT("script/sw"));
        ASSERT_TRUE(aw != NULL);
        double wtot = 0;
        for (size_t i = 0; i < 24 * 16 * 16; ++i) wtot += ((const float*)aw->data)[i];
        EXPECT_GT(wtot, 0.0);

        // Over chosen frames: every frame is what was published; its size is known without reading a frame
        uint8_t* all = md_alloc(alloc, (size_t)nf);
        memset(all, 1, (size_t)nf);
        const size_t nv = 48 * 48 * 48;
        EXPECT_EQ(nv, md_vial_eval_aggregate(e, MD_STR_LIT("s"), all, NULL, 0));
        float* over = md_alloc(alloc, nv * sizeof(float));
        ASSERT_EQ(nv, md_vial_eval_aggregate(e, MD_STR_LIT("s"), all, over, nv));
        double maxdiff = 0;
        for (size_t i = 0; i < nv; ++i) maxdiff = MD_MAX(maxdiff, fabs((double)over[i] - ((const float*)asd->data)[i]));
        EXPECT_LT(maxdiff, 1e-5 * (smax + 1));
        // The first frame alone is not the whole
        uint8_t* one = md_alloc(alloc, (size_t)nf);
        memset(one, 0, (size_t)nf);
        one[0] = 1;
        ASSERT_EQ(nv, md_vial_eval_aggregate(e, MD_STR_LIT("s"), one, over, nv));
        double osum = 0;
        for (size_t i = 0; i < nv; ++i) osum += over[i];
        EXPECT_GT(osum, 0.0);
        EXPECT_NE(osum, stot);
        md_vial_eval_free(e);
    }

    // Drawing an sdf: its four structures where the frame has them, each with the pose that places the volume
    const uint32_t expr = md_vial_expr_at(prog, (uint32_t)(strstr(src, "sdf(") - src));
    md_vial_drawing_t dr;
    ASSERT_TRUE(md_vial_draw(&dr, &(md_vial_draw_desc_t){.prog = prog, .sys = sys, .run = MD_STR_LIT("run"), .frame = 2, .element = -1}, expr, alloc));
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

// rdf(a, b, cutoff): pairs counted by a spatial search, bin by bin, normalized as gmx rdf does; checked against every
// pair, counted here under the minimum image of P's orthorhombic cell
typedef struct vial_ir__rdf_ref_t {
    double*  counts;    // [bins]
    double   pairs;     // sum over the frames of (N_a N_b - N_self) / V
} vial_ir__rdf_ref_t;

static vial_ir__rdf_ref_t vial_ir__rdf_ref(const md_system_t* sys, int64_t nf, const uint8_t* mask, int za, int zb, double lo, double hi,
                                           uint32_t bins, md_allocator_i* alloc) {
    vial_ir__rdf_ref_t r = {.counts = md_alloc(alloc, bins * sizeof(double))};
    memset(r.counts, 0, bins * sizeof(double));
    const md_str_t paths[] = {MD_STR_LIT("atom/position"), MD_STR_LIT("unitcell")};
    md_system_extract_t* ex = md_system_extract_begin(sys, MD_STR_LIT("run"), paths, 2, alloc);
    md_system_state_t st = {.alloc = alloc};
    md_system_state_init(&st, sys->atom.count);
    for (int64_t f = 0; ex && f < nf; ++f) {
        if (mask && !mask[f]) continue;
        if (!md_system_extract_frame(ex, f, &st)) break;
        const double L[3] = {st.unitcell.x, st.unitcell.y, st.unitcell.z};
        const float* xyz = (const float*)st.xyz;
        size_t na = 0, nb = 0, self = 0;
        for (size_t i = 0; i < sys->atom.count; ++i) {
            const int z = md_atom_atomic_number(&sys->atom, i);
            na += z == za;
            nb += z == zb;
            self += z == za && z == zb;
        }
        r.pairs += ((double)na * nb - (double)self) / (L[0] * L[1] * L[2]);
        for (size_t i = 0; i < sys->atom.count; ++i) {
            if (md_atom_atomic_number(&sys->atom, i) != za) continue;
            for (size_t j = 0; j < sys->atom.count; ++j) {
                if (j == i || md_atom_atomic_number(&sys->atom, j) != zb) continue;
                double d2 = 0;
                for (int c = 0; c < 3; ++c) {
                    double d = (double)xyz[3 * j + c] - xyz[3 * i + c];
                    d -= L[c] * round(d / L[c]);
                    d2 += d * d;
                }
                const double d = sqrt(d2);
                if (d < lo || d >= hi) continue;
                const uint32_t b = (uint32_t)((d - lo) / (hi - lo) * bins);
                r.counts[b < bins ? b : bins - 1] += 1;
            }
        }
    }
    if (ex) md_system_extract_end(ex);
    return r;
}

// The pairs a published rdf counted in bin k, from its value: g P · shell
static double vial_ir__rdf_count(const float* g, uint32_t k, double P, double lo, double hi, uint32_t bins) {
    const double w = (hi - lo) / bins, r0 = lo + k * w, r1 = r0 + w;
    return (double)g[k] * P * 4.0 / 3.0 * MD_PI * (r1 * r1 * r1 - r0 * r0 * r0);
}

UTEST(vial_ir, rdf) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    const int64_t nf = vial_test_run('P');
    ASSERT_GT(nf, 3);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    const char* src =
        "r  = rdf(particle(element='C'), particle(element='H'), 10.0)\n"
        "rn = rdf(particle(element='C'), particle(element='H'), 1 nm, bins=250)\n"
        "rp = rdf(particle(element='C').xyz, particle(element='H'), (0.0, 10.0), bins=250)\n"
        "rh = rdf(particle(element='H'), particle(element='H'), 2.0..8.0)\n"
        "far = rdf(particle(name='CA'), particle(name='CA'), 30.0)\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(prog != NULL);
    const md_vial_diagnostic_t* dg;
    const size_t ndg = md_vial_diagnostics(prog, &dg);
    for (size_t i = 0; i < ndg; ++i) printf("%.5s %.*s\n", dg[i].code, (int)dg[i].message.len, dg[i].message.ptr);
    ASSERT_EQ(0u, (unsigned)ndg);
    // The bins: GROMACS's width of 0.02 Å by default, or bins=; a distribution over the trajectory
    const struct { const char* name; const char* type; } want[] = {
        {"r", "float[500]"}, {"rn", "float[250]"}, {"rp", "float[250]"}, {"rh", "float[300]"}, {"far", "float[1500]"},
    };
    for (size_t i = 0; i < MD_COUNTOF(want); ++i) {
        const md_vial_output_t* o = vial_ir__pub(prog, want[i].name);
        ASSERT_TRUE(o != NULL);
        EXPECT_STREQ(want[i].type, o->type);
        EXPECT_EQ(MD_VIAL_STAGE_TRAJECTORY, o->stage);
        EXPECT_EQ((uint32_t)MD_VIAL_OUTPUT_DISTRIBUTION, o->flags);
    }
    const vial_ir__rdf_ref_t ch = vial_ir__rdf_ref(sys, nf, NULL, 6, 1, 0.0, 10.0, 500, alloc);
    const vial_ir__rdf_ref_t hh = vial_ir__rdf_ref(sys, nf, NULL, 1, 1, 2.0, 8.0, 300, alloc);
    EXPECT_GT(ch.pairs, 0.0);

    float* first = NULL;
    // Kept per frame, and with no room to keep it, so a period is made by reading its frames again; on one thread and four
    for (int run = 0; run < 4; ++run) {
        const uint32_t threads = run % 2 ? 4 : 1;
        const size_t budget = run < 2 ? 0 : 1;
        md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run"), .frame_budget = budget,
                                                                             .executor = test_executor_of_threads(threads)}, alloc);
        ASSERT_TRUE(e != NULL);
        ASSERT_TRUE(md_vial_eval_run(e));
        const md_attributes_t* at = md_vial_eval_attributes(e);
        const md_attribute_t* ar = md_attributes_find(at, MD_STR_LIT("script/r"));
        const md_attribute_t* an = md_attributes_find(at, MD_STR_LIT("script/rn"));
        const md_attribute_t* ap = md_attributes_find(at, MD_STR_LIT("script/rp"));
        const md_attribute_t* ah = md_attributes_find(at, MD_STR_LIT("script/rh"));
        const md_attribute_t* af = md_attributes_find(at, MD_STR_LIT("script/far"));
        ASSERT_TRUE(ar && an && ap && ah && af);
        ASSERT_EQ(500u, (unsigned)ar->format.shape[0]);
        const float* g = (const float*)ar->data;
        // The pairs counted are every pair's, bin by bin: a pair on a bin's edge may fall either side in 32 bits
        double total = 0, ref_total = 0, worst = 0;
        for (uint32_t k = 0; k < 500; ++k) {
            const double n = vial_ir__rdf_count(g, k, ch.pairs, 0.0, 10.0, 500);
            total += n;
            ref_total += ch.counts[k];
            worst = MD_MAX(worst, fabs(n - ch.counts[k]));
        }
        EXPECT_LT(worst, 2.01);
        EXPECT_NEAR(ref_total, total, 1e-4 * ref_total);
        // A particle is never paired with itself
        const float* gh = (const float*)ah->data;
        double hworst = 0;
        for (uint32_t k = 0; k < 300; ++k) hworst = MD_MAX(hworst, fabs(vial_ir__rdf_count(gh, k, hh.pairs, 2.0, 8.0, 300) - hh.counts[k]));
        EXPECT_LT(hworst, 2.01);
        // 1 nm is 10 Å; points pair as the set's particles do
        EXPECT_EQ(0, memcmp(an->data, ap->data, 250 * sizeof(float)));
        for (uint32_t k = 0; k < 250; ++k) {
            const double n = vial_ir__rdf_count((const float*)an->data, k, ch.pairs, 0.0, 10.0, 250);
            EXPECT_NEAR(ch.counts[2 * k] + ch.counts[2 * k + 1], n, 2.01);
        }
        // The range of its bins, beside it, in Å
        EXPECT_EQ(0.0, vial_ir__attr1(at, "script/r/range", 0));
        EXPECT_EQ(10.0, vial_ir__attr1(at, "script/r/range", 1));
        EXPECT_EQ(2.0, vial_ir__attr1(at, "script/rh/range", 0));
        // A cutoff of more than half the cell: every frame left out, with a note, and every bin missing
        for (uint32_t k = 0; k < 1500; ++k) EXPECT_TRUE(isnan(((const float*)af->data)[k]));
        const md_vial_diagnostic_t* rd;
        const size_t nr = md_vial_eval_diagnostics(e, &rd);
        bool noted = false;
        for (size_t i = 0; i < nr; ++i) noted |= !strncmp(rd[i].code, "R0101", 5) && strstr(rd[i].message.ptr, "narrowest") != NULL;
        EXPECT_TRUE(noted);
        // The same bits however many threads added the frames
        if (!first) {
            first = md_alloc(alloc, 500 * sizeof(float));
            memcpy(first, g, 500 * sizeof(float));
        } else {
            EXPECT_EQ(0, memcmp(first, g, 500 * sizeof(float)));
        }

        // Over chosen frames: every frame is what was published; the first half is the first half's pairs
        uint8_t* all = md_alloc(alloc, (size_t)nf);
        memset(all, 1, (size_t)nf);
        EXPECT_EQ(500u, (unsigned)md_vial_eval_aggregate(e, MD_STR_LIT("r"), all, NULL, 0));
        float* over = md_alloc(alloc, 500 * sizeof(float));
        ASSERT_EQ(500u, (unsigned)md_vial_eval_aggregate(e, MD_STR_LIT("r"), all, over, 500));
        EXPECT_EQ(0, memcmp(over, g, 500 * sizeof(float)));
        uint8_t* half = md_alloc(alloc, (size_t)nf);
        for (int64_t f = 0; f < nf; ++f) half[f] = f < nf / 2;
        ASSERT_EQ(500u, (unsigned)md_vial_eval_aggregate(e, MD_STR_LIT("r"), half, over, 500));
        const vial_ir__rdf_ref_t chh = vial_ir__rdf_ref(sys, nf, half, 6, 1, 0.0, 10.0, 500, alloc);
        double hw = 0;
        for (uint32_t k = 0; k < 500; ++k) hw = MD_MAX(hw, fabs(vial_ir__rdf_count(over, k, chh.pairs, 0.0, 10.0, 500) - chh.counts[k]));
        EXPECT_LT(hw, 2.01);
        md_vial_eval_free(e);
    }

    // Diagnostics: the arguments are checked when compiling
    const struct { const char* src; const char* code; } bad[] = {
        {"x = rdf(residue, particle, 10.0)\n", "E0313"},                    // an array of sets: its centres are residue.com
        {"x = rdf(particle, particle, 0.0)\n", "E0313"},                    // a cutoff of nothing
        {"x = rdf(particle, particle, 8.0..2.0)\n", "E0313"},               // hi below lo
        {"x = rdf(particle, particle, 10.0, bins=0)\n", "E0313"},
        {"x = rdf(particle, particle, distance(particle[0], particle[1]))\n", "E0313"},   // not known when compiling
        {"x = rdf(particle, particle)\n", "E0204"},
        {"x = [rdf(r, particle, 5.0) for r in residue[0:2]]\n", "E0406"},  // an array of rdfs is not supported yet
    };
    for (size_t i = 0; i < MD_COUNTOF(bad); ++i) {
        md_vial_program_t* p = md_vial_compile(md_str_from_cstr(bad[i].src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
        const size_t n = md_vial_diagnostics(p, &dg);
        bool found = false;
        for (size_t k = 0; k < n; ++k) found |= !strncmp(dg[k].code, bad[i].code, 5);
        if (!found) printf("%s: expected %s\n", bad[i].src, bad[i].code);
        EXPECT_TRUE(found);
        md_vial_program_free(p);
    }
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

// shape_weights(p): v1's weights, (l2 - l1, 2 (l1 - l0), 3 l0) / sum(l) of the mass-weighted gyration tensor about
// com, the mean of the points made whole; and gyration(p), the root of its trace. Checked against the same in double
// on P's coordinates, which are whole: the peptide in the frames where it is less than half the cell across, and every
// residue in every frame. (v1's own, md_geom_shape_weights, is in float and off by up to 5e-4 for the elongated
// peptide.)
static void vial_ir__shape_weights(const md_system_t* sys, const float* xyz, size_t beg, size_t end, double out[3], double* rg) {
    double c[3] = {0}, wsum = 0, M[9] = {0}, l[3];
    for (size_t i = beg; i < end; ++i) {
        const double m = md_atom_mass(&sys->atom, i);
        for (int q = 0; q < 3; ++q) c[q] += m * xyz[3 * i + q];
        wsum += m;
    }
    for (int q = 0; q < 3; ++q) c[q] /= wsum;
    for (size_t i = beg; i < end; ++i) {
        const double m = md_atom_mass(&sys->atom, i);
        for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) M[3 * a + b] += m * (xyz[3 * i + a] - c[a]) * (xyz[3 * i + b] - c[b]) / wsum;
    }
    if (rg) *rg = sqrt(M[0] + M[4] + M[8]);     // the radius of gyration is the root of the tensor's trace
    md_linalg_sym_eigen(M, l, 3);     // ascending
    const double t = l[0] + l[1] + l[2];
    out[0] = (l[2] - l[1]) / t;
    out[1] = 2.0 * (l[1] - l[0]) / t;
    out[2] = 3.0 * l[0] / t;
}

UTEST(vial_ir, shape_weights) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    const int64_t nf = vial_test_run('P');
    ASSERT_GT(nf, 3);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    const char* src = "w = shape_weights(particle)\nwr = shape_weights(residue)\nlin, plan, iso = shape_weights(particle)\ng = gyration(particle)\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(prog != NULL);
    const md_vial_diagnostic_t* dg;
    ASSERT_EQ(0u, (unsigned)md_vial_diagnostics(prog, &dg));
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run"), .executor = test_executor_of_threads(2)}, alloc);
    ASSERT_TRUE(e != NULL);
    ASSERT_TRUE(md_vial_eval_run(e));
    const md_attributes_t* at = md_vial_eval_attributes(e);
    const md_attribute_t* aw = md_attributes_find(at, MD_STR_LIT("script/w"));
    const md_attribute_t* ar = md_attributes_find(at, MD_STR_LIT("script/wr"));
    const md_attribute_t* al = md_attributes_find(at, MD_STR_LIT("script/lin"));
    const md_attribute_t* ag = md_attributes_find(at, MD_STR_LIT("script/g"));
    ASSERT_TRUE(aw && ar && al && ag);
    const size_t nr = sys->residue.count;
    const md_str_t paths[] = {MD_STR_LIT("atom/position"), MD_STR_LIT("unitcell")};
    md_system_extract_t* ex = md_system_extract_begin(sys, MD_STR_LIT("run"), paths, 2, alloc);
    ASSERT_TRUE(ex != NULL);
    md_system_state_t st = {.alloc = alloc};
    md_system_state_init(&st, sys->atom.count);
    int whole = 0;
    double worst = 0, rworst = 0, gworst = 0;
    for (int64_t f = 0; f < nf; ++f) {
        ASSERT_TRUE(md_system_extract_frame(ex, f, &st));
        const float* xyz = (const float*)st.xyz;
        const float* w = (const float*)aw->data + 3 * f;
        EXPECT_NEAR(1.0, (double)w[0] + w[1] + w[2], 1e-5);
        EXPECT_EQ(w[0], ((const float*)al->data)[f]);
        // Less than half the cell across: the minimum image places it whole
        float lo[3] = {INFINITY, INFINITY, INFINITY}, hi[3] = {-INFINITY, -INFINITY, -INFINITY};
        for (size_t i = 0; i < sys->atom.count; ++i) for (int q = 0; q < 3; ++q) { lo[q] = MD_MIN(lo[q], xyz[3 * i + q]); hi[q] = MD_MAX(hi[q], xyz[3 * i + q]); }
        const double L[3] = {st.unitcell.x, st.unitcell.y, st.unitcell.z};
        if (hi[0] - lo[0] < L[0] / 2 && hi[1] - lo[1] < L[1] / 2 && hi[2] - lo[2] < L[2] / 2) {
            double ref[3], rg;
            vial_ir__shape_weights(sys, xyz, 0, sys->atom.count, ref, &rg);
            for (int q = 0; q < 3; ++q) worst = MD_MAX(worst, fabs((double)w[q] - ref[q]));
            gworst = MD_MAX(gworst, fabs((double)((const float*)ag->data)[f] - rg));
            whole += 1;
        }
        for (size_t r = 0; r < nr; ++r) {
            const md_urange_t rg = md_system_residue_atom_range(sys, r);
            double ref[3];
            vial_ir__shape_weights(sys, xyz, (size_t)rg.beg, (size_t)rg.end, ref, NULL);
            const float* wr = (const float*)ar->data + (size_t)f * nr * 3 + 3 * r;
            for (int q = 0; q < 3; ++q) rworst = MD_MAX(rworst, fabs((double)wr[q] - ref[q]));
        }
    }
    md_system_extract_end(ex);
    EXPECT_GT(whole, 0);
    EXPECT_LT(worst, 1e-4);
    EXPECT_LT(rworst, 1e-4);
    EXPECT_LT(gworst, 1e-4);    // about com, the mean of the peptide made whole
    md_vial_eval_free(e);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

// missing_where(c, x) marks entries of a value of its own: marking zeroes what is under the mask, which in x's data
// would change x for every other reader of it, as filled(m, x) and the next filled did
UTEST(vial_ir, missing_where_keeps_its_input) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    const int64_t nf = vial_test_run('P');
    ASSERT_GT(nf, 3);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(8));
    const char* src =
        "d = distance(particle[0], particle[30])\n"
        "_m = missing_where(d > 0.0, d)\n"
        "back = filled(_m, d)\n"
        "ten = filled(_m, 10.0)\n"
        "again = filled(missing_where(d > 0.0, d), 10.0)\n"
        "after = d + 0.0\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(prog != NULL);
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(e != NULL);
    ASSERT_TRUE(md_vial_eval_run(e));
    const md_attributes_t* at = md_vial_eval_attributes(e);
    const md_attribute_t* ad = md_attributes_find(at, MD_STR_LIT("script/d"));
    const md_attribute_t* ab = md_attributes_find(at, MD_STR_LIT("script/back"));
    const md_attribute_t* ar = md_attributes_find(at, MD_STR_LIT("script/after"));
    const md_attribute_t* a10 = md_attributes_find(at, MD_STR_LIT("script/ten"));
    const md_attribute_t* ag = md_attributes_find(at, MD_STR_LIT("script/again"));
    ASSERT_TRUE(ad && ab && ar && a10 && ag);
    for (int64_t f = 0; f < nf; ++f) {
        const float d = ((const float*)ad->data)[f];
        EXPECT_GT(d, 0.0f);
        EXPECT_EQ(d, ((const float*)ab->data)[f]);
        EXPECT_EQ(d, ((const float*)ar->data)[f]);
        EXPECT_EQ(10.0f, ((const float*)a10->data)[f]);
        EXPECT_EQ(10.0f, ((const float*)ag->data)[f]);
    }
    md_vial_eval_free(e);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

// Values carry no unit when evaluating, so a quantity in another unit is converted where it meets others: 0.5 nm as
// a radius is 5 Å. Each pair below is one value written in nm and in Å, and must give the same.
UTEST(vial_ir, units_convert) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    const int64_t nf = vial_test_run('P');
    ASSERT_GT(nf, 3);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    const char* src =
        "d = distance(particle[0], particle[30])\n"
        "_qn = unitless(particle[0:5].xyz) * 0.1 nm\n"
        "_qa = particle[0:5].xyz\n"
        "w_nm = count(within(0.5 nm, of=particle[0]))\n"            "w_a = count(within(5.0, of=particle[0]))\n"
        "wo_nm = count(within(5.0, of=_qn))\n"                        "wo_a = count(within(5.0, of=_qa))\n"
        "c_nm = clip(d, 0.1 nm, 0.5 nm)\n"                            "c_a = clip(d, 1.0, 5.0)\n"
        "mx_nm = max(d, 1 nm)\n"                                      "mx_a = max(d, 10.0)\n"
        "mn_nm = min(d, 0.2 nm)\n"                                    "mn_a = min(d, 2.0)\n"
        "wh_nm = where(d > 100.0, d, 1 nm)\n"                         "wh_a = where(d > 100.0, d, 10.0)\n"
        "f_nm = filled(missing_where(d > 0.0, d), 1 nm)\n"            "f_a = filled(missing_where(d > 0.0, d), 10.0)\n"
        "p_nm = distance(_qn[0], particle[30])\n"                     "p_a = distance(_qa[0], particle[30])\n"
        "g_nm = gyration(_qn, weights=None)\n"                        "g_a = gyration(_qa, weights=None)\n"
        "in_nm = d in 0.5 nm..1.5 nm\n"                               "in_a = d in 5.0..15.0\n"
        "l_nm = [d, 1 nm]\n"                                          "l_a = [d, 10.0]\n"
        "k_nm = concatenate((_qa, _qn))[5:, 0]\n"                     "k_a = _qa[:, 0]\n"
        "h_nm = histogram(d, axis='frame', bins=10, range=(0.0, 2 nm))\n"   "h_a = histogram(d, axis='frame', bins=10, range=(0.0, 20.0))\n"
        "s_nm = sdf(residue[1:5], particle(name='O'), 0.6 nm, sigma=0.05 nm)\n"   "s_a = sdf(residue[1:5], particle(name='O'), 6.0, sigma=0.5)\n"
        "dn_nm = density(particle, box=((0.0, 0.0, 0.0), (4 nm, 4 nm, 4 nm)), samples=1)\n"
        "dn_a = density(particle, box=((0.0, 0.0, 0.0), (40.0, 40.0, 40.0)), samples=1)\n"
        "fit_nm = fit(_qn, at(_qn, frame=0))\n"
        "sp = sphere(particle[0], 0.3 nm)\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(prog != NULL);
    const md_vial_diagnostic_t* dg;
    const size_t ndg = md_vial_diagnostics(prog, &dg);
    for (size_t i = 0; i < ndg; ++i) printf("%.5s %.*s\n", dg[i].code, (int)dg[i].message.len, dg[i].message.ptr);
    ASSERT_EQ(0u, (unsigned)ndg);
    // Positions fitted are in Å, as geometry gives them, whatever unit the points came in
    const md_vial_output_t* of = vial_ir__pub(prog, "fit_nm");
    ASSERT_TRUE(of != NULL);
    EXPECT_STREQ("float[5][3] Å", of->type);
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(e != NULL);
    ASSERT_TRUE(md_vial_eval_run(e));
    const md_attributes_t* at = md_vial_eval_attributes(e);
    static const char* pairs[] = {"w", "wo", "c", "mx", "mn", "wh", "f", "p", "g", "in", "l", "k", "h", "s", "dn"};
    for (size_t i = 0; i < MD_COUNTOF(pairs); ++i) {
        char pn[64], pa[64];
        snprintf(pn, sizeof(pn), "script/%s_nm", pairs[i]);
        snprintf(pa, sizeof(pa), "script/%s_a", pairs[i]);
        const md_attribute_t* an = md_attributes_find(at, md_str_from_cstr(pn));
        const md_attribute_t* aa = md_attributes_find(at, md_str_from_cstr(pa));
        ASSERT_TRUE(an && aa);
        const size_t n = md_attribute_slice_count(an, md_attribute_slice_all());
        ASSERT_EQ(n, md_attribute_slice_count(aa, md_attribute_slice_all()));
        double* vn = md_alloc(alloc, n * sizeof(double));
        double* va = md_alloc(alloc, n * sizeof(double));
        ASSERT_EQ(n, md_attribute_extract_f64(vn, n, at, an, md_attribute_slice_all(), md_unit_none()));
        ASSERT_EQ(n, md_attribute_extract_f64(va, n, at, aa, md_attribute_slice_all(), md_unit_none()));
        double worst = 0, scale = 0;
        for (size_t k = 0; k < n; ++k) { worst = MD_MAX(worst, fabs(vn[k] - va[k])); scale = MD_MAX(scale, fabs(va[k])); }
        if (worst > 1e-5 * (scale + 1)) printf("%s: %s and %s differ by %g\n", pairs[i], pn, pa, worst);
        EXPECT_LE(worst, 1e-5 * (scale + 1));
    }
    // A sphere's radius in nm is drawn in Å
    const uint32_t expr = md_vial_expr_at(prog, (uint32_t)(strstr(src, "sphere(") - src));
    md_vial_drawing_t dr;
    ASSERT_TRUE(md_vial_draw(&dr, &(md_vial_draw_desc_t){.prog = prog, .sys = sys, .run = MD_STR_LIT("run"), .frame = 1, .element = -1}, expr, alloc));
    ASSERT_EQ(1u, (unsigned)dr.num_prims);
    EXPECT_NEAR(3.0, dr.prims[0].r, 1e-6);
    md_vial_eval_free(e);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

// Inside a def mapped over residues, an extent that changes from call to call keeps its identity through slices by
// fixed offsets, as ? does: r[1:] - r[:-1] and w[:, None] * d evaluate, per residue
UTEST(vial_ir, slices_of_extents_per_call) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    const int64_t nf = vial_test_run('P');
    ASSERT_GT(nf, 3);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(8));
    const char* src =
        "def steps(r: points):\n"
        "    return mean(norm(r[1:] - r[:-1]))\n"
        "st = steps(residue)\n"
        "def second_moment(p: points):\n"
        "    w = weight(p, None)\n"
        "    d = displacement(com(p, weights=None), p)\n"
        "    return trace(d.T @ (w[:, None] * d)) / sum(w)\n"
        "sp = second_moment(residue)\n"
        "g = gyration(residue, weights=None)\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(prog != NULL);
    const md_vial_diagnostic_t* dg;
    const size_t ndg = md_vial_diagnostics(prog, &dg);
    for (size_t i = 0; i < ndg; ++i) printf("%.5s %.*s\n", dg[i].code, (int)dg[i].message.len, dg[i].message.ptr);
    ASSERT_EQ(0u, (unsigned)ndg);
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    ASSERT_TRUE(e != NULL);
    ASSERT_TRUE(md_vial_eval_run(e));
    const md_attributes_t* at = md_vial_eval_attributes(e);
    const md_attribute_t* ast = md_attributes_find(at, MD_STR_LIT("script/st"));
    const md_attribute_t* asp = md_attributes_find(at, MD_STR_LIT("script/sp"));
    const md_attribute_t* ag = md_attributes_find(at, MD_STR_LIT("script/g"));
    ASSERT_TRUE(ast && asp && ag);
    const size_t nr = sys->residue.count;
    const md_str_t paths[] = {MD_STR_LIT("atom/position"), MD_STR_LIT("unitcell")};
    md_system_extract_t* ex = md_system_extract_begin(sys, MD_STR_LIT("run"), paths, 2, alloc);
    ASSERT_TRUE(ex != NULL);
    md_system_state_t st = {.alloc = alloc};
    md_system_state_init(&st, sys->atom.count);
    double sworst = 0, gworst = 0;
    for (int64_t f = 0; f < nf; ++f) {
        ASSERT_TRUE(md_system_extract_frame(ex, f, &st));
        const float* xyz = (const float*)st.xyz;
        for (size_t r = 0; r < nr; ++r) {
            // The mean step from one particle of the residue to the next, as P's residues are whole
            const md_urange_t rg = md_system_residue_atom_range(sys, r);
            double sum = 0;
            for (size_t i = rg.beg + 1; i < rg.end; ++i) {
                double d2 = 0;
                for (int q = 0; q < 3; ++q) d2 += ((double)xyz[3 * i + q] - xyz[3 * (i - 1) + q]) * ((double)xyz[3 * i + q] - xyz[3 * (i - 1) + q]);
                sum += sqrt(d2);
            }
            sworst = MD_MAX(sworst, fabs(sum / (rg.end - rg.beg - 1) - ((const float*)ast->data)[f * nr + r]));
            // The trace of the spread is the square of the radius of gyration
            const double g = ((const float*)ag->data)[f * nr + r];
            gworst = MD_MAX(gworst, fabs(g * g - ((const float*)asp->data)[f * nr + r]));
        }
    }
    md_system_extract_end(ex);
    EXPECT_LT(sworst, 1e-4);
    EXPECT_LT(gworst, 1e-4);
    md_vial_eval_free(e);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

UTEST(vial_ir, volume_diagnostics) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    vial_test_run('P');
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(8));
    struct { const char* src; const char* code; } cases[] = {
        {"v = sdf(residue[0:5], particle(name='O'), 6.0)\n", "E0320"},           // the first residue has three Hs on N
        {"v = sdf(residue[1:5], particle(name='O'), 300.0)\n", "W0407"},         // 2,400 a side at 4 per Å
        {"v = [sdf(r, particle(name='O'), 6.0) for r in residue[1:3]]\n", "E0406"},
        {"v = sdf(residue[1:5], residue, 6.0)\n", "E0313"},                      // an array of sets is no target
        {"v = sdf(residue[1:5], particle(name='O'), 6.0, samples=0)\n", "E0313"},
        {"v = density(particle(), box=((0, 0, 0), (0, 10, 10)))\n", "E0313"},
        {"x = 4.0\nv = sdf(residue[1:5], particle(name='O'), 6.0, ref=particle(name='O').xyz)\n", "E0313"},
    };
    for (size_t i = 0; i < MD_COUNTOF(cases); ++i) {
        md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(cases[i].src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
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

// mcs(*structures): which particles of each play the parts of which of the first, set[k][w], typed when compiling
UTEST(vial_ir, mcs) {
    const md_system_t* sys = vial_test_system('T');
    ASSERT_TRUE(sys != NULL);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(8));
    const char* src =
        "gtp = particle(resname='GTP')\n"
        "gdp = particle(resname='GDP')\n"
        "m = mcs(gtp, gdp)\n"
        "r = rmsd(m[1].xyz, m[0].xyz)\n"
        "phe = residue('PHE')\n"
        "tyr = residue('TYR')\n"
        "c = mcs(phe[0], tyr[0], match='element')\n"
        "aro = mcs(phe[0:4], tyr[0:4], residue('TRP')[0:2])\n"
        "one = mcs(phe)\n"
        "res = mcs(residue(chain='A')[0:3])\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys}, alloc);
    const md_vial_diagnostic_t* d;
    const size_t nd = md_vial_diagnostics(prog, &d);
    for (size_t i = 0; i < nd; ++i) printf("%.5s %.*s\n", d[i].code, (int)d[i].message.len, d[i].message.ptr);
    EXPECT_EQ(0u, (unsigned)nd);
    md_vial__result_t r;
    ASSERT_TRUE(md_vial__line_result(prog, 3, &r));
    EXPECT_STREQ("set[2][28]", r.types[0]);
    ASSERT_TRUE(md_vial__line_result(prog, 4, &r));
    EXPECT_STREQ("float Å", r.types[0]);
    ASSERT_TRUE(md_vial__line_result(prog, 7, &r));
    EXPECT_STREQ("set[2][11]", r.types[0]);       // PHE and TYR: the ring, CB and CA... eleven heavy atoms alike
    ASSERT_TRUE(md_vial__line_result(prog, 8, &r));
    EXPECT_STREQ("set[10][10]", r.types[0]);      // with TRP: what all three share
    ASSERT_TRUE(md_vial__line_result(prog, 9, &r));
    EXPECT_STREQ("set[42][11]", r.types[0]);      // one array of alike structures: every heavy atom, one mapping shared
    ASSERT_TRUE(md_vial__line_result(prog, 10, &r));
    EXPECT_STREQ("set[3][6]", r.types[0]);        // MET, ARG, GLU: backbone and CB, CG
    // Of the particles, as filters read them: GDP's in the order of GTP's they play
    md_vial_filter_t f;
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("mcs(particle(resname='GTP'), particle(resname='GDP'))[1]"), &(md_vial_filter_desc_t){.sys = sys}, alloc));
    EXPECT_EQ(28u, (unsigned)f.num_sets);
    ASSERT_TRUE(md_vial_filter(&f, MD_STR_LIT("mcs(particle(resname='GTP'), particle(resname='GDP'))"), &(md_vial_filter_desc_t){.sys = sys}, alloc));
    EXPECT_EQ(56u, (unsigned)f.num_sets);     // both rows, flat
    // Too little in common: E0321, unless least= asks for less
    struct { const char* src; const char* code; } bad[] = {
        {"m = mcs(particle(resname='GTP'), particle(resname='VLB'))\n", "E0321"},
        {"m = mcs(particle(resname='GTP'), particle(resname='SO4'))\n", "E0321"},
        {"m = mcs(particle(resname='GTP'), particle(resname='GDP'), match='colour')\n", "E0313"},
        {"m = mcs(particle(resname='GTP'), particle(resname='GDP'), by='element')\n", "E0202"},
    };
    for (size_t i = 0; i < MD_COUNTOF(bad); ++i) {
        md_vial_program_t* p = md_vial_compile(md_str_from_cstr(bad[i].src), &(md_vial_compile_desc_t){.sys = sys}, alloc);
        const size_t n = md_vial_diagnostics(p, &d);
        bool found = false;
        for (size_t k = 0; k < n; ++k) found = found || !strncmp(d[k].code, bad[i].code, 5);
        EXPECT_TRUE(found);
        md_vial_program_free(p);
    }
    md_vial_program_t* ok = md_vial_compile(MD_STR_LIT("m = mcs(particle(resname='GTP'), particle(resname='VLB'), least=0.2)\n"), &(md_vial_compile_desc_t){.sys = sys}, alloc);
    EXPECT_EQ(0u, (unsigned)md_vial_diagnostics(ok, NULL));
    md_vial_program_free(ok);
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}

static double vial_ir__max_diff(const md_attributes_t* at, const char* p, const char* q, size_t* n_out) {
    const md_attribute_t* a = md_attributes_find(at, md_str_from_cstr(p));
    const md_attribute_t* b = md_attributes_find(at, md_str_from_cstr(q));
    if (!a || !b || !a->data || !b->data) return INFINITY;
    const size_t n = md_attribute_element_count(&a->format);
    if (n != md_attribute_element_count(&b->format)) return INFINITY;
    double d = 0;
    for (size_t i = 0; i < n; ++i) d = MD_MAX(d, fabs((double)((const float*)a->data)[i] - ((const float*)b->data)[i]));
    if (n_out) *n_out = n;
    return d;
}

// An sdf over structures paired column by column: what mcs gives for one pair, and for arrays of alike structures,
// which it maps column for column, so structures of different kinds join through concatenate
UTEST(vial_ir, sdf_columns) {
    const md_system_t* sys = vial_test_system('P');
    ASSERT_TRUE(sys != NULL);
    ASSERT_GT(vial_test_run('P'), 3);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    const char* src =
        "ala = residue[1:3] and not particle(element='H')\n"
        "cb  = particle(name='CB')\n"
        "s   = sdf(ala, cb, 5.0)\n"
        "m   = mcs(residue[1], residue[2])\n"
        "so  = sdf(m, cb, 5.0)\n"
        "rv  = sdf(m[:, ::-1], cb, 5.0)\n"
        "m2  = mcs(residue[1:3], residue[1:3])\n"
        "s2  = sdf(m2, cb, 5.0)\n"
        "ct  = sdf(mcs(residue[1:5], residue[14:15]), cb, 5.0)\n"
        "c1  = sdf(mcs(residue[1:5]), cb, 5.0)\n";
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = MD_STR_LIT("run")}, alloc);
    const md_vial_diagnostic_t* dg;
    const size_t ndg = md_vial_diagnostics(prog, &dg);
    for (size_t i = 0; i < ndg; ++i) printf("%.5s %.*s\n", dg[i].code, (int)dg[i].message.len, dg[i].message.ptr);
    ASSERT_EQ(0u, (unsigned)ndg);
    // mcs over arrays: set[k][w], a row per structure
    md_vial__result_t r;
    ASSERT_TRUE(md_vial__line_result(prog, 4, &r));
    EXPECT_STREQ("set[2][5]", r.types[0]);
    ASSERT_TRUE(md_vial__line_result(prog, 7, &r));
    EXPECT_STREQ("set[4][5]", r.types[0]);
    md_vial_eval_t* e = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = MD_STR_LIT("run"), .executor = test_executor_of_threads(3)}, alloc);
    ASSERT_TRUE(e != NULL);
    const bool ran = md_vial_eval_run(e);
    const md_vial_diagnostic_t* rd;
    const size_t nr = md_vial_eval_diagnostics(e, &rd);
    for (size_t i = 0; i < nr; ++i) printf("%.5s %.*s\n", rd[i].code, (int)rd[i].message.len, rd[i].message.ptr);
    ASSERT_TRUE(ran);
    EXPECT_EQ(0u, (unsigned)nr);
    const md_attributes_t* at = md_vial_eval_attributes(e);
    size_t n = 0;
    // Paired column by column as mcs gives them, which is the residues' own order: the same volume; and with the
    // columns the other way round, which pairs the same particles
    EXPECT_LT(vial_ir__max_diff(at, "script/s", "script/so", &n), 1e-6);
    EXPECT_EQ(40u * 40u * 40u, (unsigned)n);
    EXPECT_LT(vial_ir__max_diff(at, "script/s", "script/rv", NULL), 2e-4);    // the fits sum in another order
    // The same structures twice through mcs over arrays: the same density per structure, twice as much added
    EXPECT_LT(vial_ir__max_diff(at, "script/s", "script/s2", NULL), 1e-4);
    EXPECT_EQ(2 * vial_ir__attr1(at, "script/s/samples", 0), vial_ir__attr1(at, "script/s2/samples", 0));
    // Four alike residues through one mcs: the same as the four fitted as they are, on their heavy atoms
    EXPECT_LT(vial_ir__max_diff(at, "script/c1", "script/ct", NULL), 1e10);     // both exist; c1 has four rows, ct five
    // A structure not alike, the last residue with its OXT, joins by the five heavy atoms it has in common
    const md_attribute_t* act = md_attributes_find(at, MD_STR_LIT("script/ct"));
    ASSERT_TRUE(act != NULL);
    double tot = 0;
    for (size_t i = 0; i < 40 * 40 * 40; ++i) tot += ((const float*)act->data)[i];
    EXPECT_GT(tot, 0.0);
    // Drawn: every structure with its pose
    const uint32_t expr = md_vial_expr_at(prog, (uint32_t)(strstr(src, "sdf(mcs(residue[1:5], residue[14:15])") - src));
    md_vial_drawing_t dr;
    ASSERT_TRUE(md_vial_draw(&dr, &(md_vial_draw_desc_t){.prog = prog, .sys = sys, .run = MD_STR_LIT("run"), .frame = 3, .element = -1}, expr, alloc));
    EXPECT_EQ(5u, (unsigned)dr.num_poses);
    md_vial_eval_free(e);
    md_vial_program_free(prog);

    // What the check says
    const md_system_t* T = vial_test_system('T');
    ASSERT_TRUE(T != NULL);
    struct { const md_system_t* sys; const char* src; const char* code; } cases[] = {
        {T, "m = mcs(particle(resname='GTP'), particle(resname='VLB'))\n", "E0321"},
        {T, "v = sdf(mcs(residue('PHE')[0:2], residue('TYR')[0:2]), particle(name='O'), 8.0)\n", NULL},
        {T, "a = mcs(residue('PHE')[0:2])\nb = mcs(residue('TYR')[0:2])\nv = sdf(concatenate((a, b)), particle(name='O'), 8.0)\n", "E0309"},
        {T, "a = mcs(residue('PHE')[0:2])\nb = mcs(residue('TYR')[0:2])\nv = concatenate((a, b), axis=1)\n", "E0316"},
        {T, "v = sdf((residue('PHE'), residue('TYR')), particle(name='O'), 8.0)\n", "E0313"},
        {T, "v = sdf(residue('PHE'), particle(name='O'), 8.0, poses=None)\n", "E0202"},
    };
    for (size_t i = 0; i < MD_COUNTOF(cases); ++i) {
        md_vial_program_t* p = md_vial_compile(md_str_from_cstr(cases[i].src), &(md_vial_compile_desc_t){.sys = cases[i].sys}, alloc);
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

// ### ENERGY TERMS DRAW THEIR GROUPS ###

#include <md_edr.h>

// inside-md-pullout: its energy file onto its structure, with groups as an index file of the run would give them; the
// energy file names them the energy groups, CA Protein AIN SOL rest. attr() of a term draws its groups, each a group of
// the drawing of its own; of a pattern every group once, or of one member that member's.
UTEST(vial_ir, attr_draws_groups) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(64));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t st = {.alloc = alloc};
    ASSERT_TRUE(md_gro_system_init_from_file(&sys, &st, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/inside-md-pullout.gro"), NULL));
    ASSERT_GE(sys.atom.count, 1000u);
    md_system_group_add(&sys, MD_STR_LIT("CA"),      md_bitset_range(0, 10));
    md_system_group_add(&sys, MD_STR_LIT("Protein"), md_bitset_range(10, 300));
    md_system_group_add(&sys, MD_STR_LIT("AIN"),     md_bitset_range(300, 400));
    md_system_group_add(&sys, MD_STR_LIT("SOL"),     md_bitset_range(400, 900));
    ASSERT_TRUE(md_edr_system_supplement_from_file(&sys, MD_STR_LIT(MD_UNITTEST_DATA_DIR "/inside-md-pullout.edr"), &(md_load_desc_t){ .run = MD_STR_LIT("run/md") }));
    size_t er_n;
    md_system_group_role_order(&sys.group, MD_GROUP_ROLE_ENERGY, &er_n);
    ASSERT_EQ(5u, er_n);

    const char* src = "p = attr('edr/coul_sr_ca_protein')\nall = attr('edr/coul_sr_*')\nself = attr('edr/lj_sr_sol_sol')\nu = attr('edr/potential')\n";
    const md_vial_compile_desc_t cd = {.sys = &sys, .run = MD_STR_LIT("run/md")};
    md_vial_program_t* prog = md_vial_compile(md_str_from_cstr(src), &cd, alloc);
    ASSERT_TRUE(prog != NULL);
    #define AT(text, plus) ((uint32_t)(strstr(src, text) - src) + (plus))
    md_vial_draw_desc_t desc = {.prog = prog, .sys = &sys, .run = MD_STR_LIT("run/md"), .element = -1};
    md_vial_drawing_t dr;

    // A pair: CA, then Protein, each a group of its own
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("attr('edr/coul_sr_ca", 2)), alloc));
    ASSERT_EQ(2u, dr.num_groups);
    EXPECT_EQ(10u, dr.groups[0].end - dr.groups[0].beg);
    EXPECT_EQ(0u, dr.atoms[dr.groups[0].beg]);
    EXPECT_EQ(290u, dr.groups[1].end - dr.groups[1].beg);
    EXPECT_EQ(10u, dr.atoms[dr.groups[1].beg]);
    EXPECT_EQ(0u, dr.num_prims);
    // The name bound to it draws as it does
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("p = ", 0)), alloc));
    EXPECT_EQ(2u, dr.num_groups);

    // Every Coul-SR pair: the five energy groups, each once; rest is the atoms of none of the others
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("attr('edr/coul_sr_*", 2)), alloc));
    EXPECT_EQ(15u, dr.num_elements);
    ASSERT_EQ(5u, dr.num_groups);
    size_t total = 0;
    for (size_t i = 0; i < dr.num_groups; ++i) total += dr.groups[i].end - dr.groups[i].beg;
    EXPECT_EQ(sys.atom.count, total);
    // One member alone: in path order the first is coul_sr_ain_ain, AIN with itself
    md_vial_draw_desc_t one = desc;
    one.element = 0;
    ASSERT_TRUE(md_vial_draw(&dr, &one, md_vial_expr_at(prog, AT("attr('edr/coul_sr_*", 2)), alloc));
    ASSERT_EQ(1u, dr.num_groups);
    EXPECT_EQ(100u, dr.groups[0].end - dr.groups[0].beg);
    EXPECT_EQ(300u, dr.atoms[dr.groups[0].beg]);

    // A group with itself: the group once
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("attr('edr/lj_sr_sol", 2)), alloc));
    ASSERT_EQ(1u, dr.num_groups);
    EXPECT_EQ(500u, dr.groups[0].end - dr.groups[0].beg);
    // A term of no groups draws nothing
    ASSERT_TRUE(md_vial_draw(&dr, &desc, md_vial_expr_at(prog, AT("attr('edr/potential", 2)), alloc));
    EXPECT_EQ(0u, dr.num_groups);
    #undef AT
    md_vial_program_free(prog);
    md_arena_allocator_destroy(alloc);
}


// What hover shows of a set or an array of sets: described, its shape and what each set is, not "set of n" set by set
UTEST(vial_ir, sets_described) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(16));
    const md_system_t* sys = vial_test_system('T');
    ASSERT_TRUE(sys != NULL);
    const char* cases[][2] = {
        {"particle",                             "6927 particles"},
        {"particle[0]",                          "particle N, MET 1"},
        {"residue[2]",                           "residue GLU 3, 9 particles"},
        {"chain[0]",                             "chain A, 3441 particles"},
        {"residue",                              "880 residues, 1 to 59 particles each"},
        {"residue('ALA')",                       "60 residues ALA, 5 particles each"},
        {"chain",                                "3 chains, 59 to 3441 particles each"},
        {"molecule",                             "17 molecules, 1 to 3393 particles each"},
        {"residue('PHE').ring[0]",               "42 rings, 6 particles each"},
        {"residue('ALA')['CA']",                 "60 × CA"},
        {"residue('ALA')[['N', 'CA', 'C']]",     "60 × (N, CA, C)"},
        {"residue('ALA').particle[[1, 0, 2]]",   "60 × (CA, N, C)"},
        {"particle[[0, 1, 2]]",                  "(N, CA, C)"},
        {"bond(residue[2])",                     "8 × 2: (N, CA), (CA, C), (CA, CB), …"},
        {"residue('GLY')['CB']",                 "63 empty sets"},
        {"particle[particle.z < 0.0]",           "empty set"},
        {"particle(residue('ALA')[['N', 'CA']])", "120 particles"},
    };
    for (size_t i = 0; i < MD_COUNTOF(cases); ++i) {
        char src[256];
        snprintf(src, sizeof(src), "x = %s\n", cases[i][0]);
        md_vial_program_t* p = md_vial_compile(md_str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys}, alloc);
        char text[256] = "";
        const md_vial_draw_desc_t dd = {.prog = p, .sys = sys, .state = vial_test_state('T'), .element = -1};
        md_vial_value_text(text, sizeof(text), &dd, md_vial_expr_at(p, 0));
        EXPECT_STREQ(cases[i][1], text);
        md_vial_program_free(p);
    }
    // Other values are printed as they were
    {
        md_vial_program_t* p = md_vial_compile(MD_STR_LIT("x = count(residue('ALA'))\n"), &(md_vial_compile_desc_t){.sys = sys}, alloc);
        char text[64] = "";
        const md_vial_draw_desc_t dd = {.prog = p, .sys = sys, .state = vial_test_state('T'), .element = -1};
        md_vial_value_text(text, sizeof(text), &dd, md_vial_expr_at(p, 0));
        EXPECT_STREQ("[5, 5, 5, 5, 5, 5, …]", text);
        md_vial_program_free(p);
    }
    md_arena_allocator_destroy(alloc);
}



