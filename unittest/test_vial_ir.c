#include "utest.h"

#include <core/md_common.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

#include <md_vial.h>
#include "../src/vial/vial_internal.h"
#include "../src/vial/ir.h"

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
    const vial_ir_t* ir = md_vial__ir(r->prog);
    if (md_array_size(ir->failures)) {
        printf("not lowered: %s\n", ir->failures[0].what);
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
    const vial_ir_t* ir = md_vial__ir(r->prog);
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
    const vial_ir_t* ir = md_vial__ir(r->prog);
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
