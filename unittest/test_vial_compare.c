#include "utest.h"

#include <core/md_common.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

#include <md_vial.h>
#include <md_script.h>
#include <md_attributes.h>
#include "../src/vial/vial_internal.h"
#include "../src/vial/ir.h"
#include "../src/vial/interp.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "vial_test_systems.h"

// The gate of phase 3: v2 gives v1's numbers on the same data. Each pair is a v1 script and a v2 script that should
// compute the same property over every frame of P's run; they agree within a float tolerance, since v1's formulas
// are not v2's to the last bit.

typedef struct vial_pair_t {
    const char* v1;         // the v1 script; its property is `p`
    const char* v2;         // the v2 script; its output is `p`
    double      tol;        // absolute
} vial_pair_t;

static const vial_pair_t vial_pairs[] = {
    {"p = distance(1, 2);",                                  "p = distance(particles[0], particles[1])", 1e-4},
    {"p = distance(1, 30);",                                 "p = distance(particles[0], particles[29])", 1e-4},
    {"p = angle(1, 2, 3);",                                  "p = angle(particles[0], particles[1], particles[2])", 1e-3},
    {"p = dihedral(5, 7, 9, 15);",                           "p = dihedral(particles[4], particles[6], particles[8], particles[14])", 1e-3},
    {"p = distance(residue(1), residue(5));",                "p = distance(residues[0].com, residues[4].com)", 1e-4},
    {"p = distance_min(residue(1), residue(5));",            "p = distance(residues[0], residues[4], between='closest')", 1e-4},
    {"p = distance_max(residue(1), residue(5));",            "p = distance(residues[0], residues[4], between='farthest')", 1e-4},
    {"p = angle(residue(1), residue(2), residue(3));",       "p = angle(residues[0].com, residues[1].com, residues[2].com)", 1e-3},
    {"p = count(within(5.0, residue(3)));",                  "p = len(within(5.0, of=residues[2]) and not residues[2])", 0.0},
};

// v1's property p at every frame
static bool v1_values(const char* src, const md_system_t* sys, uint32_t nf, double* out, md_allocator_i* alloc) {
    md_script_ir_t* ir = md_script_ir_create(alloc);
    if (!md_script_ir_compile_from_source(ir, str_from_cstr(src), sys, NULL)) { printf("v1 does not compile: %s\n", src); return false; }
    md_script_eval_t* eval = md_script_eval_create(nf, ir, alloc);
    if (!eval || !md_script_eval_frame_range(eval, ir, sys, STR_LIT("run"), 0, nf)) { printf("v1 does not evaluate: %s\n", src); return false; }
    const md_attribute_t* attr = md_attributes_find(md_script_eval_attributes(eval), STR_LIT("script/p"));
    if (!attr) return false;
    // v2 gives angles in degrees, v1 in radians
    md_unit_t unit = attr->unit;
    double factor;
    if (!md_unit_is_none(unit) && md_unit_conversion_factor(&factor, unit, md_unit_degree())) unit = md_unit_degree();
    for (uint32_t f = 0; f < nf; ++f) {
        float v = 0.0f;
        if (md_attribute_extract_f32(&v, 1, attr, md_attribute_slice_1(f), unit) != 1) return false;
        out[f] = v;
    }
    md_script_eval_free(eval);
    md_script_ir_free(ir);
    return true;
}

// v2's output p at every frame
static bool v2_values(const char* src, const md_system_t* sys, uint32_t nf, double* out, md_allocator_i* alloc) {
    const md_vial_compile_desc_t desc = {.sys = sys};
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &desc, alloc);
    const md_vial_diagnostic_t* d;
    const size_t nd = md_vial_diagnostics(prog, &d);
    for (size_t i = 0; i < nd; ++i) if (d[i].code[0] == 'E') { printf("v2 does not compile: %s: %.*s\n", src, (int)d[i].message.len, d[i].message.ptr); return false; }
    const vial_ir_t* ir = md_vial__ir(prog);
    vial_ssa_t p = 0;
    for (size_t i = 0; i < md_array_size(ir->outputs); ++i) if (str_eq(ir->outputs[i].name, STR_LIT("p"))) p = ir->outputs[i].value;
    if (!p) return false;
    vial_sys_t vs;
    vial_sys_init(&vs, sys, alloc);
    vial_interp_t in;
    vial_interp_init(&in, ir, &vs, sys, md_get_heap_allocator());
    md_system_state_t state = {.alloc = alloc};
    md_system_state_init(&state, sys->atom.count);
    const str_t paths[] = {STR_LIT("atom/position"), STR_LIT("unitcell")};
    md_system_extract_t* ex = md_system_extract_begin(sys, STR_LIT("run"), paths, 2, alloc);
    bool ok = ex != NULL;
    for (uint32_t f = 0; ok && f < nf; ++f) {
        ok = md_system_extract_frame(ex, f, &state);
        const vial_frame_t fr = {.xyz = (const float*)state.xyz, .num_atoms = state.num_atoms, .cell = state.unitcell, .index = f, .time = f};
        ok = ok && vial_interp_frame(&in, &fr);
        const vial_val_t* v = ok ? vial_interp_value(&in, p) : NULL;
        if (!v || !v->data || vial_val_count(v) != 1) { ok = false; break; }
        out[f] = v->dtype == VIAL_DT_FLOAT ? *(const float*)v->data : v->dtype == VIAL_DT_INT ? *(const int32_t*)v->data : *(const double*)v->data;
    }
    if (!ok && in.failed) printf("v2 fails: %s: %s\n", src, in.error.message);
    if (ex) md_system_extract_end(ex);
    vial_interp_free(&in);
    md_vial_program_free(prog);
    return ok;
}

UTEST(vial_compare, v1_and_v2_agree) {
    const md_system_t* sys = vial_test_system('P');
    const uint32_t nf = (uint32_t)vial_test_run('P');
    ASSERT_TRUE(sys != NULL);
    ASSERT_GT(nf, 1u);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
    double* a = md_alloc(alloc, nf * sizeof(double));
    double* b = md_alloc(alloc, nf * sizeof(double));
    int agree = 0;
    for (size_t i = 0; i < ARRAY_SIZE(vial_pairs); ++i) {
        const vial_pair_t* pr = &vial_pairs[i];
        if (!v1_values(pr->v1, sys, nf, a, alloc)) { EXPECT_TRUE(false); continue; }
        if (!v2_values(pr->v2, sys, nf, b, alloc)) { EXPECT_TRUE(false); continue; }
        double worst = 0.0;
        uint32_t at = 0;
        for (uint32_t f = 0; f < nf; ++f) if (fabs(a[f] - b[f]) > worst) { worst = fabs(a[f] - b[f]); at = f; }
        if (worst > pr->tol) {
            printf("v1 and v2 differ by %g at frame %u: v1 %.6f, v2 %.6f\n  %s\n  %s\n", worst, at, a[at], b[at], pr->v1, pr->v2);
            EXPECT_TRUE(false);
        } else {
            agree += 1;
        }
    }
    printf("vial compare: %d of %d pairs agree over %u frames\n", agree, (int)ARRAY_SIZE(vial_pairs), nf);
    md_arena_allocator_destroy(alloc);
}
