#include "utest.h"
#include "test_executor.h"

#include <core/md_common.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

#include <md_vial.h>
#include <md_script.h>
#include <md_attributes.h>
#include "../src/vial/vial_internal.h"
#include "../src/vial/ir.h"
#include "../src/vial/interp.h"

#include <core/md_os.h>

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
    {"p = distance(1, 2);",                                  "p = distance(particle[0], particle[1])", 1e-4},
    {"p = distance(1, 30);",                                 "p = distance(particle[0], particle[29])", 1e-4},
    {"p = angle(1, 2, 3);",                                  "p = angle(particle[0], particle[1], particle[2])", 1e-3},
    {"p = dihedral(5, 7, 9, 15);",                           "p = dihedral(particle[4], particle[6], particle[8], particle[14])", 1e-3},
    {"p = distance(residue(1), residue(5));",                "p = distance(residue[0].com, residue[4].com)", 1e-4},
    {"p = distance_min(residue(1), residue(5));",            "p = distance(residue[0], residue[4], between='closest')", 1e-4},
    {"p = distance_max(residue(1), residue(5));",            "p = distance(residue[0], residue[4], between='farthest')", 1e-4},
    {"p = angle(residue(1), residue(2), residue(3));",       "p = angle(residue[0].com, residue[1].com, residue[2].com)", 1e-3},
    {"p = count(within(5.0, residue(3)));",                  "p = len(within(5.0, of=residue[2]) and not residue[2])", 0.0},
    {"p = count(within(8.0, atom(1:10)));",                  "p = len(within(8.0, of=particle[0:10]) and not particle[0:10])", 0.0},
    {"p = angle(1, 50, 120);",                               "p = angle(particle[0], particle[49], particle[119])", 1e-3},
    {"p = dihedral(residue(1), residue(2), residue(3), residue(4));",
                                                             "p = dihedral(residue[0].com, residue[1].com, residue[2].com, residue[3].com)", 1e-3},
    {"p = distance_min(atom(1:20), atom(100:130));",         "p = distance(particle[0:20], particle[99:130], between='closest')", 1e-4},
    {"p = distance_max(atom(1:20), atom(100:130));",         "p = distance(particle[0:20], particle[99:130], between='farthest')", 1e-4},
    {"p = coord_x(atom(10));",                               "p = particle[9].x", 1e-5},
    {"p = coord_z(com(residue(7)));",                        "p = residue[6].com.z", 1e-4},
    {"p = rmsd(atom(1:40));",                                "x = particle[0:40]\np = rmsd(x, at(x.xyz, frame=0), weights='mass')", 1e-3},
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

// v2's output p at every frame, through the public API: published as script/p
static bool v2_values(const char* src, const md_system_t* sys, uint32_t nf, double* out, md_allocator_i* alloc) {
    const md_vial_compile_desc_t desc = {.sys = sys, .run = STR_LIT("run")};
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &desc, alloc);
    const md_vial_diagnostic_t* d;
    const size_t nd = md_vial_diagnostics(prog, &d);
    for (size_t i = 0; i < nd; ++i) if (d[i].code[0] == 'E') { printf("v2 does not compile: %s: %.*s\n", src, (int)d[i].message.len, d[i].message.ptr); return false; }
    md_vial_eval_t* eval = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    if (!eval) { printf("v2 does not evaluate: %s\n", src); md_vial_program_free(prog); return false; }
    bool ok = md_vial_eval_run(eval);
    const md_vial_diagnostic_t* rd;
    const size_t nr = md_vial_eval_diagnostics(eval, &rd);
    for (size_t i = 0; i < nr; ++i) printf("v2 fails: %s: %.5s %.*s\n", src, rd[i].code, (int)rd[i].message.len, rd[i].message.ptr);
    const md_attribute_t* attr = md_attributes_find(md_vial_eval_attributes(eval), STR_LIT("script/p"));
    ok = ok && attr && (attr->flags & MD_ATTRIBUTE_FLAG_TEMPORAL) && attr->format.shape[0] == nf;
    for (uint32_t f = 0; ok && f < nf; ++f) {
        double v = 0.0;
        ok = md_attribute_extract_f64(&v, 1, attr, md_attribute_slice_1(f), attr->unit) == 1;
        out[f] = v;
    }
    md_vial_eval_free(eval);
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
        const md_tick_t t0 = md_tick_now();
        if (!v1_values(pr->v1, sys, nf, a, alloc)) { EXPECT_TRUE(false); continue; }
        const md_tick_t t1 = md_tick_now();
        if (!v2_values(pr->v2, sys, nf, b, alloc)) { EXPECT_TRUE(false); continue; }
        const md_tick_t t2 = md_tick_now();
        if (getenv("VIAL_COMPARE_TIME")) {
            printf("vial compare r1: v1 %7.3f ms, v2 %7.3f ms: %s\n", md_tick_to_milliseconds(t1 - t0), md_tick_to_milliseconds(t2 - t1), pr->v2);
        }
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

// sdf: v1's volume and VIAL's are of the same particles in the same frame. They are not the same numbers: v1 counts
// each particle into one voxel and sums over the structures, VIAL spreads it over eight and gives a density per
// structure. So what they hold is compared: the particles per structure in the box, and where they are.
UTEST(vial_compare, sdf_agrees_with_v1) {
    const md_system_t* sys = vial_test_system('P');
    const uint32_t nf = (uint32_t)vial_test_run('P');
    ASSERT_TRUE(sys != NULL);
    ASSERT_GT(nf, 1u);
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(64));
    enum { D = 128 };
    const double ext = 6.0, step = 2.0 * ext / D;
    // v1
    md_script_ir_t* ir = md_script_ir_create(alloc);
    ASSERT_TRUE(md_script_ir_compile_from_source(ir, STR_LIT("p = sdf(residue(2:5), name(\"O\"), 6.0);"), sys, NULL));
    md_script_eval_t* ev1 = md_script_eval_create(nf, ir, alloc);
    ASSERT_TRUE(ev1 && md_script_eval_frame_range(ev1, ir, sys, STR_LIT("run"), 0, nf));
    const md_attribute_t* a1 = md_attributes_find(md_script_eval_attributes(ev1), STR_LIT("script/p"));
    ASSERT_TRUE(a1 && a1->data);
    // VIAL, at v1's 128 a side
    const char* src = "p = sdf(residue[1:5], particle(name='O'), 6.0, samples=128 / 12)\n";
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &(md_vial_compile_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    ASSERT_EQ(0u, (unsigned)md_vial_diagnostics(prog, NULL));
    md_vial_eval_t* ev2 = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run"), .executor = test_executor_of_threads(4)}, alloc);
    ASSERT_TRUE(ev2 && md_vial_eval_run(ev2));
    const md_attribute_t* a2 = md_attributes_find(md_vial_eval_attributes(ev2), STR_LIT("script/p"));
    ASSERT_TRUE(a2 && a2->data);
    ASSERT_EQ((uint32_t)D, a2->format.shape[0]);
    const float* v1 = (const float*)a1->data;
    const float* v2 = (const float*)a2->data;
    double n1 = 0, n2 = 0, m1[3] = {0}, m2[3] = {0};
    for (uint32_t x = 0; x < D; ++x) for (uint32_t y = 0; y < D; ++y) for (uint32_t z = 0; z < D; ++z) {
        const double c[3] = {(x + 0.5) * step - ext, (y + 0.5) * step - ext, (z + 0.5) * step - ext};
        const double w1 = v1[((size_t)z * D + y) * D + x] / 4.0;                       // x fastest, summed over 4
        const double w2 = v2[((size_t)x * D + y) * D + z] * step * step * step;         // z fastest, per Å³
        n1 += w1; n2 += w2;
        for (int a = 0; a < 3; ++a) { m1[a] += w1 * c[a]; m2[a] += w2 * c[a]; }
    }
    printf("sdf: v1 %.4f, VIAL %.4f O per structure; centres (%.3f %.3f %.3f) and (%.3f %.3f %.3f)\n", n1, n2,
           m1[0] / n1, m1[1] / n1, m1[2] / n1, m2[0] / n2, m2[1] / n2, m2[2] / n2);
    EXPECT_GT(n1, 1.0);
    EXPECT_NEAR(n1, n2, 0.02 * n1);
    for (int a = 0; a < 3; ++a) EXPECT_NEAR(m1[a] / n1, m2[a] / n2, 0.05);
    md_vial_eval_free(ev2);
    md_vial_program_free(prog);
    md_script_eval_free(ev1);
    md_script_ir_free(ir);
    md_arena_allocator_destroy(alloc);
}
