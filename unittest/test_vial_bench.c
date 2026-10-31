#include "utest.h"

// A benchmark of v2 against v1, run on demand: VIAL_BENCH=1 md_unittest --filter=vial_bench.*
//
//   VIAL_BENCH_TOP    a .gro or .pdb file (default test_data/catalyst.gro)
//   VIAL_BENCH_TRAJ   an .xtc, .trr, .dcd or multi-model .pdb file (default test_data/catalyst.xtc)
//   VIAL_BENCH_REPS   repetitions; each figure is the fastest (default 3)
//
// Per pair of scripts: v1's and v2's compile time and evaluation time over every frame. Then what the deferred
// cache across compiles would save: compiling 100 selection lines, and at()'s pre-step.

#include <core/md_common.h>
#include <core/md_arena_allocator.h>
#include <core/md_os.h>
#include <core/md_str.h>

#include <md_vial.h>
#include <md_script.h>
#include <md_system.h>
#include <md_attributes.h>
#include <md_gro.h>
#include <md_pdb.h>
#include <md_xtc.h>
#include <md_trr.h>
#include <md_dcd.h>
#include <md_util.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define VIAL_BENCH_VERSION "vial bench r6"

typedef struct vial_bench_pair_t {
    const char* v1;
    const char* v2;
} vial_bench_pair_t;

static const vial_bench_pair_t vial_bench_pairs[] = {
    {"p = distance(1, 2);",                                   "p = distance(particles[0], particles[1])"},
    {"p = distance(com(atom(1:50)), com(atom(51:100)));",     "p = distance(particles[0:50].com, particles[50:100].com)"},
    {"p = angle(1, 2, 3);",                                   "p = angle(particles[0], particles[1], particles[2])"},
    {"p = distance_min(atom(1:50), atom(51:100));",           "p = distance(particles[0:50], particles[50:100], between='closest')"},
    {"p = count(within(5.0, atom(1:50)));",                   "p = len(within(5.0, of=particles[0:50]) and not particles[0:50])"},
    {"p = distance(residue(1), residue(2));",                 "p = distance(residues[0].com, residues[1].com)"},
    // The rest of the comparison corpus (test_vial_compare.c), on what the catalyst has
    {"p = distance(1, 300);",                                 "p = distance(particles[0], particles[299])"},
    {"p = dihedral(5, 7, 9, 15);",                            "p = dihedral(particles[4], particles[6], particles[8], particles[14])"},
    {"p = angle(1, 50, 120);",                                "p = angle(particles[0], particles[49], particles[119])"},
    {"p = distance_max(atom(1:50), atom(51:100));",           "p = distance(particles[0:50], particles[50:100], between='farthest')"},
    {"p = distance_min(atom(1:20), atom(100:130));",          "p = distance(particles[0:20], particles[99:130], between='closest')"},
    {"p = count(within(8.0, atom(1:10)));",                   "p = len(within(8.0, of=particles[0:10]) and not particles[0:10])"},
    {"p = coord_x(atom(10));",                                "p = particles[9].x"},
    {"p = coord_z(com(residue(2)));",                         "p = residues[1].com.z"},
    {"p = angle(residue(1), residue(2), residue(3));",        "p = angle(residues[0].com, residues[1].com, residues[2].com)"},
    {"p = rmsd(atom(1:40));",                                 "x = particles[0:40]\np = rmsd(x, at(x.xyz, frame=0), weights='mass')"},
    {"p = distance_min(atom(1:400), atom(401:1336));",        "p = distance(particles[0:400], particles[400:1336], between='closest')"},
    // The same reductions written as a matrix of every pair, which the compiler computes by a search
    {"p = distance_min(atom(1:400), atom(401:1336));",        "p = min(distance(particles[0:400].xyz[:, None], particles[400:1336].xyz))"},
    {"p = distance_max(atom(1:50), atom(51:100));",           "p = max(distance(particles[0:50].xyz[:, None], particles[50:100].xyz))"},
};

static bool str_ends(str_t s, const char* e) {
    const size_t n = strlen(e);
    return s.len >= n && strncmp(s.ptr + s.len - n, e, n) == 0;
}

static bool vial_bench_load(md_system_t* sys, md_system_state_t* state, str_t top, str_t traj) {
    bool ok = str_ends(top, ".gro") ? md_gro_system_init_from_file(sys, state, top) :
              str_ends(top, ".pdb") ? md_pdb_system_init_from_file(sys, state, top, MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE) : false;
    if (!ok) return false;
    md_util_system_infer(sys, state, MD_UTIL_INFER_ALL);
    const uint32_t flags = MD_RUN_FLAG_DISABLE_CACHE_WRITE;
    return str_ends(traj, ".xtc") ? md_xtc_system_publish_run(sys, traj, STR_LIT("run"), flags) :
           str_ends(traj, ".trr") ? md_trr_system_publish_run(sys, traj, STR_LIT("run"), flags) :
           str_ends(traj, ".dcd") ? md_dcd_system_publish_run(sys, traj, STR_LIT("run"), flags) :
           str_ends(traj, ".pdb") ? md_pdb_system_publish_run(sys, traj, STR_LIT("run"), flags) : false;
}

static double ms_since(md_tick_t t0) { return md_tick_to_milliseconds(md_tick_now() - t0); }

static bool bench_v1(const char* src, const md_system_t* sys, uint32_t nf, md_allocator_i* alloc, double* compile_ms, double* eval_ms) {
    md_tick_t t0 = md_tick_now();
    md_script_ir_t* ir = md_script_ir_create(alloc);
    if (!md_script_ir_compile_from_source(ir, str_from_cstr(src), sys, NULL)) return false;
    *compile_ms = ms_since(t0);
    t0 = md_tick_now();
    md_script_eval_t* eval = md_script_eval_create(nf, ir, alloc);
    const bool ok = eval && md_script_eval_frame_range(eval, ir, sys, STR_LIT("run"), 0, nf);
    *eval_ms = ms_since(t0);
    if (eval) md_script_eval_free(eval);
    md_script_ir_free(ir);
    return ok;
}

static bool bench_v2(const char* src, const md_system_t* sys, md_allocator_i* alloc, double* compile_ms, double* eval_ms) {
    md_tick_t t0 = md_tick_now();
    const md_vial_compile_desc_t desc = {.sys = sys, .run = STR_LIT("run")};
    md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &desc, alloc);
    *compile_ms = ms_since(t0);
    t0 = md_tick_now();
    md_vial_eval_t* eval = md_vial_eval_create(prog, &(md_vial_eval_desc_t){.sys = sys, .run = STR_LIT("run")}, alloc);
    const bool ok = eval && md_vial_eval_run(eval);
    *eval_ms = ms_since(t0);
    if (eval) md_vial_eval_free(eval);
    md_vial_program_free(prog);
    return ok;
}

UTEST(vial_bench, against_v1) {
    if (!getenv("VIAL_BENCH")) return;
    const char* top_env = getenv("VIAL_BENCH_TOP");
    const char* traj_env = getenv("VIAL_BENCH_TRAJ");
    const str_t top = top_env ? str_from_cstr(top_env) : STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.gro");
    const str_t traj = traj_env ? str_from_cstr(traj_env) : STR_LIT(MD_UNITTEST_DATA_DIR "/catalyst.xtc");
    const int reps = getenv("VIAL_BENCH_REPS") ? atoi(getenv("VIAL_BENCH_REPS")) : 3;

    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(16));
    md_system_t sys = {.alloc = alloc};
    md_system_state_t state = {.alloc = alloc};
    ASSERT_TRUE(vial_bench_load(&sys, &state, top, traj));
    const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/time"));
    ASSERT_TRUE(time != NULL);
    const uint32_t nf = time->format.shape[0];
    printf("%s: %.*s, %zu particles, %zu residues; %.*s, %u frames; fastest of %d\n", VIAL_BENCH_VERSION, (int)top.len, top.ptr,
        sys.atom.count, md_system_residue_count(&sys), (int)traj.len, traj.ptr, nf, reps);
    printf("%s: %2s %-72s %10s %10s %10s %10s\n", VIAL_BENCH_VERSION, "", "v2 script", "v1 compile", "v2 compile", "v1 eval", "v2 eval");

    // VIAL_BENCH_ONLY=i: pair i alone, for a profiler; nothing else is run
    const int only = getenv("VIAL_BENCH_ONLY") ? atoi(getenv("VIAL_BENCH_ONLY")) : -1;
    double sum_v1 = 0.0, sum_v2 = 0.0;
    for (size_t i = 0; i < ARRAY_SIZE(vial_bench_pairs); ++i) {
        if (only >= 0 && (size_t)only != i) continue;
        double best[4] = {1e30, 1e30, 1e30, 1e30};
        bool ok = true;
        for (int r = 0; r < reps && ok; ++r) {
            md_allocator_i* a = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
            double c1, e1, c2, e2;
            // Every other rep runs v2 first, so neither gains from running second, in a warm cache
            ok = r % 2 ? bench_v2(vial_bench_pairs[i].v2, &sys, a, &c2, &e2) && bench_v1(vial_bench_pairs[i].v1, &sys, nf, a, &c1, &e1)
                       : bench_v1(vial_bench_pairs[i].v1, &sys, nf, a, &c1, &e1) && bench_v2(vial_bench_pairs[i].v2, &sys, a, &c2, &e2);
            if (ok) {
                if (c1 < best[0]) best[0] = c1;
                if (c2 < best[1]) best[1] = c2;
                if (e1 < best[2]) best[2] = e1;
                if (e2 < best[3]) best[3] = e2;
            }
            md_arena_allocator_destroy(a);
        }
        // The script on one line, its lines joined by "; "
        char shown[256];
        size_t w = 0;
        for (const char* c = vial_bench_pairs[i].v2; *c && w + 3 < sizeof(shown); ++c) {
            if (*c == '\n') { shown[w++] = ';'; shown[w++] = ' '; }
            else shown[w++] = *c;
        }
        shown[w] = '\0';
        if (!ok) { printf("%s: %2zu %-72s failed\n", VIAL_BENCH_VERSION, i, shown); continue; }
        printf("%s: %2zu %-72s %8.3f ms %8.3f ms %8.2f ms %8.2f ms  %5.2fx\n", VIAL_BENCH_VERSION, i, shown, best[0], best[1], best[2], best[3],
            best[3] / best[2]);
        sum_v1 += best[2];
        sum_v2 += best[3];
    }
    printf("%s: %2s %-72s %10s %10s %8.2f ms %8.2f ms  %5.2fx\n", VIAL_BENCH_VERSION, "", "every pair", "", "", sum_v1, sum_v2, sum_v2 / sum_v1);
    if (only >= 0) {
        md_arena_allocator_destroy(alloc);
        return;
    }

    // What a cache across compiles would save: 100 lines of selections, compiled as on every keystroke
    {
        char* src = NULL;
        md_allocator_i* a = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
        const uint32_t nres = (uint32_t)md_system_residue_count(&sys);
        for (uint32_t k = 0; k < 100; ++k) {
            char line[160];
            const uint32_t lo = nres ? (k * 7) % nres : 0;
            const int len = k % 2 ? snprintf(line, sizeof(line), "s%u = residues[%u:%u]\n", k, lo, lo + 5)
                                  : snprintf(line, sizeof(line), "s%u = particles(element='C') and not particles[%u:%u]\n", k, lo, lo + 20);
            md_array_push_array(src, line, (size_t)len, a);
        }
        md_array_push(src, '\0', a);
        double best = 1e30, empty = 1e30;
        for (int r = 0; r < reps; ++r) {
            double c, e;
            md_tick_t t0 = md_tick_now();
            md_vial_program_t* prog = md_vial_compile(str_from_cstr(src), &(md_vial_compile_desc_t){.sys = &sys}, a);
            c = ms_since(t0);
            md_vial_program_free(prog);
            if (c < best) best = c;
            t0 = md_tick_now();
            prog = md_vial_compile(STR_LIT(""), &(md_vial_compile_desc_t){.sys = &sys}, a);
            e = ms_since(t0);
            md_vial_program_free(prog);
            if (e < empty) empty = e;
        }
        printf("%s: compiling 100 selection lines: %.3f ms; an empty script: %.3f ms\n", VIAL_BENCH_VERSION, best, empty);
        md_arena_allocator_destroy(a);
    }

    // at()'s pre-step, computed at every evaluation since it is not cached
    {
        double with_at = 1e30, without = 1e30;
        for (int r = 0; r < reps; ++r) {
            md_allocator_i* a = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(4));
            double c, e;
            if (bench_v2("ref = at(particles.xyz, frame=0)\np = rmsd(particles, ref)", &sys, a, &c, &e) && e < with_at) with_at = e;
            if (bench_v2("p = rmsd(particles, particles.xyz)", &sys, a, &c, &e) && e < without) without = e;
            md_arena_allocator_destroy(a);
        }
        printf("%s: rmsd of every particle against at(frame=0): %.2f ms; against itself, no at(): %.2f ms\n", VIAL_BENCH_VERSION, with_at, without);
    }
    md_arena_allocator_destroy(alloc);
}
