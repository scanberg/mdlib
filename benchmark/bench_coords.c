#include "ubench.h"

// The kernels that read or write a state's coordinates, on a large system: what the coordinate
// layout of md_system_state_t costs or saves.

#include <md_gro.h>
#include <md_pdb.h>
#include <md_trr.h>
#include <md_system.h>
#include <md_util.h>
#include <md_script.h>
#include <core/md_coord_stream.h>
#include <core/md_spatial_acc.h>
#include <core/md_arena_allocator.h>
#include <core/md_bitfield.h>
#include <core/md_vec_math.h>

#define BIG_GRO     MD_BENCHMARK_DATA_DIR "/centered.gro"
#define PROTEIN_PDB MD_BENCHMARK_DATA_DIR "/tubulin-A-B.pdb"

static bool load_big(md_system_t* sys, md_system_state_t* st, md_allocator_i* arena) {
    sys->alloc = arena;
    st->alloc  = arena;
    if (!md_gro_system_init_from_file(sys, st, STR_LIT(BIG_GRO))) {
        fprintf(stderr, "Failed to load " BIG_GRO "\n");
        return false;
    }
    return true;
}

// A second state, every atom moved a little, some across the cell boundary
static void jitter(md_system_state_t* dst, const md_system_state_t* src, float amount) {
    md_system_state_copy(dst, src);
    for (size_t i = 0; i < src->num_atoms; ++i) {
        const float s = (float)((i * 2654435761u) % 1000) / 1000.0f - 0.5f;
        dst->xyz[i].x += amount * s;
        dst->xyz[i].y -= amount * s;
        dst->xyz[i].z += amount * 0.5f * s;
    }
}

UBENCH_EX(coords, com_pbc) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));
    md_system_t sys = {0}; md_system_state_t st = {0};
    if (load_big(&sys, &st, arena)) {
        UBENCH_SET_BYTES(st.num_atoms * 12);
        vec3_t com = {0};
        UBENCH_DO_BENCHMARK() {
            com = md_util_com_compute(st.xyz, NULL, NULL, st.num_atoms, &st.unitcell);
            UBENCH_DO_NOTHING(&com);
        }
    }
    md_vm_arena_destroy(arena);
}

UBENCH_EX(coords, com_nopbc) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));
    md_system_t sys = {0}; md_system_state_t st = {0};
    if (load_big(&sys, &st, arena)) {
        UBENCH_SET_BYTES(st.num_atoms * 12);
        vec3_t com = {0};
        UBENCH_DO_BENCHMARK() {
            com = md_util_com_compute(st.xyz, NULL, NULL, st.num_atoms, NULL);
            UBENCH_DO_NOTHING(&com);
        }
    }
    md_vm_arena_destroy(arena);
}

UBENCH_EX(coords, aabb) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));
    md_system_t sys = {0}; md_system_state_t st = {0};
    if (load_big(&sys, &st, arena)) {
        UBENCH_SET_BYTES(st.num_atoms * 12);
        float mn[3], mx[3];
        UBENCH_DO_BENCHMARK() {
            md_util_aabb_compute(mn, mx, st.xyz, NULL, NULL, st.num_atoms);
            UBENCH_DO_NOTHING(mn);
        }
    }
    md_vm_arena_destroy(arena);
}

UBENCH_EX(coords, pbc_wrap) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));
    md_system_t sys = {0}; md_system_state_t st = {0}, moved = {.alloc = arena};
    if (load_big(&sys, &st, arena)) {
        jitter(&moved, &st, 30.0f);
        UBENCH_SET_BYTES(st.num_atoms * 12);
        UBENCH_DO_BENCHMARK() {
            md_util_pbc(moved.xyz, NULL, moved.num_atoms, &st.unitcell);
            UBENCH_DO_NOTHING(moved.xyz);
        }
    }
    md_vm_arena_destroy(arena);
}

UBENCH_EX(coords, interpolate_linear) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));
    md_system_t sys = {0}; md_system_state_t st = {0}, b = {.alloc = arena}, out = {.alloc = arena};
    if (load_big(&sys, &st, arena)) {
        jitter(&b, &st, 3.0f);
        md_system_state_init(&out, st.num_atoms);
        const vec3_t* xyz[2] = {st.xyz, b.xyz};
        UBENCH_SET_BYTES(st.num_atoms * 12 * 3);
        UBENCH_DO_BENCHMARK() {
            md_util_interpolate_linear(out.xyz, xyz, st.num_atoms, &st.unitcell, 0.3f);
            UBENCH_DO_NOTHING(out.xyz);
        }
    }
    md_vm_arena_destroy(arena);
}

UBENCH_EX(coords, interpolate_cubic) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));
    md_system_t sys = {0}; md_system_state_t st = {0}, s[3] = {{.alloc = arena}, {.alloc = arena}, {.alloc = arena}}, out = {.alloc = arena};
    if (load_big(&sys, &st, arena)) {
        for (int i = 0; i < 3; ++i) jitter(&s[i], &st, 1.0f + i);
        md_system_state_init(&out, st.num_atoms);
        const vec3_t* xyz[4] = {st.xyz, s[0].xyz, s[1].xyz, s[2].xyz};
        UBENCH_SET_BYTES(st.num_atoms * 12 * 5);
        UBENCH_DO_BENCHMARK() {
            md_util_interpolate_cubic_spline(out.xyz, xyz, st.num_atoms, &st.unitcell, 0.3f, 1.0f);
            UBENCH_DO_NOTHING(out.xyz);
        }
    }
    md_vm_arena_destroy(arena);
}

UBENCH_EX(coords, rmsd) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));
    md_system_t sys = {0}; md_system_state_t st = {0}, b = {.alloc = arena};
    if (load_big(&sys, &st, arena)) {
        jitter(&b, &st, 1.0f);
        const vec3_t* xyz[2] = {st.xyz, b.xyz};
        UBENCH_SET_BYTES(st.num_atoms * 12 * 2);
        double r = 0;
        UBENCH_DO_BENCHMARK() {
            const vec3_t com[2] = {
                md_util_com_compute(st.xyz, NULL, NULL, st.num_atoms, NULL),
                md_util_com_compute(b.xyz,  NULL, NULL, b.num_atoms,  NULL),
            };
            r = md_util_rmsd_compute(xyz, NULL, NULL, st.num_atoms, com);
            UBENCH_DO_NOTHING(&r);
        }
    }
    md_vm_arena_destroy(arena);
}

UBENCH_EX(coords, sort_spatial) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));
    md_system_t sys = {0}; md_system_state_t st = {0};
    if (load_big(&sys, &st, arena)) {
        uint32_t* idx = md_alloc(arena, st.num_atoms * sizeof(uint32_t));
        UBENCH_DO_BENCHMARK() {
            md_util_sort_spatial(idx, st.xyz, st.num_atoms);
            UBENCH_DO_NOTHING(idx);
        }
    }
    md_vm_arena_destroy(arena);
}

static void count_pairs(const uint32_t* i_idx, const uint32_t* j_idx, const float* ij_dist2, size_t num_pairs, void* user_param) {
    (void)i_idx; (void)j_idx; (void)ij_dist2;
    *(size_t*)user_param += num_pairs;
}

UBENCH_EX(coords, spatial_acc_build_and_pairs) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));
    md_system_t sys = {0}; md_system_state_t st = {0};
    if (load_big(&sys, &st, arena)) {
        md_coord_stream_t stream = md_coord_stream_from_aos((const float*)st.xyz, sizeof(vec3_t), NULL, st.num_atoms);
        size_t count = 0;
        UBENCH_DO_BENCHMARK() {
            md_allocator_i* tmp = md_arena_allocator_create(md_get_heap_allocator(), MEGABYTES(64));
            md_spatial_acc_t acc = { .alloc = tmp };
            md_spatial_acc_init(&acc, &stream, 3.0, &st.unitcell, 0);
            md_spatial_acc_for_each_internal_pair_within_cutoff(&acc, 3.0, count_pairs, &count);
            md_arena_allocator_destroy(tmp);
            UBENCH_DO_NOTHING(&count);
        }
    }
    md_vm_arena_destroy(arena);
}

UBENCH_EX(coords, infer_covalent_bonds) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));
    md_system_t sys = {0}; md_system_state_t st = {0};
    if (load_big(&sys, &st, arena)) {
        UBENCH_DO_BENCHMARK() {
            md_util_system_infer_covalent_bonds(&sys, &st);
            UBENCH_DO_NOTHING(&sys.bond);
        }
    }
    md_vm_arena_destroy(arena);
}

UBENCH_EX(coords, protein_backbone_and_hbonds) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));
    md_system_t sys = {.alloc = arena};
    md_system_state_t st = {.alloc = arena};
    if (md_pdb_system_init_from_file(&sys, &st, STR_LIT(PROTEIN_PDB), MD_PDB_OPTION_DISABLE_CACHE_FILE_WRITE) &&
        md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL)) {
        const size_t n = sys.protein_backbone.segment.count;
        md_secondary_structure_t* ss = md_alloc(arena, n * sizeof(md_secondary_structure_t));
        md_backbone_angles_t* ang = md_alloc(arena, n * sizeof(md_backbone_angles_t));
        UBENCH_DO_BENCHMARK() {
            md_util_backbone_secondary_structure_infer(ss, n, st.xyz, &st.unitcell, &sys.protein_backbone);
            md_util_backbone_angles_compute(ang, n, st.xyz, &st.unitcell, &sys.protein_backbone);
            md_util_hydrogen_bond_infer(&sys.hydrogen_bond, st.xyz, &st.unitcell, 3.0, 120.0);
            UBENCH_DO_NOTHING(ss);
        }
    } else {
        fprintf(stderr, "Failed to load " PROTEIN_PDB "\n");
    }
    md_vm_arena_destroy(arena);
}

// Script evaluation over a run: extraction of each frame and the coordinate kernels behind the
// properties. The TRR is uncompressed, so decoding is a small part of it.
UBENCH_EX(coords, script_eval_trr) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));
    md_system_t sys = {.alloc = arena};
    md_system_state_t st = {.alloc = arena};
    const str_t run = STR_LIT("run/t");
    if (md_gro_system_init_from_file(&sys, &st, STR_LIT(MD_BENCHMARK_DATA_DIR "/tryptophan-md.gro")) &&
        md_util_system_infer(&sys, &st, MD_UTIL_INFER_ALL) &&
        md_trr_system_publish_run(&sys, STR_LIT(MD_BENCHMARK_DATA_DIR "/tryptophan-md.trr"), run, MD_RUN_FLAG_DISABLE_CACHE_WRITE)) {
        md_script_ir_t* ir = md_script_ir_create(arena);
        const str_t src = STR_LIT(
            "c = com(all);\n"
            "d = distance(com(residue(1)), com(all));\n"
            "g = rdf(element('C'), element('O'), 8.0);\n"
            "r = rmsd(all);\n");
        if (md_script_ir_compile_from_source(ir, src, &sys, NULL)) {
            const md_attribute_t* time = md_attributes_find(&sys.attributes, STR_LIT("run/t/time"));
            const uint32_t F = time ? (uint32_t)time->format.shape[0] : 0;
            md_script_eval_t* eval = md_script_eval_create(F, ir, arena);
            UBENCH_DO_BENCHMARK() {
                md_script_eval_frame_range(eval, ir, &sys, run, 0, F);
                UBENCH_DO_NOTHING(eval);
            }
            md_script_eval_free(eval);
        } else {
            fprintf(stderr, "script did not compile\n");
        }
    }
    md_vm_arena_destroy(arena);
}

// One frame after another into a state, positions and cell
UBENCH_EX(coords, extract_trr) {
    md_allocator_i* arena = md_vm_arena_create(GIGABYTES(4));
    md_system_t sys = {.alloc = arena};
    md_system_state_t st = {.alloc = arena};
    const str_t run = STR_LIT("run/t");
    if (md_gro_system_init_from_file(&sys, &st, STR_LIT(MD_BENCHMARK_DATA_DIR "/tryptophan-md.gro")) &&
        md_trr_system_publish_run(&sys, STR_LIT(MD_BENCHMARK_DATA_DIR "/tryptophan-md.trr"), run, MD_RUN_FLAG_DISABLE_CACHE_WRITE)) {
        const str_t paths[] = { STR_LIT("atom/position"), STR_LIT("unitcell") };
        md_system_extract_t* ex = md_system_extract_begin(&sys, run, paths, 2, md_get_heap_allocator());
        UBENCH_DO_BENCHMARK() {
            for (int64_t f = 0; f < 101; ++f) md_system_extract_frame(ex, f, &st);
            UBENCH_DO_NOTHING(&st);
        }
        md_system_extract_end(ex);
    }
    md_vm_arena_destroy(arena);
}
