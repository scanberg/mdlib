#include "utest.h"

#include "qm_test_util.h"

#include <md_molden.h>
#include <md_topo.h>

#include <core/md_allocator.h>
#include <core/md_str.h>

#include <math.h>
#include <string.h>

// Certified critical points of a GTO density (md_topo_compute_extremum_graph_gto), driven exactly the way
// a consumer drives it: a QM file is read into a system, and the basis, the total density matrix and the
// geometry come back out of the attribute table. Water, cc-pVDZ, RHF (see test_data/molden/README.md):
// the topology is three nuclear attractors and two O-H bond critical points, nothing else.

typedef struct topo_gto_input_t {
    qm_test_t      t;
    md_gto_basis_t basis;
    double*        density;
    float          xyz[3 * 3];  // Bohr
} topo_gto_input_t;

static bool topo_gto_load_water(topo_gto_input_t* in) {
    memset(in, 0, sizeof(*in));
    qm_test_init(&in->t, MEGABYTES(16));
    if (!md_molden_system_init_from_file(&in->t.sys, &in->t.state, str_from_cstr(MD_UNITTEST_DATA_DIR "/molden/h2o_ccpvdz.molden"))) return false;
    if (!qm_test_basis(&in->basis, &in->t)) return false;
    size_t dim = 0;
    in->density = qm_test_matrix(&in->t, STR_LIT("orbital/total/density"), &dim);
    if (!in->density || dim != md_gto_basis_num_ao(&in->basis)) return false;
    double xyz[9] = {0};
    const md_attribute_t* coord = qm_test_attr(&in->t, STR_LIT("qm/atom/coordinate"));   // Angstrom
    if (!coord || md_attribute_extract_f64(xyz, 9, coord, md_attribute_slice_all(), md_unit_none()) != 9) return false;
    for (int i = 0; i < 9; ++i) in->xyz[i] = (float)(xyz[i] * QM_TEST_ANGSTROM_TO_BOHR);
    return true;
}

static md_topo_gto_desc_t topo_gto_desc(const topo_gto_input_t* in) {
    return (md_topo_gto_desc_t){
        .basis = &in->basis,
        .atom_xyz = in->xyz,
        .density_matrix = in->density,
        .rho_min = 1.0e-4,
        .h_min = 1.0e-4,
        .trace_separatrices = true,
    };
}

static double topo_gto_dist(const md_topo_vert_t* v, const float* p) {
    const double dx = v->x - p[0], dy = v->y - p[1], dz = v->z - p[2];
    return sqrt(dx * dx + dy * dy + dz * dz);
}

UTEST(topo_gto, water_topology_is_certified) {
    topo_gto_input_t in;
    ASSERT_TRUE(topo_gto_load_water(&in));
    const md_topo_gto_desc_t desc = topo_gto_desc(&in);

    md_topo_extremum_graph_t graph = { .alloc = in.t.alloc };
    md_topo_gto_info_t info;
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&graph, &info, &desc));

    // Every cube was decided: nothing near-degenerate here, so the result is proven complete.
    EXPECT_TRUE(info.complete);
    EXPECT_EQ(0u, info.num_unresolved_boxes);
    EXPECT_EQ(1, info.poincare_hopf);
    EXPECT_GT(info.domain_pad, 3.0);

    uint32_t count[MD_TOPO_NUM_TYPES] = {0};
    md_topo_count_vertex_types(count, &graph);
    EXPECT_EQ(3u, count[MD_TOPO_MAXIMUM]);
    EXPECT_EQ(2u, count[MD_TOPO_SPLIT_SADDLE]);
    EXPECT_EQ(0u, count[MD_TOPO_JOIN_SADDLE]);
    EXPECT_EQ(0u, count[MD_TOPO_MINIMUM]);

    // Each maximum is a nuclear attractor. A GTO density has no cusp, so attractors are displaced from the
    // nuclei: negligibly for O, by several hundredths of a Bohr for H.
    int max_of_atom[3] = { -1, -1, -1 };
    for (uint32_t i = 0; i < graph.num_vertices; ++i) {
        if (graph.types[i] != MD_TOPO_MAXIMUM) continue;
        int best = -1; double bd = 1e9;
        for (int a = 0; a < 3; ++a) {
            const double d = topo_gto_dist(&graph.vertices[i], in.xyz + 3 * a);
            if (d < bd) { bd = d; best = a; }
        }
        ASSERT_GE(best, 0);
        EXPECT_EQ(-1, max_of_atom[best]);
        max_of_atom[best] = (int)i;
        EXPECT_LT(bd, best == 0 ? 1.0e-3 : 0.15);
    }

    // Each bond critical point lies between O and one H, and its two separatrices (the bond path) end in
    // exactly those two attractors.
    uint32_t bcp_seen = 0;
    for (uint32_t i = 0; i < graph.num_vertices; ++i) {
        if (graph.types[i] != MD_TOPO_SPLIT_SADDLE) continue;
        bcp_seen++;
        EXPECT_GT(graph.vertices[i].value, 0.2f);   // a covalent O-H bond, rho_b ~ 0.37
        int ends[2] = { -1, -1 }, n = 0;
        for (uint32_t e = 0; e < graph.num_edges; ++e) {
            if (graph.edges[e].from == i && n < 2) ends[n++] = (int)graph.edges[e].to;
        }
        ASSERT_EQ(2, n);
        const bool o_end = ends[0] == max_of_atom[0] || ends[1] == max_of_atom[0];
        const bool h_end = ends[0] == max_of_atom[1] || ends[1] == max_of_atom[1] || ends[0] == max_of_atom[2] || ends[1] == max_of_atom[2];
        EXPECT_TRUE(o_end);
        EXPECT_TRUE(h_end);
    }
    EXPECT_EQ(2u, bcp_seen);
    EXPECT_EQ(4u, graph.num_edges);

    md_topo_extremum_graph_free(&graph);
    qm_test_free(&in.t);
}

// Same input, same output, bit for bit, whatever the thread count: cubes are decided independently and
// their results merged in cube order, with no atomics and no ordering that depends on scheduling.
UTEST(topo_gto, deterministic_across_thread_counts) {
    topo_gto_input_t in;
    ASSERT_TRUE(topo_gto_load_water(&in));
    md_topo_gto_desc_t da = topo_gto_desc(&in);
    md_topo_gto_desc_t db = topo_gto_desc(&in);
    da.num_threads = 1;
    db.num_threads = 4;

    md_topo_extremum_graph_t a = { .alloc = in.t.alloc };
    md_topo_extremum_graph_t b = { .alloc = in.t.alloc };
    md_topo_gto_info_t ia, ib;
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&a, &ia, &da));
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&b, &ib, &db));
    EXPECT_EQ(1u, ia.num_threads);
    EXPECT_EQ(4u, ib.num_threads);
    ASSERT_EQ(a.num_vertices, b.num_vertices);
    ASSERT_EQ(a.num_edges, b.num_edges);
    EXPECT_EQ(ia.num_box_evals, ib.num_box_evals);
    EXPECT_EQ(0, memcmp(a.vertices, b.vertices, sizeof(md_topo_vert_t) * a.num_vertices));
    EXPECT_EQ(0, memcmp(a.types, b.types, sizeof(md_topo_critical_point_type_t) * a.num_vertices));
    EXPECT_EQ(0, memcmp(a.edges, b.edges, sizeof(md_topo_edge_t) * a.num_edges));

    md_topo_extremum_graph_free(&a);
    md_topo_extremum_graph_free(&b);
    qm_test_free(&in.t);
}

// A density threshold above the bond critical points removes them, and only them: what is reported is
// exactly the set of CPs with rho >= rho_min.
UTEST(topo_gto, rho_min_is_respected) {
    topo_gto_input_t in;
    ASSERT_TRUE(topo_gto_load_water(&in));
    md_topo_gto_desc_t desc = topo_gto_desc(&in);
    desc.rho_min = 0.5;
    desc.trace_separatrices = false;

    md_topo_extremum_graph_t graph = { .alloc = in.t.alloc };
    md_topo_gto_info_t info;
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&graph, &info, &desc));
    uint32_t count[MD_TOPO_NUM_TYPES] = {0};
    md_topo_count_vertex_types(count, &graph);
    EXPECT_EQ(1u, count[MD_TOPO_MAXIMUM]);     // only the O attractor has rho > 0.5
    EXPECT_EQ(0u, count[MD_TOPO_SPLIT_SADDLE]);
    for (uint32_t i = 0; i < graph.num_vertices; ++i) EXPECT_GE(graph.vertices[i].value, 0.5f);

    md_topo_extremum_graph_free(&graph);
    qm_test_free(&in.t);
}

// A cancel flag that is already raised stops the search before any cube is decided: the call reports
// failure, says it was cancelled, and never claims completeness.
UTEST(topo_gto, cancel) {
    topo_gto_input_t in;
    ASSERT_TRUE(topo_gto_load_water(&in));
    md_topo_gto_desc_t desc = topo_gto_desc(&in);
    volatile int32_t cancel = 1;
    desc.cancel = &cancel;

    md_topo_extremum_graph_t graph = { .alloc = in.t.alloc };
    md_topo_gto_info_t info;
    EXPECT_FALSE(md_topo_compute_extremum_graph_gto(&graph, &info, &desc));
    EXPECT_TRUE(info.cancelled);
    EXPECT_FALSE(info.complete);

    md_topo_extremum_graph_free(&graph);
    qm_test_free(&in.t);
}

static bool topo_gto_same_graph(const md_topo_extremum_graph_t* a, const md_topo_extremum_graph_t* b) {
    // Same CPs (matched by position: symmetry-equivalent CPs of equal density may be listed in another
    // order) and the same edges between them.
    if (a->num_vertices != b->num_vertices || a->num_edges != b->num_edges || a->num_vertices > 64) return false;
    uint32_t map[64];
    bool used[64] = { false };
    for (uint32_t i = 0; i < a->num_vertices; ++i) {
        double best = 1e9; uint32_t bj = 0;
        for (uint32_t j = 0; j < b->num_vertices; ++j) {
            if (a->types[i] != b->types[j]) continue;
            const double dx = a->vertices[i].x - b->vertices[j].x, dy = a->vertices[i].y - b->vertices[j].y, dz = a->vertices[i].z - b->vertices[j].z;
            const double d = sqrt(dx * dx + dy * dy + dz * dz);
            if (d < best) { best = d; bj = j; }
        }
        if (best > 1e-5 || used[bj]) return false;
        if (fabs(a->vertices[i].value - b->vertices[bj].value) > 1e-6 * (1.0 + fabs(a->vertices[i].value))) return false;
        used[bj] = true;
        map[i] = bj;
    }
    for (uint32_t e = 0; e < a->num_edges; ++e) {
        bool found = false;
        for (uint32_t f = 0; f < b->num_edges && !found; ++f) found = b->edges[f].from == map[a->edges[e].from] && b->edges[f].to == map[a->edges[e].to];
        if (!found) return false;
    }
    return true;
}

// The factored form (D = sum_k l_k c_k c_k^T; rank 5 here, the occupied orbitals of RHF water) proves
// the same topology as the matrix form, with the same CPs and edges.
UTEST(topo_gto, factored_density_matches_matrix) {
    topo_gto_input_t in;
    ASSERT_TRUE(topo_gto_load_water(&in));
    md_topo_gto_desc_t dm = topo_gto_desc(&in), df = topo_gto_desc(&in);
    dm.density_form = MD_TOPO_GTO_DENSITY_MATRIX;
    df.density_form = MD_TOPO_GTO_DENSITY_FACTORED;

    md_topo_extremum_graph_t m = { .alloc = in.t.alloc }, f = { .alloc = in.t.alloc };
    md_topo_gto_info_t im, inf;
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&m, &im, &dm));
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&f, &inf, &df));
    EXPECT_EQ(0u, im.density_rank);
    EXPECT_EQ(0u, im.num_factored_evals);
    EXPECT_EQ(5u, inf.density_rank);
    EXPECT_EQ(inf.num_box_evals, inf.num_factored_evals);
    EXPECT_TRUE(inf.complete);
    EXPECT_TRUE(topo_gto_same_graph(&m, &f));

    md_topo_extremum_graph_free(&m);
    md_topo_extremum_graph_free(&f);
    qm_test_free(&in.t);
}

// A density matrix of full rank (here water's plus 1e-3 I) is not factored by auto, where it would not
// pay; forced, its full-rank factorization proves the same topology as the matrix.
UTEST(topo_gto, full_rank_density) {
    topo_gto_input_t in;
    ASSERT_TRUE(topo_gto_load_water(&in));
    const size_t N = md_gto_basis_num_ao(&in.basis);
    for (size_t i = 0; i < N; ++i) in.density[i * N + i] += 1.0e-3;
    md_topo_gto_desc_t da = topo_gto_desc(&in), dm = topo_gto_desc(&in), df = topo_gto_desc(&in);
    da.rho_min = dm.rho_min = df.rho_min = 1.0e-3;
    dm.density_form = MD_TOPO_GTO_DENSITY_MATRIX;
    df.density_form = MD_TOPO_GTO_DENSITY_FACTORED;

    md_topo_extremum_graph_t a = { .alloc = in.t.alloc }, m = { .alloc = in.t.alloc }, f = { .alloc = in.t.alloc };
    md_topo_gto_info_t ia, im, inf;
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&a, &ia, &da));
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&m, &im, &dm));
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&f, &inf, &df));
    EXPECT_EQ(0u, ia.density_rank);
    EXPECT_EQ(ia.num_box_evals, im.num_box_evals);
    EXPECT_EQ((uint32_t)N, inf.density_rank);
    EXPECT_TRUE(inf.complete);
    EXPECT_TRUE(topo_gto_same_graph(&m, &f));

    md_topo_extremum_graph_free(&a);
    md_topo_extremum_graph_free(&m);
    md_topo_extremum_graph_free(&f);
    qm_test_free(&in.t);
}

#if MD_ENABLE_GPU
#include <core/md_gpu.h>

// The GPU sweep (fp32 with rounding margins; roots polished in double on the CPU) reaches the same
// certified topology as the CPU reference, and is itself deterministic.
UTEST(topo_gto, gpu_matches_cpu) {
    md_gpu_device_t dev = md_gpu_device_create(&(md_gpu_device_desc_t){ .label = "topo_gto unittest" });
    if (!dev) UTEST_SKIP("no GPU device");
    md_gpu_stream_t stream = md_gpu_stream_default(dev, MD_GPU_STREAM_COMPUTE);

    topo_gto_input_t in;
    ASSERT_TRUE(topo_gto_load_water(&in));
    md_topo_gto_desc_t desc = topo_gto_desc(&in);
    desc.rho_min = 1.0e-3;      // the same 5 CPs in a smaller domain: quick on software Vulkan too

    md_topo_extremum_graph_t cpu = { .alloc = in.t.alloc }, gpu = { .alloc = in.t.alloc }, gpu2 = { .alloc = in.t.alloc };
    md_topo_gto_info_t ic, ig, ig2;
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto(&cpu, &ic, &desc));
    ASSERT_TRUE(md_topo_compute_extremum_graph_gto_gpu(&gpu, &ig, &desc, stream));
    EXPECT_TRUE(ig.used_gpu);
    EXPECT_TRUE(ig.complete);
    EXPECT_GT(ig.num_gpu_box_evals, 0u);
    EXPECT_EQ(1, ig.poincare_hopf);
    EXPECT_TRUE(topo_gto_same_graph(&cpu, &gpu));

    ASSERT_TRUE(md_topo_compute_extremum_graph_gto_gpu(&gpu2, &ig2, &desc, stream));
    EXPECT_EQ(ig.num_box_evals, ig2.num_box_evals);
    ASSERT_EQ(gpu.num_vertices, gpu2.num_vertices);
    ASSERT_EQ(gpu.num_edges, gpu2.num_edges);
    EXPECT_EQ(0, memcmp(gpu.vertices, gpu2.vertices, sizeof(md_topo_vert_t) * gpu.num_vertices));
    EXPECT_EQ(0, memcmp(gpu.edges, gpu2.edges, sizeof(md_topo_edge_t) * gpu.num_edges));

    md_topo_extremum_graph_free(&cpu);
    md_topo_extremum_graph_free(&gpu);
    md_topo_extremum_graph_free(&gpu2);
    qm_test_free(&in.t);
    md_topo_gpu_shutdown();      // releases the kernel before its device
    md_gpu_device_destroy(dev);
}
#endif
