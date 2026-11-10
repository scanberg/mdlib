#include "utest.h"

#include <md_flow.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>

#include <math.h>

// A three-column graph shaped like the real thing: two groups of two atoms donating, three
// orbitals in the middle, two groups of two atoms accepting. Every band is a genuine joint
// distribution, so the marginals agree with the node weights by construction rather than by luck.
//
//   col 0                col 1          col 2
//   A { a0 .3, a1 .2 }   o0 .5          C { c0 .25, c1 .25 }
//   B { b0 .1, b1 .4 }   o1 .3          D { d0 .40, d1 .10 }
//                        o2 .2

typedef struct test_graph_t {
    md_flow_graph_t graph;
    uint32_t A, a0, a1, B, b0, b1;
    uint32_t o0, o1, o2;
    uint32_t C, c0, c1, D, d0, d1;
} test_graph_t;

static uint32_t add(md_flow_graph_t* g, uint32_t column, uint32_t parent, uint32_t level, float weight, uint64_t key) {
    md_flow_node_t node = {
        .column = column,
        .parent = parent,
        .level  = level,
        .weight = weight,
        .color  = {1,1,1,1},
        .label  = MD_STR_INIT("n"),
        .key    = key,
    };
    return md_flow_graph_add_node(g, &node);
}

static void build(test_graph_t* t, md_allocator_i* alloc) {
    md_flow_graph_t* g = &t->graph;
    md_flow_graph_init(g, 3, alloc);

    t->A  = add(g, 0, MD_FLOW_INVALID_INDEX, 0, 0.5f, 1);
    t->a0 = add(g, 0, t->A, 1, 0.3f, 2);
    t->a1 = add(g, 0, t->A, 1, 0.2f, 3);
    t->B  = add(g, 0, MD_FLOW_INVALID_INDEX, 0, 0.5f, 4);
    t->b0 = add(g, 0, t->B, 1, 0.1f, 5);
    t->b1 = add(g, 0, t->B, 1, 0.4f, 6);

    t->o0 = add(g, 1, MD_FLOW_INVALID_INDEX, 0, 0.5f, 10);
    t->o1 = add(g, 1, MD_FLOW_INVALID_INDEX, 0, 0.3f, 11);
    t->o2 = add(g, 1, MD_FLOW_INVALID_INDEX, 0, 0.2f, 12);

    t->C  = add(g, 2, MD_FLOW_INVALID_INDEX, 0, 0.5f,  20);
    t->c0 = add(g, 2, t->C, 1, 0.25f, 21);
    t->c1 = add(g, 2, t->C, 1, 0.25f, 22);
    t->D  = add(g, 2, MD_FLOW_INVALID_INDEX, 0, 0.5f,  23);
    t->d0 = add(g, 2, t->D, 1, 0.40f, 24);
    t->d1 = add(g, 2, t->D, 1, 0.10f, 25);

    // Band 0: rows are col-0 atoms, columns are orbitals.
    md_flow_graph_add_link(g, t->a0, t->o0, 0.20f);
    md_flow_graph_add_link(g, t->a0, t->o1, 0.10f);
    md_flow_graph_add_link(g, t->a1, t->o0, 0.10f);
    md_flow_graph_add_link(g, t->a1, t->o1, 0.05f);
    md_flow_graph_add_link(g, t->a1, t->o2, 0.05f);
    md_flow_graph_add_link(g, t->b0, t->o1, 0.05f);
    md_flow_graph_add_link(g, t->b0, t->o2, 0.05f);
    md_flow_graph_add_link(g, t->b1, t->o0, 0.20f);
    md_flow_graph_add_link(g, t->b1, t->o1, 0.10f);
    md_flow_graph_add_link(g, t->b1, t->o2, 0.10f);

    // Band 1: rows are orbitals, columns are col-2 atoms.
    md_flow_graph_add_link(g, t->o0, t->c0, 0.20f);
    md_flow_graph_add_link(g, t->o0, t->c1, 0.10f);
    md_flow_graph_add_link(g, t->o0, t->d0, 0.15f);
    md_flow_graph_add_link(g, t->o0, t->d1, 0.05f);
    md_flow_graph_add_link(g, t->o1, t->c0, 0.05f);
    md_flow_graph_add_link(g, t->o1, t->c1, 0.10f);
    md_flow_graph_add_link(g, t->o1, t->d0, 0.15f);
    md_flow_graph_add_link(g, t->o2, t->c1, 0.05f);
    md_flow_graph_add_link(g, t->o2, t->d0, 0.10f);
    md_flow_graph_add_link(g, t->o2, t->d1, 0.05f);
}

#define EPS 1.0e-4f

UTEST(flow, graph_validates) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    test_graph_t t = {0};
    build(&t, alloc);

    EXPECT_EQ(md_flow_graph_validate(&t.graph, EPS), MD_FLOW_OK);

    md_arena_allocator_destroy(alloc);
}

UTEST(flow, rejects_bad_hierarchy) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_flow_graph_t g = {0};
    md_flow_graph_init(&g, 2, alloc);

    const uint32_t root = add(&g, 0, MD_FLOW_INVALID_INDEX, 0, 1.0f, 1);
    // A parent in a different column must be refused at insertion time.
    EXPECT_EQ(add(&g, 1, root, 1, 1.0f, 2), MD_FLOW_INVALID_INDEX);
    // So must a column that does not exist.
    EXPECT_EQ(add(&g, 5, MD_FLOW_INVALID_INDEX, 0, 1.0f, 3), MD_FLOW_INVALID_INDEX);

    md_arena_allocator_destroy(alloc);
}

UTEST(flow, rejects_skipped_column) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_flow_graph_t g = {0};
    md_flow_graph_init(&g, 3, alloc);

    const uint32_t a = add(&g, 0, MD_FLOW_INVALID_INDEX, 0, 1.0f, 1);
    add(&g, 1, MD_FLOW_INVALID_INDEX, 0, 1.0f, 2);
    const uint32_t c = add(&g, 2, MD_FLOW_INVALID_INDEX, 0, 1.0f, 3);
    md_flow_graph_add_link(&g, a, c, 1.0f);     // column 0 -> 2, skipping the middle

    EXPECT_EQ(md_flow_graph_validate(&g, EPS), MD_FLOW_ERR_LAYERING);

    md_arena_allocator_destroy(alloc);
}

UTEST(flow, rejects_broken_partition) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    test_graph_t t = {0};
    build(&t, alloc);

    // Move the PARENT, not a child: a leaf's weight also feeds the column sum and its own link
    // totals, so shifting one would trip three checks at once and prove nothing about this one.
    t.graph.nodes[t.A].weight = 0.6f;
    EXPECT_EQ(md_flow_graph_validate(&t.graph, EPS), MD_FLOW_ERR_PARTITION);

    md_arena_allocator_destroy(alloc);
}

UTEST(flow, rejects_node_flow_mismatch) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    test_graph_t t = {0};
    build(&t, alloc);

    // A node whose ribbons do not fill it would draw flows spilling out of its own bar.
    t.graph.links[0].weight = 0.25f;
    const md_flow_result_t result = md_flow_graph_validate(&t.graph, EPS);
    EXPECT_TRUE(result == MD_FLOW_ERR_NODE_FLOW || result == MD_FLOW_ERR_BAND_SUM);

    md_arena_allocator_destroy(alloc);
}

// The invariant the whole design rests on: whatever the user expands, the diagram still adds up.
UTEST(flow, conserved_on_every_cut) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    test_graph_t t = {0};
    build(&t, alloc);

    md_flow_cut_t cut = {0};
    md_flow_cut_init(&cut, alloc);

    // Four nodes in this graph have children, so there are 16 distinct cuts. Check them all.
    const uint32_t expandable[4] = { t.A, t.B, t.C, t.D };
    for (uint32_t mask = 0; mask < 16; ++mask) {
        md_flow_cut_collapse_all(&cut);
        for (uint32_t i = 0; i < 4; ++i) {
            md_flow_cut_set_expanded(&cut, expandable[i], (mask & (1u << i)) != 0);
        }
        md_flow_cut_resolve(&cut, &t.graph);

        EXPECT_EQ(md_flow_cut_validate(&cut, &t.graph, EPS), MD_FLOW_OK);

        double total = 0.0;
        for (size_t i = 0; i < md_array_size(cut.links); ++i) {
            total += (double)cut.links[i].weight;
        }
        // Two bands, each summing to one.
        EXPECT_NEAR(total, 2.0, (double)EPS);
    }

    md_flow_cut_free(&cut);
    md_arena_allocator_destroy(alloc);
}

// Aggregation must be exact, not approximate: the collapsed A -> o0 flow is the sum of the leaf
// flows it contains, and nothing else.
UTEST(flow, aggregation_is_a_sum_of_leaves) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    test_graph_t t = {0};
    build(&t, alloc);

    md_flow_cut_t cut = {0};
    md_flow_cut_init(&cut, alloc);
    md_flow_cut_collapse_all(&cut);
    md_flow_cut_resolve(&cut, &t.graph);

    float a_to_o0 = 0.0f;
    for (size_t i = 0; i < md_array_size(cut.links); ++i) {
        if (cut.links[i].src == t.A && cut.links[i].dst == t.o0) {
            a_to_o0 = cut.links[i].weight;
        }
    }
    // a0 -> o0 is 0.20 and a1 -> o0 is 0.10.
    EXPECT_NEAR((double)a_to_o0, 0.30, (double)EPS);

    // And expanding A must split that one ribbon into exactly those two.
    md_flow_cut_set_expanded(&cut, t.A, true);
    md_flow_cut_resolve(&cut, &t.graph);

    float a0_to_o0 = 0.0f, a1_to_o0 = 0.0f;
    for (size_t i = 0; i < md_array_size(cut.links); ++i) {
        if (cut.links[i].dst != t.o0) continue;
        if (cut.links[i].src == t.a0) a0_to_o0 = cut.links[i].weight;
        if (cut.links[i].src == t.a1) a1_to_o0 = cut.links[i].weight;
    }
    EXPECT_NEAR((double)a0_to_o0, 0.20, (double)EPS);
    EXPECT_NEAR((double)a1_to_o0, 0.10, (double)EPS);
    EXPECT_NEAR((double)(a0_to_o0 + a1_to_o0), (double)a_to_o0, (double)EPS);

    md_flow_cut_free(&cut);
    md_arena_allocator_destroy(alloc);
}

UTEST(flow, threshold_folds_into_others_without_losing_weight) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    test_graph_t t = {0};
    build(&t, alloc);

    md_flow_cut_t cut = {0};
    md_flow_cut_init(&cut, alloc);
    md_flow_cut_expand_all(&cut, &t.graph);
    // Folding is per parent and never a single node: folds a0+a1 inside A, o1+o2 among the roots, c0+c1
    // inside C, and leaves b0 and d1, each alone below it in their group
    cut.threshold = 0.35f;
    md_flow_cut_resolve(&cut, &t.graph);

    EXPECT_EQ(md_array_size(cut.other), (size_t)3);
    for (size_t i = 0; i < md_array_size(cut.other); ++i) {
        // inside its parent's slot: an "Others" of atoms has the group as its parent
        if (cut.other[i].column == 0) EXPECT_EQ(cut.other[i].parent, t.A);
        if (cut.other[i].column == 1) EXPECT_EQ(cut.other[i].parent, MD_FLOW_INVALID_INDEX);
        if (cut.other[i].column == 2) EXPECT_EQ(cut.other[i].parent, t.C);
    }
    // Conservation is what makes "Others" safe rather than a rounding hole.
    EXPECT_EQ(md_flow_cut_validate(&cut, &t.graph, EPS), MD_FLOW_OK);

    md_flow_cut_free(&cut);
    md_arena_allocator_destroy(alloc);
}

UTEST(flow, layout_ribbons_match_at_both_ends) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    test_graph_t t = {0};
    build(&t, alloc);

    md_flow_cut_t cut = {0};
    md_flow_cut_init(&cut, alloc);
    md_flow_cut_expand_all(&cut, &t.graph);
    md_flow_cut_resolve(&cut, &t.graph);

    md_flow_layout_t layout = {0};
    md_flow_layout_init(&layout, alloc);
    md_flow_layout_compute(&layout, &t.graph, &cut, NULL);

    EXPECT_EQ(md_array_size(layout.nodes),   md_array_size(cut.visible));
    EXPECT_EQ(md_array_size(layout.ribbons), md_array_size(cut.links));

    for (size_t i = 0; i < md_array_size(layout.ribbons); ++i) {
        const md_flow_layout_ribbon_t* r = layout.ribbons + i;
        const float src_h = r->p1.y - r->p0.y;
        const float dst_h = r->q1.y - r->q0.y;
        // A ribbon that is not the same thickness at both ends is a ribbon that lies about how
        // much charge moved.
        EXPECT_NEAR((double)src_h, (double)dst_h, 1.0e-5);
        EXPECT_TRUE(src_h > 0.0f);
    }

    // Everything stays inside the unit square the layout promises.
    for (size_t i = 0; i < md_array_size(layout.nodes); ++i) {
        const md_flow_layout_node_t* n = layout.nodes + i;
        EXPECT_TRUE(n->min.x >= -1.0e-5f && n->max.x <= 1.0f + 1.0e-5f);
        EXPECT_TRUE(n->min.y >= -1.0e-5f && n->max.y <= 1.0f + 1.0e-5f);
        EXPECT_TRUE(n->max.y >= n->min.y);
    }

    md_flow_layout_free(&layout);
    md_flow_cut_free(&cut);
    md_arena_allocator_destroy(alloc);
}

// Ribbon stubs are sub-allocated inside their node's own span; if that ever drifts, ribbons draw
// outside the bar they belong to and the picture silently stops being a Sankey.
UTEST(flow, layout_stubs_stay_inside_their_node) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    test_graph_t t = {0};
    build(&t, alloc);

    md_flow_cut_t cut = {0};
    md_flow_cut_init(&cut, alloc);
    md_flow_cut_collapse_all(&cut);
    md_flow_cut_resolve(&cut, &t.graph);

    md_flow_layout_t layout = {0};
    md_flow_layout_init(&layout, alloc);
    md_flow_layout_compute(&layout, &t.graph, &cut, NULL);

    for (size_t i = 0; i < md_array_size(layout.ribbons); ++i) {
        const md_flow_layout_ribbon_t* r = layout.ribbons + i;
        const md_flow_link_t* link = cut.links + r->link_idx;

        const md_flow_layout_node_t* src = NULL;
        const md_flow_layout_node_t* dst = NULL;
        for (size_t j = 0; j < md_array_size(layout.nodes); ++j) {
            if (layout.nodes[j].cut_idx == link->src) src = layout.nodes + j;
            if (layout.nodes[j].cut_idx == link->dst) dst = layout.nodes + j;
        }
        ASSERT_TRUE(src != NULL);
        ASSERT_TRUE(dst != NULL);

        EXPECT_TRUE(r->p0.y >= src->min.y - 1.0e-5f && r->p1.y <= src->max.y + 1.0e-5f);
        EXPECT_TRUE(r->q0.y >= dst->min.y - 1.0e-5f && r->q1.y <= dst->max.y + 1.0e-5f);
    }

    md_flow_layout_free(&layout);
    md_flow_cut_free(&cut);
    md_arena_allocator_destroy(alloc);
}

UTEST(flow, empty_graph_is_not_a_special_case) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1));
    md_flow_graph_t g = {0};
    md_flow_graph_init(&g, 3, alloc);

    md_flow_cut_t cut = {0};
    md_flow_cut_init(&cut, alloc);
    md_flow_cut_resolve(&cut, &g);
    EXPECT_EQ(md_array_size(cut.visible), 0u);
    EXPECT_EQ(md_flow_cut_validate(&cut, &g, EPS), MD_FLOW_OK);

    md_flow_layout_t layout = {0};
    md_flow_layout_init(&layout, alloc);
    md_flow_layout_compute(&layout, &g, &cut, NULL);
    EXPECT_EQ(md_array_size(layout.nodes), 0u);

    md_flow_layout_free(&layout);
    md_flow_cut_free(&cut);
    md_arena_allocator_destroy(alloc);
}

// Expanding a group must show its atoms in the group's own slot, not scattered down the column.
// Crossing reduction orders by barycentre, and left to itself it has no reason to keep siblings
// together - so this pins the two-level ordering that makes an expansion read as a closer look at
// the same diagram rather than as a different one.
UTEST(flow, expansion_keeps_siblings_together) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    test_graph_t t = {0};
    build(&t, alloc);

    md_flow_cut_t cut = {0};
    md_flow_cut_init(&cut, alloc);
    md_flow_cut_expand_all(&cut, &t.graph);
    md_flow_cut_resolve(&cut, &t.graph);

    md_flow_layout_t layout = {0};
    md_flow_layout_init(&layout, alloc);
    md_flow_layout_compute(&layout, &t.graph, &cut, NULL);

    // Walk each column in drawn order; every time the outermost ancestor changes, that ancestor
    // must not have been seen before - which is exactly "each family occupies one run".
    for (uint32_t c = 0; c < t.graph.num_columns; ++c) {
        uint32_t seen[64] = {0};
        size_t num_seen = 0;
        uint32_t prev_root = MD_FLOW_INVALID_INDEX;

        for (size_t i = 0; i < md_array_size(layout.nodes); ++i) {
            const uint32_t n = layout.nodes[i].cut_idx;
            const md_flow_node_t* node = md_flow_cut_node(&cut, &t.graph, n);
            if (node->column != c) continue;

            uint32_t root = n;
            if (n < md_array_size(t.graph.nodes)) {
                while (t.graph.nodes[root].parent != MD_FLOW_INVALID_INDEX) {
                    root = t.graph.nodes[root].parent;
                }
            }
            if (root == prev_root) continue;

            for (size_t j = 0; j < num_seen; ++j) {
                EXPECT_NE(seen[j], root);   // this family already had its run and is starting another
            }
            seen[num_seen++] = root;
            prev_root = root;
        }
    }

    md_flow_layout_free(&layout);
    md_flow_cut_free(&cut);
    md_arena_allocator_destroy(alloc);
}

// A partially expanded cut is the normal case, not a special one: one group open, its neighbour
// closed, and the numbers still add up.
UTEST(flow, partial_expansion_is_ordinary) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    test_graph_t t = {0};
    build(&t, alloc);

    md_flow_cut_t cut = {0};
    md_flow_cut_init(&cut, alloc);
    md_flow_cut_collapse_all(&cut);
    md_flow_cut_set_expanded(&cut, t.A, true);      // A open, B closed, in the same column
    md_flow_cut_resolve(&cut, &t.graph);

    EXPECT_EQ(md_flow_cut_validate(&cut, &t.graph, EPS), MD_FLOW_OK);

    bool saw_a0 = false, saw_a1 = false, saw_A = false, saw_B = false;
    for (size_t i = 0; i < md_array_size(cut.visible); ++i) {
        const uint32_t n = cut.visible[i];
        if (n == t.a0) saw_a0 = true;
        if (n == t.a1) saw_a1 = true;
        if (n == t.A)  saw_A  = true;
        if (n == t.B)  saw_B  = true;
    }
    EXPECT_TRUE(saw_a0 && saw_a1);
    EXPECT_FALSE(saw_A);    // an expanded node is replaced by its children, not drawn behind them
    EXPECT_TRUE(saw_B);

    md_flow_layout_t layout = {0};
    md_flow_layout_init(&layout, alloc);
    md_flow_layout_compute(&layout, &t.graph, &cut, NULL);
    for (size_t i = 0; i < md_array_size(layout.ribbons); ++i) {
        const md_flow_layout_ribbon_t* r = layout.ribbons + i;
        EXPECT_NEAR((double)(r->p1.y - r->p0.y), (double)(r->q1.y - r->q0.y), 1.0e-5);
    }

    md_flow_layout_free(&layout);
    md_flow_cut_free(&cut);
    md_arena_allocator_destroy(alloc);
}

// A band brackets an expanded ancestor. Its extent is taken from where the descendants actually
// landed, so it cannot disagree with them - and because siblings are kept contiguous, the band is
// a solid run rather than something with holes in it.
UTEST(flow, bands_bracket_their_expanded_ancestors) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    test_graph_t t = {0};
    build(&t, alloc);

    md_flow_cut_t cut = {0};
    md_flow_cut_init(&cut, alloc);
    md_flow_cut_collapse_all(&cut);
    md_flow_cut_set_expanded(&cut, t.A, true);
    md_flow_cut_resolve(&cut, &t.graph);

    md_flow_layout_t layout = {0};
    md_flow_layout_init(&layout, alloc);
    md_flow_layout_compute(&layout, &t.graph, &cut, NULL);

    // Exactly one band: A is open, B is closed and is drawn as a node instead.
    ASSERT_EQ(md_array_size(layout.bands), 1u);
    EXPECT_EQ(layout.bands[0].node, t.A);
    EXPECT_EQ(layout.bands[0].depth, 0u);

    const md_flow_layout_node_t *a0 = NULL, *a1 = NULL;
    for (size_t i = 0; i < md_array_size(layout.nodes); ++i) {
        if (layout.nodes[i].cut_idx == t.a0) a0 = layout.nodes + i;
        if (layout.nodes[i].cut_idx == t.a1) a1 = layout.nodes + i;
    }
    ASSERT_TRUE(a0 && a1);

    const float lo = a0->min.y < a1->min.y ? a0->min.y : a1->min.y;
    const float hi = a0->max.y > a1->max.y ? a0->max.y : a1->max.y;
    EXPECT_NEAR((double)layout.bands[0].min.y, (double)lo, 1.0e-5);
    EXPECT_NEAR((double)layout.bands[0].max.y, (double)hi, 1.0e-5);

    // Collapsing A leaves nothing to bracket: a group drawn as a node is not also a band.
    md_flow_cut_collapse_all(&cut);
    md_flow_cut_resolve(&cut, &t.graph);
    md_flow_layout_compute(&layout, &t.graph, &cut, NULL);
    EXPECT_EQ(md_array_size(layout.bands), 0u);

    md_flow_layout_free(&layout);
    md_flow_cut_free(&cut);
    md_arena_allocator_destroy(alloc);
}

// A container - a node whose children have children, a molecule of fragments - keeps its band when it is
// collapsed, spanning its own row and flagged so, at the depth it has open: the level stays where it is. A
// node of leaves (a fragment of atoms) has none while collapsed.
UTEST(flow, a_collapsed_container_keeps_its_band) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    md_flow_graph_t g = {0};
    md_flow_graph_init(&g, 2, alloc);
    const uint32_t M  = add(&g, 0, MD_FLOW_INVALID_INDEX, 0, 1.0f, 1);
    const uint32_t F0 = add(&g, 0, M,  1, 0.6f, 2);
    const uint32_t a0 = add(&g, 0, F0, 2, 0.6f, 3);
    const uint32_t F1 = add(&g, 0, M,  1, 0.4f, 4);
    const uint32_t a1 = add(&g, 0, F1, 2, 0.4f, 5);
    const uint32_t P  = add(&g, 1, MD_FLOW_INVALID_INDEX, 0, 1.0f, 6);
    md_flow_graph_add_link(&g, a0, P, 0.6f);
    md_flow_graph_add_link(&g, a1, P, 0.4f);

    md_flow_cut_t cut = {0};
    md_flow_cut_init(&cut, alloc);
    md_flow_cut_collapse_all(&cut);
    md_flow_cut_resolve(&cut, &g);
    md_flow_layout_t layout = {0};
    md_flow_layout_init(&layout, alloc);
    md_flow_layout_compute(&layout, &g, &cut, NULL);

    ASSERT_EQ(md_array_size(layout.bands), 1u);
    EXPECT_EQ(layout.bands[0].node, M);
    EXPECT_EQ(layout.bands[0].depth, 0u);
    EXPECT_EQ(layout.bands[0].flags, (uint32_t)MD_FLOW_BAND_FLAG_COLLAPSED);
    const md_flow_layout_node_t* m = NULL;
    for (size_t i = 0; i < md_array_size(layout.nodes); ++i) if (layout.nodes[i].cut_idx == M) m = layout.nodes + i;
    ASSERT_TRUE(m);
    EXPECT_NEAR((double)layout.bands[0].min.y, (double)m->min.y, 1.0e-5);
    EXPECT_NEAR((double)layout.bands[0].max.y, (double)m->max.y, 1.0e-5);

    // Open: the band brackets the fragments, which are of leaves and have none of their own
    md_flow_cut_set_expanded(&cut, M, true);
    md_flow_cut_resolve(&cut, &g);
    md_flow_layout_compute(&layout, &g, &cut, NULL);
    ASSERT_EQ(md_array_size(layout.bands), 1u);
    EXPECT_EQ(layout.bands[0].node, M);
    EXPECT_EQ(layout.bands[0].flags, 0u);

    md_flow_layout_free(&layout);
    md_flow_cut_free(&cut);
    md_arena_allocator_destroy(alloc);
}

// Ordering by contribution has to hold at EVERY level of expansion, not just the top one: within
// a column, families are ranked by the ancestor's weight and members by their own.
UTEST(flow, weight_order_ranks_within_families) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    test_graph_t t = {0};
    build(&t, alloc);

    md_flow_cut_t cut = {0};
    md_flow_cut_init(&cut, alloc);
    md_flow_cut_expand_all(&cut, &t.graph);
    md_flow_cut_resolve(&cut, &t.graph);

    md_flow_layout_params_t params = md_flow_layout_params_default();
    params.order = MD_FLOW_ORDER_WEIGHT;

    md_flow_layout_t layout = {0};
    md_flow_layout_init(&layout, alloc);
    md_flow_layout_compute(&layout, &t.graph, &cut, &params);

    // Column 0 fully expanded: A{a0 .3, a1 .2} and B{b0 .1, b1 .4}. B outranks A only if you rank
    // by the members, which is wrong - both groups weigh .5, so the tie is broken by member
    // weight inside each family. Whatever the family order, a0 must precede a1 and b1 precede b0.
    float y_a0 = -1, y_a1 = -1, y_b0 = -1, y_b1 = -1;
    for (size_t i = 0; i < md_array_size(layout.nodes); ++i) {
        const uint32_t n = layout.nodes[i].cut_idx;
        if (n == t.a0) y_a0 = layout.nodes[i].min.y;
        if (n == t.a1) y_a1 = layout.nodes[i].min.y;
        if (n == t.b0) y_b0 = layout.nodes[i].min.y;
        if (n == t.b1) y_b1 = layout.nodes[i].min.y;
    }
    ASSERT_TRUE(y_a0 >= 0 && y_a1 >= 0 && y_b0 >= 0 && y_b1 >= 0);
    EXPECT_LT(y_a0, y_a1);      // .3 above .2
    EXPECT_LT(y_b1, y_b0);      // .4 above .1

    // The middle column is flat, so it is a plain ranking: o0 .5, o1 .3, o2 .2.
    float y_o0 = -1, y_o1 = -1, y_o2 = -1;
    for (size_t i = 0; i < md_array_size(layout.nodes); ++i) {
        const uint32_t n = layout.nodes[i].cut_idx;
        if (n == t.o0) y_o0 = layout.nodes[i].min.y;
        if (n == t.o1) y_o1 = layout.nodes[i].min.y;
        if (n == t.o2) y_o2 = layout.nodes[i].min.y;
    }
    EXPECT_LT(y_o0, y_o1);
    EXPECT_LT(y_o1, y_o2);

    md_flow_layout_free(&layout);
    md_flow_cut_free(&cut);
    md_arena_allocator_destroy(alloc);
}

// A long tail of small nodes must not let "Others" swallow a quarter of the column. Folding runs
// smallest-first and stops at the cap, so what stays hidden is the genuinely negligible part.
UTEST(flow, others_never_exceeds_its_cap) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));

    // One dominant node and forty small ones, the shape a delocalized virtual column takes.
    md_flow_graph_t g = {0};
    md_flow_graph_init(&g, 2, alloc);

    const uint32_t num_small = 40;
    const float small_w = 0.6f / (float)num_small;      // 1.5% each - all under a 2% threshold
    uint32_t src[64] = {0};

    src[0] = add(&g, 0, MD_FLOW_INVALID_INDEX, 0, 0.4f, 1);
    for (uint32_t i = 0; i < num_small; ++i) {
        src[i + 1] = add(&g, 0, MD_FLOW_INVALID_INDEX, 0, small_w, 100 + i);
    }
    const uint32_t dst = add(&g, 1, MD_FLOW_INVALID_INDEX, 0, 1.0f, 999);
    for (uint32_t i = 0; i <= num_small; ++i) {
        md_flow_graph_add_link(&g, src[i], dst, g.nodes[src[i]].weight);
    }
    ASSERT_EQ(md_flow_graph_validate(&g, EPS), MD_FLOW_OK);

    md_flow_cut_t cut = {0};
    md_flow_cut_init(&cut, alloc);
    cut.threshold = 0.02f;      // every small node qualifies for folding

    // Uncapped, the whole 60% tail disappears into one grey box.
    cut.other_max = 1.0f;
    md_flow_cut_resolve(&cut, &g);
    ASSERT_EQ(md_array_size(cut.other), 1u);
    EXPECT_NEAR((double)cut.other[0].weight, 0.60, 0.01);
    EXPECT_EQ(md_flow_cut_validate(&cut, &g, EPS), MD_FLOW_OK);

    // Capped at 10%, folding stops there and the rest are drawn individually.
    cut.other_max = 0.10f;
    md_flow_cut_resolve(&cut, &g);
    ASSERT_EQ(md_array_size(cut.other), 1u);
    EXPECT_TRUE(cut.other[0].weight <= 0.10f + 1.0e-5f);
    EXPECT_TRUE(cut.other_count[0] > 1u);
    // Conservation is what makes the cap safe rather than a second way to lose weight.
    EXPECT_EQ(md_flow_cut_validate(&cut, &g, EPS), MD_FLOW_OK);

    // And the survivors are on screen rather than merely uncounted.
    EXPECT_EQ(md_array_size(cut.visible), (size_t)(1 + num_small - cut.other_count[0] + 1 /*Others*/) + 1 /*col 1*/);

    md_flow_cut_free(&cut);
    md_arena_allocator_destroy(alloc);
}

// Folding exactly one node would replace a labelled row with a grey box saying the same thing.
UTEST(flow, others_is_not_created_for_a_single_node) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    test_graph_t t = {0};
    build(&t, alloc);

    md_flow_cut_t cut = {0};
    md_flow_cut_init(&cut, alloc);
    md_flow_cut_expand_all(&cut, &t.graph);
    cut.other_max = 1.0f;
    cut.threshold = 0.15f;      // in column 0 only b0 (.10) falls below this
    md_flow_cut_resolve(&cut, &t.graph);

    for (size_t i = 0; i < md_array_size(cut.other); ++i) {
        EXPECT_TRUE(cut.other_count[i] > 1u);
    }
    EXPECT_EQ(md_flow_cut_validate(&cut, &t.graph, EPS), MD_FLOW_OK);

    md_flow_cut_free(&cut);
    md_arena_allocator_destroy(alloc);
}

// Compressing a column changes only how its band is DIVIDED, never what the numbers mean. The
// stubs must still exactly fill each node, at both ends, or ribbons overflow the bars they belong
// to - which is the whole reason the proportional layout was safe in the first place.
UTEST(flow, compressed_column_still_fills_its_nodes) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    test_graph_t t = {0};
    build(&t, alloc);

    md_flow_cut_t cut = {0};
    md_flow_cut_init(&cut, alloc);
    md_flow_cut_collapse_all(&cut);
    md_flow_cut_resolve(&cut, &t.graph);

    md_flow_layout_params_t params = md_flow_layout_params_default();
    params.weight_exponent[1] = 0.25f;      // squash the middle column only

    md_flow_layout_t layout = {0};
    md_flow_layout_init(&layout, alloc);
    md_flow_layout_compute(&layout, &t.graph, &cut, &params);

    for (size_t i = 0; i < md_array_size(layout.ribbons); ++i) {
        const md_flow_layout_ribbon_t* r = layout.ribbons + i;
        const md_flow_link_t* link = cut.links + r->link_idx;

        const md_flow_layout_node_t *src = NULL, *dst = NULL;
        for (size_t j = 0; j < md_array_size(layout.nodes); ++j) {
            if (layout.nodes[j].cut_idx == link->src) src = layout.nodes + j;
            if (layout.nodes[j].cut_idx == link->dst) dst = layout.nodes + j;
        }
        ASSERT_TRUE(src && dst);
        EXPECT_TRUE(r->p0.y >= src->min.y - 1.0e-5f && r->p1.y <= src->max.y + 1.0e-5f);
        EXPECT_TRUE(r->q0.y >= dst->min.y - 1.0e-5f && r->q1.y <= dst->max.y + 1.0e-5f);
    }

    // The middle column's rows are far closer in size than their weights are (.5 / .3 / .2), which
    // is the point: a channel you cannot see is a channel you cannot pick.
    float lo = 1.0e30f, hi = 0.0f;
    for (size_t i = 0; i < md_array_size(layout.nodes); ++i) {
        const md_flow_node_t* n = md_flow_cut_node(&cut, &t.graph, layout.nodes[i].cut_idx);
        if (n->column != 1) continue;
        const float h = layout.nodes[i].max.y - layout.nodes[i].min.y;
        lo = h < lo ? h : lo;
        hi = h > hi ? h : hi;
    }
    EXPECT_TRUE(hi / lo < 2.0f);        // weights differ 2.5x; drawn heights much less

    md_flow_layout_free(&layout);
    md_flow_cut_free(&cut);
    md_arena_allocator_destroy(alloc);
}

// With the cap off, the threshold means exactly what it says: nothing below it is drawn. Opening
// the column's "Others" is how you look anyway - which is what lets the threshold stay honest
// instead of being quietly overridden.
UTEST(flow, threshold_is_honest_and_others_opens) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));

    md_flow_graph_t g = {0};
    md_flow_graph_init(&g, 2, alloc);

    const uint32_t num_small = 40;
    const float small_w = 0.6f / (float)num_small;      // 1.5% each
    uint32_t src[64] = {0};
    src[0] = add(&g, 0, MD_FLOW_INVALID_INDEX, 0, 0.4f, 1);
    for (uint32_t i = 0; i < num_small; ++i) src[i + 1] = add(&g, 0, MD_FLOW_INVALID_INDEX, 0, small_w, 100 + i);
    const uint32_t dst = add(&g, 1, MD_FLOW_INVALID_INDEX, 0, 1.0f, 999);
    for (uint32_t i = 0; i <= num_small; ++i) md_flow_graph_add_link(&g, src[i], dst, g.nodes[src[i]].weight);

    md_flow_cut_t cut = {0};
    md_flow_cut_init(&cut, alloc);
    cut.threshold = 0.02f;
    cut.other_max = 1.0f;           // uncapped: the honest default
    md_flow_cut_resolve(&cut, &g);

    // Column 0 shows the dominant node plus one "Others". Nothing sub-threshold leaks through.
    size_t col0 = 0;
    for (size_t i = 0; i < md_array_size(cut.visible); ++i) {
        const md_flow_node_t* n = md_flow_cut_node(&cut, &g, cut.visible[i]);
        if (n->column != 0) continue;
        col0++;
        EXPECT_TRUE(n->weight >= cut.threshold || (n->flags & MD_FLOW_NODE_FLAG_OTHER));
    }
    EXPECT_EQ(col0, 2u);

    // Opening it draws every member individually, and the diagram still adds up.
    md_flow_cut_set_other_expanded(&cut, 0, true);
    md_flow_cut_resolve(&cut, &g);
    EXPECT_EQ(md_array_size(cut.other), 0u);

    size_t col0_open = 0;
    for (size_t i = 0; i < md_array_size(cut.visible); ++i) {
        if (md_flow_cut_node(&cut, &g, cut.visible[i])->column == 0) col0_open++;
    }
    EXPECT_EQ(col0_open, (size_t)(num_small + 1));
    EXPECT_EQ(md_flow_cut_validate(&cut, &g, EPS), MD_FLOW_OK);

    // And closing it puts them back.
    md_flow_cut_set_other_expanded(&cut, 0, false);
    md_flow_cut_resolve(&cut, &g);
    EXPECT_EQ(md_array_size(cut.other), 1u);

    md_flow_cut_free(&cut);
    md_arena_allocator_destroy(alloc);
}

// Two mirrored columns of the same three fragments, added in the same order: GRAPH order keeps that
// order on both sides, and a flow from a fragment to itself is reported as a self flow - including
// after a collapse merges two leaves of one group, whose flow between them becomes the group's own.
UTEST(flow, graph_order_mirrors_and_marks_self_flow) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    md_flow_graph_t g = {0};
    md_flow_graph_init(&g, 2, alloc);

    // hole: X{x0 .3, x1 .2} Y .3 Z .2 ; particle: X{x0 .1, x1 .1} Y .1 Z .7
    uint32_t n[2][5];
    const float w[2][5] = { { 0.5f, 0.3f, 0.2f, 0.3f, 0.2f }, { 0.2f, 0.1f, 0.1f, 0.1f, 0.7f } };
    for (uint32_t c = 0; c < 2; ++c) {
        md_flow_node_t X  = { .column = c, .parent = MD_FLOW_INVALID_INDEX, .weight = w[c][0], .key = 10 + c, .entity = 1 };
        n[c][0] = md_flow_graph_add_node(&g, &X);
        md_flow_node_t x0 = { .column = c, .parent = n[c][0], .level = 1, .weight = w[c][1], .key = 20 + c, .entity = 11 };
        n[c][1] = md_flow_graph_add_node(&g, &x0);
        md_flow_node_t x1 = { .column = c, .parent = n[c][0], .level = 1, .weight = w[c][2], .key = 30 + c, .entity = 12 };
        n[c][2] = md_flow_graph_add_node(&g, &x1);
        md_flow_node_t Y  = { .column = c, .parent = MD_FLOW_INVALID_INDEX, .weight = w[c][3], .key = 40 + c, .entity = 2 };
        n[c][3] = md_flow_graph_add_node(&g, &Y);
        md_flow_node_t Z  = { .column = c, .parent = MD_FLOW_INVALID_INDEX, .weight = w[c][4], .key = 50 + c, .entity = 3 };
        n[c][4] = md_flow_graph_add_node(&g, &Z);
    }
    md_flow_graph_add_link(&g, n[0][1], n[1][2], 0.1f);   // x0 -> x1: inside X
    md_flow_graph_add_link(&g, n[0][1], n[1][4], 0.2f);   // x0 -> Z
    md_flow_graph_add_link(&g, n[0][2], n[1][1], 0.1f);   // x1 -> x0: inside X
    md_flow_graph_add_link(&g, n[0][2], n[1][4], 0.1f);   // x1 -> Z
    md_flow_graph_add_link(&g, n[0][3], n[1][3], 0.1f);   // Y -> Y
    md_flow_graph_add_link(&g, n[0][3], n[1][4], 0.2f);   // Y -> Z
    md_flow_graph_add_link(&g, n[0][4], n[1][4], 0.2f);   // Z -> Z
    ASSERT_EQ(md_flow_graph_validate(&g, 1.0e-5f), MD_FLOW_OK);

    md_flow_cut_t cut = {0};
    md_flow_cut_init(&cut, alloc);
    md_flow_cut_resolve(&cut, &g);   // X collapsed

    md_flow_layout_params_t params = md_flow_layout_params_default();
    params.order = MD_FLOW_ORDER_GRAPH;
    md_flow_layout_t layout = {0};
    md_flow_layout_init(&layout, alloc);
    md_flow_layout_compute(&layout, &g, &cut, &params);

    // X, Y, Z top to bottom in both columns, despite Z being the heaviest on the right
    for (uint32_t c = 0; c < 2; ++c) {
        float y[3] = {0};
        const uint32_t want[3] = { n[c][0], n[c][3], n[c][4] };
        for (size_t i = 0; i < md_array_size(layout.nodes); ++i) {
            for (int k = 0; k < 3; ++k) if (layout.nodes[i].cut_idx == want[k]) y[k] = layout.nodes[i].min.y;
        }
        EXPECT_LT(y[0], y[1]);
        EXPECT_LT(y[1], y[2]);
    }

    // self flows: X -> X (0.2, the two cross links merged), Y -> Y, Z -> Z; the rest is transfer
    double self = 0.0, transfer = 0.0;
    for (size_t i = 0; i < md_array_size(layout.ribbons); ++i) {
        const md_flow_link_t* l = cut.links + layout.ribbons[i].link_idx;
        if (layout.ribbons[i].flags & MD_FLOW_RIBBON_FLAG_SELF) {
            self += l->weight;
            EXPECT_EQ(md_flow_cut_node(&cut, &g, l->src)->entity, md_flow_cut_node(&cut, &g, l->dst)->entity);
        } else {
            transfer += l->weight;
        }
    }
    EXPECT_NEAR(self, 0.5, 1.0e-6);
    EXPECT_NEAR(transfer, 0.5, 1.0e-6);

    // expanded, x0 -> x1 is a transfer between atoms, not a self flow
    md_flow_cut_set_expanded(&cut, n[0][0], true);
    md_flow_cut_set_expanded(&cut, n[1][0], true);
    md_flow_cut_resolve(&cut, &g);
    md_flow_layout_compute(&layout, &g, &cut, &params);
    self = 0.0;
    for (size_t i = 0; i < md_array_size(layout.ribbons); ++i) {
        if (layout.ribbons[i].flags & MD_FLOW_RIBBON_FLAG_SELF) self += cut.links[layout.ribbons[i].link_idx].weight;
    }
    EXPECT_NEAR(self, 0.3, 1.0e-6);

    md_flow_layout_free(&layout);
    md_flow_cut_free(&cut);
    md_arena_allocator_destroy(alloc);
}

// Three groups of 4, 3 and 2 atoms, mirrored in two columns, every atom flowing to itself. 'hole' and
// 'part' are the atom weights of the two columns.
typedef struct mirror_graph_t {
    md_flow_graph_t graph;
    uint32_t group[2][3];
    uint32_t atom[2][9];
} mirror_graph_t;

static void build_mirror(mirror_graph_t* m, md_allocator_i* alloc, const float hole[9], const float part[9]) {
    md_flow_graph_t* g = &m->graph;
    md_flow_graph_init(g, 2, alloc);
    const uint32_t size[3] = { 4, 3, 2 };
    for (uint32_t c = 0; c < 2; ++c) {
        const float* w = c == 0 ? hole : part;
        uint32_t a = 0;
        for (uint32_t k = 0; k < 3; ++k) {
            float gw = 0.0f;
            for (uint32_t j = 0; j < size[k]; ++j) gw += w[a + j];
            md_flow_node_t gn = { .column = c, .parent = MD_FLOW_INVALID_INDEX, .weight = gw, .key = 100 * c + k, .entity = 1000 + k };
            m->group[c][k] = md_flow_graph_add_node(g, &gn);
            for (uint32_t j = 0; j < size[k]; ++j, ++a) {
                md_flow_node_t an = { .column = c, .parent = m->group[c][k], .level = 1, .weight = w[a], .key = 100 * c + 10 + a, .entity = 2000 + a };
                m->atom[c][a] = md_flow_graph_add_node(g, &an);
            }
        }
    }
}

static uint32_t visible_in_column(const md_flow_cut_t* cut, const md_flow_graph_t* g, uint32_t column) {
    uint32_t n = 0;
    for (size_t i = 0; i < md_array_size(cut->visible); ++i) n += md_flow_cut_node(cut, g, cut->visible[i])->column == column ? 1 : 0;
    return n;
}

UTEST(flow, budget_opens_the_heaviest_within_the_rows) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    const float w[9] = { 0.15f, 0.15f, 0.15f, 0.15f, 0.1f, 0.1f, 0.1f, 0.05f, 0.05f };
    mirror_graph_t m = {0};
    build_mirror(&m, alloc, w, w);
    for (uint32_t a = 0; a < 9; ++a) md_flow_graph_add_link(&m.graph, m.atom[0][a], m.atom[1][a], w[a]);
    ASSERT_EQ(md_flow_graph_validate(&m.graph, EPS), MD_FLOW_OK);

    md_flow_cut_t cut = {0};
    md_flow_cut_init(&cut, alloc);
    cut.mirror = true;

    // No budget: nothing opens by itself
    md_flow_cut_resolve(&cut, &m.graph);
    EXPECT_EQ(visible_in_column(&cut, &m.graph, 0), 3u);

    // 6 rows: the heaviest group opens (3 + 3), on both sides; neither of the others fits after it
    cut.max_rows = 6;
    md_flow_cut_resolve(&cut, &m.graph);
    for (uint32_t c = 0; c < 2; ++c) {
        EXPECT_TRUE (md_flow_cut_is_open(&cut, m.group[c][0]));
        EXPECT_FALSE(md_flow_cut_is_open(&cut, m.group[c][1]));
        EXPECT_FALSE(md_flow_cut_is_open(&cut, m.group[c][2]));
        EXPECT_EQ(visible_in_column(&cut, &m.graph, c), 6u);
    }
    EXPECT_EQ(md_flow_cut_validate(&cut, &m.graph, EPS), MD_FLOW_OK);

    // The user closes the heaviest: the budget goes to the next ones, and the pin holds on both sides
    md_flow_cut_pin(&cut, &m.graph, m.group[0][0], false);
    md_flow_cut_resolve(&cut, &m.graph);
    for (uint32_t c = 0; c < 2; ++c) {
        EXPECT_FALSE(md_flow_cut_is_open(&cut, m.group[c][0]));
        EXPECT_TRUE (md_flow_cut_is_open(&cut, m.group[c][1]));
        EXPECT_TRUE (md_flow_cut_is_open(&cut, m.group[c][2]));
        EXPECT_EQ(visible_in_column(&cut, &m.graph, c), 6u);
    }

    // A pin open wins over the budget
    md_flow_cut_clear_pins(&cut);
    cut.max_rows = 3;
    md_flow_cut_pin(&cut, &m.graph, m.group[1][2], true);
    md_flow_cut_resolve(&cut, &m.graph);
    for (uint32_t c = 0; c < 2; ++c) {
        EXPECT_TRUE (md_flow_cut_is_open(&cut, m.group[c][2]));
        EXPECT_FALSE(md_flow_cut_is_open(&cut, m.group[c][0]));
    }
    EXPECT_EQ(md_flow_cut_validate(&cut, &m.graph, EPS), MD_FLOW_OK);

    md_flow_cut_free(&cut);
    md_arena_allocator_destroy(alloc);
}

UTEST(flow, others_fold_per_parent_and_mirrored) {
    md_allocator_i* alloc = md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(4));
    // atom 1 is small as a hole and large as a particle: mirrored, it must not fold on either side
    const float hole[9] = { 0.20f, 0.04f, 0.04f, 0.04f, 0.2f, 0.2f, 0.1f, 0.1f, 0.08f };
    const float part[9] = { 0.04f, 0.20f, 0.04f, 0.04f, 0.2f, 0.2f, 0.1f, 0.1f, 0.08f };
    mirror_graph_t m = {0};
    build_mirror(&m, alloc, hole, part);
    // atom 0 -> atom 1 carries 0.16, the rest goes to itself
    md_flow_graph_add_link(&m.graph, m.atom[0][0], m.atom[1][0], 0.04f);
    md_flow_graph_add_link(&m.graph, m.atom[0][0], m.atom[1][1], 0.16f);
    md_flow_graph_add_link(&m.graph, m.atom[0][1], m.atom[1][1], 0.04f);
    for (uint32_t a = 2; a < 9; ++a) md_flow_graph_add_link(&m.graph, m.atom[0][a], m.atom[1][a], hole[a]);
    ASSERT_EQ(md_flow_graph_validate(&m.graph, EPS), MD_FLOW_OK);

    md_flow_cut_t cut = {0};
    md_flow_cut_init(&cut, alloc);
    cut.mirror = true;
    cut.threshold = 0.05f;
    md_flow_cut_pin(&cut, &m.graph, m.group[0][0], true);
    md_flow_cut_resolve(&cut, &m.graph);

    // In each column: atoms 2 and 3 fold into an Others inside group 0. Atom 1 stays (0.2 as a particle).
    ASSERT_EQ(md_array_size(cut.other), (size_t)2);
    for (size_t i = 0; i < 2; ++i) {
        EXPECT_EQ(cut.other[i].parent, m.group[cut.other[i].column][0]);
        EXPECT_EQ(cut.other_count[i], 2u);
    }
    for (uint32_t c = 0; c < 2; ++c) {
        bool atom1 = false;
        for (size_t v = 0; v < md_array_size(cut.visible); ++v) atom1 |= cut.visible[v] == m.atom[c][1];
        EXPECT_TRUE(atom1);
    }
    EXPECT_EQ(md_flow_cut_validate(&cut, &m.graph, EPS), MD_FLOW_OK);

    // It sits inside its parent's slot: right after the last visible sibling
    for (size_t v = 0; v + 1 < md_array_size(cut.visible); ++v) {
        const md_flow_node_t* n = md_flow_cut_node(&cut, &m.graph, cut.visible[v + 1]);
        if (n->flags & MD_FLOW_NODE_FLAG_OTHER) {
            EXPECT_EQ(md_flow_cut_node(&cut, &m.graph, cut.visible[v])->parent, n->parent);
        }
    }

    // Opening one opens its mirror image too
    md_flow_cut_open_other(&cut, &m.graph, cut.other + 0, true);
    md_flow_cut_resolve(&cut, &m.graph);
    EXPECT_EQ(md_array_size(cut.other), (size_t)0);
    EXPECT_TRUE(md_flow_cut_is_other_open_for(&cut, m.graph.nodes + m.atom[1][2]));
    EXPECT_EQ(md_flow_cut_validate(&cut, &m.graph, EPS), MD_FLOW_OK);

    md_flow_cut_free(&cut);
    md_arena_allocator_destroy(alloc);
}
