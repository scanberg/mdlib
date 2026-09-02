#include <md_flow.h>

#include <core/md_allocator.h>
#include <core/md_common.h>
#include <core/md_log.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

// Key namespace for the synthesized "Others" nodes, kept clear of any caller's key scheme so the
// cut can be keyed by node identity without colliding with the graph's own keys.
#define MD_FLOW_OTHER_KEY_TAG 0xF10A0000ULL

// ---------------------------------------------------------------------------
// Graph
// ---------------------------------------------------------------------------

void md_flow_graph_init(md_flow_graph_t* graph, uint32_t num_columns, struct md_allocator_i* alloc) {
    ASSERT(graph);
    ASSERT(alloc);
    MEMSET(graph, 0, sizeof(md_flow_graph_t));
    graph->num_columns = num_columns;
    graph->alloc = alloc;
}

void md_flow_graph_free(md_flow_graph_t* graph) {
    ASSERT(graph);
    if (graph->alloc) {
        md_array_free(graph->nodes, graph->alloc);
        md_array_free(graph->links, graph->alloc);
    }
    MEMSET(graph, 0, sizeof(md_flow_graph_t));
}

void md_flow_graph_clear(md_flow_graph_t* graph) {
    ASSERT(graph);
    md_array_shrink(graph->nodes, 0);
    md_array_shrink(graph->links, 0);
}

uint32_t md_flow_graph_add_node(md_flow_graph_t* graph, const md_flow_node_t* node) {
    ASSERT(graph);
    ASSERT(node);

    if (node->column >= graph->num_columns) {
        MD_LOG_ERROR("flow: node column %u is out of range (graph has %u columns)", node->column, graph->num_columns);
        return MD_FLOW_INVALID_INDEX;
    }

    // A parent must be an EARLIER node in the SAME column. Requiring 'earlier' is what makes the
    // hierarchy acyclic by construction rather than by a check.
    if (node->parent != MD_FLOW_INVALID_INDEX) {
        const uint32_t count = (uint32_t)md_array_size(graph->nodes);
        if (node->parent >= count) {
            MD_LOG_ERROR("flow: node parent %u is not an earlier node", node->parent);
            return MD_FLOW_INVALID_INDEX;
        }
        if (graph->nodes[node->parent].column != node->column) {
            MD_LOG_ERROR("flow: node parent %u lives in column %u, child in column %u", node->parent,
                graph->nodes[node->parent].column, node->column);
            return MD_FLOW_INVALID_INDEX;
        }
    }

    const uint32_t idx = (uint32_t)md_array_size(graph->nodes);
    md_array_push(graph->nodes, *node, graph->alloc);
    return idx;
}

bool md_flow_graph_add_link(md_flow_graph_t* graph, uint32_t src, uint32_t dst, float weight) {
    ASSERT(graph);

    if (!(weight > 0.0f)) {
        // Zero and NaN links carry no information and would only add degenerate ribbons.
        return false;
    }

    const uint32_t count = (uint32_t)md_array_size(graph->nodes);
    if (src >= count || dst >= count) {
        MD_LOG_ERROR("flow: link (%u -> %u) references a node that does not exist", src, dst);
        return false;
    }

    const md_flow_link_t link = { .src = src, .dst = dst, .weight = weight };
    md_array_push(graph->links, link, graph->alloc);
    return true;
}

bool md_flow_node_is_leaf(const md_flow_graph_t* graph, uint32_t node_idx) {
    ASSERT(graph);
    const size_t count = md_array_size(graph->nodes);
    for (size_t i = node_idx + 1; i < count; ++i) {
        // Children always come after their parent, so the scan can start just past the node.
        if (graph->nodes[i].parent == node_idx) {
            return false;
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------

const char* md_flow_result_str(md_flow_result_t result) {
    switch (result) {
    case MD_FLOW_OK:                return "ok";
    case MD_FLOW_ERR_HIERARCHY:     return "a node's parent is not an earlier node in the same column";
    case MD_FLOW_ERR_LAYERING:      return "a link does not connect two leaves in adjacent columns";
    case MD_FLOW_ERR_COLUMN_SUM:    return "the leaf weights of a column do not sum to one";
    case MD_FLOW_ERR_BAND_SUM:      return "the link weights across a band do not sum to one";
    case MD_FLOW_ERR_NODE_FLOW:     return "a node's links do not sum to its own weight";
    case MD_FLOW_ERR_PARTITION:     return "a node's children do not partition its weight";
    default:                        return "unknown";
    }
}

// Marks, for every node, whether anything claims it as a parent. One pass instead of the
// quadratic md_flow_node_is_leaf() when the answer is needed for the whole graph.
static void flow_mark_internal(bool* out_has_children, const md_flow_graph_t* graph) {
    const size_t count = md_array_size(graph->nodes);
    MEMSET(out_has_children, 0, sizeof(bool) * count);
    for (size_t i = 0; i < count; ++i) {
        const uint32_t parent = graph->nodes[i].parent;
        if (parent != MD_FLOW_INVALID_INDEX) {
            out_has_children[parent] = true;
        }
    }
}

md_flow_result_t md_flow_graph_validate(const md_flow_graph_t* graph, float eps) {
    ASSERT(graph);

    const size_t num_nodes = md_array_size(graph->nodes);
    const size_t num_links = md_array_size(graph->links);

    if (num_nodes == 0) {
        return MD_FLOW_OK;
    }

    md_temp_scope_t temp = md_temp_begin();

    bool*  has_children  = md_temp_alloc_array(temp, bool,  num_nodes);
    double* child_sum    = md_temp_alloc_zero_array(temp, double, num_nodes);
    double* out_sum      = md_temp_alloc_zero_array(temp, double, num_nodes);
    double* in_sum       = md_temp_alloc_zero_array(temp, double, num_nodes);
    double* column_sum   = md_temp_alloc_zero_array(temp, double, graph->num_columns);
    double* band_sum     = md_temp_alloc_zero_array(temp, double, graph->num_columns);

    flow_mark_internal(has_children, graph);

    md_flow_result_t result = MD_FLOW_OK;

    // (hierarchy) + accumulate the child sums and the per-column leaf sums
    for (size_t i = 0; i < num_nodes && result == MD_FLOW_OK; ++i) {
        const md_flow_node_t* node = graph->nodes + i;

        if (node->column >= graph->num_columns) {
            result = MD_FLOW_ERR_HIERARCHY;
            break;
        }
        if (node->parent != MD_FLOW_INVALID_INDEX) {
            if (node->parent >= i || graph->nodes[node->parent].column != node->column) {
                result = MD_FLOW_ERR_HIERARCHY;
                break;
            }
            child_sum[node->parent] += node->weight;
        }
        if (!has_children[i]) {
            column_sum[node->column] += node->weight;
        }
    }

    // (1) layering - links go c -> c+1, and only between leaves, because everything coarser is
    // derived by aggregation and a link on an internal node would be counted twice.
    for (size_t i = 0; i < num_links && result == MD_FLOW_OK; ++i) {
        const md_flow_link_t* link = graph->links + i;
        const uint32_t src_col = graph->nodes[link->src].column;
        const uint32_t dst_col = graph->nodes[link->dst].column;
        if (src_col + 1 != dst_col || has_children[link->src] || has_children[link->dst]) {
            result = MD_FLOW_ERR_LAYERING;
            break;
        }
        band_sum[src_col] += link->weight;
        out_sum[link->src] += link->weight;
        in_sum[link->dst]  += link->weight;
    }

    // (3) hierarchy partitions
    for (size_t i = 0; i < num_nodes && result == MD_FLOW_OK; ++i) {
        if (has_children[i] && fabs(child_sum[i] - (double)graph->nodes[i].weight) > (double)eps) {
            result = MD_FLOW_ERR_PARTITION;
        }
    }

    // A leaf's ribbons must exactly fill it. This is not decoration: the layout sub-allocates
    // ribbon stubs inside a node's own span, so a node whose links overrun its weight would draw
    // ribbons spilling out of the bar they belong to.
    for (size_t i = 0; i < num_nodes && result == MD_FLOW_OK; ++i) {
        if (has_children[i]) continue;
        const md_flow_node_t* node = graph->nodes + i;
        if (node->column + 1 < graph->num_columns && fabs(out_sum[i] - (double)node->weight) > (double)eps) {
            result = MD_FLOW_ERR_NODE_FLOW;
        }
        if (node->column > 0 && fabs(in_sum[i] - (double)node->weight) > (double)eps) {
            result = MD_FLOW_ERR_NODE_FLOW;
        }
    }

    // (2) conservation
    for (uint32_t c = 0; c < graph->num_columns && result == MD_FLOW_OK; ++c) {
        if (fabs(column_sum[c] - 1.0) > (double)eps) {
            result = MD_FLOW_ERR_COLUMN_SUM;
        }
    }
    for (uint32_t c = 0; c + 1 < graph->num_columns && result == MD_FLOW_OK; ++c) {
        if (fabs(band_sum[c] - 1.0) > (double)eps) {
            result = MD_FLOW_ERR_BAND_SUM;
        }
    }

    md_temp_end(temp);
    return result;
}

// ---------------------------------------------------------------------------
// Cut
// ---------------------------------------------------------------------------

void md_flow_cut_init(md_flow_cut_t* cut, struct md_allocator_i* alloc) {
    ASSERT(cut);
    ASSERT(alloc);
    MEMSET(cut, 0, sizeof(md_flow_cut_t));
    md_bitfield_init(&cut->expanded, alloc);
    cut->other_max = 1.0f;      // uncapped unless the caller says otherwise
    cut->alloc = alloc;
}

void md_flow_cut_free(md_flow_cut_t* cut) {
    ASSERT(cut);
    if (cut->alloc) {
        md_bitfield_free(&cut->expanded);
        md_array_free(cut->visible, cut->alloc);
        md_array_free(cut->links,   cut->alloc);
        md_array_free(cut->other,   cut->alloc);
        md_array_free(cut->other_count, cut->alloc);
    }
    MEMSET(cut, 0, sizeof(md_flow_cut_t));
}

void md_flow_cut_set_expanded(md_flow_cut_t* cut, uint32_t node_idx, bool expanded) {
    ASSERT(cut);
    if (expanded) {
        md_bitfield_set_bit(&cut->expanded, node_idx);
    } else {
        md_bitfield_clear_bit(&cut->expanded, node_idx);
    }
}

void md_flow_cut_set_other_expanded(md_flow_cut_t* cut, uint32_t column, bool expanded) {
    ASSERT(cut);
    if (column >= MD_FLOW_MAX_COLUMNS) return;
    if (expanded) cut->other_expanded |=  (1u << column);
    else          cut->other_expanded &= ~(1u << column);
}

bool md_flow_cut_is_other_expanded(const md_flow_cut_t* cut, uint32_t column) {
    ASSERT(cut);
    return column < MD_FLOW_MAX_COLUMNS && (cut->other_expanded & (1u << column)) != 0;
}

bool md_flow_cut_is_expanded(const md_flow_cut_t* cut, uint32_t node_idx) {
    ASSERT(cut);
    return md_bitfield_test_bit(&cut->expanded, node_idx);
}

void md_flow_cut_collapse_all(md_flow_cut_t* cut) {
    ASSERT(cut);
    md_bitfield_clear(&cut->expanded);
}

void md_flow_cut_expand_all(md_flow_cut_t* cut, const md_flow_graph_t* graph) {
    ASSERT(cut);
    ASSERT(graph);
    const size_t count = md_array_size(graph->nodes);
    for (size_t i = 0; i < count; ++i) {
        md_bitfield_set_bit(&cut->expanded, i);
    }
}

const md_flow_node_t* md_flow_cut_node(const md_flow_cut_t* cut, const md_flow_graph_t* graph, uint32_t cut_idx) {
    ASSERT(cut);
    ASSERT(graph);
    const uint32_t num_graph_nodes = (uint32_t)md_array_size(graph->nodes);
    if (cut_idx < num_graph_nodes) {
        return graph->nodes + cut_idx;
    }
    const uint32_t other_idx = cut_idx - num_graph_nodes;
    if (other_idx < (uint32_t)md_array_size(cut->other)) {
        return cut->other + other_idx;
    }
    return NULL;
}

static int flow_link_cmp(const void* a, const void* b) {
    const md_flow_link_t* la = (const md_flow_link_t*)a;
    const md_flow_link_t* lb = (const md_flow_link_t*)b;
    if (la->src != lb->src) return la->src < lb->src ? -1 : 1;
    if (la->dst != lb->dst) return la->dst < lb->dst ? -1 : 1;
    return 0;
}

void md_flow_cut_resolve(md_flow_cut_t* cut, const md_flow_graph_t* graph) {
    ASSERT(cut);
    ASSERT(graph);

    md_array_shrink(cut->visible, 0);
    md_array_shrink(cut->links,   0);
    md_array_shrink(cut->other,   0);
    md_array_shrink(cut->other_count, 0);

    const size_t num_nodes = md_array_size(graph->nodes);
    if (num_nodes == 0) {
        return;
    }

    md_temp_scope_t temp = md_temp_begin();

    bool* has_children = md_temp_alloc_array(temp, bool, num_nodes);
    flow_mark_internal(has_children, graph);

    // Pass 1: the antichain, before thresholding. A node is in it when it is either collapsed or
    // a leaf, and every one of its ancestors is expanded.
    bool* in_cut = md_temp_alloc_zero_array(temp, bool, num_nodes);
    for (size_t i = 0; i < num_nodes; ++i) {
        const md_flow_node_t* node = graph->nodes + i;

        bool ancestors_expanded = true;
        for (uint32_t p = node->parent; p != MD_FLOW_INVALID_INDEX; p = graph->nodes[p].parent) {
            if (!md_bitfield_test_bit(&cut->expanded, p)) {
                ancestors_expanded = false;
                break;
            }
        }
        if (!ancestors_expanded) continue;

        const bool expanded = has_children[i] && md_bitfield_test_bit(&cut->expanded, i);
        if (!expanded) {
            in_cut[i] = true;
        }
    }

    // Pass 2: threshold. Weights are already a share of their column, so the test is direct.
    // Everything folded in a column lands in that column's single synthesized "Others".
    const uint32_t num_graph_nodes = (uint32_t)num_nodes;
    uint32_t* other_of_column = md_temp_alloc_array(temp, uint32_t, graph->num_columns);
    double*   other_weight    = md_temp_alloc_zero_array(temp, double, graph->num_columns);
    for (uint32_t c = 0; c < graph->num_columns; ++c) {
        other_of_column[c] = MD_FLOW_INVALID_INDEX;
    }

    // rep[i] is the cut index that node i's flow is attributed to, for nodes in the antichain.
    uint32_t* rep = md_temp_alloc_array(temp, uint32_t, num_nodes);
    for (size_t i = 0; i < num_nodes; ++i) {
        rep[i] = MD_FLOW_INVALID_INDEX;
    }

    // Fold smallest-first and stop at the cap, so what ends up hidden is the genuinely negligible
    // tail rather than "whatever happened to be under the threshold" - which in a delocalized
    // column can be most of it.
    bool* folded = md_temp_alloc_zero_array(temp, bool, num_nodes);
    uint32_t* other_count = md_temp_alloc_zero_array(temp, uint32_t, graph->num_columns);

    if (cut->threshold > 0.0f) {
        uint32_t* candidate = md_temp_alloc_array(temp, uint32_t, num_nodes);
        const float cap = (cut->other_max > 0.0f) ? cut->other_max : 1.0f;

        for (uint32_t c = 0; c < graph->num_columns; ++c) {
            if (md_flow_cut_is_other_expanded(cut, c)) continue;   // opened: fold nothing here

            uint32_t count = 0;
            for (size_t i = 0; i < num_nodes; ++i) {
                if (!in_cut[i] || graph->nodes[i].column != c) continue;
                if (graph->nodes[i].weight < cut->threshold) candidate[count++] = (uint32_t)i;
            }
            if (count == 0) continue;

            // Ascending by weight. Insertion sort: these lists are short and already near-sorted
            // whenever the caller built the column in any sensible order.
            for (uint32_t a = 1; a < count; ++a) {
                const uint32_t v = candidate[a];
                uint32_t b = a;
                while (b > 0 && graph->nodes[candidate[b - 1]].weight > graph->nodes[v].weight) {
                    candidate[b] = candidate[b - 1];
                    b--;
                }
                candidate[b] = v;
            }

            double sum = 0.0;
            for (uint32_t a = 0; a < count; ++a) {
                const double w = (double)graph->nodes[candidate[a]].weight;
                if (sum + w > (double)cap) break;
                sum += w;
                folded[candidate[a]] = true;
                other_count[c] += 1;
            }

            // One survivor in an Others is just that node wearing a worse label.
            if (other_count[c] == 1) {
                for (size_t i = 0; i < num_nodes; ++i) {
                    if (folded[i] && graph->nodes[i].column == c) { folded[i] = false; break; }
                }
                other_count[c] = 0;
                sum = 0.0;
            }

            other_weight[c] = sum;
            if (other_count[c] > 0) other_of_column[c] = 0;  // marker: this column needs an Others
        }
    }

    for (size_t i = 0; i < num_nodes; ++i) {
        if (!in_cut[i] || folded[i]) continue;
        rep[i] = (uint32_t)i;
        md_array_push(cut->visible, (uint32_t)i, cut->alloc);
    }

    // Materialize the Others nodes now that their weights are known, then point the folded nodes
    // at them. They are appended after the graph's nodes in cut index space.
    for (uint32_t c = 0; c < graph->num_columns; ++c) {
        if (other_of_column[c] == MD_FLOW_INVALID_INDEX) continue;
        const md_flow_node_t other = {
            .column = c,
            .parent = MD_FLOW_INVALID_INDEX,
            .level  = 0,
            .weight = (float)other_weight[c],
            .color  = {0.6f, 0.6f, 0.6f, 1.0f},
            .label  = STR_LIT("Others"),
            .key    = (MD_FLOW_OTHER_KEY_TAG << 32) | (uint64_t)c,
            .flags  = MD_FLOW_NODE_FLAG_OTHER,
        };
        other_of_column[c] = num_graph_nodes + (uint32_t)md_array_size(cut->other);
        md_array_push(cut->other, other, cut->alloc);
        md_array_push(cut->other_count, other_count[c], cut->alloc);
        md_array_push(cut->visible, other_of_column[c], cut->alloc);
    }

    for (size_t i = 0; i < num_nodes; ++i) {
        if (in_cut[i] && folded[i]) {
            rep[i] = other_of_column[graph->nodes[i].column];
        }
    }

    // Pass 3: every leaf attributes to its visible ancestor-or-self, so the links can be
    // aggregated. This is the step that guarantees a collapsed view and an expanded view tell
    // the same story - both are sums over the same leaf links.
    uint32_t* leaf_rep = md_temp_alloc_array(temp, uint32_t, num_nodes);
    for (size_t i = 0; i < num_nodes; ++i) {
        uint32_t n = (uint32_t)i;
        while (n != MD_FLOW_INVALID_INDEX && !in_cut[n]) {
            n = graph->nodes[n].parent;
        }
        leaf_rep[i] = (n == MD_FLOW_INVALID_INDEX) ? MD_FLOW_INVALID_INDEX : rep[n];
    }

    const size_t num_links = md_array_size(graph->links);
    if (num_links > 0) {
        md_flow_link_t* scratch = md_temp_alloc_array(temp, md_flow_link_t, num_links);
        size_t count = 0;
        for (size_t i = 0; i < num_links; ++i) {
            const md_flow_link_t* link = graph->links + i;
            const uint32_t s = leaf_rep[link->src];
            const uint32_t d = leaf_rep[link->dst];
            if (s == MD_FLOW_INVALID_INDEX || d == MD_FLOW_INVALID_INDEX) continue;
            scratch[count].src    = s;
            scratch[count].dst    = d;
            scratch[count].weight = link->weight;
            count++;
        }

        qsort(scratch, count, sizeof(md_flow_link_t), flow_link_cmp);

        for (size_t i = 0; i < count;) {
            size_t j = i + 1;
            double sum = (double)scratch[i].weight;
            while (j < count && scratch[j].src == scratch[i].src && scratch[j].dst == scratch[i].dst) {
                sum += (double)scratch[j].weight;
                j++;
            }
            const md_flow_link_t merged = { .src = scratch[i].src, .dst = scratch[i].dst, .weight = (float)sum };
            md_array_push(cut->links, merged, cut->alloc);
            i = j;
        }
    }

    md_temp_end(temp);
}

md_flow_result_t md_flow_cut_validate(const md_flow_cut_t* cut, const md_flow_graph_t* graph, float eps) {
    ASSERT(cut);
    ASSERT(graph);

    if (md_array_size(cut->visible) == 0) {
        return MD_FLOW_OK;
    }

    md_temp_scope_t temp = md_temp_begin();
    const uint32_t cut_space = (uint32_t)(md_array_size(graph->nodes) + md_array_size(cut->other));
    double* column_sum = md_temp_alloc_zero_array(temp, double, graph->num_columns);
    double* band_sum   = md_temp_alloc_zero_array(temp, double, graph->num_columns);
    double* out_sum    = md_temp_alloc_zero_array(temp, double, cut_space);
    double* in_sum     = md_temp_alloc_zero_array(temp, double, cut_space);

    md_flow_result_t result = MD_FLOW_OK;

    const size_t num_visible = md_array_size(cut->visible);
    for (size_t i = 0; i < num_visible; ++i) {
        const md_flow_node_t* node = md_flow_cut_node(cut, graph, cut->visible[i]);
        column_sum[node->column] += (double)node->weight;
    }
    const size_t num_links = md_array_size(cut->links);
    for (size_t i = 0; i < num_links; ++i) {
        const md_flow_node_t* src = md_flow_cut_node(cut, graph, cut->links[i].src);
        band_sum[src->column] += (double)cut->links[i].weight;
        out_sum[cut->links[i].src] += (double)cut->links[i].weight;
        in_sum[cut->links[i].dst]  += (double)cut->links[i].weight;
    }

    for (size_t i = 0; i < num_visible && result == MD_FLOW_OK; ++i) {
        const uint32_t n = cut->visible[i];
        const md_flow_node_t* node = md_flow_cut_node(cut, graph, n);
        if (node->column + 1 < graph->num_columns && fabs(out_sum[n] - (double)node->weight) > (double)eps) {
            result = MD_FLOW_ERR_NODE_FLOW;
        }
        if (node->column > 0 && fabs(in_sum[n] - (double)node->weight) > (double)eps) {
            result = MD_FLOW_ERR_NODE_FLOW;
        }
    }

    for (uint32_t c = 0; c < graph->num_columns && result == MD_FLOW_OK; ++c) {
        if (fabs(column_sum[c] - 1.0) > (double)eps) {
            result = MD_FLOW_ERR_COLUMN_SUM;
        }
    }
    for (uint32_t c = 0; c + 1 < graph->num_columns && result == MD_FLOW_OK; ++c) {
        if (fabs(band_sum[c] - 1.0) > (double)eps) {
            result = MD_FLOW_ERR_BAND_SUM;
        }
    }

    md_temp_end(temp);
    return result;
}

// ---------------------------------------------------------------------------
// Layout
// ---------------------------------------------------------------------------
//
// Everything below works in graph space: the flow axis (x) runs 0..1 across the columns and the
// weight axis (y) runs 0..1 down each column. Pan, zoom and pixels are the caller's business,
// which is what lets a zoom happen without a relayout.

md_flow_layout_params_t md_flow_layout_params_default(void) {
    md_flow_layout_params_t params = {
        .node_thickness    = 0.12f,
        .node_gap          = 0.02f,
        .order             = MD_FLOW_ORDER_WEIGHT,
        .barycentre_sweeps = 2,
    };
    for (uint32_t c = 0; c < MD_FLOW_MAX_COLUMNS; ++c) params.weight_exponent[c] = 1.0f;
    return params;
}

void md_flow_layout_init(md_flow_layout_t* layout, struct md_allocator_i* alloc) {
    ASSERT(layout);
    ASSERT(alloc);
    MEMSET(layout, 0, sizeof(md_flow_layout_t));
    layout->alloc = alloc;
}

void md_flow_layout_free(md_flow_layout_t* layout) {
    ASSERT(layout);
    if (layout->alloc) {
        md_array_free(layout->nodes,   layout->alloc);
        md_array_free(layout->ribbons, layout->alloc);
        md_array_free(layout->bands,   layout->alloc);
        md_array_free(layout->order,   layout->alloc);
    }
    MEMSET(layout, 0, sizeof(md_flow_layout_t));
}

// Stable insertion sort over an index array, keyed by a parallel float array. Columns and bands
// hold tens of entries at most, so this beats reaching for qsort with no context parameter.
static void flow_sort_by_key(uint32_t* indices, const float* key, size_t count) {
    for (size_t i = 1; i < count; ++i) {
        const uint32_t v = indices[i];
        const float    k = key[v];
        size_t j = i;
        while (j > 0 && key[indices[j - 1]] > k) {
            indices[j] = indices[j - 1];
            j--;
        }
        indices[j] = v;
    }
}


// Sort keys for one column: (which family, then where inside it).
//
// The families are ranked by an INTEGER rather than by their raw score. That matters: two families
// with identical scores - two groups of equal weight, say - would interleave under a float key,
// because their members' scores then decide the global order and siblings end up separated. An
// integer rank cannot tie, so a family always occupies one contiguous run.
static void flow_column_keys(float* out_key, const uint32_t* nodes, uint32_t count, const uint32_t* root_of,
                             const float* family_score, const float* member_score, float* scratch_score,
                             uint32_t* scratch_root)
{
    uint32_t num_roots = 0;
    for (uint32_t i = 0; i < count; ++i) {
        const uint32_t r = root_of[nodes[i]];
        bool seen = false;
        for (uint32_t j = 0; j < num_roots; ++j) {
            if (scratch_root[j] == r) { seen = true; break; }
        }
        if (!seen) {
            scratch_root[num_roots] = r;
            scratch_score[num_roots] = family_score[nodes[i]];
            num_roots++;
        }
    }

    for (uint32_t i = 0; i < count; ++i) {
        const uint32_t r = root_of[nodes[i]];
        uint32_t rank = 0;
        float score = 0.0f;
        for (uint32_t j = 0; j < num_roots; ++j) {
            if (scratch_root[j] == r) { score = scratch_score[j]; break; }
        }
        // Rank counts the families that score HIGHER, so the highest-scoring family gets rank 0
        // and lands at the top. Callers that want "small value at the top" negate their score.
        for (uint32_t j = 0; j < num_roots; ++j) {
            if (scratch_score[j] > score || (scratch_score[j] == score && scratch_root[j] < r)) {
                rank++;
            }
        }
        // Member scores are in [-1, 1] for weights and [0, count) for barycentres, so the family
        // stride has to clear the larger of the two.
        out_key[nodes[i]] = (float)rank * (float)(count + 2) + member_score[nodes[i]];
    }
}

void md_flow_layout_compute(md_flow_layout_t* layout, const md_flow_graph_t* graph, const md_flow_cut_t* cut,
                            const md_flow_layout_params_t* params)
{
    ASSERT(layout);
    ASSERT(graph);
    ASSERT(cut);

    md_array_shrink(layout->nodes,   0);
    md_array_shrink(layout->ribbons, 0);
    md_array_shrink(layout->bands,   0);
    md_array_shrink(layout->order,   0);

    const size_t num_visible = md_array_size(cut->visible);
    const uint32_t num_columns = graph->num_columns;
    if (num_visible == 0 || num_columns == 0) {
        return;
    }

    const md_flow_layout_params_t p = params ? *params : md_flow_layout_params_default();

    md_temp_scope_t temp = md_temp_begin();

    // Cut index space is dense: graph nodes first, then the synthesized "Others".
    const uint32_t cut_space = (uint32_t)(md_array_size(graph->nodes) + md_array_size(cut->other));

    float*    row      = md_temp_alloc_array(temp, float,    cut_space);
    float*    bary     = md_temp_alloc_array(temp, float,    cut_space);
    uint32_t* root_of  = md_temp_alloc_array(temp, uint32_t, cut_space);
    float*    y_min    = md_temp_alloc_zero_array(temp, float, cut_space);
    float*    y_max    = md_temp_alloc_zero_array(temp, float, cut_space);
    float*    key      = md_temp_alloc_zero_array(temp, float, cut_space);
    float*    fam_key  = md_temp_alloc_zero_array(temp, float, cut_space);
    float*    mem_key  = md_temp_alloc_zero_array(temp, float, cut_space);
    float*    out_cur  = md_temp_alloc_zero_array(temp, float, cut_space);
    float*    in_cur   = md_temp_alloc_zero_array(temp, float, cut_space);
    uint32_t* col_of   = md_temp_alloc_array(temp, uint32_t, cut_space);

    for (uint32_t i = 0; i < cut_space; ++i) {
        row[i] = 0.0f;
        col_of[i] = MD_FLOW_INVALID_INDEX;
        root_of[i] = i;
    }

    // A visible node's outermost ancestor. Crossing reduction orders by this first, so that
    // expanding a group shows its members side by side in the group's own slot instead of
    // scattering them down the column - an expansion that shuffles the whole column reads as a
    // different diagram rather than as a closer look at the same one.
    {
        const uint32_t num_graph_nodes = (uint32_t)md_array_size(graph->nodes);
        for (uint32_t i = 0; i < num_graph_nodes; ++i) {
            uint32_t r = i;
            while (graph->nodes[r].parent != MD_FLOW_INVALID_INDEX) {
                r = graph->nodes[r].parent;
            }
            root_of[i] = r;
        }
    }

    // Bucket the visible nodes by column, preserving insertion order as the initial row order.
    uint32_t* col_count = md_temp_alloc_zero_array(temp, uint32_t, num_columns);
    for (size_t i = 0; i < num_visible; ++i) {
        const uint32_t cut_idx = cut->visible[i];
        const md_flow_node_t* node = md_flow_cut_node(cut, graph, cut_idx);
        if (!node) continue;
        col_of[cut_idx] = node->column;
        row[cut_idx] = (float)col_count[node->column]++;
    }

    uint32_t max_col_count = 0;
    for (uint32_t c = 0; c < num_columns; ++c) {
        max_col_count = MAX(max_col_count, col_count[c]);
    }
    if (max_col_count == 0) {
        md_temp_end(temp);
        return;
    }

    // Column-major index lists, so a sweep can reorder one column at a time.
    uint32_t** col_nodes = md_temp_alloc_array(temp, uint32_t*, num_columns);
    for (uint32_t c = 0; c < num_columns; ++c) {
        col_nodes[c] = col_count[c] ? md_temp_alloc_array(temp, uint32_t, col_count[c]) : NULL;
    }
    {
        uint32_t* fill = md_temp_alloc_zero_array(temp, uint32_t, num_columns);
        for (size_t i = 0; i < num_visible; ++i) {
            const uint32_t cut_idx = cut->visible[i];
            const uint32_t c = col_of[cut_idx];
            if (c == MD_FLOW_INVALID_INDEX) continue;
            col_nodes[c][fill[c]++] = cut_idx;
        }
    }

    // Crossing reduction: order each column by the weighted mean row of its neighbours in the
    // column the sweep came from. Two sweeps is plenty at these sizes and it is the single
    // largest readability win over a fixed order.
    const size_t num_cut_links = md_array_size(cut->links);

    // Rank by contribution: heaviest first, hierarchy respected. The ancestor's own weight ranks
    // the families and a node's weight ranks it inside its family, so a column reads as a ranking
    // at every level of expansion rather than only at the top one.
    if (p.order == MD_FLOW_ORDER_WEIGHT) {
        float*    scratch_score = md_temp_alloc_array(temp, float,    cut_space);
        uint32_t* scratch_root  = md_temp_alloc_array(temp, uint32_t, cut_space);

        for (uint32_t c = 0; c < num_columns; ++c) {
            if (col_count[c] < 2) continue;
            for (uint32_t i = 0; i < col_count[c]; ++i) {
                const uint32_t n = col_nodes[c][i];
                const uint32_t r = root_of[n];
                // A family's score is the ancestor's own weight when it has one, and the node's
                // own weight when the node is its own root.
                fam_key[n] = (r != n && r < (uint32_t)md_array_size(graph->nodes))
                           ? graph->nodes[r].weight
                           : md_flow_cut_node(cut, graph, n)->weight;
                // Negated, because the sort is ascending and the heaviest belongs at the top.
                mem_key[n] = -md_flow_cut_node(cut, graph, n)->weight;
            }
            flow_column_keys(key, col_nodes[c], col_count[c], root_of, fam_key, mem_key,
                             scratch_score, scratch_root);
            flow_sort_by_key(col_nodes[c], key, col_count[c]);
            for (uint32_t i = 0; i < col_count[c]; ++i) {
                row[col_nodes[c][i]] = (float)i;
            }
        }
    }

    float*    sweep_score = md_temp_alloc_array(temp, float,    cut_space);
    uint32_t* sweep_root  = md_temp_alloc_array(temp, uint32_t, cut_space);
    const uint32_t sweeps = (p.order == MD_FLOW_ORDER_BARYCENTRE) ? p.barycentre_sweeps : 0u;
    for (uint32_t sweep = 0; sweep < sweeps; ++sweep) {
        for (int dir = 0; dir < 2; ++dir) {
            const bool forward = (dir == 0);
            for (uint32_t step = 0; step + 1 < num_columns; ++step) {
                const uint32_t c = forward ? (step + 1) : (num_columns - 2 - step);
                if (col_count[c] < 2) continue;

                for (uint32_t i = 0; i < col_count[c]; ++i) {
                    bary[col_nodes[c][i]] = row[col_nodes[c][i]];
                }

                for (uint32_t i = 0; i < col_count[c]; ++i) {
                    const uint32_t n = col_nodes[c][i];
                    double num = 0.0;
                    double den = 0.0;
                    for (size_t l = 0; l < num_cut_links; ++l) {
                        const md_flow_link_t* link = cut->links + l;
                        const uint32_t self  = forward ? link->dst : link->src;
                        const uint32_t other = forward ? link->src : link->dst;
                        if (self != n) continue;
                        num += (double)row[other] * (double)link->weight;
                        den += (double)link->weight;
                    }
                    if (den > 0.0) {
                        bary[n] = (float)(num / den);
                    }
                }

                // Order by (which family, own barycentre), the family ranked by the weighted mean
                // barycentre of its members. Negated below because a family's rank counts how many
                // score LOWER, and here a smaller barycentre belongs nearer the top.
                for (uint32_t i = 0; i < col_count[c]; ++i) {
                    const uint32_t n = col_nodes[c][i];
                    const uint32_t r = root_of[n];
                    double num = 0.0;
                    double den = 0.0;
                    for (uint32_t j = 0; j < col_count[c]; ++j) {
                        const uint32_t m = col_nodes[c][j];
                        if (root_of[m] != r) continue;
                        const float w = md_flow_cut_node(cut, graph, m)->weight;
                        num += (double)bary[m] * (double)w;
                        den += (double)w;
                    }
                    fam_key[n] = -((den > 0.0) ? (float)(num / den) : bary[n]);  // small bary -> top
                    mem_key[n] = bary[n];
                }
                flow_column_keys(key, col_nodes[c], col_count[c], root_of, fam_key, mem_key,
                                 sweep_score, sweep_root);

                flow_sort_by_key(col_nodes[c], key, col_count[c]);
                for (uint32_t i = 0; i < col_count[c]; ++i) {
                    row[col_nodes[c][i]] = (float)i;
                }
            }
        }
    }

    // One weight scale shared by every column, so a ribbon has the same thickness at both ends.
    // The gap budget is set by the busiest column; shorter columns are centred in the area.
    float gap = p.node_gap;
    float usable = 1.0f - gap * (float)(max_col_count > 0 ? max_col_count - 1 : 0);
    if (usable < 0.1f) {
        // Too many nodes for the requested gap - shrink the gap rather than invert the layout.
        usable = 0.1f;
        gap = (max_col_count > 1) ? (1.0f - usable) / (float)(max_col_count - 1) : 0.0f;
    }

    const float thickness = CLAMP(p.node_thickness, 0.0f, 1.0f);
    const float span = (num_columns > 1) ? (1.0f - thickness) / (float)(num_columns - 1) : 0.0f;

    // A node's drawn share of its column, after the column's own exponent. Kept separately from
    // the node's weight, which stays the truth the tooltips and the invariants speak in.
    float* draw_share = md_temp_alloc_zero_array(temp, float, cut_space);

    for (uint32_t c = 0; c < num_columns; ++c) {
        if (col_count[c] == 0) continue;

        const float e = (c < MD_FLOW_MAX_COLUMNS) ? p.weight_exponent[c] : 1.0f;

        double share_sum = 0.0;
        for (uint32_t i = 0; i < col_count[c]; ++i) {
            const uint32_t n = col_nodes[c][i];
            const float w = md_flow_cut_node(cut, graph, n)->weight;
            const double sh = (e == 1.0f) ? (double)w : pow((double)MAX(w, 0.0f), (double)e);
            draw_share[n] = (float)sh;
            share_sum += sh;
        }
        if (share_sum <= 0.0) continue;
        for (uint32_t i = 0; i < col_count[c]; ++i) {
            draw_share[col_nodes[c][i]] /= (float)share_sum;
        }

        // Compressed or not, a column fills the same band - only the split between its rows moves.
        const float height = usable + gap * (float)(col_count[c] - 1);
        float y = (1.0f - height) * 0.5f;

        const float x0 = (num_columns > 1) ? (float)c * span : (1.0f - thickness) * 0.5f;
        const float x1 = x0 + thickness;

        for (uint32_t i = 0; i < col_count[c]; ++i) {
            const uint32_t n = col_nodes[c][i];
            const float h = draw_share[n] * usable;

            y_min[n] = y;
            y_max[n] = y + h;
            out_cur[n] = y;
            in_cur[n] = y;

            const md_flow_layout_node_t ln = {
                .cut_idx = n,
                .min = { x0, y },
                .max = { x1, y + h },
            };
            md_array_push(layout->nodes, ln, layout->alloc);
            md_array_push(layout->order, n, layout->alloc);

            y += h + gap;
        }
    }

    // Bands: the extent of every expanded ancestor, taken from where its descendants ended up.
    // Derived rather than laid out, so a band can never disagree with the nodes it brackets.
    {
        const uint32_t num_graph_nodes = (uint32_t)md_array_size(graph->nodes);
        float* band_min = md_temp_alloc_array(temp, float, num_graph_nodes ? num_graph_nodes : 1);
        float* band_max = md_temp_alloc_array(temp, float, num_graph_nodes ? num_graph_nodes : 1);
        bool*  band_any = md_temp_alloc_zero_array(temp, bool, num_graph_nodes ? num_graph_nodes : 1);

        for (uint32_t i = 0; i < num_graph_nodes; ++i) {
            band_min[i] =  1.0e30f;
            band_max[i] = -1.0e30f;
        }

        for (size_t i = 0; i < num_visible; ++i) {
            const uint32_t n = cut->visible[i];
            if (n >= num_graph_nodes) continue;     // an "Others" node has no ancestry
            if (col_of[n] == MD_FLOW_INVALID_INDEX) continue;

            for (uint32_t p = graph->nodes[n].parent; p != MD_FLOW_INVALID_INDEX; p = graph->nodes[p].parent) {
                band_any[p] = true;
                band_min[p] = MIN(band_min[p], y_min[n]);
                band_max[p] = MAX(band_max[p], y_max[n]);
            }
        }

        for (uint32_t i = 0; i < num_graph_nodes; ++i) {
            if (!band_any[i]) continue;

            const uint32_t c = graph->nodes[i].column;
            const float x0 = (num_columns > 1) ? (float)c * span : (1.0f - thickness) * 0.5f;

            uint32_t depth = 0;
            for (uint32_t p = graph->nodes[i].parent; p != MD_FLOW_INVALID_INDEX; p = graph->nodes[p].parent) {
                depth++;
            }

            const md_flow_layout_band_t band = {
                .node   = i,
                .column = c,
                .depth  = depth,
                .min    = { x0, band_min[i] },
                .max    = { x0 + thickness, band_max[i] },
            };
            md_array_push(layout->bands, band, layout->alloc);
        }
    }

    // Ribbon stubs. Outgoing are allocated in destination-row order and incoming in source-row
    // order; that alone removes most self-crossings and costs nothing.
    if (num_cut_links > 0) {
        uint32_t* by_out  = md_temp_alloc_array(temp, uint32_t, num_cut_links);
        uint32_t* by_in   = md_temp_alloc_array(temp, uint32_t, num_cut_links);
        float*    link_key = md_temp_alloc_array(temp, float,   num_cut_links);

        // Sort keys pack the two row indices into one float. Rows are small integers, so a
        // scale of 4096 keeps the pairs ordered without needing a comparator with context.
        for (size_t l = 0; l < num_cut_links; ++l) {
            by_out[l] = (uint32_t)l;
            by_in[l]  = (uint32_t)l;
        }

        for (size_t l = 0; l < num_cut_links; ++l) {
            link_key[l] = row[cut->links[l].src] * 4096.0f + row[cut->links[l].dst];
        }
        flow_sort_by_key(by_out, link_key, num_cut_links);

        // A stub takes the link's share OF ITS OWN NODE, scaled by that node's drawn height. With
        // every exponent at 1 this is identical to weight * usable; with a compressed column it is
        // what keeps the stubs exactly filling the node instead of overflowing it.
        float* out_span_min = md_temp_alloc_array(temp, float, num_cut_links);
        float* out_span_max = md_temp_alloc_array(temp, float, num_cut_links);
        for (size_t i = 0; i < num_cut_links; ++i) {
            const uint32_t l = by_out[i];
            const uint32_t s = cut->links[l].src;
            const md_flow_node_t* src = md_flow_cut_node(cut, graph, s);
            const float denom = (src && src->weight > 0.0f) ? src->weight : 1.0f;
            const float h = (cut->links[l].weight / denom) * (y_max[s] - y_min[s]);
            out_span_min[l] = out_cur[s];
            out_span_max[l] = out_cur[s] + h;
            out_cur[s] += h;
        }

        for (size_t l = 0; l < num_cut_links; ++l) {
            link_key[l] = row[cut->links[l].dst] * 4096.0f + row[cut->links[l].src];
        }
        flow_sort_by_key(by_in, link_key, num_cut_links);

        for (size_t i = 0; i < num_cut_links; ++i) {
            const uint32_t l = by_in[i];
            const uint32_t s = cut->links[l].src;
            const uint32_t d = cut->links[l].dst;
            const md_flow_node_t* dst = md_flow_cut_node(cut, graph, d);
            const float denom = (dst && dst->weight > 0.0f) ? dst->weight : 1.0f;
            const float h = (cut->links[l].weight / denom) * (y_max[d] - y_min[d]);

            const uint32_t src_col = col_of[s];
            const uint32_t dst_col = col_of[d];
            if (src_col == MD_FLOW_INVALID_INDEX || dst_col == MD_FLOW_INVALID_INDEX) continue;

            const float src_x = ((num_columns > 1) ? (float)src_col * span : (1.0f - thickness) * 0.5f) + thickness;
            const float dst_x = (num_columns > 1) ? (float)dst_col * span : (1.0f - thickness) * 0.5f;

            const md_flow_layout_ribbon_t ribbon = {
                .link_idx = l,
                .p0 = { src_x, out_span_min[l] },
                .p1 = { src_x, out_span_max[l] },
                .q0 = { dst_x, in_cur[d] },
                .q1 = { dst_x, in_cur[d] + h },
            };
            md_array_push(layout->ribbons, ribbon, layout->alloc);
            in_cur[d] += h;
        }
    }

    md_temp_end(temp);
}

static const md_flow_layout_node_t* flow_find_node(const md_flow_layout_t* layout, uint32_t cut_idx) {
    const size_t count = md_array_size(layout->nodes);
    for (size_t i = 0; i < count; ++i) {
        if (layout->nodes[i].cut_idx == cut_idx) return layout->nodes + i;
    }
    return NULL;
}

static const md_flow_layout_ribbon_t* flow_find_ribbon(const md_flow_layout_t* layout, uint32_t link_idx) {
    const size_t count = md_array_size(layout->ribbons);
    for (size_t i = 0; i < count; ++i) {
        if (layout->ribbons[i].link_idx == link_idx) return layout->ribbons + i;
    }
    return NULL;
}

static vec2_t flow_lerp2(vec2_t a, vec2_t b, float t) {
    vec2_t r = { a.x + (b.x - a.x) * t, a.y + (b.y - a.y) * t };
    return r;
}

void md_flow_layout_lerp(md_flow_layout_t* out, const md_flow_layout_t* a, const md_flow_layout_t* b, float t) {
    ASSERT(out);
    ASSERT(a);
    ASSERT(b);

    md_array_shrink(out->nodes,   0);
    md_array_shrink(out->ribbons, 0);
    md_array_shrink(out->bands,   0);
    md_array_shrink(out->order,   0);

    t = CLAMP(t, 0.0f, 1.0f);

    // Anything present on only one side collapses to zero extent on the other, so an expansion
    // reads as a node splitting rather than as a pop.
    const size_t num_b_nodes = md_array_size(b->nodes);
    for (size_t i = 0; i < num_b_nodes; ++i) {
        const md_flow_layout_node_t* nb = b->nodes + i;
        const md_flow_layout_node_t* na = flow_find_node(a, nb->cut_idx);
        md_flow_layout_node_t n = *nb;
        if (na) {
            n.min = flow_lerp2(na->min, nb->min, t);
            n.max = flow_lerp2(na->max, nb->max, t);
        } else {
            const vec2_t mid = { (nb->min.x + nb->max.x) * 0.5f, (nb->min.y + nb->max.y) * 0.5f };
            n.min = flow_lerp2(mid, nb->min, t);
            n.max = flow_lerp2(mid, nb->max, t);
        }
        md_array_push(out->nodes, n, out->alloc);
        md_array_push(out->order, n.cut_idx, out->alloc);
    }

    const size_t num_a_nodes = md_array_size(a->nodes);
    for (size_t i = 0; i < num_a_nodes; ++i) {
        const md_flow_layout_node_t* na = a->nodes + i;
        if (flow_find_node(b, na->cut_idx)) continue;
        const vec2_t mid = { (na->min.x + na->max.x) * 0.5f, (na->min.y + na->max.y) * 0.5f };
        md_flow_layout_node_t n = *na;
        n.min = flow_lerp2(na->min, mid, t);
        n.max = flow_lerp2(na->max, mid, t);
        md_array_push(out->nodes, n, out->alloc);
        md_array_push(out->order, n.cut_idx, out->alloc);
    }

    const size_t num_b_ribbons = md_array_size(b->ribbons);
    for (size_t i = 0; i < num_b_ribbons; ++i) {
        const md_flow_layout_ribbon_t* rb = b->ribbons + i;
        const md_flow_layout_ribbon_t* ra = flow_find_ribbon(a, rb->link_idx);
        md_flow_layout_ribbon_t r = *rb;
        if (ra) {
            r.p0 = flow_lerp2(ra->p0, rb->p0, t);
            r.p1 = flow_lerp2(ra->p1, rb->p1, t);
            r.q0 = flow_lerp2(ra->q0, rb->q0, t);
            r.q1 = flow_lerp2(ra->q1, rb->q1, t);
        } else {
            r.p1 = flow_lerp2(rb->p0, rb->p1, t);
            r.q1 = flow_lerp2(rb->q0, rb->q1, t);
        }
        md_array_push(out->ribbons, r, out->alloc);
    }
}
