#pragma once

// Layered flow graphs: the model behind VIAMD's charge-transfer / transition diagram.
//
// A flow graph is a layered DAG. Nodes live in columns; links only ever connect column c to
// column c+1. Within a column, nodes form a hierarchy (group -> atom -> orbital), and a node's
// children partition its weight exactly.
//
// The point of the split between md_flow_graph_t and md_flow_cut_t: the GRAPH is built once,
// at the finest level available, from whatever produced it. The CUT is what the user is looking
// at right now - which nodes they expanded, what they thresholded away. Every coarser view is
// derived from the graph by aggregation, never by asking the source for different numbers. That
// is what makes it impossible for the collapsed view and the expanded view to disagree.
//
// Layout is computed in GRAPH SPACE (x = column, y in [0,1] by cumulative weight). Pan and zoom
// belong to whoever draws it, so changing the zoom never triggers a relayout.
//
// See docs/transition_flow_design.md in the viamd repository for the reasoning.

#include <core/md_array.h>
#include <core/md_bitfield.h>
#include <core/md_str.h>
#include <core/md_vec_math.h>

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

struct md_allocator_i;

#define MD_FLOW_INVALID_INDEX 0xFFFFFFFFu
#define MD_FLOW_MAX_COLUMNS   8

typedef enum md_flow_node_flag_t {
    MD_FLOW_NODE_FLAG_NONE  = 0,
    // Synthesized by the cut to absorb every sibling below the threshold. Never present in a
    // graph, only ever in a cut, and it is expandable like any other node.
    MD_FLOW_NODE_FLAG_OTHER = 1,
} md_flow_node_flag_t;

typedef struct md_flow_node_t {
    uint32_t column;    // 0 .. num_columns-1
    uint32_t parent;    // coarser node in the SAME column, or MD_FLOW_INVALID_INDEX
    uint32_t level;     // 0 = coarsest. Purely descriptive; the hierarchy is 'parent'.
    float    weight;    // share of this node's column total, [0,1]
    vec4_t   color;
    str_t    label;     // NOT owned - must outlive the graph
    str_t    sublabel;  // optional second line, same ownership rule; {0} for none
    uint64_t key;       // stable identity across rebuilds; interaction state keys off this
    uint32_t flags;     // md_flow_node_flag_t
} md_flow_node_t;

typedef struct md_flow_link_t {
    uint32_t src;       // node in column c
    uint32_t dst;       // node in column c+1
    float    weight;    // share of the grand total
} md_flow_link_t;

typedef struct md_flow_graph_t {
    md_array(md_flow_node_t) nodes;
    md_array(md_flow_link_t) links;
    uint32_t num_columns;
    struct md_allocator_i* alloc;
} md_flow_graph_t;

// ---------------------------------------------------------------------------
// Building
// ---------------------------------------------------------------------------

void     md_flow_graph_init (md_flow_graph_t* graph, uint32_t num_columns, struct md_allocator_i* alloc);
void     md_flow_graph_free (md_flow_graph_t* graph);
void     md_flow_graph_clear(md_flow_graph_t* graph);   // keeps the allocation

// Returns the new node's index, or MD_FLOW_INVALID_INDEX if the node is malformed (column out of
// range, or a parent that is not an earlier node in the same column).
uint32_t md_flow_graph_add_node(md_flow_graph_t* graph, const md_flow_node_t* node);

// Links are only meaningful between LEAF nodes: coarser links are derived by the cut. A link with
// non-positive weight is dropped.
bool md_flow_graph_add_link(md_flow_graph_t* graph, uint32_t src, uint32_t dst, float weight);

bool md_flow_node_is_leaf(const md_flow_graph_t* graph, uint32_t node_idx);

// ---------------------------------------------------------------------------
// Validation - the four invariants from the design doc
// ---------------------------------------------------------------------------

typedef enum md_flow_result_t {
    MD_FLOW_OK = 0,
    MD_FLOW_ERR_HIERARCHY,      // parent in another column, forward reference, or a cycle
    MD_FLOW_ERR_LAYERING,       // (1) a link that does not connect two leaves in adjacent columns
    MD_FLOW_ERR_COLUMN_SUM,     // (2) leaf weights within a column do not sum to 1
    MD_FLOW_ERR_BAND_SUM,       // (2) link weights across a band do not sum to 1
    MD_FLOW_ERR_NODE_FLOW,      // (2) a node's links do not sum to its own weight
    MD_FLOW_ERR_PARTITION,      // (3) a node's children do not partition its weight
} md_flow_result_t;

const char* md_flow_result_str(md_flow_result_t result);

// eps is an absolute tolerance on the weight sums. 1.0e-4f is a sane default.
md_flow_result_t md_flow_graph_validate(const md_flow_graph_t* graph, float eps);

// ---------------------------------------------------------------------------
// The cut: what is actually on screen
// ---------------------------------------------------------------------------
//
// A cut is one antichain per column - every leaf of the graph has exactly one ancestor-or-self
// in the cut. Expanding one node while its siblings stay collapsed is the normal case.
//
// Visible nodes are addressed in CUT INDEX SPACE: an index below md_array_size(graph->nodes) is
// a graph node, anything at or above it is a synthesized "Others". Use md_flow_cut_node() rather
// than indexing either array directly.

typedef struct md_flow_cut_t {
    md_bitfield_t expanded;                 // graph node index -> is expanded
    float threshold;                        // fold nodes below this share of their COLUMN, [0,1)

    // An "Others" never grows past this share of its column: folding proceeds from the smallest
    // node upward and stops once taking the next would cross the cap.
    //
    // DEFAULT 1 (uncapped), and it should usually stay there. A cap makes the threshold stop
    // meaning what it says - folding halts early and the LARGEST sub-threshold nodes stay on
    // screen, so a reader who asked for "above 1%" sees 0.5% rows and has no way to tell why.
    // Since an "Others" can now be opened (below), the cap is a preference rather than the fix it
    // was introduced as.
    float other_max;

    // Bit per column: that column's "Others" is open, so nothing in it folds. This is what lets
    // the threshold be honest - hide everything below it, and let the reader look anyway.
    uint32_t other_expanded;

    // Outputs of md_flow_cut_resolve()
    md_array(uint32_t)       visible;       // cut index space, grouped by column
    md_array(md_flow_link_t) links;         // aggregated, endpoints in cut index space
    md_array(md_flow_node_t) other;         // synthesized "Others" nodes
    md_array(uint32_t)       other_count;   // how many nodes each "Others" absorbed

    struct md_allocator_i* alloc;
} md_flow_cut_t;

void md_flow_cut_init(md_flow_cut_t* cut, struct md_allocator_i* alloc);
void md_flow_cut_free(md_flow_cut_t* cut);

void md_flow_cut_set_expanded(md_flow_cut_t* cut, uint32_t node_idx, bool expanded);

// Open or close a column's "Others". An open one folds nothing, so its members are drawn
// individually however small they are.
void md_flow_cut_set_other_expanded(md_flow_cut_t* cut, uint32_t column, bool expanded);
bool md_flow_cut_is_other_expanded(const md_flow_cut_t* cut, uint32_t column);
bool md_flow_cut_is_expanded (const md_flow_cut_t* cut, uint32_t node_idx);
void md_flow_cut_collapse_all(md_flow_cut_t* cut);
void md_flow_cut_expand_all  (md_flow_cut_t* cut, const md_flow_graph_t* graph);

// Recompute 'visible', 'links' and 'other'. Call after any change to expansion, threshold or the
// graph itself. Cheap - it is an aggregation over the graph, no physics involved.
void md_flow_cut_resolve(md_flow_cut_t* cut, const md_flow_graph_t* graph);

// Resolve a cut index to a node. Never returns NULL for an index that appears in 'visible'.
const md_flow_node_t* md_flow_cut_node(const md_flow_cut_t* cut, const md_flow_graph_t* graph, uint32_t cut_idx);

// Invariant (4): weight is conserved on every valid cut - for any mix of expanded and collapsed
// nodes. This is the check that says an expansion cannot make the diagram lie.
md_flow_result_t md_flow_cut_validate(const md_flow_cut_t* cut, const md_flow_graph_t* graph, float eps);

// ---------------------------------------------------------------------------
// Layout - pure, in graph space, no drawing
// ---------------------------------------------------------------------------

// How the rows of a column are ordered. Both policies respect the hierarchy first - siblings stay
// contiguous within their ancestor's slot - and differ only in how they rank within it.
typedef enum md_flow_order_t {
    // Largest contribution at the top. Makes the column readable as a ranking, at the cost of some
    // ribbon crossings.
    MD_FLOW_ORDER_WEIGHT = 0,
    // Barycentre sweeps: order by the mean position of the neighbours a node connects to. Fewer
    // crossings, but the rows are in no order the reader can name.
    MD_FLOW_ORDER_BARYCENTRE,
} md_flow_order_t;

typedef struct md_flow_layout_params_t {
    float node_thickness;       // node bar thickness along the flow axis, graph units
    float node_gap;             // gap between adjacent nodes in a column, as a share of the column
    md_flow_order_t order;
    uint32_t barycentre_sweeps; // only for MD_FLOW_ORDER_BARYCENTRE; 2 is plenty

    // Per column, the exponent applied to a node's weight before it is turned into a height:
    // h ∝ w^e, renormalized within the column. 1 = strictly proportional, the Sankey default.
    // Below 1 compresses the range, so a column whose weight is nearly all in one node still gives
    // the small ones a body you can read, label and click. 0 would make every row identical.
    //
    // This makes a column's scale LOCAL, so a ribbon spanning two columns with different exponents
    // tapers. That is not a lie: each end is proportional within its own column's scale, and the
    // stubs still exactly fill the node at both ends. It IS a reason not to compress a column whose
    // absolute proportions are the story - the atom columns - while compressing one that is really
    // a channel picker, like the NTO pairs.
    float weight_exponent[MD_FLOW_MAX_COLUMNS];
} md_flow_layout_params_t;

md_flow_layout_params_t md_flow_layout_params_default(void);

typedef struct md_flow_layout_node_t {
    uint32_t cut_idx;
    vec2_t   min;               // graph space; x is the flow axis, y is the weight axis
    vec2_t   max;
} md_flow_layout_node_t;

// A ribbon is a quad with two curved edges. p0/p1 are the near and far edge of the band where it
// leaves the source; q0/q1 where it arrives at the destination. The curve between them is the
// drawing backend's business.
typedef struct md_flow_layout_ribbon_t {
    uint32_t link_idx;          // index into cut->links
    vec2_t   p0, p1;
    vec2_t   q0, q1;
} md_flow_layout_ribbon_t;

// An ancestor that is currently EXPANDED, and therefore not drawn as a node itself. It still has
// an extent - the span its visible descendants occupy - and a diagram that shows that extent is
// one where the user can see which atoms belong to which group while looking at the atoms.
//
// Bands are what makes a group an organizing layer OUTSIDE the column rather than an alternative
// to it. The layout reports the extent and the column's own x range; which side of the column to
// put it on, and how far out, is the drawing backend's decision.
typedef struct md_flow_layout_band_t {
    uint32_t node;      // graph node index of the ancestor
    uint32_t column;
    uint32_t depth;     // 0 = outermost ancestor, 1 = its child, ...
    vec2_t   min;       // y spans the visible descendants; x matches the column's node rects
    vec2_t   max;
} md_flow_layout_band_t;

typedef struct md_flow_layout_t {
    md_array(md_flow_layout_node_t)   nodes;
    md_array(md_flow_layout_ribbon_t) ribbons;
    md_array(md_flow_layout_band_t)   bands;
    // Row order per column, cut index space. Kept so a caller can animate between two layouts.
    md_array(uint32_t) order;
    struct md_allocator_i* alloc;
} md_flow_layout_t;

void md_flow_layout_init(md_flow_layout_t* layout, struct md_allocator_i* alloc);
void md_flow_layout_free(md_flow_layout_t* layout);

// Lays out into the unit square [0,1] x [0,1] in graph space. Pan, zoom and the mapping to pixels
// belong to the caller, so zooming never needs a relayout.
void md_flow_layout_compute(md_flow_layout_t* layout, const md_flow_graph_t* graph, const md_flow_cut_t* cut,
                            const md_flow_layout_params_t* params);

// Linear interpolation between two layouts of the SAME cut index space, for expansion animation.
// Nodes and ribbons present in only one side fade by collapsing to zero extent.
void md_flow_layout_lerp(md_flow_layout_t* out, const md_flow_layout_t* a, const md_flow_layout_t* b, float t);

#ifdef __cplusplus
}
#endif
