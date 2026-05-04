// topo_setup.slang
// Single-invocation (1,1,1) pass run after critical_point detection.
// Reads the 4 raw type counts and produces:
//   counters[4]    – write-cursor base offsets for the compaction pass
//   type_counts[5] – per-type counts + total vertex count for the extraction pass

[[vk::binding(0)]] RWStructuredBuffer<uint> counts;       // [0]=maxima [1]=split_saddles [2]=minima [3]=join_saddles
[[vk::binding(1)]] RWStructuredBuffer<uint> counters;     // [0..3] output: compaction write-cursor initial values
[[vk::binding(2)]] RWStructuredBuffer<uint> type_counts;  // [0..4] output: per-type counts + total

[shader("compute")]
[numthreads(1, 1, 1)]
void main() {
    uint n0 = counts[0]; // maxima
    uint n1 = counts[1]; // split saddles
    uint n2 = counts[2]; // minima
    uint n3 = counts[3]; // join saddles

    // Prefix-sum → write-cursor base offsets (type-sorted layout)
    counters[0] = 0;
    counters[1] = n0;
    counters[2] = n0 + n1;
    counters[3] = n0 + n1 + n2;

    // Per-type counts and total for the extraction pass
    type_counts[0] = n0;
    type_counts[1] = n1;
    type_counts[2] = n2;
    type_counts[3] = n3;
    type_counts[4] = n0 + n1 + n2 + n3;
}
