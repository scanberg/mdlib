// topo_setup.slang
// Single-invocation (1,1,1) pass run after critical_point detection.
// Reads the 4 raw type counts and produces:
//   counters[4]    – write-cursor base offsets for the compaction pass
//   type_counts[5] – per-type counts + total vertex count for the extraction pass

struct RootArgs {
    uint* counts;       // [0]=maxima [1]=split_saddles [2]=minima [3]=join_saddles
    uint* counters;     // [0..3] output: compaction write-cursor initial values
    uint* type_counts;  // [0..4] output: per-type counts + total
};

[[vk::push_constant]]
ConstantBuffer<RootArgs> args;

[shader("compute")]
[numthreads(1, 1, 1)]
void main() {
    uint n0 = args.counts[0]; // maxima
    uint n1 = args.counts[1]; // split saddles
    uint n2 = args.counts[2]; // minima
    uint n3 = args.counts[3]; // join saddles

    // Prefix-sum → write-cursor base offsets (type-sorted layout)
    args.counters[0] = 0;
    args.counters[1] = n0;
    args.counters[2] = n0 + n1;
    args.counters[3] = n0 + n1 + n2;

    // Per-type counts and total for the extraction pass
    args.type_counts[0] = n0;
    args.type_counts[1] = n1;
    args.type_counts[2] = n2;
    args.type_counts[3] = n3;
    args.type_counts[4] = n0 + n1 + n2 + n3;
}
