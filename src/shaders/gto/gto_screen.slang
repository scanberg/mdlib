// gto_screen.slang
// The screening of a workgroup's CGTOs into a list in the order of their index, the same on every GPU and every run.
//
// The list had been built in the order its waves won an atomic: where the GPU runs a group's waves at once (Metal),
// that order, and with it the order of the sums over the list's pairs, changed run to run, and their last bits with
// it. Here each round of WG_SIZE CGTOs sets one bit per CGTO that passes (an atomic OR: the same in any order), and a
// passing thread's slot is the bits below its own: independent of the scheduling and of the wave size.
//
// Include after WG_SIZE (a multiple of 32) is defined. Every thread of the group calls gto_screen_round, the same
// number of times (it holds barriers).

#define GTO_SCREEN_WORDS (WG_SIZE / 32)

groupshared uint gto_screen_mask[GTO_SCREEN_WORDS];

// The round's slot of the calling thread among the threads that pass (by thread index) in x, and how many pass in y:
// the same y in every thread.
uint2 gto_screen_round(bool passed, uint tid) {
    if (tid < GTO_SCREEN_WORDS) gto_screen_mask[tid] = 0u;
    GroupMemoryBarrierWithGroupSync();
    const uint w = tid >> 5u;
    const uint bit = tid & 31u;
    if (passed) InterlockedOr(gto_screen_mask[w], 1u << bit);
    GroupMemoryBarrierWithGroupSync();
    uint rank = countbits(gto_screen_mask[w] & ((1u << bit) - 1u));
    uint count = 0u;
    for (uint k = 0u; k < GTO_SCREEN_WORDS; ++k) {
        const uint c = countbits(gto_screen_mask[k]);
        rank  += k < w ? c : 0u;
        count += c;
    }
    GroupMemoryBarrierWithGroupSync();   // the mask is read until here: the next round clears it
    return uint2(rank, count);
}
