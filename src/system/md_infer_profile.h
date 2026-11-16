#pragma once

// The time each step of md_system_infer takes, for profiling it (benchmark/bench_infer.c). Not a public header.
//
// Null unless a profiler sets it, and then md_system_infer records into it the wall time of each of its steps, in the
// order they run. One call at a time: a profiler sets it, calls md_system_infer, reads it and unsets it. The molecules,
// the rings and the backbones, which run side by side, are each recorded with their own time, then together as
// "branches", the wall time of the three. A step of a step (the parts of the bonds) is named indented, and is part of
// the time of the step after it.

#include <stdint.h>

#define MDI_INFER_PROFILE_MAX_STEPS 32

typedef struct mdi_infer_profile_t {
    int         count;
    const char* name[MDI_INFER_PROFILE_MAX_STEPS];
    double      ms[MDI_INFER_PROFILE_MAX_STEPS];
} mdi_infer_profile_t;

extern mdi_infer_profile_t* mdi_infer_profile;

// The step named ended now, since *t (md_tick_now), which it sets to now. Nothing unless profiled.
void mdi_infer_profile_step(int64_t* t, const char* name);
