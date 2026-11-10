// md_parse_bench: what md_parse_f64 costs per number, by the shape the number is written in.
//
//   md_parse_bench [--reps n] [--count n] [--seed n]
//
//   --reps <n>   timed passes over every set, the best and the median are reported (default 9)
//   --count <n>  numbers in each set (default 1000000)
//   --seed <n>   of the generator the numbers are drawn from (default 1)
//
// Each set is a million random numbers printed the way a format mdlib reads prints them, parsed one after another from
// one buffer, so the lengths and the digits vary as they do in a file and the branch predictor cannot learn the input.
// Each line is a set: best and median ns per number for md_parse_f64 and for strtod on the same text, the numbers that
// differ from strtod in any bit (must be 0), and a hash of the bits md_parse_f64 gave. The hash depends only on the
// input, so two builds with the same seed and count must print the same hashes: a change of md_parse that keeps them is
// one that keeps every value. The header names the benchmark's revision (bumped with every change of it) and the mdlib
// commit it was built from.

#include <core/md_parse.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include "bench_rev.h"                    // MD_BENCH_SOURCE_REV, generated at build time by bench_version.cmake

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// The benchmark's own revision: bumped with every change of it, so that a result names the build it came from
#define MD_PARSE_BENCH_REV "R1"

#define MAX_REPS 64

typedef enum {
    SET_FIXED3,     // "%.3f" in [-1000, 1000): pdb and gro coordinates
    SET_FIXED6,     // "%.6f" in [-100, 100): xyz, mol2, lammps data
    SET_PRMTOP,     // "%.8E" in [-1e3, 1e3] over many decades: amber prmtop and rst7 fields
    SET_CUBE,       // "%.5E" over many decades, about zero: cube file values
    SET_ROUNDTRIP,  // "%.17g" of an arbitrary double: what a JSON writer emits
    SET_INT,        // "%d": an integer field read as a float (csv, a count)
    SET_MIXED,      // all of the above, interleaved at random
    SET_COUNT
} set_kind_t;

static const char* set_name[SET_COUNT] = { "fixed %.3f", "fixed %.6f", "prmtop %.8E", "cube %.5E", "roundtrip %.17g", "int %d", "mixed" };

static uint64_t rng_state;
static uint64_t rng_next(void) {
    // splitmix64
    uint64_t z = (rng_state += UINT64_C(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31);
}
static double rng_unit(void) { return (double)(rng_next() >> 11) * (1.0 / 9007199254740992.0); }    // [0, 1)
static double rng_range(double lo, double hi) { return lo + (hi - lo) * rng_unit(); }

static int print_number(char* buf, size_t cap, set_kind_t kind) {
    switch (kind) {
    case SET_FIXED3:    return snprintf(buf, cap, "%.3f", rng_range(-1000.0, 1000.0));
    case SET_FIXED6:    return snprintf(buf, cap, "%.6f", rng_range(-100.0, 100.0));
    case SET_PRMTOP:    return snprintf(buf, cap, "%.8E", rng_range(-1.0, 1.0) * pow(10.0, floor(rng_range(-6.0, 4.0))));
    case SET_CUBE:      return snprintf(buf, cap, "%.5E", rng_range(-1.0, 1.0) * pow(10.0, floor(rng_range(-12.0, 1.0))));
    case SET_ROUNDTRIP: return snprintf(buf, cap, "%.17g", rng_range(-1.0, 1.0) * pow(10.0, floor(rng_range(-20.0, 20.0))));
    case SET_INT:       return snprintf(buf, cap, "%d", (int)rng_range(-100000.0, 100000.0));
    default:            return print_number(buf, cap, (set_kind_t)(rng_next() % SET_MIXED));
    }
}

typedef struct {
    char*  text;    // the numbers, each terminated by a NUL so that strtod can read the same bytes
    md_str_t* tok;     // views of them, without the NUL
    size_t count;
} number_set_t;

static void set_init(number_set_t* set, set_kind_t kind, size_t count) {
    const size_t stride = 32;
    set->text  = malloc(count * stride);
    set->tok   = malloc(count * sizeof(md_str_t));
    set->count = count;
    char* c = set->text;
    for (size_t i = 0; i < count; ++i) {
        const int n = print_number(c, stride, kind);
        set->tok[i] = (md_str_t){ c, (size_t)n };
        c += n + 1;
    }
}

static void set_free(number_set_t* set) {
    free(set->text);
    free(set->tok);
}

static uint64_t bits_of(double d) {
    uint64_t b;
    memcpy(&b, &d, sizeof(b));
    return b;
}

static int cmp_double(const void* a, const void* b) {
    const double x = *(const double*)a, y = *(const double*)b;
    return (x > y) - (x < y);
}

static double parse_pass(const number_set_t* set, double* out) {
    const md_tick_t t0 = md_tick_now();
    for (size_t i = 0; i < set->count; ++i) {
        md_parse_f64(&out[i], set->tok[i]);
    }
    const md_tick_t t1 = md_tick_now();
    return md_tick_to_nanoseconds(t1 - t0) / (double)set->count;
}

static double strtod_pass(const number_set_t* set, double* out) {
    const md_tick_t t0 = md_tick_now();
    for (size_t i = 0; i < set->count; ++i) {
        out[i] = strtod(set->tok[i].ptr, NULL);
    }
    const md_tick_t t1 = md_tick_now();
    return md_tick_to_nanoseconds(t1 - t0) / (double)set->count;
}

int main(int argc, char** argv) {
    int reps = 9;
    size_t count = 1000000;
    uint64_t seed = 1;
    for (int i = 1; i < argc; ++i) {
        if      (!strcmp(argv[i], "--reps")  && i + 1 < argc) reps  = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--count") && i + 1 < argc) count = (size_t)strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--seed")  && i + 1 < argc) seed  = strtoull(argv[++i], NULL, 10);
        else {
            fprintf(stderr, "usage: md_parse_bench [--reps n] [--count n] [--seed n]\n");
            return 1;
        }
    }
    if (reps < 1) reps = 1;
    if (reps > MAX_REPS) reps = MAX_REPS;
    if (count < 1) count = 1;

    printf("md_parse_bench %s, mdlib %s\n", MD_PARSE_BENCH_REV, MD_BENCH_SOURCE_REV);
    printf("%zu numbers per set, %d reps, seed %llu; ns per number\n\n", count, reps, (unsigned long long)seed);
    printf("%-16s %10s %10s   %10s %10s   %8s   %s\n", "set", "parse best", "median", "strtod best", "median", "!=strtod", "hash");

    double* a = malloc(count * sizeof(double));
    double* b = malloc(count * sizeof(double));
    for (int k = 0; k < SET_COUNT; ++k) {
        rng_state = seed * 1000 + (uint64_t)k;
        number_set_t set;
        set_init(&set, (set_kind_t)k, count);

        double tp[MAX_REPS], ts[MAX_REPS];
        for (int r = 0; r < reps; ++r) {
            tp[r] = parse_pass(&set, a);
            ts[r] = strtod_pass(&set, b);
        }
        qsort(tp, (size_t)reps, sizeof(double), cmp_double);
        qsort(ts, (size_t)reps, sizeof(double), cmp_double);

        size_t mismatch = 0;
        uint64_t hash = UINT64_C(14695981039346656037);
        for (size_t i = 0; i < count; ++i) {
            const uint64_t x = bits_of(a[i]);
            mismatch += x != bits_of(b[i]);
            hash = (hash ^ x) * UINT64_C(1099511628211);
        }
        printf("%-16s %10.2f %10.2f   %10.2f %10.2f   %8zu   %016llx\n", set_name[k], tp[0], tp[reps / 2], ts[0], ts[reps / 2],
               mismatch, (unsigned long long)hash);
        set_free(&set);
    }
    free(a);
    free(b);
    return 0;
}
