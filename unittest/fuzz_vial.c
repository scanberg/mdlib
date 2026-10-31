// libFuzzer target for VIAL: the front end and the checker, against T. Every input must keep the promises of vial_fuzz_check.h; a broken promise
// aborts, which libFuzzer reports with the input. Built with -DMD_VIAL_FUZZ=ON (clang only), with ASan and UBSan:
//
//   cmake -S . -B build_fuzz -G Ninja -DCMAKE_C_COMPILER=clang -DMD_VIAL_FUZZ=ON
//   ninja -C build_fuzz fuzz_vial
//   python3 unittest/vial_fuzz_seeds.py build_fuzz/vial_seeds
//   build_fuzz/bin/fuzz_vial -dict=unittest/vial_fuzz.dict -max_len=4096 -timeout=2 -max_total_time=3600 \
//       build_fuzz/vial_corpus build_fuzz/vial_seeds
//
// The built-in fuzzer in test_vial_fuzz.c runs the same check in every unittest run, for a few seconds.

#include "vial_fuzz_check.h"
#include "vial_test_systems.h"

#include <core/md_arena_allocator.h>
#include <core/md_common.h>

#include <stdio.h>
#include <stdlib.h>

int LLVMFuzzerTestOneInput(const uint8_t* data, size_t size) {
    // A copy of exactly size bytes on the heap, with no NUL after it, so ASan sees any read past the end
    char* buf = malloc(size ? size : 1);
    if (!buf) return 0;
    if (size) memcpy(buf, data, size);
    md_allocator_i* arena = md_arena_allocator_create(md_get_heap_allocator(), KILOBYTES(256));
    const char* err = vial_fuzz_check((str_t){buf, size}, vial_test_system('T'), vial_test_state('T'), arena);
    if (err) {
        fprintf(stderr, "vial_fuzz_check: %s\n", err);
        abort();
    }
    md_arena_allocator_destroy(arena);
    free(buf);
    return 0;
}
