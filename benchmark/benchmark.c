#include "ubench.h"
#include "bench_rev.h"      // MD_BENCH_SOURCE_REV, generated at build time by bench_version.cmake

#include <stdio.h>

UBENCH_STATE();

int main(int argc, const char *const argv[]) {
	// The source measured, so that a result is never one of a stale build
	printf("mdlib %s\n", MD_BENCH_SOURCE_REV);
	return ubench_main(argc, argv);
}
