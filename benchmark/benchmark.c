#include "ubench.h"

UBENCH_STATE();

extern void bench_com();

int main(int argc, const char *const argv[]) {
	bench_com();
	return ubench_main(argc, argv);
}
