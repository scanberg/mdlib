cmake_minimum_required(VERSION 3.15)

project(md_benchmark VERSION 0.1 LANGUAGES C)

set(SRC_FILES
    benchmark.c
    bench_allocator.c
    bench_spatial_hash.c
    bench_xtc.c
    bench_gro.c
    bench_pdb.c
    bench_xyz.c
    bench_str.c
    bench_util.c
    bench_lammps.c
)

if (${MD_ENABLE_VLX})
    set (SRC_FILES ${SRC_FILES} bench_gto.c)
endif()

add_executable(md_benchmark ${SRC_FILES})

target_compile_definitions(md_benchmark PRIVATE MD_BENCHMARK_DATA_DIR=\"${MD_BENCHMARK_DATA_DIR}\" ${MD_DEFINES})
target_compile_options(md_benchmark PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
target_link_options(md_benchmark PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
target_link_libraries(md_benchmark PRIVATE mdlib)

target_include_directories(md_benchmark
  PRIVATE
  gen/
  ../ext/xtc/
)
