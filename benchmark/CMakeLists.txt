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
    bench_coords.c
    bench_lammps.c
)

if (MD_ENABLE_HDF5)
    set (SRC_FILES ${SRC_FILES} bench_gto.c)
endif()

add_executable(md_benchmark ${SRC_FILES})

target_compile_definitions(md_benchmark PRIVATE MD_BENCHMARK_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/../test_data\" ${MD_DEFINES})
target_compile_options(md_benchmark PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
target_link_options(md_benchmark PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
target_link_libraries(md_benchmark PRIVATE mdlib)

target_include_directories(md_benchmark
  PRIVATE
  gen/
  ../ext/xtc/
)

# Timings of the certified GTO critical point search, CPU and GPU, per phase (see bench_topo_gto.c).
# Its own executable: it prints tables rather than ubench statistics, and takes options.
add_executable(md_topo_gto_bench bench_topo_gto.c)
target_compile_definitions(md_topo_gto_bench PRIVATE MD_BENCHMARK_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/../test_data\" ${MD_DEFINES})
target_compile_options(md_topo_gto_bench PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
target_link_options(md_topo_gto_bench PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
target_link_libraries(md_topo_gto_bench PRIVATE mdlib)

# The mdlib commit the benchmark is built from ("branch@hash", plus any locally changed tracked files),
# refreshed on every build and printed in the header, so each result says which source produced it.
add_custom_target(md_topo_gto_bench_rev
    COMMAND ${CMAKE_COMMAND} -DSRC=${CMAKE_CURRENT_SOURCE_DIR}/.. -DOUT=${CMAKE_CURRENT_BINARY_DIR}/bench_rev.h
            -P ${CMAKE_CURRENT_SOURCE_DIR}/bench_version.cmake
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/bench_rev.h
    VERBATIM)
add_dependencies(md_topo_gto_bench md_topo_gto_bench_rev)
target_include_directories(md_topo_gto_bench PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
