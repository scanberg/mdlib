cmake_minimum_required(VERSION 3.15)

project(md_benchmark VERSION 0.1 LANGUAGES C)

add_executable(md_benchmark benchmark.c "bench_lammps.c")

target_compile_definitions(md_benchmark PRIVATE MD_BENCHMARK_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/../test_data\" ${MD_DEFINES})
target_compile_options(md_benchmark PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
target_link_options(md_benchmark PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
target_link_libraries(md_benchmark PRIVATE mdlib)

target_sources(md_benchmark
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/bench_allocator.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bench_spatial_hash.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bench_xtc.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bench_gro.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bench_pdb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bench_xyz.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bench_str.c
    ${CMAKE_CURRENT_SOURCE_DIR}/bench_util.c
  )