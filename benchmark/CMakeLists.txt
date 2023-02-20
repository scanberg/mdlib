cmake_minimum_required(VERSION 3.15)

project(md_benchmark VERSION 0.1 LANGUAGES C)

add_executable(md_benchmark benchmark.c)
target_link_libraries(md_benchmark PRIVATE mdlib)
target_compile_definitions(md_benchmark PRIVATE MD_BENCHMARK_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/../test_data\")

target_sources(md_benchmark
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/test_allocator.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_spatial_hash.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_xtc.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_gro.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_pdb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_xyz.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_str.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_util.c
  )