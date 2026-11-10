project(md_benchmark LANGUAGES C)

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
    set (SRC_FILES ${SRC_FILES} bench_gto.c bench_gto_int.c)
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

# How structure files (PDB, GRO) load over threads: best and median per thread count, the speedup, and a hash of the
# system each count loads (see bench_load.c). Its own executable, which takes the files to load.
add_executable(md_load_bench bench_load.c)
target_compile_options(md_load_bench PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
target_link_options(md_load_bench PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
target_link_libraries(md_load_bench PRIVATE mdlib $<$<PLATFORM_ID:Windows>:psapi>)
add_dependencies(md_load_bench md_topo_gto_bench_rev)
target_include_directories(md_load_bench PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

# What md_parse_f64 costs per number, by the shape the number is written in, against strtod on the same text (see
# bench_parse.c). Its own executable, which prints a table and takes options.
add_executable(md_parse_bench bench_parse.c)
target_compile_options(md_parse_bench PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
target_link_options(md_parse_bench PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
target_link_libraries(md_parse_bench PRIVATE mdlib)
add_dependencies(md_parse_bench md_topo_gto_bench_rev)
target_include_directories(md_parse_bench PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
