project(md_benchmark LANGUAGES C)

# bench_rev.h: the mdlib commit the benchmarks are built from ("branch@hash", plus any locally changed tracked files),
# which each prints in its header, so every result says which source produced it (see bench_version.cmake).
# A target, so that it runs on every build: at configure time it would go stale on a checkout, a commit or an edit,
# none of which reconfigures, and a custom command reruns only when its listed DEPENDS change, which git state has
# none of. The header is rewritten only when the text changes, so an unchanged checkout rebuilds nothing.
add_custom_target(md_bench_rev
    COMMAND ${CMAKE_COMMAND} -DSRC=${CMAKE_CURRENT_SOURCE_DIR}/.. -DOUT=${CMAKE_CURRENT_BINARY_DIR}/bench_rev.h
            -P ${CMAKE_CURRENT_SOURCE_DIR}/bench_version.cmake
    BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/bench_rev.h
    VERBATIM)

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
add_dependencies(md_benchmark md_bench_rev)

target_include_directories(md_benchmark
  PRIVATE
  gen/
  ${CMAKE_CURRENT_BINARY_DIR}
  ../ext/xtc/
)

# Timings of the certified GTO critical point search, CPU and GPU, per phase (see bench_topo_gto.c).
# Its own executable: it prints tables rather than ubench statistics, and takes options.
add_executable(md_topo_gto_bench bench_topo_gto.c)
target_compile_definitions(md_topo_gto_bench PRIVATE MD_BENCHMARK_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/../test_data\" ${MD_DEFINES})
target_compile_options(md_topo_gto_bench PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
target_link_options(md_topo_gto_bench PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
target_link_libraries(md_topo_gto_bench PRIVATE mdlib)
add_dependencies(md_topo_gto_bench md_bench_rev)
target_include_directories(md_topo_gto_bench PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

# The GPU evaluations on the shared GEMM core (the density, ELF and LOL, the orbitals with their shell kernel, the dmat
# build), as viamd evaluates them, each checked against the CPU at sampled voxels (see bench_gto_gemm.c). Its own
# executable, which takes options.
if (MD_ENABLE_HDF5 AND MD_ENABLE_GPU)
    add_executable(md_gto_gemm_bench bench_gto_gemm.c)
    target_compile_definitions(md_gto_gemm_bench PRIVATE MD_BENCHMARK_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/../test_data\" ${MD_DEFINES})
    target_compile_options(md_gto_gemm_bench PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
    target_link_options(md_gto_gemm_bench PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
    target_link_libraries(md_gto_gemm_bench PRIVATE mdlib)
    add_dependencies(md_gto_gemm_bench md_bench_rev)
    target_include_directories(md_gto_gemm_bench PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
endif()

# What a density matrix kept in memory in 16 x 16 tiles would hold: the patterns of its consumers (the GPU density,
# ELF/LOL, the certified topology, a D-weighted bound) and the ESP's error with D cut to them (see bench_dmat_census.c).
# CPU only. Its own executable, which takes options.
if (MD_ENABLE_HDF5)
    add_executable(md_dmat_census_bench bench_dmat_census.c)
    target_compile_definitions(md_dmat_census_bench PRIVATE MD_BENCHMARK_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/../test_data\" ${MD_DEFINES})
    target_compile_options(md_dmat_census_bench PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
    target_link_options(md_dmat_census_bench PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
    target_link_libraries(md_dmat_census_bench PRIVATE mdlib)
    add_dependencies(md_dmat_census_bench md_bench_rev)
    target_include_directories(md_dmat_census_bench PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
endif()

# How structure files (PDB, GRO) load over threads: best and median per thread count, the speedup, and a hash of the
# system each count loads (see bench_load.c). Its own executable, which takes the files to load.
add_executable(md_load_bench bench_load.c bench_common.h)
target_compile_options(md_load_bench PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
target_link_options(md_load_bench PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
target_link_libraries(md_load_bench PRIVATE mdlib $<$<PLATFORM_ID:Windows>:psapi>)
add_dependencies(md_load_bench md_bench_rev)
target_include_directories(md_load_bench PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

# What md_parse_f64 costs per number, by the shape the number is written in, against strtod on the same text (see
# bench_parse.c). Its own executable, which prints a table and takes options.
add_executable(md_parse_bench bench_parse.c)
target_compile_options(md_parse_bench PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
target_link_options(md_parse_bench PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
target_link_libraries(md_parse_bench PRIVATE mdlib)
add_dependencies(md_parse_bench md_bench_rev)
target_include_directories(md_parse_bench PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

# How the topology of a system is inferred, and over threads: md_system_infer (or md_bond_infer alone), best and median
# per thread count, the speedup, the page faults and a hash of what it inferred (see bench_infer.c). Its own executable,
# which takes the files to infer.
add_executable(md_infer_bench bench_infer.c bench_common.h)
target_compile_options(md_infer_bench PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
target_link_options(md_infer_bench PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
target_link_libraries(md_infer_bench PRIVATE mdlib $<$<PLATFORM_ID:Windows>:psapi>)
add_dependencies(md_infer_bench md_bench_rev)
target_include_directories(md_infer_bench PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_compile_definitions(md_infer_bench PRIVATE MD_BENCHMARK_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/../test_data\" ${MD_DEFINES})
