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

# GPU GTO kernel comparison (densities and molecular orbitals). Runs from
# benchmark/data/density alone; with HDF5 it also uses the SCF density and orbitals of
# test_data/vlx/mol.h5 when that file is present.
if (MD_ENABLE_GPU)
    # mdlib revision for the log (taken at configure time; the bench scripts reconfigure).
    set(MD_BENCH_GIT_REV "unknown")
    find_package(Git QUIET)
    if (GIT_FOUND)
        execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse --short=10 HEAD
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        OUTPUT_VARIABLE _rev OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET RESULT_VARIABLE _rev_res)
        if (_rev_res EQUAL 0 AND _rev)
            execute_process(COMMAND ${GIT_EXECUTABLE} status --porcelain --untracked-files=no
                            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                            OUTPUT_VARIABLE _dirty OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
            set(MD_BENCH_GIT_REV "${_rev}")
            if (_dirty)
                set(MD_BENCH_GIT_REV "${_rev} (with local changes)")
            endif()
        endif()
    endif()

    add_executable(md_bench_gto_gpu bench_gto_gpu.c)
    target_compile_definitions(md_bench_gto_gpu PRIVATE
        MD_BENCHMARK_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/../test_data\"
        MD_DENSITY_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/data/density\"
        MD_BENCH_GIT_REV=\"${MD_BENCH_GIT_REV}\"
        MD_BENCH_BUILD_TYPE=\"$<CONFIG>\"
        ${MD_DEFINES})
    target_compile_options(md_bench_gto_gpu PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
    target_link_options(md_bench_gto_gpu PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
    target_link_libraries(md_bench_gto_gpu PRIVATE mdlib)
    if (WIN32)
        target_link_libraries(md_bench_gto_gpu PRIVATE advapi32)
    endif()
endif()
