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

# md_render: how well clusters bound their atoms on real systems (CPU), and what md_render costs on
# generated scenes (GPU) (see bench_render.c). Its own executable, which prints tables and takes options.
add_executable(md_render_bench bench_render.c)
target_compile_definitions(md_render_bench PRIVATE MD_BENCHMARK_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/../test_data\" ${MD_DEFINES})
target_compile_options(md_render_bench PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
target_link_options(md_render_bench PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
target_link_libraries(md_render_bench PRIVATE mdlib)
add_dependencies(md_render_bench md_bench_rev)
target_include_directories(md_render_bench PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

# md_render_snap: a structure through md_render_mol and md_render, written as PNG images (headless; see
# render_snap.c): for looking at what the renderer draws on real data.
add_executable(md_render_snap render_snap.c render_scene.h)
target_compile_definitions(md_render_snap PRIVATE ${MD_DEFINES})
target_compile_options(md_render_snap PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
target_link_options(md_render_snap PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
target_link_libraries(md_render_snap PRIVATE mdlib)
add_dependencies(md_render_snap md_bench_rev)
target_include_directories(md_render_snap PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

# md_render_window: md_render_snap's scene in a window, presented through md_gpu's surface (see
# render_window.c): whether presentation works end to end on a machine, locally or remotely. Off by
# default, since it fetches GLFW (X11 only on Linux: what remote desktops serve, and Wayland desktops
# run it through XWayland).
option(MD_RENDER_WINDOW "Build md_render_window, a presentation test (fetches GLFW)" OFF)
if (MD_RENDER_WINDOW AND MD_ENABLE_GPU)
  include(FetchContent)
  set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
  set(GLFW_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
  set(GLFW_BUILD_DOCS     OFF CACHE BOOL "" FORCE)
  set(GLFW_INSTALL        OFF CACHE BOOL "" FORCE)
  set(GLFW_BUILD_WAYLAND  OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG        3.4
    GIT_SHALLOW    TRUE
  )
  FetchContent_MakeAvailable(glfw)

  add_executable(md_render_window render_window.c render_scene.h)
  target_compile_definitions(md_render_window PRIVATE ${MD_DEFINES})
  target_compile_options(md_render_window PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
  target_link_options(md_render_window PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
  target_link_libraries(md_render_window PRIVATE mdlib glfw)
  add_dependencies(md_render_window md_bench_rev)
  target_include_directories(md_render_window PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR}/gen)
  compile_gpu_shaders(render_window_shaders.inl
    TARGET    md_render_window
    NAMESPACE md_shader
    SOURCE    ${CMAKE_CURRENT_SOURCE_DIR}/shaders/render_window.slang
    VERTEX    full_vs
    FRAGMENT  compose_fs
  )
endif()

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
