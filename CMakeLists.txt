cmake_minimum_required(VERSION 3.20)

include(cmake/CheckVCPKG.cmake)

project(mdlib LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
if (CMAKE_COMPILER_IS_GNUCC)
  set(CMAKE_C_EXTENSIONS ON)
endif()

include(CMakeDependentOption)
include(cmake/FindAVX.cmake)
include(cmake/Common.cmake)

# Options
option(MD_UNITTEST "Enable unittests" ON)
option(MD_BENCHMARK "Enable benchmarks" ON)
option(MD_LINK_STDLIB_STATIC "Link against static stdlibs" ON)
option(MD_ENABLE_GPU "Enable GPU acceleration" OFF)
option(MD_ENABLE_AVX  "Enable AVX extensions" ${HAVE_AVX_EXTENSIONS})
option(MD_ENABLE_AVX2 "Enable AVX2 extensions" ${HAVE_AVX2_EXTENSIONS})
option(MD_ENABLE_AVX512 "Enable AVX512 extensions" ${HAVE_AVX512_EXTENSIONS})
cmake_dependent_option(MD_ENABLE_FMA "Enable FMA extensions" ON "HAVE_FMA_EXTENSIONS" OFF)
option(MD_FORCE_ASSERTIONS "Force assertions" OFF)
option(MD_ENABLE_VLX "Enable Veloxchem Support" OFF)
set(MD_GL_SPLINE_SUBDIVISION_COUNT "8" CACHE STRING "Number of subdivision segments for splines used in secondary structures")

if (MD_ENABLE_VLX)
  set(MD_ENABLE_HDF5 TRUE)
endif()

#file(GLOB CORE_FILES src/core/*.h src/core/*.c src/core/*.inl)
file(GLOB_RECURSE EXT_FILES ext/*)
file(GLOB_RECURSE PROJ_FILES cmake/*)

set(CORE_FILES
  src/core/md_allocator.c
  src/core/md_allocator.h
  src/core/md_arena_allocator.c
  src/core/md_arena_allocator.h
  src/core/md_array.h
  src/core/md_base64.c
  src/core/md_base64.h
  src/core/md_bitfield.c
  src/core/md_bitfield.h
  src/core/md_bitop.inl
  src/core/md_common.c
  src/core/md_common.h
  src/core/md_compiler.h
  src/core/md_gl_util.c
  src/core/md_gl_util.h
  src/core/md_gpu.h
  src/core/md_grid.c
  src/core/md_grid.h
  src/core/md_handle.c
  src/core/md_handle.h
  src/core/md_hash.h
  src/core/md_intrinsics.h
  src/core/md_linear_allocator.c
  src/core/md_linear_allocator.h
  src/core/md_log.c
  src/core/md_log.h
  src/core/md_lru_cache.inl
  src/core/md_os.c
  src/core/md_os.h
  src/core/md_parse.h
  src/core/md_platform.h
  src/core/md_ring_allocator.c
  src/core/md_ring_allocator.h
  src/core/md_simd.h
  src/core/md_spatial_acc.c
  src/core/md_spatial_acc.h
  src/core/md_spatial_hash.c
  src/core/md_spatial_hash.h
  src/core/md_str_builder.c
  src/core/md_str_builder.h
  src/core/md_str.c
  src/core/md_str.h
  src/core/md_tracking_allocator.c
  src/core/md_tracking_allocator.h
  src/core/md_unit.c
  src/core/md_unit.h
  src/core/md_vec_math.c
  src/core/md_vec_math.h
  src/core/md_virtual_allocator.c
  src/core/md_virtual_allocator.h
)

set(SRC_FILES
  src/md_atomic_infer.c
  src/md_csv.c
  src/md_csv.h
  src/md_cube.c
  src/md_cube.h
  src/md_edr.c
  src/md_edr.h
  src/md_filter.h
  src/md_frame_cache.c
  src/md_frame_cache.h
  src/md_gfx.c
  src/md_gfx.h
  src/md_gl.c
  src/md_gl.h
  src/md_gro.c
  src/md_gro.h
  src/md_gto.c
  src/md_gto.h
  src/md_lammps.c
  src/md_lammps.h
  src/md_mmcif.c
  src/md_mmcif.h
  src/md_system.c
  src/md_system.h
  src/md_pdb.c
  src/md_pdb.h
  src/md_script.c
  src/md_script.h
  src/md_script_functions.inl
  src/md_smiles.c
  src/md_smiles.h
  src/md_topo.c
  src/md_topo.h
  src/md_trajectory.h
  src/md_trr.c
  src/md_trr.h
  src/md_types.h
  src/md_util.c
  src/md_util.h
  src/md_xtc.c
  src/md_xtc.h
  src/md_xvg.c
  src/md_xvg.h
  src/md_xyz.c
  src/md_xyz.h
)

set(GL_SHADER_FILES 
  src/shaders/compute_spline.vert
  src/shaders/compute_spline.geom
  src/shaders/compute_velocity.vert
  src/shaders/cartoon.vert
  src/shaders/cartoon.geom
  src/shaders/cartoon.frag
  src/shaders/licorice.vert
  src/shaders/licorice.geom
  src/shaders/licorice.frag
  src/shaders/ribbons.vert
  src/shaders/ribbons.geom
  src/shaders/ribbons.frag
  src/shaders/spacefill.vert
  src/shaders/spacefill.geom
  src/shaders/spacefill.frag
)

set(GTO_SHADER_FILES 
  src/shaders/eval_alie.comp
  src/shaders/eval_gto.comp
  src/shaders/eval_gto_density.comp
  src/shaders/eval_gto_density_grad.comp
  src/shaders/voronoi_segment.comp
)

set(TOPO_SHADER_FILES
  src/shaders/topo/critical_point_compaction.comp
  src/shaders/topo/bidirectional_manifold.comp
  src/shaders/topo/path_compression.comp
  src/shaders/topo/critical_points.comp
  src/shaders/topo/vertex_edge_extraction.comp
)

set(GEN_FILES)

if (MD_ENABLE_GPU)
  if (APPLE)
    set(MD_GPU_BACKEND "METAL")
    enable_language(OBJC)
    list(APPEND CORE_FILES src/core/md_gpu_metal.m)
    set_source_files_properties(
      src/core/md_gpu_metal.m
      PROPERTIES
      LANGUAGE OBJC
    )
  else()
    set(MD_GPU_BACKEND "VULKAN")
  endif()
  message(STATUS "mdlib: GPU backend = ${MD_GPU_BACKEND}")
endif()

include(cmake/Shaders.cmake)

 # Generate gen folder if not present
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/gen)

#md_compile_shaders(gen/topo_shaders.inl
#    NAMESPACE topo
#    SOURCES ${TOPO_SHADER_FILES}
#)

list(APPEND GEN_FILES ${CMAKE_CURRENT_BINARY_DIR}/gen/topo_shaders.inl)

set(MD_DEFINES MD_GL_SPLINE_SUBDIVISION_COUNT=${MD_GL_SPLINE_SUBDIVISION_COUNT})

if (CMAKE_C_BYTE_ORDER STREQUAL BIG_ENDIAN)
  set(MD_DEFINES ${MD_DEFINES} __BIG_ENDIAN__)
else()
  set(MD_DEFINES ${MD_DEFINES} __LITTLE_ENDIAN__)
endif()

if (${MD_ENABLE_FMA})
  set(MD_DEFINES ${MD_DEFINES} __FMA__)
endif()

if (MD_FORCE_ASSERTIONS)
  set(MD_DEFINES ${MD_DEFINES} __FORCE_ASSERTIONS__=1)
else()
  set(MD_DEFINES ${MD_DEFINES} __FORCE_ASSERTIONS__=0)    
endif()

set(MD_FLAGS)
set(MD_FLAGS_EXT)
set(MD_FLAGS_REL)
set(MD_FLAGS_DEB)
set(MD_FLAGS_AVX)

set(MD_LINK_FLAGS)
set(MD_LINK_FLAGS_REL)
set(MD_LINK_FLAGS_DEB)
set(MD_LIBS)
set(MD_PROPERTIES)

if (MD_ENABLE_VLX)
  list(APPEND SRC_FILES src/md_vlx.c src/md_vlx.h)
  if(NOT DEFINED MD_VLX_BASIS_FOLDER)
    set(MD_VLX_BASIS_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}/ext/veloxchem/basis")
  endif()
  list(APPEND MD_DEFINES MD_VLX MD_VLX_BASIS_FOLDER="${MD_VLX_BASIS_FOLDER}")
endif()

if (MD_ENABLE_NETCDF)
  find_package(NetCDF REQUIRED)
  list(APPEND MD_LIBS NetCDF::NetCDF_C)
endif()

if (MD_ENABLE_HDF5)
  if (WIN32)
    set(HDF5_USE_STATIC_LIBRARIES TRUE)
  endif()
  set(HDF5_PREFER_PARALLEL FALSE)
  find_package(HDF5 REQUIRED)
  list(APPEND MD_LIBS HDF5::HDF5)
endif()

create_resources("${GL_SHADER_FILES}"  "gen/gl_shaders.inl")
create_resources("${GTO_SHADER_FILES}" "gen/gto_shaders.inl")
create_resources("${TOPO_SHADER_FILES}" "gen/topo_shaders.inl")

# Detect AppleClang vs other Clang
if (CMAKE_C_COMPILER_ID STREQUAL "Clang")
  if (CMAKE_C_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")  # clang-cl
    set(COMPILER_TYPE "MSVC-Clang")
  elseif (CMAKE_C_COMPILER_FRONTEND_VARIANT STREQUAL "AppleClang")
    set(COMPILER_TYPE "AppleClang")
  elseif (CMAKE_C_COMPILER_FRONTEND_VARIANT STREQUAL "GNU")  # native Clang
    set(COMPILER_TYPE "Clang")
  else()
    set(COMPILER_TYPE "Clang-Unknown")
  endif()
elseif (CMAKE_C_COMPILER_ID STREQUAL "AppleClang")  # CMake 3.x
  set(COMPILER_TYPE "AppleClang")
elseif (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  set(COMPILER_TYPE "GCC")
elseif (CMAKE_C_COMPILER_ID STREQUAL "MSVC")
  set(COMPILER_TYPE "MSVC")
else()
  set(COMPILER_TYPE "Unknown")
endif()

message(STATUS "mdlib: Compiler detected as ${COMPILER_TYPE}")

if (COMPILER_TYPE MATCHES "Clang")
  list(APPEND MD_FLAGS_EXT "-w")
  list(APPEND MD_FLAGS
    "-Wall"
    "-Wextra"
    "-Wpedantic"
    "-Wno-missing-field-initializers"
    "-Wno-missing-braces"
    "-Wno-unused-function"
    "-Wno-unused-parameter"
    "-Wno-gnu-zero-variadic-macro-arguments"
    "-Wno-nested-anon-types")

  if (COMPILER_TYPE STREQUAL "Clang") # native Clang
    list(APPEND MD_LIBS "m" "atomic")
  endif()

  if (COMPILER_TYPE STREQUAL "AppleClang")
    list(APPEND MD_FLAGS "-Wno-format")
  endif()

  # AVX/FMA flags
  if (MD_ENABLE_AVX512)
    list(APPEND MD_FLAGS_AVX "-mavx512f" "-mavx512dq")
  elseif (MD_ENABLE_AVX2)
    list(APPEND MD_FLAGS_AVX "-mavx2")
  elseif (MD_ENABLE_AVX)
    list(APPEND MD_FLAGS_AVX "-mavx")
  endif()

  if (MD_ENABLE_FMA)
    list(APPEND MD_FLAGS_AVX "-mfma" "-ffp-contract=fast")
  endif()

elseif (COMPILER_TYPE STREQUAL "GCC")
    list(APPEND MD_FLAGS_EXT "-w")
    list(APPEND MD_FLAGS
      "-Wall"
      "-Wextra"
      "-Wpedantic"
      "-Wno-missing-field-initializers"
      "-Wno-missing-braces"
      "-Wno-unused-function"
      "-Wno-unused-parameter"
      "-Wno-format-truncation")

    if (MD_ENABLE_AVX512)
        list(APPEND MD_FLAGS_AVX "-mavx512f" "-mavx512dq")
    elseif (MD_ENABLE_AVX2)
        list(APPEND MD_FLAGS_AVX "-mavx2")
    elseif (MD_ENABLE_AVX)
        list(APPEND MD_FLAGS_AVX "-mavx")
    endif()

    if (MD_ENABLE_FMA)
        list(APPEND MD_FLAGS_AVX "-mfma" "-ffast-math")
    endif()

    if (MD_LINK_STDLIB_STATIC)
      list(APPEND MD_LIBS
        "-static-libgcc"
        "-static-libstdc++")
    endif()
    find_package(Threads REQUIRED)
    list(APPEND MD_DEFINES "${_GNU_SOURCE}")
    list(APPEND MD_LIBS
        "m"
        "Threads::Threads")
elseif ((COMPILER_TYPE STREQUAL "MSVC"))
    list(APPEND MD_FLAGS_EXT "/W2")
    list(APPEND MD_FLAGS
      "/W4"
      "/wd4201"
      "/wd4305"
      "/wd4324"
      "/wd26451"
      "/w44255"
      "/MP"
      "/GR-"
      "/Oi"
      "/fp:except-")

    if (MD_ENABLE_AVX512)
        list(APPEND MD_FLAGS_AVX "/arch:AVX512" "/fp:fast")
    elseif (MD_ENABLE_AVX2)
        list(APPEND MD_FLAGS_AVX "/arch:AVX2" "/fp:fast")
    elseif (MD_ENABLE_AVX)
        list(APPEND MD_FLAGS_AVX "/arch:AVX" "/fp:fast")
    endif()
      
    list(APPEND MD_DEFINES "_CRT_SECURE_NO_WARNINGS")
    
    list(APPEND MD_FLAGS_DEB "/Ob1" "/ZI")
    list(APPEND MD_FLAGS_REL "/GS-" "/GL")
    list(APPEND MD_LINK_FLAGS_DEB "/INCREMENTAL")
    list(APPEND MD_LINK_FLAGS_REL "/LTCG")
    list(APPEND PROJ_FILES "mdlib.natvis")

    if (MD_LINK_STDLIB_STATIC)
        list(APPEND MD_PROPERTIES MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    else()
        list(APPEND MD_PROPERTIES MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
    endif()
endif()

set(MD_FLAGS ${MD_FLAGS} ${MD_FLAGS_AVX})

if (MD_UNITTEST)
  enable_testing()
  add_subdirectory(unittest)
endif()

if (MD_BENCHMARK)
  add_subdirectory(benchmark)
endif()

source_group("src"      FILES ${SRC_FILES}) 
source_group("src/core" FILES ${CORE_FILES})
source_group("ext"      FILES ${EXT_FILES})
source_group("gen"      FILES ${GEN_FILES})
source_group("shaders"  FILES ${SHADER_FILES})
source_group("proj"     FILES ${PROJ_FILES})

add_library(mdlib ${CORE_FILES} ${SRC_FILES} ${EXT_FILES} ${GEN_FILES} ${SHADER_FILES} ${PROJ_FILES})

target_include_directories(mdlib
  PUBLIC
  src/
  ext/gl3w/
  ext/simde/
  ext/xxhash/

  PRIVATE
  gen/
  ext/fastlz/
  ext/dcd/
  ext/xtc/
  ext/svd3/
  ext/stb/
  ext/libdivide
)

if (MD_PROPERTIES)
  set_target_properties(mdlib PROPERTIES ${MD_PROPERTIES})
endif()

target_compile_definitions(mdlib PRIVATE ${MD_DEFINES})
target_compile_options(mdlib PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
target_link_options(mdlib PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
set_source_files_properties(${EXT_FILES}, PROPERTIES COMPILE_OPTIONS ${MD_FLAGS_EXT})
target_compile_features(mdlib PRIVATE c_std_11)
target_link_libraries(mdlib PUBLIC ${MD_LIBS})

get_directory_property(has_parent PARENT_DIRECTORY)
if (has_parent)
    set(MD_DEFINES ${MD_DEFINES} PARENT_SCOPE)
    set(MD_FLAGS ${MD_FLAGS} PARENT_SCOPE)
    set(MD_FLAGS_DEB ${MD_FLAGS_DEB} PARENT_SCOPE)
    set(MD_FLAGS_REL ${MD_FLAGS_REL} PARENT_SCOPE)
    set(MD_FLAGS_EXT ${MD_FLAGS_EXT} PARENT_SCOPE)
    set(MD_FLAGS_AVX ${MD_FLAGS_AVX} PARENT_SCOPE)
    set(MD_LINK_FLAGS ${MD_LINK_FLAGS} PARENT_SCOPE)
    set(MD_LINK_FLAGS_DEB ${MD_LINK_FLAGS_DEB} PARENT_SCOPE)
    set(MD_LINK_FLAGS_REL ${MD_LINK_FLAGS_REL} PARENT_SCOPE)
    set(MD_ENABLE_VLX ${MD_ENABLE_VLX} PARENT_SCOPE)
endif()
