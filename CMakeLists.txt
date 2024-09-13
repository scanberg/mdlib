cmake_minimum_required(VERSION 3.20)

include(cmake/common.cmake)

project(mdlib VERSION 0.1.12 LANGUAGES C)

include(CMakeDependentOption)
include(cmake/FindAVX.cmake)

# Options
option(MD_UNITTEST "Enable unittests" ON)
option(MD_BENCHMARK "Enable benchmarks" ON)
option(MD_LINK_STDLIB_STATIC "Link against static stdlibs" ON)
option(MD_ENABLE_AVX  "Enable AVX extensions" ${HAVE_AVX_EXTENSIONS})
option(MD_ENABLE_AVX2 "Enable AVX2 extensions" ${HAVE_AVX2_EXTENSIONS})
option(MD_ENABLE_AVX512 "Enable AVX512 extensions" ${HAVE_AVX512_EXTENSIONS})
cmake_dependent_option(MD_ENABLE_FMA "Enable FMA extensions" ON "HAVE_FMA_EXTENSIONS" OFF)
option(MD_FORCE_ASSERTIONS "Force assertions" OFF)
option(MD_ENABLE_VLX "Enable Veloxchem Support" OFF)
set(MD_GL_SPLINE_SUBDIVISION_COUNT "8" CACHE STRING "Number of subdivision segments for splines used in secondary structures")

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
if (CMAKE_COMPILER_IS_GNUCC)
  set(CMAKE_C_EXTENSIONS ON)
endif()

if (${MD_ENABLE_VLX})
    set(MD_ENABLE_HDF5 TRUE)
endif()

if (APPLE)
    # Remove once we have proper support for Apple Silicon?
    set(CMAKE_OSX_ARCHITECTURES "x86_64")
endif()

file(GLOB CORE_FILES src/core/*.h src/core/*.c src/core/*.inl)
file(GLOB_RECURSE EXT_FILES ext/*)
file(GLOB_RECURSE PROJ_FILES cmake/*)

set(SRC_FILES
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
    src/md_molecule.c
    src/md_molecule.h
    src/md_pdb.c
    src/md_pdb.h
    src/md_script.c
    src/md_script.h
    src/md_script_functions.inl
    src/md_smiles.c
    src/md_smiles.h
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

set(SHADER_FILES 
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

if (${MD_ENABLE_VLX})
    set (SRC_FILES ${SRC_FILES} src/md_vlx.c src/md_vlx.h)
    set (MD_DEFINES ${MD_DEFINES} MD_VLX MD_VLX_BASIS_FOLDER="${CMAKE_CURRENT_SOURCE_DIR}/ext/veloxchem/basis")
endif()

if (MD_ENABLE_NETCDF)
    find_package(NetCDF REQUIRED)
    set(MD_LIBS ${MD_LIBS} NetCDF::NetCDF_C)
endif()

if (MD_ENABLE_HDF5)
    set(HDF5_USE_STATIC_LIBRARIES TRUE)
    find_package(HDF5 REQUIRED)
    set(MD_LIBS ${MD_LIBS} HDF5::HDF5)
endif()

create_resources("${SHADER_FILES}" "gen/gl_shaders.inl")
create_resources("src/shaders/eval_gto.comp" "gen/gto_shaders.inl")

if (CMAKE_C_COMPILER_ID MATCHES "Clang")
  set(MD_FLAGS_EXT ${MD_FLAGS_EXT} "-w")
  set(MD_FLAGS ${MD_FLAGS}
    "-Wall"
    "-Wextra"
    "-Wpedantic"
    "-Wno-missing-field-initializers"
    "-Wno-missing-braces"
    "-Wno-unused-function"
    "-Wno-unused-parameter"
    "-Wno-gnu-zero-variadic-macro-arguments"
    "-Wno-nested-anon-types")

    if (MD_ENABLE_AVX512)
        set(MD_FLAGS_AVX ${MD_FLAGS_AVX} "-mavx512f" "-mavx512dq")
    elseif (MD_ENABLE_AVX2)
        set(MD_FLAGS_AVX ${MD_FLAGS_AVX} "-mavx2")
    elseif (MD_ENABLE_AVX)
        set(MD_FLAGS_AVX ${MD_FLAGS_AVX} "-mavx")
    endif()

    if (MD_ENABLE_FMA)
        set(MD_FLAGS_AVX ${MD_FLAGS_AVX} "-mfma" "-ffp-contract=fast")
    endif()

  if (CMAKE_C_COMPILER_ID STREQUAL "Clang")
    if (CMAKE_C_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC") # clang-cl
    message(STATUS "mdlib: MSVC-Clang configuration was used")
    elseif (CMAKE_C_COMPILER_FRONTEND_VARIANT STREQUAL "GNU") # clang native
    message(STATUS "mdlib: Clang configuration was used")
        set(MD_LIBS ${MD_LIBS} "m" "atomic")
    endif()
  elseif (CMAKE_C_COMPILER_ID STREQUAL "AppleClang")
    message(STATUS "mdlib: AppleClang configuration was used")
  endif()
elseif (CMAKE_C_COMPILER_ID STREQUAL "GNU") # GCC
    message(STATUS "mdlib: GCC configuration was used")
    set(MD_FLAGS_EXT ${MD_FLAGS_EXT} "-w")
    set(MD_FLAGS ${MD_FLAGS}
      "-Wall"
      "-Wextra"
      "-Wpedantic"
      "-Wno-missing-field-initializers"
      "-Wno-missing-braces"
      "-Wno-unused-function"
      "-Wno-unused-parameter"
      "-Wno-format-truncation")

    if (MD_ENABLE_AVX512)
        set(MD_FLAGS_AVX ${MD_FLAGS_AVX} "-mavx512f" "-mavx512dq")
    elseif (MD_ENABLE_AVX2)
        set(MD_FLAGS_AVX ${MD_FLAGS_AVX} "-mavx2")
    elseif (MD_ENABLE_AVX)
        set(MD_FLAGS_AVX ${MD_FLAGS_AVX} "-mavx")
    endif()

    if (MD_ENABLE_FMA)
        set(MD_FLAGS_AVX ${MD_FLAGS_AVX} "-mfma" "-ffast-math")
    endif()

    if (MD_LINK_STDLIB_STATIC)
      set(MD_LIBS ${MD_LIBS}
        "-static-libgcc"
        "-static-libstdc++")
    endif()
    find_package(Threads REQUIRED)
    set(MD_DEFINES ${MD_DEFINES} "${_GNU_SOURCE}")
    set(MD_LIBS ${MD_LIBS}
        "m"
        "Threads::Threads")
elseif ((CMAKE_C_COMPILER_ID STREQUAL "MSVC")) # msvc
    message(STATUS "mdlib: MSVC configuration was used")
    set(MD_FLAGS_EXT ${MD_FLAGS_EXT} "/W2")
    set(MD_FLAGS ${MD_FLAGS}
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
        set(MD_FLAGS_AVX ${MD_FLAGS_AVX} "/arch:AVX512" "/fp:fast")
    elseif (MD_ENABLE_AVX2)
        set(MD_FLAGS_AVX ${MD_FLAGS_AVX} "/arch:AVX2" "/fp:fast")
    elseif (MD_ENABLE_AVX)
        set(MD_FLAGS_AVX ${MD_FLAGS_AVX} "/arch:AVX" "/fp:fast")
    endif()
      
    set(MD_DEFINES ${MD_DEFINES} "_CRT_SECURE_NO_WARNINGS")
    
    set(MD_FLAGS_DEB ${MD_FLAGS_DEB} "/Ob1" "/ZI")
    set(MD_FLAGS_REL ${MD_FLAGS_REL} "/GS-" "/GL")
    set(MD_LINK_FLAGS_DEB ${MD_LINK_FLAGS_DEB} "/INCREMENTAL")
    set(MD_LINK_FLAGS_REL ${MD_LINK_FLAGS_REL} "/LTCG")
    set(PROJ_FILES ${PROJ_FILES} "mdlib.natvis")
endif()

set(MD_FLAGS ${MD_FLAGS} ${MD_FLAGS_AVX})
set(MD_PROPERTIES)

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
source_group("shaders"  FILES ${SHADER_FILES})
source_group("proj"     FILES ${PROJ_FILES})

add_library(mdlib ${CORE_FILES} ${SRC_FILES} ${EXT_FILES} ${SHADER_FILES} ${PROJ_FILES})

if (MD_LINK_STDLIB_STATIC)
    set(MD_PROPERTIES ${MD_PROPERTIES} MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
else()
    set(MD_PROPERTIES ${MD_PROPERTIES} MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()

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
)

set_target_properties(mdlib PROPERTIES ${MD_PROPERTIES})

target_compile_definitions(mdlib PRIVATE ${MD_DEFINES})
target_compile_options(mdlib PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
target_link_options(mdlib PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
set_source_files_properties(${EXT_FILES}, PROPERTIES COMPILE_OPTIONS ${MD_FLAGS_EXT})
target_compile_features(mdlib PRIVATE c_std_11)
target_link_libraries(mdlib PRIVATE ${MD_LIBS})

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
