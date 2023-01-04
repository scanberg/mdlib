cmake_minimum_required(VERSION 3.15)

project(mdlib VERSION 0.1.11 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
if(CMAKE_COMPILER_IS_GNUCC)
  set(CMAKE_C_EXTENSIONS ON)
endif()

include(CMakeDependentOption)
include(cmake/common.cmake)
include(cmake/FindAVX.cmake)

# Options
option(MD_UNITTEST "Enable unittest" ON)
option(MD_LINK_STDLIB_STATIC "Link against static stdlibs" ON)
option(MD_USE_RELATIVE_RESOURCE_PATH "Use relative path for resources, this will copy the relevant files into the binary folder upon building" ON)
option(MD_ENABLE_AVX  "Enable AVX extensions" ${HAVE_AVX_EXTENSIONS})
option(MD_ENABLE_AVX2 "Enable AVX2 extensions" ${HAVE_AVX2_EXTENSIONS})
cmake_dependent_option(MD_ENABLE_FMA "Enable FMA extensions" ON "HAVE_FMA_EXTENSIONS" OFF)

if(MD_UNITTEST)
  enable_testing()
  add_subdirectory(unittest)
endif()

file(GLOB SRC_FILES src/*.h src/*.c src/*.inl)
file(GLOB CORE_FILES src/core/*.h src/core/*.c src/core/*.inl)
file(GLOB_RECURSE EXT_FILES ext/*)
file(GLOB_RECURSE SHADER_FILES src/shaders/*)
file(GLOB_RECURSE PROJ_FILES cmake/*)

set(MD_DEFINES)
if (${MD_ENABLE_FMA})
    set(MD_DEFINES ${MD_DEFINES} "__FMA__")
endif()
set(MD_FLAGS_EXT)
set(MD_FLAGS)
set(MD_FLAGS_REL)
set(MD_FLAGS_DEB)
if (MD_USE_RELATIVE_RESOURCE_PATH)
  set(MD_SHADER_DIR "shaders")
  set(MD_SHADER_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/shaders")
  set(MD_SHADER_DST_DIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/shaders")
else()
  set(MD_SHADER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/shaders")
endif()
set(MD_LIBS)

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

    if (MD_ENABLE_AVX2)
        set(MD_FLAGS ${MD_FLAGS} "-mavx2")
    elseif (MD_ENABLE_AVX)
		set(MD_FLAGS ${MD_FLAGS} "-mavx")
    endif()

    if (MD_ENABLE_FMA)
		set(MD_FLAGS ${MD_FLAGS} "-mfma" "-ffp-contract=fast")
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
      "-Wno-unused-parameter")

    if (MD_ENABLE_AVX2)
        set(MD_FLAGS ${MD_FLAGS} "-mavx2")
    elseif (MD_ENABLE_AVX)
		set(MD_FLAGS ${MD_FLAGS} "-mavx")
    endif()

    if (MD_ENABLE_FMA)
		set(MD_FLAGS ${MD_FLAGS} "-mfma" "-ffast-math")
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
    set(MD_SHADER_DST_DIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIG>/shaders")
    set(MD_FLAGS_EXT ${MD_FLAGS_EXT} "/W2")
    set(MD_FLAGS ${MD_FLAGS}
      "/W4"
      "/wd4201"
      "/wd4305"
      "/wd4324"
      "/wd26451"
      "/MP"
      "/GR-"
      "/Oi")

    if (MD_ENABLE_AVX2)
        set(MD_FLAGS ${MD_FLAGS} "/arch:AVX2" "/fp:fast")
    elseif (MD_ENABLE_AVX)
		set(MD_FLAGS ${MD_FLAGS} "/arch:AVX" "/fp:fast")
    endif()
      
    set(MD_DEFINES ${MD_DEFINES} "_CRT_SECURE_NO_WARNINGS")
    
    set(MD_FLAGS_DEB ${MD_FLAGS_DEB}
        "/ZI")
    set(MD_FLAGS_REL ${MD_FLAGS_REL}
        "/GS-")
    set(PROJ_FILES ${PROJ_FILES} "mdlib.natvis")
endif()

source_group("src" FILES ${SRC_FILES}) 
source_group("core" FILES ${CORE_FILES})
source_group("ext" FILES ${EXT_FILES})
source_group("shaders" FILES ${SHADER_FILES})
source_group("proj" FILES ${PROJ_FILES})

add_library(mdlib ${CORE_FILES} ${SRC_FILES} ${EXT_FILES} ${SHADER_FILES} ${PROJ_FILES})

if (MD_LINK_STDLIB_STATIC)
    set_property(TARGET mdlib PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
else()
    set_property(TARGET mdlib PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()

if (MD_USE_RELATIVE_RESOURCE_PATH)
  create_copy_resource_dir_target(md_copy_shaders ${MD_SHADER_SRC_DIR} ${MD_SHADER_DST_DIR})
  add_dependencies(mdlib md_copy_shaders) 
endif()

target_include_directories(mdlib PUBLIC
  src/
  ext/dcd/
  ext/fastlz/
  ext/gl3w/
  ext/svd3/
  ext/xtc/)

target_compile_definitions(mdlib PUBLIC MD_SHADER_DIR="${MD_SHADER_DIR}" PRIVATE ${MD_DEFINES})
target_compile_options(mdlib PRIVATE ${MD_FLAGS} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
set_source_files_properties(${EXT_FILES}, PROPERTIES COMPILE_OPTIONS ${MD_FLAGS_EXT})
target_compile_features(mdlib PRIVATE c_std_11)
target_link_libraries(mdlib PRIVATE ${MD_LIBS})

get_directory_property(has_parent PARENT_DIRECTORY)
if (has_parent)
    set(MD_FLAGS ${MD_FLAGS} PARENT_SCOPE)
    set(MD_FLAGS_DEB ${MD_FLAGS_DEB} PARENT_SCOPE)
    set(MD_FLAGS_REL ${MD_FLAGS_REL} PARENT_SCOPE)
endif()