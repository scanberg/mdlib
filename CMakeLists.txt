cmake_minimum_required(VERSION 3.5)

project(mdlib VERSION 0.1.2 LANGUAGES C)

INCLUDE(cmake/common.cmake)
INCLUDE(cmake/FindAVX.cmake)

# Options
option(MD_UNITTEST "Enable unittest" ON)
option(MD_LINK_STDLIB_STATIC "Link against static stdlibs" ON)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
if(CMAKE_COMPILER_IS_GNUCC)
  set(CMAKE_C_EXTENSIONS ON)
endif()

add_subdirectory(ext)

if(MD_UNITTEST)
  enable_testing()
  add_subdirectory(unittest)
endif()

file(GLOB SRC_FILES src/*.h src/*.c)
file(GLOB CORE_FILES src/core/*.h src/core/*.c)
file(GLOB_RECURSE SHADER_FILES src/shaders/*)

source_group("src" FILES ${SRC_FILES})
source_group("core" FILES ${CORE_FILES})
source_group("shaders" FILES ${SHADER_FILES})

add_library(mdlib ${CORE_FILES} ${SRC_FILES} ${EXT_FILES} ${SHADER_FILES})

set(MD_WARNING_FLAGS)
set(MD_COMMON_COMPILE_FLAGS ${AVX_FLAGS})
set(MD_RELEASE_COMPILE_FLAGS)
set(MD_DEBUG_COMPILE_FLAGS)
set(MD_COMPILE_DEFINITIONS MD_SHADER_DIR=\"shaders\")
set(MD_LIBS)
set(MD_SHADER_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/shaders")
set(MD_SHADER_DST_DIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/shaders")

if (CMAKE_C_COMPILER_ID MATCHES "Clang")
  list(APPEND MD_WARNING_FLAGS -Wall -Wextra -Wpedantic -Wno-missing-field-initializers -Wno-missing-braces -Wno-unused-function -Wno-unused-parameter -Wno-gnu-anonymous-struct -Wno-nested-anon-types)
  if (CMAKE_C_COMPILER_ID STREQUAL "Clang")
    if (CMAKE_C_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC") # clang-cl
    message(STATUS "mdlib: MSVC-Clang configuration was used")
    elseif (CMAKE_C_COMPILER_FRONTEND_VARIANT STREQUAL "GNU") # clang native
    message(STATUS "mdlib: Clang configuration was used")
    endif()
  elseif (CMAKE_C_COMPILER_ID STREQUAL "AppleClang")
    message(STATUS "mdlib: AppleClang configuration was used")
  endif()
elseif (CMAKE_C_COMPILER_ID STREQUAL "GNU") # GCC
    message(STATUS "mdlib: GCC configuration was used")
    list(APPEND MD_WARNING_FLAGS -Wall -Wextra -Wpedantic -Wno-missing-field-initializers -Wno-missing-braces -Wno-unused-function -Wno-unused-parameter)
    if (MD_LINK_STDLIB_STATIC)
      list(APPEND MD_LIBS -static-libgcc -static-libstdc++)
    endif()
    find_package(Threads REQUIRED)
    list(APPEND MD_COMPILE_DEFINITIONS _GNU_SOURCE)
    list(APPEND MD_LIBS m Threads::Threads)
elseif ((CMAKE_C_COMPILER_ID STREQUAL "MSVC")) # msvc
    message(STATUS "mdlib: MSVC configuration was used")
    set(MD_SHADER_DST_DIR "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/$<CONFIG>/shaders")
    list(APPEND MD_WARNING_FLAGS /W4 /wd4201 /wd4324 /wd26451)
    list(APPEND MD_COMPILE_FLAGS /MP /GR- /fp:fast)
    list(APPEND MD_COMPILE_DEFINITIONS _CRT_SECURE_NO_WARNINGS)
    
    if (MD_LINK_STDLIB_STATIC)
        list(APPEND MD_DEBUG_COMPILE_FLAGS /MTd)
        list(APPEND MD_RELEASE_COMPILE_FLAGS /MT /Oi /GS-)
    else()
        list(APPEND MD_DEBUG_COMPILE_FLAGS /MDd)
        list(APPEND MD_RELEASE_COMPILE_FLAGS /MD /Oi /GS-)
    endif()
endif()

create_copy_resource_dir_target(md_copy_shaders ${MD_SHADER_SRC_DIR} ${MD_SHADER_DST_DIR})

add_dependencies(mdlib md_copy_shaders)

target_include_directories(mdlib PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_compile_definitions(mdlib PRIVATE ${MD_COMPILE_DEFINITIONS})
target_compile_options(mdlib PRIVATE ${MD_WARNING_FLAGS} ${MD_COMMON_COMPILE_FLAGS})
add_compile_options_config(DEBUG ${MD_DEBUG_COMPILE_FLAGS})
add_compile_options_config(RELEASE ${MD_RELEASE_COMPILE_FLAGS})
target_compile_features(mdlib PRIVATE c_std_11)
target_link_libraries(mdlib PRIVATE mdext ${MD_LIBS})
