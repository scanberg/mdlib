cmake_minimum_required(VERSION 3.0)

project(libmold VERSION 0.1.1 LANGUAGES C)

source_group("core" FILES core/arena_alloc.h core/bitop.h core/bitop.c)
source_group("src" FILES mold_draw.h mold_draw.c mold_filter.h mold_filter.c mold_util.h mold_util.c mold_molecule.h)
source_group("ext" FILES ext/gl3w.c)

if (MSVC)
    message(STATUS "MSVC configuration was used")
    add_compile_options(/MT /W4 /MP /GR- /fp:fast /MP /arch:AVX)# /arch:AVX2
endif()

add_library(mold mold_draw.h mold_draw.c mold_filter.h mold_filter.c mold_molecule.h mold_util.h mold_util.c ext/gl3w.c)

target_include_directories(mold
	PUBLIC
		$<INSTALL_INTERFACE:
			${CMAKE_CURRENT_SOURCE_DIR}/mold_draw.h>
		$<BUILD_INTERFACE:
			${CMAKE_CURRENT_SOURCE_DIR}/mold_draw.h>
)

target_compile_definitions(mold PRIVATE
	MOLD_SHADER_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/shaders\")

