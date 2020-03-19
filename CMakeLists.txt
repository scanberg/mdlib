cmake_minimum_required(VERSION 3.0)

project(libmold VERSION 0.1.1 LANGUAGES C)

source_group("src" FILES mold.h mold.c)
add_library(mold mold.h mold.c)

target_include_directories(mold
	PUBLIC
		$<INSTALL_INTERFACE:
			${CMAKE_CURRENT_SOURCE_DIR};>
		$<BUILD_INTERFACE:
			${CMAKE_CURRENT_SOURCE_DIR};>
)

target_compile_definitions(mold PRIVATE
	MOLD_SHADER_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/shaders\")