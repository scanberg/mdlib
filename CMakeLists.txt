cmake_minimum_required(VERSION 3.0)

project(mold)

option(MOLD_UNITTEST "Enable MOLD unittest" TRUE)

add_subdirectory(src)
if (MOLD_UNITTEST)
	add_subdirectory(unittest)
endif()