cmake_minimum_required(VERSION 3.0)

option(MOLD_UNITTEST "Enable MOLD unittest" FALSE)

add_subdirectory(src)
if (MOLD_UNITTEST)
	add_subdirectory(unittest)
endif()