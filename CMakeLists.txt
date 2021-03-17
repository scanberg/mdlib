cmake_minimum_required(VERSION 3.0)

project(mold VERSION 0.1.1 LANGUAGES C)

option(MOLD_UNITTEST "Enable MOLD unittest" ON)
option(MD_USE_RELATIVE_PATHS "Use relative paths" OFF)

add_subdirectory(src)

if(MOLD_UNITTEST)
  add_subdirectory(unittest)
endif()
