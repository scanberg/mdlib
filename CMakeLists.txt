cmake_minimum_required(VERSION 3.0)

project(mold VERSION 0.1.1 LANGUAGES C)

option(MOLD_UNITTEST "Enable MOLD unittest" ON)
option(MD_USE_RELATIVE_PATHS "Use relative paths" OFF)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_subdirectory(src)

if(MOLD_UNITTEST)
  add_subdirectory(unittest)
endif()
