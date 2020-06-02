cmake_minimum_required(VERSION 3.0)

project(mold_unittest VERSION 0.1.1 LANGUAGES C)
add_executable(mold_unittest unittest.c)
target_link_libraries(mold_unittest mold)