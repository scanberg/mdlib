cmake_minimum_required(VERSION 3.0)

project(md_unittest VERSION 0.1.1 LANGUAGES C)
add_executable(md_unittest unittest.c)
target_link_libraries(md_unittest mold)
target_compile_definitions(md_unittest PRIVATE MD_UNITTEST_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/data\")

set_property(TARGET md_unittest PROPERTY C_STANDARD 99)
set_property(TARGET md_unittest PROPERTY C_STANDARD_REQUIRED ON)

add_test(NAME mytest
         COMMAND testDriver --config $<CONFIG>
                            --exe $<TARGET_FILE:md_unittest>)