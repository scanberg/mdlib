set(SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/unittest.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_allocator.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_array.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_bitop.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_intrinsics.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_script.c
)

add_executable(md_unittest ${SOURCE_FILES})
target_link_libraries(md_unittest PRIVATE mold)
target_compile_definitions(md_unittest PRIVATE MD_UNITTEST_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/data\")

add_test(NAME mytest
	COMMAND $<TARGET_FILE:md_unittest>)
