add_executable(md_unittest unittest.c)
target_link_libraries(md_unittest PRIVATE mold)
target_compile_definitions(md_unittest PRIVATE MD_UNITTEST_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/data\")

add_test(NAME mytest
	COMMAND $<TARGET_FILE:md_unittest>)
