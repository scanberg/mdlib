add_executable(md_unittest unittest.c)
target_link_libraries(md_unittest PRIVATE mold ${VIAMD_STDLIBS})

target_compile_definitions(md_unittest PRIVATE MD_UNITTEST_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/data\")

target_sources(md_unittest
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/test_allocator.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_array.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_bitop.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_bitfield.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_intrinsics.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_sync.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_os.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_pdb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_xtc.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_util.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_script.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_frame_cache.c
    ${CMAKE_CURRENT_SOURCE_DIR}/test_spatial_hash.c
    ${CMAKE_CURRENT_SOURCE_DIR}/rmsd.c
  )

# adapted from https://github.com/Kitware/CMake/blob/master/Modules/GoogleTest.cmake
# and https://github.com/Kitware/CMake/blob/master/Modules/GoogleTestAddTests.cmake
function(discover_tests TARGET)
  set(options)
  set(oneValueArgs
    TARGET
    WORKING_DIRECTORY
    DISCOVERY_TIMEOUT
  )
  set(multiValueArgs
    PROPERTIES
  )
  cmake_parse_arguments(
    ""
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
  )

  if(NOT _WORKING_DIRECTORY)
    set(_WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  endif()
  set(_TEST_LIST ${TARGET}_TESTS)
  if(NOT _DISCOVERY_TIMEOUT)
    set(_DISCOVERY_TIMEOUT 5)
  endif()

  set(_DISCOVERY_MODE POST_BUILD)

  get_property(
    has_counter
    TARGET ${TARGET}
    PROPERTY CTEST_DISCOVERED_TEST_COUNTER
    SET
  )
  if(has_counter)
    get_property(
      counter
      TARGET ${TARGET}
      PROPERTY CTEST_DISCOVERED_TEST_COUNTER
    )
    math(EXPR counter "${counter} + 1")
  else()
    set(counter 1)
  endif()
  set_property(
    TARGET ${TARGET}
    PROPERTY CTEST_DISCOVERED_TEST_COUNTER
    ${counter}
  )

  # Define rule to generate test list for aforementioned test executable
  set(ctest_file_base "${CMAKE_CURRENT_BINARY_DIR}/${TARGET}[${counter}]")
  set(ctest_include_file "${ctest_file_base}_include.cmake")
  set(ctest_tests_file "${ctest_file_base}_tests.cmake")

  add_custom_command(
    TARGET ${TARGET} POST_BUILD
    BYPRODUCTS "${ctest_tests_file}"
    COMMAND "${CMAKE_COMMAND}"
            -D "TEST_TARGET=${TARGET}"  # maybe not needed?
            -D "TEST_EXECUTABLE=$<TARGET_FILE:${TARGET}>"
            -D "TEST_WORKING_DIR=${_WORKING_DIRECTORY}"
            -D "TEST_PROPERTIES=${_PROPERTIES}"
            -D "CTEST_FILE=${ctest_tests_file}"
            -D "TEST_DISCOVERY_TIMEOUT=${_DISCOVERY_TIMEOUT}"
            -P "${CMAKE_CURRENT_LIST_DIR}/add_tests.cmake"
    VERBATIM
  )

  file(WRITE "${ctest_include_file}"
    "if(EXISTS \"${ctest_tests_file}\")\n"
    "  include(\"${ctest_tests_file}\")\n"
    "else()\n"
    "  add_test(${TARGET}_NOT_BUILT ${TARGET}_NOT_BUILT)\n"
    "endif()\n"
  )

  # Add discovered tests to directory TEST_INCLUDE_FILES
  set_property(DIRECTORY
    APPEND PROPERTY TEST_INCLUDE_FILES "${ctest_include_file}"
  )
endfunction()

discover_tests(md_unittest
  PROPERTIES
    LABELS "unit"
  )
