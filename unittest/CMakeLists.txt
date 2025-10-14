cmake_minimum_required(VERSION 3.15)

project(md_unittest VERSION 0.1 LANGUAGES C)

add_executable(md_unittest unittest.c)

target_compile_definitions(md_unittest PRIVATE MD_UNITTEST_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/../test_data\" UTEST_USE_OLD_QPC ${MD_DEFINES})
target_compile_options(md_unittest PRIVATE ${MD_FLAGS_EXT} ${MD_FLAGS_AVX} $<$<CONFIG:Debug>:${MD_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_FLAGS_REL}>)
target_link_options(md_unittest PRIVATE ${MD_LINK_FLAGS} $<$<CONFIG:Debug>:${MD_LINK_FLAGS_DEB}> $<$<CONFIG:Release>:${MD_LINK_FLAGS_REL}>)
target_link_libraries(md_unittest PRIVATE mdlib ${MD_LIBS})

set (SRC_FILES
    test_allocator.c
    test_api_equivalence.c
    test_array.c
    test_atomic.c
    test_bitop.c
    test_bitfield.c
    test_element_guess_compat.c
    test_str.c
    test_hash.c
    test_edr.c
    test_cube.c
    test_intrinsics.c
    test_simd.c
    test_vec_math.c
    test_os.c
    test_pdb.c
    test_gro.c
    test_trr.c
    test_mmcif.c
    test_util.c
    test_unit.c
    test_script.c
    test_filter.c
    test_frame_cache.c
    test_spatial_hash.c
    test_string_builder.c
    test_base64.c
    test_xtc.c
    test_xyz.c
    test_xvg.c
    rmsd.c
    test_lammps.c
)

if (${MD_ENABLE_VLX})
    set (SRC_FILES ${SRC_FILES} test_vlx.c test_gto.c)
endif()

target_sources(md_unittest PRIVATE ${SRC_FILES})

target_include_directories(md_unittest
  PRIVATE
  ../ext/xtc/
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
