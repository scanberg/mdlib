# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION ${CMAKE_VERSION})

# Overwrite possibly existing ${_CTEST_FILE} with empty file
set(flush_tests_MODE WRITE)

# Flushes script to ${_CTEST_FILE}
macro(flush_script)
  file(${flush_tests_MODE} "${_CTEST_FILE}" "${script}")
  set(flush_tests_MODE APPEND)

  set(script "")
endmacro()

# Flushes tests_buffer to tests
macro(flush_tests_buffer)
  list(APPEND tests "${tests_buffer}")
  set(tests_buffer "")
endmacro()

macro(add_command NAME)
  set(_args "")
  foreach(_arg ${ARGN})
    if(_arg MATCHES "[^-./:a-zA-Z0-9_]")
      string(APPEND _args " [==[${_arg}]==]")
    else()
      string(APPEND _args " ${_arg}")
    endif()
  endforeach()
  string(APPEND script "${NAME}(${_args})\n")
  string(LENGTH "${script}" _script_len)
  if(${_script_len} GREATER "50000")
    flush_script()
  endif()
  # Unsets macro local variables to prevent leakage outside of this macro.
  unset(_args)
  unset(_script_len)
endmacro()

function(discover_tests_impl)

  cmake_parse_arguments(
    ""
    ""
    "TEST_EXECUTABLE;TEST_WORKING_DIR;CTEST_FILE;TEST_DISCOVERY_TIMEOUT"
    "TEST_PROPERTIES"
    ${ARGN}
  )
  set(properties ${_TEST_PROPERTIES})
  set(script)
  set(suite)
  set(tests)
  set(tests_buffer)

  # Run test executable to get list of available tests
  if(NOT EXISTS "${_TEST_EXECUTABLE}")
    message(FATAL_ERROR
      "Specified test executable does not exist.\n"
      "  Path: '${_TEST_EXECUTABLE}'"
    )
  endif()
  execute_process(
    COMMAND "${_TEST_EXECUTABLE}" --list-tests
    WORKING_DIRECTORY "${_TEST_WORKING_DIR}"
    TIMEOUT ${_TEST_DISCOVERY_TIMEOUT}
    OUTPUT_VARIABLE output
    RESULT_VARIABLE result
  )
  if(NOT ${result} EQUAL 0)
    string(REPLACE "\n" "\n    " output "${output}")
    message(FATAL_ERROR
      "Error running test executable.\n"
      "  Path: '${_TEST_EXECUTABLE}'\n"
      "  Result: ${result}\n"
      "  Output:\n"
      "    ${output}\n"
    )
  endif()

  # Preserve semicolon in test-parameters
  string(REPLACE [[;]] [[\;]] output "${output}")
  string(REPLACE "\n" ";" output "${output}")

  # Parse output
  foreach(line ${output})
    # extract suite and test names
    string(REPLACE "." ";" _line "${line}")
    list(GET _line 0 suite)
    list(GET _line 1 name)
    set(testname "${suite}.${name}")
    add_command(add_test
      "${testname}"
      "${_TEST_EXECUTABLE}"
      "--filter=${suite}.${name}"
    )
    add_command(set_tests_properties
      "${testname}"
      PROPERTIES
        WORKING_DIRECTORY "${_TEST_WORKING_DIR}"
        LABELS ${suite}
        SKIP_REGULAR_EXPRESSION "\\\\[  SKIPPED \\\\]"
        ${properties}
    )
    list(APPEND tests_buffer "${testname}")
    list(LENGTH tests_buffer tests_buffer_length)
    if(${tests_buffer_length} GREATER "250")
      flush_tests_buffer()
    endif()
  endforeach()

  # Create a list of all discovered tests, which users may use to e.g. set
  # properties on the tests
  flush_tests_buffer()
  add_command(set ${_TEST_LIST} ${tests})

  # Write CTest script
  flush_script()

endfunction()

if(CMAKE_SCRIPT_MODE_FILE)
  discover_tests_impl(
    TEST_EXECUTABLE ${TEST_EXECUTABLE}
    TEST_WORKING_DIR ${TEST_WORKING_DIR}
    CTEST_FILE ${CTEST_FILE}
    TEST_DISCOVERY_TIMEOUT ${TEST_DISCOVERY_TIMEOUT}
    TEST_PROPERTIES ${TEST_PROPERTIES}
  )
endif()
