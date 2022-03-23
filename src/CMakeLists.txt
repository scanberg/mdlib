# add_compile_options_config(<CONFIG> <option> ...)
function(add_compile_options_config CONFIG)
    foreach(opt ${ARGN})
        add_compile_options("$<$<CONFIG:${CONFIG}>:${opt}>")
    endforeach()
endfunction()

file(GLOB SRC_FILES *)
file(GLOB CORE_FILES core/*)
file(GLOB SHADER_FILES shaders/*)

source_group("src" FILES ${SRC_FILES})
source_group("core" FILES ${CORE_FILES})
source_group("shaders" FILES ${SHADER_FILES})

add_library(mold ${CORE_FILES} ${SRC_FILES} ${EXT_FILES})
set(MOLD_WARNING_FLAGS "")

if (CMAKE_C_COMPILER_ID MATCHES "Clang")
  set(MOLD_WARNING_FLAGS -Wall -Wno-missing-field-initializers -Wno-missing-braces -Wno-unused-function -Wno-unused-parameter -Wextra -Wpedantic -Wno-gnu-anonymous-struct -Wno-nested-anon-types)
  target_compile_options(mold PRIVATE -march=native)
  if (MD_PROFILE_COMPILE_TIMES)
    add_compile_options(-ftime-trace)
  endif()
  if (CMAKE_C_COMPILER_ID STREQUAL "Clang")
    if (CMAKE_C_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC") # clang-cl
    message(STATUS "MOLD: MSVC-Clang configuration was used")
    elseif (CMAKE_C_COMPILER_FRONTEND_VARIANT STREQUAL "GNU") # clang native
    message(STATUS "MOLD: Clang configuration was used")
    endif()
  endif()
elseif (CMAKE_C_COMPILER_ID STREQUAL "GNU") # GCC
    message(STATUS "MOLD: GCC configuration was used")
    set(MOLD_WARNING_FLAGS -Wall -Wno-missing-field-initializers -Wno-missing-braces -Wno-unused-function -Wno-unused-parameter -Wextra -Wpedantic)
    target_compile_options(mold PRIVATE -march=native)
elseif ((CMAKE_C_COMPILER_ID STREQUAL "MSVC")) # msvc
    message(STATUS "MOLD: MSVC configuration was used")
    set(MOLD_WARNING_FLAGS "/W4")
    target_compile_options(mold PRIVATE /MP /GR- /fp:fast /arch:AVX)
    
    if (MD_PROFILE_COMPILE_TIMES)
        message(STATUS "MOLD: MSVC compile time profiling enabled")
        add_compile_options(/d2cgsummary)
    endif()
    add_compile_options_config(RELEASE /MT /Oi /GS-)
    add_compile_options_config(DEBUG /MTd)
endif()

add_subdirectory(ext)

target_include_directories(mold
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

if (MD_USE_RELATIVE_PATHS)
    target_compile_definitions(mold PRIVATE MD_SHADER_DIR=\"shaders\")

	# COPY SHADER FILES TO EXECUTABLE DIRECTORY TO KEEP RELATIVE PATH
    # CANNOT COPY ENTIRE DIRECTORY -> WILL REMOVE ANY EXISTING DIRECTORY
    file(GLOB_RECURSE COPY_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} shaders/*)

    foreach(FILENAME ${COPY_FILES})
    set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}")
    set(DST "${CMAKE_BINARY_DIR}/$<CONFIG>/${FILENAME}")

    add_custom_command(
        TARGET mold POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy ${SRC} ${DST}
        )
    endforeach(FILENAME)
else()
	target_compile_definitions(mold PRIVATE MD_SHADER_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/shaders\")
endif()

if (MSVC)
    target_compile_definitions(mold PRIVATE _CRT_SECURE_NO_WARNINGS)
endif()

if (CMAKE_COMPILER_IS_GNUCC)
find_package(Threads REQUIRED)
    target_compile_definitions(mold PRIVATE _GNU_SOURCE)
    list(APPEND MOLD_LIBS m Threads::Threads)
endif()

target_compile_options(mold PRIVATE ${MOLD_WARNING_FLAGS})

target_link_libraries(mold PRIVATE mold_ext ${MOLD_LIBS} ${VIAMD_STDLIBS})