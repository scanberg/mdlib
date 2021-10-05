# add_compile_options_config(<CONFIG> <option> ...)
function(add_compile_options_config CONFIG)
    foreach(opt ${ARGN})
        add_compile_options("$<$<CONFIG:${CONFIG}>:${opt}>")
    endforeach()
endfunction()

file(GLOB SRC_FILES *)
file(GLOB CORE_FILES core/*)
file(GLOB SHADER_FILES shaders/*)
file(GLOB_RECURSE EXT_FILES ext/*)

source_group("src" FILES ${SRC_FILES})
source_group("core" FILES ${CORE_FILES})
source_group("ext" FILES ${EXT_FILES})
source_group("shaders" FILES ${SHADER_FILES})

add_library(mold ${CORE_FILES} ${SRC_FILES} ${EXT_FILES})

if (MSVC)
    message(STATUS "MSVC configuration was used")
    target_compile_options(mold PRIVATE /W4 /MP /GR- /fp:fast /arch:AVX)
    #add_compile_options(/W4 /MP /GR- /fp:fast /arch:AVX)# /arch:AVX2

    if (MD_PROFILE_COMPILE_TIMES)
        message(STATUS "MSVC compile time profiling enabled")
        #add_compile_options(/Bt+ /d1reportTime /d2cgsummary)
        add_compile_options(/d2cgsummary)
        #add_link_options(/time+)
    endif()
    add_compile_options_config(RELEASE /MT /Oi /GS-)
    add_compile_options_config(DEBUG /MTd)
else()
    message(STATUS "Assumed GCC or CLANG configuration was used")
    add_compile_options(
        -Wall -Wno-unused-function -Wno-unused-parameter -Wextra -Wpedantic -Wno-gnu-anonymous-struct -Wno-nested-anon-types -mavx #-maes
    )
    if (MD_PROFILE_COMPILE_TIMES)
        add_compile_options(-ftime-trace)
    endif()
endif()

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


target_link_libraries(mold PRIVATE ${MOLD_LIBS} ${VIAMD_STDLIBS})