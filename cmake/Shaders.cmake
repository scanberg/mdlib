# cmake/Shaders.cmake

function(md_compile_shaders OUT_FILE)
    set(options)
    set(oneValueArgs NAMESPACE)
    set(multiValueArgs SOURCES)
    cmake_parse_arguments(SHADER "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    set(GENERATED_C "${CMAKE_CURRENT_BINARY_DIR}/${OUT_FILE}")

    add_custom_command(
        OUTPUT ${GENERATED_C}
        COMMAND ${CMAKE_COMMAND}
                -DMD_GPU_BACKEND=${MD_GPU_BACKEND}
                -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
                -DSHADER_SOURCES=${SHADER_SOURCES}
                -DSHADER_NAMESPACE=${SHADER_NAMESPACE}
                -DOUTPUT_FILE=${GENERATED_C}
                -DGLSLC_EXECUTABLE=${GLSLC_EXECUTABLE}
                -DXCRUN_EXECUTABLE=${XCRUN_EXECUTABLE}
                -DSPIRV_CROSS_EXECUTABLE=${SPIRV_CROSS_EXECUTABLE}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/CompileShaders.cmake
        DEPENDS ${SHADER_SOURCES}
        COMMENT "Compiling shaders (${MD_GPU_BACKEND}) → ${OUT_FILE}"
        VERBATIM
    )

    set(${OUT_FILE} ${GENERATED_C} PARENT_SCOPE)
endfunction()

if (DEFINED MD_GPU_BACKEND)
    # Required for GLSL -> SPIR-V compilation
    find_program(GLSLC_EXECUTABLE glslc REQUIRED)
    message(STATUS "Found glslc: ${GLSLC_EXECUTABLE}")
endif()

if (MD_GPU_BACKEND STREQUAL "METAL")
    find_program(XCRUN_EXECUTABLE xcrun REQUIRED)
    find_program(SPIRV_CROSS_EXECUTABLE spirv-cross REQUIRED)
    message(STATUS "Found xcrun: ${XCRUN_EXECUTABLE}")
    message(STATUS "Found spirv-cross: ${SPIRV_CROSS_EXECUTABLE}")
endif()
