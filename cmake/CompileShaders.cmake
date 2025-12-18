cmake_minimum_required(VERSION 3.20)

if (NOT DEFINED SHADER_SOURCES)
    message(FATAL_ERROR "CompileShaders.cmake: SHADER_SOURCES not set")
endif()

file(WRITE "${OUTPUT_FILE}" "// Generated shader resources\n\n")

foreach(SHADER ${SHADER_SOURCES})
    if(NOT IS_ABSOLUTE "${SHADER}")
        set(SHADER_ABS "${SOURCE_DIR}/${SHADER}")
    else()
        set(SHADER_ABS "${SHADER}")
    endif()
    get_filename_component(SHADER_NAME ${SHADER} NAME)
    get_filename_component(SHADER_BASE ${SHADER} NAME_WE)

    message(STATUS "Processing shader: ${SHADER_NAME}")

    if (MD_GPU_BACKEND STREQUAL "VULKAN" OR MD_GPU_BACKEND STREQUAL "METAL")
        set(SPV_FILE "${SHADER_BASE}.spv")

        execute_process(
            COMMAND ${GLSLC_EXECUTABLE} ${SHADER_ABS} -o ${SPV_FILE}
            RESULT_VARIABLE RC
            ERROR_VARIABLE  GLSLC_STDERR
        )
        if (RC)
            message(STATUS "COMMAND: ${GLSLC_EXECUTABLE} ${SHADER_ABS} -o ${SPV_FILE}")
            message(FATAL_ERROR "glslc failed for ${SHADER_NAME}")
            message(STATUS "STDERR: ${GLSLC_STDERR}")
        endif()
    endif()

    # Set empty BIN_FILE variable (the compiled binary file to write to output)
    set(BIN_FILE "")

    if (MD_GPU_BACKEND STREQUAL "VULKAN")
        # Assert that SPV_FILE exists
        if (NOT EXISTS "${SPV_FILE}")
            message(FATAL_ERROR "SPIR-V file not found: ${SPV_FILE}")
        endif()

        set(BIN_FILE "${SPV_FILE}")

    elseif (MD_GPU_BACKEND STREQUAL "METAL")
        # Assert that SPV_FILE exists
        if (NOT EXISTS "${SPV_FILE}")
            message(FATAL_ERROR "SPIR-V file not found: ${SPV_FILE}")
        endif()

        set(MSL_FILE "${SHADER_BASE}.metal")
        set(LIB_FILE "${SHADER_BASE}.metallib")

        execute_process(
            COMMAND ${SPIRV_CROSS_EXECUTABLE} ${SPV_FILE} --msl --output ${MSL_FILE}
            RESULT_VARIABLE RC
            ERROR_VARIABLE  SPIRV_CROSS_STDERR
        )
        if (RC)
            message(STATUS "COMMAND: ${SPIRV_CROSS_EXECUTABLE} ${SPV_FILE} --msl --output ${MSL_FILE}")
            message(STATUS "STDERR: ${SPIRV_CROSS_STDERR}")
            message(FATAL_ERROR "spirv-cross failed for ${SHADER_NAME}")
        endif()

        execute_process(
            COMMAND ${XCRUN_EXECUTABLE} -sdk macosx metal
                    ${MSL_FILE} -o ${LIB_FILE}
            RESULT_VARIABLE RC
            OUTPUT_VARIABLE METAL_STDOUT
            ERROR_VARIABLE  METAL_STDERR
        )
        if (RC)
            message(STATUS "COMMAND: ${XCRUN_EXECUTABLE} -sdk macosx metal ${MSL_FILE} -o ${LIB_FILE}")
            message(STATUS "STDOUT: ${METAL_STDOUT}")
            message(STATUS "STDERR: ${METAL_STDERR}")
            message(FATAL_ERROR "metal failed for ${SHADER_NAME}")
        endif()

        set(BIN_FILE "${LIB_FILE}")
    else()
        message(FATAL_ERROR "Unknown MD_GPU_BACKEND: ${MD_GPU_BACKEND}")
    endif()

    # Read binary → hex
    file(READ "${BIN_FILE}" BIN HEX)

    # Emit C array
    file(APPEND "${OUTPUT_FILE}"
        "// ${SHADER_NAME}\n"
        "static const unsigned char ${SHADER_NAMESPACE}_${SHADER_BASE}[] = {"
    )

    string(LENGTH "${BIN}" BIN_LEN)
    math(EXPR BYTE_COUNT "${BIN_LEN} / 2")

    foreach(I RANGE 0 ${BYTE_COUNT}-1)
        math(EXPR OFFSET "${I} * 2")
        string(SUBSTRING "${BIN}" ${OFFSET} 2 BYTE)
        file(APPEND "${OUTPUT_FILE}" "0x${BYTE},")
    endforeach()

    file(APPEND "${OUTPUT_FILE}" "};\n\n")
endforeach()
