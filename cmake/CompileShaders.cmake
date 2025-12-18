cmake_minimum_required(VERSION 3.20)

if (NOT DEFINED SHADER_SOURCES)
    message(FATAL_ERROR "CompileShaders.cmake: SHADER_SOURCES not set")
endif()

file(WRITE "${OUTPUT_FILE}" "// Generated shader resources\n\n")

foreach(SHADER ${SHADER_SOURCES})
    get_filename_component(SHADER_NAME ${SHADER} NAME)
    get_filename_component(SHADER_BASE ${SHADER} NAME_WE)

    message(STATUS "Processing shader: ${SHADER_NAME}")

    if (MD_GPU_BACKEND STREQUAL "VULKAN")

        set(SPV_FILE "${SHADER_BASE}.spv")

        execute_process(
            COMMAND ${GLSLC_EXECUTABLE} ${SHADER} -o ${SPV_FILE}
            RESULT_VARIABLE RC
        )
        if (RC)
            message(FATAL_ERROR "glslc failed for ${SHADER}")
        endif()

        # Read binary → hex
        file(READ "${SPV_FILE}" BIN HEX)

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

    elseif (MD_GPU_BACKEND STREQUAL "METAL")

        set(SPV_FILE "${SHADER_BASE}.spv")
        set(MSL_FILE "${SHADER_BASE}.metal")
        set(AIR_FILE "${SHADER_BASE}.air")
        set(LIB_FILE "${SHADER_BASE}.metallib")

        execute_process(
            COMMAND ${GLSLC_EXECUTABLE} ${SHADER} -o ${SPV_FILE}
            RESULT_VARIABLE RC
        )
        if (RC)
            message(FATAL_ERROR "glslc failed for ${SHADER}")
        endif()

        execute_process(
            COMMAND ${SPIRV_CROSS_EXECUTABLE} ${SPV_FILE} --msl --output ${MSL_FILE}
            RESULT_VARIABLE RC
        )
        if (RC)
            message(FATAL_ERROR "spirv-cross failed for ${SHADER}")
        endif()

        execute_process(
            COMMAND ${XCRUN_EXECUTABLE} -sdk macosx metal
                    ${MSL_FILE} -o ${AIR_FILE}
            RESULT_VARIABLE RC
        )
        if (RC)
            message(FATAL_ERROR "metal failed for ${SHADER}")
        endif()

        execute_process(
            COMMAND ${XCRUN_EXECUTABLE} -sdk macosx metallib
                    ${AIR_FILE} -o ${LIB_FILE}
            RESULT_VARIABLE RC
        )
        if (RC)
            message(FATAL_ERROR "metallib failed for ${SHADER}")
        endif()

        file(READ "${LIB_FILE}" BIN HEX)

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

    else()
        message(FATAL_ERROR "Unknown MD_GPU_BACKEND: ${MD_GPU_BACKEND}")
    endif()
endforeach()
