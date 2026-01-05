include(${CMAKE_CURRENT_LIST_DIR}/EmbedBinaryFiles.cmake)

function(compile_shaders OUT_FILE)
    set(options)
    set(oneValueArgs NAMESPACE TARGET)
    set(multiValueArgs SOURCES)
    cmake_parse_arguments(SHADER "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if (NOT SHADER_TARGET OR NOT SHADER_NAMESPACE)
        message(FATAL_ERROR "compile_shaders: TARGET and NAMESPACE required")
    endif()

    set(GEN_DIR ${CMAKE_CURRENT_BINARY_DIR}/gen)
    file(MAKE_DIRECTORY ${GEN_DIR})

    set(COMPILED_SHADERS "")

    if (MD_GPU_BACKEND STREQUAL "VULKAN" OR MD_GPU_BACKEND STREQUAL "METAL")
        if (NOT DEFINED GLSLC_EXECUTABLE)
            message(FATAL_ERROR "GLSLC_EXECUTABLE not defined but required for shader compilation")
        endif()
    endif()

    foreach(SRC ${SHADER_SOURCES})
        get_filename_component(NAME ${SRC} NAME_WE)
        get_filename_component(ABS_SRC ${SRC} ABSOLUTE)

        if (MD_GPU_BACKEND STREQUAL "VULKAN" OR MD_GPU_BACKEND STREQUAL "METAL")
            set(SPV_FILE "${GEN_DIR}/${NAME}.spv")

            # Compile GLSL to SPIR-V
            add_custom_command(
                OUTPUT ${SPV_FILE}
                COMMAND ${GLSLC_EXECUTABLE} ${ABS_SRC} -o ${SPV_FILE}
                DEPENDS ${SRC}
            )

            set(BIN_FILE "${SPV_FILE}")

            if (MD_GPU_BACKEND STREQUAL "METAL")
                set(MSL_FILE "${GEN_DIR}/${NAME}.metal")
                set(LIB_FILE "${GEN_DIR}/${NAME}.metallib")

                # Cross-compile SPIR-V to MSL
                add_custom_command(
                    OUTPUT ${MSL_FILE}
                    COMMAND ${SPIRV_CROSS_EXECUTABLE} --msl --output ${MSL_FILE} ${SPV_FILE}
                    DEPENDS ${SPV_FILE}
                )

                # Compile MSL to metallib
                add_custom_command(
                    OUTPUT ${LIB_FILE}
                    COMMAND ${XCRUN_EXECUTABLE} -sdk macosx metal ${MSL_FILE} -o ${LIB_FILE}
                    DEPENDS ${MSL_FILE}
                )

                set(BIN_FILE "${LIB_FILE}")
            endif()
        else()
            # OpenGL backend, no compilation needed
            set(BIN_FILE "${ABS_SRC}")
        endif()

        list(APPEND COMPILED_SHADERS ${BIN_FILE})
    endforeach()

    add_custom_target(${SHADER_NAMESPACE}_shaders DEPENDS ${COMPILED_SHADERS})
    add_dependencies(${SHADER_TARGET} ${SHADER_NAMESPACE}_shaders)

    embed_binary_files(
        TARGET ${SHADER_TARGET}
        NAMESPACE ${SHADER_NAMESPACE}
        OUTPUT ${OUT_FILE}
        FILES ${COMPILED_SHADERS}
    )
endfunction()
