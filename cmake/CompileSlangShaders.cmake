include(${CMAKE_CURRENT_LIST_DIR}/EmbedBinaryFiles.cmake)

function(compile_slang_shaders OUT_FILE)
    set(options)
    set(oneValueArgs NAMESPACE TARGET)
    set(multiValueArgs SOURCES)
    cmake_parse_arguments(SHADER "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if (NOT SHADER_TARGET OR NOT SHADER_NAMESPACE)
        message(FATAL_ERROR "compile_slang_shaders: TARGET and NAMESPACE required")
    endif()

    set(GEN_DIR ${CMAKE_CURRENT_BINARY_DIR}/gen)
    file(MAKE_DIRECTORY ${GEN_DIR})

    set(COMPILED_SHADERS "")

    if (NOT DEFINED SLANGC_EXECUTABLE)
        message(FATAL_ERROR "SLANGC_EXECUTABLE not defined but required for shader compilation")
    endif()

    foreach(SRC ${SHADER_SOURCES})
        get_filename_component(NAME ${SRC} NAME_WE)
        get_filename_component(ABS_SRC ${SRC} ABSOLUTE)

        if (NOT EXISTS ${ABS_SRC})
            message(FATAL_ERROR "Shader source file not found: ${ABS_SRC}")
        endif()

        if (MD_GPU_BACKEND STREQUAL "VULKAN")
            set(SPV_FILE "${GEN_DIR}/${NAME}.spv")
            message(STATUS "Compiling Slang shader ${ABS_SRC} -> ${SPV_FILE}")

            add_custom_command(
                OUTPUT ${SPV_FILE}
                COMMAND ${SLANGC_EXECUTABLE}
                    ${ABS_SRC}
                    -target spirv
                    -emit-spirv-directly
                    -force-glsl-scalar-layout
                    -o ${SPV_FILE}
                DEPENDS ${ABS_SRC}
                COMMENT "slangc: ${NAME}.slang -> ${NAME}.spv"
            )

            set(BIN_FILE "${SPV_FILE}")

        elseif (MD_GPU_BACKEND STREQUAL "METAL")
            set(LIB_FILE "${GEN_DIR}/${NAME}.metallib")
            message(STATUS "Compiling Slang shader ${ABS_SRC} -> ${LIB_FILE}")

            add_custom_command(
                OUTPUT ${LIB_FILE}
                COMMAND ${SLANGC_EXECUTABLE}
                    ${ABS_SRC}
                    -target metallib
                    -o ${LIB_FILE}
                DEPENDS ${ABS_SRC}
                COMMENT "slangc: ${NAME}.slang -> ${NAME}.metallib"
            )

            set(BIN_FILE "${LIB_FILE}")
        else()
            message(FATAL_ERROR "compile_slang_shaders: unsupported backend '${MD_GPU_BACKEND}' (Slang shaders require VULKAN or METAL)")
        endif()

        list(APPEND COMPILED_SHADERS ${BIN_FILE})
    endforeach()

    add_custom_target(${SHADER_NAMESPACE}_slang_shaders DEPENDS ${COMPILED_SHADERS})
    add_dependencies(${SHADER_TARGET} ${SHADER_NAMESPACE}_slang_shaders)

    embed_binary_files(
        TARGET ${SHADER_TARGET}
        NAMESPACE ${SHADER_NAMESPACE}
        OUTPUT ${OUT_FILE}
        FILES ${COMPILED_SHADERS}
    )
endfunction()
