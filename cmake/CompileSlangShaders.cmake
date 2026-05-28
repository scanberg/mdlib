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
    set(SLANG_REFLECTION_TO_HEADER_SCRIPT "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/SlangReflectionToHeader.cmake")

    set(COMPILED_SHADERS "")
    set(REFLECTION_HEADERS "")

    if (NOT DEFINED SLANG_EXECUTABLE)
        message(FATAL_ERROR "SLANG_EXECUTABLE not defined but required for shader compilation")
    endif()

    # Metal reserves 'main', so Slang renames shader entry points to 'main_0'.
    # Runtime Metal loading explicitly requests that generated entry point name.
    set(SLANG_FLAGS "-Wno-40100")

    foreach(SRC ${SHADER_SOURCES})
        get_filename_component(NAME ${SRC} NAME_WE)
        get_filename_component(NAME_WITH_EXT ${SRC} NAME)
        get_filename_component(ABS_SRC ${SRC} ABSOLUTE)
        string(REPLACE "." "_" SYMBOL_NAME ${NAME})
        set(SYMBOL "${SHADER_NAMESPACE}_${SYMBOL_NAME}")
        set(REFLECTION_JSON "${GEN_DIR}/${NAME}.reflection.json")
        set(LOGICAL_REFLECTION_JSON "${GEN_DIR}/${NAME}.logical.reflection.json")
        set(REFLECTION_HEADER "${GEN_DIR}/${SYMBOL}.reflection.inl")

        if (NOT EXISTS ${ABS_SRC})
            message(FATAL_ERROR "Shader source file not found: ${ABS_SRC}")
        endif()

        if (MD_GPU_BACKEND STREQUAL "VULKAN")
            set(SPV_FILE "${GEN_DIR}/${NAME}.spv")
            message(STATUS "Compiling Slang shader ${ABS_SRC} -> ${SPV_FILE}")

            add_custom_command(
                OUTPUT ${SPV_FILE}
                BYPRODUCTS ${REFLECTION_JSON}
                COMMAND ${SLANG_EXECUTABLE}
                    ${ABS_SRC}
                    ${SLANG_FLAGS}
                    -target spirv
                    -emit-spirv-directly
                    -force-glsl-scalar-layout
                    -reflection-json ${REFLECTION_JSON}
                    -o ${SPV_FILE}
                DEPENDS ${ABS_SRC}
                COMMENT "slangc: ${NAME}.slang -> ${NAME}.spv"
            )

            set(BIN_FILE "${SPV_FILE}")
            set(LOGICAL_REFLECTION_JSON "${REFLECTION_JSON}")

        elseif (MD_GPU_BACKEND STREQUAL "METAL")
            set(MSL_FILE "${GEN_DIR}/${NAME}.metal")
            set(LIB_FILE "${GEN_DIR}/${NAME}.metallib")
            set(LOGICAL_SPV_FILE "${GEN_DIR}/${NAME}.logical.spv")
            message(STATUS "Compiling Slang shader ${ABS_SRC} -> ${LIB_FILE}")

            # Step 1: slangc -> MSL source
            add_custom_command(
                OUTPUT ${MSL_FILE}
                BYPRODUCTS ${REFLECTION_JSON}
                COMMAND ${SLANG_EXECUTABLE}
                    ${ABS_SRC}
                    ${SLANG_FLAGS}
                    -target metal
                    -reflection-json ${REFLECTION_JSON}
                    -o ${MSL_FILE}
                DEPENDS ${ABS_SRC}
                COMMENT "slangc: ${NAME}.slang -> ${NAME}.metal"
            )

            # Step 2: xcrun metal -> metallib, with explicit deployment target so
            # the metallib doesn't claim a higher OS version than CMAKE_OSX_DEPLOYMENT_TARGET.
            find_program(XCRUN_EXECUTABLE xcrun REQUIRED)
            # Slang lowers two-argument Interlocked* calls to Metal atomic_fetch_*
            # expressions whose returned old value is intentionally ignored.
            set(_METAL_FLAGS "-Wno-unused-variable")
            if (CMAKE_OSX_DEPLOYMENT_TARGET)
                list(APPEND _METAL_FLAGS "-mmacosx-version-min=${CMAKE_OSX_DEPLOYMENT_TARGET}")
            endif()

            add_custom_command(
                OUTPUT ${LIB_FILE}
                COMMAND ${XCRUN_EXECUTABLE} -sdk macosx metal ${_METAL_FLAGS} ${MSL_FILE} -o ${LIB_FILE}
                DEPENDS ${MSL_FILE}
                COMMENT "metal: ${NAME}.metal -> ${NAME}.metallib"
            )

            # Generate a logical SPIR-V reflection alongside the Metal-specific reflection.
            # The SPIR-V reflection preserves descriptor set space information, which Metal's
            # target-flattened reflection drops.
            add_custom_command(
                OUTPUT ${LOGICAL_REFLECTION_JSON}
                BYPRODUCTS ${LOGICAL_SPV_FILE}
                COMMAND ${SLANG_EXECUTABLE}
                    ${ABS_SRC}
                    ${SLANG_FLAGS}
                    -target spirv
                    -emit-spirv-directly
                    -force-glsl-scalar-layout
                    -reflection-json ${LOGICAL_REFLECTION_JSON}
                    -o ${LOGICAL_SPV_FILE}
                DEPENDS ${ABS_SRC}
                COMMENT "slangc: ${NAME}.slang -> ${NAME}.logical.reflection.json"
            )

            set(BIN_FILE "${LIB_FILE}")
        else()
            message(FATAL_ERROR "compile_slang_shaders: unsupported backend '${MD_GPU_BACKEND}' (Slang shaders require VULKAN or METAL)")
        endif()

        list(APPEND COMPILED_SHADERS ${BIN_FILE})
        add_custom_command(
            OUTPUT ${REFLECTION_HEADER}
            COMMAND ${CMAKE_COMMAND}
                -DINPUT=${REFLECTION_JSON}
                -DLOGICAL_INPUT=${LOGICAL_REFLECTION_JSON}
                -DOUTPUT=${REFLECTION_HEADER}
                -DSYMBOL=${SYMBOL}
                -DSOURCE=${NAME_WITH_EXT}
                -P ${SLANG_REFLECTION_TO_HEADER_SCRIPT}
            DEPENDS ${REFLECTION_JSON} ${LOGICAL_REFLECTION_JSON} ${SLANG_REFLECTION_TO_HEADER_SCRIPT}
            COMMENT "slang-reflect: ${NAME}.reflection.json -> ${SYMBOL}.reflection.inl"
            VERBATIM
        )
        list(APPEND REFLECTION_HEADERS ${REFLECTION_HEADER})
    endforeach()

    add_custom_target(${SHADER_NAMESPACE}_slang_shaders DEPENDS ${COMPILED_SHADERS})
    add_dependencies(${SHADER_TARGET} ${SHADER_NAMESPACE}_slang_shaders)

    get_filename_component(REFLECTION_OUT_NAME ${OUT_FILE} NAME_WE)
    set(REFLECTION_OUT_FILE "${GEN_DIR}/${REFLECTION_OUT_NAME}_reflection.inl")
    file(WRITE ${REFLECTION_OUT_FILE} "#pragma once\n\n")
    foreach(REFLECTION_HEADER ${REFLECTION_HEADERS})
        get_filename_component(REFLECTION_HEADER_NAME ${REFLECTION_HEADER} NAME)
        file(APPEND ${REFLECTION_OUT_FILE} "#include \"${REFLECTION_HEADER_NAME}\"\n")
    endforeach()

    add_custom_target(${SHADER_NAMESPACE}_slang_reflection DEPENDS ${REFLECTION_HEADERS})
    add_dependencies(${SHADER_TARGET} ${SHADER_NAMESPACE}_slang_reflection)

    embed_binary_files(
        TARGET ${SHADER_TARGET}
        NAMESPACE ${SHADER_NAMESPACE}
        OUTPUT ${OUT_FILE}
        FILES ${COMPILED_SHADERS}
    )
endfunction()
