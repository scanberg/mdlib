# CompileGpuShaders.cmake
#
# Compiles Slang compute kernels for md_gpu and embeds the resulting binaries.
#
# Deliberately much simpler than CompileSlangShaders.cmake: md_gpu has no
# descriptor sets and no per-dispatch resource declarations, so there is no
# reflection step and no generated binding table. All the host side needs is
# the raw SPIR-V (Vulkan) or metallib (Metal) bytes.
#
#   compile_gpu_shaders(<out_header>
#       TARGET   <target>
#       NAMESPACE <prefix>
#       SOURCE   <file.slang>
#       ENTRIES  <entry> [<entry> ...]
#   )
#
# Produces, for each entry point, symbols
#
#   extern const uint8_t <prefix>_<stem>_<entry>_start[];
#   extern const size_t  <prefix>_<stem>_<entry>_byte_size;

include_guard(GLOBAL)
include(${CMAKE_CURRENT_LIST_DIR}/EmbedBinaryFiles.cmake)

# Descriptor set that Slang's DescriptorHandle<T> heap is placed in. Slang then
# uses (space, binding 0) for samplers, 2 for sampled images and 3 for storage
# images. Must match MD_VK_BINDLESS_SPACE in src/core/md_gpu_vulkan.c.
set(MD_GPU_BINDLESS_SPACE 0 CACHE STRING "Descriptor space for the md_gpu bindless heap")

# Every md_gpu kernel #includes the prelude, and the prelude decides the
# bindless binding assignment -- editing it changes the SPIR-V of every shader.
# slangc does not report its includes to the build system, so name the prelude
# explicitly or a stale binary survives a prelude change and the descriptors
# silently stop matching the set layout.
set(MD_GPU_SHADER_DEPS "${CMAKE_CURRENT_LIST_DIR}/../src/shaders/md_gpu.slang")

function(compile_gpu_shaders OUT_HEADER)
    set(oneValueArgs TARGET NAMESPACE SOURCE)
    set(multiValueArgs ENTRIES)
    cmake_parse_arguments(G2 "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if (NOT G2_TARGET OR NOT G2_NAMESPACE OR NOT G2_SOURCE)
        message(FATAL_ERROR "compile_gpu_shaders: TARGET, NAMESPACE and SOURCE are required")
    endif()
    if (NOT DEFINED SLANG_EXECUTABLE)
        message(FATAL_ERROR "compile_gpu_shaders: SLANG_EXECUTABLE not defined")
    endif()

    get_filename_component(ABS_SRC ${G2_SOURCE} ABSOLUTE)
    get_filename_component(STEM ${G2_SOURCE} NAME_WE)
    if (NOT EXISTS ${ABS_SRC})
        message(FATAL_ERROR "compile_gpu_shaders: source not found: ${ABS_SRC}")
    endif()

    set(GEN_DIR ${CMAKE_CURRENT_BINARY_DIR}/gen)
    file(MAKE_DIRECTORY ${GEN_DIR})

    # Metal reserves 'main', so Slang renames entry points; silence that note.
    set(SLANG_FLAGS "-Wno-40100")

    # Reject argument structs whose layout differs between SPIR-V and MSL.
    # Vectors and bindless handles are the constructs that diverge, and they
    # diverge silently, so this is checked at build time rather than trusted.
    #
    # The stamp is wired in as a file dependency of the shader binaries rather
    # than as a custom target. A target per kernel would put ten pseudo-targets
    # in every IDE's target list to run one script; a stamp the compile step
    # already depends on gives the same ordering and the same incremental
    # behaviour with nothing to look at. Same pattern EmbedBinaryFiles.cmake
    # uses for its generated sources.
    find_package(Python3 COMPONENTS Interpreter QUIET)
    set(LINT_STAMP "")
    if (Python3_Interpreter_FOUND)
        set(LINT_SCRIPT ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../tools/check_gpu_arg_layout.py)
        set(LINT_STAMP  ${GEN_DIR}/${STEM}.arglayout.stamp)
        add_custom_command(
            OUTPUT ${LINT_STAMP}
            COMMAND ${Python3_EXECUTABLE} ${LINT_SCRIPT}
                --slangc ${SLANG_EXECUTABLE}
                --bindless-space ${MD_GPU_BINDLESS_SPACE}
                ${ABS_SRC} ${G2_ENTRIES}
            COMMAND ${CMAKE_COMMAND} -E touch ${LINT_STAMP}
            DEPENDS ${ABS_SRC} ${MD_GPU_SHADER_DEPS} ${LINT_SCRIPT}
            COMMENT "md_gpu: checking ${STEM}.slang argument-struct portability"
            VERBATIM
        )
    else()
        message(WARNING "md_gpu: Python3 not found, skipping argument-struct portability check")
    endif()

    set(BIN_FILES "")
    foreach(ENTRY ${G2_ENTRIES})
        if (MD_GPU_BACKEND STREQUAL "VULKAN")
            set(BIN "${GEN_DIR}/${STEM}_${ENTRY}.spv")
            add_custom_command(
                OUTPUT ${BIN}
                COMMAND ${SLANG_EXECUTABLE}
                    ${ABS_SRC} ${SLANG_FLAGS}
                    -target spirv
                    -emit-spirv-directly
                    -profile glsl_450
                    -bindless-space-index ${MD_GPU_BINDLESS_SPACE}
                    -entry ${ENTRY}
                    -o ${BIN}
                DEPENDS ${ABS_SRC} ${MD_GPU_SHADER_DEPS} ${LINT_STAMP}
                COMMENT "slangc: ${STEM}.slang [${ENTRY}] -> ${STEM}_${ENTRY}.spv"
            )
        elseif (MD_GPU_BACKEND STREQUAL "METAL")
            set(MSL "${GEN_DIR}/${STEM}_${ENTRY}.metal")
            set(AIR "${GEN_DIR}/${STEM}_${ENTRY}.air")
            set(BIN "${GEN_DIR}/${STEM}_${ENTRY}.metallib")
            add_custom_command(
                OUTPUT ${MSL}
                COMMAND ${SLANG_EXECUTABLE}
                    ${ABS_SRC} ${SLANG_FLAGS}
                    -target metal
                    -entry ${ENTRY}
                    -o ${MSL}
                DEPENDS ${ABS_SRC} ${MD_GPU_SHADER_DEPS} ${LINT_STAMP}
                COMMENT "slangc: ${STEM}.slang [${ENTRY}] -> ${STEM}_${ENTRY}.metal"
            )
            add_custom_command(
                OUTPUT ${BIN}
                COMMAND xcrun -sdk macosx metal -c ${MSL} -o ${AIR}
                COMMAND xcrun -sdk macosx metallib ${AIR} -o ${BIN}
                DEPENDS ${MSL}
                COMMENT "metallib: ${STEM}_${ENTRY}.metal -> ${STEM}_${ENTRY}.metallib"
            )
        else()
            message(FATAL_ERROR "compile_gpu_shaders: unknown MD_GPU_BACKEND '${MD_GPU_BACKEND}'")
        endif()
        list(APPEND BIN_FILES ${BIN})
    endforeach()

    embed_binary_files(
        TARGET     ${G2_TARGET}
        NAMESPACE  ${G2_NAMESPACE}
        OUTPUT     ${OUT_HEADER}
        FILES      ${BIN_FILES}
    )
endfunction()
