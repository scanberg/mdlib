# CompileGpuShaders.cmake
#
# Compiles Slang kernels and vertex/fragment shaders for md_gpu and embeds the
# resulting binaries.
#
# md_gpu has no descriptor sets and no per-dispatch resource declarations, so
# there is no generated binding table. What the host side needs is the raw
# SPIR-V (Vulkan) or metallib/MSL (Metal) bytes, plus each kernel's group size
# and argument-struct size, which tools/check_gpu_arg_layout.py reads out of the
# compiled shader while it checks argument-struct portability.
#
#   compile_gpu_shaders(<out_header>
#       TARGET   <target>
#       NAMESPACE <prefix>
#       SOURCE   <file.slang>
#       ENTRIES  <entry> [<entry> ...]     # compute kernels
#       VERTEX   <entry> [<entry> ...]     # vertex stages
#       FRAGMENT <entry> [<entry> ...]     # fragment stages
#       DEPENDS  <file> [<file> ...]       # optional: files the source #includes
#   )
#
# At least one of ENTRIES / VERTEX / FRAGMENT is required. Produces, for each
# entry point, symbols
#
#   extern const uint8_t <prefix>_<stem>_<entry>_start[];
#   extern const size_t  <prefix>_<stem>_<entry>_byte_size;
#
# and, for a kernel or a raster stage respectively,
#
#   static inline md_gpu_kernel_desc_t <prefix>_<stem>_<entry>_kernel(void);
#   static inline md_gpu_shader_t      <prefix>_<stem>_<entry>_shader(void);
#
# all reachable by including <out_header>. The last two are what call sites use:
#
#   md_gpu_kernel_desc_t d = md_shader_topo_critical_points_main_kernel();
#   md_gpu_kernel_t k = md_gpu_kernel_create(device, &d);
#
#   md_gpu_pipeline_desc_t pd = {0};
#   pd.vertex   = md_shader_spheres_vs_main_shader();
#   pd.fragment = md_shader_spheres_fs_main_shader();

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

# Offline Metal compilation needs Apple's `metal` and `metallib`, which ship
# with Xcode's Metal toolchain and NOT with the macOS Command Line Tools. They
# are an optimisation, not a requirement: without them the MSL that slangc
# emits is embedded verbatim and compiled by the Metal framework at runtime
# (see md_mtl_library_from_blob in src/core/md_gpu_metal.m). Configuration must
# therefore never fail over a missing `metal`.
#
# `xcrun --find metal` is not a usable test: since Xcode 15 it resolves a stub
# that reports "Metal toolchain not installed" only when it is actually run.
# So the probe compiles a three-line kernel and believes the exit code.
#
# The result is cached, so -DMD_GPU_METAL_COMPILER=OFF on a machine that has
# Xcode forces the runtime path -- which is how configuration A of the build
# matrix gets tested without uninstalling anything.
function(md_gpu_probe_metal_compiler)
    # Unconditional: the metallib rules below need it even when the probe is
    # skipped because the user pinned MD_GPU_METAL_COMPILER on the command line.
    find_program(MD_GPU_XCRUN_EXECUTABLE xcrun)
    if (NOT MD_GPU_XCRUN_EXECUTABLE)
        set(MD_GPU_XCRUN_EXECUTABLE "xcrun" CACHE FILEPATH "xcrun" FORCE)
    endif()

    if (DEFINED CACHE{MD_GPU_METAL_COMPILER})
        return()
    endif()

    set(_ok OFF)
    if (APPLE)
        set(_probe_dir ${CMAKE_CURRENT_BINARY_DIR}/gen/metal_probe)
        file(MAKE_DIRECTORY ${_probe_dir})
        file(WRITE ${_probe_dir}/probe.metal
            "#include <metal_stdlib>\n"
            "using namespace metal;\n"
            "kernel void md_gpu_probe(device uint* p [[buffer(0)]]) { p[0] = 1u; }\n")
        execute_process(
            COMMAND ${MD_GPU_XCRUN_EXECUTABLE} -sdk macosx metal
                    -c ${_probe_dir}/probe.metal -o ${_probe_dir}/probe.air
            RESULT_VARIABLE _rc OUTPUT_QUIET ERROR_QUIET)
        if (_rc EQUAL 0)
            execute_process(
                COMMAND ${MD_GPU_XCRUN_EXECUTABLE} -sdk macosx metallib
                        ${_probe_dir}/probe.air -o ${_probe_dir}/probe.metallib
                RESULT_VARIABLE _rc OUTPUT_QUIET ERROR_QUIET)
            if (_rc EQUAL 0)
                set(_ok ON)
            endif()
        endif()
    endif()

    set(MD_GPU_METAL_COMPILER ${_ok} CACHE BOOL
        "Apple's offline Metal compiler (xcrun metal/metallib) is present and usable")
endfunction()

if (MD_GPU_BACKEND STREQUAL "METAL")
    md_gpu_probe_metal_compiler()
    # Printed on every configure, including when the cache value was pinned by
    # hand, so the build log always says which of the two paths is in effect.
    if (MD_GPU_METAL_COMPILER)
        message(STATUS "md_gpu: offline Metal compiler available -- kernels compiled to .metallib at build time")
    else()
        message(STATUS "md_gpu: no offline Metal compiler -- embedding MSL source, kernels compiled at runtime")
    endif()
endif()

function(compile_gpu_shaders OUT_HEADER)
    set(oneValueArgs TARGET NAMESPACE SOURCE)
    set(multiValueArgs ENTRIES VERTEX FRAGMENT DEPENDS)
    cmake_parse_arguments(G2 "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if (NOT G2_TARGET OR NOT G2_NAMESPACE OR NOT G2_SOURCE)
        message(FATAL_ERROR "compile_gpu_shaders: TARGET, NAMESPACE and SOURCE are required")
    endif()
    set(ALL_ENTRIES ${G2_ENTRIES} ${G2_VERTEX} ${G2_FRAGMENT})
    if (NOT ALL_ENTRIES)
        message(FATAL_ERROR "compile_gpu_shaders: at least one of ENTRIES, VERTEX, FRAGMENT is required")
    endif()
    set(STAGE_ARGS "")
    if (G2_VERTEX)
        list(APPEND STAGE_ARGS --vertex ${G2_VERTEX})
    endif()
    if (G2_FRAGMENT)
        list(APPEND STAGE_ARGS --fragment ${G2_FRAGMENT})
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
    if (NOT Python3_Interpreter_FOUND)
        message(FATAL_ERROR "compile_gpu_shaders: Python 3 is required to check argument "
                            "structs and to generate kernel descriptors")
    endif()
    set(LINT_SCRIPT ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../tools/check_gpu_arg_layout.py)
    set(LINT_STAMP  ${GEN_DIR}/${STEM}.arglayout.stamp)
    # Kernel descriptors: one md_gpu_kernel_desc_t per entry point, with the
    # group size and argument-struct size read from the compiled shader.
    set(KERNELS_INL ${GEN_DIR}/${STEM}_kernels.inl)
    add_custom_command(
        OUTPUT ${LINT_STAMP} ${KERNELS_INL}
        COMMAND ${Python3_EXECUTABLE} ${LINT_SCRIPT}
            --slangc ${SLANG_EXECUTABLE}
            --bindless-space ${MD_GPU_BINDLESS_SPACE}
            --emit ${KERNELS_INL}
            --namespace ${G2_NAMESPACE}
            ${ABS_SRC} ${G2_ENTRIES} ${STAGE_ARGS}
        COMMAND ${CMAKE_COMMAND} -E touch ${LINT_STAMP}
        DEPENDS ${ABS_SRC} ${MD_GPU_SHADER_DEPS} ${G2_DEPENDS} ${LINT_SCRIPT}
        COMMENT "md_gpu: checking ${STEM}.slang argument-struct portability"
        VERBATIM
    )
    target_sources(${G2_TARGET} PRIVATE ${KERNELS_INL})

    set(BIN_FILES "")
    foreach(ENTRY ${ALL_ENTRIES})
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
                DEPENDS ${ABS_SRC} ${MD_GPU_SHADER_DEPS} ${G2_DEPENDS} ${LINT_STAMP}
                COMMENT "slangc: ${STEM}.slang [${ENTRY}] -> ${STEM}_${ENTRY}.spv"
            )
        elseif (MD_GPU_BACKEND STREQUAL "METAL")
            set(MSL "${GEN_DIR}/${STEM}_${ENTRY}.metal")
            set(AIR "${GEN_DIR}/${STEM}_${ENTRY}.air")
            add_custom_command(
                OUTPUT ${MSL}
                COMMAND ${SLANG_EXECUTABLE}
                    ${ABS_SRC} ${SLANG_FLAGS}
                    -target metal
                    -entry ${ENTRY}
                    -o ${MSL}
                DEPENDS ${ABS_SRC} ${MD_GPU_SHADER_DEPS} ${G2_DEPENDS} ${LINT_STAMP}
                COMMENT "slangc: ${STEM}.slang [${ENTRY}] -> ${STEM}_${ENTRY}.metal"
            )
            if (MD_GPU_METAL_COMPILER)
                set(BIN "${GEN_DIR}/${STEM}_${ENTRY}.metallib")
                add_custom_command(
                    OUTPUT ${BIN}
                    COMMAND ${MD_GPU_XCRUN_EXECUTABLE} -sdk macosx metal -c ${MSL} -o ${AIR}
                    COMMAND ${MD_GPU_XCRUN_EXECUTABLE} -sdk macosx metallib ${AIR} -o ${BIN}
                    DEPENDS ${MSL}
                    COMMENT "metallib: ${STEM}_${ENTRY}.metal -> ${STEM}_${ENTRY}.metallib"
                )
            else()
                # Embed the MSL itself. embed_binary_files() names its symbols
                # from NAME_WE, so <stem>_<entry>.metal and <stem>_<entry>.metallib
                # produce the *same* md_shader_<stem>_<entry>_start symbol and no
                # consumer changes. The backend tells the two apart by their bytes.
                set(BIN ${MSL})
            endif()
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
    # The embed header declares the code symbols; the descriptors built on
    # them come after, so including <OUT_HEADER> gives both.
    file(APPEND ${GEN_DIR}/${OUT_HEADER} "#include \"${STEM}_kernels.inl\"\n")
endfunction()
