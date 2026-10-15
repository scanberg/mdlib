# CompileGpuShaders.cmake
#
# Compiles Slang kernels and vertex/fragment shaders for md_gpu and embeds the
# resulting binaries.
#
# md_gpu has no descriptor sets and no per-dispatch resource declarations, so
# there is no generated binding table. What the host side needs is the raw
# SPIR-V (Vulkan) or metallib/MSL (Metal) bytes, plus each kernel's group size
# and argument-struct size, which are read out of the compiled shader while its
# argument structs are checked for portability (tools/check_gpu_arg_layout.py).
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
#
# Each entry point is built by a single command, tools/md_gpu_build_shader.py:
# one slangc run for both SPIR-V and MSL, the argument-struct check on the two,
# the SPIR-V post-pass, the descriptor, and the binary embedded as C. Nothing
# else in the target waits on that command, so every generator can run all of
# them in parallel -- MSBuild included, since CMake 3.27 (policy CMP0147).

include_guard(GLOBAL)

# Every shader build step is a Python script. Found once, here: inside
# compile_gpu_shaders() the result would not outlive the call, and FindPython
# re-validates the interpreter every time it runs.
find_package(Python3 COMPONENTS Interpreter QUIET)

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
    if (NOT DEFINED SLANG_EXECUTABLE)
        message(FATAL_ERROR "compile_gpu_shaders: SLANG_EXECUTABLE not defined")
    endif()
    if (NOT Python3_Interpreter_FOUND)
        message(FATAL_ERROR "compile_gpu_shaders: Python 3 is required to build md_gpu shaders")
    endif()
    if (NOT MD_GPU_BACKEND STREQUAL "VULKAN" AND NOT MD_GPU_BACKEND STREQUAL "METAL")
        message(FATAL_ERROR "compile_gpu_shaders: unknown MD_GPU_BACKEND '${MD_GPU_BACKEND}'")
    endif()

    get_filename_component(ABS_SRC ${G2_SOURCE} ABSOLUTE)
    get_filename_component(STEM ${G2_SOURCE} NAME_WE)
    if (NOT EXISTS ${ABS_SRC})
        message(FATAL_ERROR "compile_gpu_shaders: source not found: ${ABS_SRC}")
    endif()

    set(GEN_DIR ${CMAKE_CURRENT_BINARY_DIR}/gen)
    file(MAKE_DIRECTORY ${GEN_DIR})

    set(TOOLS_DIR    ${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../tools)
    set(BUILD_SCRIPT ${TOOLS_DIR}/md_gpu_build_shader.py)
    set(BIN2C_SCRIPT ${TOOLS_DIR}/bin2c.py)
    # The build script imports the others, so editing any of them rebuilds.
    set(SCRIPT_DEPS ${BUILD_SCRIPT} ${BIN2C_SCRIPT}
                    ${TOOLS_DIR}/check_gpu_arg_layout.py ${TOOLS_DIR}/spirv_deforward.py)

    # With Apple's offline compiler the embedded binary is a metallib, built
    # from the MSL by a second command. Otherwise the build script embeds the
    # binary itself: the SPIR-V, or on Metal the MSL, which the Metal framework
    # then compiles at runtime (see md_mtl_library_from_blob).
    set(USE_METALLIB OFF)
    if (MD_GPU_BACKEND STREQUAL "METAL" AND MD_GPU_METAL_COMPILER)
        set(USE_METALLIB ON)
    endif()

    set(DECLS "")
    set(DESCS "")
    foreach(ENTRY ${ALL_ENTRIES})
        if (ENTRY IN_LIST G2_ENTRIES)
            set(STAGE compute)
        elseif (ENTRY IN_LIST G2_VERTEX)
            set(STAGE vertex)
        else()
            set(STAGE fragment)
        endif()

        # Must match the symbol md_gpu_build_shader.py and bin2c.py define.
        set(SYMBOL ${G2_NAMESPACE}_${STEM}_${ENTRY})
        set(SPV    ${GEN_DIR}/${STEM}_${ENTRY}.spv)
        set(MSL    ${GEN_DIR}/${STEM}_${ENTRY}.metal)
        set(DESC   ${GEN_DIR}/${STEM}_${ENTRY}_desc.inl)
        set(C_FILE ${GEN_DIR}/${SYMBOL}.c)

        set(BUILD_CMD ${Python3_EXECUTABLE} ${BUILD_SCRIPT}
            --slangc         ${SLANG_EXECUTABLE}
            --bindless-space ${MD_GPU_BINDLESS_SPACE}
            --backend        ${MD_GPU_BACKEND}
            --stage          ${STAGE}
            --namespace      ${G2_NAMESPACE}
            --spv ${SPV} --msl ${MSL} --desc ${DESC})
        set(BUILD_DEPS ${ABS_SRC} ${MD_GPU_SHADER_DEPS} ${G2_DEPENDS} ${SCRIPT_DEPS})
        # DEPENDS_EXPLICIT_ONLY below: a shader build reads only its sources and
        # the scripts, never the output of another target. Without it Ninja would
        # hold every command back until all of the target's dependencies are
        # built -- for md_unittest, all of mdlib.

        if (USE_METALLIB)
            set(AIR      ${GEN_DIR}/${STEM}_${ENTRY}.air)
            set(METALLIB ${GEN_DIR}/${STEM}_${ENTRY}.metallib)
            add_custom_command(
                OUTPUT  ${SPV} ${MSL} ${DESC}
                COMMAND ${BUILD_CMD} ${ABS_SRC} ${ENTRY}
                DEPENDS ${BUILD_DEPS}
                DEPENDS_EXPLICIT_ONLY
                COMMENT "slangc: ${STEM}.slang [${ENTRY}]"
                VERBATIM
            )
            add_custom_command(
                OUTPUT  ${METALLIB} ${C_FILE}
                COMMAND ${MD_GPU_XCRUN_EXECUTABLE} -sdk macosx metal -c ${MSL} -o ${AIR}
                COMMAND ${MD_GPU_XCRUN_EXECUTABLE} -sdk macosx metallib ${AIR} -o ${METALLIB}
                COMMAND ${Python3_EXECUTABLE} ${BIN2C_SCRIPT} ${METALLIB} ${C_FILE} ${SYMBOL}
                DEPENDS ${MSL} ${BIN2C_SCRIPT}
                DEPENDS_EXPLICIT_ONLY
                COMMENT "metallib: ${STEM}_${ENTRY}.metal -> ${STEM}_${ENTRY}.metallib"
                VERBATIM
            )
        else()
            add_custom_command(
                OUTPUT  ${SPV} ${MSL} ${DESC} ${C_FILE}
                COMMAND ${BUILD_CMD} --embed ${C_FILE} ${ABS_SRC} ${ENTRY}
                DEPENDS ${BUILD_DEPS}
                DEPENDS_EXPLICIT_ONLY
                COMMENT "slangc: ${STEM}.slang [${ENTRY}]"
                VERBATIM
            )
        endif()
        target_sources(${G2_TARGET} PRIVATE ${C_FILE} ${DESC})

        string(APPEND DECLS
"// ${STEM}.slang [${ENTRY}]
extern const uint8_t ${SYMBOL}_start[];
extern const size_t ${SYMBOL}_byte_size;
static inline size_t ${SYMBOL}_size(void) {
    return ${SYMBOL}_byte_size;
}

")
        string(APPEND DESCS "#include \"${STEM}_${ENTRY}_desc.inl\"\n")
    endforeach()

    # The code symbols, then the descriptors built on them. Written at configure
    # time, and only when its content changes, so that reconfiguring does not
    # recompile everything that includes it.
    file(CONFIGURE OUTPUT ${GEN_DIR}/${OUT_HEADER} @ONLY CONTENT
"/* Generated by compile_gpu_shaders() from ${STEM}.slang -- do not edit. */
#pragma once
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern \"C\" {
#endif

${DECLS}#ifdef __cplusplus
}
#endif

${DESCS}")
    target_include_directories(${G2_TARGET} PRIVATE ${GEN_DIR})
endfunction()
