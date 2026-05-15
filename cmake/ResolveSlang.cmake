# ========================================================================
# ResolveSlang.cmake
#
# Self-contained Slang/slangc resolver for CMake projects.
#
# Features:
#   - Downloads correct Slang binary release for host platform
#   - Supports:
#       * Linux x86_64
#       * Linux aarch64
#       * macOS x86_64
#       * macOS arm64 (Apple Silicon)
#       * Windows x64
#   - Version pinning
#   - SHA256 verification support (optional)
#   - Automatic extraction
#   - Cache reuse
#   - Exposes:
#
#         SLANG_ROOT
#         SLANG_EXECUTABLE
#
#   - Creates imported executable target:
#
#         slangc
#
# Usage:
#
#   include(cmake/ResolveSlang.cmake)
#
#   resolve_slang(
#       VERSION 2025.8.1
#   )
#
#   message(STATUS "slangc: ${SLANG_EXECUTABLE}")
#
# Example:
#
#   add_custom_command(
#       OUTPUT shader.spv
#       COMMAND ${SLANG_EXECUTABLE}
#           shaders/test.slang
#           -target spirv
#           -entry computeMain
#           -o shader.spv
#       DEPENDS shaders/test.slang
#   )
#
# ========================================================================

include_guard(GLOBAL)

function(resolve_slang)

    set(options)
    set(oneValueArgs
        VERSION
        CACHE_DIR
        SHA256
    )

    cmake_parse_arguments(SLANG
        "${options}"
        "${oneValueArgs}"
        ""
        ${ARGN}
    )

    if(NOT SLANG_VERSION)
        message(FATAL_ERROR
            "resolve_slang() requires VERSION")
    endif()

    # ================================================================
    # Cache location
    # ================================================================

    if(SLANG_CACHE_DIR)
        set(_slang_cache_dir "${SLANG_CACHE_DIR}")
    else()
        set(_slang_cache_dir
            "${CMAKE_BINARY_DIR}/third_party/slang")
    endif()

    file(MAKE_DIRECTORY "${_slang_cache_dir}")

    # ================================================================
    # Detect platform
    # ================================================================

    set(_platform "")
    set(_archive_ext "")
    set(_slang_exe "slangc")

    if(WIN32)

        set(_platform "windows-x86_64")
        set(_archive_ext "zip")
        set(_slang_exe "slangc.exe")

    elseif(APPLE)

        if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")

            set(_platform "macos-aarch64")

        else()

            set(_platform "macos-x86_64")

        endif()

        set(_archive_ext "tar.gz")

    elseif(UNIX)

        if(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")

            set(_platform "linux-aarch64")

        else()

            set(_platform "linux-x86_64")

        endif()

        set(_archive_ext "tar.gz")

    else()

        message(FATAL_ERROR
            "Unsupported platform")

    endif()

    # ================================================================
    # Release URL
    # ================================================================

    set(_archive_name
        "slang-${SLANG_VERSION}-${_platform}.${_archive_ext}")

    set(_download_url
        "https://github.com/shader-slang/slang/releases/download/v${SLANG_VERSION}/${_archive_name}")

    # ================================================================
    # Paths
    # ================================================================

    set(_archive_path
        "${_slang_cache_dir}/${_archive_name}")

    set(_extract_dir
        "${_slang_cache_dir}/${SLANG_VERSION}-${_platform}")

    set(_bin_dir
        "${_extract_dir}/bin")

    set(_slang_executable
        "${_bin_dir}/${_slang_exe}")

    # ================================================================
    # Download if missing
    # ================================================================

    if(NOT EXISTS "${_slang_executable}")

        message(STATUS "")
        message(STATUS "Resolving Slang ${SLANG_VERSION}")
        message(STATUS "Platform: ${_platform}")
        message(STATUS "Downloading: ${_download_url}")
        message(STATUS "")

        if(NOT EXISTS "${_archive_path}")

            if(SLANG_SHA256)

                file(
                    DOWNLOAD
                    "${_download_url}"
                    "${_archive_path}"
                    EXPECTED_HASH SHA256=${SLANG_SHA256}
                    SHOW_PROGRESS
                    STATUS _download_status
                )

            else()

                file(
                    DOWNLOAD
                    "${_download_url}"
                    "${_archive_path}"
                    SHOW_PROGRESS
                    STATUS _download_status
                )

            endif()

            list(GET _download_status 0 _download_code)

            if(NOT _download_code EQUAL 0)

                list(GET _download_status 1 _download_error)

                message(FATAL_ERROR
                    "Failed downloading Slang:\n"
                    "${_download_error}")

            endif()

        endif()

        # ============================================================
        # Clean extraction dir
        # ============================================================

        file(REMOVE_RECURSE "${_extract_dir}")
        file(MAKE_DIRECTORY "${_extract_dir}")

        # ============================================================
        # Extract
        # ============================================================

        message(STATUS "Extracting Slang...")

        file(
            ARCHIVE_EXTRACT
            INPUT "${_archive_path}"
            DESTINATION "${_extract_dir}"
        )

        # ============================================================
        # Some releases contain nested top-level directories.
        # Normalize layout if needed.
        # ============================================================

        if(NOT EXISTS "${_slang_executable}")

            file(GLOB _children
                RELATIVE "${_extract_dir}"
                "${_extract_dir}/*")

            list(LENGTH _children _child_count)

            if(_child_count EQUAL 1)

                list(GET _children 0 _child)

                if(EXISTS "${_extract_dir}/${_child}/bin/${_slang_exe}")

                    set(_extract_dir
                        "${_extract_dir}/${_child}")

                    set(_bin_dir
                        "${_extract_dir}/bin")

                    set(_slang_executable
                        "${_bin_dir}/${_slang_exe}")

                endif()

            endif()

        endif()

        # ============================================================
        # Validation
        # ============================================================

        if(NOT EXISTS "${_slang_executable}")

            message(FATAL_ERROR
                "Could not locate slang executable after extraction:\n"
                "${_slang_executable}")

        endif()

        # ============================================================
        # Ensure executable permissions
        # ============================================================

        if(UNIX)

            execute_process(
                COMMAND chmod +x "${_slang_executable}"
            )

        endif()

    endif()

    # ================================================================
    # Export variables
    # ================================================================

    set(SLANG_ROOT
        "${_extract_dir}"
        CACHE PATH
        "Slang root directory"
        FORCE)

    set(SLANG_EXECUTABLE
        "${_slang_executable}"
        CACHE FILEPATH
        "Path to slangc executable"
        FORCE)

    # ================================================================
    # Imported executable target
    # ================================================================

    if(NOT TARGET slangc)

        add_executable(slangc IMPORTED GLOBAL)

        set_target_properties(slangc PROPERTIES
            IMPORTED_LOCATION "${SLANG_EXECUTABLE}"
        )

    endif()

    # ================================================================
    # Status
    # ================================================================

    message(STATUS "")
    message(STATUS "Slang resolved successfully")
    message(STATUS "SLANG_ROOT       = ${SLANG_ROOT}")
    message(STATUS "SLANG_EXECUTABLE = ${SLANG_EXECUTABLE}")
    message(STATUS "")

endfunction()