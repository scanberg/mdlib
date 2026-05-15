# ========================================================================
# ResolveVulkan.cmake
#
# Hybrid Vulkan resolver:
#
#   Tier 1: Vulkan SDK / system Vulkan (preferred)
#   Tier 2: system libvulkan + vendored headers/volk
#   Tier 3: fully vendored (headers + volk via FetchContent)
#
# Exposes:
#
#   vulkan_backend (INTERFACE target)
#
# Optional variables:
#
#   VULKAN_ENABLE_VALIDATION (default ON if SDK found)
#
# Behavior:
#   - Uses Vulkan SDK if available (find_package)
#   - Otherwise fetches Vulkan-Headers + volk
#   - Always provides a consistent interface target
#
# ========================================================================

include_guard(GLOBAL)

include(FetchContent)

# ------------------------------------------------------------------------
# Main resolver function
# ------------------------------------------------------------------------

function(resolve_vulkan)

    message(STATUS "")
    message(STATUS "Resolving Vulkan dependency...")
    message(STATUS "")

    # --------------------------------------------------------------------
    # 1. Try system Vulkan / SDK
    # --------------------------------------------------------------------

    find_package(Vulkan QUIET)

    if(Vulkan_FOUND)

        set(_use_system_vulkan ON PARENT_SCOPE)

        message(STATUS "Vulkan: Found system/Vulkan SDK")
        message(STATUS "Vulkan include: ${Vulkan_INCLUDE_DIRS}")
        message(STATUS "Vulkan library: ${Vulkan_LIBRARIES}")

    else()

        set(_use_system_vulkan OFF PARENT_SCOPE)

        message(STATUS "Vulkan: Not found system-wide, using vendored fallback")

    endif()

    # --------------------------------------------------------------------
    # 2. Create unified interface target
    # --------------------------------------------------------------------

    if(NOT TARGET vulkan_backend)
        add_library(vulkan_backend INTERFACE)
    endif()

    # --------------------------------------------------------------------
    # 3. System Vulkan path
    # --------------------------------------------------------------------

    if(Vulkan_FOUND)

        target_link_libraries(vulkan_backend INTERFACE Vulkan::Vulkan)

        # Validation layers are only meaningful when SDK is present
        option(Vulkan_ENABLE_VALIDATION "Enable Vulkan validation layers" ON)

        message(STATUS "Vulkan validation layers: ${Vulkan_ENABLE_VALIDATION}")

    else()

        set(Vulkan_ENABLE_VALIDATION OFF CACHE BOOL "" FORCE)

        # ----------------------------------------------------------------
        # 4. Fetch minimal Vulkan stack
        # ----------------------------------------------------------------

        message(STATUS "Fetching Vulkan-Headers + volk...")

        FetchContent_Declare(
            VulkanHeaders
            GIT_REPOSITORY https://github.com/KhronosGroup/Vulkan-Headers.git
            GIT_TAG vulkan-sdk-1.4.350.0
        )

        FetchContent_Declare(
            volk
            GIT_REPOSITORY https://github.com/zeux/volk.git
            GIT_TAG vulkan-sdk-1.4.350.0
        )

        FetchContent_MakeAvailable(VulkanHeaders volk)

        # ----------------------------------------------------------------
        # 5. Configure volk usage
        # ----------------------------------------------------------------

        target_compile_definitions(vulkan_backend INTERFACE
            VK_NO_PROTOTYPES
        )

        target_link_libraries(vulkan_backend INTERFACE volk)

        target_include_directories(vulkan_backend INTERFACE
            ${vulkanheaders_SOURCE_DIR}/include
        )

    endif()

    # --------------------------------------------------------------------
    # 6. Debug output
    # --------------------------------------------------------------------

    message(STATUS "")
    message(STATUS "Vulkan resolution complete")
    message(STATUS "System Vulkan: ${Vulkan_FOUND}")
    message(STATUS "Validation enabled: ${Vulkan_ENABLE_VALIDATION}")
    message(STATUS "")

endfunction()