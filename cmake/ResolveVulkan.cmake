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

function(resolve_vulkan)

    message(STATUS "Resolving Vulkan dependency (volk loader)...")

    # 1. Headers: prefer system/SDK, fall back to FetchContent
    find_package(Vulkan QUIET)

    if(NOT TARGET vulkan_backend)
        add_library(vulkan_backend INTERFACE)
    endif()

    if(Vulkan_FOUND)
        message(STATUS "Vulkan headers: ${Vulkan_INCLUDE_DIRS} (system/SDK)")
        target_include_directories(vulkan_backend INTERFACE ${Vulkan_INCLUDE_DIRS})
    else()
        message(STATUS "Vulkan headers: fetching Vulkan-Headers")
        FetchContent_Declare(
            VulkanHeaders
            GIT_REPOSITORY https://github.com/KhronosGroup/Vulkan-Headers.git
            GIT_TAG vulkan-sdk-1.4.350.0
        )
        FetchContent_MakeAvailable(VulkanHeaders)
        target_include_directories(vulkan_backend INTERFACE
            ${vulkanheaders_SOURCE_DIR}/include)
    endif()

    # 2. Loader: always volk, fetched and compiled in
    FetchContent_Declare(
        volk
        GIT_REPOSITORY https://github.com/zeux/volk.git
        GIT_TAG vulkan-sdk-1.4.350.0
    )
    FetchContent_MakeAvailable(volk)

    # 3. WSI platform defines must be set when volk.c is compiled
    set(_volk_platform_defs)
    if(WIN32)
        list(APPEND _volk_platform_defs VK_USE_PLATFORM_WIN32_KHR)
    elseif(APPLE)
        list(APPEND _volk_platform_defs VK_USE_PLATFORM_METAL_EXT)
    elseif(UNIX)
        # viamd uses GLFW X11 only (GLFW_BUILD_WAYLAND OFF)
        list(APPEND _volk_platform_defs VK_USE_PLATFORM_XLIB_KHR)
        # Add VK_USE_PLATFORM_WAYLAND_KHR here if you ever enable Wayland.
    endif()
    target_compile_definitions(volk PRIVATE ${_volk_platform_defs})

    # 4. Consumers must NOT include vulkan.h directly — only volk.h
    target_compile_definitions(vulkan_backend INTERFACE
        VK_NO_PROTOTYPES
        ${_volk_platform_defs})

    target_link_libraries(vulkan_backend INTERFACE volk)

    # 5. Validation layers: only meaningful if the SDK/loader provides them
    option(Vulkan_ENABLE_VALIDATION "Enable Vulkan validation layers" ${Vulkan_FOUND})
    message(STATUS "Vulkan validation layers: ${Vulkan_ENABLE_VALIDATION}")

endfunction()