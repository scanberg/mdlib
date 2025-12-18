# We do not want to have a hard dependency on VCPKG
# As it is only used for big external libs
# So if VCPKG_ROOT is defined in the environment, we use the cmake toolchain of vcpkg
if (NOT DEFINED CMAKE_TOOLCHAIN_FILE AND DEFINED ENV{VCPKG_ROOT})
    if(WIN32)
        set(VCPKG_TARGET_TRIPLET "x64-windows-static" CACHE STRING "VCPKG Target Triplet to use")
    endif()

    set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake")
    message(STATUS "VCPKG_ROOT is defined, using vcpkg-toolchain: " ${CMAKE_TOOLCHAIN_FILE})
endif()