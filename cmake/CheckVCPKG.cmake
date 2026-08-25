# We do not want to have a hard dependency on VCPKG
# As it is only used for big external libs
# So if VCPKG_ROOT is defined in the environment, we use the cmake toolchain of vcpkg
if (NOT DEFINED CMAKE_TOOLCHAIN_FILE AND DEFINED ENV{VCPKG_ROOT})
    if(WIN32)
        set(VCPKG_TARGET_TRIPLET "x64-windows-static" CACHE STRING "VCPKG Target Triplet to use")
    endif()

    # vcpkg only switches to manifest mode when it finds a vcpkg.json next to the top level
    # CMakeLists.txt. That holds when mdlib is configured directly, but not when it is consumed
    # as a subdirectory, where the manifest still lives here in the mdlib tree. Point vcpkg at it
    # explicitly so both configurations resolve the same dependencies - in particular hdf5 without
    # its default szip feature, which pulls in a libaec source fetch from a host that is regularly
    # unreachable from CI runners.
    #
    # Gated on HDF5 actually being requested: manifest mode installs the whole manifest at
    # configure time, and a build without Veloxchem support has no use for it.
    if (VIAMD_ENABLE_VELOXCHEM OR MD_ENABLE_VLX OR MD_ENABLE_HDF5)
        get_filename_component(MD_VCPKG_MANIFEST_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
        if (NOT DEFINED VCPKG_MANIFEST_DIR
            AND NOT EXISTS "${CMAKE_SOURCE_DIR}/vcpkg.json"
            AND EXISTS "${MD_VCPKG_MANIFEST_DIR}/vcpkg.json")
            set(VCPKG_MANIFEST_DIR "${MD_VCPKG_MANIFEST_DIR}")
            message(STATUS "Using mdlib vcpkg manifest: ${VCPKG_MANIFEST_DIR}/vcpkg.json")
        endif()
    endif()

    set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake")
    message(STATUS "VCPKG_ROOT is defined, using vcpkg-toolchain: " ${CMAKE_TOOLCHAIN_FILE})
endif()