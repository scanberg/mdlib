enable_language(C CXX)
include(GNUInstallDirs)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

# Mirrors SRC_DIR into DST_DIR at build time, for executables that find their resources
# relative to themselves. Added, edited and removed files all reach DST_DIR, and nothing
# else is written there, so the folder can ship as it is.
function(create_copy_resource_dir_target target_name SRC_DIR DST_DIR)
    # DST_DIR is deleted and copied anew, so refuse anything outside the build tree.
    cmake_path(IS_PREFIX CMAKE_BINARY_DIR "${DST_DIR}" NORMALIZE inside)
    if (NOT inside)
        message(FATAL_ERROR "create_copy_resource_dir_target: ${DST_DIR} is outside the build tree")
    endif()

    file(GLOB_RECURSE files CONFIGURE_DEPENDS LIST_DIRECTORIES false RELATIVE "${SRC_DIR}" "${SRC_DIR}/*")
    if (NOT files)
        message(FATAL_ERROR "create_copy_resource_dir_target: ${SRC_DIR} has no files")
    endif()
    list(TRANSFORM files PREPEND "${SRC_DIR}/" OUTPUT_VARIABLE src_files)
    list(TRANSFORM files PREPEND "${DST_DIR}/" OUTPUT_VARIABLE dst_files)

    # The file list, rewritten only when a file is added or removed: a removal leaves every
    # remaining output up to date, and this is what makes the copy run and drop the file.
    # Rewriting it on every configure would make every reconfigure copy again.
    set(manifest "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${target_name}.files")
    list(JOIN files "\n" content)
    file(CONFIGURE OUTPUT "${manifest}" CONTENT "${content}\n" @ONLY)

    add_custom_command(
        OUTPUT  ${dst_files}
        COMMAND ${CMAKE_COMMAND} -E rm -rf "${DST_DIR}"
        COMMAND ${CMAKE_COMMAND} -E copy_directory "${SRC_DIR}" "${DST_DIR}"
        DEPENDS ${src_files} "${manifest}"
        COMMENT "Copying ${target_name} resources"
        VERBATIM
    )
    add_custom_target(${target_name} DEPENDS ${dst_files})
    set_target_properties(${target_name} PROPERTIES FOLDER "resources")
endfunction()

# https://stackoverflow.com/questions/11813271/embed-resources-eg-shader-code-images-into-executable-library-with-cmake
# Creates C resources file from files in given directory
function(create_resources input_list output_file)
    # The resources are embedded when CMake configures: a change to one of the files has to configure again, or a build
    # keeps the old contents without a word (the shaders, which are edited without touching any CMake file)
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${input_list})
    # Create empty output file
    file(WRITE ${output_file} "")
    message(STATUS "Writing resources to ${output_file}")
    # Iterate through input files
    foreach(bin ${input_list})
        # Get short filename
        string(REGEX MATCH "([^/]+)$" filename ${bin})
        # Replace filename spaces & extension separator for C compatibility
        string(REGEX REPLACE "\\.| |-" "_" ident ${filename})
        # Read hex data from file
        file(READ "${bin}" filedata HEX)
        # Convert hex data for C compatibility
        string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," filedata ${filedata})
        # Append data to output file
        file(APPEND ${output_file} "const unsigned char ${ident}[] = {${filedata}0x0};\nconst unsigned int ${ident}_size = sizeof(${ident});\n")
    endforeach()
endfunction()
