enable_language(C CXX)
include(GNUInstallDirs)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

function(create_copy_resource_dir_target target_name SRC_DIR DST_DIR)
    if (CMAKE_VERSION VERSION_LESS 3.26)
      message(FATAL_ERROR "create_copy_resource_dir_target requires CMake 3.26 or newer")
    endif()

    file(GLOB_RECURSE COPY_FILES CONFIGURE_DEPENDS LIST_DIRECTORIES false RELATIVE "${SRC_DIR}" "${SRC_DIR}/*")

    set(SRC_FILES)
    set(DST_FILES)
    foreach(FILE IN LISTS COPY_FILES)
      list(APPEND SRC_FILES "${SRC_DIR}/${FILE}")
      list(APPEND DST_FILES "${DST_DIR}/${FILE}")
    endforeach()

    set(MANIFEST "${CMAKE_CURRENT_BINARY_DIR}/${target_name}_resources.txt")
    string(REPLACE ";" "\n" MANIFEST_CONTENT "${COPY_FILES}")
    file(WRITE "${MANIFEST}" "${MANIFEST_CONTENT}\n")

    set(STAMP "${DST_DIR}/.${target_name}.stamp")
    add_custom_command(
      OUTPUT ${DST_FILES} "${STAMP}"
      COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different "${SRC_DIR}" "${DST_DIR}"
      COMMAND ${CMAKE_COMMAND} -E touch ${DST_FILES} "${STAMP}"
      DEPENDS ${SRC_FILES} "${MANIFEST}"
      COMMENT "Copying ${target_name} resources"
      VERBATIM
    )
    add_custom_target(${target_name} DEPENDS "${STAMP}" ${DST_FILES})
endfunction()

# https://stackoverflow.com/questions/11813271/embed-resources-eg-shader-code-images-into-executable-library-with-cmake
# Creates C resources file from files in given directory
function(create_resources input_list output_file)
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
