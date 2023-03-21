string(REGEX REPLACE "( -DNDEBUG$|-DNDEBUG )" "" CMAKE_CXX_FLAGS_BETATEST "${CMAKE_CXX_FLAGS_RELEASE}" ) 
string(REGEX REPLACE "( -DNDEBUG$|-DNDEBUG )" "" CMAKE_C_FLAGS_BETATEST "${CMAKE_C_FLAGS_RELEASE}" )    
string(REGEX REPLACE "-DNDEBUG " "" CMAKE_CXX_FLAGS_RELWITHDEBUG "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -DDEBUG" )
string(REGEX REPLACE "-DNDEBUG " "" CMAKE_C_FLAGS_RELWITHDEBUG "${CMAKE_C_FLAGS_RELWITHDEBINFO} -DDEBUG" )

string(APPEND CMAKE_CXX_FLAGS_BETATEST " -g")
string(APPEND CMAKE_C_FLAGS_BETATEST " -g")

include(GNUInstallDirs)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

function(create_copy_resource_dir_target target_name SRC_DIR DST_DIR)
    set(DST_FILES)
    file(GLOB_RECURSE COPY_FILES RELATIVE "${SRC_DIR}/" "${SRC_DIR}/*")
    foreach(FILE ${COPY_FILES})
      set(SRC "${SRC_DIR}/${FILE}")
      set(DST "${DST_DIR}/${FILE}")
      add_custom_command(
        OUTPUT  ${DST}
        COMMAND ${CMAKE_COMMAND} -E make_directory ${DST_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy ${SRC} ${DST}
        DEPENDS ${SRC}
      )
      list(APPEND DST_FILES ${DST})
    endforeach(FILE)
    add_custom_target(${target_name} DEPENDS ${DST_FILES})
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
        file(APPEND ${output_file} "const unsigned char ${ident}[] = {${filedata}};\nconst unsigned ${ident}_size = sizeof(${ident});\n")
    endforeach()
endfunction()