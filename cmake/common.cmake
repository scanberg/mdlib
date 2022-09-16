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

# add_compile_options_config(<CONFIG> <option> ...)
function(add_compile_options_config CONFIG)
    foreach(opt ${ARGN})
        add_compile_options("$<$<CONFIG:${CONFIG}>:${opt}>")
    endforeach()
endfunction()

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