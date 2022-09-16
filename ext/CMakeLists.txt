add_library(mdext
    dcd/endianswap.h
    dcd/fastio.h
    dcd/readdcd.h
    fastlz/fastlz.h
    fastlz/fastlz.c
    gl3w/gl3w.h
    gl3w/gl3w.c
    gl3w/glcorearb.h
    svd3/svd3.h
    svd3/svd3.c
    xtc/xdrfile.c
    xtc/xdrfile.h
    xtc/xdrfile_trr.c
    xtc/xdrfile_trr.h
    xtc/xdrfile_xtc.c
    xtc/xdrfile_xtc.h
)

target_include_directories(mdext
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

if (MSVC)
    target_compile_definitions(mdext PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(mdext PRIVATE "/W1")
endif()

if (CMAKE_COMPILER_IS_GNUCC)
    find_package(Threads REQUIRED)
    target_compile_definitions(mdext PRIVATE _GNU_SOURCE)
endif()

target_link_libraries(mdext PRIVATE ${MOLD_LIBS})
target_compile_features(mdext PRIVATE c_std_11)