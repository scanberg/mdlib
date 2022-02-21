add_library(mold_ext
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

target_include_directories(mold_ext
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

if (MSVC)
    target_compile_definitions(mold_ext PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(mold_ext PRIVATE "/W1")
endif()

if (CMAKE_COMPILER_IS_GNUCC)
find_package(Threads REQUIRED)
    target_compile_definitions(mold_ext PRIVATE _GNU_SOURCE)
    list(APPEND MOLD_LIBS m Threads::Threads)
endif()


target_link_libraries(mold_ext PRIVATE ${MOLD_LIBS} ${VIAMD_STDLIBS})
target_compile_features(imgui PRIVATE c_std_11)