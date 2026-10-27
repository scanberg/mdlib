# Writes OUT (a header defining MD_BENCH_SOURCE_REV) from the git state of the mdlib checkout at SRC:
# "<branch>@<short hash>", plus " +local changes: <files>" when tracked files differ from that commit.
# Run at build time (see CMakeLists.txt), so every benchmark run reports the source it was built from.
# The file is rewritten only when the text changes, so an unchanged checkout triggers no rebuild.

set(rev "unknown (no git)")
find_package(Git QUIET)
if (GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} -C ${SRC} rev-parse --short HEAD
                    OUTPUT_VARIABLE hash OUTPUT_STRIP_TRAILING_WHITESPACE RESULT_VARIABLE res ERROR_QUIET)
    if (res EQUAL 0)
        execute_process(COMMAND ${GIT_EXECUTABLE} -C ${SRC} rev-parse --abbrev-ref HEAD
                        OUTPUT_VARIABLE branch OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
        execute_process(COMMAND ${GIT_EXECUTABLE} -C ${SRC} --no-optional-locks status --porcelain
                                --untracked-files=no --ignore-submodules=all
                        OUTPUT_VARIABLE dirty OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
        set(rev "${branch}@${hash}")
        if (NOT dirty STREQUAL "")
            string(REPLACE "\n" ";" lines "${dirty}")
            set(files "")
            foreach (line IN LISTS lines)
                string(SUBSTRING "${line}" 3 -1 path)
                list(APPEND files "${path}")
            endforeach()
            list(LENGTH files nfiles)
            if (nfiles GREATER 4)
                list(SUBLIST files 0 4 files)
                list(APPEND files "...")
            endif()
            list(JOIN files ", " files)
            set(rev "${rev} +local changes: ${files}")
        endif()
    endif()
endif()

string(REPLACE "\\" "\\\\" rev "${rev}")
string(REPLACE "\"" "\\\"" rev "${rev}")
set(text "#define MD_BENCH_SOURCE_REV \"${rev}\"\n")
if (EXISTS ${OUT})
    file(READ ${OUT} old)
endif()
if (NOT "${old}" STREQUAL "${text}")
    file(WRITE ${OUT} "${text}")
endif()
