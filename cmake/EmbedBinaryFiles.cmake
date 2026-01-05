function(embed_binary_files)
    set(options)
    set(oneValueArgs TARGET NAMESPACE OUTPUT)
    set(multiValueArgs FILES)
    cmake_parse_arguments(EMBED "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if (NOT EMBED_TARGET OR NOT EMBED_NAMESPACE OR NOT EMBED_OUTPUT)
        message(FATAL_ERROR "embed_binary_files: TARGET, NAMESPACE, OUTPUT required")
    endif()

    set(GEN_DIR ${CMAKE_CURRENT_BINARY_DIR}/gen)
    file(MAKE_DIRECTORY ${GEN_DIR})

    set(OUT_HEADER ${GEN_DIR}/${EMBED_OUTPUT})

    file(WRITE ${OUT_HEADER}
"#pragma once
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern \"C\" {
#endif

")

    foreach(FILE ${EMBED_FILES})
        get_filename_component(NAME_WITH_EXT ${FILE} NAME)
        get_filename_component(NAME ${FILE} NAME_WE)
        string(REPLACE "." "_" SYMBOL ${NAME})
        set(SYMBOL ${EMBED_NAMESPACE}_${SYMBOL})

        get_filename_component(ABS ${FILE} ABSOLUTE)

        file(APPEND ${OUT_HEADER}
"// From file: ${NAME_WITH_EXT}
extern const uint8_t ${SYMBOL}_start[];
extern const uint8_t ${SYMBOL}_end[];
static inline size_t ${SYMBOL}_size(void) {
    return (size_t)(${SYMBOL}_end - ${SYMBOL}_start);
}

")

        if (MSVC)
            set(C_FILE ${GEN_DIR}/${SYMBOL}.c)

            file(READ ${ABS} HEX_CONTENT HEX)
            
            string(LENGTH "${HEX_CONTENT}" HEX_LENGTH)
            math(EXPR BYTE_COUNT "${HEX_LENGTH} / 2")
            
            # Convert hex string to C array format
            string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," ARRAY_DATA "${HEX_CONTENT}")
            string(REGEX REPLACE ",$" "" ARRAY_DATA "${ARRAY_DATA}")
            
            file(WRITE ${C_FILE}
"#include <stdint.h>

const uint8_t ${SYMBOL}_start[] = {
    ${ARRAY_DATA}
};

const uint8_t ${SYMBOL}_end[] = {};
")

            target_sources(${EMBED_TARGET} PRIVATE ${C_FILE})

        elseif (APPLE)
            target_link_options(${EMBED_TARGET} PRIVATE
                "-Wl,-sectcreate,__DATA,__${SYMBOL},${ABS}"
            )

        else()
            set(OBJ ${GEN_DIR}/${SYMBOL}.o)

            add_custom_command(
                OUTPUT ${OBJ}
                COMMAND ${CMAKE_OBJCOPY}
                    --input binary
                    --output elf64-x86-64
                    --binary-architecture i386:x86-64
                    ${ABS} ${OBJ}
                DEPENDS ${ABS}
            )

            target_sources(${EMBED_TARGET} PRIVATE ${OBJ})
        endif()
    endforeach()

    file(APPEND ${OUT_HEADER}
"#ifdef __cplusplus
}
#endif
")

    target_include_directories(${EMBED_TARGET} PRIVATE ${GEN_DIR})
endfunction()
