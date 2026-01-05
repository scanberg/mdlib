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
            enable_language(ASM_MASM)
            set(ASM ${GEN_DIR}/${SYMBOL}.asm)

            file(WRITE ${ASM}
".data
PUBLIC ${SYMBOL}_start
PUBLIC ${SYMBOL}_end

${SYMBOL}_start LABEL BYTE
    INCBIN \"${ABS}\"
${SYMBOL}_end LABEL BYTE
")

            target_sources(${EMBED_TARGET} PRIVATE ${ASM})

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
