# Usage:
# cmake -DINPUT=... -DOUTPUT=... -DSYMBOL=... -P BinToC.cmake

if (NOT DEFINED INPUT OR NOT DEFINED OUTPUT OR NOT DEFINED SYMBOL)
  message(FATAL_ERROR "BinToC.cmake requires INPUT, OUTPUT, SYMBOL")
endif()

if (NOT EXISTS "${INPUT}")
  message(FATAL_ERROR "BinToC.cmake: input file does not exist: ${INPUT}")
endif()

file(SIZE "${INPUT}" INPUT_SIZE)
if (INPUT_SIZE EQUAL 0)
  message(FATAL_ERROR "BinToC.cmake: input file is empty: ${INPUT}")
endif()

# Read binary as hex string (2 chars per byte)
file(READ "${INPUT}" HEX_CONTENT HEX)
string(TOLOWER "${HEX_CONTENT}" HEX_CONTENT)

# Build array data by iterating 2 hex chars (one uint8_t) at a time
string(LENGTH "${HEX_CONTENT}" HEX_LEN)
math(EXPR NUM_BYTES "${HEX_LEN} / 2")
set(ARRAY_DATA "")
set(POS 0)
foreach(I RANGE 1 ${NUM_BYTES})
  string(SUBSTRING "${HEX_CONTENT}" ${POS} 2 BYTE)
  if (ARRAY_DATA)
    string(APPEND ARRAY_DATA ",0x${BYTE}")
  else()
    set(ARRAY_DATA "0x${BYTE}")
  endif()
  math(EXPR POS "${POS} + 2")
endforeach()

file(WRITE "${OUTPUT}"
"#include <stdint.h>\n#include <stddef.h>\n\n"
"const uint8_t ${SYMBOL}_start[] = {\n    ${ARRAY_DATA}\n};\n\n"
"const size_t ${SYMBOL}_byte_size = sizeof(${SYMBOL}_start);\n\n"
"_Static_assert(sizeof(${SYMBOL}_start) > 0, \"${SYMBOL}_start is empty - check the source binary file\");\n"
)