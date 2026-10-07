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

# "0xab,0xcd,..." in one regex pass (appending byte by byte is quadratic in CMake: minutes for a large
# SPIR-V module), without the trailing comma.
string(REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," ARRAY_DATA "${HEX_CONTENT}")
string(REGEX REPLACE ",$" "" ARRAY_DATA "${ARRAY_DATA}")

file(WRITE "${OUTPUT}"
"#include <stdint.h>\n#include <stddef.h>\n\n"
"const uint8_t ${SYMBOL}_start[] = {\n    ${ARRAY_DATA}\n};\n\n"
"const size_t ${SYMBOL}_byte_size = sizeof(${SYMBOL}_start);\n\n"
"_Static_assert(sizeof(${SYMBOL}_start) > 0, \"${SYMBOL}_start is empty - check the source binary file\");\n"
)