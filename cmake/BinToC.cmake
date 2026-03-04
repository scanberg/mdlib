# Usage:
# cmake -DINPUT=... -DOUTPUT=... -DSYMBOL=... -P BinToC.cmake

if (NOT DEFINED INPUT OR NOT DEFINED OUTPUT OR NOT DEFINED SYMBOL)
  message(FATAL_ERROR "BinToC.cmake requires INPUT, OUTPUT, SYMBOL")
endif()

# Read binary as hex string (2 chars per byte)
file(READ "${INPUT}" HEX_CONTENT HEX)
string(TOLOWER "${HEX_CONTENT}" HEX_CONTENT)

# Pad to 32-bit (4 bytes => 8 hex chars). If not aligned, append '0's.
string(LENGTH "${HEX_CONTENT}" HEX_LEN)
math(EXPR REM "${HEX_LEN} % 8")
if (NOT REM EQUAL 0)
  math(EXPR PAD "8 - ${REM}")
  string(RANDOM LENGTH ${PAD} ALPHABET "0" PAD_ZEROS) # produces all '0'
  string(APPEND HEX_CONTENT "${PAD_ZEROS}")
endif()

# Convert each 32-bit word (8 hex digits) to "0x........,"
string(REGEX REPLACE "([0-9a-f]{8})" "0x\\1," ARRAY_DATA "${HEX_CONTENT}")
string(REGEX REPLACE ",$" "" ARRAY_DATA "${ARRAY_DATA}")

file(WRITE "${OUTPUT}"
"#include <stdint.h>\n\n"
"const uint32_t ${SYMBOL}_start[] = {\n    ${ARRAY_DATA}\n};\n\n"
"const uint32_t * const ${SYMBOL}_end = ${SYMBOL}_start + (sizeof(${SYMBOL}_start) / sizeof(${SYMBOL}_start[0]));\n"
)