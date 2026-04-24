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
  string(REPEAT "0" ${PAD} PAD_ZEROS)
  string(APPEND HEX_CONTENT "${PAD_ZEROS}")
  string(LENGTH "${HEX_CONTENT}" HEX_LEN)
endif()

# Build array data by iterating 8 hex chars (one uint32_t) at a time
math(EXPR NUM_WORDS "${HEX_LEN} / 8")
set(ARRAY_DATA "")
set(POS 0)
foreach(I RANGE 1 ${NUM_WORDS})
  string(SUBSTRING "${HEX_CONTENT}" ${POS} 8 WORD)
  if (ARRAY_DATA)
    string(APPEND ARRAY_DATA ",0x${WORD}")
  else()
    set(ARRAY_DATA "0x${WORD}")
  endif()
  math(EXPR POS "${POS} + 8")
endforeach()

file(WRITE "${OUTPUT}"
"#include <stdint.h>\n\n"
"const uint32_t ${SYMBOL}_start[] = {\n    ${ARRAY_DATA}\n};\n\n"
"const uint32_t * const ${SYMBOL}_end = ${SYMBOL}_start + (sizeof(${SYMBOL}_start) / sizeof(${SYMBOL}_start[0]));\n"
)