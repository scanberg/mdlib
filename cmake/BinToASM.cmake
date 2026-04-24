# Usage:
# cmake -DINPUT=... -DOUTPUT=... -DSYMBOL=... [-DAPPLE_ABI=ON] -P BinToASM.cmake
#
# Generates a .S assembly file that embeds a binary file using .incbin.
# On Apple (Mach-O), C symbols are prefixed with underscore in assembly.

if (NOT DEFINED INPUT OR NOT DEFINED OUTPUT OR NOT DEFINED SYMBOL)
  message(FATAL_ERROR "BinToASM.cmake requires INPUT, OUTPUT, SYMBOL")
endif()

if (APPLE_ABI)
    set(SECTION ".section __DATA,__const")
    set(PREFIX "_")
else()
    set(SECTION ".section .rodata")
    set(PREFIX "")
endif()

file(WRITE "${OUTPUT}"
"${SECTION}\n"
".p2align 2\n"
".global ${PREFIX}${SYMBOL}_start\n"
"${PREFIX}${SYMBOL}_start:\n"
".incbin \"${INPUT}\"\n"
".p2align 2\n"
".global ${PREFIX}${SYMBOL}_end\n"
"${PREFIX}${SYMBOL}_end:\n"
)
