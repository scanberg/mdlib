# Usage:
# cmake -DINPUT=... -DOUTPUT=... -DSYMBOL=... -DSOURCE=... -P SlangReflectionToHeader.cmake

if (NOT DEFINED INPUT OR NOT DEFINED OUTPUT OR NOT DEFINED SYMBOL)
  message(FATAL_ERROR "SlangReflectionToHeader.cmake requires INPUT, OUTPUT, SYMBOL")
endif()

if (NOT DEFINED SOURCE)
  set(SOURCE "${INPUT}")
endif()

if (NOT EXISTS "${INPUT}")
  message(FATAL_ERROR "SlangReflectionToHeader.cmake: input file does not exist: ${INPUT}")
endif()

function(slang_reflect_sanitize_identifier OUT IN)
  string(REGEX REPLACE "[^A-Za-z0-9_]" "_" IDENT "${IN}")
  if (IDENT MATCHES "^[0-9]")
    set(IDENT "_${IDENT}")
  endif()
  set(${OUT} "${IDENT}" PARENT_SCOPE)
endfunction()

function(slang_reflect_binding_kind OUT KIND)
  if (KIND STREQUAL "constantBuffer")
    set(VALUE 1)
  elseif (KIND STREQUAL "pushConstantBuffer")
    set(VALUE 2)
  elseif (KIND STREQUAL "shaderResource")
    set(VALUE 3)
  elseif (KIND STREQUAL "descriptorTableSlot")
    set(VALUE 4)
  elseif (KIND STREQUAL "sampler")
    set(VALUE 5)
  elseif (KIND STREQUAL "uniform")
    set(VALUE 6)
  else()
    set(VALUE 0)
  endif()
  set(${OUT} ${VALUE} PARENT_SCOPE)
endfunction()

function(slang_reflect_type_kind OUT KIND)
  if (KIND STREQUAL "scalar")
    set(VALUE 1)
  elseif (KIND STREQUAL "vector")
    set(VALUE 2)
  elseif (KIND STREQUAL "matrix")
    set(VALUE 3)
  elseif (KIND STREQUAL "struct")
    set(VALUE 4)
  elseif (KIND STREQUAL "pointer")
    set(VALUE 5)
  elseif (KIND STREQUAL "resource")
    set(VALUE 6)
  elseif (KIND STREQUAL "constantBuffer")
    set(VALUE 7)
  else()
    set(VALUE 0)
  endif()
  set(${OUT} ${VALUE} PARENT_SCOPE)
endfunction()

function(slang_reflect_scalar_c_type OUT SCALAR_TYPE)
  if (SCALAR_TYPE STREQUAL "float32" OR SCALAR_TYPE STREQUAL "float")
    set(VALUE "float")
  elseif (SCALAR_TYPE STREQUAL "float64" OR SCALAR_TYPE STREQUAL "double")
    set(VALUE "double")
  elseif (SCALAR_TYPE STREQUAL "uint32" OR SCALAR_TYPE STREQUAL "uint")
    set(VALUE "uint32_t")
  elseif (SCALAR_TYPE STREQUAL "int32" OR SCALAR_TYPE STREQUAL "int")
    set(VALUE "int32_t")
  elseif (SCALAR_TYPE STREQUAL "uint64")
    set(VALUE "uint64_t")
  elseif (SCALAR_TYPE STREQUAL "int64")
    set(VALUE "int64_t")
  elseif (SCALAR_TYPE STREQUAL "bool")
    set(VALUE "uint32_t")
  else()
    set(VALUE "uint8_t")
  endif()
  set(${OUT} "${VALUE}" PARENT_SCOPE)
endfunction()

function(slang_reflect_field_decl OUT PARAM_INDEX FIELD_INDEX FIELD_SYMBOL FIELD_SIZE)
  string(JSON FIELD_KIND ERROR_VARIABLE FIELD_KIND_ERROR GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementType fields ${FIELD_INDEX} type kind)
  if (FIELD_KIND_ERROR)
    set(${OUT} "uint8_t ${FIELD_SYMBOL}[${FIELD_SIZE}];" PARENT_SCOPE)
    return()
  endif()

  if (FIELD_KIND STREQUAL "scalar")
    string(JSON SCALAR_TYPE ERROR_VARIABLE SCALAR_TYPE_ERROR GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementType fields ${FIELD_INDEX} type scalarType)
    if (SCALAR_TYPE_ERROR)
      set(SCALAR_TYPE "unknown")
    endif()
    slang_reflect_scalar_c_type(C_TYPE "${SCALAR_TYPE}")
    set(${OUT} "${C_TYPE} ${FIELD_SYMBOL};" PARENT_SCOPE)
  elseif (FIELD_KIND STREQUAL "vector")
    string(JSON ELEMENT_COUNT ERROR_VARIABLE ELEMENT_COUNT_ERROR GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementType fields ${FIELD_INDEX} type elementCount)
    string(JSON SCALAR_TYPE ERROR_VARIABLE SCALAR_TYPE_ERROR GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementType fields ${FIELD_INDEX} type elementType scalarType)
    if (ELEMENT_COUNT_ERROR)
      set(ELEMENT_COUNT 1)
    endif()
    if (SCALAR_TYPE_ERROR)
      set(SCALAR_TYPE "unknown")
    endif()
    slang_reflect_scalar_c_type(C_TYPE "${SCALAR_TYPE}")
    set(${OUT} "${C_TYPE} ${FIELD_SYMBOL}[${ELEMENT_COUNT}];" PARENT_SCOPE)
  elseif (FIELD_KIND STREQUAL "matrix")
    string(JSON ROW_COUNT ERROR_VARIABLE ROW_COUNT_ERROR GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementType fields ${FIELD_INDEX} type rowCount)
    string(JSON COLUMN_COUNT ERROR_VARIABLE COLUMN_COUNT_ERROR GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementType fields ${FIELD_INDEX} type columnCount)
    string(JSON SCALAR_TYPE ERROR_VARIABLE SCALAR_TYPE_ERROR GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementType fields ${FIELD_INDEX} type elementType scalarType)
    if (ROW_COUNT_ERROR)
      set(ROW_COUNT 1)
    endif()
    if (COLUMN_COUNT_ERROR)
      set(COLUMN_COUNT 1)
    endif()
    if (SCALAR_TYPE_ERROR)
      set(SCALAR_TYPE "unknown")
    endif()
    slang_reflect_scalar_c_type(C_TYPE "${SCALAR_TYPE}")
    set(${OUT} "${C_TYPE} ${FIELD_SYMBOL}[${ROW_COUNT}][${COLUMN_COUNT}];" PARENT_SCOPE)
  elseif (FIELD_KIND STREQUAL "pointer")
    set(${OUT} "uint64_t ${FIELD_SYMBOL};" PARENT_SCOPE)
  else()
    set(${OUT} "uint8_t ${FIELD_SYMBOL}[${FIELD_SIZE}];" PARENT_SCOPE)
  endif()
endfunction()

file(READ "${INPUT}" REFLECTION_JSON)

file(WRITE "${OUTPUT}"
"#pragma once\n"
"#include <stdint.h>\n"
"#include <stddef.h>\n\n"
"// Generated from Slang reflection: ${SOURCE}\n\n"
"#ifndef MD_GPU_SHADER_REFLECTION_CONSTANTS_DEFINED\n"
"#define MD_GPU_SHADER_REFLECTION_CONSTANTS_DEFINED\n"
"#define MD_GPU_SHADER_BINDING_KIND_UNKNOWN 0u\n"
"#define MD_GPU_SHADER_BINDING_KIND_CONSTANT_BUFFER 1u\n"
"#define MD_GPU_SHADER_BINDING_KIND_PUSH_CONSTANT_BUFFER 2u\n"
"#define MD_GPU_SHADER_BINDING_KIND_SHADER_RESOURCE 3u\n"
"#define MD_GPU_SHADER_BINDING_KIND_DESCRIPTOR_TABLE_SLOT 4u\n"
"#define MD_GPU_SHADER_BINDING_KIND_SAMPLER 5u\n"
"#define MD_GPU_SHADER_BINDING_KIND_UNIFORM 6u\n"
"#define MD_GPU_SHADER_TYPE_KIND_UNKNOWN 0u\n"
"#define MD_GPU_SHADER_TYPE_KIND_SCALAR 1u\n"
"#define MD_GPU_SHADER_TYPE_KIND_VECTOR 2u\n"
"#define MD_GPU_SHADER_TYPE_KIND_MATRIX 3u\n"
"#define MD_GPU_SHADER_TYPE_KIND_STRUCT 4u\n"
"#define MD_GPU_SHADER_TYPE_KIND_POINTER 5u\n"
"#define MD_GPU_SHADER_TYPE_KIND_RESOURCE 6u\n"
"#define MD_GPU_SHADER_TYPE_KIND_CONSTANT_BUFFER 7u\n"
"#endif\n\n")

file(APPEND "${OUTPUT}"
"#ifndef MD_GPU_SHADER_STATIC_ASSERT\n"
"#if defined(__cplusplus)\n"
"#define MD_GPU_SHADER_STATIC_ASSERT static_assert\n"
"#else\n"
"#define MD_GPU_SHADER_STATIC_ASSERT _Static_assert\n"
"#endif\n"
"#endif\n\n")

string(JSON ENTRY_COUNT ERROR_VARIABLE ENTRY_COUNT_ERROR LENGTH "${REFLECTION_JSON}" entryPoints)
if (NOT ENTRY_COUNT_ERROR AND ENTRY_COUNT GREATER 0)
  math(EXPR ENTRY_LAST "${ENTRY_COUNT} - 1")
  foreach(ENTRY_INDEX RANGE 0 ${ENTRY_LAST})
    string(JSON ENTRY_STAGE ERROR_VARIABLE ENTRY_STAGE_ERROR GET "${REFLECTION_JSON}" entryPoints ${ENTRY_INDEX} stage)
    if (NOT ENTRY_STAGE_ERROR AND ENTRY_STAGE STREQUAL "compute")
      string(JSON TG_X ERROR_VARIABLE TG_X_ERROR GET "${REFLECTION_JSON}" entryPoints ${ENTRY_INDEX} threadGroupSize 0)
      string(JSON TG_Y ERROR_VARIABLE TG_Y_ERROR GET "${REFLECTION_JSON}" entryPoints ${ENTRY_INDEX} threadGroupSize 1)
      string(JSON TG_Z ERROR_VARIABLE TG_Z_ERROR GET "${REFLECTION_JSON}" entryPoints ${ENTRY_INDEX} threadGroupSize 2)
      if (TG_X_ERROR)
        set(TG_X 0)
      endif()
      if (TG_Y_ERROR)
        set(TG_Y 0)
      endif()
      if (TG_Z_ERROR)
        set(TG_Z 0)
      endif()
      file(APPEND "${OUTPUT}"
"enum {\n"
"    ${SYMBOL}_thread_group_size_x = ${TG_X}u,\n"
"    ${SYMBOL}_thread_group_size_y = ${TG_Y}u,\n"
"    ${SYMBOL}_thread_group_size_z = ${TG_Z}u,\n"
"};\n\n")
      break()
    endif()
  endforeach()
endif()

string(JSON PARAM_COUNT ERROR_VARIABLE PARAM_COUNT_ERROR LENGTH "${REFLECTION_JSON}" parameters)
if (PARAM_COUNT_ERROR OR PARAM_COUNT EQUAL 0)
  return()
endif()

math(EXPR PARAM_LAST "${PARAM_COUNT} - 1")
foreach(PARAM_INDEX RANGE 0 ${PARAM_LAST})
  string(JSON PARAM_NAME GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} name)
  slang_reflect_sanitize_identifier(PARAM_SYMBOL "${PARAM_NAME}")

  string(JSON BINDING_KIND ERROR_VARIABLE BINDING_KIND_ERROR GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} binding kind)
  string(JSON BINDING_INDEX ERROR_VARIABLE BINDING_INDEX_ERROR GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} binding index)
  if (BINDING_KIND_ERROR)
    set(BINDING_KIND "unknown")
  endif()
  if (BINDING_INDEX_ERROR)
    set(BINDING_INDEX 0)
  endif()
  slang_reflect_binding_kind(BINDING_KIND_VALUE "${BINDING_KIND}")

  string(JSON TYPE_KIND ERROR_VARIABLE TYPE_KIND_ERROR GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type kind)
  if (TYPE_KIND_ERROR)
    set(TYPE_KIND "unknown")
  endif()
  slang_reflect_type_kind(TYPE_KIND_VALUE "${TYPE_KIND}")

  if (TYPE_KIND STREQUAL "constantBuffer")
    string(JSON ROOT_SIZE ERROR_VARIABLE ROOT_STRUCT_SIZE_ERROR GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementVarLayout binding size)
    if (ROOT_STRUCT_SIZE_ERROR)
      set(ROOT_SIZE 0)
    endif()

    string(JSON FIELD_COUNT ERROR_VARIABLE ROOT_FIELD_COUNT_ERROR LENGTH "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementType fields)
    if (ROOT_FIELD_COUNT_ERROR)
      set(FIELD_COUNT 0)
    endif()

    set(STRUCT_NAME "${SYMBOL}_${PARAM_SYMBOL}_t")
    file(APPEND "${OUTPUT}" "typedef struct ${STRUCT_NAME} {\n")

    set(CURRENT_OFFSET 0)
    set(PADDING_INDEX 0)
    if (FIELD_COUNT GREATER 0)
      math(EXPR FIELD_LAST "${FIELD_COUNT} - 1")
      foreach(FIELD_INDEX RANGE 0 ${FIELD_LAST})
        string(JSON FIELD_NAME GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementType fields ${FIELD_INDEX} name)
        slang_reflect_sanitize_identifier(FIELD_SYMBOL "${FIELD_NAME}")
        string(JSON FIELD_OFFSET ERROR_VARIABLE FIELD_OFFSET_ERROR GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementType fields ${FIELD_INDEX} binding offset)
        string(JSON FIELD_SIZE ERROR_VARIABLE FIELD_SIZE_ERROR GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementType fields ${FIELD_INDEX} binding size)
        if (FIELD_OFFSET_ERROR)
          set(FIELD_OFFSET ${CURRENT_OFFSET})
        endif()
        if (FIELD_SIZE_ERROR)
          set(FIELD_SIZE 0)
        endif()
        if (FIELD_OFFSET GREATER CURRENT_OFFSET)
          math(EXPR PAD_SIZE "${FIELD_OFFSET} - ${CURRENT_OFFSET}")
          file(APPEND "${OUTPUT}" "    uint8_t _padding${PADDING_INDEX}[${PAD_SIZE}];\n")
          math(EXPR PADDING_INDEX "${PADDING_INDEX} + 1")
        endif()
        slang_reflect_field_decl(FIELD_DECL ${PARAM_INDEX} ${FIELD_INDEX} "${FIELD_SYMBOL}" ${FIELD_SIZE})
        file(APPEND "${OUTPUT}" "    ${FIELD_DECL}\n")
        math(EXPR CURRENT_OFFSET "${FIELD_OFFSET} + ${FIELD_SIZE}")
      endforeach()
    endif()
    if (ROOT_SIZE GREATER CURRENT_OFFSET)
      math(EXPR PAD_SIZE "${ROOT_SIZE} - ${CURRENT_OFFSET}")
      file(APPEND "${OUTPUT}" "    uint8_t _padding${PADDING_INDEX}[${PAD_SIZE}];\n")
    endif()
    file(APPEND "${OUTPUT}" "} ${STRUCT_NAME};\n")
    file(APPEND "${OUTPUT}" "MD_GPU_SHADER_STATIC_ASSERT(sizeof(${STRUCT_NAME}) == ${ROOT_SIZE}u, \"${STRUCT_NAME} size mismatch\");\n")

    if (FIELD_COUNT GREATER 0)
      math(EXPR FIELD_LAST "${FIELD_COUNT} - 1")
      foreach(FIELD_INDEX RANGE 0 ${FIELD_LAST})
        string(JSON FIELD_NAME GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementType fields ${FIELD_INDEX} name)
        slang_reflect_sanitize_identifier(FIELD_SYMBOL "${FIELD_NAME}")
        string(JSON FIELD_OFFSET ERROR_VARIABLE FIELD_OFFSET_ERROR GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementType fields ${FIELD_INDEX} binding offset)
        if (FIELD_OFFSET_ERROR)
          set(FIELD_OFFSET 0)
        endif()
        file(APPEND "${OUTPUT}" "MD_GPU_SHADER_STATIC_ASSERT(offsetof(${STRUCT_NAME}, ${FIELD_SYMBOL}) == ${FIELD_OFFSET}u, \"${STRUCT_NAME}.${FIELD_SYMBOL} offset mismatch\");\n")
      endforeach()
    endif()
    file(APPEND "${OUTPUT}" "\n")
  endif()

  file(APPEND "${OUTPUT}"
"enum {\n"
"    ${SYMBOL}_${PARAM_SYMBOL}_binding_kind = ${BINDING_KIND_VALUE}u,\n"
"    ${SYMBOL}_${PARAM_SYMBOL}_binding_index = ${BINDING_INDEX}u,\n"
"    ${SYMBOL}_${PARAM_SYMBOL}_type_kind = ${TYPE_KIND_VALUE}u,\n")

  if (TYPE_KIND STREQUAL "constantBuffer")
    string(JSON ROOT_SIZE ERROR_VARIABLE ROOT_SIZE_ERROR GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementVarLayout binding size)
    if (ROOT_SIZE_ERROR)
      set(ROOT_SIZE 0)
    endif()

    string(JSON FIELD_COUNT ERROR_VARIABLE FIELD_COUNT_ERROR LENGTH "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementType fields)
    if (FIELD_COUNT_ERROR)
      set(FIELD_COUNT 0)
    endif()

    file(APPEND "${OUTPUT}"
  "    ${SYMBOL}_${PARAM_SYMBOL}_field_count = ${FIELD_COUNT}u,\n")

    if (FIELD_COUNT GREATER 0)
      math(EXPR FIELD_LAST "${FIELD_COUNT} - 1")
      foreach(FIELD_INDEX RANGE 0 ${FIELD_LAST})
        string(JSON FIELD_NAME GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementType fields ${FIELD_INDEX} name)
        slang_reflect_sanitize_identifier(FIELD_SYMBOL "${FIELD_NAME}")
        string(JSON FIELD_OFFSET ERROR_VARIABLE FIELD_OFFSET_ERROR GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementType fields ${FIELD_INDEX} binding offset)
        string(JSON FIELD_SIZE ERROR_VARIABLE FIELD_SIZE_ERROR GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementType fields ${FIELD_INDEX} binding size)
        string(JSON FIELD_TYPE_KIND ERROR_VARIABLE FIELD_TYPE_KIND_ERROR GET "${REFLECTION_JSON}" parameters ${PARAM_INDEX} type elementType fields ${FIELD_INDEX} type kind)
        if (FIELD_OFFSET_ERROR)
          set(FIELD_OFFSET 0)
        endif()
        if (FIELD_SIZE_ERROR)
          set(FIELD_SIZE 0)
        endif()
        if (FIELD_TYPE_KIND_ERROR)
          set(FIELD_TYPE_KIND "unknown")
        endif()
        slang_reflect_type_kind(FIELD_TYPE_KIND_VALUE "${FIELD_TYPE_KIND}")
        file(APPEND "${OUTPUT}"
"    ${SYMBOL}_${PARAM_SYMBOL}_${FIELD_SYMBOL}_type_kind = ${FIELD_TYPE_KIND_VALUE}u,\n")
      endforeach()
    endif()
  endif()

  file(APPEND "${OUTPUT}" "};\n\n")
endforeach()