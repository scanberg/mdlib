#pragma once

#include <core/md_str.h>

typedef struct {
	str_t src;
	str_t marker;
} md_gl_shader_src_injection_t;

#ifdef __cplusplus
extern "C" {
#endif

// Add markers for profiling
void md_gl_debug_push(const char* lbl);
void md_gl_debug_pop(void);

bool md_gl_shader_compile(uint32_t shader, str_t src, const md_gl_shader_src_injection_t injections[], size_t num_injections);

bool md_gl_program_attach_and_link(uint32_t program, const uint32_t shaders[], size_t shader_count);
bool md_gl_program_attach_and_link_transform_feedback(uint32_t program, const uint32_t shader[], size_t shader_count, const char* varying[], size_t varying_count, uint32_t capture_mode);

#ifdef __cplusplus
}
#endif