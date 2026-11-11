#include "md_gl_util.h"
#include <core/md_gl_util.h>
#include <core/md_parse.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_writer.h>
#include <core/md_log.h>

#include <GL/gl3w.h>

void md_gl_debug_push(const char* lbl) {
    if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl);
}

void md_gl_debug_pop(void) {
    if (glPopDebugGroup) glPopDebugGroup();
}

bool md_gl_shader_compile(uint32_t shader, md_str_t src, const md_gl_shader_src_injection_t injections[], size_t num_injections) {
    md_temp_scope_t temp = md_temp_begin();
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    md_str_writer_t sb = md_str_writer_create(temp_alloc);

    bool result = false;

    if (!md_str_begins_with(src, MD_STR_LIT("#version "))) {
        MD_LOG_ERROR("Missing version string as first line in shader!");
        goto done;
    }

    if (injections && num_injections > 0) {
        // Insert potential empty markers, they go after #version line
        md_str_t version_line;
        md_parse_next_line(&version_line, &src);
        md_str_writer_push_str(&sb, version_line);
        md_str_writer_push_char(&sb, '\n');
        for (size_t i = 0; i < num_injections; ++i) {
            if (md_str_empty(injections[i].marker)) {
                md_str_writer_push_str(&sb, injections[i].src);
                md_str_writer_push_char(&sb, '\n');
            }
        }

        md_str_t line;
        while (md_parse_next_line(&line, &src)) {
            if (md_str_begins_with(line, MD_STR_LIT("#pragma"))) {
                md_str_t marker = md_str_trim(md_str_substr(line, sizeof("#pragma"), SIZE_MAX));
                for (size_t i = 0; i < num_injections; ++i) {
                    if (md_str_eq(injections[i].marker, marker)) {
                        md_str_writer_push_str(&sb, injections[i].src);
                        md_str_writer_push_char(&sb, '\n');
                    }
                }
            } else {
                md_str_writer_push_str(&sb, line);
                md_str_writer_push_char(&sb, '\n');
            }
        }
    } else {
        md_str_writer_push_str(&sb, src);
    }

    const char* csrc = md_str_writer_to_cstr(&sb);

    glShaderSource(shader, 1, &csrc, 0);
    glCompileShader(shader);

    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);

    if (!success) {
        char err_buf[1024];
        glGetShaderInfoLog(shader, MD_COUNTOF(err_buf), NULL, err_buf);
        MD_LOG_ERROR("Shader compile error:\n%s\n", err_buf);
        goto done;
    }

    result = true;
done:
    md_temp_end(temp);
    return result;
}

bool md_gl_program_attach_and_link(uint32_t program, const uint32_t shaders[], size_t shader_count) {
    MD_ASSERT(program);

    for (size_t i = 0; i < shader_count; i++) {
        GLint compile_status;
        glGetShaderiv(shaders[i], GL_COMPILE_STATUS, &compile_status);
        if (glIsShader(shaders[i]) && compile_status) {
            glAttachShader(program, shaders[i]);
        } else {
            MD_LOG_ERROR("Program link error: One or more shaders are invalid\n");
            return false;
        }
    }

    GLint success;
    glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (!success) {
        char err_buf[256];
        glGetProgramInfoLog(program, MD_COUNTOF(err_buf), NULL, err_buf);
        MD_LOG_ERROR("Program link error:\n%s\n", err_buf);
        return false;
    }

    for (size_t i = 0; i < shader_count; i++) {
        glDetachShader(program, shaders[i]);
    }

    return true;
}

bool md_gl_program_attach_and_link_transform_feedback(uint32_t program, const uint32_t shader[], size_t shader_count, const char* varying[], size_t varying_count, uint32_t capture_mode) {
    MD_ASSERT(program);
    MD_ASSERT(capture_mode == GL_INTERLEAVED_ATTRIBS || capture_mode == GL_SEPARATE_ATTRIBS);

    for (size_t i = 0; i < shader_count; i++) {
        GLint compile_status;
        glGetShaderiv(shader[i], GL_COMPILE_STATUS, &compile_status);
        if (glIsShader(shader[i]) && compile_status) {
            glAttachShader(program, shader[i]);
        } else {
            MD_LOG_ERROR("One or more shaders are invalid");
            return false;
        }
    }

    GLint link_status;
    glTransformFeedbackVaryings(program, (GLsizei)varying_count, varying, capture_mode);
    glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &link_status);
    if (!link_status) {
        char err_buf[256];
        glGetProgramInfoLog(program, MD_COUNTOF(err_buf), NULL, err_buf);
        MD_LOG_ERROR("Program link error:\n%s\n", err_buf);
        return false;
    }

    for (uint32_t i = 0; i < shader_count; i++) {
        glDetachShader(program, shader[i]);
    }

    return true;
}
