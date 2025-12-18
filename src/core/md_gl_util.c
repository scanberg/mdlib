#include "md_gl_util.h"
#include <core/md_gl_util.h>
#include <core/md_str_builder.h>
#include <core/md_log.h>

#include <GL/gl3w.h>

void md_gl_debug_push(const char* lbl) {
    if (glPushDebugGroup) glPushDebugGroup(GL_DEBUG_SOURCE_APPLICATION, GL_KHR_debug, -1, lbl);
}

void md_gl_debug_pop(void) {
    if (glPopDebugGroup) glPopDebugGroup();
}

bool md_gl_shader_compile(uint32_t shader, str_t src, const md_gl_shader_src_injection_t injections[], size_t num_injections) {
    size_t temp_pos = md_temp_get_pos();
    md_strb_t sb = md_strb_create(md_get_temp_allocator());

    bool result = false;

    if (!str_eq_cstr_n(src, "#version ", 9)) {
        MD_LOG_ERROR("Missing version string as first line in shader!");
        goto done;
    }

    if (injections && num_injections > 0) {
        // Insert potential empty markers, they go after #version line
        str_t version_line;
        str_extract_line(&version_line, &src);
        md_strb_push_str(&sb, version_line);
        md_strb_push_char(&sb, '\n');
        for (size_t i = 0; i < num_injections; ++i) {
            if (str_empty(injections[i].marker)) {
                md_strb_push_str (&sb, injections[i].src);
                md_strb_push_char(&sb, '\n');
            }
        }

        str_t line;
        while (str_extract_line(&line, &src)) {
            if (str_begins_with(line, STR_LIT("#pragma"))) {
                str_t marker = str_trim(str_substr(line, sizeof("#pragma"), SIZE_MAX));
                for (size_t i = 0; i < num_injections; ++i) {
                    if (str_eq(injections[i].marker, marker)) {
                        md_strb_push_str (&sb, injections[i].src);
                        md_strb_push_char(&sb, '\n');
                    }
                }
            } else {
                md_strb_push_str(&sb, line);
                md_strb_push_char(&sb, '\n');
            }
        }
    } else {
        md_strb_push_str(&sb, src);
    }

    const char* csrc = md_strb_to_cstr(sb);

    glShaderSource(shader, 1, &csrc, 0);
    glCompileShader(shader);

    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);

    if (!success) {
        char err_buf[1024];
        glGetShaderInfoLog(shader, ARRAY_SIZE(err_buf), NULL, err_buf);
        MD_LOG_ERROR("Shader compile error:\n%s\n", err_buf);
        goto done;
    }

    result = true;
done:
    md_temp_set_pos_back(temp_pos);
    return result;
}

bool md_gl_program_attach_and_link(uint32_t program, const uint32_t shaders[], size_t shader_count) {
    ASSERT(program);

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
        glGetProgramInfoLog(program, ARRAY_SIZE(err_buf), NULL, err_buf);
        MD_LOG_ERROR("Program link error:\n%s\n", err_buf);
        return false;
    }

    for (size_t i = 0; i < shader_count; i++) {
        glDetachShader(program, shaders[i]);
    }

    return true;
}

bool md_gl_program_attach_and_link_transform_feedback(uint32_t program, const uint32_t shader[], size_t shader_count, const char* varying[], size_t varying_count, uint32_t capture_mode) {
    ASSERT(program);
    ASSERT(capture_mode == GL_INTERLEAVED_ATTRIBS || capture_mode == GL_SEPARATE_ATTRIBS);

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
        glGetProgramInfoLog(program, ARRAY_SIZE(err_buf), NULL, err_buf);
        MD_LOG_ERROR("Program link error:\n%s\n", err_buf);
        return false;
    }

    for (uint32_t i = 0; i < shader_count; i++) {
        glDetachShader(program, shader[i]);
    }

    return true;
}
