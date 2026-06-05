#include <md_vector_graphics.h>

#include <core/md_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_str_builder.h>

#include <math.h>
#include <string.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

static md_vg_style_t safe_style(const md_vg_style_t* style) {
    return style ? *style : md_vg_style_default();
}

static bool paint_is_none(md_vg_paint_t paint) {
    return paint.type == MD_VG_PAINT_NONE;
}

static bool paint_is_gradient(md_vg_paint_t paint) {
    return paint.type == MD_VG_PAINT_LINEAR_GRADIENT;
}

static uint32_t push_points(md_vg_scene_t* scene, const vec2_t* points, size_t count) {
    ASSERT(scene);
    ASSERT(scene->alloc);
    const uint32_t offset = (uint32_t)md_array_size(scene->points);
    md_array_push_array(scene->points, points, count, scene->alloc);
    return offset;
}

static uint32_t push_verbs(md_vg_scene_t* scene, const md_vg_path_verb_t* verbs, size_t count) {
    ASSERT(scene);
    ASSERT(scene->alloc);
    const uint32_t offset = (uint32_t)md_array_size(scene->verbs);
    md_array_push_array(scene->verbs, verbs, count, scene->alloc);
    return offset;
}

static uint32_t push_string(md_vg_scene_t* scene, str_t str) {
    ASSERT(scene);
    ASSERT(scene->alloc);
    const uint32_t offset = (uint32_t)md_array_size(scene->strings);
    md_array_push_array(scene->strings, str.ptr, str.len, scene->alloc);
    md_array_push(scene->strings, '\0', scene->alloc);
    return offset;
}

static void push_cmd(md_vg_scene_t* scene, const md_vg_cmd_t* cmd) {
    ASSERT(scene);
    ASSERT(scene->alloc);
    md_array_push(scene->commands, *cmd, scene->alloc);
}

static void emit_xml_escaped(md_strb_t* sb, str_t text) {
    for (size_t i = 0; i < text.len; ++i) {
        switch (text.ptr[i]) {
        case '&':
            md_strb_push_cstr(sb, "&amp;");
            break;
        case '<':
            md_strb_push_cstr(sb, "&lt;");
            break;
        case '>':
            md_strb_push_cstr(sb, "&gt;");
            break;
        case '\"':
            md_strb_push_cstr(sb, "&quot;");
            break;
        case '\'':
            md_strb_push_cstr(sb, "&apos;");
            break;
        default:
            md_strb_push_char(sb, text.ptr[i]);
            break;
        }
    }
}

static void color_unpack(md_vg_color_t color, uint8_t rgba[4]) {
    rgba[0] = (uint8_t)((color >> 24u) & 0xFFu);
    rgba[1] = (uint8_t)((color >> 16u) & 0xFFu);
    rgba[2] = (uint8_t)((color >> 8u) & 0xFFu);
    rgba[3] = (uint8_t)(color & 0xFFu);
}

static void emit_color_attr(md_strb_t* sb, const char* attr_name, md_vg_color_t color) {
    uint8_t rgba[4];
    color_unpack(color, rgba);
    md_strb_fmt(sb, " %s=\"rgb(%u,%u,%u)\"", attr_name, rgba[0], rgba[1], rgba[2]);
    if (rgba[3] != 255) {
        md_strb_fmt(sb, " %s-opacity=\"%.9g\"", attr_name, rgba[3] / 255.0);
    }
}

static void emit_paint_attr(md_strb_t* sb, const char* attr_name, md_vg_paint_t paint, size_t cmd_idx) {
    switch (paint.type) {
    default:
    case MD_VG_PAINT_NONE:
        md_strb_fmt(sb, " %s=\"none\"", attr_name);
        break;
    case MD_VG_PAINT_COLOR:
        emit_color_attr(sb, attr_name, paint.data.color);
        break;
    case MD_VG_PAINT_LINEAR_GRADIENT:
        md_strb_fmt(sb, " %s=\"url(#mdvg-grad-%zu-%s)\"", attr_name, cmd_idx, attr_name);
        break;
    }
}

static const char* line_cap_str(md_vg_line_cap_t cap) {
    switch (cap) {
    default:
    case MD_VG_LINE_CAP_BUTT:   return "butt";
    case MD_VG_LINE_CAP_ROUND:  return "round";
    case MD_VG_LINE_CAP_SQUARE: return "square";
    }
}

static const char* line_join_str(md_vg_line_join_t join) {
    switch (join) {
    default:
    case MD_VG_LINE_JOIN_MITER: return "miter";
    case MD_VG_LINE_JOIN_ROUND: return "round";
    case MD_VG_LINE_JOIN_BEVEL: return "bevel";
    }
}

static const char* fill_rule_str(md_vg_fill_rule_t fill_rule) {
    switch (fill_rule) {
    default:
    case MD_VG_FILL_RULE_NON_ZERO: return "nonzero";
    case MD_VG_FILL_RULE_EVEN_ODD: return "evenodd";
    }
}

static void emit_style(md_strb_t* sb, const md_vg_style_t* style, size_t cmd_idx) {
    emit_paint_attr(sb, "fill", style->fill, cmd_idx);
    emit_paint_attr(sb, "stroke", style->stroke, cmd_idx);

    if (!paint_is_none(style->stroke)) {
        md_strb_fmt(sb, " stroke-width=\"%.9g\"", style->stroke_width);
        md_strb_fmt(sb, " stroke-linecap=\"%s\"", line_cap_str(style->line_cap));
        md_strb_fmt(sb, " stroke-linejoin=\"%s\"", line_join_str(style->line_join));
        md_strb_fmt(sb, " stroke-miterlimit=\"%.9g\"", style->miter_limit);
    }
    if (!paint_is_none(style->fill)) {
        md_strb_fmt(sb, " fill-rule=\"%s\"", fill_rule_str(style->fill_rule));
    }
}

static void emit_gradient_def(md_strb_t* sb, size_t cmd_idx, const char* slot, const md_vg_linear_gradient_t* grad) {
    uint8_t rgba0[4];
    uint8_t rgba1[4];
    color_unpack(grad->col0, rgba0);
    color_unpack(grad->col1, rgba1);

    md_strb_fmt(sb, "    <linearGradient id=\"mdvg-grad-%zu-%s\" gradientUnits=\"userSpaceOnUse\" x1=\"%.9g\" y1=\"%.9g\" x2=\"%.9g\" y2=\"%.9g\">\n",
        cmd_idx, slot,
        grad->p0.x, grad->p0.y,
        grad->p1.x, grad->p1.y);
    md_strb_fmt(sb, "      <stop offset=\"0%%\" stop-color=\"rgb(%u,%u,%u)\"", rgba0[0], rgba0[1], rgba0[2]);
    if (rgba0[3] != 255) {
        md_strb_fmt(sb, " stop-opacity=\"%.9g\"", rgba0[3] / 255.0);
    }
    md_strb_push_cstr(sb, "/>\n");
    md_strb_fmt(sb, "      <stop offset=\"100%%\" stop-color=\"rgb(%u,%u,%u)\"", rgba1[0], rgba1[1], rgba1[2]);
    if (rgba1[3] != 255) {
        md_strb_fmt(sb, " stop-opacity=\"%.9g\"", rgba1[3] / 255.0);
    }
    md_strb_push_cstr(sb, "/>\n");
    md_strb_push_cstr(sb, "    </linearGradient>\n");
}

static void emit_defs(md_strb_t* sb, const md_vg_scene_t* scene) {
    bool has_gradients = false;
    for (size_t i = 0; i < md_array_size(scene->commands); ++i) {
        const md_vg_style_t* style = &scene->commands[i].style;
        has_gradients |= paint_is_gradient(style->fill) || paint_is_gradient(style->stroke);
    }
    if (!has_gradients) return;

    md_strb_push_cstr(sb, "  <defs>\n");
    for (size_t i = 0; i < md_array_size(scene->commands); ++i) {
        const md_vg_style_t* style = &scene->commands[i].style;
        if (paint_is_gradient(style->fill)) {
            emit_gradient_def(sb, i, "fill", &style->fill.data.linear_gradient);
        }
        if (paint_is_gradient(style->stroke)) {
            emit_gradient_def(sb, i, "stroke", &style->stroke.data.linear_gradient);
        }
    }
    md_strb_push_cstr(sb, "  </defs>\n");
}

static void emit_points(md_strb_t* sb, const vec2_t* points, size_t count) {
    for (size_t i = 0; i < count; ++i) {
        if (i) md_strb_push_char(sb, ' ');
        md_strb_fmt(sb, "%.9g,%.9g", points[i].x, points[i].y);
    }
}

static void emit_path_data(md_strb_t* sb, const md_vg_scene_t* scene, const md_vg_cmd_t* cmd) {
    const md_vg_path_verb_t* verbs = scene->verbs + cmd->data.path.verb_offset;
    const vec2_t* points = scene->points + cmd->data.path.point_offset;
    size_t point_idx = 0;

    for (size_t i = 0; i < cmd->data.path.verb_count; ++i) {
        switch (verbs[i]) {
        case MD_VG_PATH_MOVE_TO:
            md_strb_fmt(sb, "M %.9g %.9g ", points[point_idx + 0].x, points[point_idx + 0].y);
            point_idx += 1;
            break;
        case MD_VG_PATH_LINE_TO:
            md_strb_fmt(sb, "L %.9g %.9g ", points[point_idx + 0].x, points[point_idx + 0].y);
            point_idx += 1;
            break;
        case MD_VG_PATH_QUADRATIC_TO:
            md_strb_fmt(sb, "Q %.9g %.9g %.9g %.9g ",
                points[point_idx + 0].x, points[point_idx + 0].y,
                points[point_idx + 1].x, points[point_idx + 1].y);
            point_idx += 2;
            break;
        case MD_VG_PATH_CUBIC_TO:
            md_strb_fmt(sb, "C %.9g %.9g %.9g %.9g %.9g %.9g ",
                points[point_idx + 0].x, points[point_idx + 0].y,
                points[point_idx + 1].x, points[point_idx + 1].y,
                points[point_idx + 2].x, points[point_idx + 2].y);
            point_idx += 3;
            break;
        case MD_VG_PATH_CLOSE:
            md_strb_push_cstr(sb, "Z ");
            break;
        default:
            break;
        }
    }
}

static void emit_svg_scene(md_strb_t* sb, const md_vg_scene_t* scene) {
    const char* font_family = scene->default_font_family.ptr ? scene->default_font_family.ptr : "sans-serif";
    const int font_family_len = scene->default_font_family.ptr ? (int)scene->default_font_family.len : (int)strlen(font_family);

    md_strb_push_cstr(sb, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    md_strb_fmt(sb,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" width=\"%.9g\" height=\"%.9g\" viewBox=\"0 0 %.9g %.9g\">\n",
        scene->size.x, scene->size.y, scene->size.x, scene->size.y);

    emit_defs(sb, scene);

    for (size_t i = 0; i < md_array_size(scene->commands); ++i) {
        const md_vg_cmd_t* cmd = &scene->commands[i];
        switch (cmd->type) {
        case MD_VG_CMD_LINE:
            md_strb_fmt(sb, "  <line x1=\"%.9g\" y1=\"%.9g\" x2=\"%.9g\" y2=\"%.9g\"",
                cmd->data.line.p0.x, cmd->data.line.p0.y,
                cmd->data.line.p1.x, cmd->data.line.p1.y);
            emit_style(sb, &cmd->style, i);
            md_strb_push_cstr(sb, " />\n");
            break;
        case MD_VG_CMD_RECT: {
            const vec2_t size = { cmd->data.rect.max.x - cmd->data.rect.min.x, cmd->data.rect.max.y - cmd->data.rect.min.y };
            md_strb_fmt(sb, "  <rect x=\"%.9g\" y=\"%.9g\" width=\"%.9g\" height=\"%.9g\"",
                cmd->data.rect.min.x, cmd->data.rect.min.y, size.x, size.y);
            if (cmd->data.rect.rounding > 0.0f) {
                md_strb_fmt(sb, " rx=\"%.9g\" ry=\"%.9g\"", cmd->data.rect.rounding, cmd->data.rect.rounding);
            }
            emit_style(sb, &cmd->style, i);
            md_strb_push_cstr(sb, " />\n");
            break;
        }
        case MD_VG_CMD_CIRCLE:
            md_strb_fmt(sb, "  <circle cx=\"%.9g\" cy=\"%.9g\" r=\"%.9g\"",
                cmd->data.circle.center.x, cmd->data.circle.center.y, cmd->data.circle.radius);
            emit_style(sb, &cmd->style, i);
            md_strb_push_cstr(sb, " />\n");
            break;
        case MD_VG_CMD_ELLIPSE:
            md_strb_fmt(sb, "  <ellipse cx=\"%.9g\" cy=\"%.9g\" rx=\"%.9g\" ry=\"%.9g\"",
                cmd->data.ellipse.center.x, cmd->data.ellipse.center.y,
                cmd->data.ellipse.radius.x, cmd->data.ellipse.radius.y);
            if (cmd->data.ellipse.rotation != 0.0f) {
                md_strb_fmt(sb, " transform=\"rotate(%.9g %.9g %.9g)\"",
                    cmd->data.ellipse.rotation * (180.0 / PI),
                    cmd->data.ellipse.center.x,
                    cmd->data.ellipse.center.y);
            }
            emit_style(sb, &cmd->style, i);
            md_strb_push_cstr(sb, " />\n");
            break;
        case MD_VG_CMD_POLYLINE: {
            const vec2_t* points = scene->points + cmd->data.polyline.point_offset;
            if (cmd->data.polyline.closed) {
                md_strb_push_cstr(sb, "  <polygon points=\"");
            } else {
                md_strb_push_cstr(sb, "  <polyline points=\"");
            }
            emit_points(sb, points, cmd->data.polyline.point_count);
            md_strb_push_char(sb, '\"');
            emit_style(sb, &cmd->style, i);
            md_strb_push_cstr(sb, " />\n");
            break;
        }
        case MD_VG_CMD_POLYGON: {
            const vec2_t* points = scene->points + cmd->data.polygon.point_offset;
            md_strb_push_cstr(sb, "  <polygon points=\"");
            emit_points(sb, points, cmd->data.polygon.point_count);
            md_strb_push_char(sb, '\"');
            emit_style(sb, &cmd->style, i);
            md_strb_push_cstr(sb, " />\n");
            break;
        }
        case MD_VG_CMD_BEZIER_QUADRATIC:
            md_strb_fmt(sb, "  <path d=\"M %.9g %.9g Q %.9g %.9g %.9g %.9g\"",
                cmd->data.bezier_quadratic.p1.x, cmd->data.bezier_quadratic.p1.y,
                cmd->data.bezier_quadratic.p2.x, cmd->data.bezier_quadratic.p2.y,
                cmd->data.bezier_quadratic.p3.x, cmd->data.bezier_quadratic.p3.y);
            emit_style(sb, &cmd->style, i);
            md_strb_push_cstr(sb, " />\n");
            break;
        case MD_VG_CMD_BEZIER_CUBIC:
            md_strb_fmt(sb, "  <path d=\"M %.9g %.9g C %.9g %.9g %.9g %.9g %.9g %.9g\"",
                cmd->data.bezier_cubic.p1.x, cmd->data.bezier_cubic.p1.y,
                cmd->data.bezier_cubic.p2.x, cmd->data.bezier_cubic.p2.y,
                cmd->data.bezier_cubic.p3.x, cmd->data.bezier_cubic.p3.y,
                cmd->data.bezier_cubic.p4.x, cmd->data.bezier_cubic.p4.y);
            emit_style(sb, &cmd->style, i);
            md_strb_push_cstr(sb, " />\n");
            break;
        case MD_VG_CMD_PATH:
            md_strb_push_cstr(sb, "  <path d=\"");
            emit_path_data(sb, scene, cmd);
            md_strb_push_char(sb, '\"');
            emit_style(sb, &cmd->style, i);
            md_strb_push_cstr(sb, " />\n");
            break;
        case MD_VG_CMD_TEXT: {
            const str_t text = {
                .ptr = scene->strings + cmd->data.text.text_offset,
                .len = cmd->data.text.text_len,
            };
            md_strb_fmt(sb, "  <text x=\"%.9g\" y=\"%.9g\" font-size=\"%.9g\" font-family=\"%.*s\" dominant-baseline=\"hanging\" xml:space=\"preserve\"",
                cmd->data.text.pos.x,
                cmd->data.text.pos.y,
                cmd->data.text.font_size,
                font_family_len,
                font_family);
            emit_style(sb, &cmd->style, i);
            md_strb_push_cstr(sb, ">");
            emit_xml_escaped(sb, text);
            md_strb_push_cstr(sb, "</text>\n");
            break;
        }
        default:
            break;
        }
    }

    md_strb_push_cstr(sb, "</svg>\n");
}

md_vg_color_t md_vg_color_from_vec4(vec4_t color) {
    uint32_t r = (uint32_t)CLAMP(color.x * 255.0f + 0.5f, 0.0f, 255.0f);
    uint32_t g = (uint32_t)CLAMP(color.y * 255.0f + 0.5f, 0.0f, 255.0f);
    uint32_t b = (uint32_t)CLAMP(color.z * 255.0f + 0.5f, 0.0f, 255.0f);
    uint32_t a = (uint32_t)CLAMP(color.w * 255.0f + 0.5f, 0.0f, 255.0f);
    return MD_VG_COLOR_RGBA(r, g, b, a);
}

vec4_t md_vg_color_to_vec4(md_vg_color_t color) {
    uint8_t rgba[4];
    color_unpack(color, rgba);
    return (vec4_t){
        .x = rgba[0] / 255.0f,
        .y = rgba[1] / 255.0f,
        .z = rgba[2] / 255.0f,
        .w = rgba[3] / 255.0f,
    };
}

void md_vg_scene_init(md_vg_scene_t* scene, vec2_t size, md_allocator_i* alloc) {
    ASSERT(scene);
    ASSERT(alloc);
    *scene = (md_vg_scene_t){0};
    scene->size = size;
    scene->alloc = alloc;
    scene->default_font_family = str_copy(STR_LIT("sans-serif"), alloc);
}

void md_vg_scene_free(md_vg_scene_t* scene) {
    if (!scene || !scene->alloc) return;
    md_array_free(scene->commands, scene->alloc);
    md_array_free(scene->points, scene->alloc);
    md_array_free(scene->verbs, scene->alloc);
    md_array_free(scene->strings, scene->alloc);
    md_array_free(scene->path_points, scene->alloc);
    md_array_free(scene->path_verbs, scene->alloc);
    str_free(scene->default_font_family, scene->alloc);
    *scene = (md_vg_scene_t){0};
}

void md_vg_scene_clear(md_vg_scene_t* scene) {
    ASSERT(scene);
    md_array_shrink(scene->commands, 0);
    md_array_shrink(scene->points, 0);
    md_array_shrink(scene->verbs, 0);
    md_array_shrink(scene->strings, 0);
    md_array_shrink(scene->path_points, 0);
    md_array_shrink(scene->path_verbs, 0);
}

void md_vg_scene_set_size(md_vg_scene_t* scene, vec2_t size) {
    ASSERT(scene);
    scene->size = size;
}

void md_vg_scene_set_default_font_family(md_vg_scene_t* scene, str_t font_family) {
    ASSERT(scene);
    ASSERT(scene->alloc);
    str_free(scene->default_font_family, scene->alloc);
    scene->default_font_family = str_copy(font_family, scene->alloc);
}

str_t md_vg_scene_to_svg(const md_vg_scene_t* scene, md_allocator_i* alloc) {
    ASSERT(scene);
    ASSERT(alloc);

    md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);
    md_strb_t sb = md_strb_create(temp_alloc);
    emit_svg_scene(&sb, scene);
    str_t result = str_copy(md_strb_to_str(sb), alloc);
    md_temp_end(temp_scope);
    return result;
}

bool md_vg_scene_write_svg_file(const md_vg_scene_t* scene, str_t path) {
    ASSERT(scene);

    md_temp_scope_t temp_scope = md_temp_begin_avoid(scene->alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp_scope);
    md_strb_t sb = md_strb_create(temp_alloc);
    emit_svg_scene(&sb, scene);
    str_t svg = md_strb_to_str(sb);

    md_file_t file = {0};
    bool result = false;
    if (md_file_open(&file, path, MD_FILE_WRITE | MD_FILE_CREATE | MD_FILE_TRUNCATE)) {
        result = md_file_write(file, svg.ptr, svg.len) == svg.len;
        md_file_close(&file);
    }

    md_temp_end(temp_scope);
    return result;
}

void md_vg_add_line(md_vg_scene_t* scene, vec2_t p0, vec2_t p1, md_vg_color_t color, float thickness) {
    md_vg_add_line_styled(scene, p0, p1, &(md_vg_style_t){
        .fill = { .type = MD_VG_PAINT_NONE },
        .stroke = { .type = MD_VG_PAINT_COLOR, .data.color = color },
        .stroke_width = thickness,
        .miter_limit = 4.0f,
        .line_cap = MD_VG_LINE_CAP_BUTT,
        .line_join = MD_VG_LINE_JOIN_MITER,
        .fill_rule = MD_VG_FILL_RULE_NON_ZERO,
    });
}

void md_vg_add_line_styled(md_vg_scene_t* scene, vec2_t p0, vec2_t p1, const md_vg_style_t* style) {
    md_vg_cmd_t cmd = {
        .type = MD_VG_CMD_LINE,
        .style = safe_style(style),
        .data.line = { .p0 = p0, .p1 = p1 },
    };
    cmd.style.fill = md_vg_paint_none();
    push_cmd(scene, &cmd);
}

void md_vg_add_rect(md_vg_scene_t* scene, vec2_t min, vec2_t max, md_vg_color_t color, float rounding, float thickness) {
    md_vg_add_rect_styled(scene, min, max, rounding, &(md_vg_style_t){
        .fill = { .type = MD_VG_PAINT_NONE },
        .stroke = { .type = MD_VG_PAINT_COLOR, .data.color = color },
        .stroke_width = thickness,
        .miter_limit = 4.0f,
        .line_cap = MD_VG_LINE_CAP_BUTT,
        .line_join = MD_VG_LINE_JOIN_MITER,
        .fill_rule = MD_VG_FILL_RULE_NON_ZERO,
    });
}

void md_vg_add_rect_filled(md_vg_scene_t* scene, vec2_t min, vec2_t max, md_vg_color_t color, float rounding) {
    md_vg_add_rect_styled(scene, min, max, rounding, &(md_vg_style_t){
        .fill = { .type = MD_VG_PAINT_COLOR, .data.color = color },
        .stroke = { .type = MD_VG_PAINT_NONE },
        .stroke_width = 1.0f,
        .miter_limit = 4.0f,
        .line_cap = MD_VG_LINE_CAP_BUTT,
        .line_join = MD_VG_LINE_JOIN_MITER,
        .fill_rule = MD_VG_FILL_RULE_NON_ZERO,
    });
}

void md_vg_add_rect_styled(md_vg_scene_t* scene, vec2_t min, vec2_t max, float rounding, const md_vg_style_t* style) {
    md_vg_cmd_t cmd = {
        .type = MD_VG_CMD_RECT,
        .style = safe_style(style),
        .data.rect = {
            .min = min,
            .max = max,
            .rounding = rounding,
        },
    };
    push_cmd(scene, &cmd);
}

void md_vg_add_circle(md_vg_scene_t* scene, vec2_t center, float radius, md_vg_color_t color, int num_segments, float thickness) {
    (void)num_segments;
    md_vg_add_circle_styled(scene, center, radius, num_segments, &(md_vg_style_t){
        .fill = { .type = MD_VG_PAINT_NONE },
        .stroke = { .type = MD_VG_PAINT_COLOR, .data.color = color },
        .stroke_width = thickness,
        .miter_limit = 4.0f,
        .line_cap = MD_VG_LINE_CAP_BUTT,
        .line_join = MD_VG_LINE_JOIN_MITER,
        .fill_rule = MD_VG_FILL_RULE_NON_ZERO,
    });
}

void md_vg_add_circle_filled(md_vg_scene_t* scene, vec2_t center, float radius, md_vg_color_t color, int num_segments) {
    (void)num_segments;
    md_vg_add_circle_styled(scene, center, radius, num_segments, &(md_vg_style_t){
        .fill = { .type = MD_VG_PAINT_COLOR, .data.color = color },
        .stroke = { .type = MD_VG_PAINT_NONE },
        .stroke_width = 1.0f,
        .miter_limit = 4.0f,
        .line_cap = MD_VG_LINE_CAP_BUTT,
        .line_join = MD_VG_LINE_JOIN_MITER,
        .fill_rule = MD_VG_FILL_RULE_NON_ZERO,
    });
}

void md_vg_add_circle_styled(md_vg_scene_t* scene, vec2_t center, float radius, int num_segments, const md_vg_style_t* style) {
    (void)num_segments;
    md_vg_cmd_t cmd = {
        .type = MD_VG_CMD_CIRCLE,
        .style = safe_style(style),
        .data.circle = {
            .center = center,
            .radius = radius,
        },
    };
    push_cmd(scene, &cmd);
}

void md_vg_add_ellipse(md_vg_scene_t* scene, vec2_t center, vec2_t radius, md_vg_color_t color, float rotation, int num_segments, float thickness) {
    (void)num_segments;
    md_vg_add_ellipse_styled(scene, center, radius, rotation, num_segments, &(md_vg_style_t){
        .fill = { .type = MD_VG_PAINT_NONE },
        .stroke = { .type = MD_VG_PAINT_COLOR, .data.color = color },
        .stroke_width = thickness,
        .miter_limit = 4.0f,
        .line_cap = MD_VG_LINE_CAP_BUTT,
        .line_join = MD_VG_LINE_JOIN_MITER,
        .fill_rule = MD_VG_FILL_RULE_NON_ZERO,
    });
}

void md_vg_add_ellipse_filled(md_vg_scene_t* scene, vec2_t center, vec2_t radius, md_vg_color_t color, float rotation, int num_segments) {
    (void)num_segments;
    md_vg_add_ellipse_styled(scene, center, radius, rotation, num_segments, &(md_vg_style_t){
        .fill = { .type = MD_VG_PAINT_COLOR, .data.color = color },
        .stroke = { .type = MD_VG_PAINT_NONE },
        .stroke_width = 1.0f,
        .miter_limit = 4.0f,
        .line_cap = MD_VG_LINE_CAP_BUTT,
        .line_join = MD_VG_LINE_JOIN_MITER,
        .fill_rule = MD_VG_FILL_RULE_NON_ZERO,
    });
}

void md_vg_add_ellipse_styled(md_vg_scene_t* scene, vec2_t center, vec2_t radius, float rotation, int num_segments, const md_vg_style_t* style) {
    (void)num_segments;
    md_vg_cmd_t cmd = {
        .type = MD_VG_CMD_ELLIPSE,
        .style = safe_style(style),
        .data.ellipse = {
            .center = center,
            .radius = radius,
            .rotation = rotation,
        },
    };
    push_cmd(scene, &cmd);
}

void md_vg_add_triangle(md_vg_scene_t* scene, vec2_t p0, vec2_t p1, vec2_t p2, md_vg_color_t color, float thickness) {
    const vec2_t points[] = { p0, p1, p2 };
    md_vg_add_polyline(scene, points, ARRAY_SIZE(points), color, true, thickness);
}

void md_vg_add_triangle_filled(md_vg_scene_t* scene, vec2_t p0, vec2_t p1, vec2_t p2, md_vg_color_t color) {
    const vec2_t points[] = { p0, p1, p2 };
    md_vg_add_convex_poly_filled(scene, points, ARRAY_SIZE(points), color);
}

void md_vg_add_triangle_styled(md_vg_scene_t* scene, vec2_t p0, vec2_t p1, vec2_t p2, const md_vg_style_t* style) {
    const vec2_t points[] = { p0, p1, p2 };
    md_vg_add_polygon_styled(scene, points, ARRAY_SIZE(points), style);
}

void md_vg_add_polyline(md_vg_scene_t* scene, const vec2_t* points, size_t count, md_vg_color_t color, bool closed, float thickness) {
    md_vg_add_polyline_styled(scene, points, count, closed, &(md_vg_style_t){
        .fill = { .type = MD_VG_PAINT_NONE },
        .stroke = { .type = MD_VG_PAINT_COLOR, .data.color = color },
        .stroke_width = thickness,
        .miter_limit = 4.0f,
        .line_cap = MD_VG_LINE_CAP_BUTT,
        .line_join = MD_VG_LINE_JOIN_MITER,
        .fill_rule = MD_VG_FILL_RULE_NON_ZERO,
    });
}

void md_vg_add_polyline_styled(md_vg_scene_t* scene, const vec2_t* points, size_t count, bool closed, const md_vg_style_t* style) {
    if (!points || count == 0) return;
    md_vg_cmd_t cmd = {
        .type = MD_VG_CMD_POLYLINE,
        .style = safe_style(style),
        .data.polyline = {
            .point_offset = push_points(scene, points, count),
            .point_count = (uint32_t)count,
            .closed = closed,
        },
    };
    if (!closed) {
        cmd.style.fill = md_vg_paint_none();
    }
    push_cmd(scene, &cmd);
}

void md_vg_add_convex_poly_filled(md_vg_scene_t* scene, const vec2_t* points, size_t count, md_vg_color_t color) {
    md_vg_add_polygon_styled(scene, points, count, &(md_vg_style_t){
        .fill = { .type = MD_VG_PAINT_COLOR, .data.color = color },
        .stroke = { .type = MD_VG_PAINT_NONE },
        .stroke_width = 1.0f,
        .miter_limit = 4.0f,
        .line_cap = MD_VG_LINE_CAP_BUTT,
        .line_join = MD_VG_LINE_JOIN_MITER,
        .fill_rule = MD_VG_FILL_RULE_NON_ZERO,
    });
}

void md_vg_add_polygon_styled(md_vg_scene_t* scene, const vec2_t* points, size_t count, const md_vg_style_t* style) {
    if (!points || count == 0) return;
    md_vg_cmd_t cmd = {
        .type = MD_VG_CMD_POLYGON,
        .style = safe_style(style),
        .data.polygon = {
            .point_offset = push_points(scene, points, count),
            .point_count = (uint32_t)count,
        },
    };
    push_cmd(scene, &cmd);
}

void md_vg_add_bezier_quadratic(md_vg_scene_t* scene, vec2_t p0, vec2_t p1, vec2_t p2, md_vg_color_t color, float thickness) {
    md_vg_add_bezier_quadratic_styled(scene, p0, p1, p2, &(md_vg_style_t){
        .fill = { .type = MD_VG_PAINT_NONE },
        .stroke = { .type = MD_VG_PAINT_COLOR, .data.color = color },
        .stroke_width = thickness,
        .miter_limit = 4.0f,
        .line_cap = MD_VG_LINE_CAP_BUTT,
        .line_join = MD_VG_LINE_JOIN_MITER,
        .fill_rule = MD_VG_FILL_RULE_NON_ZERO,
    });
}

void md_vg_add_bezier_quadratic_styled(md_vg_scene_t* scene, vec2_t p0, vec2_t p1, vec2_t p2, const md_vg_style_t* style) {
    md_vg_cmd_t cmd = {
        .type = MD_VG_CMD_BEZIER_QUADRATIC,
        .style = safe_style(style),
        .data.bezier_quadratic = {
            .p1 = p0,
            .p2 = p1,
            .p3 = p2,
        },
    };
    cmd.style.fill = md_vg_paint_none();
    push_cmd(scene, &cmd);
}

void md_vg_add_bezier_cubic(md_vg_scene_t* scene, vec2_t p0, vec2_t p1, vec2_t p2, vec2_t p3, md_vg_color_t color, float thickness) {
    md_vg_add_bezier_cubic_styled(scene, p0, p1, p2, p3, &(md_vg_style_t){
        .fill = { .type = MD_VG_PAINT_NONE },
        .stroke = { .type = MD_VG_PAINT_COLOR, .data.color = color },
        .stroke_width = thickness,
        .miter_limit = 4.0f,
        .line_cap = MD_VG_LINE_CAP_BUTT,
        .line_join = MD_VG_LINE_JOIN_MITER,
        .fill_rule = MD_VG_FILL_RULE_NON_ZERO,
    });
}

void md_vg_add_bezier_cubic_styled(md_vg_scene_t* scene, vec2_t p0, vec2_t p1, vec2_t p2, vec2_t p3, const md_vg_style_t* style) {
    md_vg_cmd_t cmd = {
        .type = MD_VG_CMD_BEZIER_CUBIC,
        .style = safe_style(style),
        .data.bezier_cubic = {
            .p1 = p0,
            .p2 = p1,
            .p3 = p2,
            .p4 = p3,
        },
    };
    cmd.style.fill = md_vg_paint_none();
    push_cmd(scene, &cmd);
}

void md_vg_add_text(md_vg_scene_t* scene, vec2_t pos, float font_size, md_vg_color_t color, str_t text) {
    md_vg_add_text_styled(scene, pos, font_size, &(md_vg_style_t){
        .fill = { .type = MD_VG_PAINT_COLOR, .data.color = color },
        .stroke = { .type = MD_VG_PAINT_NONE },
        .stroke_width = 1.0f,
        .miter_limit = 4.0f,
        .line_cap = MD_VG_LINE_CAP_BUTT,
        .line_join = MD_VG_LINE_JOIN_MITER,
        .fill_rule = MD_VG_FILL_RULE_NON_ZERO,
    }, text);
}

void md_vg_add_text_styled(md_vg_scene_t* scene, vec2_t pos, float font_size, const md_vg_style_t* style, str_t text) {
    if (!text.ptr || text.len == 0) return;
    md_vg_cmd_t cmd = {
        .type = MD_VG_CMD_TEXT,
        .style = safe_style(style),
        .data.text = {
            .pos = pos,
            .font_size = font_size,
            .text_offset = push_string(scene, text),
            .text_len = (uint32_t)text.len,
        },
    };
    push_cmd(scene, &cmd);
}

void md_vg_path_clear(md_vg_scene_t* scene) {
    ASSERT(scene);
    md_array_shrink(scene->path_points, 0);
    md_array_shrink(scene->path_verbs, 0);
}

void md_vg_path_move_to(md_vg_scene_t* scene, vec2_t pos) {
    ASSERT(scene);
    md_array_push(scene->path_verbs, MD_VG_PATH_MOVE_TO, scene->alloc);
    md_array_push(scene->path_points, pos, scene->alloc);
}

void md_vg_path_line_to(md_vg_scene_t* scene, vec2_t pos) {
    ASSERT(scene);
    if (md_array_size(scene->path_verbs) == 0) {
        md_vg_path_move_to(scene, pos);
        return;
    }
    md_array_push(scene->path_verbs, MD_VG_PATH_LINE_TO, scene->alloc);
    md_array_push(scene->path_points, pos, scene->alloc);
}

void md_vg_path_bezier_quadratic_curve_to(md_vg_scene_t* scene, vec2_t p1, vec2_t p2) {
    ASSERT(scene);
    md_array_push(scene->path_verbs, MD_VG_PATH_QUADRATIC_TO, scene->alloc);
    md_array_push(scene->path_points, p1, scene->alloc);
    md_array_push(scene->path_points, p2, scene->alloc);
}

void md_vg_path_bezier_cubic_curve_to(md_vg_scene_t* scene, vec2_t p1, vec2_t p2, vec2_t p3) {
    ASSERT(scene);
    md_array_push(scene->path_verbs, MD_VG_PATH_CUBIC_TO, scene->alloc);
    md_array_push(scene->path_points, p1, scene->alloc);
    md_array_push(scene->path_points, p2, scene->alloc);
    md_array_push(scene->path_points, p3, scene->alloc);
}

void md_vg_path_close(md_vg_scene_t* scene) {
    ASSERT(scene);
    if (md_array_size(scene->path_verbs) == 0) return;
    md_array_push(scene->path_verbs, MD_VG_PATH_CLOSE, scene->alloc);
}

static void commit_path(md_vg_scene_t* scene, const md_vg_style_t* style) {
    if (!scene || md_array_size(scene->path_verbs) == 0) return;
    md_vg_cmd_t cmd = {
        .type = MD_VG_CMD_PATH,
        .style = safe_style(style),
        .data.path = {
            .verb_offset = push_verbs(scene, scene->path_verbs, md_array_size(scene->path_verbs)),
            .verb_count = (uint32_t)md_array_size(scene->path_verbs),
            .point_offset = push_points(scene, scene->path_points, md_array_size(scene->path_points)),
            .point_count = (uint32_t)md_array_size(scene->path_points),
        },
    };
    push_cmd(scene, &cmd);
    md_vg_path_clear(scene);
}

void md_vg_path_stroke(md_vg_scene_t* scene, md_vg_color_t color, float thickness) {
    md_vg_path_stroke_styled(scene, &(md_vg_style_t){
        .fill = { .type = MD_VG_PAINT_NONE },
        .stroke = { .type = MD_VG_PAINT_COLOR, .data.color = color },
        .stroke_width = thickness,
        .miter_limit = 4.0f,
        .line_cap = MD_VG_LINE_CAP_BUTT,
        .line_join = MD_VG_LINE_JOIN_MITER,
        .fill_rule = MD_VG_FILL_RULE_NON_ZERO,
    });
}

void md_vg_path_stroke_styled(md_vg_scene_t* scene, const md_vg_style_t* style) {
    md_vg_style_t cmd_style = safe_style(style);
    cmd_style.fill = md_vg_paint_none();
    commit_path(scene, &cmd_style);
}

void md_vg_path_fill(md_vg_scene_t* scene, md_vg_color_t color) {
    md_vg_path_fill_styled(scene, &(md_vg_style_t){
        .fill = { .type = MD_VG_PAINT_COLOR, .data.color = color },
        .stroke = { .type = MD_VG_PAINT_NONE },
        .stroke_width = 1.0f,
        .miter_limit = 4.0f,
        .line_cap = MD_VG_LINE_CAP_BUTT,
        .line_join = MD_VG_LINE_JOIN_MITER,
        .fill_rule = MD_VG_FILL_RULE_NON_ZERO,
    });
}

void md_vg_path_fill_styled(md_vg_scene_t* scene, const md_vg_style_t* style) {
    commit_path(scene, style);
}
