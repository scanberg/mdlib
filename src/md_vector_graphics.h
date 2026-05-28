#pragma once

#include <core/md_array.h>
#include <core/md_str.h>
#include <core/md_vec_math.h>

#include <stdbool.h>
#include <stdint.h>

struct md_allocator_i;

typedef uint32_t md_vg_color_t;

#define MD_VG_COLOR_RGBA(r, g, b, a) ((((uint32_t)(r) & 0xFFu) << 24u) | (((uint32_t)(g) & 0xFFu) << 16u) | (((uint32_t)(b) & 0xFFu) << 8u) | ((uint32_t)(a) & 0xFFu))
#define MD_VG_COLOR_RGB(r, g, b)      MD_VG_COLOR_RGBA((r), (g), (b), 255)

typedef enum md_vg_paint_type_t {
    MD_VG_PAINT_NONE = 0,
    MD_VG_PAINT_COLOR,
    MD_VG_PAINT_LINEAR_GRADIENT,
} md_vg_paint_type_t;

typedef struct md_vg_linear_gradient_t {
    vec2_t p0;
    vec2_t p1;
    md_vg_color_t col0;
    md_vg_color_t col1;
} md_vg_linear_gradient_t;

typedef struct md_vg_paint_t {
    md_vg_paint_type_t type;
    union {
        md_vg_color_t color;
        md_vg_linear_gradient_t linear_gradient;
    } data;
} md_vg_paint_t;

typedef enum md_vg_line_cap_t {
    MD_VG_LINE_CAP_BUTT = 0,
    MD_VG_LINE_CAP_ROUND,
    MD_VG_LINE_CAP_SQUARE,
} md_vg_line_cap_t;

typedef enum md_vg_line_join_t {
    MD_VG_LINE_JOIN_MITER = 0,
    MD_VG_LINE_JOIN_ROUND,
    MD_VG_LINE_JOIN_BEVEL,
} md_vg_line_join_t;

typedef enum md_vg_fill_rule_t {
    MD_VG_FILL_RULE_NON_ZERO = 0,
    MD_VG_FILL_RULE_EVEN_ODD,
} md_vg_fill_rule_t;

typedef struct md_vg_style_t {
    md_vg_paint_t fill;
    md_vg_paint_t stroke;
    float stroke_width;
    float miter_limit;
    md_vg_line_cap_t line_cap;
    md_vg_line_join_t line_join;
    md_vg_fill_rule_t fill_rule;
} md_vg_style_t;

typedef enum md_vg_path_verb_t {
    MD_VG_PATH_MOVE_TO = 0,
    MD_VG_PATH_LINE_TO,
    MD_VG_PATH_QUADRATIC_TO,
    MD_VG_PATH_CUBIC_TO,
    MD_VG_PATH_CLOSE,
} md_vg_path_verb_t;

typedef enum md_vg_cmd_type_t {
    MD_VG_CMD_LINE = 0,
    MD_VG_CMD_RECT,
    MD_VG_CMD_CIRCLE,
    MD_VG_CMD_ELLIPSE,
    MD_VG_CMD_POLYLINE,
    MD_VG_CMD_POLYGON,
    MD_VG_CMD_BEZIER_QUADRATIC,
    MD_VG_CMD_BEZIER_CUBIC,
    MD_VG_CMD_PATH,
    MD_VG_CMD_TEXT,
} md_vg_cmd_type_t;

typedef struct md_vg_cmd_t {
    md_vg_cmd_type_t type;
    md_vg_style_t style;
    union {
        struct {
            vec2_t p0;
            vec2_t p1;
        } line;
        struct {
            vec2_t min;
            vec2_t max;
            float rounding;
        } rect;
        struct {
            vec2_t center;
            float radius;
        } circle;
        struct {
            vec2_t center;
            vec2_t radius;
            float rotation;
        } ellipse;
        struct {
            uint32_t point_offset;
            uint32_t point_count;
            bool closed;
        } polyline;
        struct {
            uint32_t point_offset;
            uint32_t point_count;
        } polygon;
        struct {
            vec2_t p1;
            vec2_t p2;
            vec2_t p3;
        } bezier_quadratic;
        struct {
            vec2_t p1;
            vec2_t p2;
            vec2_t p3;
            vec2_t p4;
        } bezier_cubic;
        struct {
            uint32_t verb_offset;
            uint32_t verb_count;
            uint32_t point_offset;
            uint32_t point_count;
        } path;
        struct {
            vec2_t pos;
            float font_size;
            uint32_t text_offset;
            uint32_t text_len;
        } text;
    } data;
} md_vg_cmd_t;

typedef struct md_vg_scene_t {
    struct md_allocator_i* alloc;
    vec2_t size;
    str_t default_font_family;
    md_array(md_vg_cmd_t) commands;
    md_array(vec2_t) points;
    md_array(md_vg_path_verb_t) verbs;
    md_array(char) strings;
    md_array(vec2_t) path_points;
    md_array(md_vg_path_verb_t) path_verbs;
} md_vg_scene_t;

#ifdef __cplusplus
extern "C" {
#endif

static inline md_vg_paint_t md_vg_paint_none(void) {
    md_vg_paint_t paint;
    MEMSET(&paint, 0, sizeof(paint));
        paint.type = MD_VG_PAINT_NONE;
    return paint;
}

static inline md_vg_paint_t md_vg_paint_color(md_vg_color_t color) {
    md_vg_paint_t paint;
    MEMSET(&paint, 0, sizeof(paint));
    paint.type = MD_VG_PAINT_COLOR;
    paint.data.color = color;
    return paint;
}

static inline md_vg_paint_t md_vg_paint_linear_gradient(vec2_t p0, vec2_t p1, md_vg_color_t col0, md_vg_color_t col1) {
    md_vg_paint_t paint;
    MEMSET(&paint, 0, sizeof(paint));
    paint.type = MD_VG_PAINT_LINEAR_GRADIENT;
    paint.data.linear_gradient.p0 = p0;
    paint.data.linear_gradient.p1 = p1;
    paint.data.linear_gradient.col0 = col0;
    paint.data.linear_gradient.col1 = col1;
    return paint;
}

static inline md_vg_style_t md_vg_style_default(void) {
    md_vg_style_t style;
    MEMSET(&style, 0, sizeof(style));
    style.fill.type = MD_VG_PAINT_NONE;
    style.stroke.type = MD_VG_PAINT_NONE;
    style.stroke_width = 1.0f;
    style.miter_limit = 4.0f;
    style.line_cap = MD_VG_LINE_CAP_BUTT;
    style.line_join = MD_VG_LINE_JOIN_MITER;
    style.fill_rule = MD_VG_FILL_RULE_NON_ZERO;
    return style;
}

static inline md_vg_style_t md_vg_style_fill(md_vg_color_t color) {
    md_vg_style_t style = md_vg_style_default();
    style.fill = md_vg_paint_color(color);
    return style;
}

static inline md_vg_style_t md_vg_style_stroke(md_vg_color_t color, float width) {
    md_vg_style_t style = md_vg_style_default();
    style.stroke = md_vg_paint_color(color);
    style.stroke_width = width;
    return style;
}

static inline md_vg_style_t md_vg_style_fill_and_stroke(md_vg_color_t fill, md_vg_color_t stroke, float width) {
    md_vg_style_t style = md_vg_style_default();
    style.fill = md_vg_paint_color(fill);
    style.stroke = md_vg_paint_color(stroke);
    style.stroke_width = width;
    return style;
}

md_vg_color_t md_vg_color_from_vec4(vec4_t color);
vec4_t md_vg_color_to_vec4(md_vg_color_t color);

void md_vg_scene_init(md_vg_scene_t* scene, vec2_t size, struct md_allocator_i* alloc);
void md_vg_scene_free(md_vg_scene_t* scene);
void md_vg_scene_clear(md_vg_scene_t* scene);
void md_vg_scene_set_size(md_vg_scene_t* scene, vec2_t size);
void md_vg_scene_set_default_font_family(md_vg_scene_t* scene, str_t font_family);

str_t md_vg_scene_to_svg(const md_vg_scene_t* scene, struct md_allocator_i* alloc);
bool md_vg_scene_write_svg_file(const md_vg_scene_t* scene, str_t path);

void md_vg_add_line(md_vg_scene_t* scene, vec2_t p0, vec2_t p1, md_vg_color_t color, float thickness);
void md_vg_add_line_styled(md_vg_scene_t* scene, vec2_t p0, vec2_t p1, const md_vg_style_t* style);

void md_vg_add_rect(md_vg_scene_t* scene, vec2_t min, vec2_t max, md_vg_color_t color, float rounding, float thickness);
void md_vg_add_rect_filled(md_vg_scene_t* scene, vec2_t min, vec2_t max, md_vg_color_t color, float rounding);
void md_vg_add_rect_styled(md_vg_scene_t* scene, vec2_t min, vec2_t max, float rounding, const md_vg_style_t* style);

void md_vg_add_circle(md_vg_scene_t* scene, vec2_t center, float radius, md_vg_color_t color, int num_segments, float thickness);
void md_vg_add_circle_filled(md_vg_scene_t* scene, vec2_t center, float radius, md_vg_color_t color, int num_segments);
void md_vg_add_circle_styled(md_vg_scene_t* scene, vec2_t center, float radius, int num_segments, const md_vg_style_t* style);

void md_vg_add_ellipse(md_vg_scene_t* scene, vec2_t center, vec2_t radius, md_vg_color_t color, float rotation, int num_segments, float thickness);
void md_vg_add_ellipse_filled(md_vg_scene_t* scene, vec2_t center, vec2_t radius, md_vg_color_t color, float rotation, int num_segments);
void md_vg_add_ellipse_styled(md_vg_scene_t* scene, vec2_t center, vec2_t radius, float rotation, int num_segments, const md_vg_style_t* style);

void md_vg_add_triangle(md_vg_scene_t* scene, vec2_t p0, vec2_t p1, vec2_t p2, md_vg_color_t color, float thickness);
void md_vg_add_triangle_filled(md_vg_scene_t* scene, vec2_t p0, vec2_t p1, vec2_t p2, md_vg_color_t color);
void md_vg_add_triangle_styled(md_vg_scene_t* scene, vec2_t p0, vec2_t p1, vec2_t p2, const md_vg_style_t* style);

void md_vg_add_polyline(md_vg_scene_t* scene, const vec2_t* points, size_t count, md_vg_color_t color, bool closed, float thickness);
void md_vg_add_polyline_styled(md_vg_scene_t* scene, const vec2_t* points, size_t count, bool closed, const md_vg_style_t* style);
void md_vg_add_convex_poly_filled(md_vg_scene_t* scene, const vec2_t* points, size_t count, md_vg_color_t color);
void md_vg_add_polygon_styled(md_vg_scene_t* scene, const vec2_t* points, size_t count, const md_vg_style_t* style);

void md_vg_add_bezier_quadratic(md_vg_scene_t* scene, vec2_t p0, vec2_t p1, vec2_t p2, md_vg_color_t color, float thickness);
void md_vg_add_bezier_quadratic_styled(md_vg_scene_t* scene, vec2_t p0, vec2_t p1, vec2_t p2, const md_vg_style_t* style);

void md_vg_add_bezier_cubic(md_vg_scene_t* scene, vec2_t p0, vec2_t p1, vec2_t p2, vec2_t p3, md_vg_color_t color, float thickness);
void md_vg_add_bezier_cubic_styled(md_vg_scene_t* scene, vec2_t p0, vec2_t p1, vec2_t p2, vec2_t p3, const md_vg_style_t* style);

void md_vg_add_text(md_vg_scene_t* scene, vec2_t pos, float font_size, md_vg_color_t color, str_t text);
void md_vg_add_text_styled(md_vg_scene_t* scene, vec2_t pos, float font_size, const md_vg_style_t* style, str_t text);

void md_vg_path_clear(md_vg_scene_t* scene);
void md_vg_path_move_to(md_vg_scene_t* scene, vec2_t pos);
void md_vg_path_line_to(md_vg_scene_t* scene, vec2_t pos);
void md_vg_path_bezier_quadratic_curve_to(md_vg_scene_t* scene, vec2_t p1, vec2_t p2);
void md_vg_path_bezier_cubic_curve_to(md_vg_scene_t* scene, vec2_t p1, vec2_t p2, vec2_t p3);
void md_vg_path_close(md_vg_scene_t* scene);
void md_vg_path_stroke(md_vg_scene_t* scene, md_vg_color_t color, float thickness);
void md_vg_path_stroke_styled(md_vg_scene_t* scene, const md_vg_style_t* style);
void md_vg_path_fill(md_vg_scene_t* scene, md_vg_color_t color);
void md_vg_path_fill_styled(md_vg_scene_t* scene, const md_vg_style_t* style);

#ifdef __cplusplus
}
#endif
