#include "utest.h"

#include <md_vector_graphics.h>

#include <core/md_allocator.h>
#include <core/md_str.h>

#include <string.h>

UTEST(vector_graphics, svg_export) {
    md_allocator_i* alloc = md_get_heap_allocator();

    md_vg_scene_t scene = {0};
    md_vg_scene_init(&scene, (vec2_t){ 320.0f, 180.0f }, alloc);
    md_vg_scene_set_default_font_family(&scene, STR_LIT("Helvetica"));

    md_vg_add_line(&scene, (vec2_t){ 8.0f, 8.0f }, (vec2_t){ 64.0f, 24.0f }, MD_VG_COLOR_RGB(255, 0, 0), 2.0f);
    md_vg_add_rect_filled(&scene, (vec2_t){ 12.0f, 32.0f }, (vec2_t){ 88.0f, 72.0f }, MD_VG_COLOR_RGBA(0, 128, 255, 200), 6.0f);
    md_vg_add_circle(&scene, (vec2_t){ 128.0f, 56.0f }, 18.0f, MD_VG_COLOR_RGB(0, 0, 0), 0, 1.5f);
    md_vg_add_ellipse_filled(&scene, (vec2_t){ 184.0f, 56.0f }, (vec2_t){ 28.0f, 16.0f }, MD_VG_COLOR_RGB(0, 200, 120), 0.2f, 0);
    md_vg_add_triangle_filled(&scene, (vec2_t){ 228.0f, 72.0f }, (vec2_t){ 250.0f, 32.0f }, (vec2_t){ 272.0f, 72.0f }, MD_VG_COLOR_RGB(220, 180, 0));

    md_vg_style_t gradient_style = md_vg_style_stroke(MD_VG_COLOR_RGB(0, 0, 0), 6.0f);
    gradient_style.stroke = md_vg_paint_linear_gradient((vec2_t){ 40.0f, 120.0f }, (vec2_t){ 164.0f, 84.0f }, MD_VG_COLOR_RGB(0, 96, 180), MD_VG_COLOR_RGB(255, 180, 0));
    gradient_style.line_cap = MD_VG_LINE_CAP_ROUND;
    md_vg_add_bezier_cubic_styled(&scene,
        (vec2_t){ 40.0f, 120.0f },
        (vec2_t){ 80.0f, 84.0f },
        (vec2_t){ 124.0f, 156.0f },
        (vec2_t){ 164.0f, 112.0f },
        &gradient_style);

    md_vg_path_move_to(&scene, (vec2_t){ 196.0f, 112.0f });
    md_vg_path_line_to(&scene, (vec2_t){ 236.0f, 112.0f });
    md_vg_path_line_to(&scene, (vec2_t){ 252.0f, 148.0f });
    md_vg_path_close(&scene);
    md_vg_path_fill(&scene, MD_VG_COLOR_RGBA(40, 40, 40, 180));

    md_vg_add_text(&scene, (vec2_t){ 12.0f, 144.0f }, 14.0f, MD_VG_COLOR_RGB(20, 20, 20), STR_LIT("Vector export"));

    str_t svg = md_vg_scene_to_svg(&scene, alloc);

    ASSERT_TRUE(svg.ptr != 0);
    ASSERT_TRUE(strstr(svg.ptr, "<svg") != NULL);
    ASSERT_TRUE(strstr(svg.ptr, "<line") != NULL);
    ASSERT_TRUE(strstr(svg.ptr, "<rect") != NULL);
    ASSERT_TRUE(strstr(svg.ptr, "<circle") != NULL);
    ASSERT_TRUE(strstr(svg.ptr, "<ellipse") != NULL);
    ASSERT_TRUE(strstr(svg.ptr, "<polygon") != NULL);
    ASSERT_TRUE(strstr(svg.ptr, "<path") != NULL);
    ASSERT_TRUE(strstr(svg.ptr, "<text") != NULL);
    ASSERT_TRUE(strstr(svg.ptr, "<linearGradient") != NULL);
    ASSERT_TRUE(strstr(svg.ptr, "font-family=\"Helvetica\"") != NULL);
    ASSERT_TRUE(strstr(svg.ptr, "url(#mdvg-grad-") != NULL);
    ASSERT_TRUE(strstr(svg.ptr, "Vector export") != NULL);

    str_free(svg, alloc);
    md_vg_scene_free(&scene);
}