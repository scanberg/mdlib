#version 460

#include "common.h"

layout(location = 0) in vec2 uv;

layout(binding = DEPTH_TEX_BINDING)  uniform sampler2D in_depth_tex;
layout(binding = COLOR_TEX_BINDING)  uniform sampler2D in_color_tex;
layout(binding = NORMAL_TEX_BINDING) uniform sampler2D in_normal_tex;
layout(binding = SS_VEL_TEX_BINDING) uniform sampler2D in_ss_vel_tex;
layout(binding = INDEX_TEX_BINDING)  uniform sampler2D in_index_tex;

// These are the outputs of the external FBO, this needs to be configured in some way
layout(location = 0) out vec4 out_color;
layout(location = 1) out vec2 out_normal;
layout(location = 2) out vec2 out_velocity;
layout(location = 3) out vec4 out_index;

void main()
{
	gl_FragDepth = texture(in_depth_tex, uv).x;
	out_color    = texture(in_color_tex, uv);
	out_normal   = texture(in_normal_tex, uv).xy;
	out_velocity = texture(in_ss_vel_tex, uv).xy;
	out_index	 = texture(in_index_tex, uv);
}