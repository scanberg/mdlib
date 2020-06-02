#version 330 core

layout (std140) uniform ubo {
	mat4 u_mvp_mat;
};
uniform sampler2D u_depth_tex;

layout (location = 0) in vec3 in_aabb_min;
layout (location = 1) in vec3 in_aabb_max;

out uint out_visible;

bool in_frustum(in vec3 aabb_min, in vec3 aabb_max) {
    //@TODO: Implement
    return true;
}

uint get_cull_bits(vec4 hpos) {
	uint bits = 0U;
	bits |= hpos.x < -hpos.w ?  1U : 0u;
	bits |= hpos.x >  hpos.w ?  2U : 0u;
	bits |= hpos.y < -hpos.w ?  4U : 0u;
	bits |= hpos.y >  hpos.w ?  8U : 0u;
	bits |= hpos.z < -hpos.w ? 16U : 0u;
	bits |= hpos.z >  hpos.w ? 32U : 0u;
	bits |= hpos.w <= 0      ? 64U : 0u;
	return bits;
}

void compute_clip_extent(out vec3 clip_min, out vec3 clip_max, out uint clip_bits, in vec3 aabb_min, in vec3 aabb_max) {
    vec4 v = u_mvp_mat * vec4(aabb_min, 1);
    v.xyz /= v.w;
    clip_min = v.xyz;
    clip_max = v.xyz;
    clip_bits = get_cull_bits(v);
    for (int i = 1; i < 8; i++) {
        vec3 ext = aabb_max * vec3(float(i & 1), float((i >> 1) & 1), float((i >> 2) & 1));
        v = u_mvp_mat * vec4(aabb_min + ext, 1);
        v.xyz /= v.w;
        clip_min = min(clip_min, v.xyz);
        clip_max = max(clip_max, v.xyz);
        clip_bits &= get_cull_bits(v);
    }
}

void main() {
    uint visible = 0U;

    if (in_frustum(in_aabb_min, in_aabb_max)) {
        vec3 clip_min;
        vec3 clip_max;
        uint clip_bits;
        compute_clip_extent(clip_min, clip_max, clip_bits, in_aabb_min, in_aabb_max);

        visible = uint(clip_bits == 0U);
        if (visible != 0U) {
	        vec2 uv_min = clip_min.xy * 0.5 + 0.5;
	        vec2 uv_max = clip_max.xy * 0.5 + 0.5;
	        vec2 uv_ext = (clip_max.xy - clip_min.xy);
	        ivec2 tex_size = textureSize(u_depth_tex, 0);
	        float max_size = max(uv_ext.x, uv_ext.y) * float(max(tex_size.x, tex_size.y));
	        float mip_level = ceil(log2(max_size));
	        
	        float depth = 0;
	        float a = textureLod(u_depth_tex, clip_min.xy, mip_level).r;
	        float b = textureLod(u_depth_tex, vec2(clip_max.x, clip_min.y), mip_level).r;
	        float c = textureLod(u_depth_tex, clip_max.xy, mip_level).r;
	        float d = textureLod(u_depth_tex, vec2(clip_min.x, clip_max.y), mip_level).r;
	        depth = max(depth,max(max(max(a,b),c),d));
	        visible = uint(clip_min.z <= depth);
    	}
    }

    out_visible = visible;
}