void project_aabb(vec3 aabb_min, vec3 aabb_max, mat4 MVP, out vec3 out_aabb_min, out vec3 out_aabb_max) {
	vec3 min_xyz = vec3(POS_INF,POS_INF,POS_INF);
	vec3 max_xyz = vec3(NEG_INF,NEG_INF,NEG_INF);
	for (int i = 0; i < 8; ++i) {
		vec3 v = vec3(
			bool((i & 1)) ? aabb_min.x : aabb_max.x,
			bool((i & 2)) ? aabb_min.y : aabb_max.y,
			bool((i & 4)) ? aabb_min.z : aabb_max.z);
		vec4 p = MVP * vec4(v, 1);
		p.xyz = p.xyz / p.w;
		min_xyz = min(min_xyz, p.xyz);
		max_xyz = max(max_xyz, p.xyz);
	}
	// [-1,1] to [0,1]
	out_aabb_min = min_xyz * 0.5 + 0.5;
	out_aabb_max = max_xyz * 0.5 + 0.5;
}

bool cull_frustum(vec3 uv_aabb_min, vec3 uv_aabb_max) {
	return (0 < uv_aabb_max.x && uv_aabb_min.x < 1) && (0 < uv_aabb_max.y && uv_aabb_min.y < 1) && (0 < uv_aabb_max.z);
}

bool cull_hiz(vec3 uv_aabb_min, vec3 uv_aabb_max, sampler2D hiz_tex, uint hiz_width, uint hiz_height) {
	float width  = (uv_aabb_max.x - uv_aabb_min.x) * hiz_width;
	float height = (uv_aabb_max.y - uv_aabb_min.y) * hiz_height;
	float max_size = max(width, height);
	float level = floor(log2(max_size));

	// Sampler is set up to do max reduction, so this computes the maximum depth of a 2x2 texel quad
	float depth_ref = textureLod(hiz_tex, (uv_aabb_min.xy + uv_aabb_max.xy) * 0.5, level).x;
	float depth_min = uv_aabb_min.z;

	return (depth_min <= depth_ref);
}
