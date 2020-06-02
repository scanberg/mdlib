#version 330 core

uniform int       u_depth_lod;
uniform int       u_even_lod;
uniform sampler2D u_depth_tex;

in vec2 uv;

float sample_depth(ivec2 coord, ivec2 tex_size) {
	return texelFetch(u_depth_tex, clamp(coord, ivec2(0), tex_size - 1), u_depth_lod).r;
}

void main() {
    ivec2 lod_size = textureSize(u_depth_tex, u_depth_lod);
    float depth = 0;
  
    if (u_even_lod == 1) {
        ivec2 offsets[] = ivec2[](
			ivec2(0,0),
			ivec2(0,1),
			ivec2(1,1),
			ivec2(1,0)
    	);
	    ivec2 coord = ivec2(gl_FragCoord.xy);
	    coord *= 2;
	    
	    for (int i = 0; i < 4; i++){
	    	depth = max(depth, sample_depth(coord + offsets[i], lod_size));
	    }
	}
  	else {
	    // need this to handle non-power of two
	    // very conservative
	    
	    vec2 offsets[] = vec2[](
	    	vec2(-1,-1),
	      	vec2( 0,-1),
	      	vec2( 1,-1),
	      	vec2(-1, 0),
	      	vec2( 0, 0),
	      	vec2( 1, 0),
	      	vec2(-1, 1),
	      	vec2( 0, 1),
	      	vec2( 1, 1)
	    );
	    vec2 coord = uv;
	    vec2 texel = 1.0 / (vec2(lod_size));
	    
	    for (int i = 0; i < 9; i++){
	    	vec2 pos = coord + offsets[i] * texel;
	    	ivec2 coord = ivec2(pos * lod_size);
			depth = max(depth, sample_depth(coord, lod_size));
	    }
	}

	gl_FragDepth = depth;
}