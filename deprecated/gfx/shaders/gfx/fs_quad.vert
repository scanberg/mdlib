#version 460

layout(location = 0) out vec2 uv;
void main() {
	vec4 pos =  vec4(
		(float( gl_VertexID       & 1)) * 4.0 - 1.0,
		(float((gl_VertexID >> 1) & 1)) * 4.0 - 1.0,
		0, 1.0);
	uv = pos.xy * 0.5 + 0.5;
	gl_Position = pos;
}
