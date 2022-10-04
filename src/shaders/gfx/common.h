#define REP_TYPE_SPACEFILL	0
#define REP_TYPE_LICORICE	1
#define REP_TYPE_VDW		2
#define REP_TYPE_CARTOON	3
#define REP_TYPE_RIBBONS	4
#define REP_TYPE_GAUSSIAN	5
#define REP_TYPE_SES		6

// Buffer bindings
#define POSITION_BINDING  0
#define RADIUS_BINDING    1
#define GROUP_BINDING     2
#define COLOR_BINDING     3
#define TRANSFORM_BINDING 4
#define DRAW_OP_BINDING	5
#define DRAW_INDIRECT_BINDING 6
#define DRAW_SPHERE_INDEX_BINDING 7

#define DEBUG_BINDING 8

// Texture bind slots
#define DEPTH_TEX_BINDING     0
#define COLOR_TEX_BINDING     1
#define NORMAL_TEX_BINDING    2
#define SS_VEL_TEX_BINDING    3
#define INDEX_TEX_BINDING     4

// Output locations in drawbuffer
#define COLOR_TEX_ATTACHMENT   0
#define NORMAL_TEX_ATTACHMENT  1
#define SS_VEL_TEX_ATTACHMENT  2
#define INDEX_TEX_ATTACHMENT   3

#define RASTER_SPHERE_GROUP_SIZE 32
#define CULL_GROUP_SIZE 1024
#define DEPTH_REDUCE_GROUP_SIZE 32

#define POS_INF uintBitsToFloat(0x7F800000)
#define NEG_INF uintBitsToFloat(0xFF800000)

struct UniformData {
	// Single transforms
	mat4 world_to_view;
	mat4 world_to_view_normal;
	mat4 view_to_clip;
	mat4 clip_to_view;
	// Composed transforms
	mat4 world_to_clip;
	mat4 curr_clip_to_prev_clip;	// Current frame clip to previous frame clip
	mat4 prev_world_to_clip;		// Previous frames world to clip
	// Misc
	vec4 frustum_plane[4];
	uint frustum_culling;
	uint depth_culling;
	uint depth_pyramid_width;
	uint depth_pyramid_height;
	float znear;
	float zfar;
	float _pad[2];
};

struct DebugData {
	uint total_groups_processed;
	uint total_groups_drawn;
	float aabb_depth;
	float read_depth;
	float tex_min[2];
	float tex_max[2];
	float aabb_min[2];
	float aabb_max[2];
	float aabb_width;
	float aabb_height;
	float lod;
	float pick_depth;
	uint  pick_index;
	uint _pad[3];
};

struct DrawOp {
	uint  group_offset;
	uint  group_count;
	uint  transform_idx;
	uint  atom_offset;
	uint  color_offset;
	uint  rep_type;
	float rep_args[2];
};

struct Group {
	uint atom_offset;
	uint atom_count;
};

struct Position {
	float x, y, z;
};

struct DrawArraysIndirectCommand {
	uint  count;
	uint  instance_count;
	uint  first;
	uint  base_instance;
};

struct DrawElementsIndirectCommand {
	uint  count;
	uint  instance_count;
	uint  first_index;
	uint  base_vertex;
	uint  base_instance;
};

struct DrawSpheresIndirectCommand {
	//DrawArraysIndirectCommand cmd;
	DrawElementsIndirectCommand cmd;
	//Payload
	uint transform_idx;
	uint color_offset;
	float radius_scale;
};

struct DispatchIndirectCommand {
	uint  num_groups_x;
	uint  num_groups_y;
	uint  num_groups_z;
	uint  _pad;
};

struct DrawIndirect {
	uint draw_sphere_cmd_count;
	uint raster_sphere_cmd_count;
	uint sphere_idx_count;
	uint _pad;
	DispatchIndirectCommand    raster_sphere_cmd;
	DrawSpheresIndirectCommand draw_sphere_cmd[128];
};
