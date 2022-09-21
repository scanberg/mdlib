#define REP_TYPE_SPACEFILL	0
#define REP_TYPE_LICORICE	1
#define REP_TYPE_VDW		2
#define REP_TYPE_CARTOON	3
#define REP_TYPE_RIBBONS	4
#define REP_TYPE_GAUSSIAN	5
#define REP_TYPE_SES		6

#define POSITION_BINDING  0
#define RADIUS_BINDING    1
#define GROUP_BINDING     2
#define COLOR_BINDING     3
#define TRANSFORM_BINDING 4

#define DRAW_OP_BINDING	5
#define DRAW_INDIRECT_BINDING 6
#define DRAW_SPHERE_INDEX_BINDING 7

#define RASTER_SPHERE_GROUP_SIZE 32
#define CULL_GROUP_SIZE 1024
#define DEPTH_REDUCE_GROUP_SIZE 32

#define INF (1.0/0.0)

struct UniformData {
	mat4 world_to_view;
	mat4 world_to_view_normal;
	mat4 view_to_clip;
	mat4 clip_to_view;
	vec4 frustum_plane[4];
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
	uint  instanceCount;
	uint  firstIndex;
	uint  baseVertex;
	uint  baseInstance;
	uint  _pad[3];
};

struct DrawSpheresIndirectCommand {
	DrawArraysIndirectCommand cmd;
	//Payload
	uint transform_idx;
	uint atom_offset;
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
