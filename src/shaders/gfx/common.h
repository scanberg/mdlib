#define REP_TYPE_SPACEFILL	0
#define REP_TYPE_LICORICE	1
#define REP_TYPE_VDW		2
#define REP_TYPE_CARTOON	3
#define REP_TYPE_RIBBONS	4
#define REP_TYPE_GAUSSIAN	5
#define REP_TYPE_SES		6

// Buffer bindings
#define POSITION_BINDING			0
#define RADIUS_BINDING				1
#define CLUSTER_RANGE_BINDING		2
#define CLUSTER_BOUNDS_BINDING		3
#define CLUSTER_POS_RAD_BINDING		4
#define CLUSTER_ATOM_INDEX_BINDING	5
#define COLOR_BINDING				6
#define TRANSFORM_BINDING			7

#define GPU_DRAW_OP_BINDING			9
#define DRAW_INDIRECT_PARAM_BINDING	10
#define DRAW_INDIRECT_CMD_BINDING	11
#define DRAW_SPHERE_INDEX_BINDING	12

#define CLUSTER_INST_DATA_BINDING		14
#define CLUSTER_INST_VIS_IDX_BINDING	15
#define CLUSTER_INST_OCC_IDX_BINDING	16
#define CLUSTER_WRITE_ELEM_BINDING		17	// Bind an array of indices to cluster instances to generate element indices
#define CLUSTER_INST_RAST_IDX_BINDING	18

#define INSTANCE_DATA_BINDING		20
#define INSTANCE_VIS_IDX_BINDING	21
#define INSTANCE_OCC_IDX_BINDING	22

#define DEBUG_BINDING				25
#define VIS_BINDING					26

// Texture bind slots
#define DEPTH_TEX_BINDING	0
#define COLOR_TEX_BINDING	1
#define NORMAL_TEX_BINDING	2
#define SS_VEL_TEX_BINDING	3
#define INDEX_TEX_BINDING	4

// Output locations in drawbuffer
#define COLOR_TEX_ATTACHMENT	0
#define NORMAL_TEX_ATTACHMENT	1
#define SS_VEL_TEX_ATTACHMENT	2
#define INDEX_TEX_ATTACHMENT	3

#define INSTANCE_CULL_GROUP_SIZE 1024
#define CLUSTER_CULL_GROUP_SIZE 128
#define CLUSTER_CULL_LATE_GROUP_SIZE 1024
#define DEPTH_REDUCE_GROUP_SIZE 32
#define CLUSTER_RASTERIZE_GROUP_SIZE 128
#define CLUSTER_COMPUTE_DATA_GROUP_SIZE 128
#define CLUSTER_WRITE_ELEM_IDX_SIZE 128

#define CLUSTER_RASTER_LIMIT_PIXEL_EXTENT 50

#define POS_INF uintBitsToFloat(0x7F800000)
#define NEG_INF uintBitsToFloat(0xFF800000)

// Round up division
#ifndef DIV_UP
#define DIV_UP(x, y) ((x + (y-1)) / y)
#endif

// Rounds up x to nearest multiple of y. No affiliation to the weed-killer
#ifndef ROUND_UP
#define ROUND_UP(x, y) (y * DIV_UP(x,y))
#endif

struct UniformData {
	// Single transforms
	mat4 world_to_view;
	mat4 world_to_view_normal;
	mat4 view_to_clip;
	mat4 clip_to_view;
	// Composed transforms
	mat4 world_to_clip;
	mat4 curr_clip_to_prev_clip;	// Current frame clip to previous frame clip
	mat4 curr_view_to_prev_clip;	// Current frame view to previous frame clip
	mat4 prev_world_to_clip;		// Previous frames world to clip
	mat4 prev_world_to_view;
	mat4 prev_view_to_clip;
	// Misc
	vec4 frustum_plane[4];
	uint frustum_culling;
	uint depth_culling;
	uint depth_pyramid_width;
	uint depth_pyramid_height;
	uint render_width;
	uint render_height;
	float znear;
	float zfar;
	uint late;
	uint _pad[3];
};

struct DebugData {
	uint instance_occ_count;
	uint instance_late_count;
	uint instance_clust_late_count;
	float pick_depth;
	uint  pick_index;
};

// For now this is only a spacefill version for now
struct InstanceData {
	vec4 min_xyzr;
	vec4 max_xyzr;
	uint cluster_offset;
	uint cluster_count;
	uint transform_idx;
	uint _pad;
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

struct DispatchIndirectCommand {
	uint num_groups_x;
	uint num_groups_y;
	uint num_groups_z;
	uint _pad;
};

struct DrawParameters {
	DispatchIndirectCommand write_vis_elem_cmd;
	DispatchIndirectCommand clust_raster_cmd;
	DrawElementsIndirectCommand draw_sphere_cmd;

	DispatchIndirectCommand inst_clust_cull_cmd;
	DispatchIndirectCommand inst_cull_late_cmd;
	DispatchIndirectCommand clust_cull_late_cmd;

	uint cluster_inst_count;
	uint cluster_vis_count;
	uint cluster_occ_count;
	uint cluster_rast_count;

	uint instance_count;
	uint instance_vis_count;
	uint instance_occ_count;
};

// xyz + radius min and max
struct ClusterBounds {
	vec4 cluster_min;
	vec4 cluster_max;
};

struct CompressedPosRad {
	uint xy;
	uint zr;
};

struct ClusterInstance {
	uint  cluster_range;
	uint  cluster_idx;
	uint  transform_idx;
	float rep_scale;
};
