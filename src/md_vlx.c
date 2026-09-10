#include <md_vlx.h>

#include <core/md_os.h>
#include <core/md_log.h>
#include <core/md_parse.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str_builder.h>
#include <core/md_hash.h>

#include <md_types.h>
#include <md_system.h>
#include <md_gto.h>
#include <md_util.h>

#include <hdf5.h>

// ---------------------------------------------------------------------------
// Reader-internal types.
//
// These used to live in md_vlx.h. Nothing outside this file may name them: what a VeloxChem file
// carries reaches a consumer as system attributes (see the table in md_vlx.h), not as a struct or
// an enum declared here. Keeping them file local is what makes that true rather than aspirational.
// ---------------------------------------------------------------------------

// XPS: one entry per computed core-hole state.
// VeloxChem performs the calculation per element, removing one core electron from each atom of that
// element in turn and taking the total energy difference against the ground state (delta-SCF).
// Field order is deliberate: 'ionization_energy' and 'contribution' are both double and adjacent so
// a plotting consumer can take them as x/y with stride = sizeof(vlx_xps_entry_t) - which is exactly
// what vlx_publish_column does when it publishes them as two contiguous attributes.
typedef struct vlx_xps_entry_t {
	double       ionization_energy;	// unit: eV
	double       contribution;		// Fraction of the core MO localized on 'atom_index'
	int32_t      atom_index;		// Index into the molecular structure
	int32_t      mo_index;			// Index into the MO (orbital) arrays
	md_element_t element;			// Atomic number, matches the owning group
	bool         is_delocalized;	// Core hole is spread over several symmetry equivalent atoms
} vlx_xps_entry_t;


typedef enum {
	VLX_SPIN_ALPHA = 0,
	VLX_SPIN_BETA  = 1,
} vlx_spin_t;

typedef enum {
	VLX_NTO_PARTICLE = 0,
	VLX_NTO_HOLE = 1,
} vlx_nto_type_t;

typedef enum {
	VLX_TRANSITION_ATTACHMENT = 0,
	VLX_TRANSITION_DETACHMENT = 1,
	VLX_TRANSITION_DIFFERENCE = 2,
} vlx_transition_type_t;

typedef enum {
	VLX_SCF_UNKNOWN = 0,
	VLX_SCF_RESTRICTED,
	VLX_SCF_RESTRICTED_OPENSHELL,
	VLX_SCF_UNRESTRICTED,
} vlx_scf_type_t;

typedef enum {
	VLX_RSP_UNKNOWN = 0,
	VLX_RSP_LINEAR,			// Linear response (frequency-dependent polarizabilities)
	VLX_RSP_CPP,				// Complex Polarization Propagator
	VLX_RSP_C6,				// Homomolecular C_6 value (in a.u.)
	VLX_RSP_TPA,		        // Two-photon Absorption (TPA) cross-sections
	VLX_RSP_TPA_TRANSITION,  // Two-Photon Absorption transition properties
	VLX_RSP_RIXS,			// Resonant Inelastic X-ray Scattering (RIXS) cross-sections
} vlx_rsp_type_t;

// PES (Potential Energy Surface) operations (Geometry Optimizations)
typedef enum {
	VLX_OPT_UNKNOWN = 0,
	VLX_OPT_GEOMETRY,			// Local minimum optimization (Ground state or excited state)
	VLX_OPT_CONSTRAINED,			// Optimization with geometric constraints
	VLX_OPT_TS,					// First-order saddle point optimization (Transition State)
	VLX_OPT_IRC,					// Reaction path optimization (Intrinsic Reaction Coordinate)
	VLX_OPT_COUNT,
} vlx_opt_type_t;

// A practical upper bound on the natural transition orbital pairs one excited state can carry.
// Generous for any real calculation and small enough that a row sits on the stack; the extract
// truncates at whatever capacity it is given, so this is a convenience and not a limit of the data.
#define VLX_NTO_MAX_LAMBDAS 32



// Internal only, and for the same reason as vlx_atomic_property_t above: what a density
// property looks like between reading it out of the h5 file and publishing it as an attribute.
typedef struct vlx_density_property_t {
	str_t    label;     // Display text, as authored in the file
	str_t    name;      // Dataset name in the h5 file. The identity, and what the attribute path is built from
	uint64_t key;
	size_t 	 dim[2];
	double*  data;      // dim[0] * dim[1] values, row major
} vlx_density_property_t;

#include <float.h>
#include <math.h>
#include <stdlib.h>	// qsort

#define ANGSTROM_TO_BOHR 1.8897261246257702
#define BOHR_TO_ANGSTROM 0.5291772109029999

#define VLX_NTO_POWER_ITERATIONS 256
#define VLX_NTO_EIGENVALUE_EPSILON 1.0e-14
#define VLX_NTO_CONVERGENCE_EPSILON 1.0e-10


/*

COMMENTS (Robin):

This file is meant to cover the VeloxChem file format.
Alot of the functionality for constructing and extracting basis functions and gaussian type orbitals is extracted from the VeloxChem source code.
https://github.com/VeloxChem/VeloxChem which is released under the LGPL-3.0 license.

*/

// Single contracted basis function
typedef struct basis_set_func_t {
	uint8_t  type; // Azimuthal Quantum Number
	uint8_t  param_count;
	uint16_t param_offset;
} basis_set_func_t;

typedef struct basis_set_basis_t {
	uint8_t  max_type;
	uint8_t  basis_func_count;
	uint16_t basis_func_offset;
} basis_set_basis_t;

typedef struct basis_set_t {
	str_t identifier;
	struct {
		size_t count;
		double* exponents;
		double* normalization_coefficients;
	} param;

	struct {
		size_t count;
		basis_set_func_t* data;
	} basis_func;

	// The atom basis entries are implicitly stored in the order of atomic numbers
	// 0 is a NULL entry, 1 = Hydrogen, 2 = Helium etc.
	struct {
		size_t count;
		basis_set_basis_t* data;
	} atom_basis;
} basis_set_t;

// New format

typedef struct vlx_1d_data_t {
	size_t  size;
	double* data;
} vlx_1d_data_t;

typedef struct vlx_2d_data_t {
	size_t  size[2];
	double* data;
} vlx_2d_data_t;

typedef struct vlx_orbital_t {
	vlx_2d_data_t coefficients;
    vlx_2d_data_t density;
	vlx_1d_data_t energy;
	vlx_1d_data_t occupancy;
	size_t homo_idx;
	size_t lumo_idx;
} vlx_orbital_t;


// Self Consistent Field
typedef struct vlx_scf_t {
	vlx_scf_type_t type;

	double energy;
	dvec3_t ground_state_dipole_moment;

	vlx_orbital_t alpha;
	vlx_orbital_t beta;

	vlx_2d_data_t S;
} vlx_scf_t;


typedef struct vlx_rsp_t {
	vlx_rsp_type_t type;

	size_t   number_of_frequencies;
	size_t   num_core;
	size_t   num_valence;
	size_t   num_virtual;

	dvec3_t* electric_transition_dipoles;
	dvec3_t* magnetic_transition_dipoles;
	dvec3_t* velocity_transition_dipoles;


	// Linear only, Should have dimensions [number_of_frequencies][num_occ * num_vir * 2]
	vlx_2d_data_t solution_matrix;
} vlx_rsp_t;



typedef struct vlx_t {
	basis_set_t basis_set;

	str_t  basis_set_ident;
	str_t  dft_func_label;
	str_t  potfile_text;

	size_t number_of_atoms;
	size_t number_of_alpha_electrons;
	size_t number_of_beta_electrons;

	double molecular_charge;
	double nuclear_repulsion_energy;
	size_t spin_multiplicity;

	// Arrays (length = number_of_atoms)
	dvec3_t* atom_coordinates;

	vlx_density_property_t* density_properties; // Optional data, may be NULL, length is number of density properties

	// Data blocks
	vlx_scf_t scf;
	vlx_rsp_t rsp;

	md_element_t* atomic_numbers;
	int* local_to_global_atom_idx; // Maps local atom indices to global system indices for subsystems. NULL if not a subsystem.
	// ao_remap[shell_ao_idx] = vlx_ao_idx
	// Maps from shell order (angl→atom→func→isph) to VeloxChem matrix row order (angl→isph→atom→func).
	// Built once after the basis set is parsed; used to permute C, D, S matrices into shell order.
	int* ao_remap;

	// The system every block is read INTO. A reader's destination is the attribute table, not this
	// struct: what stays here is only what a later step still has to look at - the atom list the
	// basis is built over, and the AO matrices which are permuted and converted before they can be
	// published. Everything else goes straight into 'sys' as it is read.
	struct md_system_t* sys;

	struct md_allocator_i* arena;
} vlx_t;

// ---------------------------------------------------------------------------
// ATTRIBUTE PUBLISHING
//
// Everything a file carries reaches a consumer as an attribute on the system, so these sit ahead of
// the readers rather than after them: a reader's destination IS the table, and the shortest path
// from an HDF5 dataset to it is vlx_publish_h5_* below, which reads into the storage the table just
// reserved. Nothing here names vlx_t - the table is the output, not the reader.
// ---------------------------------------------------------------------------


// Units the QM blocks are stated in. Not in md_unit.h because nothing outside quantum chemistry asks
// for them, and a unit only ever constructed in one place is better constructed there than named
// globally - but they are needed by the readers now that a reader publishes what it reads, so they
// sit here rather than inside one function.
static inline md_unit_t vlx_unit_hartree(void)			{ return md_unit_hartree(); }
static inline md_unit_t vlx_unit_e_bohr(void)			{ return md_unit_elementary_charge_bohr(); }
static inline md_unit_t vlx_unit_bohr_magneton(void)	{ return md_unit_bohr_magneton(); }
static inline md_unit_t vlx_unit_bohr_velocity(void)	{ return md_unit_bohr_velocity(); }
static inline md_unit_t vlx_unit_angstrom(void)			{ return md_unit_angstrom(); }
static inline md_unit_t vlx_unit_wavenumber(void)		{ return md_unit_pow(md_unit_scl(md_unit_meter(), 1.0e-2), -1); }              // cm^-1
static inline md_unit_t vlx_unit_km_per_mol(void)		{ return md_unit_div(md_unit_scl(md_unit_meter(), 1.0e3), md_unit_mole()); }
static inline md_unit_t vlx_unit_amu(void)				{ return md_unit_scl(md_unit_kilogram(), 1.66053906660e-27); }

// The two calculation-kind enums as TEXT, for the attribute table.
//
// Text and not the enum's integer: a number in the table is a contract on THIS header's ordering,
// which nothing else can read and which an inserted enumerator silently breaks. A second QM reader
// naming its own run "rixs" says the same thing without ever having heard of vlx_rsp_type_t,
// which is the whole reason these leave the reader at all. The strings are the enumerator names
// lowercased, so nothing is needed to read the mapping by.
//
// UNKNOWN publishes nothing rather than the word "unknown": an absent path already means "this file
// does not say", and a consumer which has to handle the absent case gains nothing from a second
// spelling of it.
static str_t vlx_rsp_type_str(vlx_rsp_type_t type) {
	switch (type) {
	case VLX_RSP_LINEAR:			return STR_LIT("linear");
	case VLX_RSP_CPP:			return STR_LIT("cpp");
	case VLX_RSP_C6:				return STR_LIT("c6");
	case VLX_RSP_TPA:			return STR_LIT("tpa");
	case VLX_RSP_TPA_TRANSITION:	return STR_LIT("tpa_transition");
	case VLX_RSP_RIXS:			return STR_LIT("rixs");
	case VLX_RSP_UNKNOWN:
	default:						return (str_t){0};
	}
}

static str_t vlx_scf_type_str(vlx_scf_type_t type) {
	switch (type) {
	case VLX_SCF_RESTRICTED:				return STR_LIT("restricted");
	case VLX_SCF_RESTRICTED_OPENSHELL:	return STR_LIT("restricted_openshell");
	case VLX_SCF_UNRESTRICTED:			return STR_LIT("unrestricted");
	case VLX_SCF_UNKNOWN:
	default:								return (str_t){0};
	}
}

static str_t vlx_opt_type_str(vlx_opt_type_t type) {
	switch (type) {
	case VLX_OPT_GEOMETRY:		return STR_LIT("geometry");
	case VLX_OPT_CONSTRAINED:	return STR_LIT("constrained");
	case VLX_OPT_TS:				return STR_LIT("transition_state");
	case VLX_OPT_IRC:			return STR_LIT("irc");
	case VLX_OPT_UNKNOWN:
	case VLX_OPT_COUNT:
	default:						return (str_t){0};
	}
}

// Builds an attribute path from a fixed group prefix and a name taken from the file. A '/' inside
// the name would silently introduce a group level in a namespace where the separator is the only
// structure there is, so it is folded to '_'. Returns an empty str_t when the name does not fit,
// which the callers treat as "skip this one" rather than as a reason to stop publishing.
static str_t vlx_attribute_path(char* buf, size_t cap, str_t group, str_t name) {
	int len = snprintf(buf, cap, STR_FMT "/" STR_FMT, STR_ARG(group), STR_ARG(name));
	if (len <= 0 || (size_t)len >= cap) {
		MD_LOG_ERROR("Attribute path '" STR_FMT "/" STR_FMT "' does not fit in %zu characters", STR_ARG(group), STR_ARG(name), cap - 1);
		return (str_t){0};
	}
	for (int c = (int)group.len + 1; c < len; ++c) {
		if (buf[c] == '/') buf[c] = '_';
	}
	return str_from_cstrn(buf, (size_t)len);
}

// Publishes one attribute under a path this publisher owns, replacing whatever was there - see
// md_attributes_replace on why that is what a producer wants.
static md_attribute_id_t vlx_publish(md_system_t* sys, str_t path, str_t label, md_unit_t unit, md_attribute_format_t format, const void* data, size_t byte_size) {
	return md_attributes_replace(&sys->attributes, &(md_attribute_desc_t){
		.path      = path,
		.format    = format,
		.unit      = unit,
		.label     = label,
		.data      = data,
		.byte_size = byte_size,
	});
}

// The same, for an attribute computed through a provider instead of one copied in. The provider's
// user_data is 'sys' itself (see the transition density providers below), a borrowed pointer that
// needs no bookkeeping and outlives 'vlx'.
static md_attribute_id_t vlx_publish_virtual(md_system_t* sys, str_t path, str_t label, md_unit_t unit, md_attribute_format_t format, const md_attribute_virtual_t* virt) {
	return md_attributes_replace(&sys->attributes, &(md_attribute_desc_t){
		.path   = path,
		.format = format,
		.unit   = unit,
		.label  = label,
		.virt   = virt,
	});
}

// rank 1 {N}, one scalar per element.
static md_attribute_id_t vlx_publish_series(md_system_t* sys, str_t path, str_t label, md_unit_t unit, const double* values, size_t count) {
	if (!values || count == 0) {
		return MD_ATTRIBUTE_INVALID;
	}
	md_attribute_format_t format = {
		.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)count },
	};
	return vlx_publish(sys, path, label, unit, format, values, count * sizeof(double));
}

// Gives an already published attribute a second, format neutral name. Both names then read one
// datum - no copy, and a consumer of either is unaffected when the other appears or goes.
static md_attribute_id_t vlx_alias(md_system_t* sys, md_attribute_id_t target, str_t path) {
	if (target != MD_ATTRIBUTE_INVALID) {
		const md_attribute_t* existing = md_attributes_find(&sys->attributes, path);
		if (existing) {
			md_attributes_remove(&sys->attributes, existing->id);
		}
		return md_attributes_alias(&sys->attributes, target, path, (str_t){0}, (str_t){0});
	}
	return MD_ATTRIBUTE_INVALID;
}

// Publishes beta's copy of a per orbital series, or a SECOND NAME for alpha's when the two share
// storage. The reader shallow copies beta from alpha for anything but an unrestricted calculation,
// so comparing the POINTERS is what tells the cases apart - and it is the only test that gets the
// restricted open shell case right, where the orbitals are shared but the occupations are read
// separately. Switching on vlx->scf.type instead would alias an occupation array that differs.
static md_attribute_id_t vlx_publish_or_alias(md_system_t* sys, md_attribute_id_t alpha_id, str_t path, str_t label, md_unit_t unit,
                                              const double* alpha_values, const double* beta_values, size_t count) {
	if (beta_values && beta_values == alpha_values) {
		return vlx_alias(sys, alpha_id, path);
	}
	return vlx_publish_series(sys, path, label, unit, beta_values, count);
}

// rank 0, a single scalar. The value is copied, so a local is fine.
static void vlx_publish_scalar(md_system_t* sys, str_t path, str_t label, md_unit_t unit, double value) {
	md_attribute_format_t format = {
		.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 0,
	};
	vlx_publish(sys, path, label, unit, format, &value, sizeof(double));
}

// rank 1 {N} of 3 component values. dvec3_t is three contiguous doubles, so the source array is
// already the interleaved layout an attribute stores and this is a straight copy.
static void vlx_publish_vec3_series(md_system_t* sys, str_t path, str_t label, md_unit_t unit, const dvec3_t* values, size_t count) {
	if (!values || count == 0) {
		return;
	}
	md_attribute_format_t format = {
		.type = MD_ATTRIBUTE_TYPE_F64, .components = 3, .rank = 1, .shape = { (uint32_t)count },
	};
	vlx_publish(sys, path, label, unit, format, values, count * 3 * sizeof(double));
}

// A single string is rank 1 {1}, by the same rule that makes a single 3-vector rank 2 {1,3}. The
// descriptor carries the TEXT and the table stores a handle - see the STRINGS note in md_system.h.
static void vlx_publish_str(md_system_t* sys, str_t path, str_t label, str_t value) {
	if (str_empty(value)) {
		return;
	}
	md_attribute_format_t format = {
		.type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 1, .shape = { 1 },
	};
	vlx_publish(sys, path, label, md_unit_none(), format, &value, sizeof(str_t));
}

// rank 2 {A,B}, one scalar per (a,b), row major with b fastest - which is the layout md_vlx.h
// documents for the 2D response quantities, so no rearrangement happens here.
static void vlx_publish_matrix(md_system_t* sys, str_t path, str_t label, md_unit_t unit, const double* values, size_t rows, size_t cols) {
	if (!values || rows == 0 || cols == 0) {
		return;
	}
	md_attribute_format_t format = {
		.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2, .shape = { (uint32_t)rows, (uint32_t)cols },
	};
	vlx_publish(sys, path, label, unit, format, values, rows * cols * sizeof(double));
}

// Publishes one COLUMN of an array of structs. The source is strided and an attribute is
// contiguous, so the values are gathered into the table's own storage through md_attributes_data
// rather than through a temporary which is then copied again.
//
// A struct of mixed types is not one attribute - a value has ONE type - so a record with six fields
// becomes six sibling paths over the same index space. That is what the ATTRIBUTES note means by
// independent quantities being sibling paths: the transposition from the file's row layout is the
// whole cost, and it happens once, here.
static void vlx_publish_column(md_system_t* sys, str_t path, str_t label, md_unit_t unit, md_attribute_type_t type, const void* base, size_t stride, size_t count) {
	if (!base || count == 0) {
		return;
	}

	md_attribute_format_t format = {
		.type = type, .components = 1, .rank = 1, .shape = { (uint32_t)count },
	};
	md_attribute_id_t id = vlx_publish(sys, path, label, unit, format, NULL, 0);
	if (id == MD_ATTRIBUTE_INVALID) {
		return;
	}

	uint8_t* dst = (uint8_t*)md_attributes_data(&sys->attributes, id, type);
	if (!dst) {
		md_attributes_remove(&sys->attributes, id);
		return;
	}

	const size_t   elem_size = md_attribute_type_size(type);
	const uint8_t* src       = (const uint8_t*)base;
	for (size_t i = 0; i < count; ++i) {
		MEMCPY(dst + i * elem_size, src + i * stride, elem_size);
	}
}


// The anchor of a dipole group: rank 0, one 3 component value, constant over whatever index space
// the group's vector has. Angstrom because it is a point in system space, unlike the moment.
//
// It is NOT replicated to match the vector's shape. Group members share an index space, not a
// shape; storing the same three numbers once per excited state would be N copies with nothing
// keeping them equal, to save a consumer one line.
static void vlx_publish_origin(md_system_t* sys, str_t path, dvec3_t origin) {
	md_attribute_format_t format = {
		.type = MD_ATTRIBUTE_TYPE_F64, .components = 3, .rank = 0,
	};
	vlx_publish(sys, path, (str_t){0}, md_unit_angstrom(), format, &origin, 3 * sizeof(double));
}

// Forward declarations for the reader-internal accessors below. All static: the vlx object is a
// parse-time scratch representation and nothing outside this file names it.
static vlx_t* vlx_create(struct md_allocator_i* backing, struct md_system_t* sys);
static size_t vlx_number_of_atoms(const vlx_t* vlx);
static size_t vlx_number_of_electrons(const vlx_t* vlx, vlx_spin_t spin);
static const dvec3_t* vlx_atom_coordinates(const vlx_t* vlx);
static const uint8_t* vlx_atomic_numbers(const vlx_t* vlx);
static const int* vlx_local_to_global_atom_idx(const vlx_t* vlx);
static dvec3_t vlx_scf_ground_state_dipole_moment(const vlx_t* vlx);
static size_t vlx_scf_number_of_atomic_orbitals (const vlx_t* vlx);
static size_t vlx_scf_number_of_molecular_orbitals(const vlx_t* vlx);
static const double* vlx_scf_mo_occupancy(const vlx_t* vlx, vlx_spin_t spin);
static const double* vlx_scf_mo_energy(const vlx_t* vlx, vlx_spin_t spin);
static bool vlx_gto_basis_extract(md_gto_basis_t* out, const vlx_t* vlx, struct md_allocator_i* alloc);
static const double* vlx_scf_mo_coefficients(const vlx_t* vlx, size_t mo_idx, vlx_spin_t spin);
static size_t vlx_rsp_nto_coefficients_extract(double* out_coefficients, double* out_lambdas, const vlx_t* vlx, size_t state_idx, vlx_nto_type_t type, size_t lambda_count);
static size_t vlx_scf_overlap_matrix_size(const vlx_t* vlx);
static const double* vlx_scf_overlap_matrix_data(const vlx_t* vlx);
static size_t vlx_rsp_number_of_excited_states(const vlx_t* vlx);
static const dvec3_t* vlx_rsp_electric_transition_dipole_moments(const vlx_t* vlx);
static const dvec3_t* vlx_rsp_magnetic_transition_dipole_moments(const vlx_t* vlx);
static const dvec3_t* vlx_rsp_velocity_transition_dipole_moments(const vlx_t* vlx);
static bool vlx_rsp_has_nto(const vlx_t* vlx);
static size_t vlx_rsp_nto_lambdas_extract(double* out_lambdas, const vlx_t* vlx, size_t state_idx, size_t lambda_count);
static bool vlx_system_begin(vlx_t* vlx, md_system_state_t* state);
// Publishes what could NOT be published as it was read.
//
// Every block that is a straight pass-through - the SCF history, the response series, the
// vibrational and optimisation blocks, XPS, the per atom properties - is published by its own
// reader, straight into the storage the attribute table reserved. What is left here is everything
// whose shape or content is not known until the whole file is in hand: the AO data, which has to be
// permuted into shell order and converted out of the spherical basis before it means anything; the
// densities and transition densities computed from it; the NTOs; and the dipole groups, whose origin
// is a property of the molecule rather than of any one block.
static void vlx_publish_whole_file_attributes(struct md_system_t* sys, const vlx_t* vlx);

static int char_to_angular_momentum_type(int c) {
	switch (c) {
	case 'S': return 0;
	case 'P': return 1;
	case 'D': return 2;
	case 'F': return 3;
	case 'G': return 4;
	default: return -1;
	}
}

static inline basis_set_basis_t* basis_set_get_atom_basis(const basis_set_t* basis_set, int atomic_number) {
	if (atomic_number < (int)basis_set->atom_basis.count) {
		return basis_set->atom_basis.data + atomic_number;
	}
	return NULL;
}

static inline int compute_max_angular_momentum(const basis_set_t* basis_set, const md_element_t* atomic_numbers, size_t count) {
	ASSERT(basis_set);
	ASSERT(atomic_numbers);
	int max_angl = 0;
	for (size_t i = 0; i < count; ++i) {
		const basis_set_basis_t* atom_basis = basis_set_get_atom_basis(basis_set, atomic_numbers[i]);
		max_angl = MAX(max_angl, (int)atom_basis->max_type);
	}
	return max_angl;
}

#define d3  3.464101615137754587
#define f5  1.581138830084189666
#define f15 7.745966692414833770
#define f3  1.224744871391589049
#define g35 4.0 * 5.916079783099616042
#define g17 4.0 * 4.183300132670377739
#define g5  4.0 * 2.236067977499789696
#define g2  4.0 * 1.581138830084189666

static const double		S_factors[] = {1.0};
static const uint8_t    S_indices[] = {0};
static const uint8_t	S_num_fac[] = {1};

static const double		P_factors[] = {1.0, 1.0, 1.0};
static const uint8_t	P_indices[] = {1, 2, 0};
static const uint8_t	P_offsets[] = {0, 1, 2};
static const uint8_t	P_num_fac[] = {1, 1, 1};

static const double		D_factors[] = {d3, d3, -1.0, -1.0, 2.0, d3, 0.5 * d3, -0.5 * d3};
static const uint8_t    D_indices[] = {1, 4, 0, 3, 5, 2, 0, 3};
static const uint8_t    D_offsets[] = {0, 1, 2, 5, 6};
static const uint8_t	D_num_fac[] = {1, 1, 3, 1, 2};

static const double		F_factors[] = {3.0 * f5, -f5, f15, 4.0 * f3, -f3, -f3, 2.0, -3.0, -3.0, 4.0 * f3, -f3, -f3, 0.5 * f15, -0.5 * f15, f5, -3.0 * f5};
static const uint8_t    F_indices[] = {1, 6, 4, 8, 1, 6, 9, 2, 7, 5, 0, 3, 2, 7, 0, 3};
static const uint8_t	F_offsets[] = {0, 2, 3, 6, 9, 12, 14};
static const uint8_t	F_num_fac[] = {2, 1, 3, 3, 3, 2, 2};

static const double		G_factors[] = {
	g35, -g35, 3.0 * g17, -g17, 6.0 * g5, -g5, -g5, 4.0 * g2, -3.0 * g2, -3.0 * g2,
	8.0, 3.0, 3.0, 6.0, -24.0, -24.0, 4.0 * g2, -3.0 * g2, -3.0 * g2, 3.0 * g5,
	-3.0 * g5, -0.5 * g5, 0.5 * g5,  g17,  -3.0 * g17, 0.25 * g35, 0.25 * g35, -1.50 * g35};
static const uint8_t    G_indices[] = {1, 6, 4, 11, 8, 1, 6, 13, 4, 11, 14, 0, 10, 3, 5, 12, 9, 2, 7, 5, 12, 0, 10, 2, 7, 0, 10, 3};
static const uint8_t	G_offsets[] = {0, 2, 4, 7, 10, 16, 19, 23, 25};
static const uint8_t	G_num_fac[] = {2, 2, 3, 3, 6, 3, 4, 2, 3};

#undef d3
#undef f5
#undef f15
#undef f3
#undef g35
#undef g17
#undef g5 
#undef g2 

static inline int spherical_momentum_num_components(int angl) {
	return angl * 2 + 1;
}

static inline int spherical_momentum_num_factors(int angl, int isph) {
	switch(angl) {
	case 0:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(S_num_fac));
		return S_num_fac[isph];
	case 1:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(P_num_fac));
		return P_num_fac[isph];
	case 2:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(D_num_fac));
		return D_num_fac[isph];
	case 3:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(F_num_fac));
		return F_num_fac[isph];
	case 4:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(G_num_fac));
		return G_num_fac[isph];
	default:
		ASSERT(false);
		return 0;
	}
}

static inline const double* spherical_momentum_factors(int angl, int isph) {
	switch(angl) {
	case 0:
		ASSERT(isph == 0);
		return S_factors;
	case 1:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(P_offsets));
		return P_factors + P_offsets[isph];
	case 2:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(D_offsets));
		return D_factors + D_offsets[isph];
	case 3:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(F_offsets));
		return F_factors + F_offsets[isph];
	case 4:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(G_offsets));
		return G_factors + G_offsets[isph];
	default:
		ASSERT(false);
		return NULL;
	}
}

static inline const uint8_t* spherical_momentum_indices(int angl, int isph) {
	switch(angl) {
	case 0:
		ASSERT(isph == 0);
		return S_indices;
	case 1:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(P_offsets));
		return P_indices + P_offsets[isph];
	case 2:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(D_offsets));
		return D_indices + D_offsets[isph];
	case 3:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(F_offsets));
		return F_indices + F_offsets[isph];
	case 4:
		ASSERT(0 <= isph && (size_t)isph < ARRAY_SIZE(G_offsets));
		return G_indices + G_offsets[isph];
	default:
		ASSERT(false);
		return NULL;
	}
}

typedef uint8_t lmn_t[3];

// S: 0
static const lmn_t S_lmn[1] = {{0,0,0}};
// P: x y z
static const lmn_t P_lmn[3] = {{1,0,0}, {0,1,0}, {0,0,1}};
// D: xx xy xz yy yz zz
static const lmn_t D_lmn[6] = {{2,0,0}, {1,1,0}, {1,0,1}, {0,2,0}, {0,1,1}, {0,0,2}};
// F: xxx xxy xxz xyy xyz xzz yyy yyz yzz zzz
static const lmn_t F_lmn[10] = {{3,0,0}, {2,1,0}, {2,0,1}, {1,2,0}, {1,1,1}, {1,0,2}, {0,3,0}, {0,2,1}, {0,1,2}, {0,0,3}};
// G: xxxx xxxy xxxz xxyy xxyz xxzz xyyy xyyz xyzz xzzz yyyy yyyz yyzz yzzz zzzz
static const lmn_t G_lmn[15] = {{4,0,0}, {3,1,0}, {3,0,1}, {2,2,0}, {2,1,1}, {2,0,2}, {1,3,0}, {1,2,1}, {1,1,2}, {1,0,3}, {0,4,0}, {0,3,1}, {0,2,2}, {0,1,3}, {0,0,4}};

static inline const lmn_t* cartesian_angular_momentum(int angl) {
	switch (angl) {
	case 0: return S_lmn;
	case 1: return P_lmn;
	case 2: return D_lmn;
	case 3: return F_lmn;
	case 4: return G_lmn;
	default: ASSERT(false); return NULL;
	}
}

typedef struct basis_func_t {
    int type;
    int count;
    double* exponents;
    double* normalization_coefficients;
} basis_func_t;

static inline basis_func_t get_basis_func(const basis_set_t* basis_set, int basis_func_idx) {
	basis_set_func_t func = basis_set->basis_func.data[basis_func_idx];
	return (basis_func_t) {
		.type = func.type,
		.count = func.param_count,
		.exponents = basis_set->param.exponents + func.param_offset,
		.normalization_coefficients = basis_set->param.normalization_coefficients + func.param_offset,
	};
}

static size_t basis_set_extract_atomic_basis_func_angl(basis_func_t* out_funcs, size_t cap_funcs, const basis_set_t* basis_set, int atomic_number, int angl) {
    size_t count = 0;

    basis_set_basis_t* atom_basis = basis_set_get_atom_basis(basis_set, atomic_number);
    if (atom_basis) {
        int beg = atom_basis->basis_func_offset;
        int end = atom_basis->basis_func_offset + atom_basis->basis_func_count;
        for (int i = beg; i < end; ++i) {
            if (count == cap_funcs) break;
            if (basis_set->basis_func.data[i].type == angl) {
                out_funcs[count++] = get_basis_func(basis_set, i);
            }
        }
    }

    return count;
}

static size_t basis_set_count_atomic_basis_func(const basis_set_t* basis_set, int atomic_number, int angl) {
    size_t count = 0;
    basis_set_basis_t* atom_basis = basis_set_get_atom_basis(basis_set, atomic_number);
    if (atom_basis) {
        int beg = atom_basis->basis_func_offset;
        int end = atom_basis->basis_func_offset + atom_basis->basis_func_count;
        for (int i = beg; i < end; ++i) {
            if (basis_set->basis_func.data[i].type == angl) {
                count++;
            }
        }
    }
    return count;
}

static size_t extract_ao_to_atom_idx(int* out_ao_to_atom, const md_atomic_number_t* atomic_numbers, size_t number_of_atoms, const basis_set_t* basis_set) {
	int natoms = (int)number_of_atoms;
	int max_angl = compute_max_angular_momentum(basis_set, atomic_numbers, number_of_atoms);

	size_t count = 0;

	basis_func_t basis_funcs[128];

	// azimuthal quantum number: s,p,d,f,...
	for (int angl = 0; angl <= max_angl; angl++) {
		//CSphericalMomentum sphmom(angl);
		int nsph = spherical_momentum_num_components(angl);
		// magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...
		for (int isph = 0; isph < nsph; isph++) {
			// int	ncomp = spherical_momentum_num_factors(angl, isph);

			// go through atoms
			for (int atomidx = 0; atomidx < natoms; atomidx++) {
				int idelem = atomic_numbers[atomidx];
				size_t num_ao = basis_set_extract_atomic_basis_func_angl(basis_funcs, ARRAY_SIZE(basis_funcs), basis_set, idelem, angl);

				for (size_t iao = 0; iao < num_ao; iao++) {
					if (out_ao_to_atom) {
						out_ao_to_atom[count] = atomidx;
					}
					count += 1;
				}
			}
		}
	}
	return count;
}

static size_t compute_basis_num_atomic_orbitals(const vlx_t* vlx) {
	ASSERT(vlx);
	if (!vlx->basis_set.atom_basis.count || !vlx->atomic_numbers || vlx->number_of_atoms == 0) {
		return 0;
	}
	return extract_ao_to_atom_idx(NULL, vlx->atomic_numbers, vlx->number_of_atoms, &vlx->basis_set);
}

// This is a ported reference implementation from VeloxChem found in VisualizationDriver.cpp
static size_t compPhiAtomicOrbitals(double* out_phi, size_t phi_cap,
	const dvec3_t* atom_coordinates, const md_element_t* atomic_numbers, size_t num_atoms,
	const basis_set_t* basis_set,
	double xp,
	double yp,
	double zp)
{
	int natoms = (int)num_atoms;
	int max_angl = compute_max_angular_momentum(basis_set, atomic_numbers, num_atoms);

	size_t count = 0;

	basis_func_t basis_funcs[128];

	// azimuthal quantum number: s,p,d,f,...
	for (int aoidx = 0, angl = 0; angl <= max_angl; angl++) {
		//CSphericalMomentum sphmom(angl);
		int nsph = spherical_momentum_num_components(angl);
		const lmn_t* lmn = cartesian_angular_momentum(angl);
		// magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...
		for (int isph = 0; isph < nsph; isph++) {
			// prepare Cartesian components (Maximum number of components should be 6 here for the currently supported basis set)
			double lx[8];
			double ly[8];
			double lz[8];
			int			      ncomp = spherical_momentum_num_factors(angl, isph);
			const double*	fcarts  = spherical_momentum_factors(angl, isph);
			const uint8_t*	indices = spherical_momentum_indices(angl, isph);

			for (int icomp = 0; icomp < ncomp; icomp++) {
				int cartind = indices[icomp];

				lx[icomp] = lmn[cartind][0];
				ly[icomp] = lmn[cartind][1];
				lz[icomp] = lmn[cartind][2];
			}

			// go through atoms

			for (int atomidx = 0; atomidx < natoms; atomidx++) {
				// process coordinates
				// Conversion from Ångström to Bohr
				double rx = (xp - atom_coordinates[atomidx].x) * ANGSTROM_TO_BOHR;
				double ry = (yp - atom_coordinates[atomidx].y) * ANGSTROM_TO_BOHR;
				double rz = (zp - atom_coordinates[atomidx].z) * ANGSTROM_TO_BOHR;
				double r2 = rx*rx + ry*ry + rz*rz;

				// process atomic orbitals
				int idelem = atomic_numbers[atomidx];

				size_t num_basis_funcs = basis_set_extract_atomic_basis_func_angl(basis_funcs, ARRAY_SIZE(basis_funcs), basis_set, idelem, angl);
				for (size_t funcidx = 0; funcidx < num_basis_funcs; funcidx++, aoidx++) {
					double phiao = 0.0;

					basis_func_t bf = basis_funcs[funcidx];

					// process primitives
					for (int iprim = 0; iprim < bf.count; iprim++) {
						double expon = exp(-bf.exponents[iprim] * r2);
						double coef1 = bf.normalization_coefficients[iprim];

						// transform from Cartesian to spherical harmonics
						for (int icomp = 0; icomp < ncomp; icomp++) {
							double coef2 = coef1 * fcarts[icomp];
							double powxyz = pow(rx, lx[icomp]) * pow(ry, ly[icomp]) * pow(rz, lz[icomp]);
							phiao += coef2 * powxyz * expon;
						}
					}

					out_phi[count++] = phiao;
					if (count == phi_cap) {
						return count;
					}
				}
			}
		}
	}

	return count;
}

// Build a permutation table that maps from shell order (angl→atom→func→isph)
// to VeloxChem matrix row order (angl→isph→atom→func).
// ao_remap[shell_ao_idx] = vlx_ao_idx.
// Returns the total number of AOs (length of the table), or 0 on failure.
static bool build_ao_remap(int* out_remap, size_t capacity, const vlx_t* vlx) {
	ASSERT(out_remap);
	ASSERT(vlx);

	int natoms   = (int)vlx->number_of_atoms;
	int max_angl = compute_max_angular_momentum(&vlx->basis_set, vlx->atomic_numbers, vlx->number_of_atoms);

	// First, count total AOs
	size_t num_aos = 0;
	for (int angl = 0; angl <= max_angl; angl++) {
		int nsph = spherical_momentum_num_components(angl);
		for (int atomidx = 0; atomidx < natoms; atomidx++) {
			int idelem = vlx->atomic_numbers[atomidx];
			size_t num_funcs = basis_set_count_atomic_basis_func(&vlx->basis_set, idelem, angl);
			num_aos += (size_t)nsph * num_funcs;
		}
	}

	if (num_aos != capacity) {
		MD_LOG_ERROR("Capacity of remap table did not match the number of AOs found in the basis set");
		return false;
	}

	MEMSET(out_remap, 0, sizeof(int) * num_aos);

	// VeloxChem matrix AO index: angl → isph → atom → func
	// Shell AO index:            angl → atom → func → isph
	//
	// We walk shell order (outer loop) and record where each entry maps in VLX order.
	// vlx_ao_start[angl][isph][atom] is needed. We pre-compute the vlx base offset per
	// (angl, isph, atom) by walking the vlx ordering once.

	// Compute vlx_base[angl][isph][atomidx] = starting VLX AO index for that group.
	// Flat: vlx_base[angl * (max_angl+1) * natoms + isph * natoms + atomidx]
	// But isph is up to 2*max_angl+1. Use a temp VLA-style allocation.
	int max_nsph = 2 * max_angl + 1;
	md_temp_scope_t temp = md_temp_begin();
	int* vlx_base = (int*)md_temp_alloc(temp, sizeof(int) * (max_angl + 1) * max_nsph * natoms);
	int* vlx_num  = (int*)md_temp_alloc(temp, sizeof(int) * (max_angl + 1) * natoms);
	MEMSET(vlx_base, 0, sizeof(int) * (max_angl + 1) * max_nsph * natoms);
	MEMSET(vlx_num,  0, sizeof(int) * (max_angl + 1) * natoms);

	// Count funcs per (angl, atom) so we know stride for vlx_base
	for (int angl = 0; angl <= max_angl; angl++) {
		for (int atomidx = 0; atomidx < natoms; atomidx++) {
			int idelem = vlx->atomic_numbers[atomidx];
			basis_func_t bf[128];
			size_t nf = basis_set_extract_atomic_basis_func_angl(bf, ARRAY_SIZE(bf), &vlx->basis_set, idelem, angl);
			vlx_num[angl * natoms + atomidx] = (int)nf;
		}
	}

	// Walk VLX order to fill vlx_base
	int vlx_idx = 0;
	for (int angl = 0; angl <= max_angl; angl++) {
		int nsph = spherical_momentum_num_components(angl);
		for (int isph = 0; isph < nsph; isph++) {
			for (int atomidx = 0; atomidx < natoms; atomidx++) {
				vlx_base[angl * max_nsph * natoms + isph * natoms + atomidx] = vlx_idx;
				vlx_idx += vlx_num[angl * natoms + atomidx];
			}
		}
	}

	// Now walk shell order and fill remap
	int shell_idx = 0;
	for (int angl = 0; angl <= max_angl; angl++) {
		int nsph = spherical_momentum_num_components(angl);
		for (int atomidx = 0; atomidx < natoms; atomidx++) {
			int nfuncs = vlx_num[angl * natoms + atomidx];
			for (int funcidx = 0; funcidx < nfuncs; funcidx++) {
				for (int isph = 0; isph < nsph; isph++, shell_idx++) {
					int base = vlx_base[angl * max_nsph * natoms + isph * natoms + atomidx];
					out_remap[shell_idx] = base + funcidx;
				}
			}
		}
	}

	md_temp_end(temp);
	return true;
}

static size_t vlx_pgto_count(const vlx_t* vlx) {
	int natoms = (int)vlx->number_of_atoms;
	int max_angl = compute_max_angular_momentum(&vlx->basis_set, vlx->atomic_numbers, vlx->number_of_atoms);

	size_t count = 0;

	basis_func_t basis_funcs[256];

	// azimuthal quantum number: s,p,d,f,...
	for (int angl = 0; angl <= max_angl; angl++) {
		//CSphericalMomentum sphmom(angl);
		int nsph = spherical_momentum_num_components(angl);
		// magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...
		for (int isph = 0; isph < nsph; isph++) {
			int	ncomp = spherical_momentum_num_factors(angl, isph);
			// go through atoms
			for (int atomidx = 0; atomidx < natoms; atomidx++) {
				int idelem = vlx->atomic_numbers[atomidx];

				// process atomic orbitals
                size_t num_basis_funcs = basis_set_extract_atomic_basis_func_angl(basis_funcs, ARRAY_SIZE(basis_funcs), &vlx->basis_set, idelem, angl);
				for (size_t funcidx = 0; funcidx < num_basis_funcs; funcidx++) {
					// process primitives
					count += basis_funcs[funcidx].count * ncomp;
				}
			}
		}
	}

	return count;
}

static inline double compute_overlap(basis_func_t func, int i, int j) {
	const double fab  = 1.0 / (func.exponents[i] + func.exponents[j]);
	const double fab2 = fab * fab;
	const double ovl = func.normalization_coefficients[i] * func.normalization_coefficients[j] * pow(PI * fab, 1.5);

	switch (func.type) {
	case 0: return ovl;
	case 1: return 0.5 * fab * ovl;
	case 2: return 3.0 * fab2 * ovl;
	case 3: return 7.5 * fab2 * fab * ovl;
	case 4: return 420.0 * fab2 * fab2 * ovl;
	case 5: return 1890.0 * fab2 * fab2 * fab * ovl;
	case 6: return 41580.0 * fab2 * fab2 * fab2 * ovl;
	default:
		ASSERT(false);
		return 0;
	}
}

static void rescale_basis_func(basis_func_t func) {
	const double fpi = 2.0 / PI;

	for (int i = 0; i < func.count; i++) {
		func.normalization_coefficients[i] *= pow(func.exponents[i] * fpi, 0.75);
	}

	if (func.type < 0 || 6 < func.type) {
		MD_LOG_DEBUG("Invalid basis function type supplied in rescaling");
		return;
	}

	static const double f_table[] = {
		0,
		2.0,
		1.15470053837925152902, // 2.0 / sqrt(3.0)
		1.03279555898864450271, // 4.0 / sqrt(15.0)
		0.19518001458970663587, // 2.0 / sqrt(105.0)
		0.13012000972647109058, // 4.0 / sqrt(945.0)
		0.03923265908909997910, // 4.0 / sqrt(10395.0)
	};

	double f = f_table[func.type];
	double e = (double)func.type * 0.5;

	for (int i = 0; i < func.count; i++) {
		func.normalization_coefficients[i] *= pow(f * func.exponents[i], e);
	}
}

static void normalize_basis_set(basis_set_t* basis_set) {
	for (size_t func_idx = 0; func_idx < basis_set->basis_func.count; ++func_idx) {
		basis_func_t func = get_basis_func(basis_set, (int)func_idx);
		// uncontracted basis, set expansion coeficient to 1.0
		if (func.count == 1) func.normalization_coefficients[0] = 1.0;

		// normalize primitive GBFs
		rescale_basis_func(func);

		// compute overlap
		double ovl = 0.0;
		for (int i = 0; i < func.count; i++) {
			ovl += compute_overlap(func, i, i);
			for (int j = i + 1; j < func.count; j++) {
				ovl += 2.0 * compute_overlap(func, i, j);
			}
		}

		// renormalize primitive BFs
		ovl = 1.0 / sqrt(ovl);
		for (int i = 0; i < func.count; i++) {
			func.normalization_coefficients[i] *= ovl;
		}
	}
}

static bool parse_basis_set(basis_set_t* basis_set, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	MEMSET(basis_set, 0, sizeof(basis_set_t));

	str_t line;
	str_t tok[4];
	size_t line_count = 0;

	// Insert null_basis element for index 0
	const basis_set_basis_t null_basis = {0};

	basis_set_basis_t* curr_atom_basis = NULL;
	while (md_buffered_reader_extract_line(&line, reader)) {
		line_count += 1;
		str_t line_original = line;
		size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (!num_tok) continue;

		if ((num_tok == 2 || num_tok == 3) && str_eq(tok[0], STR_LIT("@BASIS_SET"))) {
			// In the case of renamed identifiers using the alias table, there is an additional token which provides the original identifier (prepended with an !)
			str_t ident = {0};
			if (num_tok == 2)
				ident = tok[1];
			else {
				ident = str_substr(tok[2], 1, SIZE_MAX);
			}
			MD_LOG_DEBUG("Parsing Basis Set with identifier: '" STR_FMT "'", STR_ARG(ident));
			basis_set->identifier = str_copy(ident, alloc);
		}
		else if (num_tok == 2 && str_eq(tok[0], STR_LIT("@ATOMBASIS"))) {
			int atomic_number = md_atomic_number_from_symbol(tok[1], true);
			if (atomic_number == 0) {
				MD_LOG_ERROR("Unrecognized element '" STR_FMT "' in basis set", STR_ARG(tok[1]));
				return false;
			}
			basis_set_basis_t atom_basis = {
				.max_type = 0,
				.basis_func_count = 0,
				.basis_func_offset = (uint16_t)basis_set->basis_func.count,
			};

			// Grow the array and fill in slots with null_basis
			while ((int)md_array_size(basis_set->atom_basis.data) < atomic_number) {
				md_array_push(basis_set->atom_basis.data, null_basis, alloc);
			}

			md_array_push(basis_set->atom_basis.data, atom_basis, alloc);
			curr_atom_basis = md_array_last(basis_set->atom_basis.data);

			basis_set->atom_basis.count = md_array_size(basis_set->atom_basis.data);
		}
		else if (num_tok == 1 && str_eq(tok[0], STR_LIT("@END"))) {
			curr_atom_basis = NULL;
		}
		else if (num_tok == 3) {
			int type = char_to_angular_momentum_type(tok[0].ptr[0]);
			if (type == -1) {
				MD_LOG_ERROR("Unrecognized angular momentum type '" STR_FMT "' in basis set", STR_ARG(tok[0]));
				MD_LOG_ERROR("This occured on line %zu: '" STR_FMT "'", line_count, STR_ARG(line_original));

				return false;
			}
			int count = (int)parse_int(tok[1]);
			if (count == 0 || count > 255) {
				MD_LOG_ERROR("Invalid number of coefficients in atom basis in basis set");
				return false;
			}

			if (curr_atom_basis == NULL) {
				MD_LOG_ERROR("No atom basis has been defined for supplied coefficients");
				return false;
			}

			// We have a new basis function for the current atom basis
			basis_set_func_t basis_func = {
				.type = (uint8_t)type,
				.param_count = (uint8_t)count,
				.param_offset = (uint16_t)basis_set->param.count,
			};
			md_array_push(basis_set->basis_func.data, basis_func, alloc);
			basis_set->basis_func.count += 1;

			for (int i = 0; i < count; ++i) {
				if (!md_buffered_reader_extract_line(&line, reader)) {
					MD_LOG_ERROR("Failed to parse coefficients in atom basis function");
					return false;
				}
				num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
				if (num_tok != 2) {
					MD_LOG_ERROR("Unexpected number of parameters in atom basis function coefficients");
					return false;
				}
				char buf[64];
				str_copy_to_char_buf(buf, sizeof(buf), tok[0]);

				double exponent = parse_float(tok[0]);
				double coeff    = parse_float(tok[1]);
				
				md_array_push(basis_set->param.exponents, exponent, alloc);
				md_array_push(basis_set->param.normalization_coefficients, coeff, alloc);
				basis_set->param.count += 1;
			}

			curr_atom_basis->basis_func_count += 1;
			curr_atom_basis->max_type = MAX(curr_atom_basis->max_type, (uint8_t)type);
		}
	}

	return true;
}

// HDF5 prints a full diagnostic stack to stderr for every failed call. This reader
// probes for optional data constantly and reports its own failures through MD_LOG_*
// with the offending field name, so the automatic handler is pure noise. Silence it
// for the duration of a parse and restore whatever the host application had set --
// mdlib must not leave global HDF5 state modified.
typedef struct h5_error_scope_t {
	H5E_auto2_t func;
	void*       client_data;
} h5_error_scope_t;

static h5_error_scope_t h5_error_scope_begin(void) {
	h5_error_scope_t scope = {0};
	H5Eget_auto2(H5E_DEFAULT, &scope.func, &scope.client_data);
	H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
	return scope;
}

static void h5_error_scope_end(h5_error_scope_t scope) {
	H5Eset_auto2(H5E_DEFAULT, scope.func, scope.client_data);
}

// H5Lexists() only tolerates a missing *final* path component. If an intermediate
// group is absent it fails outright and, with the default error handler installed,
// dumps a full diagnostic stack to stderr. Nearly every probe in this file is for
// optional data (a file without TPA has no "tpa_strengths" group at all), so walk
// the path one component at a time and stop at the first miss.
//
// Returns false both for "not present" and for a genuine error; callers here only
// care whether the data is usable.
static bool h5_link_exists(hid_t loc_id, const char* path) {
	if (!path || !*path) return false;

	const size_t len = strlen(path);

	char buf[512];
	if (len >= sizeof(buf)) {
		MD_LOG_ERROR("HDF5 path exceeds %zu characters: '%s'", sizeof(buf) - 1, path);
		return false;
	}
	MEMCPY(buf, path, len + 1);

	for (size_t i = 0; i < len; ++i) {
		if (buf[i] != '/') continue;
		if (i == 0 || buf[i - 1] == '/') continue;  // leading or repeated separator

		buf[i] = '\0';
		const htri_t exists = H5Lexists(loc_id, buf, H5P_DEFAULT);
		buf[i] = '/';

		if (exists <= 0) return false;
	}

	return H5Lexists(loc_id, buf, H5P_DEFAULT) > 0;
}

static H5I_type_t h5_get_object_type(hid_t loc_id, const char* name) {
    if (!h5_link_exists(loc_id, name)) {
        return H5I_UNINIT;
    }

    hid_t obj_id = H5Oopen(loc_id, name, H5P_DEFAULT);
    if (obj_id < 0) {
        return H5I_BADID;
    }

    H5I_type_t type = H5Iget_type(obj_id);
    H5Oclose(obj_id);
    return type;
}

// Number of elements in a dataset's dataspace. Returns false on failure.
//
// Every read below uses H5S_ALL for both the memory and file selection, which makes
// H5Dread write *the whole dataset* into the caller's buffer. Any read into a fixed
// size destination must therefore check this first, or a file with an unexpected
// shape silently overruns the destination.
static bool h5_dataset_num_elements(hsize_t* out_count, hid_t dataset_id, const char* field_name) {
	hid_t space_id = H5Dget_space(dataset_id);
	if (space_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to query H5 dataspace for dataset: '%s'", field_name);
		return false;
	}

	const hssize_t npoints = H5Sget_simple_extent_npoints(space_id);
	H5Sclose(space_id);

	if (npoints < 0) {
		MD_LOG_ERROR("Failed to query element count for H5 dataset: '%s'", field_name);
		return false;
	}

	*out_count = (hsize_t)npoints;
	return true;
}

// Reads a dataset expected to hold exactly one element. Returns false if it holds
// anything else, rather than overrunning 'buf'.
static bool h5_read_scalar(void* buf, hid_t file_id, hid_t mem_type_id, const char* field_name) {
	if (!h5_link_exists(file_id, field_name)) {
		return false;
	}

	// Open the dataset containing the double value
	hid_t dataset_id = H5Dopen(file_id, field_name, H5P_DEFAULT);
	if (dataset_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 dataset: '%s'", field_name);
		return false;
	}

	bool result = false;

	hsize_t num_elem = 0;
	if (!h5_dataset_num_elements(&num_elem, dataset_id, field_name)) {
		goto done;
	}
	if (num_elem != 1) {
		MD_LOG_ERROR("Expected a scalar H5 dataset for '%s', got %llu elements", field_name, (unsigned long long)num_elem);
		goto done;
	}

	// Read the dataset into the 'value' variable
	herr_t status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf);
	if (status < 0) {
		MD_LOG_ERROR("Failed to read data for H5 dataset: '%s'", field_name);
		goto done;
	}

	result = true;
done:
	H5Dclose(dataset_id);

	return result;
}

static bool h5_read_str(str_t* str, hid_t file_id, const char* field_name, md_allocator_i* alloc) {
	bool result = false;

	if (!h5_link_exists(file_id, field_name)) {
		return false;
	}

	hid_t dataset_id = H5Dopen(file_id, field_name, H5P_DEFAULT);
	if (dataset_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 dataset: '%s'", field_name);
		return false;
	}

	hid_t datatype_id = H5Dget_type(dataset_id);
	hid_t space_id = H5Dget_space(dataset_id);

	if (datatype_id == H5I_INVALID_HID || space_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to query H5 datatype/space for dataset: '%s'", field_name);
		goto done;
	}

	if (H5Tget_class(datatype_id) != H5T_STRING) {
		MD_LOG_ERROR("H5 dataset is not a string: '%s'", field_name);
		goto done;
	}

	{
		// Both branches below read with H5S_ALL into storage for a single string.
		hsize_t num_elem = 0;
		if (!h5_dataset_num_elements(&num_elem, dataset_id, field_name)) {
			goto done;
		}
		if (num_elem != 1) {
			MD_LOG_ERROR("Expected a single string in H5 dataset '%s', got %llu", field_name, (unsigned long long)num_elem);
			goto done;
		}
	}

	if (H5Tis_variable_str(datatype_id)) {
		char* tmp = NULL;
		herr_t status = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmp);
		if (status < 0) {
			MD_LOG_ERROR("Failed to read variable-length string for H5 dataset: '%s'", field_name);
			goto done;
		}

		size_t len = tmp ? strlen(tmp) : 0;
		str_t data = str_alloc(len, alloc);
		if (len > 0) {
			MEMCPY((char*)data.ptr, tmp, len);
		}
		*str = data;

		if (tmp) {
			H5free_memory(tmp);
		}
	} else {
		bool fixed_result = false;
		const size_t raw_len = H5Tget_size(datatype_id);
		md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
		char* tmp = md_temp_alloc(temp_scope, raw_len + 1);
		if (!tmp) {
			MD_LOG_ERROR("Failed to allocate temporary buffer for H5 dataset: '%s'", field_name);
			goto fixed_done;
		}
		MEMSET(tmp, 0, raw_len + 1);

		herr_t status = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
		if (status < 0) {
			MD_LOG_ERROR("Failed to read fixed-length string for H5 dataset: '%s'", field_name);
			goto fixed_done;
		}

		size_t len = strnlen(tmp, raw_len);
		str_t data = str_alloc(len, alloc);
		if (len > 0) {
			MEMCPY((char*)data.ptr, tmp, len);
		}
		*str = data;
		fixed_result = true;

	fixed_done:
		md_temp_end(temp_scope);
		if (!fixed_result) goto done;
	}

	result = true;

done:
	if (datatype_id != H5I_INVALID_HID) H5Tclose(datatype_id);
	if (space_id != H5I_INVALID_HID) H5Sclose(space_id);
	H5Dclose(dataset_id);

	return result;
}

static size_t h5_read_cstr(char* out_str, size_t str_cap, hid_t file_id, const char* field_name) {
    ASSERT(out_str);
	ASSERT(str_cap > 0);
	out_str[0] = '\0';

	if (!h5_link_exists(file_id, field_name)) {
		return 0;
	}

	// Open the dataset
	hid_t dataset_id = H5Dopen(file_id, field_name, H5P_DEFAULT);
	if (dataset_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 dataset: '%s'", field_name);
		return 0;
	}

	size_t result = 0;

	// Get the datatype and space
	hid_t datatype_id = H5Dget_type(dataset_id);  // Get datatype
	if (datatype_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to query H5 datatype for dataset: '%s'", field_name);
		goto done;
	}

	if (H5Tget_class(datatype_id) != H5T_STRING) {
		MD_LOG_ERROR("H5 dataset is not a string: '%s'", field_name);
		goto done;
	}

	{
		// Both branches below read with H5S_ALL into storage for a single string.
		hsize_t num_elem = 0;
		if (!h5_dataset_num_elements(&num_elem, dataset_id, field_name)) {
			goto done;
		}
		if (num_elem != 1) {
			MD_LOG_ERROR("Expected a single string in H5 dataset '%s', got %llu", field_name, (unsigned long long)num_elem);
			goto done;
		}
	}

	if (H5Tis_variable_str(datatype_id)) {
		// Variable-length string
		char* tmp = NULL;
		herr_t status = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &tmp);
		if (status < 0) {
			MD_LOG_ERROR("Failed to read variable-length string for H5 dataset: '%s'", field_name);
			goto done;
		}

		if (tmp) {
			size_t len = strlen(tmp);
			result = MIN(len, str_cap - 1);
			MEMCPY(out_str, tmp, result);
			out_str[result] = '\0';
			H5free_memory(tmp);
		}
	} else {
      // Fixed-length strings may not be null-terminated, so read into a temporary buffer first.
		bool fixed_result = false;
		size_t size = H5Tget_size(datatype_id);
		md_temp_scope_t temp = md_temp_begin();

		char* tmp = md_temp_alloc(temp, size + 1);
		if (!tmp) {
			MD_LOG_ERROR("Failed to allocate temporary buffer for H5 dataset: '%s'", field_name);
			goto fixed_done;
		}
		MEMSET(tmp, 0, size + 1);

		herr_t status = H5Dread(dataset_id, datatype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);
		if (status < 0) {
			MD_LOG_ERROR("Failed to read fixed-length string for H5 dataset: '%s'", field_name);
			goto fixed_done;
		}

		size_t len = strnlen(tmp, size);
		result = MIN(len, str_cap - 1);
		MEMCPY(out_str, tmp, result);
		out_str[result] = '\0';
		fixed_result = true;

	fixed_done:
		md_temp_end(temp);
		if (!fixed_result) goto done;
	}

done:
	// Close HDF5 resources
  if (datatype_id != H5I_INVALID_HID) H5Tclose(datatype_id);
	H5Dclose(dataset_id);

	return result;
}

#define H5_MAX_RANK 32

// Writes the extent of each dimension into dims[0 .. rank-1] and returns the rank.
//
// Returns 0 for "not present" and for every failure. Never returns a negative value:
// most call sites test this truthily, and a -1 would read as success.
static int h5_read_dataset_dims(size_t* dims, int max_dims, hid_t file_id, const char* field_name) {
	ASSERT(dims);

	if (max_dims <= 0 || max_dims > H5_MAX_RANK) {
		MD_LOG_ERROR("Invalid max_dims (%i) requested for H5 dataset: '%s'", max_dims, field_name);
		return 0;
	}

	if (!h5_link_exists(file_id, field_name)) {
		return 0;
	}

	// Open the dataset
	hid_t dataset_id = H5Dopen(file_id, field_name, H5P_DEFAULT);
	if (dataset_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 dataset: '%s'", field_name);
		return 0;
	}

	int   result   = 0;
	hid_t space_id = H5Dget_space(dataset_id);
	if (space_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to query H5 dataspace for dataset: '%s'", field_name);
		goto done;
	}

	int ndim = H5Sget_simple_extent_ndims(space_id);
	if (ndim < 0) {
		MD_LOG_ERROR("Failed to get number of dimensions for H5 dataset: '%s'", field_name);
		goto done;
	}

	if (ndim > max_dims) {
		MD_LOG_ERROR("H5 dataset '%s' has rank %i, caller supplied room for %i", field_name, ndim, max_dims);
		goto done;
	}

	// hsize_t is always 64-bit while size_t is not, so read into a correctly typed
	// buffer and widen/narrow explicitly rather than aliasing the caller's array.
	hsize_t extent[H5_MAX_RANK] = {0};
	ndim = H5Sget_simple_extent_dims(space_id, extent, NULL);
	if (ndim < 0) {
		MD_LOG_ERROR("Failed to get dimensions for H5 dataset: '%s'", field_name);
		goto done;
	}

	for (int i = 0; i < ndim; ++i) {
		if (extent[i] > (hsize_t)SIZE_MAX) {
			MD_LOG_ERROR("H5 dataset '%s' dimension %i does not fit in size_t", field_name, i);
			goto done;
		}
		dims[i] = (size_t)extent[i];
	}

	result = ndim;

done:
	if (space_id != H5I_INVALID_HID) H5Sclose(space_id);
	H5Dclose(dataset_id);

	return result;
}

static bool h5_check_dataset_exists(hid_t file_id, const char* field_name) {
	return h5_link_exists(file_id, field_name);
}

// Reads a string attribute into a caller-supplied buffer, always null terminated
// and always truncated to fit. Returns false if the attribute is absent, is not a
// string, or cannot be read; every handle opened here is released on every path.
//
// The size check on the fixed-length branch matters: the stored size is whatever
// the writer chose, so reading it straight into the destination overflows for any
// attribute longer than the caller's buffer.
static bool h5_read_string_attribute(char* out_buf, size_t out_cap, hid_t obj_id, const char* attr_name) {
	ASSERT(out_buf);
	ASSERT(out_cap > 0);
	out_buf[0] = '\0';

	// htri_t: negative on error, and negative is truthy.
	if (H5Aexists(obj_id, attr_name) <= 0) {
		return false;
	}

	hid_t attr_id = H5Aopen(obj_id, attr_name, H5P_DEFAULT);
	if (attr_id == H5I_INVALID_HID) {
		return false;
	}

	bool  result    = false;
	hid_t attr_type = H5Aget_type(attr_id);
	if (attr_type == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to query type of attribute '%s'", attr_name);
		goto done;
	}

	if (H5Tget_class(attr_type) != H5T_STRING) {
		MD_LOG_ERROR("Attribute '%s' is not a string", attr_name);
		goto done;
	}

	if (H5Tis_variable_str(attr_type)) {
		char* var_str = NULL;
		if (H5Aread(attr_id, attr_type, &var_str) < 0) {
			MD_LOG_ERROR("Failed to read variable-length string attribute '%s'", attr_name);
			goto done;
		}
		// On success HDF5 may still hand back NULL for an empty string.
		if (var_str) {
			const size_t len = MIN(strlen(var_str), out_cap - 1);
			MEMCPY(out_buf, var_str, len);
			out_buf[len] = '\0';
			H5free_memory(var_str);
		}
		result = true;
	} else {
		const size_t size = H5Tget_size(attr_type);
		if (size == 0) {
			MD_LOG_ERROR("Attribute '%s' has zero size", attr_name);
			goto done;
		}

		// Fixed-length strings are not necessarily null terminated.
		md_temp_scope_t temp = md_temp_begin();
		char* tmp = md_temp_alloc(temp, size + 1);
		if (tmp) {
			MEMSET(tmp, 0, size + 1);
			if (H5Aread(attr_id, attr_type, tmp) >= 0) {
				const size_t len = MIN(strnlen(tmp, size), out_cap - 1);
				MEMCPY(out_buf, tmp, len);
				out_buf[len] = '\0';
				result = true;
			} else {
				MD_LOG_ERROR("Failed to read fixed-length string attribute '%s'", attr_name);
			}
		}
		md_temp_end(temp);
	}

done:
	if (attr_type != H5I_INVALID_HID) H5Tclose(attr_type);
	H5Aclose(attr_id);
	return result;
}

// Checks a square AO-basis matrix for symmetry and forces it if it is close but not
// exact. Returns true if it was already symmetric to tolerance.
//
// Consumers of AO density matrices in this codebase read only the upper triangle, so
// an asymmetric matrix is not merely inaccurate -- half of it is discarded without a
// trace. Anything beyond rounding is reported with the offending magnitude so it is
// visible rather than absorbed.
static bool vlx_report_and_enforce_symmetry(double* mat, size_t dim, const char* label) {
	ASSERT(mat);

	double max_asym = 0.0;
	double max_abs  = 0.0;
	for (size_t i = 0; i < dim; ++i) {
		for (size_t j = i + 1; j < dim; ++j) {
			const double a = mat[i * dim + j];
			const double b = mat[j * dim + i];
			max_asym = MAX(max_asym, fabs(a - b));
			max_abs  = MAX(max_abs, MAX(fabs(a), fabs(b)));
		}
	}

	// Scale-relative, so this does not fire on accumulated rounding in a large matrix.
	const double tolerance = 1.0e-10 * (max_abs > 0.0 ? max_abs : 1.0);
	if (max_asym <= tolerance) {
		return true;
	}

	MD_LOG_INFO("Density property '%s' is not symmetric (max deviation %g, largest element %g). "
				"Symmetrizing: the density evaluation path only reads the upper triangle.",
				label, max_asym, max_abs);

	for (size_t i = 0; i < dim; ++i) {
		for (size_t j = i + 1; j < dim; ++j) {
			const double value = 0.5 * (mat[i * dim + j] + mat[j * dim + i]);
			mat[i * dim + j] = value;
			mat[j * dim + i] = value;
		}
	}
	return false;
}

typedef bool (*h5_group_visit_cb_t)(vlx_t* vlx, hid_t group_handle, const char* group_path, void* user_data);

// Guards against pathological or cyclic (soft/external link) hierarchies. VeloxChem
// files nest a handful of levels; anything deeper is not something we should follow.
#define H5_MAX_GROUP_DEPTH 32

static bool h5_visit_groups_recursive_impl(vlx_t* vlx, hid_t group_handle, const char* group_path, h5_group_visit_cb_t callback, void* user_data, int depth) {
	ASSERT(vlx);
	ASSERT(group_path);
	ASSERT(callback);

	if (!callback(vlx, group_handle, group_path, user_data)) {
		return false;
	}

	if (depth >= H5_MAX_GROUP_DEPTH) {
		MD_LOG_ERROR("HDF5 group nesting exceeds %i levels at '%s', not descending further", H5_MAX_GROUP_DEPTH, group_path);
		return true;
	}

	H5G_info_t info = { 0 };
	if (H5Gget_info(group_handle, &info) < 0) {
		MD_LOG_ERROR("Failed to get group info when traversing HDF5 groups");
		return false;
	}

	char name_buf[256];
	for (hsize_t i = 0; i < info.nlinks; ++i) {
		// H5Gget_objname_by_idx / H5Gget_objtype_by_idx are the deprecated 1.6 API and
		// are compiled out when HDF5 is built without the compatibility layer.
		// H5Lget_name_by_idx returns the length excluding the terminator, or negative.
		const ssize_t name_len = H5Lget_name_by_idx(group_handle, ".", H5_INDEX_NAME, H5_ITER_INC, i, name_buf, sizeof(name_buf), H5P_DEFAULT);
		if (name_len < 0) {
			continue;
		}
		if ((size_t)name_len >= sizeof(name_buf)) {
			// Truncated: opening the shortened name would resolve to the wrong object.
			MD_LOG_ERROR("Skipping HDF5 link under '%s' whose name exceeds %zu characters", group_path, sizeof(name_buf) - 1);
			continue;
		}

		if (h5_get_object_type(group_handle, name_buf) != H5I_GROUP) {
			continue;
		}

		hid_t child_group = H5Gopen(group_handle, name_buf, H5P_DEFAULT);
		if (child_group == H5I_INVALID_HID) {
			continue;
		}

		char child_path[512];
		if (strcmp(group_path, "/") == 0) {
			snprintf(child_path, sizeof(child_path), "/%s", name_buf);
		} else {
			snprintf(child_path, sizeof(child_path), "%s/%s", group_path, name_buf);
		}

		bool result = h5_visit_groups_recursive_impl(vlx, child_group, child_path, callback, user_data, depth + 1);
		H5Gclose(child_group);
		if (!result) {
			return false;
		}
	}

	return true;
}

static bool h5_visit_groups_recursive(vlx_t* vlx, hid_t group_handle, const char* group_path, h5_group_visit_cb_t callback, void* user_data) {
	return h5_visit_groups_recursive_impl(vlx, group_handle, group_path, callback, user_data, 0);
}

static bool h5_read_dataset_data(void* out_data, size_t num_samples, hid_t file_id, hid_t mem_type_id, const char* field_name) {
	ASSERT(out_data);

	if (!h5_link_exists(file_id, field_name)) {
		return false;
	}

	hid_t dataset_id = H5Dopen(file_id, field_name, H5P_DEFAULT);
	if (dataset_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 dataset: '%s'", field_name);
		return false;
	}

	bool result = false;
	hid_t space_id = H5Dget_space(dataset_id);
	if (space_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 space");
		goto done;
	}

	const hssize_t num_points = H5Sget_simple_extent_npoints(space_id);
	if (num_points < 0) {
		MD_LOG_ERROR("Failed to query element count for H5 dataset: '%s'", field_name);
		goto done;
	}

	if ((hsize_t)num_points != (hsize_t)num_samples) {
		MD_LOG_ERROR("Unexpected number of points reading H5 dataset '%s', got %llu, expected %zu",
			field_name, (unsigned long long)num_points, num_samples);
		goto done;
	}

	herr_t status = H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, out_data);

	if (status < 0) {
		MD_LOG_ERROR("Failed to read H5 dataset: '%s'", field_name);
		goto done;
	}

	result = true;
done:
	H5Sclose(space_id);
	H5Dclose(dataset_id);

    return result;
}

// Reads an HDF5 dataset STRAIGHT INTO the storage the attribute table just reserved. No staging
// array in between: the table is where the values are going, so it is where they are read.
//
// Publishes nothing when the dataset is absent or unreadable, and leaves nothing behind when a read
// fails halfway - an absent path is how a consumer learns a file did not carry something, so a
// half filled attribute would be worse than none.
static bool vlx_publish_h5(md_system_t* sys, hid_t handle, const char* dataset, str_t path, str_t label,
						   md_unit_t unit, md_attribute_format_t format) {
	ASSERT(sys);

	const size_t count = md_attribute_element_count(&format);
	if (count == 0 || !h5_check_dataset_exists(handle, dataset)) {
		return false;
	}

	md_attribute_id_t id = vlx_publish(sys, path, label, unit, format, NULL, 0);
	if (id == MD_ATTRIBUTE_INVALID) {
		return false;
	}

	void* dst = md_attributes_data(&sys->attributes, id, format.type);
	if (!dst || !h5_read_dataset_data(dst, count, handle, H5T_NATIVE_DOUBLE, dataset)) {
		md_attributes_remove(&sys->attributes, id);
		return false;
	}
	return true;
}

// rank 1 {count}, one f64 per entry.
static bool vlx_publish_h5_series(md_system_t* sys, hid_t handle, const char* dataset, str_t path, str_t label,
								  md_unit_t unit, size_t count) {
	md_attribute_format_t format = {
		.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)count },
	};
	return vlx_publish_h5(sys, handle, dataset, path, label, unit, format);
}

// rank 2 {rows,cols}, row major with cols fastest - the layout the datasets are stored in, so this
// is a straight read and never a rearrangement.
static bool vlx_publish_h5_matrix(md_system_t* sys, hid_t handle, const char* dataset, str_t path, str_t label,
								  md_unit_t unit, size_t rows, size_t cols) {
	md_attribute_format_t format = {
		.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2, .shape = { (uint32_t)rows, (uint32_t)cols },
	};
	return vlx_publish_h5(sys, handle, dataset, path, label, unit, format);
}

// A single f64 from a scalar dataset, published as rank 0.
static bool vlx_publish_h5_scalar(md_system_t* sys, hid_t handle, const char* dataset, str_t path, str_t label, md_unit_t unit) {
	double value = 0.0;
	if (!h5_check_dataset_exists(handle, dataset) || !h5_read_scalar(&value, handle, H5T_NATIVE_DOUBLE, dataset)) {
		return false;
	}
	vlx_publish_scalar(sys, path, label, unit, value);
	return true;
}

// HDF5 has no native complex type. h5py stores numpy complex arrays as a compound type with two
// floating point members, named 'r' and 'i' by default (configurable via h5py's 'complex_names'),
// laid out as interleaved (re, im) pairs on disk.
//
// This reads such a dataset and splits it into two separate, tightly packed arrays of doubles.
// HDF5 does the de-interleaving for us: a memory datatype that declares only one of the two members
// makes H5Dread gather just that component, so no interleaved temporary buffer is needed.
//
// out_real and out_imag must each hold num_samples doubles. Either may be NULL to skip that
// component. Plain (non-compound) real datasets are also accepted, in which case out_imag is zeroed.
static bool h5_read_complex_dataset_split(double* out_real, double* out_imag, size_t num_samples, hid_t file_id, const char* field_name) {
	if (!h5_link_exists(file_id, field_name)) {
		return false;
	}

	hid_t dataset_id = H5Dopen(file_id, field_name, H5P_DEFAULT);
	if (dataset_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 dataset: '%s'", field_name);
		return false;
	}

	bool  result       = false;
	hid_t space_id     = H5I_INVALID_HID;
	hid_t file_type_id = H5I_INVALID_HID;
	hid_t real_type_id = H5I_INVALID_HID;
	hid_t imag_type_id = H5I_INVALID_HID;
	char* member_name[2] = { NULL, NULL };

	space_id = H5Dget_space(dataset_id);
	if (space_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open H5 space for dataset: '%s'", field_name);
		goto done;
	}

	hsize_t num_points = H5Sget_simple_extent_npoints(space_id);
	if (num_points != (hsize_t)num_samples) {
		MD_LOG_ERROR("Unexpected number of points when reading dataset '%s', got %i, expected %i", field_name, (int)num_points, (int)num_samples);
		goto done;
	}

	file_type_id = H5Dget_type(dataset_id);
	if (file_type_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to get H5 datatype for dataset: '%s'", field_name);
		goto done;
	}

	if (H5Tget_class(file_type_id) != H5T_COMPOUND) {
		// Not stored as complex, treat the data as purely real.
		if (out_real) {
			if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, out_real) < 0) {
				MD_LOG_ERROR("An error occured when reading H5 dataset: '%s'", field_name);
				goto done;
			}
		}
		if (out_imag) {
			MEMSET(out_imag, 0, num_samples * sizeof(double));
		}
		result = true;
		goto done;
	}

	if (H5Tget_nmembers(file_type_id) != 2) {
		MD_LOG_ERROR("Expected 2 members in the compound H5 datatype of dataset '%s'", field_name);
		goto done;
	}

	// The member names depend on how the file was written, so query them instead of assuming 'r'/'i'.
	member_name[0] = H5Tget_member_name(file_type_id, 0);
	member_name[1] = H5Tget_member_name(file_type_id, 1);
	if (!member_name[0] || !member_name[1]) {
		MD_LOG_ERROR("Failed to get the compound member names of H5 dataset: '%s'", field_name);
		goto done;
	}

	// Members are conventionally ordered (real, imaginary), but do not rely on it.
	int real_idx = (member_name[0][0] == 'i' || member_name[0][0] == 'I') ? 1 : 0;
	int imag_idx = 1 - real_idx;

	if (out_real) {
		real_type_id = H5Tcreate(H5T_COMPOUND, sizeof(double));
		if (real_type_id == H5I_INVALID_HID || H5Tinsert(real_type_id, member_name[real_idx], 0, H5T_NATIVE_DOUBLE) < 0) {
			MD_LOG_ERROR("Failed to create a memory datatype for the real part of H5 dataset: '%s'", field_name);
			goto done;
		}
		if (H5Dread(dataset_id, real_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, out_real) < 0) {
			MD_LOG_ERROR("An error occured when reading the real part of H5 dataset: '%s'", field_name);
			goto done;
		}
	}

	if (out_imag) {
		imag_type_id = H5Tcreate(H5T_COMPOUND, sizeof(double));
		if (imag_type_id == H5I_INVALID_HID || H5Tinsert(imag_type_id, member_name[imag_idx], 0, H5T_NATIVE_DOUBLE) < 0) {
			MD_LOG_ERROR("Failed to create a memory datatype for the imaginary part of H5 dataset: '%s'", field_name);
			goto done;
		}
		if (H5Dread(dataset_id, imag_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, out_imag) < 0) {
			MD_LOG_ERROR("An error occured when reading the imaginary part of H5 dataset: '%s'", field_name);
			goto done;
		}
	}

	result = true;
done:
	if (member_name[0]) H5free_memory(member_name[0]);
	if (member_name[1]) H5free_memory(member_name[1]);
	if (real_type_id != H5I_INVALID_HID) H5Tclose(real_type_id);
	if (imag_type_id != H5I_INVALID_HID) H5Tclose(imag_type_id);
	if (file_type_id != H5I_INVALID_HID) H5Tclose(file_type_id);
	if (space_id     != H5I_INVALID_HID) H5Sclose(space_id);
	H5Dclose(dataset_id);

	return result;
}

static bool h5_read_atomic_properties_in_group(vlx_t* vlx, hid_t group_handle, const char* group_path, void* user_data) {
	(void)group_path;
	(void)user_data;

	H5G_info_t info = { 0 };
	if (H5Gget_info(group_handle, &info) < 0) {
		MD_LOG_ERROR("Failed to get group info when reading atomic properties");
		return false;
	}

	char name_buf[256];
	for (hsize_t i = 0; i < info.nlinks; ++i) {
		ssize_t size = H5Gget_objname_by_idx(group_handle, i, name_buf, sizeof(name_buf));
		if (size < 0) {
			continue;
		}
		H5G_obj_t type = H5Gget_objtype_by_idx(group_handle, i);

		// Ensure that the type is a dataset, if not we skip
		if (type != H5G_DATASET) {
			continue;
		}
		hid_t dataset_id = H5Dopen(group_handle, name_buf, H5P_DEFAULT);
		if (dataset_id == H5I_INVALID_HID) {
			continue;
		}

		// The presence of this attribute is what marks a dataset as an atomic property.
		// A dataset that lacks it, or whose attribute is unreadable, is simply skipped
		// rather than aborting the whole traversal.
		char property_label[256] = { 0 };
		if (!h5_read_string_attribute(property_label, sizeof(property_label), dataset_id, "atomic_property")) {
			H5Dclose(dataset_id);
			continue;
		}

		if (property_label[0] == '\0') {
			// Not an error, just use the field name as the property label
			snprintf(property_label, sizeof(property_label), "%s", name_buf);
		}

		// We have a property label, we attempt to read the dataset as an array of doubles with the length of number of atoms
		hid_t space_id = H5Dget_space(dataset_id);
		if (space_id == H5I_INVALID_HID) {
			MD_LOG_ERROR("Failed to get dataspace for dataset '%s'", name_buf);
			goto done;
		}

		int num_dims = H5Sget_simple_extent_ndims(space_id);
		if (num_dims < 0) {
			MD_LOG_ERROR("Failed to get number of dimensions for dataset '%s'", name_buf);
			H5Sclose(space_id);
			goto done;
		}
		if (num_dims < 1 || num_dims > 2) {
			MD_LOG_ERROR("Unsupported rank for atomic property dataset '%s', expected 1 or 2, got %i", name_buf, num_dims);
			H5Sclose(space_id);
			goto done;
		}

		size_t dims[2] = { 0 };
		H5Sget_simple_extent_dims(space_id, (hsize_t*)dims, 0);

		// We expect the inner most dimension to be the number of atoms. Otherwise we skip.
		if (dims[num_dims - 1] != vlx->number_of_atoms) {
			MD_LOG_ERROR("Unexpected size of innermost dimension for atomic property dataset '%s', expected %zu, got %zu", name_buf, vlx->number_of_atoms, dims[num_dims - 1]);
			H5Sclose(space_id);
			continue;
		}

		size_t num_points = H5Sget_simple_extent_npoints(space_id);
		
		H5Sclose(space_id);

		// Published straight into the system's attribute table, with the dataset read into the
		// storage the table just reserved. The dimensions are kept exactly as the dataspace reported
		// them: this reader only accepts a dataset whose INNERMOST dimension is the atom count,
		// which is the attribute convention that the atom axis is the last index axis. So a plain
		// per atom property is rank 1 {N} and one with variants (excited states, spins, whatever the
		// file meant) is rank 2 {S,N}, with scalar values in both cases.
		//
		// The path is built from the DATASET NAME, not from the label: a label is display text that
		// two datasets are free to share, and the path is the property's identity.
		char path_buf[256];
		str_t name = str_from_cstr(name_buf);
		str_t path = vlx_attribute_path(path_buf, sizeof(path_buf), STR_LIT("atom"), name);
		if (str_empty(path)) {
			H5Dclose(dataset_id);
			continue;
		}

		md_attribute_format_t format = {
			.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = (uint32_t)num_dims,
		};
		for (int d = 0; d < num_dims; ++d) {
			format.shape[d] = (uint32_t)dims[d];
		}

		// The label is what the file called it for a human; when it is the dataset name there is
		// nothing for it to add, and an absent label is a valid state.
		str_t label = str_eq(str_from_cstr(property_label), name) ? (str_t){0} : str_copy_cstr(property_label, vlx->arena);

		md_attribute_id_t id = vlx_publish(vlx->sys, path, label, md_unit_none(), format, NULL, 0);
		double* dst = id != MD_ATTRIBUTE_INVALID ? (double*)md_attributes_data(&vlx->sys->attributes, id, MD_ATTRIBUTE_TYPE_F64) : NULL;
		if (!dst) {
			if (id != MD_ATTRIBUTE_INVALID) md_attributes_remove(&vlx->sys->attributes, id);
			H5Dclose(dataset_id);
			continue;
		}

		herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dst);
		if (status < 0) {
			MD_LOG_ERROR("Failed to read data for atomic property dataset '%s'", name_buf);
			md_attributes_remove(&vlx->sys->attributes, id);
			goto done;
		}
		(void)num_points;
	done:
		// The attribute handles are owned and released by h5_read_string_attribute().
		H5Dclose(dataset_id);
	}
	return true;
}

static bool h5_read_density_properties_in_group(vlx_t* vlx, hid_t group_handle, const char* group_path, void* user_data) {
	(void)group_path;
	(void)user_data;

	H5G_info_t info = { 0 };
	if (H5Gget_info(group_handle, &info) < 0) {
		MD_LOG_ERROR("Failed to get group info when reading density properties");
		return false;
	}

	char name_buf[256];
	for (hsize_t i = 0; i < info.nlinks; ++i) {
		ssize_t size = H5Gget_objname_by_idx(group_handle, i, name_buf, sizeof(name_buf));
		if (size < 0) {
			continue;
		}
		H5G_obj_t type = H5Gget_objtype_by_idx(group_handle, i);

		// Ensure that the type is a dataset, if not we skip
		if (type != H5G_DATASET) {
			continue;
		}
		hid_t dataset_id = H5Dopen(group_handle, name_buf, H5P_DEFAULT);
		if (dataset_id == H5I_INVALID_HID) {
			continue;
		}

		// As above: the attribute is the marker, so skip rather than abort.
		char property_label[256] = { 0 };
		if (!h5_read_string_attribute(property_label, sizeof(property_label), dataset_id, "density_property")) {
			H5Dclose(dataset_id);
			continue;
		}

		if (property_label[0] == '\0') {
			snprintf(property_label, sizeof(property_label), "%s", name_buf);
		}

		hid_t space_id = H5Dget_space(dataset_id);
		if (space_id == H5I_INVALID_HID) {
			MD_LOG_ERROR("Failed to get dataspace for dataset '%s'", name_buf);
			goto done;
		}

		int num_dims = H5Sget_simple_extent_ndims(space_id);
		if (num_dims < 0) {
			MD_LOG_ERROR("Failed to get number of dimensions for dataset '%s'", name_buf);
			H5Sclose(space_id);
			goto done;
		}

		if (num_dims != 2) {
			MD_LOG_ERROR("Too many dimensions for density property dataset '%s', expected 2, got %i", name_buf, num_dims);
			H5Sclose(space_id);
			goto done;
		}

		size_t dims[2] = { 0 };
		H5Sget_simple_extent_dims(space_id, (hsize_t*)dims, 0);

		size_t num_aos = vlx_scf_number_of_atomic_orbitals(vlx);
		if (dims[0] != num_aos || dims[1] != num_aos) {
			MD_LOG_ERROR("Unexpected dimensions for density property dataset '%s', expected [%zu x %zu], got [%zu x %zu]", name_buf, num_aos, num_aos, dims[0], dims[1]);
			H5Sclose(space_id);
			continue;
		}

		size_t num_points = H5Sget_simple_extent_npoints(space_id);
		
		H5Sclose(space_id);

		// Construct a unique uint64_t key for this property.
		uint64_t key = md_hash64(name_buf, sizeof(name_buf), 0);

		vlx_density_property_t property = {
			 .label = str_copy_cstr(property_label, vlx->arena),
			 .name = str_copy_cstr(name_buf, vlx->arena),
			 .key = key,
			 .dim[0] = dims[0],
			 .dim[1] = dims[1],
			 .data = NULL,
		};

		md_array_resize(property.data, num_points, vlx->arena);
		herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, property.data);
		if (status < 0) {
			MD_LOG_ERROR("Failed to read data for density property dataset '%s'", name_buf);
			goto done;
		}

		// Symmetry is a load-bearing assumption downstream: the GL/GPU density path packs only the
		// upper triangle (density_matrix_upper_tri_extract_float in md_gto.c) and the lower half is
		// never read. The SCF and transition densities are symmetric by construction -- the latter is
		// explicitly symmetrized -- but density properties are a generic pass-through from the file,
		// so nothing has checked them until here. Report and enforce rather than letting half the
		// matrix be silently dropped.
		//
		// Staged rather than published here, unlike the atomic properties above: a density property
		// is an AO matrix, so it goes through the spherical to Cartesian conversion with the rest of
		// the AO data and comes out a DIFFERENT SIZE. It can only be published once that has run.
		if (dims[0] == dims[1] && dims[0] > 1) {
			vlx_report_and_enforce_symmetry(property.data, dims[0], property_label);
		}

		md_array_push(vlx->density_properties, property, vlx->arena);
		MD_LOG_DEBUG("Read density property '%s' with dimensions [%zu x %zu]", property_label, dims[0], dims[1]);
	done:
		// The attribute handles are owned and released by h5_read_string_attribute().
		H5Dclose(dataset_id);
	}
	return true;
}

static bool h5_read_atomic_properties(vlx_t* vlx, hid_t group_handle) {
	return h5_visit_groups_recursive(vlx, group_handle, "/", h5_read_atomic_properties_in_group, NULL);
}

static bool h5_read_density_properties(vlx_t* vlx, hid_t group_handle) {
	return h5_visit_groups_recursive(vlx, group_handle, "/", h5_read_density_properties_in_group, NULL);
}

// ---------------------------------------------------------------------------
// AO permutation helpers
// Reorder AO-indexed matrices from VeloxChem order into shell order in-place.
// ---------------------------------------------------------------------------

// Permute rows of an AO×MO matrix (num_ao rows, num_mo columns, row-major).
// remap[shell_ao] = vlx_ao  =>  dst_row[shell_ao] = src_row[vlx_ao]
// Permute rows of a [num_ao x num_mo] matrix according to remap and transpose to [num_mo x num_ao].
// On entry  mat is [num_ao][num_mo] in VeloxChem AO order.
// On return mat is [num_mo][num_ao] in shell order — each MO is a contiguous row.
static void ao_permute(double* mat, size_t num_ao, size_t num_mo, const int* remap) {
	md_temp_scope_t temp = md_temp_begin();
	double* tmp = md_temp_alloc_array(temp, double, num_ao * num_mo);
    if (!tmp) {
        MD_LOG_ERROR("Failed to allocate temporary buffer for AO permutation");
		goto done;
    }
	MEMCPY(tmp, mat, sizeof(double) * num_ao * num_mo);
	for (size_t mo = 0; mo < num_mo; mo++) {
		for (size_t ao = 0; ao < num_ao; ao++) {
			mat[mo * num_ao + ao] = tmp[(size_t)remap[ao] * num_mo + mo];
		}
	}
done:
	md_temp_end(temp);
}


// Permute both rows and columns of a square AO×AO matrix (num_ao × num_ao, row-major).
static void ao_permute_square(double* mat, size_t num_ao, const int* remap) {
	md_temp_scope_t temp = md_temp_begin();
	double* tmp = md_temp_alloc_array(temp, double, num_ao * num_ao);
    if (!tmp) {
        MD_LOG_ERROR("Failed to allocate temporary buffer for AO permutation");
		goto done;
    }
	MEMCPY(tmp, mat, sizeof(double) * num_ao * num_ao);
	for (size_t i = 0; i < num_ao; i++) {
		size_t si = (size_t)remap[i];
		for (size_t j = 0; j < num_ao; j++) {
			size_t sj = (size_t)remap[j];
			mat[i * num_ao + j] = tmp[si * num_ao + sj];
		}
	}
done:
	md_temp_end(temp);
}

// Permute columns of a [num_mo x num_ao] matrix according to remap.
static void ao_permute_cols(double* mat, size_t num_mo, size_t num_ao, const int* remap) {
	md_temp_scope_t temp = md_temp_begin();
	double* tmp = md_temp_alloc_array(temp, double, num_mo * num_ao);
	if (!tmp) {
		MD_LOG_ERROR("Failed to allocate temporary buffer for AO permutation");
		goto done;
	}
	MEMCPY(tmp, mat, sizeof(double) * num_mo * num_ao);
	for (size_t mo = 0; mo < num_mo; ++mo) {
		for (size_t ao = 0; ao < num_ao; ++ao) {
			mat[mo * num_ao + ao] = tmp[mo * num_ao + (size_t)remap[ao]];
		}
	}
done:
	md_temp_end(temp);
}

// ---------------------------------------------------------------------------
// Pure/spherical -> Cartesian AO conversion
// ---------------------------------------------------------------------------
// VeloxChem stores AO data in the pure (2l+1) basis. md_gto_basis_t consumes the
// Cartesian ((l+1)(l+2)/2) basis -- see the AO CONVENTION block in md_gto.h.
// Everything AO-indexed is converted once, here, at the end of parsing.
//
// ORDER MATTERS. This must run *after* ao_permute*() has put the matrices in
// shell order, and after every piece of format-internal math that touches the
// AO basis. The Cartesian embedding is rank deficient (the s-type contaminant of
// a d shell is unoccupied), so a converted overlap matrix is singular: anything
// that inverts or Lowdin-orthogonalizes S must have run already.

// [num_mo][n_sph] -> [num_mo][n_cart], reallocated from 'arena'.
static bool vlx_cart_convert_coeff(vlx_2d_data_t* mat, const md_gto_basis_t* basis,
	size_t n_sph, size_t n_cart, md_allocator_i* arena, const char* label)
{
	if (!mat->data) return true;

	if (mat->size[1] != n_sph) {
		MD_LOG_ERROR("%s: expected %zu spherical AO columns, got %zu", label, n_sph, mat->size[1]);
		return false;
	}

	const size_t num_mo = mat->size[0];
	double* dst = (double*)md_alloc(arena, sizeof(double) * num_mo * n_cart);
	if (!dst) {
		MD_LOG_ERROR("%s: failed to allocate Cartesian coefficient matrix", label);
		return false;
	}

	for (size_t mo = 0; mo < num_mo; ++mo) {
		if (md_gto_sph_to_cart_vector(dst + mo * n_cart, mat->data + mo * n_sph, basis) != n_cart) {
			MD_LOG_ERROR("%s: spherical to Cartesian conversion failed for MO %zu", label, mo);
			return false;
		}
	}

	mat->data    = dst;
	mat->size[1] = n_cart;
	return true;
}

// [n_sph][n_sph] -> [n_cart][n_cart], reallocated from 'arena'.
static bool vlx_cart_convert_square(vlx_2d_data_t* mat, const md_gto_basis_t* basis,
	size_t n_sph, size_t n_cart, md_allocator_i* arena, const char* label)
{
	if (!mat->data) return true;

	if (mat->size[0] != n_sph || mat->size[1] != n_sph) {
		MD_LOG_ERROR("%s: expected [%zu x %zu] spherical AO matrix, got [%zu x %zu]",
			label, n_sph, n_sph, mat->size[0], mat->size[1]);
		return false;
	}

	double* dst = (double*)md_alloc(arena, sizeof(double) * n_cart * n_cart);
	if (!dst) {
		MD_LOG_ERROR("%s: failed to allocate Cartesian matrix", label);
		return false;
	}

	if (md_gto_sph_to_cart_matrix(dst, mat->data, basis) != n_cart) {
		MD_LOG_ERROR("%s: spherical to Cartesian matrix conversion failed", label);
		return false;
	}

	mat->data    = dst;
	mat->size[0] = n_cart;
	mat->size[1] = n_cart;
	return true;
}

static bool vlx_convert_ao_data_to_cartesian(vlx_t* vlx) {
	md_temp_scope_t temp = md_temp_begin();
	md_allocator_i* temp_alloc = md_temp_allocator(temp);
	bool result = false;

	md_gto_basis_t basis = {0};
	if (!vlx_gto_basis_extract(&basis, vlx, temp_alloc)) {
		MD_LOG_ERROR("Failed to extract GTO basis for Cartesian AO conversion");
		goto done;
	}

	const size_t n_sph  = md_gto_basis_num_sph_ao(&basis);
	const size_t n_cart = md_gto_basis_num_ao(&basis);
	if (n_sph == 0 || n_cart == 0) {
		MD_LOG_ERROR("Cartesian AO conversion: empty basis (n_sph=%zu n_cart=%zu)", n_sph, n_cart);
		goto done;
	}

	// In the restricted case the beta orbital shares alpha's buffers (see the
	// struct memcpy in the parse path). Detect that so we convert once and
	// re-alias rather than converting the same memory twice.
	const bool beta_aliases_alpha =
		(vlx->scf.beta.coefficients.data != NULL &&
		 vlx->scf.beta.coefficients.data == vlx->scf.alpha.coefficients.data);

	if (!vlx_cart_convert_coeff(&vlx->scf.alpha.coefficients, &basis, n_sph, n_cart, vlx->arena, "Alpha orbital coefficients")) goto done;
	if (!vlx_cart_convert_square(&vlx->scf.alpha.density,     &basis, n_sph, n_cart, vlx->arena, "Alpha density")) goto done;

	if (beta_aliases_alpha) {
		MEMCPY(&vlx->scf.beta.coefficients, &vlx->scf.alpha.coefficients, sizeof(vlx_2d_data_t));
		MEMCPY(&vlx->scf.beta.density,      &vlx->scf.alpha.density,      sizeof(vlx_2d_data_t));
	} else {
		if (!vlx_cart_convert_coeff(&vlx->scf.beta.coefficients, &basis, n_sph, n_cart, vlx->arena, "Beta orbital coefficients")) goto done;
		if (!vlx_cart_convert_square(&vlx->scf.beta.density,     &basis, n_sph, n_cart, vlx->arena, "Beta density")) goto done;
	}

	if (!vlx_cart_convert_square(&vlx->scf.S, &basis, n_sph, n_cart, vlx->arena, "SCF overlap")) goto done;

	// Density properties are AO-basis [N][N] matrices read straight from the file.
	for (size_t i = 0; i < md_array_size(vlx->density_properties); ++i) {
		vlx_density_property_t* prop = &vlx->density_properties[i];
		if (!prop->data) continue;

		vlx_2d_data_t view = { .size = { prop->dim[0], prop->dim[1] }, .data = prop->data };
		if (!vlx_cart_convert_square(&view, &basis, n_sph, n_cart, vlx->arena, "Density property")) goto done;

		prop->data   = view.data;
		prop->dim[0] = view.size[0];
		prop->dim[1] = view.size[1];
	}

	// Derive the AO -> atom map from the shell list, so it cannot drift out of
	// step with the AO ordering the evaluator walks.
	result = true;
done:
	md_temp_end(temp);
	return result;
}

static bool validate_square_matrix_dims(const vlx_2d_data_t* data, const char* label) {
	ASSERT(data);
	ASSERT(label);

	if (data->size[0] == 0 || data->size[1] == 0) {
		MD_LOG_ERROR("%s matrix has invalid dimensions [%zu x %zu]", label, data->size[0], data->size[1]);
		return false;
	}
	if (data->size[0] != data->size[1]) {
		MD_LOG_ERROR("%s matrix must be square, got [%zu x %zu]", label, data->size[0], data->size[1]);
		return false;
	}
	return true;
}

static bool infer_num_mo_from_coeff_dims(size_t* num_mo, const size_t coeff_dim[2], size_t num_ao, const char* label) {
	ASSERT(num_mo);
	ASSERT(coeff_dim);
	ASSERT(label);

	if (coeff_dim[0] == 0 || coeff_dim[1] == 0) {
		MD_LOG_ERROR("%s coefficient matrix has invalid dimensions [%zu x %zu]", label, coeff_dim[0], coeff_dim[1]);
		return false;
	}

	if (coeff_dim[0] == num_ao) {
		*num_mo = coeff_dim[1];
		return true;
	}
	if (coeff_dim[1] == num_ao) {
		*num_mo = coeff_dim[0];
		return true;
	}

	MD_LOG_ERROR("%s coefficient matrix [%zu x %zu] does not contain AO dimension %zu", label, coeff_dim[0], coeff_dim[1], num_ao);
	return false;
}

static bool validate_orbital_canonical_layout(const vlx_orbital_t* orb, size_t num_ao, const char* label) {
	ASSERT(orb);
	ASSERT(label);

	if (orb->coefficients.data) {
		const size_t num_mo = orb->coefficients.size[0];
		if (num_mo == 0 || orb->coefficients.size[1] != num_ao) {
			MD_LOG_ERROR("%s coefficient matrix is not in canonical [MO x AO] layout, got [%zu x %zu], expected [num_mo x %zu]", label, orb->coefficients.size[0], orb->coefficients.size[1], num_ao);
			return false;
		}
		if (orb->energy.data && orb->energy.size != num_mo) {
			MD_LOG_ERROR("%s energy vector length mismatch, expected %zu, got %zu", label, num_mo, orb->energy.size);
			return false;
		}
		if (orb->occupancy.data && orb->occupancy.size != num_mo) {
			MD_LOG_ERROR("%s occupancy vector length mismatch, expected %zu, got %zu", label, num_mo, orb->occupancy.size);
			return false;
		}
	}

	if (orb->density.data) {
		if (!validate_square_matrix_dims(&orb->density, label)) {
			return false;
		}
		if (orb->density.size[0] != num_ao) {
			MD_LOG_ERROR("%s density matrix dimension mismatch, expected %zu, got %zu", label, num_ao, orb->density.size[0]);
			return false;
		}
	}

	return true;
}

static bool normalize_orbital_coefficients(vlx_orbital_t* orb, size_t num_ao, const int* remap, const char* label) {
	ASSERT(orb);
	ASSERT(label);

	if (!orb->coefficients.data) {
		return true;
	}

	const size_t rows = orb->coefficients.size[0];
	const size_t cols = orb->coefficients.size[1];
	if (rows == 0 || cols == 0) {
		MD_LOG_ERROR("%s coefficient matrix has invalid dimensions [%zu x %zu]", label, rows, cols);
		return false;
	}

	if (rows == num_ao) {
		const size_t num_mo = cols;
		if (remap) {
			ao_permute(orb->coefficients.data, num_ao, num_mo, remap);
		}
		orb->coefficients.size[0] = num_mo;
		orb->coefficients.size[1] = num_ao;
		return true;
	}

	if (cols == num_ao) {
		const size_t num_mo = rows;
		if (remap) {
			ao_permute_cols(orb->coefficients.data, num_mo, num_ao, remap);
		}
		return true;
	}

	MD_LOG_ERROR("%s coefficient matrix [%zu x %zu] does not contain AO dimension %zu", label, rows, cols, num_ao);
	return false;
}

static bool validate_scf_canonical_layout(const vlx_t* vlx) {
	ASSERT(vlx);

	if (vlx->scf.S.data) {
		if (!validate_square_matrix_dims(&vlx->scf.S, "SCF overlap")) {
			return false;
		}
	}

	if (vlx->scf.alpha.density.data) {
		if (!validate_square_matrix_dims(&vlx->scf.alpha.density, "Alpha density")) {
			return false;
		}
		if (!validate_orbital_canonical_layout(&vlx->scf.alpha, vlx->scf.alpha.density.size[0], "Alpha orbital")) {
			return false;
		}
	}

	if (vlx->scf.beta.density.data) {
		if (!validate_square_matrix_dims(&vlx->scf.beta.density, "Beta density")) {
			return false;
		}
		if (!validate_orbital_canonical_layout(&vlx->scf.beta, vlx->scf.beta.density.size[0], "Beta orbital")) {
			return false;
		}
	}

	if (vlx->scf.alpha.density.data && vlx->scf.S.data && vlx->scf.alpha.density.size[0] != vlx->scf.S.size[0]) {
		MD_LOG_ERROR("SCF overlap/AO dimension mismatch, alpha density is %zu and overlap is %zu", vlx->scf.alpha.density.size[0], vlx->scf.S.size[0]);
		return false;
	}
	if (vlx->scf.alpha.density.data && vlx->scf.beta.density.data && vlx->scf.alpha.density.size[0] != vlx->scf.beta.density.size[0]) {
		MD_LOG_ERROR("SCF alpha/beta AO dimension mismatch, alpha density is %zu and beta density is %zu", vlx->scf.alpha.density.size[0], vlx->scf.beta.density.size[0]);
		return false;
	}

	return true;
}


// Data extraction procedures
static bool h5_read_scf_data(vlx_t* vlx, hid_t handle) {
	md_system_t* sys = vlx->sys;
	char scf_type[64] = {0};
	if (!h5_read_cstr(scf_type, sizeof(scf_type), handle, "scf_type")) {
		return false;
	}

	if (str_eq_cstr(STR_LIT("restricted"), scf_type)) {
		vlx->scf.type = VLX_SCF_RESTRICTED;
	} else if (str_eq_cstr(STR_LIT("restricted_openshell"), scf_type)) {
		vlx->scf.type = VLX_SCF_RESTRICTED_OPENSHELL;
	} else if (str_eq_cstr(STR_LIT("unrestricted"), scf_type)) {
		vlx->scf.type = VLX_SCF_UNRESTRICTED;
	} else {
		vlx->scf.type = VLX_SCF_UNKNOWN;
		MD_LOG_ERROR("Unrecognized scf type present in h5 scf section: '%s'", scf_type);
		return false;
	}

	if (!h5_read_scalar(&vlx->scf.energy, handle, H5T_NATIVE_DOUBLE, "scf_energy")) {
		return false;
	}

	size_t dim[2];
	h5_read_dataset_dims(dim, 2, handle, "C_alpha");

	// Density dimensions (May differ from dim is always square)
	size_t den_dim[2];
    h5_read_dataset_dims(den_dim, 2, handle, "D_alpha");
	if (!validate_square_matrix_dims(&(vlx_2d_data_t){ .size = {den_dim[0], den_dim[1]}, .data = NULL }, "Alpha density")) {
		return false;
	}

	const size_t num_ao = den_dim[0];
	size_t num_mo = 0;
	if (!infer_num_mo_from_coeff_dims(&num_mo, dim, num_ao, "Alpha coefficient")) {
		return false;
	}

	md_array_resize(vlx->scf.alpha.coefficients.data, dim[0] * dim[1], vlx->arena);
	MEMCPY(vlx->scf.alpha.coefficients.size, dim, sizeof(dim));

	md_array_resize(vlx->scf.alpha.energy.data, num_mo, vlx->arena);
	vlx->scf.alpha.energy.size = num_mo;

	md_array_resize(vlx->scf.alpha.occupancy.data, num_mo, vlx->arena);
	vlx->scf.alpha.occupancy.size = num_mo;

	md_array_resize(vlx->scf.alpha.density.data, den_dim[0] * den_dim[1], vlx->arena);
    MEMCPY(vlx->scf.alpha.density.size, den_dim, sizeof(den_dim));

	// Extract alpha data
	if (!h5_read_dataset_data(vlx->scf.alpha.coefficients.data, md_array_size(vlx->scf.alpha.coefficients.data), handle, H5T_NATIVE_DOUBLE, "C_alpha")) {
		return false;
	}
	if (!h5_read_dataset_data(vlx->scf.alpha.energy.data, md_array_size(vlx->scf.alpha.energy.data), handle, H5T_NATIVE_DOUBLE, "E_alpha")) {
		return false;
	}
	if (!h5_read_dataset_data(vlx->scf.alpha.occupancy.data, md_array_size(vlx->scf.alpha.occupancy.data), handle, H5T_NATIVE_DOUBLE, "occ_alpha")) {
		return false;
	}
    if (!h5_read_dataset_data(vlx->scf.alpha.density.data, md_array_size(vlx->scf.alpha.density.data), handle, H5T_NATIVE_DOUBLE, "D_alpha")) {
        return false;
    }

	if (vlx->scf.type == VLX_SCF_UNRESTRICTED) {
		size_t beta_dim[2];
		h5_read_dataset_dims(beta_dim, 2, handle, "C_beta");
		size_t beta_num_mo = 0;
		if (!infer_num_mo_from_coeff_dims(&beta_num_mo, beta_dim, num_ao, "Beta coefficient")) {
			return false;
		}

		size_t beta_den_dim[2];
		h5_read_dataset_dims(beta_den_dim, 2, handle, "D_beta");
		if (!validate_square_matrix_dims(&(vlx_2d_data_t){ .size = {beta_den_dim[0], beta_den_dim[1]}, .data = NULL }, "Beta density")) {
			return false;
		}
		if (beta_den_dim[0] != num_ao) {
			MD_LOG_ERROR("Alpha/Beta AO dimension mismatch, alpha density is %zu and beta density is %zu", num_ao, beta_den_dim[0]);
			return false;
		}

		md_array_resize(vlx->scf.beta.coefficients.data, beta_dim[0] * beta_dim[1], vlx->arena);
		MEMCPY(vlx->scf.beta.coefficients.size, beta_dim, sizeof(beta_dim));

		md_array_resize(vlx->scf.beta.energy.data, beta_num_mo, vlx->arena);
		vlx->scf.beta.energy.size = beta_num_mo;

		md_array_resize(vlx->scf.beta.occupancy.data, beta_num_mo, vlx->arena);
		vlx->scf.beta.occupancy.size = beta_num_mo;

		md_array_resize(vlx->scf.beta.density.data, beta_den_dim[0] * beta_den_dim[1], vlx->arena);
			MEMCPY(vlx->scf.beta.density.size, beta_den_dim, sizeof(beta_den_dim));

		// Extract beta data
		if (!h5_read_dataset_data(vlx->scf.beta.coefficients.data, md_array_size(vlx->scf.beta.coefficients.data), handle, H5T_NATIVE_DOUBLE, "C_beta")) {
			return false;
		}
		if (!h5_read_dataset_data(vlx->scf.beta.energy.data, md_array_size(vlx->scf.beta.energy.data), handle, H5T_NATIVE_DOUBLE, "E_beta")) {
			return false;
		}
		if (!h5_read_dataset_data(vlx->scf.beta.occupancy.data, md_array_size(vlx->scf.beta.occupancy.data), handle, H5T_NATIVE_DOUBLE, "occ_beta")) {
			return false;
		}
        if (!h5_read_dataset_data(vlx->scf.beta.density.data, md_array_size(vlx->scf.beta.density.data), handle, H5T_NATIVE_DOUBLE, "D_beta")) {
            return false;
        }
	} else {
		// Shallow copy fields from Alpha
		MEMCPY(&vlx->scf.beta, &vlx->scf.alpha, sizeof(vlx_orbital_t));
		if (vlx->scf.type == VLX_SCF_RESTRICTED_OPENSHELL) {
			vlx->scf.beta.occupancy.data = 0;
			md_array_resize(vlx->scf.beta.occupancy.data, vlx->scf.beta.occupancy.size, vlx->arena);
			if (!h5_read_dataset_data(vlx->scf.beta.occupancy.data, md_array_size(vlx->scf.beta.occupancy.data), handle, H5T_NATIVE_DOUBLE, "occ_beta")) {
				return false;
			}
		}
	}

	// S matrix is overlap (notice dimension is the same as D)
	md_array_resize(vlx->scf.S.data, den_dim[0] * den_dim[1], vlx->arena);
    MEMCPY(vlx->scf.S.size, den_dim, sizeof(den_dim));

	if (!h5_read_dataset_data(vlx->scf.S.data, md_array_size(vlx->scf.S.data), handle, H5T_NATIVE_DOUBLE, "S")) {
		return false;
	}

	// The ground state dipole moment is not present in all versions
	if (!h5_read_dataset_data(&vlx->scf.ground_state_dipole_moment, 3, handle, H5T_NATIVE_DOUBLE, "dipole_moment")) {
		//return false;
	}

	// NOTE: H5Lexists returns htri_t -- negative on error, which is truthy. Test
	// explicitly, or a failed lookup is taken as "present" and the read proceeds.
	// The convergence history: one value per iteration, per quantity, straight into the table. Two
	// layouts in the wild - a group per iteration holding five scalars, or five flat datasets - and
	// both end up as the same five {I} attributes.
	static const struct { const char* group_field; const char* flat_field; const char* path; const char* label; bool hartree; } history[] = {
		{ "energy",        "scf_history_energy",        "vlx/scf/history/energy",        "Energy",              true  },
		{ "diff_energy",   "scf_history_diff_energy",   "vlx/scf/history/energy_diff",   "Energy Difference",   true  },
		{ "diff_density",  "scf_history_diff_density",  "vlx/scf/history/density_diff",  "Density Difference",  false },
		{ "gradient_norm", "scf_history_gradient_norm", "vlx/scf/history/gradient_norm", "Gradient Norm",       false },
		{ "max_gradient",  "scf_history_max_gradient",  "vlx/scf/history/max_gradient",  "Max Gradient",        false },
	};

	// NOTE: H5Lexists returns htri_t -- negative on error, which is truthy. Test
	// explicitly, or a failed lookup is taken as "present" and the read proceeds.
	if (h5_link_exists(handle, "scf_history")) {
		// One group per iteration, labelled '0' ... 'N', each holding the five scalars named above.
		hid_t scf_history = H5Gopen(handle, "scf_history", H5P_DEFAULT);
		if (scf_history < 0) {
			return false;
		}

		hsize_t h5_num_links = 0;
		if (H5Gget_num_objs(scf_history, &h5_num_links) < 0) {
			H5Gclose(scf_history);
			return false;
		}
		const size_t num_links = (size_t)h5_num_links;

		// Counted before anything is created, because the attribute's length is its shape and a
		// series cannot be grown after the fact. Groups are counted rather than assumed present, so
		// a file that skips one publishes a shorter history rather than a run of zeros.
		size_t num_iter = 0;
		for (size_t i = 0; i < num_links; ++i) {
			char name_buf[64];
			snprintf(name_buf, sizeof(name_buf), "%zu", i);
			if (H5Lexists(scf_history, name_buf, H5P_DEFAULT) > 0) {
				num_iter += 1;
			}
		}

		if (num_iter > 0) {
			md_attribute_format_t format = {
				.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)num_iter },
			};
			double* dst[ARRAY_SIZE(history)] = {0};
			for (size_t f = 0; f < ARRAY_SIZE(history); ++f) {
				md_attribute_id_t id = vlx_publish(sys, str_from_cstr(history[f].path), str_from_cstr(history[f].label),
												   history[f].hartree ? vlx_unit_hartree() : md_unit_none(), format, NULL, 0);
				dst[f] = id != MD_ATTRIBUTE_INVALID ? (double*)md_attributes_data(&sys->attributes, id, MD_ATTRIBUTE_TYPE_F64) : NULL;
				if (!dst[f]) {
					H5Gclose(scf_history);
					return false;
				}
			}

			size_t iter = 0;
			for (size_t i = 0; i < num_links && iter < num_iter; ++i) {
				char name_buf[64];
				snprintf(name_buf, sizeof(name_buf), "%zu", i);
				if (H5Lexists(scf_history, name_buf, H5P_DEFAULT) <= 0) {
					continue;
				}
				hid_t iter_group = H5Gopen(scf_history, name_buf, H5P_DEFAULT);
				if (iter_group < 0) {
					H5Gclose(scf_history);
					return false;
				}
				for (size_t f = 0; f < ARRAY_SIZE(history); ++f) {
					if (!h5_read_dataset_data(dst[f] + iter, 1, iter_group, H5T_NATIVE_DOUBLE, history[f].group_field)) {
						H5Gclose(iter_group);
						H5Gclose(scf_history);
						return false;
					}
				}
				H5Gclose(iter_group);
				iter += 1;
			}
		}
		H5Gclose(scf_history);
	} else if (h5_link_exists(handle, "scf_history_energy")) {
		size_t scf_hist_len = 0;
		if (!h5_read_dataset_dims(&scf_hist_len, 1, handle, "scf_history_energy")) {
			return false;
		}

		for (size_t f = 0; f < ARRAY_SIZE(history); ++f) {
			if (!vlx_publish_h5_series(sys, handle, history[f].flat_field, str_from_cstr(history[f].path), str_from_cstr(history[f].label),
									   history[f].hartree ? vlx_unit_hartree() : md_unit_none(), scf_hist_len)) {
				return false;
			}
		}
	}

	// WHICH SCF this was. A consumer can guess from whether the two spin channels share their
	// coefficients and their occupations, and that guess is right until it meets a file where one of
	// the two is missing. The reader knows; this is it saying so, on the same terms as the response
	// and optimisation types. Text rather than an integer so nothing stored depends on this file's
	// internal enum ordering.
	vlx_publish_str(sys, STR_LIT("vlx/scf/type"), STR_LIT("SCF Type"), vlx_scf_type_str(vlx->scf.type));

	return true;
}

static bool h5_read_optional_1d_data(vlx_1d_data_t* out_data, hid_t handle, const char* field_name, md_allocator_i* arena) {
	ASSERT(out_data);
	ASSERT(field_name);
	ASSERT(arena);

	if (!h5_check_dataset_exists(handle, field_name)) {
		return true;
	}

	size_t dim[2] = {0};
	int num_dims = h5_read_dataset_dims(dim, 2, handle, field_name);
	if (num_dims <= 0) {
		MD_LOG_ERROR("Invalid dimensions in H5 vector dataset '%s'", field_name);
		return false;
	}

	size_t sample_count = 1;
	for (int dim_idx = 0; dim_idx < num_dims; ++dim_idx) {
		sample_count *= dim[dim_idx];
	}

	if (sample_count == 0) {
		MD_LOG_ERROR("Empty H5 vector dataset '%s'", field_name);
		return false;
	}

	md_array_resize(out_data->data, sample_count, arena);
	MEMSET(out_data->data, 0, md_array_bytes(out_data->data));
	out_data->size = sample_count;

	if (!h5_read_dataset_data(out_data->data, sample_count, handle, H5T_NATIVE_DOUBLE, field_name)) {
		return false;
	}

	return true;
}

// The occupied/virtual split the response solution vectors are indexed by.
//
// VeloxChem writes num_core/num_valence/num_virtual only for some calculations - none of the files
// in test_data carries them - and without that split a solution vector is an undifferentiated run
// of amplitudes that nothing can be reconstructed from. So when the file is silent, derive it: an
// ordinary valence excitation spans every occupied orbital and every virtual one, which the SCF
// occupations already say.
//
// Derived, never assumed: the split has to reproduce the solution vector's own length (amp_count,
// or twice it when the vector carries both X and Y), and it is adopted only when it does. That
// check is what makes this safe for a core excitation, where the file DOES name num_core and the
// derived valence split would be wrong - there, the stored values are kept and this does nothing.
static void vlx_rsp_infer_occupied_virtual_split(vlx_t* vlx) {
	ASSERT(vlx);

	if (vlx->rsp.num_core > 0 || vlx->rsp.num_valence > 0) {
		return;     // the file said so
	}
	const size_t vec_len = vlx->rsp.solution_matrix.size[1];
	if (!vlx->rsp.solution_matrix.data || vec_len == 0) {
		return;
	}

	const double* occ    = vlx->scf.alpha.occupancy.data;
	const size_t  num_mo = vlx->scf.alpha.occupancy.size;
	if (!occ || num_mo == 0) {
		MD_LOG_ERROR("Response data has no occupied/virtual split and no SCF occupations to derive one from");
		return;
	}

	size_t nocc = 0;
	for (size_t i = 0; i < num_mo; ++i) {
		if (occ[i] > 0.0) nocc += 1;
	}
	const size_t nvir = num_mo - nocc;

	if (nocc == 0 || nvir == 0) {
		MD_LOG_ERROR("Cannot derive an occupied/virtual split: %zu of %zu orbitals are occupied", nocc, num_mo);
		return;
	}

	const size_t amp_count = nocc * nvir;
	if (vec_len != amp_count && vec_len != 2 * amp_count) {
		MD_LOG_ERROR("Derived %zu occupied x %zu virtual orbitals, which does not explain a solution vector of %zu values (expected %zu or %zu)",
			nocc, nvir, vec_len, amp_count, 2 * amp_count);
		return;
	}

	vlx->rsp.num_valence = nocc;
	vlx->rsp.num_virtual = nvir;
	MD_LOG_DEBUG("Derived the response occupied/virtual split from the SCF occupations: %zu x %zu", nocc, nvir);
}

static bool h5_read_rsp_data(vlx_t* vlx, hid_t handle) {
	md_system_t* sys = vlx->sys;

	h5_read_scalar(&vlx->rsp.number_of_frequencies, handle, H5T_NATIVE_HSIZE, "number_of_states");
	if (vlx->rsp.number_of_frequencies > 0) {
		// Standard Linear Response data
		vlx->rsp.type = VLX_RSP_LINEAR;

		if (h5_check_dataset_exists(handle, "num_core")) {
			h5_read_scalar(&vlx->rsp.num_core, handle, H5T_NATIVE_INT64, "num_core");
		}

		if (h5_check_dataset_exists(handle, "num_valence")) {
			h5_read_scalar(&vlx->rsp.num_valence, handle, H5T_NATIVE_INT64, "num_valence");
		}

		if (h5_check_dataset_exists(handle, "num_virtual")) {
			h5_read_scalar(&vlx->rsp.num_virtual, handle, H5T_NATIVE_INT64, "num_virtual");
		}

		// The excitation energies. What the frequency axis MEANS depends on the response type - these
		// for a linear response, a sampled grid for a complex polarisation propagator - which is why
		// it is one path filled from a different dataset per branch below.
		vlx_publish_h5_series(sys, handle, "eigenvalues", STR_LIT("vlx/rsp/frequency"), STR_LIT("Response Frequency"),
							  vlx_unit_hartree(), vlx->rsp.number_of_frequencies);

		// The response eigenvectors. STAGED rather than published here: the occupied/virtual split
		// they are indexed by is derived below from the SCF occupations, and they are republished
		// together in vlx_publish_whole_file_attributes once it is known.
		if (h5_check_dataset_exists(handle, "S1")) {
			size_t dims[2] = { vlx->rsp.number_of_frequencies, 0 };
			h5_read_dataset_dims(&dims[1], 1, handle, "S1");

			size_t len = dims[0] * dims[1];
			if (len > 0) {
				md_array_resize(vlx->rsp.solution_matrix.data, len, vlx->arena);
				vlx->rsp.solution_matrix.size[0] = dims[0];
				vlx->rsp.solution_matrix.size[1] = dims[1];
				char field_name[16];
				for (size_t state_idx = 0; state_idx < vlx->rsp.number_of_frequencies; ++state_idx) {
					snprintf(field_name, sizeof(field_name), "S%zu", state_idx + 1);
					double* dst = vlx->rsp.solution_matrix.data + state_idx * vlx->rsp.solution_matrix.size[1];
					if (!h5_read_dataset_data(dst, vlx->rsp.solution_matrix.size[1], handle, H5T_NATIVE_DOUBLE, field_name)) {
						return false;
					}
				}
			}
		} else if (h5_check_dataset_exists(handle, "full_solutions_matrix")) {
			size_t dims[2] = { 0 };
			h5_read_dataset_dims(dims, 2, handle, "full_solutions_matrix");
			size_t len = dims[0] * dims[1];
			if (len > 0) {
				md_array_resize(vlx->rsp.solution_matrix.data, len, vlx->arena);
				vlx->rsp.solution_matrix.size[0] = dims[0];
				vlx->rsp.solution_matrix.size[1] = dims[1];
				if (!h5_read_dataset_data(vlx->rsp.solution_matrix.data, len, handle, H5T_NATIVE_DOUBLE, "full_solutions_matrix")) {
					return false;
				}
			}
		}
	}

	// After the solution matrix is in hand, so its length is available to check a derived split
	// against, and after the SCF block was read - which vlx_read_h5_file guarantees.
	vlx_rsp_infer_occupied_virtual_split(vlx);

	// P and F for a RIXS run. Locals: nothing after this function needs them, because everything
	// they size is published here.
	size_t num_photons = 0;
	size_t num_final_states = 0;

	if (vlx->rsp.type == VLX_RSP_UNKNOWN) {
		// No standard response data, check for other types of response data by looking for type field
		if (h5_check_dataset_exists(handle, "rsp_type")) {
			char type_buf[32] = { 0 };
			h5_read_cstr(type_buf, sizeof(type_buf), handle, "rsp_type");
			if (strncmp(type_buf, "cpp", sizeof(type_buf)) == 0) {
				vlx->rsp.type = VLX_RSP_CPP;
			} else if (strncmp(type_buf, "c6", sizeof(type_buf)) == 0) {
				vlx->rsp.type = VLX_RSP_C6;
			} else if (strncmp(type_buf, "tpa_transition", sizeof(type_buf)) == 0) {
				vlx->rsp.type = VLX_RSP_TPA_TRANSITION;
			} else if (strncmp(type_buf, "tpa", sizeof(type_buf)) == 0) {
				vlx->rsp.type = VLX_RSP_TPA;
			} else if (strncmp(type_buf, "rixs", sizeof(type_buf)) == 0) {
				vlx->rsp.type = VLX_RSP_RIXS;
			}
		}

		if (vlx->rsp.type == VLX_RSP_C6) {
			// A homomolecular C6 value in a.u., and the whole of what this kind of run produces.
			double c6 = 0.0;
			if (h5_check_dataset_exists(handle, "c6")) {
				if (!h5_read_dataset_data(&c6, 1, handle, H5T_NATIVE_DOUBLE, "c6")) {
					MD_LOG_ERROR("Could not read c6 dataset");
					return false;
				}
				vlx_publish_scalar(sys, STR_LIT("vlx/rsp/c6"), STR_LIT("C6 Coefficient"), md_unit_none(), c6);
			}
		} else if (vlx->rsp.type == VLX_RSP_CPP || vlx->rsp.type == VLX_RSP_TPA) {
			size_t dim;
			if (h5_read_dataset_dims(&dim, 1, handle, "frequencies")) {
				vlx->rsp.number_of_frequencies = dim;
				if (!vlx_publish_h5_series(sys, handle, "frequencies", STR_LIT("vlx/rsp/frequency"), STR_LIT("Response Frequency"), vlx_unit_hartree(), dim)) {
					return false;
				}
			}
		} else if (vlx->rsp.type == VLX_RSP_TPA_TRANSITION) {
			size_t dim;
			if (h5_read_dataset_dims(&dim, 1, handle, "photon_energies")) {
				vlx->rsp.number_of_frequencies = dim;
				if (!vlx_publish_h5_series(sys, handle, "photon_energies", STR_LIT("vlx/rsp/frequency"), STR_LIT("Response Frequency"), vlx_unit_hartree(), dim)) {
					return false;
				}
			}
		} else if (vlx->rsp.type == VLX_RSP_RIXS) {
			// The core-excited states are the frequency axis for a RIXS run, and they are also the
			// XAS side panel's own axis - so the same values land on both paths.
			size_t dim;
			if (h5_read_dataset_dims(&dim, 1, handle, "core_eigenvalues")) {
				vlx->rsp.number_of_frequencies = dim;
				if (!vlx_publish_h5_series(sys, handle, "core_eigenvalues", STR_LIT("vlx/rsp/frequency"), STR_LIT("Response Frequency"), vlx_unit_hartree(), dim)) {
					return false;
				}
				vlx_publish_h5_series(sys, handle, "core_eigenvalues", STR_LIT("vlx/rsp/rixs/core_energy"), STR_LIT("Core Energy"), vlx_unit_hartree(), dim);
			}

			if (h5_read_dataset_dims(&dim, 1, handle, "photon_energies")) {
				num_photons = dim;
			}
		}
	}

	// WHICH response calculation this was. Not derivable from the columns: a linear response and a
	// two-photon transition run both publish peaks over the same frequency axis, and telling them
	// apart by which optional sibling happens to be present is a guess that a file carrying partial
	// data gets wrong. The reader knows which it read; this is it saying so. Text rather than an
	// integer so nothing stored depends on this file's internal enum ordering.
	vlx_publish_str(sys, STR_LIT("vlx/rsp/type"), STR_LIT("Response Type"), vlx_rsp_type_str(vlx->rsp.type));

	if (vlx->rsp.number_of_frequencies > 0) {
		const size_t num_freqs = vlx->rsp.number_of_frequencies;

		// Transition dipoles are STAGED, unlike everything else here. A dipole is published as a
		// group - a vector AND the origin it is drawn from - and the centre of charge that anchors
		// it cannot be computed until the whole file has been read. Half a group is not a dipole
		// anyone can draw, so they are published together in vlx_publish_whole_file_attributes.
		const size_t num_dipole_points = num_freqs * 3;
		if (h5_check_dataset_exists(handle, "electric_transition_dipoles")) {
			md_array_resize(vlx->rsp.electric_transition_dipoles, num_freqs, vlx->arena);
			MEMSET(vlx->rsp.electric_transition_dipoles, 0, md_array_bytes(vlx->rsp.electric_transition_dipoles));
			if (!h5_read_dataset_data(vlx->rsp.electric_transition_dipoles, num_dipole_points, handle, H5T_NATIVE_DOUBLE, "electric_transition_dipoles")) {
				md_array_free(vlx->rsp.electric_transition_dipoles, vlx->arena);
				vlx->rsp.electric_transition_dipoles = NULL;
			}
		}

		if (h5_check_dataset_exists(handle, "magnetic_transition_dipoles")) {
			md_array_resize(vlx->rsp.magnetic_transition_dipoles, num_freqs, vlx->arena);
			MEMSET(vlx->rsp.magnetic_transition_dipoles, 0, md_array_bytes(vlx->rsp.magnetic_transition_dipoles));
			if (!h5_read_dataset_data(vlx->rsp.magnetic_transition_dipoles, num_dipole_points, handle, H5T_NATIVE_DOUBLE, "magnetic_transition_dipoles")) {
				md_array_free(vlx->rsp.magnetic_transition_dipoles, vlx->arena);
				vlx->rsp.magnetic_transition_dipoles = NULL;
			}
		}

		if (h5_check_dataset_exists(handle, "velocity_transition_dipoles")) {
			md_array_resize(vlx->rsp.velocity_transition_dipoles, num_freqs, vlx->arena);
			MEMSET(vlx->rsp.velocity_transition_dipoles, 0, md_array_bytes(vlx->rsp.velocity_transition_dipoles));
			if (!h5_read_dataset_data(vlx->rsp.velocity_transition_dipoles, num_dipole_points, handle, H5T_NATIVE_DOUBLE, "velocity_transition_dipoles")) {
				md_array_free(vlx->rsp.velocity_transition_dipoles, vlx->arena);
				vlx->rsp.velocity_transition_dipoles = NULL;
			}
		}

		// Two photon absorption, over the same frequency axis.
		vlx_publish_h5_series(sys, handle, "tpa_strengths/circular", STR_LIT("vlx/rsp/tpa/circular"), STR_LIT("Circular Polarisation"), md_unit_none(), num_freqs);
		vlx_publish_h5_series(sys, handle, "tpa_strengths/linear",   STR_LIT("vlx/rsp/tpa/linear"),   STR_LIT("Linear Polarisation"),   md_unit_none(), num_freqs);

		if (vlx->rsp.type == VLX_RSP_RIXS) {
			// RIXS involves three independent dimensions:
			//   C = number of core-excited (intermediate) states  -> rsp.number_of_frequencies
			//   F = number of final (valence-excited) states       -> num_final_states
			//   P = number of incoming photon energies            -> num_photons
			// The core states are summed over coherently inside the scattering amplitude and never
			// appear as an output dimension, so F is completely unrelated to C. The 2D datasets are
			// stored row-major as [F][P]. Derive F from a representative dataset rather than assuming.
			static const char* rixs_2d_fields[] = { "cross_sections", "energy_losses", "emission_energies" };
			for (size_t i = 0; i < ARRAY_SIZE(rixs_2d_fields); ++i) {
				size_t dim[2] = { 0 };
				if (h5_read_dataset_dims(dim, (int)ARRAY_SIZE(dim), handle, rixs_2d_fields[i]) == 2 && dim[0] > 0 && dim[1] > 0) {
					num_final_states = dim[0];

					if (num_photons == 0) {
						// 'photon_energies' was missing or unreadable, recover P from the column count.
						num_photons = dim[1];
					} else if (dim[1] != num_photons) {
						MD_LOG_ERROR("RIXS: dataset '%s' has %i columns, expected %i incoming photon energies",
							rixs_2d_fields[i], (int)dim[1], (int)num_photons);
						num_final_states = 0;
					}
					break;
				}
			}

			vlx_publish_h5_matrix(sys, handle, "cross_sections",    STR_LIT("vlx/rsp/rixs/cross_section"),    STR_LIT("Cross Section"),    md_unit_none(),     num_final_states, num_photons);
			vlx_publish_h5_matrix(sys, handle, "emission_energies", STR_LIT("vlx/rsp/rixs/emission_energy"), STR_LIT("Emission Energy"), vlx_unit_hartree(), num_final_states, num_photons);
			vlx_publish_h5_matrix(sys, handle, "energy_losses",     STR_LIT("vlx/rsp/rixs/energy_loss"),     STR_LIT("Energy Loss"),     vlx_unit_hartree(), num_final_states, num_photons);

			vlx_publish_h5_series(sys, handle, "core_osc_strengths",     STR_LIT("vlx/rsp/rixs/core_oscillator_strength"), STR_LIT("Core Oscillator Strength"), md_unit_none(), num_freqs);
			vlx_publish_h5_series(sys, handle, "elastic_cross_sections", STR_LIT("vlx/rsp/rixs/elastic_cross_section"),    STR_LIT("Elastic Cross Section"),    md_unit_none(), num_photons);
			vlx_publish_h5_series(sys, handle, "photon_energies",        STR_LIT("vlx/rsp/rixs/photon_energy"),            STR_LIT("Photon Energy"),           vlx_unit_hartree(), num_photons);

			if (num_freqs > 0) {
				vlx_publish_h5_scalar(sys, handle, "gamma_fwhm_ev", STR_LIT("vlx/rsp/rixs/gamma_fwhm"), STR_LIT("Core-hole Lifetime Broadening"), md_unit_electronvolt());
			}

			if (h5_check_dataset_exists(handle, "scattering_amplitudes")) {
				// Complex, shaped [F][P][3][3]: the Cartesian scattering amplitude tensor per final
				// state and photon energy. Published as two real attributes rather than one
				// interleaved buffer, because a value has one type and a consumer wanting the
				// modulus should not have to know the storage convention to get it.
				// Derive the element count from the dataset itself rather than assuming a rank.
				size_t dim[4] = {0};
				int ndim = h5_read_dataset_dims(dim, (int)ARRAY_SIZE(dim), handle, "scattering_amplitudes");

				// h5_read_dataset_dims reports the rank of the dataset, which may exceed the number of
				// entries it actually wrote, so clamp before iterating.
				if (ndim > (int)ARRAY_SIZE(dim)) {
					MD_LOG_ERROR("Unexpected rank of H5 dataset 'scattering_amplitudes'");
					ndim = 0;
				}

				size_t num_elem = (ndim > 0) ? 1 : 0;
				for (int i = 0; i < ndim; ++i) {
					num_elem *= dim[i];
				}

				// Sanity check against the dimensions derived above: [F][P][3][3] == F * P * 9.
				const size_t num_2d_elem = num_final_states * num_photons;
				if (num_elem > 0 && num_2d_elem > 0 && num_elem != num_2d_elem * 9) {
					MD_LOG_ERROR("RIXS: 'scattering_amplitudes' holds %i elements, expected %i (%i final states x %i photon energies x 3 x 3)",
						(int)num_elem, (int)(num_2d_elem * 9), (int)num_final_states, (int)num_photons);
				}

				if (num_elem > 0 && ndim == 4) {
					md_attribute_format_t format = {
						.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 4,
						.shape = { (uint32_t)dim[0], (uint32_t)dim[1], (uint32_t)dim[2], (uint32_t)dim[3] },
					};
					md_attribute_id_t re_id = vlx_publish(sys, STR_LIT("vlx/rsp/rixs/scattering_amplitude_re"), STR_LIT("Scattering Amplitude (Re)"), md_unit_none(), format, NULL, 0);
					md_attribute_id_t im_id = vlx_publish(sys, STR_LIT("vlx/rsp/rixs/scattering_amplitude_im"), STR_LIT("Scattering Amplitude (Im)"), md_unit_none(), format, NULL, 0);
					double* re = re_id != MD_ATTRIBUTE_INVALID ? (double*)md_attributes_data(&sys->attributes, re_id, MD_ATTRIBUTE_TYPE_F64) : NULL;
					double* im = im_id != MD_ATTRIBUTE_INVALID ? (double*)md_attributes_data(&sys->attributes, im_id, MD_ATTRIBUTE_TYPE_F64) : NULL;

					// Splits the interleaved complex data on disk into the two attributes.
					if (!re || !im || !h5_read_complex_dataset_split(re, im, num_elem, handle, "scattering_amplitudes")) {
						if (re_id != MD_ATTRIBUTE_INVALID) md_attributes_remove(&sys->attributes, re_id);
						if (im_id != MD_ATTRIBUTE_INVALID) md_attributes_remove(&sys->attributes, im_id);
					}
				}
			}
		} else {
			vlx_publish_h5_series(sys, handle, "cross_sections", STR_LIT("vlx/rsp/tpa/cross_section"), STR_LIT("Cross Section"), md_unit_none(), num_freqs);
		}

		// The complex polarisation propagator outputs, sampled over the same frequency axis.
		vlx_publish_h5_series(sys, handle, "sigma",            STR_LIT("vlx/rsp/cpp/sigma"),            STR_LIT("Absorption Cross Section"), md_unit_none(), num_freqs);
		vlx_publish_h5_series(sys, handle, "optical-rotation", STR_LIT("vlx/rsp/cpp/optical_rotation"), STR_LIT("Optical Rotation"),         md_unit_none(), num_freqs);
		vlx_publish_h5_series(sys, handle, "delta-epsilon",    STR_LIT("vlx/rsp/cpp/delta_epsilon"),    STR_LIT(u8"Δε"),                     md_unit_none(), num_freqs);

		// One value per EXCITED STATE, which only a linear response has: the other response types
		// sample a frequency grid rather than resolving states, so a peak list read from one of them
		// would be indexed by something it does not have.
		if (vlx->rsp.type == VLX_RSP_LINEAR) {
			vlx_publish_h5_series(sys, handle, "oscillator_strengths", STR_LIT("vlx/rsp/oscillator_strength"), STR_LIT("Oscillator Strength"), md_unit_none(), num_freqs);
			vlx_publish_h5_series(sys, handle, "rotatory_strengths",   STR_LIT("vlx/rsp/rotatory_strength"),   STR_LIT("Rotatory Strength"),   md_unit_none(), num_freqs);
		}
	}

	return true;
}

static bool h5_read_vib_data(vlx_t* vlx, hid_t handle) {
	size_t number_of_modes = 0;

	// Attempt to read number_of_modes (Available in new format)
	if (!h5_read_scalar(&number_of_modes, handle, H5T_NATIVE_HSIZE, "number_of_modes")) {
		// Fallback (Old format, read force_constant dims to get number of modes)
		size_t dim[2];
		int num_dim = h5_read_dataset_dims(dim, 2, handle, "force_constants");
		if (num_dim <= 0) {
			return false;
		}
		// This is a fix because the input data in one version is supplied as a 2D object
		number_of_modes = (num_dim == 1) ? dim[0] : dim[1];
	}

	if (number_of_modes == 0) {
		return false;
	}

	md_system_t* sys = vlx->sys;

	// Read straight into the attribute table. There is nothing a later step has to look at here -
	// a vibrational spectrum is not indexed by anything the AO pipeline touches - so nothing is
	// staged and the reader IS the publisher.
	if (h5_check_dataset_exists(handle, "force_constants") &&
		!vlx_publish_h5_series(sys, handle, "force_constants", STR_LIT("vlx/vib/force_constant"), STR_LIT("Force Constant"), md_unit_none(), number_of_modes)) {
		return false;
	}

	if (h5_check_dataset_exists(handle, "ir_intensities") &&
		!vlx_publish_h5_series(sys, handle, "ir_intensities", STR_LIT("vlx/vib/ir_intensity"), STR_LIT("IR Intensity"), vlx_unit_km_per_mol(), number_of_modes)) {
		return false;
	}

	if (h5_check_dataset_exists(handle, "vib_frequencies") &&
		!vlx_publish_h5_series(sys, handle, "vib_frequencies", STR_LIT("vlx/vib/frequency"), STR_LIT("Frequency"), vlx_unit_wavenumber(), number_of_modes)) {
		return false;
	}

	if (h5_check_dataset_exists(handle, "reduced_masses") &&
		!vlx_publish_h5_series(sys, handle, "reduced_masses", STR_LIT("vlx/vib/reduced_mass"), STR_LIT("Reduced Mass"), vlx_unit_amu(), number_of_modes)) {
		return false;
	}

	// The displacements, {M,N} of 3 component values: the mode axis leads so one mode's
	// displacements are contiguous, which is the case the ATTRIBUTES note in md_system.h uses as its
	// example of why an atom axis is not always shape[0].
	//
	// Two layouts in the wild - one [M][N][3] dataset, or a group of per mode datasets - and both
	// are read into the same reserved storage, one plane at a time in the group case.
	if (vlx->number_of_atoms == 0) {
		MD_LOG_ERROR("Missing number of atoms, is required for normal modes");
		return false;
	}

	hid_t obj_info = H5Oopen(handle, "normal_modes", H5P_DEFAULT);
	if (obj_info == H5I_INVALID_HID) {
		MD_LOG_ERROR("Failed to open 'normal_modes' object");
		return false;
	}

	const size_t mode_len = vlx->number_of_atoms * 3;
	md_attribute_format_t mode_format = {
		.type = MD_ATTRIBUTE_TYPE_F64, .components = 3, .rank = 2,
		.shape = { (uint32_t)number_of_modes, (uint32_t)vlx->number_of_atoms },
	};

	H5I_type_t obj_type = H5Iget_type(obj_info);
	if (obj_type == H5I_GROUP) {
		hid_t normal_modes_id = H5Gopen(handle, "normal_modes", H5P_DEFAULT);
		if (normal_modes_id != H5I_INVALID_HID) {
			md_attribute_id_t id = vlx_publish(sys, STR_LIT("qm/atom/normal_mode"), STR_LIT("Normal Mode"), md_unit_none(), mode_format, NULL, 0);
			double* dst = id != MD_ATTRIBUTE_INVALID ? (double*)md_attributes_data(&sys->attributes, id, MD_ATTRIBUTE_TYPE_F64) : NULL;
			if (!dst) {
				if (id != MD_ATTRIBUTE_INVALID) md_attributes_remove(&sys->attributes, id);
				H5Gclose(normal_modes_id);
				H5Oclose(obj_info);
				return false;
			}

			char lbl[32];
			for (size_t i = 0; i < number_of_modes; ++i) {
				snprintf(lbl, sizeof(lbl), "%zu", i + 1);
				if (!h5_read_dataset_data(dst + i * mode_len, mode_len, normal_modes_id, H5T_NATIVE_DOUBLE, lbl)) {
					MD_LOG_ERROR("Failed to extract dataset in '%s' normal mode", lbl);
					md_attributes_remove(&sys->attributes, id);
					H5Gclose(normal_modes_id);
					H5Oclose(obj_info);
					return false;
				}
			}
			H5Gclose(normal_modes_id);
		}
	} else if (obj_type == H5I_DATASET) {
		size_t data_dim[3];
		int num_dim = h5_read_dataset_dims(data_dim, 3, handle, "normal_modes");

		if (num_dim != 3 || data_dim[0] != number_of_modes || data_dim[1] != vlx->number_of_atoms || data_dim[2] != 3) {
			MD_LOG_ERROR("Unexpected dimensions in normal_modes dataset");
			H5Oclose(obj_info);
			return false;
		}

		if (!vlx_publish_h5(sys, handle, "normal_modes", STR_LIT("qm/atom/normal_mode"), STR_LIT("Normal Mode"), md_unit_none(), mode_format)) {
			MD_LOG_ERROR("Failed to read normal_modes dataset");
			H5Oclose(obj_info);
			return false;
		}
	} else {
		MD_LOG_ERROR("Unrecognized object type for 'normal_modes'");
		H5Oclose(obj_info);
		return false;
	}
	H5Oclose(obj_info);

	size_t number_of_external_frequencies = 0;
	if (h5_read_scalar(&number_of_external_frequencies, handle, H5T_NATIVE_HSIZE, "number_of_external_frequencies")) {
		if (h5_check_dataset_exists(handle, "external_frequencies") &&
			!vlx_publish_h5_series(sys, handle, "external_frequencies", STR_LIT("vlx/vib/external_frequency"), STR_LIT("External Frequency"), vlx_unit_hartree(), number_of_external_frequencies)) {
			return false;
		}

		// {E,M}: one row of per mode activities per external frequency, exactly how the dataset is
		// stored. The external frequency axis leads for the same reason the mode axis leads the
		// displacements - one row is contiguous, so a consumer plotting the spectrum at one
		// frequency hands a plotting library a pointer rather than a stride.
		if (h5_check_dataset_exists(handle, "raman_activities") &&
			!vlx_publish_h5_matrix(sys, handle, "raman_activities", STR_LIT("vlx/vib/raman_activity"), STR_LIT("Raman Activity"), md_unit_none(), number_of_external_frequencies, number_of_modes)) {
			return false;
		}
	}

	return true;
}

// This procedure is an abstraction to help read illformed groups containing a collection of individual datasets containing one value each.
// The expected names of the datasets are '0', '1', ... up to the expected count. The dataset prefix is used for error messages only, it does not have to be present in the actual dataset names.
static bool h5_extract_group_as_array_double(md_array(double)* out_data, hid_t group_handle, md_allocator_i* arena) {
	ASSERT(out_data);
	ASSERT(arena);

	hsize_t num_links = 0;
	if (H5Gget_num_objs(group_handle, &num_links) < 0) {
		MD_LOG_ERROR("Failed to get number of links in group");
		return false;
	}

	size_t expected_count = (size_t)num_links;
	md_array_ensure(*out_data, expected_count, arena);

	for (size_t i = 0; i < expected_count; ++i) {
		char dataset_name[64];
		snprintf(dataset_name, sizeof(dataset_name), "%zu", i);
		if (!h5_check_dataset_exists(group_handle, dataset_name)) {
			MD_LOG_ERROR("Expected dataset '%s' in group not found", dataset_name);
			return false;
		}

		double value;
		if (!h5_read_dataset_data(&value, 1, group_handle, H5T_NATIVE_DOUBLE, dataset_name)) {
			MD_LOG_ERROR("Failed to read dataset '%s' in group", dataset_name);
			return false;
		}
		md_array_push_no_grow(*out_data, value);
	}
	return true;
}

static bool h5_extract_as_array_double(md_array(double)* out_data, hid_t handle, const char* name, md_allocator_i* arena) {
	ASSERT(out_data);
	ASSERT(name);
	ASSERT(arena);

	bool result = false;
	
	H5I_type_t obj_type = h5_get_object_type(handle, name);
	if (obj_type == H5I_GROUP) {
		hid_t obj_handle = H5Oopen(handle, name, H5P_DEFAULT);
		if (obj_handle < 0) {
			MD_LOG_ERROR("Failed to open object '%s'", name);
			return false;
		}
		result = h5_extract_group_as_array_double(out_data, obj_handle, arena);
		H5Oclose(obj_handle);
	} else if (obj_type == H5I_DATASET) {
		size_t dim;
		int num_dim = h5_read_dataset_dims(&dim, 1, handle, name);
		if (num_dim <= 0) {
			MD_LOG_ERROR("Invalid dimensions in dataset '%s'", name);
			goto done;
		}
		md_array_resize(*out_data, dim, arena);
		MEMSET(*out_data, 0, md_array_bytes(*out_data));
		result = h5_read_dataset_data(*out_data, dim, handle, H5T_NATIVE_DOUBLE, name);
	} else {
		MD_LOG_ERROR("Unrecognized object type for '%s'", name);
	}

done:
	return result;
}

static bool h5_read_opt_data(vlx_t* vlx, hid_t handle) {
	const vlx_opt_type_t opt_types[] = { VLX_OPT_GEOMETRY, VLX_OPT_CONSTRAINED, VLX_OPT_IRC };
	const char* valid_prefixes[] = { "opt", "scan", "irc" };
	const char* energy_ident = NULL;
	const char* coord_ident = NULL;
	H5I_type_t  coord_type = H5I_BADID;
    vlx_opt_type_t opt_type = VLX_OPT_UNKNOWN;

	for (size_t i = 0; i < ARRAY_SIZE(valid_prefixes); ++i) {
		H5I_type_t type = -1;
		char energy_name[32];
		char coord_name[32];

		snprintf(energy_name, sizeof(energy_name), "%s_energies", valid_prefixes[i]);
		snprintf(coord_name, sizeof(coord_name), "%s_coordinates_au", valid_prefixes[i]);

		type = h5_get_object_type(handle, coord_name);
		if (type == H5I_DATASET || type == H5I_GROUP) {
			coord_ident = coord_name;
			coord_type = type;
		}

		type = h5_get_object_type(handle, energy_name);
		if (type == H5I_DATASET || type == H5I_GROUP) {
			energy_ident = energy_name;
			opt_type = opt_types[i];
			break;
		}
	}

	if (energy_ident) {
		md_system_t* sys = vlx->sys;

		// The energies come as either one dataset or a group of one value datasets, so they are
		// gathered locally first - the only thing here that cannot be read straight into the table.
		md_array(double) energies = 0;
		if (!h5_extract_as_array_double(&energies, handle, energy_ident, vlx->arena)) {
			return false;
		}

		const size_t len = md_array_size(energies);
		if (len == 0) {
			return false;
		}

		md_attribute_format_t energy_format = {
			.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)len },
		};
		vlx_publish(sys, STR_LIT("vlx/opt/energy"), STR_LIT("Energy"), vlx_unit_hartree(), energy_format, energies, len * sizeof(double));
		md_array_free(energies, vlx->arena);

		// WHICH optimisation. A geometry optimisation and an IRC scan are the same energy-per-step
		// column, and only the reader knows that the middle of one of them is a transition state
		// rather than a step on the way down.
		vlx_publish_str(sys, STR_LIT("vlx/opt/type"), STR_LIT("Optimization Type"), vlx_opt_type_str(opt_type));

		// WHICH electronic state was optimised - 0 is the ground state. Only a geometry optimisation
		// has one; a transition state search or an IRC scan walks a path rather than relaxing a
		// state, and publishing 0 for those would name a state they never chose.
		size_t state_index = 0;
		if (h5_check_dataset_exists(handle, "state_index") && !h5_read_scalar(&state_index, handle, H5T_NATIVE_INT64, "state_index")) {
			return false;
		}
		if (opt_type == VLX_OPT_GEOMETRY) {
			// 0 when the file does not name one, which is the ground state - the state a geometry
			// optimisation relaxes unless it says otherwise.
			vlx_publish_scalar(sys, STR_LIT("vlx/opt/state_index"), STR_LIT("Optimized State"), md_unit_none(), (double)state_index);
		}

		// Only for an IRC, where it names the step the path was walked out from in both directions -
		// the energy every other step is measured against. Any other run has no such step, and
		// publishing 0 there would name one.
		size_t ts_index = 0;
		if (h5_check_dataset_exists(handle, "ts_index") && !h5_read_scalar(&ts_index, handle, H5T_NATIVE_INT64, "ts_index")) {
			return false;
		}
		if (opt_type == VLX_OPT_IRC) {
			vlx_publish_scalar(sys, STR_LIT("vlx/opt/irc_ts_index"), STR_LIT("Transition State Step"), md_unit_none(), (double)ts_index);
		}

		if (coord_ident && coord_type == H5I_DATASET) {
			size_t dim[4];
			int num_dim = h5_read_dataset_dims(dim, ARRAY_SIZE(dim), handle, coord_ident);
			if (num_dim <= 0) {
				MD_LOG_ERROR("Invalid dimensions in '%s'", coord_ident);
				return false;
			}

			if (dim[1] != vlx->number_of_atoms || dim[2] != 3) {
				MD_LOG_ERROR("Unexpected dimensions in '%s'", coord_ident);
				return false;
			}

			if (dim[0] != len) {
				MD_LOG_ERROR("Energy/coordinate step count mismatch between '%s' and '%s'", energy_ident, coord_ident);
				return false;
			}

			// {P,N} of 3 component values, read straight into the table and converted to Angstrom in
			// place - the same unit the system's own coordinates are in, so a consumer comparing a
			// step against the loaded geometry does not have to convert first.
			md_attribute_format_t coord_format = {
				.type = MD_ATTRIBUTE_TYPE_F64, .components = 3, .rank = 2,
				.shape = { (uint32_t)dim[0], (uint32_t)dim[1] },
			};
			if (!vlx_publish_h5(sys, handle, coord_ident, STR_LIT("vlx/opt/coordinate"), STR_LIT("Coordinate"), vlx_unit_angstrom(), coord_format)) {
				return false;
			}

			const md_attribute_t* attr = md_attributes_find(&sys->attributes, STR_LIT("vlx/opt/coordinate"));
			double* coord = attr ? (double*)md_attributes_data(&sys->attributes, attr->id, MD_ATTRIBUTE_TYPE_F64) : NULL;
			if (coord) {
				for (size_t i = 0; i < dim[0] * dim[1] * 3; ++i) {
					coord[i] *= BOHR_TO_ANGSTROM;
				}
			}
		}

		return true;
	}

	return false;
}

static bool h5_read_core_data(vlx_t* vlx, hid_t handle) {
	ASSERT(vlx);

	if (!h5_read_str(&vlx->basis_set_ident, handle, "basis_set", vlx->arena)) {
		return false;
	}

	if (!h5_read_str(&vlx->dft_func_label, handle, "dft_func_label", vlx->arena)) {
		return false;
	}

	if (!h5_read_scalar(&vlx->molecular_charge, handle, H5T_NATIVE_DOUBLE, "molecular_charge")) {
		return false;
	}

	if (!h5_read_scalar(&vlx->nuclear_repulsion_energy, handle, H5T_NATIVE_DOUBLE, "nuclear_repulsion")) {
		return false;
	}

	if (!h5_read_scalar(&vlx->number_of_alpha_electrons, handle, H5T_NATIVE_INT64, "number_of_alpha_electrons")) {
		return false;
	}

	if (!h5_read_scalar(&vlx->number_of_atoms, handle, H5T_NATIVE_INT64, "number_of_atoms")) {
		return false;
	}

	if (!h5_read_scalar(&vlx->number_of_beta_electrons, handle, H5T_NATIVE_INT64, "number_of_beta_electrons")) {
		return false;
	}

	if (!h5_read_str(&vlx->potfile_text, handle, "potfile_text", vlx->arena)) {
		return false;
	}

	if (!h5_read_scalar(&vlx->spin_multiplicity, handle, H5T_NATIVE_INT64, "spin_multiplicity")) {
		return false;
	}

	if (vlx->number_of_atoms == 0) {
		MD_LOG_ERROR("Number of atoms is zero");
		return false;
	}

	md_array_resize(vlx->atom_coordinates, vlx->number_of_atoms, vlx->arena);
	MEMSET(vlx->atom_coordinates, 0, md_array_bytes(vlx->atom_coordinates));
	if (!h5_read_dataset_data(vlx->atom_coordinates, md_array_size(vlx->atom_coordinates) * 3, handle, H5T_NATIVE_DOUBLE, "atom_coordinates")) {
		return false;
	}

	md_array_resize(vlx->atomic_numbers, vlx->number_of_atoms, vlx->arena);
	MEMSET(vlx->atomic_numbers, 0, md_array_bytes(vlx->atomic_numbers));
	if (!h5_read_dataset_data(vlx->atomic_numbers, md_array_size(vlx->atomic_numbers), handle, H5T_NATIVE_UINT8, "nuclear_charges")) {
		return false;
	}

	// Convert Atomic units to Ångström
	if (vlx->atom_coordinates) {
		for (size_t i = 0; i < vlx->number_of_atoms; ++i) {
			vlx->atom_coordinates[i].x *= BOHR_TO_ANGSTROM;
			vlx->atom_coordinates[i].y *= BOHR_TO_ANGSTROM;
			vlx->atom_coordinates[i].z *= BOHR_TO_ANGSTROM;
		}
	}

	if (h5_check_dataset_exists(handle, "qm_atom_indices")) {
        md_array_resize(vlx->local_to_global_atom_idx, vlx->number_of_atoms, vlx->arena);
		if (!h5_read_dataset_data(vlx->local_to_global_atom_idx, vlx->number_of_atoms, handle, H5T_NATIVE_INT32, "qm_atom_indices")) {
			return false;
		}
	}

	return true;
}

static const double* vlx_rsp_get_solution_vector(const vlx_t* vlx, size_t state_idx, size_t* out_nocc, size_t* out_nvir, size_t* out_amp_count, bool* out_has_y) {
	ASSERT(vlx);

	if (!vlx->rsp.solution_matrix.data || state_idx >= vlx->rsp.solution_matrix.size[0]) {
		return NULL;
	}

	const size_t nocc = vlx->rsp.num_core > 0 ? vlx->rsp.num_core : vlx->rsp.num_valence;
	if (nocc == 0) {
		return NULL;
	}

	const size_t nvir = vlx->rsp.num_virtual;
	if (nvir == 0) {
		return NULL;
	}

	size_t vec_size = vlx->rsp.solution_matrix.size[1];

	const size_t amp_count = nocc * nvir;
	bool has_y = false;
	if (vec_size == amp_count) {
		has_y = false;
	} else if (vec_size == 2 * amp_count) {
		has_y = true;
	} else {
		MD_LOG_ERROR("Unexpected response eigenvector length for state %zu: got %zu, expected %zu or %zu", state_idx + 1, vec_size, amp_count, 2 * amp_count);
		return NULL;
	}

	if (out_nocc) *out_nocc = nocc;
	if (out_nvir) *out_nvir = nvir;
	if (out_amp_count) *out_amp_count = amp_count;
	if (out_has_y) *out_has_y = has_y;

    return vlx->rsp.solution_matrix.data + state_idx * vec_size;
}

static void vlx_symmetrize_square(double* mat, size_t dim) {
	ASSERT(mat);
	for (size_t i = 0; i < dim; ++i) {
		for (size_t j = i + 1; j < dim; ++j) {
			const double value = 0.5 * (mat[i * dim + j] + mat[j * dim + i]);
			mat[i * dim + j] = value;
			mat[j * dim + i] = value;
		}
	}
}

static void vlx_transform_mo_density_to_ao(double* out_ao, const double* mo_density, const double* coeff, size_t mo_offset, size_t mo_count, size_t num_ao) {
	ASSERT(out_ao);
	ASSERT(mo_density);
	ASSERT(coeff);

	md_temp_scope_t temp = md_temp_begin();
	double* work = md_temp_alloc_array(temp, double, mo_count * num_ao);

	for (size_t i = 0; i < mo_count; ++i) {
		for (size_t ao = 0; ao < num_ao; ++ao) {
			double sum = 0.0;
			for (size_t j = 0; j < mo_count; ++j) {
				sum += mo_density[i * mo_count + j] * coeff[(mo_offset + j) * num_ao + ao];
			}
			work[i * num_ao + ao] = sum;
		}
	}

	for (size_t ao_i = 0; ao_i < num_ao; ++ao_i) {
		for (size_t ao_j = 0; ao_j < num_ao; ++ao_j) {
			double sum = 0.0;
			for (size_t i = 0; i < mo_count; ++i) {
				sum += coeff[(mo_offset + i) * num_ao + ao_i] * work[i * num_ao + ao_j];
			}
			out_ao[ao_i * num_ao + ao_j] = sum;
		}
	}

	md_temp_end(temp);
}

// The whole of the transition density reconstruction, in terms of a plain solution vector and a
// plain AO coefficient matrix rather than a vlx_t. This is what lets the attribute provider
// below reconstruct the same matrix with no vlx object in reach: everything it needs is one row of
// a response solution and the (already resident) alpha MO coefficients.
static bool vlx_build_transition_density_matrix(double* out_matrix, const double* solution_vector, size_t vec_len,
	size_t nocc, size_t nvir, const double* coeff, size_t num_ao, vlx_transition_type_t type)
{
	ASSERT(out_matrix);
	ASSERT(solution_vector);
	ASSERT(coeff);

	if (nocc == 0 || nvir == 0 || num_ao == 0) {
		return false;
	}

	const size_t amp_count = nocc * nvir;
	bool has_y;
	if (vec_len == amp_count) {
		has_y = false;
	} else if (vec_len == 2 * amp_count) {
		has_y = true;
	} else {
		MD_LOG_ERROR("Response solution vector holds %zu values, expected %zu or %zu (%zu occupied x %zu virtual)", vec_len, amp_count, 2 * amp_count, nocc, nvir);
		return false;
	}

	md_temp_scope_t temp = md_temp_begin();
	double* detach_mo = md_temp_alloc_array(temp, double, nocc * nocc);
	double* attach_mo = md_temp_alloc_array(temp, double, nvir * nvir);
	double* detach_ao = NULL;
	MEMSET(detach_mo, 0, sizeof(double) * nocc * nocc);
	MEMSET(attach_mo, 0, sizeof(double) * nvir * nvir);

	for (size_t i = 0; i < nocc; ++i) {
		for (size_t j = i; j < nocc; ++j) {
			double value = 0.0;
			for (size_t a = 0; a < nvir; ++a) {
				const size_t ia = i * nvir + a;
				const size_t ja = j * nvir + a;
                const double z_i = solution_vector[ia];
                const double z_j = solution_vector[ja];
                const double y_i = has_y ? solution_vector[amp_count + ia] : 0.0;
                const double y_j = has_y ? solution_vector[amp_count + ja] : 0.0;
				const double t_i = z_i - y_i;
				const double t_j = z_j - y_j;
				value += t_i * t_j;
			}
			detach_mo[i * nocc + j] = value;
			detach_mo[j * nocc + i] = value;
		}
	}

	for (size_t a = 0; a < nvir; ++a) {
		for (size_t b = a; b < nvir; ++b) {
			double value = 0.0;
			for (size_t i = 0; i < nocc; ++i) {
				const size_t ia = i * nvir + a;
				const size_t ib = i * nvir + b;
                const double z_i = solution_vector[ia];
                const double z_j = solution_vector[ib];
                const double y_i = has_y ? solution_vector[amp_count + ia] : 0.0;
                const double y_j = has_y ? solution_vector[amp_count + ib] : 0.0;
				const double t_a = z_i - y_i;
				const double t_b = z_j - y_j;
				value += t_a * t_b;
			}
			attach_mo[a * nvir + b] = value;
			attach_mo[b * nvir + a] = value;
		}
	}

	if (type == VLX_TRANSITION_DETACHMENT) {
		vlx_transform_mo_density_to_ao(out_matrix, detach_mo, coeff, 0, nocc, num_ao);
	} else {
		vlx_transform_mo_density_to_ao(out_matrix, attach_mo, coeff, nocc, nvir, num_ao);
		if (type == VLX_TRANSITION_DIFFERENCE) {
			detach_ao = md_temp_alloc_array(temp, double, num_ao * num_ao);
			vlx_transform_mo_density_to_ao(detach_ao, detach_mo, coeff, 0, nocc, num_ao);
			for (size_t i = 0; i < num_ao * num_ao; ++i) {
				out_matrix[i] -= detach_ao[i];
			}
		}
	}

	vlx_symmetrize_square(out_matrix, num_ao);
	md_temp_end(temp);
	return true;
}

static double vlx_dot(const double* a, const double* b, size_t count) {
	double sum = 0.0;
	for (size_t i = 0; i < count; ++i) {
		sum += a[i] * b[i];
	}
	return sum;
}

static double vlx_normalize(double* vec, size_t count) {
	double norm = sqrt(vlx_dot(vec, vec, count));
	if (norm > 0.0) {
		const double inv_norm = 1.0 / norm;
		for (size_t i = 0; i < count; ++i) {
			vec[i] *= inv_norm;
		}
	}
	return norm;
}

static void vlx_orthogonalize(double* vec, const double* basis, size_t basis_count, size_t dim) {
	for (size_t basis_idx = 0; basis_idx < basis_count; ++basis_idx) {
		const double* b = basis + basis_idx * dim;
		const double projection = vlx_dot(vec, b, dim);
		for (size_t i = 0; i < dim; ++i) {
			vec[i] -= projection * b[i];
		}
	}
}

static bool vlx_symmetric_top_eigenpairs(double* out_values, double* out_vectors, const double* matrix, size_t dim, size_t pair_count) {
	ASSERT(out_values);
	ASSERT(out_vectors);
	ASSERT(matrix);

	if (dim == 0 || pair_count == 0) {
		return false;
	}

	md_temp_scope_t temp = md_temp_begin();
	double* work = md_temp_alloc_array(temp, double, dim * dim);
	double* vec  = md_temp_alloc_array(temp, double, dim);
	double* next = md_temp_alloc_array(temp, double, dim);
	MEMCPY(work, matrix, sizeof(double) * dim * dim);
	MEMSET(out_values, 0, sizeof(double) * pair_count);
	MEMSET(out_vectors, 0, sizeof(double) * pair_count * dim);

	for (size_t pair_idx = 0; pair_idx < pair_count; ++pair_idx) {
		for (size_t i = 0; i < dim; ++i) {
			vec[i] = 1.0 + 0.013 * (double)(((i + 1) * (pair_idx + 3)) % 17);
		}
		vlx_orthogonalize(vec, out_vectors, pair_idx, dim);
		if (vlx_normalize(vec, dim) <= DBL_EPSILON) {
			vec[pair_idx % dim] = 1.0;
			vlx_orthogonalize(vec, out_vectors, pair_idx, dim);
			vlx_normalize(vec, dim);
		}

		for (size_t iter = 0; iter < VLX_NTO_POWER_ITERATIONS; ++iter) {
			for (size_t i = 0; i < dim; ++i) {
				double sum = 0.0;
				for (size_t j = 0; j < dim; ++j) {
					sum += work[i * dim + j] * vec[j];
				}
				next[i] = sum;
			}

			vlx_orthogonalize(next, out_vectors, pair_idx, dim);
			if (vlx_normalize(next, dim) <= DBL_EPSILON) {
				break;
			}

			double diff = 0.0;
			double neg_diff = 0.0;
			for (size_t i = 0; i < dim; ++i) {
				const double d = next[i] - vec[i];
				const double nd = next[i] + vec[i];
				diff += d * d;
				neg_diff += nd * nd;
				vec[i] = next[i];
			}
			if (sqrt(MIN(diff, neg_diff)) < VLX_NTO_CONVERGENCE_EPSILON) {
				break;
			}
		}

		double eigenvalue = 0.0;
		for (size_t i = 0; i < dim; ++i) {
			double row_sum = 0.0;
			for (size_t j = 0; j < dim; ++j) {
				row_sum += work[i * dim + j] * vec[j];
			}
			eigenvalue += vec[i] * row_sum;
		}

		if (eigenvalue < VLX_NTO_EIGENVALUE_EPSILON) {
			break;
		}

		out_values[pair_idx] = eigenvalue;
		MEMCPY(out_vectors + pair_idx * dim, vec, sizeof(double) * dim);

		for (size_t i = 0; i < dim; ++i) {
			for (size_t j = 0; j < dim; ++j) {
				work[i * dim + j] -= eigenvalue * vec[i] * vec[j];
			}
		}
	}

	md_temp_end(temp);
	return true;
}

static size_t vlx_rsp_extract_nto_from_solution(double* out_coefficients, double* out_lambdas, const vlx_t* vlx, size_t state_idx, vlx_nto_type_t type, size_t lambda_count) {
	ASSERT(vlx);

	if (lambda_count == 0 || (!out_coefficients && !out_lambdas)) {
		return 0;
	}

	size_t nocc = 0;
	size_t nvir = 0;
	size_t amp_count = 0;
	bool has_y = false;
    const double* solution_vector = vlx_rsp_get_solution_vector(vlx, state_idx, &nocc, &nvir, &amp_count, &has_y);
    if (!solution_vector) {
        return 0;
    }

	const size_t pair_count = MIN(MIN(nocc, nvir), lambda_count);
	if (pair_count == 0) {
		return 0;
	}

	const vlx_2d_data_t* scf_coeff = &vlx->scf.alpha.coefficients;
	const size_t num_ao = scf_coeff->size[1];
	const double* coeff = scf_coeff->data;

	if (out_coefficients) {
		MEMSET(out_coefficients, 0, sizeof(double) * lambda_count * num_ao);
	}
	if (out_lambdas) {
		MEMSET(out_lambdas, 0, sizeof(double) * lambda_count);
	}

	const bool use_left_gram = nocc <= nvir;
	const size_t small_dim = use_left_gram ? nocc : nvir;
	const size_t large_dim = use_left_gram ? nvir : nocc;

	size_t written_count = 0;
	md_temp_scope_t temp = md_temp_begin();
	double* transition = md_temp_alloc_array(temp, double, nocc * nvir);
	double* gram = md_temp_alloc_array(temp, double, small_dim * small_dim);
	double* eigenvalues = md_temp_alloc_array(temp, double, pair_count);
	double* small_vectors = md_temp_alloc_array(temp, double, pair_count * small_dim);
	double* large_vectors = md_temp_alloc_array(temp, double, pair_count * large_dim);
	MEMSET(gram, 0, sizeof(double) * small_dim * small_dim);
	MEMSET(large_vectors, 0, sizeof(double) * pair_count * large_dim);

	for (size_t i = 0; i < nocc; ++i) {
		for (size_t a = 0; a < nvir; ++a) {
			const size_t idx = i * nvir + a;
            const double z = solution_vector[idx];
            const double y = has_y ? solution_vector[amp_count + idx] : 0.0;
			transition[idx] = z - y;
		}
	}

	if (use_left_gram) {
		for (size_t i = 0; i < nocc; ++i) {
			for (size_t j = i; j < nocc; ++j) {
				double value = 0.0;
				for (size_t a = 0; a < nvir; ++a) {
					value += transition[i * nvir + a] * transition[j * nvir + a];
				}
				gram[i * nocc + j] = value;
				gram[j * nocc + i] = value;
			}
		}
	} else {
		for (size_t a = 0; a < nvir; ++a) {
			for (size_t b = a; b < nvir; ++b) {
				double value = 0.0;
				for (size_t i = 0; i < nocc; ++i) {
					value += transition[i * nvir + a] * transition[i * nvir + b];
				}
				gram[a * nvir + b] = value;
				gram[b * nvir + a] = value;
			}
		}
	}

	if (!vlx_symmetric_top_eigenpairs(eigenvalues, small_vectors, gram, small_dim, pair_count)) {
		goto done;
	}

	for (size_t pair_idx = 0; pair_idx < pair_count; ++pair_idx) {
		const double lambda = eigenvalues[pair_idx];
		if (lambda < VLX_NTO_EIGENVALUE_EPSILON) {
			break;
		}

		const double sigma = sqrt(lambda);
		if (sigma <= DBL_EPSILON) {
			break;
		}

		if (out_lambdas) {
			out_lambdas[pair_idx] = lambda;
		}

		if (out_coefficients) {
			const double* small = small_vectors + pair_idx * small_dim;
			double* large = large_vectors + pair_idx * large_dim;

			if (use_left_gram) {
				for (size_t a = 0; a < nvir; ++a) {
					double value = 0.0;
					for (size_t i = 0; i < nocc; ++i) {
						value += transition[i * nvir + a] * small[i];
					}
					large[a] = value / sigma;
				}
			} else {
				for (size_t i = 0; i < nocc; ++i) {
					double value = 0.0;
					for (size_t a = 0; a < nvir; ++a) {
						value += transition[i * nvir + a] * small[a];
					}
					large[i] = value / sigma;
				}
			}
			vlx_normalize(large, large_dim);

			const double* u = use_left_gram ? small : large;
			const double* v = use_left_gram ? large : small;
			double* out_coeff = out_coefficients + pair_idx * num_ao;

			if (type == VLX_NTO_PARTICLE) {
				for (size_t ao = 0; ao < num_ao; ++ao) {
					double value = 0.0;
					for (size_t a = 0; a < nvir; ++a) {
						value += coeff[(nocc + a) * num_ao + ao] * v[a];
					}
					out_coeff[ao] = value;
				}
			} else if (type == VLX_NTO_HOLE) {
				for (size_t ao = 0; ao < num_ao; ++ao) {
					double value = 0.0;
					for (size_t i = 0; i < nocc; ++i) {
						value += coeff[i * num_ao + ao] * u[i];
					}
					out_coeff[ao] = value;
				}
			} else {
				break;
			}
		}

		written_count++;
	}

done:
	md_temp_end(temp);
	return written_count;
}

static size_t vlx_rsp_extract_nto(double* out_coefficients, double* out_lambdas, const vlx_t* vlx, size_t state_idx, vlx_nto_type_t type, size_t lambda_count) {
	if (!vlx || state_idx >= vlx->rsp.number_of_frequencies || lambda_count == 0) {
		return 0;
	}

	if (vlx->rsp.solution_matrix.data && state_idx < vlx->rsp.solution_matrix.size[0]) {
		return vlx_rsp_extract_nto_from_solution(out_coefficients, out_lambdas, vlx, state_idx, type, lambda_count);
	}

	return 0;
}
static bool vlx_read_scf_results(vlx_t* vlx, str_t filename, md_system_state_t* state) {
	ASSERT(vlx);

	// Ensure a zero terminated string for interfacing to HDF5
	char buf[2048];
	str_copy_to_char_buf(buf, sizeof(buf), filename);

	h5_error_scope_t h5_err_scope = h5_error_scope_begin();

	// Open an existing file
	hid_t file_id = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Could not open HDF5 file: '"STR_FMT"'", STR_ARG(filename));
		h5_error_scope_end(h5_err_scope);
		return false;
	}

	bool result = false;

	if (!h5_read_core_data(vlx, file_id) || !vlx_system_begin(vlx, state)) {
		goto done;
	}

	if (!h5_read_scf_data(vlx, file_id)) {
		goto done;
	}

	result = true;
done:
	H5Fclose(file_id);
	h5_error_scope_end(h5_err_scope);

	return result;
}

// This is the newest version of the file format where everything is contained within a single h5 file
static bool h5_read_xps_data(vlx_t* vlx, hid_t handle);

static bool vlx_read_h5_file(vlx_t* vlx, str_t filename, md_system_state_t* state) {
	ASSERT(vlx);

	// Ensure a zero terminated string for interfacing to HDF5
	char buf[2048];
	str_copy_to_char_buf(buf, sizeof(buf), filename);

	h5_error_scope_t h5_err_scope = h5_error_scope_begin();

	// Open an existing file
	hid_t file_id = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Could not open HDF5 file: '"STR_FMT"'", STR_ARG(filename));
		h5_error_scope_end(h5_err_scope);
		return false;
	}

	bool result = false;

	// The core block first, and the system built from it before anything else is read: from here on
	// every reader publishes into vlx->sys as it goes, so the table has to exist and must not be
	// reset again afterwards.
	if (!h5_read_core_data(vlx, file_id) || !vlx_system_begin(vlx, state)) {
		goto done;
	}

	// SCF
	{
		if (H5Lexists(file_id, "scf", H5P_DEFAULT) > 0) {
			hid_t scf_id = H5Gopen(file_id, "scf", H5P_DEFAULT);
			if (scf_id != H5I_INVALID_HID) {
				result = h5_read_scf_data(vlx, scf_id);
				H5Gclose(scf_id);
				if (!result) goto done;
			}
		}
	}

	// VIB
    {
        if (H5Lexists(file_id, "vib", H5P_DEFAULT) > 0) {
            hid_t vib_id = H5Gopen(file_id, "vib", H5P_DEFAULT);
            if (vib_id != H5I_INVALID_HID) {
                result = h5_read_vib_data(vlx, vib_id);
                H5Gclose(vib_id);
                if (!result) goto done;
            }
        }
    }

	// OPT
	{
		if (H5Lexists(file_id, "opt", H5P_DEFAULT) > 0) {
			hid_t opt_id = H5Gopen(file_id, "opt", H5P_DEFAULT);
			if (opt_id != H5I_INVALID_HID) {
				result = h5_read_opt_data(vlx, opt_id);
				H5Gclose(opt_id);
				if (!result) goto done;
			}
		}
	}

	// RSP
	{
		if (H5Lexists(file_id, "rsp", H5P_DEFAULT) > 0) {
			hid_t rsp_id = H5Gopen(file_id, "rsp", H5P_DEFAULT);
			if (rsp_id != H5I_INVALID_HID) {
				result = h5_read_rsp_data(vlx, rsp_id);
				H5Gclose(rsp_id);
				if (!result) goto done;
			}
		}
	}

	// XPS. Optional top level group, independent of the response block: XPS is delta-SCF, so it may
	// appear alongside any vlx_rsp_type_t or with no response data at all.
	{
		if (H5Lexists(file_id, "xps", H5P_DEFAULT) > 0) {
			hid_t xps_id = H5Gopen(file_id, "xps", H5P_DEFAULT);
			if (xps_id != H5I_INVALID_HID) {
				result = h5_read_xps_data(vlx, xps_id);
				H5Gclose(xps_id);
				if (!result) goto done;
			}
		}
	}

	if (!h5_read_atomic_properties(vlx, file_id)) {
		goto done;
	}

	if (!h5_read_density_properties(vlx, file_id)) {
		goto done;
	}

	result = true;
done:
	H5Fclose(file_id);
	h5_error_scope_end(h5_err_scope);

	return result;
}

#define BAKE_STR(str) {str "", sizeof(str) - 1}

static inline str_t resolve_basis_set_ident(str_t input) {
	struct map_t {
        str_t in;
        str_t out;
	};

	static const struct map_t alias_table[] = {
        {BAKE_STR("6-31G*"),			BAKE_STR("6-31G_D_")},
		{BAKE_STR("6-31G**"),			BAKE_STR("6-31G_D,P_")},
		{BAKE_STR("6-31+G*"),			BAKE_STR("6-31+G_D_")},
        {BAKE_STR("6-31+G**"),			BAKE_STR("6-31+G_D,P_")},
		{BAKE_STR("6-31++G*"),			BAKE_STR("6-31++G_D_")},
		{BAKE_STR("6-31++G**"),			BAKE_STR("6-31++G_D,P_")},
        {BAKE_STR("6-311G*"),			BAKE_STR("6-311G_D_")},
		{BAKE_STR("6-311G**"),			BAKE_STR("6-311G_D,P_")},
		{BAKE_STR("6-311+G*"),			BAKE_STR("6-311+G_D_")},
        {BAKE_STR("6-311+G**"),			BAKE_STR("6-311+G_D,P_")},
		{BAKE_STR("6-311++G*"),			BAKE_STR("6-311++G_D_")},
		{BAKE_STR("6-311++G**"),		BAKE_STR("6-311++G_D,P_")},
        {BAKE_STR("6-31G(2DF,P)"),		BAKE_STR("6-31G_2DF,P_")},
		{BAKE_STR("6-31G(3DF,3PD)"),	BAKE_STR("6-31G_3DF,3PD_")},
		{BAKE_STR("6-311G(2DF,2PD)"),	BAKE_STR("6-311G_2DF,2PD_")},
        {BAKE_STR("6-311+G(2D,P)"),		BAKE_STR("6-311+G_2D,P_")},
		{BAKE_STR("6-311++G(2D,2P)"),	BAKE_STR("6-311++G_2D,2P_")},
		{BAKE_STR("6-311++G(3DF,3PD)"),	BAKE_STR("6-311++G_3DF,3PD_")},
        {BAKE_STR("DEF2-SV(P)"),		BAKE_STR("DEF2-SV_P_")},
    };

	for (size_t i = 0; i < ARRAY_SIZE(alias_table); ++i) {
        if (str_eq(input, alias_table[i].in)) {
            return alias_table[i].out;
		}
	}

	return input;
}

#undef BAKE_STR

static int vlx_xps_compare_entry(const void* a, const void* b) {
	const vlx_xps_entry_t* ea = (const vlx_xps_entry_t*)a;
	const vlx_xps_entry_t* eb = (const vlx_xps_entry_t*)b;
	if (ea->element != eb->element) {
		return (ea->element < eb->element) ? -1 : 1;
	}
	if (ea->ionization_energy != eb->ionization_energy) {
		return (ea->ionization_energy < eb->ionization_energy) ? -1 : 1;
	}
	return 0;
}

// XPS
//
// Layout of the optional '/xps' group:
//
//   /xps/<element symbol>/<n>/atom_index            scalar, integer
//                            /contribution          scalar, double
//                            /ionization_energy_ev  scalar, double, unit eV
//                            /is_delocalized        scalar, bool (h5py enum over int8)
//                            /mo_index              scalar, integer
//
// <n> is a flat '0'..'N-1' enumeration with no meaning beyond ordering, so it is discarded --
// entries are re-sorted by (element, ionization_energy) in vlx_xps_finalize().
//
// This pushes into vlx->xps.entries but does NOT finalize. vlx_parse_file() does that once, after
// every push, because finalizing takes pointers into an array that md_array_push may still realloc.
static bool h5_read_xps_data(vlx_t* vlx, hid_t handle) {
	md_array(vlx_xps_entry_t) entries = 0;
	ASSERT(vlx);

	H5G_info_t info = { 0 };
	if (H5Gget_info(handle, &info) < 0) {
		MD_LOG_ERROR("XPS: failed to get group info");
		return false;
	}

	char name_buf[256];
	for (hsize_t i = 0; i < info.nlinks; ++i) {
		if (H5Gget_objname_by_idx(handle, i, name_buf, sizeof(name_buf)) < 0) {
			continue;
		}
		if (H5Gget_objtype_by_idx(handle, i) != H5G_GROUP) {
			continue;
		}

		const md_element_t element = md_util_element_lookup(str_from_cstr(name_buf), true);
		if (element == 0) {
			MD_LOG_INFO("XPS: skipping group '%s', not a recognized element symbol", name_buf);
			continue;
		}

		hid_t elem_group = H5Gopen(handle, name_buf, H5P_DEFAULT);
		if (elem_group == H5I_INVALID_HID) {
			MD_LOG_ERROR("XPS: failed to open element group '%s'", name_buf);
			continue;
		}

		hsize_t num_links = 0;
		if (H5Gget_num_objs(elem_group, &num_links) < 0) {
			MD_LOG_ERROR("XPS: failed to count entries in element group '%s'", name_buf);
			H5Gclose(elem_group);
			continue;
		}

		// Indexed by name rather than by link index, so entries are visited in numeric order
		// regardless of how HDF5 chose to order the links.
		for (hsize_t j = 0; j < num_links; ++j) {
			char idx_buf[32];
			snprintf(idx_buf, sizeof(idx_buf), "%i", (int)j);
			if (H5Lexists(elem_group, idx_buf, H5P_DEFAULT) <= 0) {
				continue;
			}

			hid_t entry_group = H5Gopen(elem_group, idx_buf, H5P_DEFAULT);
			if (entry_group == H5I_INVALID_HID) {
				continue;
			}

			vlx_xps_entry_t entry = {
				.atom_index = -1,
				.mo_index   = -1,
				.element    = element,
			};

			// The only field without a sensible default: an entry with no energy has nothing to plot.
			if (!h5_read_dataset_data(&entry.ionization_energy, 1, entry_group, H5T_NATIVE_DOUBLE, "ionization_energy_ev")) {
				MD_LOG_ERROR("XPS: '%s/%s' has no 'ionization_energy_ev', skipping entry", name_buf, idx_buf);
				H5Gclose(entry_group);
				continue;
			}

			// Remaining fields are optional and keep their defaults if absent.
			h5_read_dataset_data(&entry.contribution, 1, entry_group, H5T_NATIVE_DOUBLE, "contribution");
			h5_read_dataset_data(&entry.atom_index,   1, entry_group, H5T_NATIVE_INT32,  "atom_index");
			h5_read_dataset_data(&entry.mo_index,     1, entry_group, H5T_NATIVE_INT32,  "mo_index");

			// h5py writes Python bools as an HDF5 enum with an int8 base; H5Dread converts
			// enum -> integer for us, so reading it as int8 works for both that and a plain int.
			int8_t is_delocalized = 0;
			if (h5_read_dataset_data(&is_delocalized, 1, entry_group, H5T_NATIVE_INT8, "is_delocalized")) {
				entry.is_delocalized = (is_delocalized != 0);
			}

			if (entry.atom_index >= 0 && (size_t)entry.atom_index >= vlx->number_of_atoms) {
				MD_LOG_INFO("XPS: '%s/%s' has out of range atom_index %i (%zu atoms), clearing it",
					name_buf, idx_buf, (int)entry.atom_index, vlx->number_of_atoms);
				entry.atom_index = -1;
			}

			md_array_push(entries, entry, vlx->arena);
			H5Gclose(entry_group);
		}

		H5Gclose(elem_group);
	}

	// Gathered before anything is published, because the entries are SORTED by (element, ionization
	// energy) and a sort needs them all in hand. That ordering is the whole reason the per element
	// grouping does not have to be published: entries end up as contiguous runs of equal element, so
	// a consumer wanting one element's states scans vlx/xps/element for its run. Publishing the runs
	// as well would be two representations of one fact, with nothing keeping them in agreement.
	const size_t num_entries = md_array_size(entries);
	if (num_entries > 0) {
		qsort(entries, num_entries, sizeof(vlx_xps_entry_t), vlx_xps_compare_entry);

		// The file's record has six fields of four different types, which is six sibling paths over
		// one {C} index space rather than one attribute - a value has one type. Note what this buys:
		// a consumer plotting ionization energy against contribution now hands a plotting library
		// two CONTIGUOUS arrays instead of one base pointer and a struct stride.
		//
		// The bool field is copied a byte at a time as U8; a wider bool would take the wrong byte on
		// a big endian target, so it is worth failing the build rather than the render.
		STATIC_ASSERT(sizeof(bool) == 1, "XPS is_delocalized is published as a single byte");

		md_system_t* sys = vlx->sys;
		const size_t stride = sizeof(vlx_xps_entry_t);
		vlx_publish_column(sys, STR_LIT("vlx/xps/ionization_energy"), STR_LIT("Ionization Energy"),   md_unit_electronvolt(), MD_ATTRIBUTE_TYPE_F64, &entries->ionization_energy, stride, num_entries);
		vlx_publish_column(sys, STR_LIT("vlx/xps/contribution"),      STR_LIT("Core MO Contribution"), md_unit_none(),        MD_ATTRIBUTE_TYPE_F64, &entries->contribution,      stride, num_entries);
		vlx_publish_column(sys, STR_LIT("vlx/xps/atom_index"),        STR_LIT("Atom Index"),           md_unit_none(),        MD_ATTRIBUTE_TYPE_I32, &entries->atom_index,        stride, num_entries);
		vlx_publish_column(sys, STR_LIT("vlx/xps/mo_index"),          STR_LIT("MO Index"),             md_unit_none(),        MD_ATTRIBUTE_TYPE_I32, &entries->mo_index,          stride, num_entries);
		vlx_publish_column(sys, STR_LIT("vlx/xps/element"),           STR_LIT("Atomic Number"),        md_unit_none(),        MD_ATTRIBUTE_TYPE_U8,  &entries->element,           stride, num_entries);
		vlx_publish_column(sys, STR_LIT("vlx/xps/is_delocalized"),    STR_LIT("Is Delocalized"),       md_unit_none(),        MD_ATTRIBUTE_TYPE_U8,  &entries->is_delocalized,    stride, num_entries);
	}

	md_array_free(entries, vlx->arena);
	return true;
}

// Reads a file into vlx->sys. The blocks publish as they are read, so what is left to do here is
// everything that needs the whole file: resolving the basis set, moving the AO matrices into shell
// order and out of the spherical basis, and publishing what falls out of that.
static bool vlx_parse_file(vlx_t* vlx, str_t filename, md_system_state_t* state) {
	md_temp_scope_t temp = md_temp_begin();
	md_allocator_i* temp_alloc = md_temp_allocator(temp);

	bool result = false;

	if (str_ends_with(filename, STR_LIT(".scf.results.h5"))) {
		if (!vlx_read_scf_results(vlx, filename, state)) {
			goto done;
		}
	} else if (str_ends_with(filename, STR_LIT(".h5"))) {
		if (!vlx_read_h5_file(vlx, filename, state)) {
			goto done;
		}
	} else {
		MD_LOG_DEBUG("Unsupported file format");
		goto done;
	}

	if (!str_empty(vlx->basis_set_ident)) {
		size_t cap = KILOBYTES(16);
		char*  buf = md_temp_alloc_array(temp, char, cap);
		md_strb_t sb = md_strb_create(temp_alloc);

		str_t ident = resolve_basis_set_ident(vlx->basis_set_ident);
		MD_LOG_DEBUG("Basis set ident: '" STR_FMT "'", STR_ARG(ident));

		char exe_buf[1024];
		str_t exe_path = {exe_buf, md_path_write_exe(exe_buf, sizeof(exe_buf))};

		str_t exe_dir = {0};
		if (!extract_folder_path(&exe_dir, exe_path)) {
			MD_LOG_ERROR("Failed to extract executable directory");
		}

		md_strb_fmt(&sb, STR_FMT "%s/" STR_FMT, STR_ARG(exe_dir), MD_VLX_BASIS_FOLDER, STR_ARG(ident));
		str_t basis_filepath = md_strb_to_str(sb);
		md_file_t basis_file = {0};
		if (md_file_open(&basis_file, basis_filepath, MD_FILE_READ)) {
			MD_LOG_DEBUG("Attempting to parse VLX basis set from file: '" STR_FMT "'", STR_ARG(basis_filepath));
			md_buffered_reader_t basis_reader = md_buffered_reader_from_file(buf, cap, basis_file);
			bool parse_result = parse_basis_set(&vlx->basis_set, &basis_reader, vlx->arena);
			md_file_close(&basis_file);
			if (!parse_result) {
				MD_LOG_ERROR("An error occured when parsing the basis set for veloxchem data");
				goto done;
			}
			normalize_basis_set(&vlx->basis_set);
		} else {
			// Attempt to read basis set file from same folder as file
			str_t folder = { 0 };
			if (!extract_folder_path(&folder, filename)) {
				MD_LOG_ERROR("An error occured when extracting the path to supplied file");
				goto done;
			}
			md_strb_reset(&sb);
			md_strb_push_str(&sb, folder);
			md_strb_push_str(&sb, ident);
			basis_filepath = md_strb_to_str(sb);
			if (md_file_open(&basis_file, basis_filepath, MD_FILE_READ)) {
				MD_LOG_DEBUG("Attempting to parse VLX basis set from file: '" STR_FMT "'", STR_ARG(basis_filepath));
				md_buffered_reader_t basis_reader = md_buffered_reader_from_file(buf, cap, basis_file);
				bool parse_result = parse_basis_set(&vlx->basis_set, &basis_reader, vlx->arena);
				md_file_close(&basis_file);
				if (!parse_result) {
					MD_LOG_ERROR("An error occured when parsing the basis set for veloxchem data");
					goto done;
				}
				normalize_basis_set(&vlx->basis_set);
			}
			else {
                MD_LOG_ERROR("Could not find basis file corresponding to identifier: '" STR_FMT "'", STR_ARG(ident));
				goto done;
			}
		}
	}

	// Build the AO remap table and apply it to all loaded matrices.
	// This must happen after the basis set has been successfully resolved,
	// since build_ao_remap() requires basis topology to be valid.
	if (vlx->basis_set.atom_basis.count > 0) {
		size_t num_ao = 0;
		if (vlx->scf.alpha.density.data) {
			num_ao = vlx->scf.alpha.density.size[0];
		} else if (vlx->scf.S.data) {
			num_ao = vlx->scf.S.size[0];
		} else if (vlx->scf.alpha.coefficients.data) {
			num_ao = compute_basis_num_atomic_orbitals(vlx);
			if (num_ao == 0) {
				MD_LOG_ERROR("Unable to infer AO dimension for SCF coefficient normalization");
				goto done;
			}
		}
		// No AO indexed data in the file - a vib or opt only run, say - so there is nothing to
		// permute and no remap to build. Building one anyway used to fail the whole load, because
		// a zero length table can never match the basis set's AO count.
		if (num_ao > 0) {
			md_array_resize(vlx->ao_remap, num_ao, vlx->arena);
			if (!build_ao_remap(vlx->ao_remap, num_ao, vlx)) {
				MD_LOG_ERROR("Failed to build AO remap table");
				goto done;
			}
		}

		if (num_ao > 0 && vlx->ao_remap) {
			// Normalize SCF coefficients to canonical [num_mo x num_ao] in shell order.
			if (!normalize_orbital_coefficients(&vlx->scf.alpha, num_ao, vlx->ao_remap, "Alpha orbital")) {
				goto done;
			}
			if (vlx->scf.alpha.density.data && num_ao == vlx->scf.alpha.density.size[0]) {
				ao_permute_square(vlx->scf.alpha.density.data, num_ao, vlx->ao_remap);
			}
			if (vlx->scf.type == VLX_SCF_UNRESTRICTED) {
				if (!normalize_orbital_coefficients(&vlx->scf.beta, num_ao, vlx->ao_remap, "Beta orbital")) {
					goto done;
				}
				if (vlx->scf.beta.density.data && num_ao == vlx->scf.beta.density.size[0]) {
					ao_permute_square(vlx->scf.beta.density.data, num_ao, vlx->ao_remap);
				}
			}
			else {
				// memcpy again from alpha into beta as dims may have changed.
				MEMCPY(&vlx->scf.beta.coefficients, &vlx->scf.alpha.coefficients, sizeof(vlx_2d_data_t));
				MEMCPY(&vlx->scf.beta.density, &vlx->scf.alpha.density, sizeof(vlx_2d_data_t));
			}
			if (vlx->scf.S.data && num_ao == vlx->scf.S.size[0]) {
				ao_permute_square(vlx->scf.S.data, num_ao, vlx->ao_remap);
			}

			// Everything is now in shell order and no further AO-basis math is
			// performed, so this is the point to leave VeloxChem's pure/spherical
			// basis for the Cartesian one that md_gto_basis_t requires.
			if (!vlx_convert_ao_data_to_cartesian(vlx)) {
				MD_LOG_ERROR("Failed to convert AO data to the Cartesian basis");
				goto done;
			}
		}
	}

	if (!validate_scf_canonical_layout(vlx)) {
		goto done;
	}

	// Identify homo and lumo
	if (vlx->scf.alpha.occupancy.data) {
		for (size_t i = 0; i < vlx->scf.alpha.occupancy.size; ++i) {
			if (vlx->scf.alpha.occupancy.data[i] == 0.0) {
				vlx->scf.alpha.homo_idx = (size_t)MAX(0, (int64_t)i - 1);
				vlx->scf.alpha.lumo_idx = i;
				break;
			}
		}
	}

	if (vlx->scf.beta.occupancy.data) {
		for (size_t i = 0; i < vlx->scf.beta.occupancy.size; ++i) {
			if (vlx->scf.beta.occupancy.data[i] == 0.0) {
				vlx->scf.beta.homo_idx = (size_t)MAX(0, (int64_t)i - 1);
				vlx->scf.beta.lumo_idx = i;
				break;
			}
		}
	}

	// NOTE: the AO to atom map is not built or kept here. It is a pure function of the shell list,
	// so a consumer derives it from the published basis/shell attributes with
	// md_gto_basis_ao_to_atom - the same call this file used to make on its behalf.

	result = true;
done:
	md_temp_end(temp);

	return result;
}

static inline void extract_row(double* dst, const vlx_2d_data_t* data, size_t row_idx) {
        ASSERT(dst);
        ASSERT(data);
	size_t num_cols = data->size[1];
	for (size_t i = 0; i < num_cols; ++i) {
		dst[i] = data->data[row_idx * num_cols + i];
	}
}

static inline void extract_col(double* dst, const vlx_2d_data_t* data, size_t col_idx) {
	ASSERT(dst);
	ASSERT(data);
	ASSERT(col_idx < data->size[1]);

	for (size_t i = 0; i < data->size[0]; ++i) {
		dst[i] = data->data[i * data->size[1] + col_idx];
	}
}

static inline size_t number_of_molecular_orbitals(const vlx_orbital_t* orb) {
	ASSERT(orb);
	return orb->coefficients.size[0];
}

static inline size_t number_of_atomic_orbitals(const vlx_orbital_t* orb) {
	ASSERT(orb);
	return orb->coefficients.size[1];
}

static inline size_t number_of_ao_coefficients(const vlx_orbital_t* orb) {
	ASSERT(orb);
	return orb->coefficients.size[1];
}

static inline void extract_ao_coefficients(double* out_coeff, const vlx_orbital_t* orb, size_t ao_idx) {
	ASSERT(out_coeff);
	ASSERT(orb);
	ASSERT(ao_idx < number_of_atomic_orbitals(orb));

	extract_col(out_coeff, &orb->coefficients, ao_idx);
}
size_t vlx_rsp_number_of_excited_states(const vlx_t* vlx) {
	if (vlx) {
		if (vlx->rsp.type == VLX_RSP_LINEAR) {
			return vlx->rsp.number_of_frequencies;
		}
	}
	return 0;
}

const dvec3_t* vlx_rsp_electric_transition_dipole_moments(const vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.electric_transition_dipoles;
	}
	return NULL;
}

const dvec3_t* vlx_rsp_magnetic_transition_dipole_moments(const vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.magnetic_transition_dipoles;
	}
	return NULL;
}

const dvec3_t* vlx_rsp_velocity_transition_dipole_moments(const vlx_t* vlx) {
	if (vlx) {
		return vlx->rsp.velocity_transition_dipoles;
	}
	return NULL;
}

// RIXS
// @TODO: The RIXS datasets 'photon_energies', 'elastic_cross_sections', 'emission_energies',
// 'energy_losses' and the 'gamma_fwhm_ev' attribute are not yet read in vlx_parse_rsp().
// 'num_incomming_photons' and 'num_final_states' also need to be assigned from the dataset dims.
// The accessors below are thin and will start returning valid data as soon as that is in place.

bool vlx_rsp_has_nto(const vlx_t* vlx) {
	if (!vlx) return false;
	if (vlx->rsp.solution_matrix.data && vlx->rsp.solution_matrix.size[0] == vlx->rsp.number_of_frequencies) {
		return true;
	}
	return false;
}

size_t vlx_rsp_nto_lambdas_extract(double* out_lambdas, const vlx_t* vlx, size_t state_idx, size_t lambda_count) {
	return vlx_rsp_extract_nto(NULL, out_lambdas, vlx, state_idx, VLX_NTO_PARTICLE, lambda_count);
}

size_t vlx_rsp_nto_coefficients_extract(double* out_coefficients, double* out_lambdas, const vlx_t* vlx, size_t state_idx, vlx_nto_type_t type, size_t lambda_count) {
	return vlx_rsp_extract_nto(out_coefficients, out_lambdas, vlx, state_idx, type, lambda_count);
}

// OPT

vlx_t* vlx_create(md_allocator_i* backing, md_system_t* sys) {
	ASSERT(backing);
	md_allocator_i* arena = md_arena_allocator_create(backing, MEGABYTES(1));
	ASSERT(arena);
	vlx_t* vlx = md_alloc(arena, sizeof(vlx_t));
	if (!vlx) {
		MD_LOG_ERROR("Failed to allocate memory for veloxchem object");
		vlx->sys = sys;
	return vlx;
	}
	MEMSET(vlx, 0, sizeof(vlx_t));
	vlx->arena = arena;
	vlx->sys = sys;
	return vlx;
}

// XPS


// The centre of charge, in Angstrom: where a dipole moment is drawn from. Nuclear charge weighted
// positions less the electronic contribution, per electron. Returns false when the file does not
// carry what it takes to compute one, in which case no dipole group is published at all - half a
// group is not a dipole anyone can draw.
static bool vlx_centre_of_charge(dvec3_t* out_angstrom, const vlx_t* vlx) {
	const size_t   num_atoms     = vlx_number_of_atoms(vlx);
	const dvec3_t* atom_coord    = vlx_atom_coordinates(vlx);
	const uint8_t* atomic_number = vlx_atomic_numbers(vlx);

	if (num_atoms == 0 || !atom_coord || !atomic_number) {
		return false;
	}

	const size_t num_electrons = vlx_number_of_electrons(vlx, VLX_SPIN_ALPHA) + vlx_number_of_electrons(vlx, VLX_SPIN_BETA);
	if (num_electrons == 0) {
		return false;
	}

	// Coordinates are Angstrom while the moment is atomic units, so the nuclear term is taken to
	// bohr first and the result taken back at the end.
	double nx = 0.0, ny = 0.0, nz = 0.0;
	for (size_t i = 0; i < num_atoms; ++i) {
		const double z = (double)atomic_number[i];
		nx += atom_coord[i].x * ANGSTROM_TO_BOHR * z;
		ny += atom_coord[i].y * ANGSTROM_TO_BOHR * z;
		nz += atom_coord[i].z * ANGSTROM_TO_BOHR * z;
	}

	const dvec3_t moment = vlx_scf_ground_state_dipole_moment(vlx);
	const double  inv_ne = 1.0 / (double)num_electrons;

	out_angstrom->x = (nx - moment.x) * inv_ne * BOHR_TO_ANGSTROM;
	out_angstrom->y = (ny - moment.y) * inv_ne * BOHR_TO_ANGSTROM;
	out_angstrom->z = (nz - moment.z) * inv_ne * BOHR_TO_ANGSTROM;
	return true;
}

// Reconstructs one excited state's AO-basis transition density purely from attributes already on
// 'sys' - the response solution vectors and the alpha MO coefficients - so this runs with no vlx
// object in reach, and keeps working after one is torn down. 'type' is baked in by the three thin
// wrappers below, one per sibling path, the same shape as the MO coefficient split above.
//
// A whole (unsliced) request reconstructs every state one after another into 'dst'; expensive, but
// no more so than the vlx-based accessor doing the same loop, and slicing by state is how a caller
// avoids paying for states it does not need.
static size_t vlx_transition_density_provide(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data, vlx_transition_type_t type) {
	md_system_t* sys = (md_system_t*)user_data;

	const md_attribute_t* sol   = md_attributes_find(&sys->attributes, STR_LIT("vlx/rsp/solution_matrix"));
	const md_attribute_t* coeff = md_attributes_find(&sys->attributes, STR_LIT("orbital/alpha/coefficient"));
	const md_attribute_t* core  = md_attributes_find(&sys->attributes, STR_LIT("vlx/rsp/num_core"));
	const md_attribute_t* val   = md_attributes_find(&sys->attributes, STR_LIT("vlx/rsp/num_valence"));
	const md_attribute_t* vir   = md_attributes_find(&sys->attributes, STR_LIT("vlx/rsp/num_virtual"));
	if (!sol || !coeff || !core || !val || !vir) {
		MD_LOG_ERROR("'" STR_FMT "' is missing the response data it reconstructs from", STR_ARG(attr->path));
		return 0;
	}

	double num_core = 0.0, num_valence = 0.0, num_virtual = 0.0;
	md_attribute_extract_f64(&num_core,    1, core, md_unit_none());
	md_attribute_extract_f64(&num_valence, 1, val,  md_unit_none());
	md_attribute_extract_f64(&num_virtual, 1, vir,  md_unit_none());

	const size_t nocc       = num_core > 0.0 ? (size_t)num_core : (size_t)num_valence;
	const size_t nvir       = (size_t)num_virtual;
	const size_t num_states = sol->format.shape[0];
	const size_t vec_len    = sol->format.shape[1];
	const size_t num_ao     = coeff->format.shape[1];
	const size_t num_mo     = coeff->format.shape[0];

	// A provider that returns 0 is reported by the extract as "wrote 0 of N", which says nothing
	// about WHICH input was missing. Each of these is a distinct, actionable state, so each says so.
	if (nocc == 0 || nvir == 0 || num_ao == 0 || num_states == 0) {
		MD_LOG_ERROR("'" STR_FMT "': %zu occupied, %zu virtual, %zu atomic orbitals, %zu states - none of these may be zero",
			STR_ARG(attr->path), nocc, nvir, num_ao, num_states);
		return 0;
	}

	// The destination is sized by the CALLER's slice, and this writes num_ao^2 per state. Nothing
	// upstream guarantees the two agree - the shape was published from
	// vlx_scf_number_of_atomic_orbitals() while num_ao here comes off the coefficient matrix -
	// so disagreeing is an overrun, not a wrong picture. Check before writing a single value.
	const size_t states_written = (slice && slice->num_idx > 0) ? 1 : num_states;
	const size_t needed = states_written * num_ao * num_ao;
	if (needed != cap) {
		MD_LOG_ERROR("'" STR_FMT "': asked for %zu values, would write %zu (%zu state(s) x %zu x %zu atomic orbitals)",
			STR_ARG(attr->path), cap, needed, states_written, num_ao, num_ao);
		return 0;
	}

	md_temp_scope_t temp = md_temp_begin();
	size_t result = 0;

	double* coeff_data = md_temp_alloc_array(temp, double, num_mo * num_ao);
	double* row        = md_temp_alloc_array(temp, double, vec_len);
	if (!coeff_data || !row) {
		MD_LOG_ERROR("'" STR_FMT "': failed to allocate %zu doubles of scratch", STR_ARG(attr->path), num_mo * num_ao + vec_len);
		goto done;
	}

	// Every step below gets its own report. An '&&' chain here costs nothing to write and tells a
	// reader of the log only that the provider produced nothing, which is exactly the position this
	// was debugged from - and there are five separate attributes it can be let down by.
	if (md_attribute_extract_f64(coeff_data, num_mo * num_ao, coeff, md_unit_none()) != num_mo * num_ao) {
		MD_LOG_ERROR("'" STR_FMT "': could not read the %zu x %zu coefficients from '" STR_FMT "'",
			STR_ARG(attr->path), num_mo, num_ao, STR_ARG(coeff->path));
		goto done;
	}

	for (size_t s = 0; s < states_written; ++s) {
		const uint32_t state_idx = (slice && slice->num_idx > 0) ? slice->idx[0] : (uint32_t)s;
		const md_attribute_slice_t row_slice = md_attribute_slice_1(state_idx);

		if (md_attribute_extract_slice_f64(row, vec_len, sol, &row_slice, md_unit_none()) != vec_len) {
			MD_LOG_ERROR("'" STR_FMT "': could not read state %u's %zu element solution vector from '" STR_FMT "' (%u x %u)",
				STR_ARG(attr->path), state_idx, vec_len, STR_ARG(sol->path), sol->format.shape[0], sol->format.shape[1]);
			goto done;
		}
		if (!vlx_build_transition_density_matrix((double*)dst + s * num_ao * num_ao, row, vec_len, nocc, nvir, coeff_data, num_ao, type)) {
			MD_LOG_ERROR("'" STR_FMT "': could not build state %u's density from %zu occupied x %zu virtual over %zu atomic orbitals",
				STR_ARG(attr->path), state_idx, nocc, nvir, num_ao);
			goto done;
		}
	}
	result = cap;

done:
	md_temp_end(temp);
	return result;
}

static size_t vlx_transition_density_attachment_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
	return vlx_transition_density_provide(dst, cap, attr, slice, user_data, VLX_TRANSITION_ATTACHMENT);
}

static size_t vlx_transition_density_detachment_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
	return vlx_transition_density_provide(dst, cap, attr, slice, user_data, VLX_TRANSITION_DETACHMENT);
}

static size_t vlx_transition_density_difference_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
	return vlx_transition_density_provide(dst, cap, attr, slice, user_data, VLX_TRANSITION_DIFFERENCE);
}

// D[ao_i][ao_j] = sum_mo occ[mo] * C[mo][ao_i] * C[mo][ao_j] - the definition of the one particle
// density matrix in terms of the (possibly fractionally occupied) orbitals that produced it, so
// this reconstructs the ground state density EXACTLY from data already resident as attributes,
// with no separate density read needed. A zero occupied MO contributes nothing, so it is skipped
// rather than paid for.
static bool vlx_build_occupation_density_matrix(double* out_matrix, const double* coeff, const double* occ, size_t num_mo, size_t num_ao) {
	ASSERT(out_matrix);
	ASSERT(coeff);
	ASSERT(occ);

	if (num_mo == 0 || num_ao == 0) {
		return false;
	}

	MEMSET(out_matrix, 0, sizeof(double) * num_ao * num_ao);
	for (size_t mo = 0; mo < num_mo; ++mo) {
		const double w = occ[mo];
		if (w == 0.0) continue;
		const double* c = coeff + mo * num_ao;
		for (size_t ao_i = 0; ao_i < num_ao; ++ao_i) {
			const double wci = w * c[ao_i];
			if (wci == 0.0) continue;
			for (size_t ao_j = 0; ao_j < num_ao; ++ao_j) {
				out_matrix[ao_i * num_ao + ao_j] += wci * c[ao_j];
			}
		}
	}
	return true;
}

// Reconstructs one spin's ground state AO density from the MO coefficients and occupations this
// system already carries as attributes - the same 'no vlx object in reach' shape as the transition
// density providers above, sharing the same rationale: recomputing on demand costs one pass over
// the occupied orbitals instead of holding a second, redundant [A][A] copy alongside them.
static size_t vlx_scf_density_provide(void* dst, size_t cap, const md_attribute_t* attr, void* user_data, vlx_spin_t spin) {
	md_system_t* sys = (md_system_t*)user_data;

	str_t coeff_path = spin == VLX_SPIN_ALPHA ? STR_LIT("orbital/alpha/coefficient")        : STR_LIT("orbital/beta/coefficient");
	str_t occ_path   = spin == VLX_SPIN_ALPHA ? STR_LIT("vlx/scf/orbital/alpha/occupation") : STR_LIT("vlx/scf/orbital/beta/occupation");

	const md_attribute_t* coeff = md_attributes_find(&sys->attributes, coeff_path);
	const md_attribute_t* occ   = md_attributes_find(&sys->attributes, occ_path);
	if (!coeff || !occ) {
		MD_LOG_ERROR("'" STR_FMT "' is missing the orbital data it reconstructs from", STR_ARG(attr->path));
		return 0;
	}

	const size_t num_mo = coeff->format.shape[0];
	const size_t num_ao = coeff->format.shape[1];
	if (occ->format.shape[0] != num_mo) {
		MD_LOG_ERROR("'" STR_FMT "': occupation holds %u values, coefficients hold %zu orbitals", STR_ARG(attr->path), occ->format.shape[0], num_mo);
		return 0;
	}

	if (num_ao * num_ao != cap) {
		MD_LOG_ERROR("'" STR_FMT "' was asked for %zu values and its coefficients span %zu atomic orbitals", STR_ARG(attr->path), cap, num_ao);
		return 0;
	}

	md_temp_scope_t temp = md_temp_begin();
	double* coeff_data = md_temp_alloc_array(temp, double, num_mo * num_ao);
	double* occ_data   = md_temp_alloc_array(temp, double, num_mo);

	bool ok = coeff_data && occ_data;
	if (!ok) {
		MD_LOG_ERROR("'" STR_FMT "': failed to allocate %zu doubles of scratch", STR_ARG(attr->path), num_mo * num_ao + num_mo);
	}
	ok = ok && md_attribute_extract_f64(coeff_data, num_mo * num_ao, coeff, md_unit_none()) == num_mo * num_ao
	        && md_attribute_extract_f64(occ_data, num_mo, occ, md_unit_none()) == num_mo
	        && vlx_build_occupation_density_matrix((double*)dst, coeff_data, occ_data, num_mo, num_ao);

	md_temp_end(temp);
	return ok ? cap : 0;
}

static size_t vlx_scf_alpha_density_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
	(void)slice; // one whole {A,A} matrix, not indexed by anything a slice could fix
	return vlx_scf_density_provide(dst, cap, attr, user_data, VLX_SPIN_ALPHA);
}

static size_t vlx_scf_beta_density_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
	(void)slice;
	return vlx_scf_density_provide(dst, cap, attr, user_data, VLX_SPIN_BETA);
}

// alpha +/- beta. Both halves are themselves virtual, so this is a derivation over derivations -
// legal because the graph stays acyclic, and the case the acyclicity rule in md_system.h exists
// for. A restricted calculation gets this for free: beta is an ALIAS of alpha there, so the total
// comes out as twice alpha and the difference as zero without a special case anywhere.
//
// Each read rebuilds both matrices. That is the caller's to cache, and a representation already
// keys one on the volume hash it computes from its own settings.
static size_t vlx_scf_density_combine(void* dst, size_t cap, const md_attribute_t* attr, void* user_data, double beta_scale) {
	md_system_t* sys = (md_system_t*)user_data;

	const md_attribute_t* alpha = md_attributes_find(&sys->attributes, STR_LIT("orbital/alpha/density"));
	const md_attribute_t* beta  = md_attributes_find(&sys->attributes, STR_LIT("orbital/beta/density"));
	if (!alpha || !beta) {
		MD_LOG_ERROR("'" STR_FMT "' is missing a spin density it combines", STR_ARG(attr->path));
		return 0;
	}

	md_temp_scope_t temp = md_temp_begin();
	double* beta_data = md_temp_alloc_array(temp, double, cap);

	bool ok = beta_data != NULL;
	if (!ok) {
		MD_LOG_ERROR("'" STR_FMT "': failed to allocate %zu doubles of scratch", STR_ARG(attr->path), cap);
	}
	if (ok && md_attribute_extract_f64((double*)dst, cap, alpha, md_unit_none()) != cap) {
		MD_LOG_ERROR("'" STR_FMT "': could not read '" STR_FMT "'", STR_ARG(attr->path), STR_ARG(alpha->path));
		ok = false;
	}
	if (ok && md_attribute_extract_f64(beta_data, cap, beta, md_unit_none()) != cap) {
		MD_LOG_ERROR("'" STR_FMT "': could not read '" STR_FMT "'", STR_ARG(attr->path), STR_ARG(beta->path));
		ok = false;
	}

	if (ok) {
		double* out = (double*)dst;
		for (size_t i = 0; i < cap; ++i) {
			out[i] += beta_scale * beta_data[i];
		}
	}

	md_temp_end(temp);
	return ok ? cap : 0;
}

static size_t vlx_scf_total_density_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
	(void)slice; // one whole {A,A} matrix, not indexed by anything a slice could fix
	return vlx_scf_density_combine(dst, cap, attr, user_data, 1.0);
}

static size_t vlx_scf_difference_density_provider(void* dst, size_t cap, const md_attribute_t* attr, const md_attribute_slice_t* slice, void* user_data) {
	(void)slice;
	return vlx_scf_density_combine(dst, cap, attr, user_data, -1.0);
}

void vlx_publish_whole_file_attributes(md_system_t* sys, const vlx_t* vlx) {
	ASSERT(sys);

	if (!vlx) {
		return;
	}
	if (!sys->attributes.alloc) {
		MD_LOG_ERROR("Attribute table allocator not set; the system has not been initialised");
		return;
	}

	const md_unit_t hartree  = vlx_unit_hartree();
	const md_unit_t e_bohr   = vlx_unit_e_bohr();
	const md_unit_t bohr_magneton = vlx_unit_bohr_magneton();
	const md_unit_t bohr_velocity = vlx_unit_bohr_velocity();

	// A label is only carried where the leaf cannot spell it: a consumer prettifies the last path
	// segment when there is none, so "gradient_norm" needs no help and "ir_intensity" does.

	// ---- SCF: one value per molecular orbital, per spin. Sibling paths rather than one {2,M}
	// attribute, because a beta set is either present or absent and never an index to loop over.
	const size_t num_mos = vlx_scf_number_of_molecular_orbitals(vlx);
	const double* alpha_energy_data = vlx_scf_mo_energy(vlx, VLX_SPIN_ALPHA);
	const double* beta_energy_data  = vlx_scf_mo_energy(vlx, VLX_SPIN_BETA);
	const double* alpha_occ_data    = vlx_scf_mo_occupancy(vlx, VLX_SPIN_ALPHA);
	const double* beta_occ_data     = vlx_scf_mo_occupancy(vlx, VLX_SPIN_BETA);

	const md_attribute_id_t alpha_energy = vlx_publish_series(sys, STR_LIT("vlx/scf/orbital/alpha/energy"),     STR_LIT("Energy"),		hartree,        alpha_energy_data,	num_mos);
	const md_attribute_id_t alpha_occ    = vlx_publish_series(sys, STR_LIT("vlx/scf/orbital/alpha/occupation"), STR_LIT("Occupation"), md_unit_none(), alpha_occ_data,		num_mos);

	// A restricted calculation has one set of orbitals; beta is a second name for it, not a second
	// copy. The open shell case shares the energies and not the occupations, and falls out of the
	// same pointer test without being special cased.
	const md_attribute_id_t beta_energy  = vlx_publish_or_alias(sys, alpha_energy, STR_LIT("vlx/scf/orbital/beta/energy"),      STR_LIT("Energy"),		hartree,        alpha_energy_data, beta_energy_data, num_mos);
	const md_attribute_id_t beta_occ     = vlx_publish_or_alias(sys, alpha_occ,    STR_LIT("vlx/scf/orbital/beta/occupation"),  STR_LIT("Occupation"), md_unit_none(), alpha_occ_data,    beta_occ_data,    num_mos);

	// An orbital's coefficients and density are already published under the neutral orbital/ tree,
	// because md_gto_basis_t is mdlib's representation and not this program's. Its energy and
	// occupation are the same objects described the same way by every QM code, so leaving them
	// only under vlx/ splits one set of orbitals across two namespaces.
	//
	// Aliasing rather than moving them: a consumer already reading the vlx/ path keeps working,
	// and one written against the neutral name works too. That is the whole reason a second name
	// is cheaper than a rename.
	vlx_alias(sys, alpha_energy, STR_LIT("orbital/alpha/energy"));
	vlx_alias(sys, alpha_occ,    STR_LIT("orbital/alpha/occupation"));
	vlx_alias(sys, beta_energy,  STR_LIT("orbital/beta/energy"));
	vlx_alias(sys, beta_occ,     STR_LIT("orbital/beta/occupation"));

	// ---- The GTO basis and the MO coefficients: everything an evaluator needs to turn an orbital
	// into samples on a grid, with no reader in the loop.
	//
	// These land on FORMAT NEUTRAL paths, unlike the vlx/ tree above, because they are not this
	// program's output: md_gto_basis_t is mdlib's own normalised Cartesian representation, and
	// vlx_gto_basis_extract is the conversion into it. Another QM reader fills the same paths
	// with the same meaning, which is the test a path has to pass to lose its prefix.
	//
	// The shell list is published as four columns rather than as md_gto_shell_t records. A record
	// in the table would be a struct layout contract, additive changes to it would silently
	// invalidate anything stored, and nothing here would be self describing. The price is one
	// interleave in md_gto_basis_extract_attributes, paid when a consumer builds a basis rather
	// than per evaluation.
	{
		md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
		md_gto_basis_t basis = {0};

		if (vlx_gto_basis_extract(&basis, vlx, md_temp_allocator(temp))) {
			const size_t num_shells     = basis.num_shells;
			const size_t num_primitives = basis.num_primitives;
			const size_t shell_stride   = sizeof(md_gto_shell_t);

			vlx_publish_column(sys, STR_LIT("basis/shell/atom_index"),       STR_LIT("Atom Index"),       md_unit_none(), MD_ATTRIBUTE_TYPE_U32, &basis.shells->atom_idx,         shell_stride, num_shells);
			vlx_publish_column(sys, STR_LIT("basis/shell/primitive_offset"), STR_LIT("Primitive Offset"), md_unit_none(), MD_ATTRIBUTE_TYPE_U32, &basis.shells->primitive_offset, shell_stride, num_shells);
			vlx_publish_column(sys, STR_LIT("basis/shell/primitive_count"),  STR_LIT("Primitive Count"),  md_unit_none(), MD_ATTRIBUTE_TYPE_U32, &basis.shells->num_primitives,   shell_stride, num_shells);
			vlx_publish_column(sys, STR_LIT("basis/shell/angular_momentum"), STR_LIT("Angular Momentum"), md_unit_none(), MD_ATTRIBUTE_TYPE_U32, &basis.shells->l,                shell_stride, num_shells);

			// Exponents are bohr^-2 and the contraction coefficients carry the shell's radial
			// normalisation; the per monomial factor is applied at evaluation. See the AO
			// CONVENTION block in md_gto.h - these values only mean anything against it.
			const md_unit_t inv_bohr_sq = md_unit_pow(md_unit_bohr_radius(), -2);
			vlx_publish_column(sys, STR_LIT("basis/primitive/exponent"),    STR_LIT("Exponent"),	inv_bohr_sq,    MD_ATTRIBUTE_TYPE_F32, basis.alpha, sizeof(float), num_primitives);
			vlx_publish_column(sys, STR_LIT("basis/primitive/coefficient"), STR_LIT("Coefficient"), md_unit_none(), MD_ATTRIBUTE_TYPE_F32, basis.coeff, sizeof(float), num_primitives);

			// The MO matrix is stored [M][A] contiguously and is already Cartesian - the spherical
			// to Cartesian conversion happens once at parse time - so this is a straight copy into
			// a rank 2 attribute, and a consumer takes one orbital with a slice extract.
			//
			// f64 and not f32: md_gto takes AO coefficients as double to keep the QM code's
			// precision at the boundary, and there is no point publishing them already narrowed.
			const size_t num_ao = vlx_scf_number_of_atomic_orbitals(vlx);
			if (num_ao > 0 && num_mos > 0) {
				if (num_ao != md_gto_basis_num_ao(&basis)) {
					MD_LOG_ERROR("MO coefficients span %zu AOs but the basis has %zu; not publishing them", num_ao, md_gto_basis_num_ao(&basis));
				} else {
					md_attribute_format_t format = {
						.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2,
						.shape = { (uint32_t)num_mos, (uint32_t)num_ao },
					};
					const size_t byte_size = num_mos * num_ao * sizeof(double);

					// Row zero is the base of the whole matrix.
					const double* alpha_coeff = vlx_scf_mo_coefficients(vlx, 0, VLX_SPIN_ALPHA);
					const double* beta_coeff  = vlx_scf_mo_coefficients(vlx, 0, VLX_SPIN_BETA);

					md_attribute_id_t alpha_coeff_id = MD_ATTRIBUTE_INVALID;
					if (alpha_coeff) {
						alpha_coeff_id = vlx_publish(sys, STR_LIT("orbital/alpha/coefficient"),	STR_LIT("Alpha Coefficient"),	md_unit_none(), format, alpha_coeff, byte_size);
					}
					// Beta previously went unpublished whenever it shared alpha's buffer, so a
					// restricted calculation simply had no orbital/beta/coefficient for anyone to
					// read. A second name costs nothing and makes it present and correct.
					if (beta_coeff) {
						if (beta_coeff == alpha_coeff) {
							vlx_alias(sys, alpha_coeff_id, STR_LIT("orbital/beta/coefficient"));
						} else {
							vlx_publish(sys, STR_LIT("orbital/beta/coefficient"),	STR_LIT("Beta Coefficient"),	md_unit_none(), format, beta_coeff, byte_size);
						}
					}

					// The AO overlap S, {A,A} and symmetric. It belongs to the BASIS and not to a spin
					// channel: both channels and every density in this table are expressed against this
					// one metric, so basis/ and not orbital/.
					//
					// It is in the same AO order and convention as the coefficients above - converted by
					// vlx_cart_convert_square and permuted by ao_permute_square, exactly as they are - so
					// the two can be used together without further ceremony. Read the AO CONVENTION block
					// in md_gto.h first: this is the CARTESIAN overlap, and the Cartesian embedding of a
					// spherical basis is rank deficient, so for any file that stored spherical data this
					// matrix is SINGULAR. That is fine for what it is wanted for - Mulliken partitioning
					// and tr(DS) - and fatal for anything that inverts or factorises it.
					const double* overlap = vlx_scf_overlap_matrix_data(vlx);
					if (overlap && vlx_scf_overlap_matrix_size(vlx) == num_ao) {
						md_attribute_format_t overlap_format = {
							.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2,
							.shape = { (uint32_t)num_ao, (uint32_t)num_ao },
						};
						vlx_publish(sys, STR_LIT("basis/overlap"), STR_LIT("AO Overlap"), md_unit_none(),
									overlap_format, overlap, num_ao * num_ao * sizeof(double));
					}

					// Ground state densities: computed on demand from the coefficients and occupations
					// just published above rather than kept as a second resident [A][A] copy. See
					// vlx_build_occupation_density_matrix() for why this is exact, not approximate.
					md_attribute_format_t density_format = {
						.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2,
						.shape = { (uint32_t)num_ao, (uint32_t)num_ao },
					};
					md_attribute_id_t alpha_density_id = MD_ATTRIBUTE_INVALID;
					if (alpha_coeff) {
						md_attribute_virtual_t alpha_density_virt = { .provider = vlx_scf_alpha_density_provider, .user_data = sys };
						alpha_density_id = vlx_publish_virtual(sys, STR_LIT("orbital/alpha/density"), STR_LIT("Alpha Density"), md_unit_none(), density_format, &alpha_density_virt);
					}
					if (beta_coeff) {
						// A density is built from coefficients AND occupations, so it is the same
						// density only when BOTH are shared. Restricted open shell shares the
						// orbitals and not the occupations, so beta gets its own provider there -
						// which now works, because the coefficients it reads exist as an alias.
						const bool same_density = (beta_coeff == alpha_coeff) &&
							(vlx_scf_mo_occupancy(vlx, VLX_SPIN_BETA) == vlx_scf_mo_occupancy(vlx, VLX_SPIN_ALPHA));

						if (same_density) {
							vlx_alias(sys, alpha_density_id, STR_LIT("orbital/beta/density"));
						} else {
							md_attribute_virtual_t beta_density_virt = { .provider = vlx_scf_beta_density_provider, .user_data = sys };
							vlx_publish_virtual(sys, STR_LIT("orbital/beta/density"), STR_LIT("Beta Density"), md_unit_none(), density_format, &beta_density_virt);
						}

						// The two the UI offers besides the spins. Published rather than combined at
						// the point of use, so that "which densities does this system have" is one
						// question asked of the table and not four cases in a consumer.
						md_attribute_virtual_t total_density_virt = { .provider = vlx_scf_total_density_provider,      .user_data = sys };
						md_attribute_virtual_t diff_density_virt  = { .provider = vlx_scf_difference_density_provider, .user_data = sys };
						vlx_publish_virtual(sys, STR_LIT("orbital/total/density"),      STR_LIT("Total Density"),       md_unit_none(), density_format, &total_density_virt);
						vlx_publish_virtual(sys, STR_LIT("orbital/difference/density"), STR_LIT("Spin Difference Density"), md_unit_none(), density_format, &diff_density_virt);
					}
				}
			}
		}

		md_temp_end(temp);
	}

	// ---- Density properties: AO basis {A,A} matrices carried through from the file as they were
	// found. Unlike the SCF densities above these are not derived from anything else the table
	// holds - there is nothing to compute them from - so they are resident, one path per property.
	//
	// The path is built from the DATASET NAME, not from the label and not from the index. An index
	// silently re-points at a different property whenever the set changes across a reload, and a
	// label is display text that two datasets are free to share.
	for (size_t i = 0; i < md_array_size(vlx->density_properties); ++i) {
		const vlx_density_property_t* prop = &vlx->density_properties[i];
		if (!prop || !prop->data || prop->dim[0] == 0 || prop->dim[1] == 0) {
			continue;
		}

		str_t name = str_empty(prop->name) ? prop->label : prop->name;
		if (str_empty(name)) {
			continue;
		}

		char path_buf[256];
		str_t path = vlx_attribute_path(path_buf, sizeof(path_buf), STR_LIT("vlx/density_property"), name);
		if (str_empty(path)) {
			continue;
		}

		// The label is what the file called it for a human, the path is its identity. When they are
		// the same word there is nothing for the label to add, and an absent one is a valid state.
		str_t label = str_eq(prop->label, name) ? (str_t){0} : prop->label;

		md_attribute_format_t format = {
			.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2,
			.shape = { (uint32_t)prop->dim[0], (uint32_t)prop->dim[1] },
		};
		vlx_publish(sys, path, label, md_unit_none(), format, prop->data, prop->dim[0] * prop->dim[1] * sizeof(double));
	}

	// The number of EXCITED STATES, which only a linear response resolves.
	const size_t num_states = vlx_rsp_number_of_excited_states(vlx);

	// ---- RSP: the raw solution vectors and the occupied/virtual split they are indexed by. These
	// are not meant for direct consumption - a consumer wants the reconstructed density, not the
	// eigenvector it came from - they exist purely so vlx_transition_density_provide() below can
	// rebuild a transition density with no 'vlx' object in reach, which is the whole point of
	// storing them as attributes rather than reaching back into 'vlx' from the provider.
	if (vlx->rsp.solution_matrix.data && num_states > 0) {
		vlx_publish_matrix(sys, STR_LIT("vlx/rsp/solution_matrix"), (str_t){0}, md_unit_none(),
			vlx->rsp.solution_matrix.data, vlx->rsp.solution_matrix.size[0], vlx->rsp.solution_matrix.size[1]);
		vlx_publish_scalar(sys, STR_LIT("vlx/rsp/num_core"),    (str_t){0}, md_unit_none(), (double)vlx->rsp.num_core);
		vlx_publish_scalar(sys, STR_LIT("vlx/rsp/num_valence"), (str_t){0}, md_unit_none(), (double)vlx->rsp.num_valence);
		vlx_publish_scalar(sys, STR_LIT("vlx/rsp/num_virtual"), (str_t){0}, md_unit_none(), (double)vlx->rsp.num_virtual);

		// One virtual attribute per sibling density, each computed on demand from the solution
		// vectors and MO coefficients above. rank {S,A,A}: slice by state for one density, or take
		// the whole thing and pay for reconstructing every state - see the caveat on
		// vlx_transition_density_provide() about that cost.
		const size_t num_ao = vlx_scf_number_of_atomic_orbitals(vlx);
		if (num_ao > 0) {
			md_attribute_format_t density_format = {
				.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 3,
				.shape = { (uint32_t)num_states, (uint32_t)num_ao, (uint32_t)num_ao },
			};
			md_attribute_virtual_t attach_virt = { .provider = vlx_transition_density_attachment_provider, .user_data = sys };
			md_attribute_virtual_t detach_virt = { .provider = vlx_transition_density_detachment_provider, .user_data = sys };
			md_attribute_virtual_t diff_virt   = { .provider = vlx_transition_density_difference_provider, .user_data = sys };
			vlx_publish_virtual(sys, STR_LIT("vlx/rsp/transition_density/attachment"), STR_LIT("Attachment Density"), md_unit_none(), density_format, &attach_virt);
			vlx_publish_virtual(sys, STR_LIT("vlx/rsp/transition_density/detachment"), STR_LIT("Detachment Density"), md_unit_none(), density_format, &detach_virt);
			vlx_publish_virtual(sys, STR_LIT("vlx/rsp/transition_density/difference"), STR_LIT("Difference Density"), md_unit_none(), density_format, &diff_virt);
		}
	}

	// ---- NTO lambdas: the weights of the natural transition orbital pairs, one row per excited
	// state. The rows are RAGGED - a state has as many pairs as it has - so they are padded to the
	// widest row and the shape is {S,Lmax}. Zero is the honest pad here rather than a sentinel: a
	// lambda IS a weight, an absent pair carries none, and the consumer which already stops at a
	// 1e-3 cutoff stops in exactly the same place. Anything ragged enough that zero would be a real
	// value does not belong in a rectangular attribute at all.
	if (vlx_rsp_has_nto(vlx) && num_states > 0) {
		double row[VLX_NTO_MAX_LAMBDAS];
		size_t max_lambdas = 0;
		for (size_t s = 0; s < num_states; ++s) {
			const size_t n = vlx_rsp_nto_lambdas_extract(row, vlx, s, ARRAY_SIZE(row));
			max_lambdas = MAX(max_lambdas, n);
		}

		if (max_lambdas > 0) {
			md_attribute_format_t format = {
				.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 2,
				.shape = { (uint32_t)num_states, (uint32_t)max_lambdas },
			};
			md_attribute_id_t id = vlx_publish(sys, STR_LIT("vlx/rsp/nto/lambda"), STR_LIT("NTO Weight"), md_unit_none(), format, NULL, 0);
			double* dst = id != MD_ATTRIBUTE_INVALID ? (double*)md_attributes_data(&sys->attributes, id, MD_ATTRIBUTE_TYPE_F64) : NULL;

			if (dst) {
				// Created zeroed, so a short row needs nothing written past its own length.
				for (size_t s = 0; s < num_states; ++s) {
					const size_t n = vlx_rsp_nto_lambdas_extract(row, vlx, s, ARRAY_SIZE(row));
					MEMCPY(dst + s * max_lambdas, row, MIN(n, max_lambdas) * sizeof(double));
				}
			} else if (id != MD_ATTRIBUTE_INVALID) {
				md_attributes_remove(&sys->attributes, id);
			}
		}

		// The NTO coefficient vectors themselves, {S,Lmax,A}, padded on the lambda axis exactly as
		// the weights above so the two index the same space and a consumer moves one slider.
		//
		// RESIDENT rather than computed on demand, unlike the transition densities. The reason is
		// size, not principle: S x Lmax x A doubles is under a megabyte for a realistic file, while
		// a transition density is A^2 PER STATE. And it is what keeps this self contained - the NTO
		// math reads the vlx object, so a provider would have to close over the reader, which is
		// the one thing the port exists to avoid. Paying it once at load buys that outright.
		const size_t num_ao_nto = vlx_scf_number_of_atomic_orbitals(vlx);
		if (max_lambdas > 0 && num_ao_nto > 0) {
			const vlx_nto_type_t types[2] = { VLX_NTO_PARTICLE, VLX_NTO_HOLE };
			str_t paths[2] = { STR_LIT("vlx/rsp/nto/particle/coefficient"), STR_LIT("vlx/rsp/nto/hole/coefficient") };
			str_t labels[2] = { STR_LIT("Particle"), STR_LIT("Hole") };

			md_attribute_format_t format = {
				.type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 3,
				.shape = { (uint32_t)num_states, (uint32_t)max_lambdas, (uint32_t)num_ao_nto },
			};

			for (int t = 0; t < 2; ++t) {
				md_attribute_id_t id = vlx_publish(sys, paths[t], labels[t], md_unit_none(), format, NULL, 0);
				double* dst = id != MD_ATTRIBUTE_INVALID ? (double*)md_attributes_data(&sys->attributes, id, MD_ATTRIBUTE_TYPE_F64) : NULL;
				if (!dst) {
					if (id != MD_ATTRIBUTE_INVALID) md_attributes_remove(&sys->attributes, id);
					continue;
				}

				// One state at a time, straight into its own plane of the table's storage. The
				// extract writes [lambda_count][num_ao] which is exactly the plane's layout, and a
				// state with fewer pairs leaves the rest of its plane at the zero it was created
				// with.
				const size_t plane = max_lambdas * num_ao_nto;
				for (size_t s = 0; s < num_states; ++s) {
					vlx_rsp_nto_coefficients_extract(dst + s * plane, NULL, vlx, s, types[t], max_lambdas);
				}
			}
		}
	}

	// ---- Dipole groups. Every vector gets an origin beside it in the SAME shape, because a dipole
	// has no index space of its own to anchor it - see the ANCHORING note in md_system.h. Publishing
	// a vector without its origin would leave a group nothing can draw, so all of this is skipped
	// together when the centre of charge cannot be computed.
	dvec3_t origin = {0, 0, 0};
	if (vlx_centre_of_charge(&origin, vlx)) {
		const dvec3_t ground_state = vlx_scf_ground_state_dipole_moment(vlx);
		vlx_publish_vec3_series(sys, STR_LIT("dipole/ground_state/vector"), STR_LIT("Ground State"), e_bohr, &ground_state, 1);
		vlx_publish_origin(sys, STR_LIT("dipole/ground_state/origin"), origin);

		// One per excited state. The magnetic and velocity forms are different quantities in
		// different units, which is exactly why unit sits on the attribute and not on the group.
		if (num_states > 0) {
			const dvec3_t* electric = vlx_rsp_electric_transition_dipole_moments(vlx);
			const dvec3_t* magnetic = vlx_rsp_magnetic_transition_dipole_moments(vlx);
			const dvec3_t* velocity = vlx_rsp_velocity_transition_dipole_moments(vlx);

			if (electric) {
				vlx_publish_vec3_series(sys, STR_LIT("dipole/electric_transition/vector"), STR_LIT("Electric Transition"), e_bohr, electric, num_states);
				vlx_publish_origin(sys, STR_LIT("dipole/electric_transition/origin"), origin);
			}
			if (magnetic) {
				vlx_publish_vec3_series(sys, STR_LIT("dipole/magnetic_transition/vector"), STR_LIT("Magnetic Transition"), bohr_magneton, magnetic, num_states);
				vlx_publish_origin(sys, STR_LIT("dipole/magnetic_transition/origin"), origin);
			}
			if (velocity) {
				vlx_publish_vec3_series(sys, STR_LIT("dipole/velocity_transition/vector"), STR_LIT("Velocity Transition"), bohr_velocity, velocity, num_states);
				vlx_publish_origin(sys, STR_LIT("dipole/velocity_transition/origin"), origin);
			}
		}
	}
}

// The core block, published the moment the system exists to hold it. Everything here is a straight
// pass-through from the file - a number or a piece of text that no later step looks at again - so
// the struct fields behind it exist only to carry the values across md_system_reset(), which clears
// the table and therefore has to happen between reading them and publishing them.
static bool vlx_publish_core(const vlx_t* vlx) {
	md_system_t* sys = vlx->sys;
	ASSERT(sys);

	if (!sys->attributes.alloc) {
		MD_LOG_ERROR("Attribute table allocator not set; the system has not been initialised");
		return false;
	}

	// ---- Molecule level scalars. rank 0 is a single value, not an array of one. ----
	vlx_publish_scalar(sys, STR_LIT("vlx/molecular_charge"),          STR_LIT("Molecular Charge"),			md_unit_none(),		vlx->molecular_charge);
	vlx_publish_scalar(sys, STR_LIT("vlx/nuclear_repulsion_energy"),  STR_LIT("Nuclear Repulsion Energy"),	vlx_unit_hartree(),        vlx->nuclear_repulsion_energy);
	// The two facts about a calculation that are TEXT and nothing else - no consumer can derive them
	// from the columns, the way it can derive the SCF type from whether the spin channels share data.
	vlx_publish_scalar(sys, STR_LIT("vlx/spin_multiplicity"),         STR_LIT("Spin Multiplicity"),		md_unit_none(), (double)vlx->spin_multiplicity);

	// The electron counts per spin channel. Not derivable from the occupations with the precision
	// this states them: a fractional occupation sums to something a consumer then has to round, and
	// which way to round is exactly the question this answers. They are also what tr(D S) is checked
	// against, which is the one cheap test that the density and the overlap agree.
	vlx_publish_scalar(sys, STR_LIT("vlx/electron_count/alpha"),      STR_LIT("Alpha Electrons"),		md_unit_none(), (double)vlx->number_of_alpha_electrons);
	vlx_publish_scalar(sys, STR_LIT("vlx/electron_count/beta"),       STR_LIT("Beta Electrons"),		md_unit_none(), (double)vlx->number_of_beta_electrons);

	vlx_publish_str(sys, STR_LIT("vlx/basis_set"),      STR_LIT("Basis Set"),      vlx->basis_set_ident);
	vlx_publish_str(sys, STR_LIT("vlx/dft_functional"), STR_LIT("DFT Functional"), vlx->dft_func_label);

	// The embedding potential the run was given, verbatim. Published because it is part of what the
	// calculation WAS and nothing else in the table records it; absent for a run without one, which
	// vlx_publish_str already handles by publishing nothing.
	vlx_publish_str(sys, STR_LIT("vlx/potfile"), STR_LIT("Potential File"), vlx->potfile_text);

	// WHICH SCF this was. A consumer can guess from whether the two spin channels share their
	// coefficients and their occupations, and that guess is right until it meets a file where one of
	// the two is missing. The reader knows; this is it saying so, on the same terms as the response
	// and optimisation types below.
	vlx_publish_str(sys, STR_LIT("vlx/scf/type"), STR_LIT("SCF Type"), vlx_scf_type_str(vlx->scf.type));

	// ---- The QM ATOM DOMAIN. ----
	//
	// The atoms this calculation covered, in ITS order and at ITS geometry. That is not the
	// system's atom set: a calculation can cover part of a loaded system - a chromophore inside a
	// protein - so the two spaces differ in length AND in order, and qm/atom/system_index is the
	// only bridge between them.
	//
	// Its own prefix, and NOT atom/*, which is the SYSTEM's atom domain: a consumer walking atom/
	// has no idea quantum chemistry exists and would index a QM column by system atom - silently,
	// and wrongly, on exactly the subset case this domain exists for. atom_property_query already
	// filters by extent and component count, but that filter passes a scalar QM column whenever the
	// two atom counts happen to agree, which is not a distinction worth resting on.
	//
	// Not under basis/ either, though the shell list indexes this space: the basis is one thing
	// defined OVER these atoms, and so are the normal modes below. A nuclear coordinate is not a
	// property of a basis set, and a file can carry a geometry without carrying a basis at all.
	// basis/shell/atom_index is an index INTO qm/atom, which is the relationship stated plainly.
	{
		const size_t num_qm_atoms = vlx->number_of_atoms;
		const md_element_t* atomic_number = vlx->atomic_numbers;

		if (num_qm_atoms > 0 && atomic_number) {
			md_attribute_format_t format = {
				.type = MD_ATTRIBUTE_TYPE_U8, .components = 1, .rank = 1, .shape = { (uint32_t)num_qm_atoms },
			};
			vlx_publish(sys, STR_LIT("qm/atom/atomic_number"), STR_LIT("Atomic Number"), md_unit_none(),
						format, atomic_number, num_qm_atoms * sizeof(md_element_t));
		}

		// Angstrom, matching the system's own coordinates rather than the bohr the evaluator works
		// in: a consumer comparing this geometry against md_system_state_t should not have to convert
		// first. This is the geometry the CALCULATION was run at, which is not necessarily where the
		// system's atoms are now - a trajectory frame or an optimisation step moves them.
		vlx_publish_vec3_series(sys, STR_LIT("qm/atom/coordinate"), STR_LIT("Coordinate"), vlx_unit_angstrom(),
								vlx->atom_coordinates, num_qm_atoms);
	}

	return true;
}

// Builds the system's atoms and state from the core block, and from here on every reader publishes
// straight into sys->attributes. It has to happen HERE, between the core block and everything else:
// md_system_reset() clears the attribute table, so a system built after the blocks were read would
// throw away everything they published.
//
// A NULL state means this file is SUPPLEMENTING a system somebody else loaded - its atoms and its
// state belong to that loader and are left exactly as they are, and only the table grows.
static bool vlx_system_begin(vlx_t* vlx, md_system_state_t* state) {
	ASSERT(vlx);
	md_system_t* sys = vlx->sys;
	ASSERT(sys);

	if (vlx->number_of_atoms == 0) {
		MD_LOG_ERROR("The veloxchem file contains no atoms");
		return false;
	}

	if (!sys->alloc) {
		MD_LOG_ERROR("System allocator not set");
		return false;
	}

	if (!state) {
		// Supplementing: no reset, no atoms, but the file's own core values still belong in the
		// table beside whatever the first loader put there.
		return vlx_publish_core(vlx);
	}

	if (!state->alloc) {
		MD_LOG_ERROR("State allocator not set");
		return false;
	}
	md_system_reset(sys);
	md_system_state_init(state, vlx->number_of_atoms);

	size_t capacity = ROUND_UP(vlx->number_of_atoms, 16);

    md_array_resize(sys->atom.type_idx, capacity, sys->alloc);
    md_array_resize(sys->atom.flags,    capacity, sys->alloc);

	MEMSET(sys->atom.type_idx,  0, md_array_bytes(sys->atom.type_idx));
    MEMSET(sys->atom.flags,		0, md_array_bytes(sys->atom.flags));

    md_atom_type_find_or_add(&sys->atom.type, STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, sys->alloc);

	for (size_t i = 0; i < vlx->number_of_atoms; ++i) {
		state->x[i] = (float)vlx->atom_coordinates[i].x;
		state->y[i] = (float)vlx->atom_coordinates[i].y;
		state->z[i] = (float)vlx->atom_coordinates[i].z;
		
		md_atomic_number_t z = vlx->atomic_numbers[i];
		str_t sym  = md_atomic_number_symbol(z);
        float mass = md_atomic_number_mass(z);
		float radius = md_atomic_number_vdw_radius(z);
		uint32_t color = md_atomic_number_cpk_color(z);

		md_atom_type_idx_t type_idx = md_atom_type_find_or_add(&sys->atom.type, sym, z, mass, radius, color, 0, sys->alloc);
		sys->atom.type_idx[i] = type_idx;
	}

	sys->atom.count = vlx->number_of_atoms;
    state->num_atoms = sys->atom.count;

	return vlx_publish_core(vlx);
}


// Which system atom each QM atom is, or nothing at all when the two spaces coincide.
//
// The file cannot decide this for itself: the same h5 carries a local-to-global map whether it is
// opened standalone - where the system IS the QM atoms and the map must NOT be applied - or against
// a larger system, where it must. What resolves it is WHICH ENTRY POINT WAS CALLED, which is why
// this takes the answer as an argument rather than trying to work it out.
//
// Publishing nothing is not the same as leaving it alone: a stale map from a previous load would
// send every evaluation to the wrong atoms, so the standalone case actively removes it.
static void vlx_publish_atom_system_index(md_system_t* sys, const vlx_t* vlx, bool supplemental) {
	ASSERT(sys);

	const str_t path = STR_LIT("qm/atom/system_index");
	// Declared before the if: this is C, not the C++ side of the tree.
	const md_attribute_t* existing = md_attributes_find(&sys->attributes, path);
	if (existing) {
		md_attributes_remove(&sys->attributes, existing->id);
	}
	if (!supplemental) {
		return;
	}

	const size_t num_qm_atoms = vlx_number_of_atoms(vlx);
	const int* local_to_global = vlx_local_to_global_atom_idx(vlx);
	if (num_qm_atoms == 0 || !local_to_global) {
		return;
	}

	md_attribute_format_t format = {
		.type = MD_ATTRIBUTE_TYPE_U32, .components = 1, .rank = 1, .shape = { (uint32_t)num_qm_atoms },
	};
	md_attribute_id_t id = md_attributes_create(&sys->attributes, &(md_attribute_desc_t){
		.path   = path,
			.format = format,
			.unit   = md_unit_none(),
			.label  = STR_LIT("System Atom Index"),
	});

	uint32_t* dst = (uint32_t*)md_attributes_data(&sys->attributes, id, MD_ATTRIBUTE_TYPE_U32);
	if (!dst) {
		if (id != MD_ATTRIBUTE_INVALID) md_attributes_remove(&sys->attributes, id);
		return;
	}
	for (size_t i = 0; i < num_qm_atoms; ++i) {
		dst[i] = (uint32_t)local_to_global[i];
	}
}

bool md_vlx_system_init_from_file(md_system_t* sys, struct md_system_state_t* state, str_t filename) {
	ASSERT(sys);

	md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);
	md_allocator_i* temp_arena = md_temp_allocator(temp_scope);
	vlx_t* vlx = vlx_create(temp_arena, sys);

	// The system is built from the core block and everything after it publishes into the table as it
	// is read - see vlx_system_begin. That is what makes the table a property of the LOAD: a system
	// carries its data the moment it is loaded, and a consumer finds it by asking the system. It
	// used to be published by the veloxchem UI component, which meant the data existed only because
	// that component was compiled in and had parsed the same file a second time.
	bool success = vlx_parse_file(vlx, filename, state);
	if (success) {
		vlx_publish_whole_file_attributes(sys, vlx);
		// Standalone: the system IS the QM atoms, so the map is cleared rather than written.
		vlx_publish_atom_system_index(sys, vlx, false);
	}

	md_temp_end(temp_scope);
	return success;
}

bool md_vlx_system_supplement_from_file(md_system_t* sys, str_t filename) {
	ASSERT(sys);

	if (!sys->attributes.alloc) {
		MD_LOG_ERROR("Cannot supplement a system which has not been initialised");
		return false;
	}

	md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);
	md_allocator_i* temp_arena = md_temp_allocator(temp_scope);
	vlx_t* vlx = vlx_create(temp_arena, sys);

	// A NULL state is what tells vlx_system_begin this file is supplementing: the atoms and the
	// state belong to whatever loaded the system first and are left alone, and this file only adds
	// to its table. Whether the file actually belongs to this system is
	// md_vlx_system_is_file_supplemental's question and the caller has already asked it.
	bool success = vlx_parse_file(vlx, filename, NULL);
	if (success) {
		vlx_publish_whole_file_attributes(sys, vlx);
		vlx_publish_atom_system_index(sys, vlx, true);
	}

	md_temp_end(temp_scope);
	return success;
}

// Attempt to open the file and check if it can supplement the existing system with QM data
bool md_vlx_system_is_file_supplemental(const md_system_t* sys, str_t filename) {
	ASSERT(sys);

	// Simple check here, we just check for the existence of a couple of fields in h5 file.

	str_t ext = {0};
	if (extract_ext(&ext, filename)) {
        if (!str_eq_ignore_case(ext, STR_LIT("h5")) && !str_eq_ignore_case(ext, STR_LIT("hdf5"))) {
			// Unsupported file extension
			return false;
		}
	}

	// Ensure a zero terminated string for interfacing to HDF5
	char buf[2048];
	str_copy_to_char_buf(buf, sizeof(buf), filename);

	h5_error_scope_t h5_err_scope = h5_error_scope_begin();

	// Open an existing file
	hid_t file_id = H5Fopen(buf, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id == H5I_INVALID_HID) {
		MD_LOG_ERROR("Could not open HDF5 file: '"STR_FMT"'", STR_ARG(filename));
		h5_error_scope_end(h5_err_scope);
		return false;
	}

	bool result = false;
	
    if (h5_check_dataset_exists(file_id, "qm_atom_indices") &&
		h5_check_dataset_exists(file_id, "nuclear_charges"))
	{
		md_temp_scope_t temp_scope = md_temp_begin();
		md_allocator_i* temp_arena = md_temp_allocator(temp_scope);
		// Derive atomic numbers from nuclear charges and verify with qm_atom_indices that they are consistent with the loaded system.
        // This is a heuristic check, but it should be sufficient to determine if the file contains data that can supplement the existing system.
		
		size_t vlx_num_atoms = 0;
		if (!h5_read_scalar(&vlx_num_atoms, file_id, H5T_NATIVE_UINT64, "number_of_atoms")) {
			goto temp_done;
		}

		md_array(int) qm_atom_indices = md_array_create(int, vlx_num_atoms, temp_arena);
		md_array(int) nuclear_charges = md_array_create(int, vlx_num_atoms, temp_arena);
		
        if (h5_read_dataset_data(qm_atom_indices, vlx_num_atoms, file_id, H5T_NATIVE_INT32, "qm_atom_indices") &&
			h5_read_dataset_data(nuclear_charges, vlx_num_atoms, file_id, H5T_NATIVE_INT32, "nuclear_charges"))
		{
			bool match = true;
			for (size_t i = 0; i < vlx_num_atoms; ++i) {
                int idx = qm_atom_indices[i];
				int z = md_atom_atomic_number(&sys->atom, idx);
                if (z != nuclear_charges[i]) {
					match = false;
					break;
				}
			}
			result = match;
		}
	temp_done:
		md_temp_end(temp_scope);
	}

	H5Fclose(file_id);
	h5_error_scope_end(h5_err_scope);

	return result;
}

// Externally visible procedures

size_t vlx_number_of_atoms(const vlx_t* vlx) {
	if (vlx) return vlx->number_of_atoms;
	return 0;
}

size_t vlx_number_of_electrons(const vlx_t* vlx, vlx_spin_t spin) {
	if (vlx) {
		if (spin == VLX_SPIN_ALPHA) {
			return vlx->number_of_alpha_electrons;
		} else if (spin == VLX_SPIN_BETA) {
			return vlx->number_of_beta_electrons;
		}
	}
	return 0;
}

const dvec3_t* vlx_atom_coordinates(const vlx_t* vlx) {
	if (vlx) return vlx->atom_coordinates;
	return NULL;
}

const uint8_t* vlx_atomic_numbers(const vlx_t* vlx) {
	if (vlx) return vlx->atomic_numbers;
	return NULL;
}

const int* vlx_local_to_global_atom_idx(const vlx_t* vlx) {
	if (vlx) return vlx->local_to_global_atom_idx;
	return NULL;
}

dvec3_t vlx_scf_ground_state_dipole_moment(const vlx_t* vlx) {
	if (vlx) return vlx->scf.ground_state_dipole_moment;
	return (dvec3_t){0};
}

size_t vlx_scf_number_of_atomic_orbitals(const vlx_t* vlx) {
	if (vlx) {
		return number_of_atomic_orbitals(&vlx->scf.alpha);
	}
	return 0;
}

size_t vlx_scf_number_of_molecular_orbitals(const vlx_t* vlx) {
	if (vlx) {
		return number_of_molecular_orbitals(&vlx->scf.alpha);
	}
	return 0;
}

const double* vlx_scf_mo_occupancy(const vlx_t* vlx, vlx_spin_t type) {
	if (vlx) {
		if (type == VLX_SPIN_ALPHA) {
			return vlx->scf.alpha.occupancy.data;
		} 
		else if (type == VLX_SPIN_BETA) {
			return vlx->scf.beta.occupancy.data;
		}
	}
	return NULL;
}

const double* vlx_scf_mo_energy(const vlx_t* vlx, vlx_spin_t type) {
	if (vlx) {
		if (type == VLX_SPIN_ALPHA) {
			return vlx->scf.alpha.energy.data;
		}
		else if (type == VLX_SPIN_BETA) {
			return vlx->scf.beta.energy.data;
		}
	}
	return NULL;
}

bool vlx_gto_basis_extract(md_gto_basis_t* out, const vlx_t* vlx, md_allocator_i* alloc) {
	if (!vlx || !out) return false;
	ASSERT(alloc);

	MEMSET(out, 0, sizeof(*out));

	int natoms   = (int)vlx->number_of_atoms;
	int max_angl = compute_max_angular_momentum(&vlx->basis_set, vlx->atomic_numbers, vlx->number_of_atoms);

	// Emit one shell per contracted function, ordered angl -> atom -> func.
	// MO coefficient vectors stored in vlx are already permuted to this order
	// by build_ao_remap() applied at parse time.
	for (int angl = 0; angl <= max_angl; angl++) {
		for (int atomidx = 0; atomidx < natoms; atomidx++) {
			int idelem = vlx->atomic_numbers[atomidx];
			basis_func_t basis_funcs[128];
			size_t num_basis_funcs = basis_set_extract_atomic_basis_func_angl(
				basis_funcs, ARRAY_SIZE(basis_funcs), &vlx->basis_set, idelem, angl);

			for (size_t funcidx = 0; funcidx < num_basis_funcs; funcidx++) {
				const basis_func_t* bf = &basis_funcs[funcidx];
				md_gto_shell_t shell = {
					.atom_idx         = (uint32_t)atomidx,
					.primitive_offset = out->num_primitives,
					.num_primitives   = (uint32_t)bf->count,
					.l                = (uint32_t)angl,
				};
				md_array_push(out->shells, shell, alloc);
				out->num_shells++;

				for (int ip = 0; ip < bf->count; ip++) {
					md_array_push(out->alpha, (float)bf->exponents[ip], alloc);
					md_array_push(out->coeff, (float)bf->normalization_coefficients[ip], alloc);
					out->num_primitives++;
				}
			}
		}
	}
	return out->num_shells > 0;
}

// Returns a direct pointer to the AO coefficient vector for MO mo_idx.
// The matrix is stored [num_mo][num_ao] after permutation and transpose at load time,
// so each MO's coefficients are contiguous and in shell order.
const double* vlx_scf_mo_coefficients(const vlx_t* vlx, size_t mo_idx, vlx_spin_t spin) {
	if (!vlx) return NULL;
	const vlx_orbital_t* orb = (spin == VLX_SPIN_ALPHA) ? &vlx->scf.alpha :
								  (spin == VLX_SPIN_BETA)  ? &vlx->scf.beta  : NULL;
	if (!orb || !orb->coefficients.data) return NULL;
	size_t num_mo = orb->coefficients.size[0];
	size_t num_ao = orb->coefficients.size[1];
	if (mo_idx >= num_mo || num_ao == 0) return NULL;
	return orb->coefficients.data + mo_idx * num_ao;
}

static inline size_t get_matrix_index(size_t i, size_t j, size_t N) {
	size_t row = (i < j) ? i : j;
	size_t col = (i < j) ? j : i;
	size_t row_offset = row * (2 * N - row + 1) / 2;
	return row_offset + (col - row);
}

// The overlap matrix is a square, symmetric matrix [N][N], this returns the length N
size_t  vlx_scf_overlap_matrix_size(const struct vlx_t* vlx) {
	if (vlx) {
		return vlx->scf.S.size[0];
	}
	return 0;
}

const double* vlx_scf_overlap_matrix_data(const struct vlx_t* vlx) {
	if (vlx) {
		return vlx->scf.S.data;
	}
	return NULL;
}
