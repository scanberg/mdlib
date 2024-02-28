#include <md_vlx.h>

#include <core/md_os.h>
#include <core/md_log.h>
#include <core/md_parse.h>
#include <core/md_arena_allocator.h>
#include <core/md_str_builder.h>
#include <md_util.h>
#include <md_molecule.h>

#include <hdf5.h>
#include <hdf5_hl.h>

enum {
	VLX_FLAG_GEOM  = 1,
	VLX_FLAG_BASIS = 2,
	VLX_FLAG_SCF   = 4,
	VLX_FLAG_RSP   = 8,
	VLX_FLAG_ALL   = (1+2+4+8),
};

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
	str_t label;
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
	// 0 would be a NULL entry, 1 corresponds to Hydrogen, 2 corresponds to Helium etc.
	struct {
		size_t count;
		basis_set_basis_t* data;
	} atom_basis;
} basis_set_t;

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

typedef struct vlx_molecule_t {
	size_t num_atoms;
	md_element_t* atomic_number;
	double* coord_x;
	double* coord_y;
	double* coord_z;
} vlx_molecule_t;

static basis_set_basis_t* basis_set_get_atom_basis(const basis_set_t* basis_set, int atomic_number) {
	if (atomic_number < basis_set->atom_basis.count) {
		return basis_set->atom_basis.data + atomic_number;
	}
	return NULL;
}

static int compute_max_angular_momentum(const basis_set_t* basis_set, const md_element_t* atomic_numbers, size_t count) {
	int max_angl = 0;
	for (size_t i = 0; i < count; ++i) {
		const basis_set_basis_t* atom_basis = basis_set_get_atom_basis(basis_set, atomic_numbers[i]);
		max_angl = MAX(max_angl, (int)atom_basis->max_type);
	}
	return max_angl;
}

static inline int spherical_momentum_num_components(int angl) {
	return angl * 2 + 1;
}

static inline int spherical_momentum_num_factors(int angl, int isph) {
	static const uint8_t table[5][9] = {
		{1},
		{1,1,1},
		{1,1,3,1,2},
		{2,1,3,3,3,2,2},
		{2,2,3,3,6,3,4,2,3},
	};
	return table[angl][isph];
}

#define d3  3.464101615137754587
#define f5  1.581138830084189666
#define f15 3.872983346207416885
#define f3  1.224744871391589049
#define g35 4.0 * 5.916079783099616042
#define g17 4.0 * 4.183300132670377739
#define g5  4.0 * 2.236067977499789696
#define g2  4.0 * 1.581138830084189666

static const double		S_factors[] = {1.0};
static const uint8_t    S_indices[] = {0};

static const double		P_factors[] = {1.0, 1.0, 1.0};
static const uint8_t	P_indices[] = {1, 2, 0};
static const uint8_t	P_offsets[] = {0, 1, 2};

static const double		D_factors[] = {d3, d3, -1.0, -1.0, 2.0, d3, 0.5 * d3, -0.5 * d3};
static const uint8_t    D_indices[] = {1, 4, 0, 3, 5, 2, 0, 3};
static const uint8_t    D_offsets[] = {0, 1, 2, 5, 6};

static const double		F_factors[] = {3.0 * f5, -f5, 2.0 * f15, 4.0 * f3, -f3, -f3, 2.0, -3.0, -3.0, 4.0 * f3, -f3, -f3, f15, -f15, f5, -3.0 * f5};
static const uint8_t    F_indices[] = {1, 6, 4, 8, 1, 6, 9, 2, 7, 5, 0, 3, 2, 7, 0, 3};
static const uint8_t	F_offsets[] = {0, 2, 3, 6, 9, 12, 14};

static const double		G_factors[] = {
	g35, -g35, 3.0 * g17, -g17, 6.0 * g5, -g5, -g5, 4.0 * g2, -3.0 * g2, -3.0 * g2,
	8.0, 3.0, 3.0, 6.0, -24.0, -24.0, 4.0 * g2, -3.0 * g2, -3.0 * g2, 3.0 * g5,
	-3.0 * g5, -0.5 * g5, 0.5 * g5,  g17,  -3.0 * g17, 0.25 * g35, 0.25 * g35, -1.50 * g35};
static const uint8_t    G_indices[] = {1, 6, 4, 11, 8, 1, 6, 13, 4, 11, 14, 0, 10, 3, 5, 12, 9, 2, 7, 5, 12, 0, 10, 2, 7, 0, 10, 3};
static const uint8_t	G_offsets[] = {0, 2, 4, 7, 10, 16, 19, 23, 25};

#undef d3
#undef f5
#undef f15
#undef f3
#undef g35
#undef g17
#undef g5 
#undef g2 

static inline const double* spherical_momentum_factors(int angl, int isph) {
	switch(angl) {
	case 0:	return S_factors;
	case 1:	return P_factors + P_offsets[isph];
	case 2: return D_factors + D_offsets[isph];
	case 3:	return F_factors + F_offsets[isph];
	case 4: return G_factors + G_offsets[isph];
	default:
		ASSERT(false);
		return NULL;
	}
}

static inline const uint8_t* spherical_momentum_indices(int angl, int isph) {
	switch(angl) {
	case 0: return S_indices;
	case 1:	return P_indices + P_offsets[isph];
	case 2: return D_indices + D_offsets[isph];
	case 3:	return F_indices + F_offsets[isph];
	case 4: return G_indices + G_offsets[isph];
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

typedef struct basis_func_range_t {
	int beg;
	int end;
} basis_func_range_t;

static basis_func_range_t basis_get_atomic_angl_basis_func_range(const basis_set_t* basis_set, int atomic_number, int angl) {
	basis_set_basis_t* atom_basis = basis_set_get_atom_basis(basis_set, atomic_number);
	ASSERT(atom_basis);

	basis_func_range_t range = {0};

	int beg = atom_basis->basis_func_offset;
	int end = atom_basis->basis_func_offset + atom_basis->basis_func_count;
	for (int i = beg; i < end; ++i) {
		basis_set_func_t basis_func = basis_set->basis_func.data[i];
		if (basis_func.type == angl) {
			range.beg = (range.end == 0) ? i : range.beg;
			range.end = i + 1;
		}
	}

	return range;
}

typedef struct basis_func_t {
	int type;
	int count;
	const double* exponents;
	const double* normalization_coefficients;
} basis_func_t;

static basis_func_t get_basis_func(const basis_set_t* basis_set, int basis_func_idx) {
	basis_set_func_t func = basis_set->basis_func.data[basis_func_idx];
	return (basis_func_t) {
		.type = func.type,
		.count = func.param_count,
		.exponents = basis_set->param.exponents + func.param_offset,
		.normalization_coefficients = basis_set->param.normalization_coefficients + func.param_offset,
	};
}

typedef struct pgto_t {
	float alpha, c_coeff, mo_coeff;
	float Nijk;
	int i, j, k, l;
	vec3_t coord;
	float cutoff;
} pgto_t;

size_t extract_pgtos(pgto_t* out_pgtos, size_t cap, const vlx_molecule_t* mol, const basis_set_t* basis_set, const double* mo_coeff, size_t num_mo_coeff) {

	size_t count_outer = 0;
	size_t count_inner = 0;

	for (size_t atom_idx = 0; atom_idx < mol->num_atoms; ++atom_idx) {
		int atomic_number = mol->atomic_number[atom_idx];
		basis_set_basis_t* atom_basis = basis_set_get_atom_basis(basis_set, atomic_number);
		
		size_t func_beg = atom_basis->basis_func_offset;
		size_t func_end = atom_basis->basis_func_offset + atom_basis->basis_func_count;
		for (size_t func_idx = func_beg; func_idx < func_end; ++func_idx) {
			basis_func_t basis_func = get_basis_func(basis_set, func_idx);

			int angl = basis_func.type;
			int nsph = spherical_momentum_num_components(angl);
			const lmn_t* lmn = cartesian_angular_momentum(angl);

			for (int isph = 0; isph < nsph; isph++) {
				// prepare Cartesian components (Maximum number of components should be 6 here for the currently supported basis set)
				double fcarts[8];
				double lx[8];
				double ly[8];
				double lz[8];
				int			      ncomp = spherical_momentum_num_factors(angl, isph);
				const double*	factors = spherical_momentum_factors(angl, isph);
				const uint8_t*	indices = spherical_momentum_indices(angl, isph);

				for (int icomp = 0; icomp < ncomp; icomp++) {
					fcarts[icomp] = factors[icomp];
					int cartind = indices[icomp];

					lx[icomp] = lmn[cartind][0];
					ly[icomp] = lmn[cartind][1];
					lz[icomp] = lmn[cartind][2];
				}

				// process primitives
				int nprims = (int)basis_func.count;
				const double* exponents = basis_func.exponents;
				const double* normcoefs = basis_func.normalization_coefficients;

				for (int iprim = 0; iprim < nprims; iprim++) {
					double expon = exponents[iprim];
					double coef1 = normcoefs[iprim];

					//if (expon < 1.0e-10) continue;
					count_outer += 1;

					// transform from Cartesian to spherical harmonics
					for (int icomp = 0; icomp < ncomp; icomp++) {

						pgto_t pgto;
						pgto.i = lx[icomp];
						pgto.j = ly[icomp];
						pgto.k = lz[icomp];
						pgto.l = angl;

						pgto.alpha = expon;
						pgto.c_coeff = coef1;
						pgto.mo_coeff = 0;

						count_inner += 1;
						//double coef2 = coef1 * fcarts[icomp];
						//double powxyz = pow(rx, lx[icomp]) * pow(ry, ly[icomp]) * pow(rz, lz[icomp]);
						//phiao += coef2 * powxyz * expon;
					}
				}
			}
		}
	}

	while(0) {};
}

static size_t compPhiAtomicOrbitals(double* out_phi, size_t phi_cap, const vlx_molecule_t* molecule,
	const basis_set_t* basis_set,
	const double xp,
	const double yp,
	const double zp)
{
	int natoms = (int)molecule->num_atoms;
	int max_angl = compute_max_angular_momentum(basis_set, molecule->atomic_number, molecule->num_atoms);

	size_t count = 0;

	// azimuthal quantum number: s,p,d,f,...
	for (int aoidx = 0, angl = 0; angl <= max_angl; angl++) {
		//CSphericalMomentum sphmom(angl);
		int nsph = spherical_momentum_num_components(angl);
		const lmn_t* lmn = cartesian_angular_momentum(angl);
		// magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...
		for (int isph = 0; isph < nsph; isph++) {
			// prepare Cartesian components (Maximum number of components should be 6 here for the currently supported basis set)
			double fcarts[8];
			double lx[8];
			double ly[8];
			double lz[8];
			int			      ncomp = spherical_momentum_num_factors(angl, isph);
			const double*	factors = spherical_momentum_factors(angl, isph);
			const uint8_t*	indices = spherical_momentum_indices(angl, isph);

			for (int icomp = 0; icomp < ncomp; icomp++) {
				fcarts[icomp] = factors[icomp];
				int cartind = indices[icomp];

				lx[icomp] = lmn[cartind][0];
				ly[icomp] = lmn[cartind][1];
				lz[icomp] = lmn[cartind][2];
			}

			// go through atoms

			for (int atomidx = 0; atomidx < natoms; atomidx++) {
				// process coordinates
				double rx = xp - molecule->coord_x[atomidx];
				double ry = yp - molecule->coord_y[atomidx];
				double rz = zp - molecule->coord_y[atomidx];
				double r2 = rx * rx + ry * ry + rz * rz;

				// process atomic orbitals
				int idelem = molecule->atomic_number[atomidx];

				basis_func_range_t range = basis_get_atomic_angl_basis_func_range(basis_set, idelem, angl);
				for (int funcidx = range.beg; funcidx < range.end; funcidx++, aoidx++) {
					double phiao = 0.0;

					basis_func_t basis_func = get_basis_func(basis_set, funcidx);

					// process primitives
					int nprims = (int)basis_func.count;
					const double* exponents = basis_func.exponents;
					const double* normcoefs = basis_func.normalization_coefficients;

					for (int iprim = 0; iprim < nprims; iprim++) {
						double expon = exp(-exponents[iprim] * r2);
						double coef1 = normcoefs[iprim];

						if (expon < 1.0e-10) continue;

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

static bool parse_basis_set(basis_set_t* basis_set, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	MEMSET(basis_set, 0, sizeof(basis_set_t));

	str_t line;
	str_t tok[4];

	// Insert null_basis element for index 0
	const basis_set_basis_t null_basis = {0};

	basis_set_basis_t* curr_atom_basis = NULL;
	while (md_buffered_reader_extract_line(&line, reader)) {
		size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (!num_tok) continue;

		if (num_tok == 2 && str_eq(tok[0], STR_LIT("@BASIS_SET"))) {
			basis_set->label = str_copy(tok[1], alloc);
		}
		else if (num_tok == 2 && str_eq(tok[0], STR_LIT("@ATOMBASIS"))) {
			int atomic_number = md_util_element_lookup_ignore_case(tok[1]);
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
			while (md_array_size(basis_set->atom_basis.data) < atomic_number) {
				md_array_push(basis_set->atom_basis.data, null_basis, alloc);
			}

			curr_atom_basis = md_array_push(basis_set->atom_basis.data, atom_basis, alloc);
			basis_set->atom_basis.count = md_array_size(basis_set->atom_basis.data);
		}
		else if (num_tok == 1 && str_eq(tok[0], STR_LIT("@END"))) {
			curr_atom_basis = NULL;
		}
		else if (num_tok == 3) {
			int type = char_to_angular_momentum_type(tok[0].ptr[0]);
			if (type == -1) {
				MD_LOG_ERROR("Unrecognized angular momentum type '" STR_FMT "' in basis set", STR_ARG(tok[0]));
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
				float exponent = (float)parse_float(tok[0]);
				float coeff    = (float)parse_float(tok[1]);
				
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

static bool parse_vlx_geom(md_vlx_geom_t* geom, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	md_buffered_reader_skip_line(reader); // ====...
	md_buffered_reader_skip_line(reader); // *empty*
	md_buffered_reader_skip_line(reader); // Atom  Coordinate X  Coordinate Y  Coordinate Z
	md_buffered_reader_skip_line(reader); // *empty* 

	str_t line;
	str_t tok[8];
	size_t count = 0;
	while (md_buffered_reader_extract_line(&line, reader)) {
		size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok == 4) {
			md_label_t sym = make_label(tok[0]);
			md_element_t nr = md_util_element_lookup_ignore_case(tok[0]);
			double x = parse_float(tok[1]);
			double y = parse_float(tok[2]);
			double z = parse_float(tok[3]);

			md_array_push(geom->atom_symbol,  sym, alloc);
			md_array_push(geom->atomic_number, nr, alloc);
			md_array_push(geom->coord_x,		x, alloc);
			md_array_push(geom->coord_y,		y, alloc);
			md_array_push(geom->coord_z,		z, alloc);
			count += 1;
		} else if (num_tok == 0) {
			// Assume valid end here upon empty line
			break;
		} else {
			MD_LOG_ERROR("Unexpected number of tokens in geometry section, expected 4, got (%zu)");
			return false;
		}
	}

	if (count == 0) {
		MD_LOG_ERROR("No atomic coordinates found");
		return false;
	}

	// If we end up here, we expect to read the following lines in order
	// Molecular charge            : int                                                                 
	// Spin multiplicity           : int                                                               
	// Number of atoms             : int                                                               
	// Number of alpha electrons   : int                                                               
	// Number of beta  electrons   : int   
	str_t field_ident[] = {
		STR_LIT("Molecular charge            :"),                                                              
		STR_LIT("Spin multiplicity           :"),                                                            
		STR_LIT("Number of atoms             :"),                                                            
		STR_LIT("Number of alpha electrons   :"),                                                            
		STR_LIT("Number of beta  electrons   :"),
	};
	int64_t field_vals[ARRAY_SIZE(field_ident)];

	size_t loc;
	for (size_t i = 0; i < ARRAY_SIZE(field_ident); ++i) {
		str_t ident = field_ident[i];
		if (!md_buffered_reader_extract_line(&line, reader) || !str_find_str(&loc, line, ident)) {
			MD_LOG_ERROR("Failed to parse line: '"STR_FMT"'", STR_ARG(field_ident[i]));
			return false;
		}
		str_t value_str = str_trim(str_substr(line, loc + str_len(ident), SIZE_MAX));
		field_vals[i] = parse_int(value_str);
	}

	geom->molecular_charge		= (int)field_vals[0];
	geom->spin_multiplicity		= (int)field_vals[1];
	geom->num_atoms				= (size_t)field_vals[2];
	geom->num_alpha_electrons	= (size_t)field_vals[3];
	geom->num_beta_electrons	= (size_t)field_vals[4];

	if (geom->num_atoms != count) {
		MD_LOG_ERROR("Incorrect number of atoms parsed, expected (%zu) entries, parsed (%zu).", geom->num_atoms, count);
		return false;
	}

	return true;
}

static bool parse_vlx_basis(md_vlx_basis_t* basis, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	md_buffered_reader_skip_line(reader); // ====...
	md_buffered_reader_skip_line(reader); // *empty*

	str_t line;
	str_t tok[8];
	int mask = 0;
	while (md_buffered_reader_extract_line(&line, reader)) {
		size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
		if (num_tok == 2 && str_eq(tok[0], STR_LIT("Basis:"))) {
			basis->ident = str_copy(tok[1], alloc);
			mask |= 1;
		} else if (num_tok == 5 && str_eq(tok[3], STR_LIT(":"))) {
			str_t first = str_join(tok[0], tok[2]);
			if (str_eq(first, STR_LIT("Contracted Basis Functions"))) {
				basis->num_contracted_basis_functions = parse_int(tok[4]);
				mask |= 2;
			} else if (str_eq(first, STR_LIT("Primitive Basis Functions"))) {
				basis->num_primitive_basis_functions = parse_int(tok[4]);
				mask |= 4;
				break;
			}
		} else if (num_tok == 1 && str_begins_with(tok[0], STR_LIT("====="))) {
			// Parsed into next section >.<
			return false;
		}
	}

	return mask == 7;
}

static bool parse_vlx_scf(md_vlx_scf_t* scf, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	str_t tok[8];
	str_t line;

	md_buffered_reader_skip_line(reader); // =====

	int mask = 0;
	while (md_buffered_reader_extract_line(&line, reader)) {
		line = str_trim(line);
		if (str_begins_with(line, STR_LIT("Iter. |")) && str_ends_with(line, STR_LIT("| Density Change"))) {
			md_buffered_reader_skip_line(reader); // -----
			// Parse table, start tokenization
			while (md_buffered_reader_extract_line(&line, reader)) {
				size_t num_tok = extract_tokens(tok, ARRAY_SIZE(tok), &line);
				if (num_tok == 6) {
					int    iteration		= (int)parse_int(tok[0]);
					double energy_tot		= parse_float(tok[1]);
					double energy_change	= parse_float(tok[2]);
					double gradient_norm	= parse_float(tok[3]);
					double max_gradient		= parse_float(tok[4]);
					double density_change	= parse_float(tok[5]);

					md_array_push(scf->iter.iteration, iteration, alloc);
					md_array_push(scf->iter.energy_total,  energy_tot, alloc);
					md_array_push(scf->iter.energy_change, energy_change, alloc);
					md_array_push(scf->iter.gradient_norm, gradient_norm, alloc);
					md_array_push(scf->iter.max_gradient, max_gradient, alloc);
					md_array_push(scf->iter.density_change, density_change, alloc);
					scf->iter.count += 1;

					mask |= 1;
				} else if (num_tok == 0) {
					// Assume valid end here upon empty line
					break;
				} else {
					MD_LOG_ERROR("Unexpected number of tokens in scf energy iteration section, expected 6, got (%zu)");
					return false;
				}
			}
		} else if (str_begins_with(line, STR_LIT("Total Energy")) && extract_tokens(tok, ARRAY_SIZE(tok), &line) == 5) {
			scf->total_energy = parse_float(tok[3]);
			mask |= 2;
		} else if (str_begins_with(line, STR_LIT("Electronic Energy")) && extract_tokens(tok, ARRAY_SIZE(tok), &line) == 5) {
			scf->electronic_energy = parse_float(tok[3]);
			mask |= 4;
		} else if (str_begins_with(line, STR_LIT("Nuclear Repulsion Energy")) && extract_tokens(tok, ARRAY_SIZE(tok), &line) == 6) {
			scf->nuclear_repulsion_energy = parse_float(tok[4]);
			mask |= 8;
		} else if (str_begins_with(line, STR_LIT("Gradient Norm")) && extract_tokens(tok, ARRAY_SIZE(tok), &line) == 5) {
			scf->gradient_norm = parse_float(tok[3]);
			mask |= 16;
		} else if (str_eq(line, STR_LIT("Ground State Dipole Moment"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_buffered_reader_skip_line(reader); // *empty*
			double vec[3];
			for (int i = 0; i < 3; ++i) {
				if (!md_buffered_reader_extract_line(&line, reader) || extract_tokens(tok, ARRAY_SIZE(tok), &line) != 6) {
					MD_LOG_ERROR("Failed to parse SCF Ground State Dipole Moment, incomplete fields!");
					return false;
				}
				vec[i] = parse_float(tok[2]);
			}
			scf->ground_state_dipole_moment.ident = STR_LIT("Ground State");
			scf->ground_state_dipole_moment.x = vec[0];
			scf->ground_state_dipole_moment.y = vec[1];
			scf->ground_state_dipole_moment.z = vec[2];
			mask |= 32;
		} else if (str_begins_with(line, STR_LIT("====="))) {
			// We've read too far and into the next section
			MD_LOG_ERROR("Failed to parse SCF section, some fields are missing");
			return false;
		}
		if (mask == 63) {
			return true;
		}
	}

	return false;
}

static bool parse_vlx_rsp_dipole_moments(md_vlx_dipole_moment_t* moments, size_t num_excited_states, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	str_t tok[8];
	str_t line;

	for (size_t i = 0; i < num_excited_states; ++i) {
		if (!md_buffered_reader_extract_line(&line, reader) || extract_tokens(tok, ARRAY_SIZE(tok), &line) != 6) {
			MD_LOG_ERROR("Unexpected number of tokens in rsp dipole moment table");
			return false;
		}
		
		str_t first = str_join(tok[0], tok[1]);
		if (!str_eq(first, STR_LIT("Excited State"))) {
			MD_LOG_ERROR("Unexpected entry in rsp dipole moment table, expected 'Excited State', got '"STR_FMT"'", first);
			return false;
		}

		str_t ident = tok[2];
		if (ident.len > 0 && ident.ptr[ident.len-1] == ':') {
			ident.len -= 1;
		}

		moments[i].ident = str_copy(ident, alloc);
		moments[i].x = parse_float(tok[3]);
		moments[i].y = parse_float(tok[4]);
		moments[i].z = parse_float(tok[5]);
	}

	return true;
}

static bool parse_vlx_rsp(md_vlx_rsp_t* rsp, md_buffered_reader_t* reader, md_allocator_i* alloc) {
	str_t tok[10];
	str_t line;

	md_buffered_reader_skip_line(reader); // =====

	int mask = 0;
	while ( md_buffered_reader_extract_line(&line, reader)) {
		line = str_trim(line);
		if (str_empty(line)) continue;

		if (str_begins_with(line, STR_LIT("Number of States")) && extract_tokens(tok, ARRAY_SIZE(tok), &line) == 5) {
			rsp->num_excited_states = parse_int(tok[4]);
			mask |= 1;
		} else if (str_eq(line, STR_LIT("Electric Transition Dipole Moments (dipole length, a.u.)"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_buffered_reader_skip_line(reader); // X Y Z
			md_array_resize(rsp->electronic_transition_length, rsp->num_excited_states, alloc);
			if (!parse_vlx_rsp_dipole_moments(rsp->electronic_transition_length, rsp->num_excited_states, reader, alloc)) {
				return false;
			}
			mask |= 2;
		} else if (str_eq(line, STR_LIT("Electric Transition Dipole Moments (dipole velocity, a.u.)"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_buffered_reader_skip_line(reader); // X Y Z
			md_array_resize(rsp->electronic_transition_velocity, rsp->num_excited_states, alloc);
			if (!parse_vlx_rsp_dipole_moments(rsp->electronic_transition_velocity, rsp->num_excited_states, reader, alloc)) {
				return false;
			}
			mask |= 4;
		} else if (str_eq(line, STR_LIT("Magnetic Transition Dipole Moments (a.u.)"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_buffered_reader_skip_line(reader); // X Y Z
			md_array_resize(rsp->magnetic_transition, rsp->num_excited_states, alloc);
			if (!parse_vlx_rsp_dipole_moments(rsp->magnetic_transition, rsp->num_excited_states, reader, alloc)) {
				return false;
			}
			mask |= 8;
		} else if (str_eq(line, STR_LIT("One-Photon Absorption"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_array_resize(rsp->absorption_ev, rsp->num_excited_states, alloc);
			md_array_resize(rsp->absorption_osc_str, rsp->num_excited_states, alloc);
			for (size_t i = 0; i < rsp->num_excited_states; ++i) {
				if (!md_buffered_reader_extract_line(&line, reader) || !extract_tokens(tok, ARRAY_SIZE(tok), &line) == 9) {
					MD_LOG_ERROR("Unexpected number of tokens in entry when parsing One-Photon Absorption");
					return false;
				}
				rsp->absorption_ev[i] = parse_float(tok[5]);
				rsp->absorption_osc_str[i] = parse_float(tok[8]);
			}
			mask |= 16;
		} else if (str_eq(line, STR_LIT("Electronic Circular Dichroism"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_array_resize(rsp->electronic_circular_dichroism_cgs, rsp->num_excited_states, alloc);
			for (size_t i = 0; i < rsp->num_excited_states; ++i) {
				if (!md_buffered_reader_extract_line(&line, reader) || !extract_tokens(tok, ARRAY_SIZE(tok), &line) == 9) {
					MD_LOG_ERROR("Unexpected number of tokens in entry when parsing Electronic Circular Dichroism");
					return false;
				}
				rsp->electronic_circular_dichroism_cgs[i] = parse_float(tok[7]);
			}
			mask |= 32;
		} else if (str_begins_with(line, STR_LIT("===="))) {
			// Parsed into next section
			return false;
		}
		if (mask == 63) {
			return true;
		}
	}
	return false;
}

bool vlx_parse(md_vlx_data_t* vlx, md_buffered_reader_t* reader, uint32_t flags) {
	str_t line;
	uint32_t mask = 0;
	while (md_buffered_reader_extract_line(&line, reader)) {
		str_t str = str_trim(line);
		if ((flags & VLX_FLAG_GEOM) && str_eq(str, STR_LIT("Molecular Geometry (Angstroms)"))) {
			if (!parse_vlx_geom(&vlx->geom, reader, vlx->alloc)) {
				MD_LOG_ERROR("Failed to parse geometry");
				return false;
			}
			mask |= VLX_FLAG_GEOM;
		} else if ((flags & VLX_FLAG_BASIS) && str_eq(str, STR_LIT("Molecular Basis (Atomic Basis)"))) {
			if (!parse_vlx_basis(&vlx->basis, reader, vlx->alloc)) {
				MD_LOG_ERROR("Failed to parse basis");
				return false;
			}
			mask |= VLX_FLAG_BASIS;
		} else if ((flags & VLX_FLAG_SCF) && str_eq(str, STR_LIT("Self Consistent Field Driver Setup"))) {
			if (!parse_vlx_scf(&vlx->scf, reader, vlx->alloc)) {
				MD_LOG_ERROR("Failed to parse SCF section");
				return false;
			}
			mask |= VLX_FLAG_SCF;
		} else if ((flags & VLX_FLAG_RSP) && str_eq(str, STR_LIT("Linear Response EigenSolver Setup"))) {
			if (!parse_vlx_rsp(&vlx->rsp, reader, vlx->alloc)) {
				MD_LOG_ERROR("Failed to parse RSP section");
				return false;
			}
			mask |= VLX_FLAG_RSP;
		}
	}

	return true;
}

bool vlx_load_scf_h5_data(md_vlx_data_t* vlx, str_t filename) {
	hid_t  file_id = 0, dataset_id = 0, space_id = 0; /* identifiers */
	herr_t status = 0;

	bool result = false;

	/* Open an existing file. */
	file_id = H5Fopen(str_beg(filename), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id > 0) {
		dataset_id = H5Dopen(file_id, "/alpha_orbitals", H5P_DEFAULT);
		if (dataset_id > 0) {
			space_id = H5Dget_space(dataset_id);
			if (space_id > 0) {
				int ndim = H5Sget_simple_extent_ndims(space_id);
				if (ndim != 2) {
					MD_LOG_ERROR("Unexpected number of dimensions in scf orbitals");
					goto done_orb;
				}

				hsize_t dims[2];
				ndim = H5Sget_simple_extent_dims(space_id, dims, 0);
				if (ndim != 2) {
					MD_LOG_ERROR("Unexpected number of dimensions in scf orbitals");
					goto done_orb;
				}		

				double* data = md_alloc(vlx->alloc, sizeof(double) * dims[0] * dims[1]);
				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

				if (status != 0) {
					MD_LOG_ERROR("Could not read scf orbitals");
					goto done_orb;
				}

				vlx->scf.alpha.orbitals.dim[0] = dims[0];
				vlx->scf.alpha.orbitals.dim[1] = dims[1];
				vlx->scf.alpha.orbitals.data = data;

			done_orb:
				H5Sclose(space_id);
			}
			H5Dclose(dataset_id);
		}

		dataset_id = H5Dopen(file_id, "/alpha_energies", H5P_DEFAULT);
		if (dataset_id > 0) {
			space_id = H5Dget_space(dataset_id);
			if (space_id > 0) {
				int ndim = H5Sget_simple_extent_ndims(space_id);
				if (ndim != 1) {
					MD_LOG_ERROR("Unexpected number of dimensions in scf energies");
					goto done_ener;
				}

				hsize_t dim;
				ndim = H5Sget_simple_extent_dims(space_id, &dim, 0);
				if (ndim != 1) {
					MD_LOG_ERROR("Unexpected number of dimensions in scf energies");
					goto done_ener;
				}		

				double* data = md_alloc(vlx->alloc, sizeof(double) * dim);
				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

				if (status != 0) {
					MD_LOG_ERROR("Could not read scf energies");
					goto done_ener;
				}

				vlx->scf.alpha.energies.count = dim;
				vlx->scf.alpha.energies.data = data;

				done_ener:
				H5Sclose(space_id);
			}
			H5Dclose(dataset_id);
		}

		dataset_id = H5Dopen(file_id, "/alpha_occupations", H5P_DEFAULT);
		if (dataset_id > 0) {
			space_id = H5Dget_space(dataset_id);
			if (space_id > 0) {
				int ndim = H5Sget_simple_extent_ndims(space_id);
				if (ndim != 1) {
					MD_LOG_ERROR("Unexpected number of dimensions in scf occupations");
					goto done_occ;
				}

				hsize_t dim;
				ndim = H5Sget_simple_extent_dims(space_id, &dim, 0);
				if (ndim != 1) {
					MD_LOG_ERROR("Unexpected number of dimensions in scf occupations");
					goto done_occ;
				}		

				double* data = md_alloc(vlx->alloc, sizeof(double) * dim);
				status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

				if (status != 0) {
					MD_LOG_ERROR("Could not read scf occupations");
					goto done_occ;
				}

				vlx->scf.alpha.occupations.count = dim;
				vlx->scf.alpha.occupations.data = data;
			done_occ:
				H5Sclose(space_id);
			}
		}
		status = H5Fclose(file_id);
	}

	result = vlx->scf.alpha.orbitals.dim[0] != 0 && vlx->scf.alpha.energies.count != 0 && vlx->scf.alpha.occupations.count != 0;

	return result;
}

void md_vlx_data_free(md_vlx_data_t* data) {
	ASSERT(data);
	if (data->alloc) {
		md_arena_allocator_destroy(data->alloc);
	}
	MEMSET(data, 0, sizeof(md_vlx_data_t));
}

static bool vlx_data_parse_str(md_vlx_data_t* vlx, str_t str, md_allocator_i* alloc, uint32_t flags) {
	MEMSET(vlx, 0, sizeof(md_vlx_data_t));
	vlx->alloc = md_arena_allocator_create(alloc, MEGABYTES(1));
	md_buffered_reader_t reader = md_buffered_reader_from_str(str);
	return vlx_parse(vlx, &reader, flags);
}

bool md_vlx_data_parse_str(md_vlx_data_t* vlx, str_t str, md_allocator_i* alloc) {
	return vlx_data_parse_str(vlx, str, alloc, VLX_FLAG_ALL);
}

static bool vlx_data_parse_file(md_vlx_data_t* vlx, str_t filename, md_allocator_i* alloc, uint32_t flags) {
	md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
	if (!file) {
		MD_LOG_ERROR("Failed to open file: '"STR_FMT"'", STR_ARG(filename));
		return false;
	}

	MEMSET(vlx, 0, sizeof(md_vlx_data_t));
	vlx->alloc = md_arena_allocator_create(alloc, MEGABYTES(1));

	size_t cap = KILOBYTES(16);
	char* buf = md_temp_push(cap);
	md_buffered_reader_t reader = md_buffered_reader_from_file(buf, cap, file);

	bool result = vlx_parse(vlx, &reader, flags);
	md_file_close(file);

	if (flags & VLX_FLAG_BASIS) {
		if (!str_empty(vlx->basis.ident)) {
			md_strb_t sb = md_strb_create(md_temp_allocator);
			md_strb_fmt(&sb, "%s/" STR_FMT, MD_VLX_BASIS_FOLDER, STR_ARG(vlx->basis.ident));
			str_t basis_path = md_strb_to_str(sb);
			md_file_o* basis_file = md_file_open(basis_path, MD_FILE_READ | MD_FILE_BINARY);
			if (basis_file) {
				md_buffered_reader_t basis_reader = md_buffered_reader_from_file(buf, cap, file);
				vlx->basis.basis_set = md_alloc(vlx->alloc, sizeof(basis_set_t));
				if (!parse_basis_set(vlx->basis.basis_set, &basis_reader, vlx->alloc)) {
					MD_LOG_ERROR("An error occured when parsing the basis set for veloxchem data");
					return false;
				}

				md_file_close(basis_file);
			}
		}


		vlx_molecule_t mol = {
			.num_atoms = vlx->geom.num_atoms,
			.atomic_number = vlx->geom.atomic_number,
			.coord_x = vlx->geom.coord_x,
			.coord_y = vlx->geom.coord_y,
			.coord_z = vlx->geom.coord_z,
		};

		size_t num_pgtos = extract_pgtos(NULL, 0, &mol, vlx->basis.basis_set, 0, 0);
		
#if 0
		size_t phi_cap = vlx->basis.num_contracted_basis_functions;
		double* phis = md_alloc(md_temp_allocator, phi_cap * sizeof(double));
		md_timestamp_t t0 = md_time_current();
		const size_t dim = 80;
		for (size_t iz = 0; iz < dim; iz++) {
			double z = iz * 0.025;
			for (size_t iy = 0; iy < dim; ++iy) {
				double y = iy * 0.025;
				for (size_t ix = 0; ix < dim; ++ix) {
					double x = ix * 0.025;
					size_t num_phi = compPhiAtomicOrbitals(phis, phi_cap, &mol, vlx->basis.basis_set, x, y, z);
				}
			}
		}
		md_timestamp_t t1 = md_time_current();
		double dt = md_time_as_milliseconds(t1-t0);
		MD_LOG_INFO("Time taken to compute phi atomic orbitals for grid size [%zu][%zu][%zu]: %f", dim, dim, dim, dt);
#endif
	}

	if (flags & VLX_FLAG_SCF) {
		str_t file_path;
		if (result == true && extract_file_path_without_ext(&file_path, filename)) {
			md_strb_t sb = md_strb_create(md_temp_allocator);
			md_strb_push_str(&sb, file_path);
			md_strb_push_str(&sb, STR_LIT(".scf.h5"));
			str_t scf_path = md_strb_to_str(sb);
			if (md_path_is_valid(scf_path)) {
				if (!vlx_load_scf_h5_data(vlx, scf_path)) {
					MD_LOG_ERROR("Failed to load scf h5 parameters from file '" STR_FMT "'", STR_ARG(scf_path));
					return false;
				}
			}
		}
	}

	return result;
}

bool md_vlx_data_parse_file(md_vlx_data_t* vlx, str_t filename, md_allocator_i* alloc) {
	return vlx_data_parse_file(vlx, filename, alloc, VLX_FLAG_ALL);
}

bool md_vlx_molecule_init(md_molecule_t* mol, const md_vlx_data_t* vlx, md_allocator_i* alloc) {
	ASSERT(mol);
	ASSERT(vlx);

	if (vlx->geom.num_atoms == 0) {
		MD_LOG_ERROR("The veloxchem data object contains no atoms");
		return false;
	}

	mol->atom.count = vlx->geom.num_atoms;
	md_array_resize(mol->atom.type, mol->atom.count, alloc);
	md_array_resize(mol->atom.x, mol->atom.count, alloc);
	md_array_resize(mol->atom.y, mol->atom.count, alloc);
	md_array_resize(mol->atom.z, mol->atom.count, alloc);

	for (size_t i = 0; i < vlx->geom.num_atoms; ++i) {
		mol->atom.type[i] = vlx->geom.atom_symbol[i];
		mol->atom.x[i] = (float)vlx->geom.coord_x[i];
		mol->atom.y[i] = (float)vlx->geom.coord_y[i];
		mol->atom.z[i] = (float)vlx->geom.coord_z[i];
	}

	return true;
}

static bool vlx_mol_init_from_str(md_molecule_t* mol, str_t str, const void* arg, md_allocator_i* alloc) {
	(void)arg;
	md_vlx_data_t vlx = {0};
	bool success = false;
	if (vlx_data_parse_str(&vlx, str, md_heap_allocator, VLX_FLAG_GEOM)) {
		success = md_vlx_molecule_init(mol, &vlx, alloc);
	}
	md_vlx_data_free(&vlx);

	return success;
}

static bool vlx_mol_init_from_file(md_molecule_t* mol, str_t filename, const void* arg, md_allocator_i* alloc) {
	(void)arg;
	md_vlx_data_t vlx = {0};
	bool success = false;
	if (vlx_data_parse_file(&vlx, filename, md_heap_allocator, VLX_FLAG_GEOM)) {
		success = md_vlx_molecule_init(mol, &vlx, alloc);
	}
	md_vlx_data_free(&vlx);

	return success;
}

static md_molecule_loader_i vlx_loader = {
	vlx_mol_init_from_str,
	vlx_mol_init_from_file
};

md_molecule_loader_i* md_vlx_molecule_api() {
	return &vlx_loader;
}
