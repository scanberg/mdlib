#include <md_vlx.h>

#include <core/md_os.h>
#include <core/md_log.h>
#include <core/md_parse.h>
#include <core/md_arena_allocator.h>
#include <core/md_str_builder.h>

#include <hdf5.h>
#include <hdf5_hl.h>

static bool parse_geom(md_vlx_geom_t* geom, md_buffered_reader_t* reader, md_allocator_i* alloc) {
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
			double x = parse_float(tok[1]);
			double y = parse_float(tok[2]);
			double z = parse_float(tok[3]);

			md_array_push(geom->atom_symbol,  sym, alloc);
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

static bool parse_basis(md_vlx_basis_t* basis, md_buffered_reader_t* reader, md_allocator_i* alloc) {
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

static bool parse_scf(md_vlx_scf_t* scf, md_buffered_reader_t* reader, md_allocator_i* alloc) {
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

static bool parse_rsp_dipole_moments(md_vlx_dipole_moment_t* moments, size_t num_excited_states, md_buffered_reader_t* reader, md_allocator_i* alloc) {
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

static bool parse_rsp(md_vlx_rsp_t* rsp, md_buffered_reader_t* reader, md_allocator_i* alloc) {
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
			if (!parse_rsp_dipole_moments(rsp->electronic_transition_length, rsp->num_excited_states, reader, alloc)) {
				return false;
			}
			mask |= 2;
		} else if (str_eq(line, STR_LIT("Electric Transition Dipole Moments (dipole velocity, a.u.)"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_buffered_reader_skip_line(reader); // X Y Z
			md_array_resize(rsp->electronic_transition_velocity, rsp->num_excited_states, alloc);
			if (!parse_rsp_dipole_moments(rsp->electronic_transition_velocity, rsp->num_excited_states, reader, alloc)) {
				return false;
			}
			mask |= 4;
		} else if (str_eq(line, STR_LIT("Magnetic Transition Dipole Moments (a.u.)"))) {
			md_buffered_reader_skip_line(reader); // -----
			md_buffered_reader_skip_line(reader); // X Y Z
			md_array_resize(rsp->magnetic_transition, rsp->num_excited_states, alloc);
			if (!parse_rsp_dipole_moments(rsp->magnetic_transition, rsp->num_excited_states, reader, alloc)) {
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

bool md_vlx_data_parse_file(md_vlx_data_t* vlx, str_t filename, struct md_allocator_i* alloc) {
	md_file_o* file = md_file_open(filename, MD_FILE_READ | MD_FILE_BINARY);
	if (!file) {
		MD_LOG_ERROR("Failed to open file: '"STR_FMT"'", STR_ARG(filename));
		return false;
	}
	bool result = false;

	MEMSET(vlx, 0, sizeof(md_vlx_data_t));

	char buf[KILOBYTES(16)];
	md_buffered_reader_t reader = md_buffered_reader_from_file(buf, sizeof(buf), file);

	vlx->alloc = md_arena_allocator_create(alloc, MEGABYTES(1));

	str_t line;
	while (md_buffered_reader_extract_line(&line, &reader)) {
		str_t str = str_trim(line);
		if (str_eq(str, STR_LIT("Molecular Geometry (Angstroms)"))) {
			if (!parse_geom(&vlx->geom, &reader, vlx->alloc)) {
				MD_LOG_ERROR("Failed to parse geometry");
				goto done;
			}
		} else if (str_eq(str, STR_LIT("Molecular Basis (Atomic Basis)"))) {
			if (!parse_basis(&vlx->basis, &reader, vlx->alloc)) {
				MD_LOG_ERROR("Failed to parse basis");
				goto done;
			}
		} else if (str_eq(str, STR_LIT("Self Consistent Field Driver Setup"))) {
			if (!parse_scf(&vlx->scf, &reader, vlx->alloc)) {
				MD_LOG_ERROR("Failed to parse SCF section");
				goto done;
			}
		} else if (str_eq(str, STR_LIT("Lhd5 read dataset type and dataspaceinear Response EigenSolver Setup"))) {
			if (!parse_rsp(&vlx->rsp, &reader, vlx->alloc)) {
				MD_LOG_ERROR("Failed to parse SCF section");
				goto done;
			}
		}
	}

	// parse scf coefficients
	md_strb_t sb = md_strb_create(md_temp_allocator);
	size_t loc;
	if (str_rfind_char(&loc, filename, '.')) {
		md_strb_push_str(&sb, str_substr(filename, 0, loc));
		md_strb_push_str(&sb, STR_LIT(".scf.h5"));
	}

	hid_t  file_id, dataset_id, space_id; /* identifiers */
	herr_t status;

	/* Open an existing file. */
	file_id = H5Fopen(md_strb_to_cstr(sb), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id < 0) {
		return false;
	}

	/* Open an existing dataset. */
	dataset_id = H5Dopen(file_id, "/alpha_orbitals", H5P_DEFAULT);
	if (dataset_id < 0) {
		return false;
	}

	space_id = H5Dget_space(dataset_id);
	if (space_id < 0) {
		return false;
	}

	int ndim = H5Sget_simple_extent_ndims(space_id);
	if (ndim < 0) {
		return false;
	}
	
	hsize_t* dims = md_temp_push(sizeof(hsize_t) * ndim);
	ndim = H5Sget_simple_extent_dims(space_id, dims, 0);

	double* data = md_alloc(vlx->alloc, sizeof(double) * 229 * 229);
	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);



	status = H5Sclose(space_id);

	/* Close the dataset. */
	status = H5Dclose(dataset_id);

	/* Close the file. */
	status = H5Fclose(file_id);

	result = true;
done:
	md_file_close(file);
	return result;
}

void md_vlx_data_free(md_vlx_data_t* data) {
	ASSERT(data);
	if (data->alloc) {
		md_arena_allocator_destroy(data->alloc);
	}
	MEMSET(data, 0, sizeof(md_vlx_data_t));
}
