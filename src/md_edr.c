#include <md_edr.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <md_xdr.h>

#include <string.h>
#include <stdio.h>


#define ENX_STRING_MAGIC -55555
#define ENX_HEADER_MAGIC -7777777
#define ENX_VERSION 5

enum
{
	enxOR,    /* Time and ensemble averaged data for orientation restraints */
	enxORI,   /* Instantaneous data for orientation restraints              */
	enxORT,   /* Order tensor(s) for orientation restraints                 */
	enxDISRE, /* Distance restraint blocks                                  */

	enxDHCOLL, /* Data about the free energy blocks in this frame.           */
	enxDHHIST, /* BAR histogram                                              */
	enxDH,     /* BAR raw delta H data                                       */

	enxAWH, /* AWH data */

	enxNR /* Total number of extra blocks in the current code,
		  * note that the enxio code can read files written by
		  * future code which contain more blocks.
		  */
};

typedef enum {
	MD_ENX_DATATYPE_INT32,
	MD_ENX_DATATYPE_FLOAT,
	MD_ENX_DATATYPE_DOUBLE,
	MD_ENX_DATATYPE_INT64,
	MD_ENX_DATATYPE_CHAR,
	MD_ENX_DATATYPE_STRING,

	MD_ENX_DATATYPE_COUNT
} md_enx_datatype_t;

typedef struct md_enxsubblock_t {
	int					nr;		// Number of items
	md_enx_datatype_t	type;	// Block data type

	union {
		float*			fval;
		double*			dval;
		int32_t*		ival;
		int64_t*		lval;
		unsigned char*	cval;
		char**			sval;
	} value;					// union of data arrays
} md_enxsubblock_t;

typedef struct md_enxblock_t {
	int					id;		// block id, from enx enums
	int					nsub;	// number of subblocks
	md_enxsubblock_t*	sub;	// subblock array
} md_enxblock_t;

typedef struct md_energy_t {
	double e;		// Current
	double eav;		// Running average
	double esum;	// Sum
} md_energy_t;

typedef struct md_enxframe_t {
	double			t;			// Timestamp of this frame
	int64_t			step;		// MD step
	int64_t			nsteps;		// Number of steps between frames
	double			dt;			// MD time step

	int				nsum;		// Number of terms for the sums in energyGroupPairTerms
	int				nre;		// Number of energies
	int				e_size;		// Size in bytes of energies
	md_energy_t*	ener;		// Energy array
	str_t*			e_names;
	str_t*			e_units;

	int				nblock;		// Number of blocks
	md_enxblock_t*	block;		// Block array
} md_enxframe_t;

typedef struct {
	md_xdr_t xdr;			// Over the whole file, which is read into memory up front
	void*    data;
	size_t   data_size;
	bool double_precision; // Are we reading double precision reals?

	struct {
		bool			old_file_open;		/* Is this an open old file? */
		bool			read_first_step;	/* Did we read the first step? */
		int				first_step;			/* First step in the energy file */
		int				step_prev;			/* Previous step */
		int				nsum_prev;			/* Previous step sum length */
		md_energy_t*	ener_prev;			/* Previous energy sums */
		md_allocator_i* alloc;
	} old;
} edr_fp_t;

static bool read_int(int* ptr, edr_fp_t* fp) {
	int32_t v;
	const bool ok = md_xdr_read_i32(&fp->xdr, &v);
	*ptr = v;
	return ok;
}

static bool read_real(double* ptr, int ndata, edr_fp_t* fp) {
	for (int i = 0; i < ndata; ++i) {
		if (fp->double_precision) {
			if (!md_xdr_read_f64(&fp->xdr, &ptr[i])) return false;
		} else {
			float flt;
			if (!md_xdr_read_f32(&fp->xdr, &flt)) return false;
			ptr[i] = flt;
		}
	}
	return true;
}

static void resize_data_subblock(md_enxsubblock_t* subblock, size_t new_size, md_allocator_i* alloc) {
	switch (subblock->type) {
	case MD_ENX_DATATYPE_INT32:  md_array_resize(subblock->value.ival, new_size, alloc); break;
	case MD_ENX_DATATYPE_FLOAT:  md_array_resize(subblock->value.fval, new_size, alloc); break;
	case MD_ENX_DATATYPE_DOUBLE: md_array_resize(subblock->value.dval, new_size, alloc); break;
	case MD_ENX_DATATYPE_INT64:  md_array_resize(subblock->value.lval, new_size, alloc); break;
	case MD_ENX_DATATYPE_CHAR:   md_array_resize(subblock->value.cval, new_size, alloc); break;
	case MD_ENX_DATATYPE_STRING: md_array_resize(subblock->value.sval, new_size, alloc); break;
	default:
		ASSERT(false);
	}
}

static void add_subblock_enxblock(md_enxblock_t* block, size_t new_size, md_allocator_i* alloc) {
	const size_t old_size = md_array_size(block->sub);
	md_array_resize(block->sub, new_size, alloc);
	if (new_size > old_size) {
		MEMSET(block->sub + old_size, 0, (new_size - old_size) * sizeof(md_enxsubblock_t));
	}
}

static void add_blocks_enxframe(md_enxframe_t* frame, size_t new_size, md_allocator_i* alloc) {
	const size_t old_size = md_array_size(frame->block);
	md_array_resize(frame->block, new_size, alloc);
	if (new_size > old_size) {
		MEMSET(frame->block + old_size, 0, (new_size - old_size) * sizeof(md_enxblock_t));
	}
}

static bool read_header(edr_fp_t* fp, md_enxframe_t* frame, int nre_test, bool* wrong_precision, md_allocator_i* alloc) {
	int magic = 0;
	int dum = 0;
	int ndisre = 0;
	int startb = 0;
	int file_version = 0;
	int nre = 0;

	double t = 0;
	double dt = 0;
	int64_t step = 0;
	int64_t nsteps = 0;
	int nsum = 0;
	int nblock = 0;

	if (wrong_precision) {
		*wrong_precision = false;
	}

	// Default real type
	const md_enx_datatype_t real_type = MD_ENX_DATATYPE_DOUBLE;

	double first_real_to_check = -2e10;
	if (!read_real(&first_real_to_check, 1, fp)) {
		MD_LOG_ERROR("Failed to read header");
		return false;
	}

	if (first_real_to_check > -1e10) {
		file_version = 1;
		if (!read_int(&dum, fp)) {
			MD_LOG_ERROR("Failed to read header step");
			return false;
		}
		t = first_real_to_check;
		step = dum;
	} else {
		if (!read_int(&magic, fp)) {
			MD_LOG_ERROR("Failed to read header magic");
			return false;
		}
		if (magic != ENX_HEADER_MAGIC) {
			MD_LOG_ERROR("Magic number mismatch in header, invalid edr file");
			return false;
		}
		if (!read_int(&file_version, fp)) {
			MD_LOG_ERROR("Falied to read header file version");
			return false;
		}
		if (file_version > ENX_VERSION) {
			MD_LOG_ERROR("Unsuported edr file version");
			return false;
		}
		if (!md_xdr_read_f64(&fp->xdr, &t)) {
			MD_LOG_ERROR("Falied to read time");
			return false;
		}
		if (!md_xdr_read_i64(&fp->xdr, &step)) {
			MD_LOG_ERROR("Falied to read step");
			return false;
		}
		if (!read_int(&nsum, fp)) {
			MD_LOG_ERROR("Falied to read nsum");
			return false;
		}
		if (file_version >= 3) {
			if (!md_xdr_read_i64(&fp->xdr, &nsteps)) {
				MD_LOG_ERROR("Falied to read nsteps");
				return false;
			}
		} else {
			nsteps = MAX(1, nsum);
		}
		if (file_version >= 5) {
			if (!md_xdr_read_f64(&fp->xdr, &dt)) {
				MD_LOG_ERROR("Falied to read dt");
				return false;
			}
		} else {
			dt = 0;
		}
	}
	if (!read_int(&nre, fp)) {
		MD_LOG_ERROR("Falied to read nre");
		return false;
	}
	if (file_version < 4) {
		if (!read_int(&ndisre, fp)) {
			MD_LOG_ERROR("Falied to read ndisre");
			return false;
		}
	}
	else {
		// Reserved for possible future use
		if (!read_int(&dum, fp)) {
			MD_LOG_ERROR("Falied to read data");
			return false;
		}
	}

	if (!read_int(&nblock, fp)) {
		MD_LOG_ERROR("Falied to read nblock");
		return false;
	}
	if (nblock < 0) {
		MD_LOG_ERROR("Negative number of blocks, edr is corrupt");
		return false;
	}

	if (ndisre != 0) {
		if (file_version >= 4) {
			MD_LOG_ERROR("Distance restraint blocks in old style in new style file");
			return false;
		}
		nblock += 1;
	}

	// Frames could have nre=0, so we can not rely only on the fr->nre check
	if (nre_test >= 0 && ((nre > 0 && nre != nre_test) || nre < 0 || ndisre < 0 || nblock < 0)) {
		if (wrong_precision) {
			*wrong_precision = true;
		}
		return false;
	}

	// we now know what these should be, or we've already bailed out because of wrong precision
	if (file_version == 1 && (t < 0 || t > 1e20 || step < 0)) {
		MD_LOG_ERROR("edr file with negative step number or unreasonable time (and without version number).");
		return false;
	}

	if (frame) {
		frame->t = t;
		frame->dt = dt;
		frame->step = step;
		frame->nsteps = nsteps;
		frame->nsum = nsum;
		frame->nre = nre;
		frame->nblock = nblock;
		add_blocks_enxframe(frame, frame->nblock, alloc);
	}

	startb = 0;
	if (ndisre > 0) {
		// sub[0] is the instantaneous data, sub[1] is time averaged
		if (frame) {
			add_subblock_enxblock(&frame->block[0], 2, alloc);
			frame->block[0].id = enxDISRE;

			frame->block[0].sub[0].type = real_type;
			frame->block[0].sub[0].nr = ndisre;

			frame->block[0].sub[1].type = real_type;
			frame->block[0].sub[1].nr = ndisre;
		}

		startb += 1;
	}

	for (int b = startb; b < nblock; ++b) {
		if (file_version < 4) {
			// blocks in old version files always have 1 subblock that consists of reals.
			int nrint;
			if (!read_int(&nrint, fp)) {
				return false;
			}
			if (frame) {
				add_subblock_enxblock(&frame->block[b], 1, alloc);
				frame->block[b].id = b - startb;
				frame->block[b].nsub = nrint;
				frame->block[b].sub[0].type = real_type;
				resize_data_subblock(&frame->block[b].sub[0], nrint, alloc);
			}
		} else {
			int id, nsub;
			if (!read_int(&id, fp)) {
				return false;
			}
			if (!read_int(&nsub, fp)) {
				return false;
			}

			if (frame) {
				frame->block[b].id = id;
				frame->block[b].nsub = nsub;
				add_subblock_enxblock(&frame->block[b], nsub, alloc);
			}

			for (int i = 0; i < nsub; ++i) {
				int typenr, nr;
				if (!read_int(&typenr, fp)) {
					return false;
				}
				if (!read_int(&nr, fp)) {
					return false;
				}
				if (typenr < 0 || typenr >= MD_ENX_DATATYPE_COUNT) {
					MD_LOG_ERROR("EDR: Invalid typenr inside subblock");
					return false;
				}
				/*
				// Can this actually happen?
				if (nr <= 0) {
					MD_LOG_ERROR("EDR: Invalid nr inside subblock");
					return false;
				}
				*/

				if (frame) {
					frame->block[b].sub[i].type = typenr;
					frame->block[b].sub[i].nr = nr;
					resize_data_subblock(&frame->block[b].sub[i], nr, alloc);
				}
			}
		}
	}

	int e_size;
	if (!read_int(&e_size, fp)) {
		return false;
	}

	if (frame) {
		frame->e_size = e_size;
	}

	// Reserved for future use
	if (!read_int(&dum, fp)) {
		return false;
	}
	if (!read_int(&dum, fp)) {
		return false;
	}

	if (file_version == 1 && nre_test < 0) {
		if (!fp->old.read_first_step)
		{
			fp->old.read_first_step = true;
			fp->old.first_step		= (int)step;
			fp->old.step_prev		= (int)step;
			fp->old.nsum_prev		= 0;
		}

		if (frame) {
			frame->nsum   = (int)frame->step - fp->old.first_step + 1;
			frame->nsteps = (int)frame->step - fp->old.step_prev;
			frame->dt     = 0;
		}
	}

	return true;
}

static inline double square(double x) { return x * x; }

static void convert_full_sums(edr_fp_t* fp, md_enxframe_t* fr) {
	int    nstep_all;
	int    ne, ns, i;
	double esum_all, eav_all;

	if (fr->nsum > 0) {
		ne = 0;
		ns = 0;
		for (i = 0; i < fr->nre; i++) {
			if (fr->ener[i].e != 0) {
				ne++;
			}
			if (fr->ener[i].esum != 0) {
				ns++;
			}
		}
		if (ne > 0 && ns == 0) {
			/* We do not have all energy sums */
			fr->nsum = 0;
		}
	}

	/* Convert old full simulation sums to sums between energy frames */
	nstep_all = (int)fr->step - fp->old.first_step + 1;
	if (fr->nsum > 1 && fr->nsum == nstep_all && fp->old.nsum_prev > 0) {
		/* Set the new sum length: the frame step difference */
		fr->nsum = (int)fr->step - fp->old.step_prev;
		for (i = 0; i < fr->nre; i++) {
			esum_all         = fr->ener[i].esum;
			eav_all          = fr->ener[i].eav;
			fr->ener[i].esum = esum_all - fp->old.ener_prev[i].esum;
			fr->ener[i].eav =
				eav_all - fp->old.ener_prev[i].eav
				- square(fp->old.ener_prev[i].esum / (nstep_all - fr->nsum) - esum_all / nstep_all)
				* (nstep_all - fr->nsum) * nstep_all / (double)(fr->nsum);
			fp->old.ener_prev[i].esum = esum_all;
			fp->old.ener_prev[i].eav  = eav_all;
		}
		fp->old.nsum_prev = nstep_all;
	}
	else if (fr->nsum > 0) {
		if (fr->nsum != nstep_all) {
			MD_LOG_ERROR("Something is wrong with the energy sums, will not use exact averages");
			fp->old.nsum_prev = 0;
		}
		else {
			fp->old.nsum_prev = nstep_all;
		}
		/* Copy all sums to ener_prev */
		for (i = 0; i < fr->nre; i++) {
			fp->old.ener_prev[i].esum = fr->ener[i].esum;
			fp->old.ener_prev[i].eav  = fr->ener[i].eav;
		}
	}

	fp->old.step_prev = (int)fr->step;
}

static bool read_frame(edr_fp_t* fp, md_enxframe_t* frame, int file_version, md_allocator_i* alloc) {
	int64_t	frame_num  = 0;
	double  frame_time = 0;
	char	buf[1024];

	if (!read_header(fp, frame, -1, NULL, alloc)) {
		MD_LOG_ERROR("Last energy frame read %d time %8.3f", (int)frame_num - 1, frame_time);
		return false;
	}

	if (frame->step < 0) {
		MD_LOG_ERROR("Something went wrong when reading header");
		return false;
	}

	frame_num += 1;
	frame_time = frame->t;

	bool sane = (frame->nre > 0);
	for (int i = 0; i < frame->nblock; ++i) {
		sane = sane || (frame->block[i].nsub > 0);
	}

	if (!(frame->step >= 0 && sane)) {
		MD_LOG_ERROR("EDR: Something went wrong when reading frame header");
        return false;
    }

	if (frame->nre > (int)md_array_size(frame->ener)) {
		const size_t new_size = (size_t)frame->nre;
		const size_t old_size = md_array_size(frame->ener);
		md_array_resize(frame->ener, (size_t)frame->nre, alloc);
		MEMSET(frame->ener + old_size, 0, (new_size - old_size) * sizeof(md_energy_t));
	}

	for (int i = 0; i < frame->nre; ++i) {
		if (!read_real(&frame->ener[i].e, 1, fp)) {
			return false;
		}
		if (file_version == 1 || frame->nsum > 0) {
			if (!read_real(&frame->ener[i].eav, 1, fp)) {
				return false;
			}
			if (!read_real(&frame->ener[i].esum, 1, fp)) {
				return false;
			}
			if (file_version == 1) {
				double dummy;
				// Old, unused real
				if (!read_real(&dummy, 1, fp)) {
					return false;
				}
			}
		}
	}

	if (fp->old.old_file_open) {
		convert_full_sums(fp, frame);
	}

	bool ok = true;
	for (int b = 0; b < frame->nblock; ++b) {
		/* now read the subblocks. */
        int nsub = frame->block[b].nsub;

		for (int i = 0; i < nsub; i++) {
			md_enxsubblock_t* sub = &(frame->block[b].sub[i]);

			/* read data */
			switch (sub->type)
			{
			case MD_ENX_DATATYPE_FLOAT:
				ok = ok && md_xdr_read_f32_array(&fp->xdr, sub->value.fval, (size_t)sub->nr);
				break;
			case MD_ENX_DATATYPE_DOUBLE:
				ok = ok && md_xdr_read_f64_array(&fp->xdr, sub->value.dval, (size_t)sub->nr);
				break;
			case MD_ENX_DATATYPE_INT32:
				ok = ok && md_xdr_read_i32_array(&fp->xdr, sub->value.ival, (size_t)sub->nr);
				break;
			case MD_ENX_DATATYPE_INT64:
				for (int j = 0; j < sub->nr && ok; ++j) {
					ok = md_xdr_read_i64(&fp->xdr, &sub->value.lval[j]);
				}
				break;
			case MD_ENX_DATATYPE_CHAR:
				// XDR u_char: one per 4 bytes
				for (int j = 0; j < sub->nr && ok; ++j) {
					uint32_t c;
					ok = md_xdr_read_u32(&fp->xdr, &c);
					sub->value.cval[j] = (unsigned char)c;
				}
				break;
			case MD_ENX_DATATYPE_STRING:
				for (int j = 0; j < (int)md_array_size(sub->value.sval) && ok; ++j) {
					str_t str;
					ok = md_xdr_read_string(&fp->xdr, &str, sizeof(buf) - 1);
					if (ok) {
						sub->value.sval[j] = md_alloc(alloc, str.len + 1);
						MEMCPY(sub->value.sval[j], str.ptr, str.len);
						sub->value.sval[j][str.len] = '\0';
					}
				}
				break;
			default:
				md_log(MD_LOG_TYPE_DEBUG, "Reading unknown block data type: this file is corrupted or from a future version");
				ok = false;
			}
			if (!ok) {
				goto done;
			}
		}
	}

done:
	if (!ok) {
		MD_LOG_ERROR("\nLast energy frame read %d", (int)frame_num - 1);
		MD_LOG_ERROR("\nWARNING: Incomplete energy frame: nr %d time %8.3f\n", (int)frame_num, frame->t);
		return false;
	}

	return true;
}

static bool read_strings(edr_fp_t* fp, md_enxframe_t* frame, int* file_version, md_allocator_i* alloc) {
	// Read Strings
	bool result = false;

	int magic;
	if (!read_int(&magic, fp)) {
		MD_LOG_ERROR("Failed to read magic from edr file");
		goto done;
	}

	int nre;
	if (magic > 0) {
		// Assume old format
		*file_version = 1;
		nre = magic;
		fp->old.old_file_open = true;
		fp->old.read_first_step = false;
	} else {
		if (magic != ENX_STRING_MAGIC) {
			MD_LOG_ERROR("Failed to open edr file: magic number mismatch when reading edr file strings");
			goto done;
		}
		*file_version = ENX_VERSION;
		read_int(file_version, fp);
		if (*file_version > ENX_VERSION) {
			MD_LOG_ERROR("Failed to open edr file: unsupported file version in edr file");
			goto done;
		}
		read_int(&nre, fp);
	}

	if (frame) {
		frame->nre = nre;
	}

	for (int i = 0; i < nre; ++i) {
		// A string ends at its first zero, whatever length the file gives it
		str_t name, unit;
		if (!md_xdr_read_string(&fp->xdr, &name, 1023)) {
			MD_LOG_ERROR("Failed to read expected number of strings within edr file");
			goto done;
		}
		name = str_from_cstrn(name.ptr, name.len);

		if (*file_version > 1) {
			if (!md_xdr_read_string(&fp->xdr, &unit, 1023)) {
				MD_LOG_ERROR("Failed to read expected number of strings within edr file");
				goto done;
			}
			unit = str_from_cstrn(unit.ptr, unit.len);
		} else {
			unit = STR_LIT("kJ/mol");
		}

		if (frame) {
			md_array_push(frame->e_names, str_copy(name, alloc), alloc);
			md_array_push(frame->e_units, str_copy(unit, alloc), alloc);
		}
	}

	result = true;
done:
	return result;
}

static bool edr_file_open(edr_fp_t* fp, str_t filename) {
	md_temp_scope_t temp_scope = md_temp_begin();
	md_allocator_i* arena = md_temp_allocator(temp_scope);
	str_t path = str_copy(filename, arena);
	bool result = false;

	{
		md_file_t file = {0};
		if (!md_file_open(&file, filename, MD_FILE_READ)) {
			MD_LOG_ERROR("Failed to open file '%.*s'", (int)path.len, path.ptr);
			md_temp_end(temp_scope);
			return false;
		}
		const int64_t size = md_file_size(file);
		fp->data = size > 0 ? md_alloc(md_get_heap_allocator(), (size_t)size) : NULL;
		fp->data_size = fp->data ? (size_t)size : 0;
		const bool read_ok = fp->data && md_file_read(file, fp->data, fp->data_size) == fp->data_size;
		md_file_close(&file);
		if (!read_ok) {
			MD_LOG_ERROR("Failed to read file '%.*s'", (int)path.len, path.ptr);
			md_temp_end(temp_scope);
			return false;
		}
		fp->xdr = md_xdr_init(fp->data, fp->data_size);
	}

	bool wrong_precision = false;
	md_enxframe_t frame = {0};
	int file_version;
	
	if (read_strings(fp, &frame, &file_version, arena) &&
		read_header(fp, &frame, frame.nre, &wrong_precision, arena) && !wrong_precision &&
		frame.e_size && (frame.nre * 4 * (int64_t)(sizeof(float)) == frame.e_size))
	{
		md_logf(MD_LOG_TYPE_INFO, "Opened '%.*s' as single precision energy file", (int)path.len, path.ptr);
	} else {
		fp->xdr = md_xdr_init(fp->data, fp->data_size);
		frame = (md_enxframe_t){0};
		fp->double_precision = true;

		if (read_strings(fp, &frame, &file_version, arena) &&
			read_header(fp, &frame, frame.nre, NULL, arena) &&
			frame.e_size && (frame.nre * 4 * (int64_t)(sizeof(double)) == frame.e_size))
		{
			md_logf(MD_LOG_TYPE_INFO, "Opened '%.*s' as double precision energy file", (int)path.len, path.ptr);
		} else {
			MD_LOG_ERROR("Failed to open edr file: format was not recognized");
			goto done;
		}
	}

	if (fp->old.old_file_open) {
		md_array_resize(fp->old.ener_prev, (size_t)frame.nre, md_get_heap_allocator());
	}

	fp->xdr = md_xdr_init(fp->data, fp->data_size);
	result = true;
done:
    md_temp_end(temp_scope);
	return result;
}

static void edr_file_close(edr_fp_t* fp) {
	if (fp->data) {
		md_free(md_get_heap_allocator(), fp->data, fp->data_size);
		fp->data = NULL;
	}
	if (fp->old.step_prev) {
		md_array_free(fp->old.ener_prev, md_get_heap_allocator());
	}
}

static md_unit_t unit_from_str(str_t str) {
	// Yields 'none' if the string is empty or could not be parsed
	md_unit_t unit;
	md_unit_parse(&unit, str);
	return unit;
}

bool md_edr_energies_parse_file(md_edr_energies_t* energies, str_t filename, struct md_allocator_i* alloc) {
	bool result = false;
	edr_fp_t fp = {0};
	if (!edr_file_open(&fp, filename)) {
		return false;
	}

	md_temp_scope_t temp_scope = md_temp_begin_avoid(alloc);
	md_allocator_i* temp = md_temp_allocator(temp_scope);
	md_enxframe_t frame = {0};

	int file_version = 0;
	if (!read_strings(&fp, &frame, &file_version, temp)) {
		goto done;
	}

	if (energies->energy != NULL && energies->alloc != 0) {
		md_log(MD_LOG_TYPE_DEBUG, "Reading energies into non-zero energy structure, potential memory leak here");
	}
	*energies = (md_edr_energies_t){0};
	energies->alloc = md_arena_allocator_create(alloc, MEGABYTES(1));

	energies->num_frames = 0;
	energies->frame_time = NULL;

	energies->num_energies = frame.nre;
	md_array_resize(energies->energy, (size_t)frame.nre, energies->alloc);
	
	for (int i = 0; i < frame.nre; ++i) {
        energies->energy[i].name	 = str_copy(frame.e_names[i], energies->alloc);
		energies->energy[i].unit_str = str_copy(frame.e_units[i], energies->alloc);
        energies->energy[i].unit	 = unit_from_str(frame.e_units[i]);
		energies->energy[i].values   = NULL;
	}

	while (true) {
		if (!read_frame(&fp, &frame, file_version, temp)) {
			MD_LOG_ERROR("Failed to read complete edr file!");
			md_edr_energies_free(energies);
			goto done;
		}

		if (frame.nre > 0) {
			// Only export frames which contain energies
			md_array_push(energies->frame_time, frame.t, energies->alloc);

			for (int i = 0; i < frame.nre; ++i) {
				md_array_push(energies->energy[i].values, (float)frame.ener[i].e, energies->alloc);
			}
			energies->num_frames += 1;
		}

		if (md_xdr_remaining(&fp.xdr) == 0) {
			break;
		}
	}
	
	result = true;
done:
    md_temp_end(temp_scope);
	edr_file_close(&fp);
	return result;
}

void md_edr_energies_free(md_edr_energies_t* energies) {
	if (energies->alloc) md_arena_allocator_destroy(energies->alloc);
	*energies = (md_edr_energies_t){0};
}
