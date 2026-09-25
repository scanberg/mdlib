#include <md_edr.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <md_xdr.h>
#include <md_system.h>

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
				md_array_push(energies->energy[i].values, frame.ener[i].e, energies->alloc);
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

// ### ATTRIBUTES ###

// "Kinetic En." -> "kinetic_en". Lowercase ASCII letters and digits are kept, every run of anything
// else becomes one '_', and none is left at either end. Returns the length, 0 if nothing survived.
static size_t edr_slug(char* buf, size_t cap, str_t name) {
	size_t len = 0;
	bool pending_sep = false;
	for (size_t i = 0; i < name.len && len + 2 < cap; ++i) {
		char c = name.ptr[i];
		if (c >= 'A' && c <= 'Z') c = (char)(c - 'A' + 'a');
		const bool keep = (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9');
		if (!keep) {
			pending_sep = (len > 0);
			continue;
		}
		if (pending_sep) {
			buf[len++] = '_';
			pending_sep = false;
		}
		buf[len++] = c;
	}
	buf[len] = '\0';
	return len;
}

// One published attribute: a lone term, or a group of terms which are the components of one value.
typedef struct edr_output_t {
	str_t       label;        // the GROMACS name, or the prefix the members share
	uint32_t    tensor;       // 1: {R,3,3}, the members row major
	uint32_t    components;   // 1, or 3 for a vector whose members are x y z
	uint32_t    count;        // members
	uint32_t    member[9];
	const char* description;
} edr_output_t;

// Index of the term called <prefix>-<suffix>, or -1.
static int edr_find_term(const md_edr_energies_t* e, str_t prefix, const char* suffix) {
	const size_t slen = strlen(suffix);
	for (size_t i = 0; i < e->num_energies; ++i) {
		str_t name = e->energy[i].name;
		if (name.len == prefix.len + 1 + slen && str_begins_with(name, prefix) && name.ptr[prefix.len] == '-' &&
			MEMCMP(name.ptr + prefix.len + 1, suffix, slen) == 0) {
			return (int)i;
		}
	}
	return -1;
}

// Every suffix present under prefix, and all in one unit - or it is not a group and the terms stay
// single. Units are compared rather than assumed: a group whose members disagree would have to pick
// one of them, and whichever it picked would be wrong for the rest.
static bool edr_try_group(uint32_t out_members[], const md_edr_energies_t* e, str_t prefix, const char* const suffixes[], uint32_t count) {
	for (uint32_t k = 0; k < count; ++k) {
		const int idx = edr_find_term(e, prefix, suffixes[k]);
		if (idx < 0) {
			return false;
		}
		if (k > 0 && !md_unit_equal(e->energy[out_members[0]].unit, e->energy[idx].unit)) {
			return false;
		}
		out_members[k] = (uint32_t)idx;
	}
	return true;
}

static const char* const edr_tensor_suffix[9] = { "XX", "XY", "XZ", "YX", "YY", "YZ", "ZX", "ZY", "ZZ" };
static const char* const edr_vector_suffix[3] = { "X", "Y", "Z" };
static const char* const edr_diag_suffix[3]   = { "XX", "YY", "ZZ" };

// Decides what gets published. Each term ends up in exactly one output, in file order of its first
// member, so the table reads the way 'gmx energy' lists the file.
static size_t edr_plan(edr_output_t* out, const md_edr_energies_t* e, bool* used) {
	size_t count = 0;
	for (size_t i = 0; i < e->num_energies; ++i) {
		if (used[i]) continue;
		edr_output_t o = {0};

		size_t loc;
		str_t name = e->energy[i].name;
		if (str_rfind_char(&loc, name, '-') && loc > 0 && loc + 1 < name.len) {
			str_t prefix = str_substr(name, 0, loc);
			if (edr_try_group(o.member, e, prefix, edr_tensor_suffix, 9)) {
				o.tensor = 1; o.components = 1; o.count = 9;
				o.description = "3x3 tensor, row major, from <name>-XX .. <name>-ZZ";
			} else if (edr_try_group(o.member, e, prefix, edr_vector_suffix, 3)) {
				o.components = 3; o.count = 3;
				o.description = "x, y and z from <name>-X, <name>-Y and <name>-Z";
			} else if (edr_try_group(o.member, e, prefix, edr_diag_suffix, 3)) {
				o.components = 3; o.count = 3;
				o.description = "the diagonal of a tensor, from <name>-XX, <name>-YY and <name>-ZZ";
			}
			if (o.count) {
				// A member already taken means some other reading of the names claimed it first.
				bool free = true;
				for (uint32_t k = 0; k < o.count; ++k) free = free && !used[o.member[k]];
				if (free) {
					o.label = prefix;
				} else {
					o = (edr_output_t){0};
				}
			}
		}
		if (!o.count) {
			o.label = name;
			o.components = 1;
			o.count = 1;
			o.member[0] = (uint32_t)i;
		}
		for (uint32_t k = 0; k < o.count; ++k) used[o.member[k]] = true;
		out[count++] = o;
	}
	return count;
}

static str_t edr_join(char* buf, size_t cap, str_t a, const char* b, size_t b_len) {
	if (a.len + 1 + b_len + 1 > cap) {
		return (str_t){0};
	}
	MEMCPY(buf, a.ptr, a.len);
	buf[a.len] = '/';
	MEMCPY(buf + a.len + 1, b, b_len);
	buf[a.len + 1 + b_len] = '\0';
	return (str_t){buf, a.len + 1 + b_len};
}

bool md_edr_system_supplement(md_system_t* sys, const md_edr_energies_t* energies, str_t run) {
	ASSERT(sys);
	ASSERT(energies);

	md_attributes_t* attributes = &sys->attributes;
	if (!attributes->alloc) {
		attributes->alloc = sys->alloc;
	}
	if (str_empty(run)) {
		MD_LOG_ERROR("EDR: no run to publish the energies under");
		return false;
	}
	if (energies->num_frames == 0 || energies->num_energies == 0) {
		MD_LOG_ERROR("EDR: the file holds no energies");
		return false;
	}
	const size_t R = energies->num_frames;
	if (R > UINT32_MAX) {
		MD_LOG_ERROR("EDR: too many frames");
		return false;
	}
	for (size_t i = 1; i < R; ++i) {
		if (energies->frame_time[i] < energies->frame_time[i - 1]) {
			MD_LOG_ERROR("EDR: frame times decrease at frame %zu (%g ps after %g ps); concatenate restarted runs without overlap first",
				i, energies->frame_time[i], energies->frame_time[i - 1]);
			return false;
		}
	}

	char group_buf[512];
	char path_buf[512];
	const str_t group = edr_join(group_buf, sizeof(group_buf), run, "edr", 3);
	if (str_empty(group)) {
		MD_LOG_ERROR("EDR: run path '" STR_FMT "' is too long", STR_ARG(run));
		return false;
	}

	md_temp_scope_t temp = md_temp_begin();
	md_allocator_i* temp_alloc = md_temp_allocator(temp);
	bool result = false;

	// The file's own axis goes in first, against a scratch table: whether it covers the run is a
	// question about the values and has to be answered before the real table is touched.
	md_attributes_t probe = { .alloc = temp_alloc };
	const md_attribute_desc_t time_desc = {
		.path   = STR_LIT("time"),
		.format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = 1, .rank = 1, .shape = { (uint32_t)R } },
		.flags  = MD_ATTRIBUTE_FLAG_TEMPORAL,
		.unit   = md_unit_picosecond(),
		.label  = STR_LIT("Time"),
		.data   = energies->frame_time,
		.byte_size = R * sizeof(double),
	};
	const md_attribute_id_t probe_id = md_attributes_create(&probe, &time_desc);
	if (!probe_id) {
		goto done;
	}

	const md_attribute_t* run_axis = md_attributes_find(attributes, edr_join(path_buf, sizeof(path_buf), run, "time", 4));
	if (run_axis && md_attributes_axis(attributes, run_axis) == run_axis) {
		const md_attribute_t* edr_axis = md_attributes_get(&probe, probe_id);
		for (size_t f = 0; f < run_axis->format.shape[0]; ++f) {
			size_t row;
			if (!md_attribute_axis_map(&row, run_axis, f, edr_axis)) {
				double t = 0;
				const md_attribute_slice_t s = md_attribute_slice_1((uint32_t)f);
				md_attribute_extract_slice_f64(&t, 1, run_axis, &s, md_unit_picosecond());
				MD_LOG_ERROR("EDR: frame %zu of '" STR_FMT "' (%g ps) has no matching time in the energy file (%g - %g ps)",
					f, STR_ARG(run), t, energies->frame_time[0], energies->frame_time[R - 1]);
				goto done;
			}
		}
	}

	// Replacing, not merging: a second file's terms would otherwise sit beside the first's along an
	// axis that only describes one of them.
	md_attributes_remove_prefix(attributes, group);

	{
		md_attribute_desc_t desc = time_desc;
		desc.path = edr_join(path_buf, sizeof(path_buf), group, "time", 4);
		if (!md_attributes_create(attributes, &desc)) {
			goto fail;
		}
	}

	{
		bool* used = md_temp_alloc_array(temp, bool, energies->num_energies);
		MEMSET(used, 0, energies->num_energies * sizeof(bool));
		edr_output_t* plan = md_temp_alloc_array(temp, edr_output_t, energies->num_energies);
		const size_t num_out = edr_plan(plan, energies, used);

		// Slugs already handed out, to keep two terms that fold to the same name apart. 'time' and
		// 'source' are taken before the first term is looked at.
		md_array(str_t) taken = 0;
		md_array_push(taken, STR_LIT("time"), temp_alloc);
		md_array_push(taken, STR_LIT("source"), temp_alloc);

		for (size_t o = 0; o < num_out; ++o) {
			const edr_output_t* out = &plan[o];

			char slug[128];
			size_t slug_len = edr_slug(slug, sizeof(slug) - 8, out->label);
			if (slug_len == 0) {
				slug_len = (size_t)snprintf(slug, sizeof(slug), "term%zu", o);
			}
			const size_t base_len = slug_len;
			for (int n = 2;; ++n) {
				bool clash = false;
				for (size_t k = 0; k < md_array_size(taken); ++k) {
					if (str_eq(taken[k], (str_t){slug, slug_len})) { clash = true; break; }
				}
				if (!clash) break;
				slug_len = base_len + (size_t)snprintf(slug + base_len, sizeof(slug) - base_len, "_%d", n);
			}
			if (slug_len != base_len) {
				MD_LOG_INFO("EDR: '" STR_FMT "' published as '%s', its plain name was taken", STR_ARG(out->label), slug);
			}
			md_array_push(taken, str_copy((str_t){slug, slug_len}, temp_alloc), temp_alloc);

			const size_t per_row = out->count;
			double* values = md_temp_alloc_array(temp, double, R * per_row);
			for (size_t r = 0; r < R; ++r) {
				for (uint32_t k = 0; k < out->count; ++k) {
					values[r * per_row + k] = energies->energy[out->member[k]].values[r];
				}
			}

			md_attribute_desc_t desc = {
				.path   = edr_join(path_buf, sizeof(path_buf), group, slug, slug_len),
				.format = { .type = MD_ATTRIBUTE_TYPE_F64, .components = out->components, .rank = 1, .shape = { (uint32_t)R } },
				.flags  = MD_ATTRIBUTE_FLAG_TEMPORAL,
				.unit   = energies->energy[out->member[0]].unit,
				.label  = out->label,
				.description = out->description ? str_from_cstr(out->description) : (str_t){0},
				.data   = values,
				.byte_size = R * per_row * sizeof(double),
			};
			if (out->tensor) {
				desc.format.rank = 3;
				desc.format.shape[1] = 3;
				desc.format.shape[2] = 3;
			}
			if (str_empty(desc.path) || !md_attributes_create(attributes, &desc)) {
				goto fail;
			}
		}
	}

	result = true;
	goto done;

fail:
	// All or nothing: half an energy file is a table that answers some questions with the new file
	// and the rest with nothing, and nothing says which.
	MD_LOG_ERROR("EDR: failed to publish the energies under '" STR_FMT "'", STR_ARG(group));
	md_attributes_remove_prefix(attributes, group);

done:
	md_attributes_free(&probe);
	md_temp_end(temp);
	return result;
}

bool md_edr_system_supplement_from_file(md_system_t* sys, str_t filename, str_t run) {
	ASSERT(sys);

	md_edr_energies_t energies = {0};
	if (!md_edr_energies_parse_file(&energies, filename, md_get_heap_allocator())) {
		return false;
	}
	bool result = md_edr_system_supplement(sys, &energies, run);
	md_edr_energies_free(&energies);

	if (result) {
		char group_buf[512];
		char path_buf[512];
		const str_t group = edr_join(group_buf, sizeof(group_buf), run, "edr", 3);
		const md_attribute_desc_t desc = {
			.path   = edr_join(path_buf, sizeof(path_buf), group, "source", 6),
			.format = { .type = MD_ATTRIBUTE_TYPE_STR, .components = 1, .rank = 0 },
			.unit   = md_unit_none(),
			.label  = STR_LIT("Source"),
			.data   = &filename,
			.byte_size = sizeof(str_t),
		};
		md_attributes_create(&sys->attributes, &desc);
	}
	return result;
}
