#include <md_trexio.h>

#include <md_system.h>
#include <md_types.h>
#include <md_gto.h>
#include <md_qm.h>

#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_common.h>
#include <core/md_log.h>
#include <core/md_str.h>
#include <core/md_unit.h>
#include <core/md_vec_math.h>

#include <hdf5.h>

#include <math.h>
#include <string.h>

#define TREXIO_BOHR_TO_ANGSTROM 0.529177210903

// ---------------------------------------------------------------------------
// HDF5 helpers
// ---------------------------------------------------------------------------
// TREXIO's HDF5 back end is regular enough that these six helpers cover the whole file: a scalar is
// an attribute on a section group, an array is a dataset under it, and everything is either an
// int64, a double or a string.

typedef struct h5_error_scope_t {
    H5E_auto2_t func;
    void*       client_data;
} h5_error_scope_t;

// HDF5 prints its own stack to stderr on every failed probe, and probing for absent optional
// datasets is most of what this reader does.
static h5_error_scope_t h5_error_scope_begin(void) {
    h5_error_scope_t scope = {0};
    H5Eget_auto2(H5E_DEFAULT, &scope.func, &scope.client_data);
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
    return scope;
}

static void h5_error_scope_end(h5_error_scope_t scope) {
    H5Eset_auto2(H5E_DEFAULT, scope.func, scope.client_data);
}

static bool trexio_has(hid_t loc, const char* name) {
    return H5Lexists(loc, name, H5P_DEFAULT) > 0;
}

// A TREXIO scalar. Absent is a normal answer - most sections are optional - so it is reported
// through the return value and never logged.
static bool trexio_read_i64_attr(int64_t* out, hid_t loc, const char* name) {
    if (H5Aexists(loc, name) <= 0) {
        return false;
    }
    hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
    if (attr == H5I_INVALID_HID) {
        return false;
    }
    const bool ok = H5Aread(attr, H5T_NATIVE_INT64, out) >= 0;
    H5Aclose(attr);
    return ok;
}

static bool trexio_read_f64_attr(double* out, hid_t loc, const char* name) {
    if (H5Aexists(loc, name) <= 0) {
        return false;
    }
    hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
    if (attr == H5I_INVALID_HID) {
        return false;
    }
    const bool ok = H5Aread(attr, H5T_NATIVE_DOUBLE, out) >= 0;
    H5Aclose(attr);
    return ok;
}

// A string attribute, fixed or variable length. Returns a copy in 'alloc', empty when absent.
static str_t trexio_read_str_attr(hid_t loc, const char* name, md_allocator_i* alloc) {
    str_t result = {0};
    if (H5Aexists(loc, name) <= 0) {
        return result;
    }
    hid_t attr = H5Aopen(loc, name, H5P_DEFAULT);
    if (attr == H5I_INVALID_HID) {
        return result;
    }
    hid_t type = H5Aget_type(attr);
    if (type == H5I_INVALID_HID || H5Tget_class(type) != H5T_STRING) {
        if (type != H5I_INVALID_HID) H5Tclose(type);
        H5Aclose(attr);
        return result;
    }

    if (H5Tis_variable_str(type)) {
        char* var = NULL;
        if (H5Aread(attr, type, &var) >= 0 && var) {
            result = str_copy_cstr(var, alloc);
            H5free_memory(var);
        }
    } else {
        const size_t size = H5Tget_size(type);
        if (size > 0) {
            char* buf = md_alloc(alloc, size + 1);
            if (buf) {
                MEMSET(buf, 0, size + 1);
                if (H5Aread(attr, type, buf) >= 0) {
                    result = str_copy_cstrn(buf, strnlen(buf, size), alloc);
                }
            }
        }
    }

    H5Tclose(type);
    H5Aclose(attr);
    return result;
}

// Total element count of a dataset, 0 when it is absent or unreadable.
static size_t trexio_dataset_count(hid_t loc, const char* name) {
    if (!trexio_has(loc, name)) {
        return 0;
    }
    hid_t dset = H5Dopen(loc, name, H5P_DEFAULT);
    if (dset == H5I_INVALID_HID) {
        return 0;
    }
    size_t count = 0;
    hid_t space = H5Dget_space(dset);
    if (space != H5I_INVALID_HID) {
        const int rank = H5Sget_simple_extent_ndims(space);
        if (rank > 0 && rank <= 4) {
            hsize_t dims[4] = {0};
            if (H5Sget_simple_extent_dims(space, dims, NULL) == rank) {
                count = 1;
                for (int i = 0; i < rank; ++i) count *= (size_t)dims[i];
            }
        }
        H5Sclose(space);
    }
    H5Dclose(dset);
    return count;
}

// Reads a whole numeric dataset into 'dst', which must hold trexio_dataset_count elements.
static bool trexio_read_dataset(void* dst, hid_t loc, const char* name, hid_t mem_type) {
    if (!trexio_has(loc, name)) {
        return false;
    }
    hid_t dset = H5Dopen(loc, name, H5P_DEFAULT);
    if (dset == H5I_INVALID_HID) {
        return false;
    }
    const bool ok = H5Dread(dset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, dst) >= 0;
    H5Dclose(dset);
    if (!ok) {
        MD_LOG_ERROR("TREXIO: failed to read dataset '%s'", name);
    }
    return ok;
}

// A dataset of strings, one per element. Returns the count written, 0 on absence or failure.
static size_t trexio_read_str_dataset(str_t dst[], size_t cap, hid_t loc, const char* name, md_allocator_i* alloc) {
    if (!trexio_has(loc, name)) {
        return 0;
    }
    hid_t dset = H5Dopen(loc, name, H5P_DEFAULT);
    if (dset == H5I_INVALID_HID) {
        return 0;
    }

    size_t written = 0;
    hid_t file_type = H5Dget_type(dset);
    hid_t mem_type  = H5I_INVALID_HID;

    if (file_type != H5I_INVALID_HID && H5Tget_class(file_type) == H5T_STRING) {
        if (H5Tis_variable_str(file_type)) {
            mem_type = H5Tcopy(H5T_C_S1);
            H5Tset_size(mem_type, H5T_VARIABLE);
            char** raw = (char**)md_alloc(alloc, sizeof(char*) * cap);
            if (raw) {
                MEMSET(raw, 0, sizeof(char*) * cap);
                if (H5Dread(dset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, raw) >= 0) {
                    for (size_t i = 0; i < cap; ++i) {
                        dst[i] = raw[i] ? str_copy_cstr(raw[i], alloc) : (str_t){0};
                    }
                    written = cap;
                    H5Dvlen_reclaim(mem_type, H5Dget_space(dset), H5P_DEFAULT, raw);
                }
            }
        } else {
            const size_t size = H5Tget_size(file_type);
            char* raw = (size > 0) ? (char*)md_alloc(alloc, size * cap + 1) : NULL;
            if (raw) {
                MEMSET(raw, 0, size * cap + 1);
                if (H5Dread(dset, file_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, raw) >= 0) {
                    for (size_t i = 0; i < cap; ++i) {
                        dst[i] = str_copy_cstrn(raw + i * size, strnlen(raw + i * size, size), alloc);
                    }
                    written = cap;
                }
            }
        }
    }

    if (mem_type != H5I_INVALID_HID) H5Tclose(mem_type);
    if (file_type != H5I_INVALID_HID) H5Tclose(file_type);
    H5Dclose(dset);
    return written;
}

// ---------------------------------------------------------------------------
// The reader's own representation
// ---------------------------------------------------------------------------

typedef struct trexio_t {
    md_allocator_i* alloc;

    str_t code;
    str_t package_version;
    str_t basis_type;
    str_t mo_type;

    int64_t electron_num;
    int64_t electron_up_num;
    int64_t electron_dn_num;
    bool    has_electron_count;

    double nuclear_repulsion;
    bool   has_nuclear_repulsion;

    size_t             num_atoms;
    md_array(uint8_t)  atomic_number;
    md_array(dvec3_t)  coord;      // Angstrom

    bool               cartesian;  // ao_cartesian != 0
    size_t             num_shells;
    size_t             num_prims;
    md_array(uint32_t) shell_atom;
    md_array(uint32_t) shell_l;
    md_array(uint32_t) shell_prim_offset;
    md_array(uint32_t) shell_prim_count;
    md_array(double)   alpha;
    md_array(double)   coeff;

    size_t             num_ao;
    md_array(double)   ao_scale;   // file AO -> md_gto unit normalised AO, per AO

    size_t             num_mo;
    md_array(double)   mo_energy;
    md_array(double)   mo_occupation;
    md_array(int64_t)  mo_spin;
    md_array(str_t)    mo_symmetry;
    md_array(double)   mo_coefficient;  // [num_mo][num_ao], in the FILE's AO order

} trexio_t;

// TREXIO's spherical order is m = 0, +1, -1, +2, -2, ... - the same as Molden's. md_gto's is m
// ascending, -l .. +l. This is md_gto's index for the file's function 'file_idx'.
static uint32_t trexio_sph_to_md_index(uint32_t l, uint32_t file_idx) {
    const int m = (file_idx == 0) ? 0 : (((file_idx & 1) ? 1 : -1) * (int)((file_idx + 1) / 2));
    return (uint32_t)(m + (int)l);
}

// ---------------------------------------------------------------------------
// Reading
// ---------------------------------------------------------------------------

static bool trexio_read_nucleus(trexio_t* trexio, hid_t file) {
    hid_t group = H5Gopen(file, "nucleus", H5P_DEFAULT);
    if (group == H5I_INVALID_HID) {
        MD_LOG_ERROR("TREXIO: file has no nucleus group");
        return false;
    }

    bool result = false;
    int64_t num = 0;
    if (!trexio_read_i64_attr(&num, group, "nucleus_num") || num <= 0) {
        MD_LOG_ERROR("TREXIO: nucleus_num is missing or not positive");
        goto done;
    }

    trexio->num_atoms = (size_t)num;

    double* charge = (double*)md_alloc(trexio->alloc, sizeof(double) * trexio->num_atoms);
    double* coord  = (double*)md_alloc(trexio->alloc, sizeof(double) * trexio->num_atoms * 3);
    str_t*  label  = (str_t*) md_alloc(trexio->alloc, sizeof(str_t)  * trexio->num_atoms);

    if (charge && coord && label
        && trexio_read_dataset(charge, group, "nucleus_charge", H5T_NATIVE_DOUBLE)
        && trexio_read_dataset(coord,  group, "nucleus_coord",  H5T_NATIVE_DOUBLE)) {

        MEMSET(label, 0, sizeof(str_t) * trexio->num_atoms);
        trexio_read_str_dataset(label, trexio->num_atoms, group, "nucleus_label", trexio->alloc);

        for (size_t i = 0; i < trexio->num_atoms; ++i) {
            // The label is authoritative where there is one: an ECP calculation states the charge
            // the electrons see, which is not the atomic number.
            md_atomic_number_t z = str_empty(label[i]) ? 0 : md_atomic_number_from_symbol(label[i], true);
            if (z == 0) {
                const double rounded = round(charge[i]);
                z = (md_atomic_number_t)CLAMP((int)rounded, 0, 118);
            }
            const dvec3_t xyz = {
                coord[i * 3 + 0] * TREXIO_BOHR_TO_ANGSTROM,
                coord[i * 3 + 1] * TREXIO_BOHR_TO_ANGSTROM,
                coord[i * 3 + 2] * TREXIO_BOHR_TO_ANGSTROM,
            };
            md_array_push(trexio->atomic_number, (uint8_t)z, trexio->alloc);
            md_array_push(trexio->coord, xyz, trexio->alloc);
        }
        result = true;
    }

    trexio->has_nuclear_repulsion = trexio_read_f64_attr(&trexio->nuclear_repulsion, group, "nucleus_repulsion");

done:
    H5Gclose(group);
    return result;
}

static void trexio_read_electron(trexio_t* trexio, hid_t file) {
    if (!trexio_has(file, "electron")) {
        return;
    }
    hid_t group = H5Gopen(file, "electron", H5P_DEFAULT);
    if (group == H5I_INVALID_HID) {
        return;
    }
    trexio->has_electron_count = trexio_read_i64_attr(&trexio->electron_num, group, "electron_num");
    trexio_read_i64_attr(&trexio->electron_up_num, group, "electron_up_num");
    trexio_read_i64_attr(&trexio->electron_dn_num, group, "electron_dn_num");
    H5Gclose(group);
}

static void trexio_read_metadata(trexio_t* trexio, hid_t file) {
    if (!trexio_has(file, "metadata")) {
        return;
    }
    hid_t group = H5Gopen(file, "metadata", H5P_DEFAULT);
    if (group == H5I_INVALID_HID) {
        return;
    }
    trexio->code            = trexio_read_str_attr(group, "metadata_code", trexio->alloc);
    trexio->package_version = trexio_read_str_attr(group, "metadata_package_version", trexio->alloc);
    H5Gclose(group);
}

static bool trexio_read_basis(trexio_t* trexio, hid_t file) {
    if (!trexio_has(file, "basis")) {
        return true;  // A geometry-only file is legal; it simply publishes no basis.
    }
    hid_t group = H5Gopen(file, "basis", H5P_DEFAULT);
    if (group == H5I_INVALID_HID) {
        return true;
    }

    bool result = true;
    int64_t shell_num = 0;
    int64_t prim_num  = 0;

    trexio->basis_type = trexio_read_str_attr(group, "basis_type", trexio->alloc);
    if (!str_empty(trexio->basis_type) && !str_eq_cstr_ignore_case(trexio->basis_type, "gaussian")) {
        MD_LOG_ERROR("TREXIO: basis_type is '" STR_FMT "'; this library evaluates Gaussian basis sets only", STR_ARG(trexio->basis_type));
        result = false;
        goto done;
    }

    if (!trexio_read_i64_attr(&shell_num, group, "basis_shell_num") ||
        !trexio_read_i64_attr(&prim_num,  group, "basis_prim_num")  ||
        shell_num <= 0 || prim_num <= 0) {
        goto done;  // No basis, which is not an error.
    }

    trexio->num_shells = (size_t)shell_num;
    trexio->num_prims  = (size_t)prim_num;

    int64_t* nucleus_index = (int64_t*)md_alloc(trexio->alloc, sizeof(int64_t) * trexio->num_shells);
    int64_t* ang_mom       = (int64_t*)md_alloc(trexio->alloc, sizeof(int64_t) * trexio->num_shells);
    int64_t* r_power       = (int64_t*)md_alloc(trexio->alloc, sizeof(int64_t) * trexio->num_shells);
    double*  shell_factor  = (double*) md_alloc(trexio->alloc, sizeof(double)  * trexio->num_shells);
    int64_t* shell_index   = (int64_t*)md_alloc(trexio->alloc, sizeof(int64_t) * trexio->num_prims);
    double*  exponent      = (double*) md_alloc(trexio->alloc, sizeof(double)  * trexio->num_prims);
    double*  coefficient   = (double*) md_alloc(trexio->alloc, sizeof(double)  * trexio->num_prims);
    double*  prim_factor   = (double*) md_alloc(trexio->alloc, sizeof(double)  * trexio->num_prims);

    result = nucleus_index && ang_mom && r_power && shell_factor && shell_index && exponent && coefficient && prim_factor;
    if (result) {
        for (size_t i = 0; i < trexio->num_shells; ++i) { shell_factor[i] = 1.0; r_power[i] = 0; }
        for (size_t i = 0; i < trexio->num_prims;  ++i) { prim_factor[i]  = 1.0; }

        result = trexio_read_dataset(nucleus_index, group, "basis_nucleus_index", H5T_NATIVE_INT64)
              && trexio_read_dataset(ang_mom,       group, "basis_shell_ang_mom", H5T_NATIVE_INT64)
              && trexio_read_dataset(shell_index,   group, "basis_shell_index",   H5T_NATIVE_INT64)
              && trexio_read_dataset(exponent,      group, "basis_exponent",      H5T_NATIVE_DOUBLE)
              && trexio_read_dataset(coefficient,   group, "basis_coefficient",   H5T_NATIVE_DOUBLE);
    }
    if (result) {
        if (trexio_has(group, "basis_shell_factor")) {
            trexio_read_dataset(shell_factor, group, "basis_shell_factor", H5T_NATIVE_DOUBLE);
        }
        if (trexio_has(group, "basis_prim_factor")) {
            trexio_read_dataset(prim_factor, group, "basis_prim_factor", H5T_NATIVE_DOUBLE);
        } else {
            // TREXIO defines it as 1 when absent, which means the contraction is stated over
            // UNNORMALISED primitives. Every writer this reader has met emits it, so say so rather
            // than silently reading a basis with the wrong radial shape.
            MD_LOG_INFO("TREXIO: basis_prim_factor is absent; the contraction is read as stated, over unnormalised primitives");
        }
        if (trexio_has(group, "basis_r_power")) {
            trexio_read_dataset(r_power, group, "basis_r_power", H5T_NATIVE_INT64);
            for (size_t s = 0; s < trexio->num_shells; ++s) {
                if (r_power[s] != 0) {
                    MD_LOG_ERROR("TREXIO: shell %zu has basis_r_power %d; this library evaluates r^l Gaussians only", s, (int)r_power[s]);
                    result = false;
                    break;
                }
            }
        }
    }

    if (result) {
        // TREXIO states the primitives as a flat list tagged with the shell they belong to, in no
        // promised order. md_gto wants each shell's primitives contiguous, so they are counted
        // first and then placed - one pass each, instead of an insert per primitive.
        md_array_resize(trexio->shell_atom,        trexio->num_shells, trexio->alloc);
        md_array_resize(trexio->shell_l,           trexio->num_shells, trexio->alloc);
        md_array_resize(trexio->shell_prim_offset, trexio->num_shells, trexio->alloc);
        md_array_resize(trexio->shell_prim_count,  trexio->num_shells, trexio->alloc);
        md_array_resize(trexio->alpha,             trexio->num_prims,  trexio->alloc);
        md_array_resize(trexio->coeff,             trexio->num_prims,  trexio->alloc);

        for (size_t s = 0; s < trexio->num_shells; ++s) {
            if (nucleus_index[s] < 0 || (size_t)nucleus_index[s] >= trexio->num_atoms) {
                MD_LOG_ERROR("TREXIO: shell %zu is centred on atom %d, the file has %zu", s, (int)nucleus_index[s], trexio->num_atoms);
                result = false;
                break;
            }
            if (ang_mom[s] < 0 || ang_mom[s] > MD_GTO_MAX_ANGULAR_MOMENTUM) {
                MD_LOG_ERROR("TREXIO: shell %zu has angular momentum %d, above the highest this library evaluates (%d)",
                             s, (int)ang_mom[s], MD_GTO_MAX_ANGULAR_MOMENTUM);
                result = false;
                break;
            }
            trexio->shell_atom[s]       = (uint32_t)nucleus_index[s];
            trexio->shell_l[s]          = (uint32_t)ang_mom[s];
            trexio->shell_prim_count[s] = 0;
        }

        for (size_t p = 0; result && p < trexio->num_prims; ++p) {
            if (shell_index[p] < 0 || (size_t)shell_index[p] >= trexio->num_shells) {
                MD_LOG_ERROR("TREXIO: primitive %zu belongs to shell %d, the file has %zu", p, (int)shell_index[p], trexio->num_shells);
                result = false;
                break;
            }
            trexio->shell_prim_count[shell_index[p]] += 1;
        }

        if (result) {
            uint32_t offset = 0;
            for (size_t s = 0; s < trexio->num_shells; ++s) {
                trexio->shell_prim_offset[s] = offset;
                offset += trexio->shell_prim_count[s];
            }
            uint32_t* cursor = (uint32_t*)md_alloc(trexio->alloc, sizeof(uint32_t) * trexio->num_shells);
            MEMSET(cursor, 0, sizeof(uint32_t) * trexio->num_shells);
            for (size_t p = 0; p < trexio->num_prims; ++p) {
                const size_t s = (size_t)shell_index[p];
                const size_t k = trexio->shell_prim_offset[s] + cursor[s]++;
                trexio->alpha[k] = exponent[p];
                // shell_factor is a property of the shell and prim_factor of the primitive; both are
                // TREXIO's own normalisation and both belong in the radial coefficient.
                trexio->coeff[k] = coefficient[p] * prim_factor[p] * shell_factor[s];
            }
        }
    }

done:
    H5Gclose(group);
    return result;
}

static bool trexio_read_ao(trexio_t* trexio, hid_t file) {
    if (!trexio_has(file, "ao")) {
        return true;
    }
    hid_t group = H5Gopen(file, "ao", H5P_DEFAULT);
    if (group == H5I_INVALID_HID) {
        return true;
    }

    int64_t ao_num = 0;
    int64_t cart   = 0;
    trexio_read_i64_attr(&ao_num, group, "ao_num");
    if (trexio_read_i64_attr(&cart, group, "ao_cartesian")) {
        trexio->cartesian = (cart != 0);
    }
    trexio->num_ao = (ao_num > 0) ? (size_t)ao_num : 0;

    if (trexio->num_ao > 0) {
        md_array_resize(trexio->ao_scale, trexio->num_ao, trexio->alloc);
        for (size_t i = 0; i < trexio->num_ao; ++i) {
            trexio->ao_scale[i] = 1.0;
        }
        if (trexio_has(group, "ao_normalization") && trexio_dataset_count(group, "ao_normalization") == trexio->num_ao) {
            trexio_read_dataset(trexio->ao_scale, group, "ao_normalization", H5T_NATIVE_DOUBLE);
        }
    }

    H5Gclose(group);
    return true;
}

static void trexio_read_mo(trexio_t* trexio, hid_t file) {
    if (!trexio_has(file, "mo") || trexio->num_ao == 0) {
        return;
    }
    hid_t group = H5Gopen(file, "mo", H5P_DEFAULT);
    if (group == H5I_INVALID_HID) {
        return;
    }

    int64_t mo_num = 0;
    trexio->mo_type = trexio_read_str_attr(group, "mo_type", trexio->alloc);
    trexio_read_i64_attr(&mo_num, group, "mo_num");
    if (mo_num <= 0 || !trexio_has(group, "mo_coefficient")) {
        H5Gclose(group);
        return;
    }
    trexio->num_mo = (size_t)mo_num;

    md_array_resize(trexio->mo_coefficient, trexio->num_mo * trexio->num_ao, trexio->alloc);
    if (!trexio_read_dataset(trexio->mo_coefficient, group, "mo_coefficient", H5T_NATIVE_DOUBLE)) {
        md_array_free(trexio->mo_coefficient, trexio->alloc);
        trexio->mo_coefficient = NULL;
        trexio->num_mo = 0;
        H5Gclose(group);
        return;
    }

    md_array_resize(trexio->mo_occupation, trexio->num_mo, trexio->alloc);
    md_array_resize(trexio->mo_energy,     trexio->num_mo, trexio->alloc);
    md_array_resize(trexio->mo_spin,       trexio->num_mo, trexio->alloc);
    MEMSET(trexio->mo_occupation, 0, sizeof(double)  * trexio->num_mo);
    MEMSET(trexio->mo_energy,     0, sizeof(double)  * trexio->num_mo);
    MEMSET(trexio->mo_spin,       0, sizeof(int64_t) * trexio->num_mo);

    if (!trexio_has(group, "mo_occupation") || !trexio_read_dataset(trexio->mo_occupation, group, "mo_occupation", H5T_NATIVE_DOUBLE)) {
        md_array_free(trexio->mo_occupation, trexio->alloc);
        trexio->mo_occupation = NULL;
    }
    if (!trexio_has(group, "mo_energy") || !trexio_read_dataset(trexio->mo_energy, group, "mo_energy", H5T_NATIVE_DOUBLE)) {
        md_array_free(trexio->mo_energy, trexio->alloc);
        trexio->mo_energy = NULL;
    }
    if (!trexio_has(group, "mo_spin") || !trexio_read_dataset(trexio->mo_spin, group, "mo_spin", H5T_NATIVE_INT64)) {
        md_array_free(trexio->mo_spin, trexio->alloc);
        trexio->mo_spin = NULL;
    }

    if (trexio_has(group, "mo_symmetry")) {
        md_array_resize(trexio->mo_symmetry, trexio->num_mo, trexio->alloc);
        MEMSET(trexio->mo_symmetry, 0, sizeof(str_t) * trexio->num_mo);
        if (trexio_read_str_dataset(trexio->mo_symmetry, trexio->num_mo, group, "mo_symmetry", trexio->alloc) != trexio->num_mo) {
            md_array_free(trexio->mo_symmetry, trexio->alloc);
            trexio->mo_symmetry = NULL;
        }
    }

    H5Gclose(group);
}

// ---------------------------------------------------------------------------
// Into md_gto's convention
// ---------------------------------------------------------------------------

static bool trexio_normalise_basis(trexio_t* trexio) {
    for (size_t s = 0; s < trexio->num_shells; ++s) {
        const uint32_t l = trexio->shell_l[s];
        double* a = trexio->alpha + trexio->shell_prim_offset[s];
        double* c = trexio->coeff + trexio->shell_prim_offset[s];
        const uint32_t n = trexio->shell_prim_count[s];

        if (n == 0) {
            MD_LOG_ERROR("TREXIO: shell %zu has no primitives", s);
            return false;
        }
        for (uint32_t p = 0; p < n; ++p) {
            if (!(a[p] > 0.0)) {
                MD_LOG_ERROR("TREXIO: shell %zu has a non positive exponent", s);
                return false;
            }
        }
        // NOT md_qm_primitive_norm_factor here, unlike the Molden reader: basis_prim_factor IS the
        // primitive normalisation and is already folded into coeff. Applying a second one would
        // change the RELATIVE weights within the contraction, which the shell normalisation below
        // cannot undo and which looks like a slightly wrong radial shape rather than like an error.
        const double norm = md_qm_shell_norm_factor(l, a, c, n, !trexio->cartesian);
        if (!(norm > 0.0)) {
            MD_LOG_ERROR("TREXIO: shell %zu has a contraction that does not normalise", s);
            return false;
        }
        for (uint32_t p = 0; p < n; ++p) {
            c[p] /= norm;
        }
    }
    return true;
}

// The factor that takes a coefficient stated against the FILE's AO onto md_gto's unit normalised
// one. ao_normalization is TREXIO's per function factor; the Cartesian part is the ratio between
// normalising a shell once, against its (l,0,0) monomial, and normalising every monomial in it -
// see md_qm_cart_coeff_factor. The two compose, and their product is 1 for a file that already
// normalised each Cartesian function individually, which is exactly what ao_normalization is for.
static bool trexio_ao_factors(double* dst, const trexio_t* trexio) {
    size_t ao = 0;
    for (size_t s = 0; s < trexio->num_shells; ++s) {
        const uint32_t l = trexio->shell_l[s];
        const uint32_t n = trexio->cartesian ? md_gto_num_cart_ao(l) : md_gto_num_sph_ao(l);
        for (uint32_t k = 0; k < n; ++k) {
            if (ao >= trexio->num_ao) {
                return false;
            }
            const double scale = md_array_size(trexio->ao_scale) ? trexio->ao_scale[ao] : 1.0;
            dst[ao] = trexio->cartesian ? scale * md_qm_cart_coeff_factor(l, k) : scale;
            ao += 1;
        }
    }
    return ao == trexio->num_ao;
}

static size_t trexio_num_cart_ao(const trexio_t* trexio) {
    size_t n = 0;
    for (size_t s = 0; s < trexio->num_shells; ++s) {
        n += md_gto_num_cart_ao(trexio->shell_l[s]);
    }
    return n;
}

static bool trexio_gto_basis_extract(md_gto_basis_t* out, const trexio_t* trexio, md_allocator_i* alloc) {
    MEMSET(out, 0, sizeof(*out));
    if (trexio->num_shells == 0) {
        return false;
    }
    md_array_resize(out->shells, trexio->num_shells, alloc);
    md_array_resize(out->alpha,  trexio->num_prims,  alloc);
    md_array_resize(out->coeff,  trexio->num_prims,  alloc);
    for (size_t s = 0; s < trexio->num_shells; ++s) {
        out->shells[s] = (md_gto_shell_t){
            .atom_idx         = trexio->shell_atom[s],
            .primitive_offset = trexio->shell_prim_offset[s],
            .num_primitives   = trexio->shell_prim_count[s],
            .l                = trexio->shell_l[s],
        };
    }
    for (size_t p = 0; p < trexio->num_prims; ++p) {
        out->alpha[p] = (float)trexio->alpha[p];
        out->coeff[p] = (float)trexio->coeff[p];
    }
    out->num_shells     = (uint32_t)trexio->num_shells;
    out->num_primitives = (uint32_t)trexio->num_prims;
    return true;
}

// md_gto's AO index for each of the file's, as one array over the whole basis. For a spherical file
// that is a reordering within each shell; for a Cartesian one the two orders already agree -
// TREXIO's canonical Cartesian ordering IS md_gto's - and it comes out as the identity.
static bool trexio_ao_permutation(uint32_t* dst, const trexio_t* trexio) {
    size_t off = 0;
    for (size_t s = 0; s < trexio->num_shells; ++s) {
        const uint32_t l = trexio->shell_l[s];
        const uint32_t n = trexio->cartesian ? md_gto_num_cart_ao(l) : md_gto_num_sph_ao(l);
        for (uint32_t k = 0; k < n; ++k) {
            if (off + k >= trexio->num_ao) {
                return false;
            }
            dst[off + k] = (uint32_t)off + (trexio->cartesian ? k : trexio_sph_to_md_index(l, k));
        }
        off += n;
    }
    return off == trexio->num_ao;
}

// ---------------------------------------------------------------------------
// Publishing
// ---------------------------------------------------------------------------

static bool trexio_publish_orbitals(md_system_t* sys, const trexio_t* trexio, md_allocator_i* temp) {
    if (trexio->num_mo == 0 || !trexio->mo_coefficient) {
        return true;
    }

    const size_t num_ao   = trexio->num_ao;
    const size_t num_cart = trexio_num_cart_ao(trexio);

    double*   factor = md_alloc(temp, sizeof(double)   * num_ao);
    uint32_t* perm   = md_alloc(temp, sizeof(uint32_t) * num_ao);
    if (!factor || !perm || !trexio_ao_factors(factor, trexio) || !trexio_ao_permutation(perm, trexio)) {
        MD_LOG_ERROR("TREXIO: the basis spans a different number of atomic orbitals than ao_num says (%zu)", num_ao);
        return false;
    }

    // Is this restricted? TREXIO says so through mo_spin, and its absence means the same thing as
    // an all-zero column: one set of orbitals, with the occupation stated over BOTH spins.
    bool unrestricted = false;
    for (size_t i = 0; trexio->mo_spin && i < trexio->num_mo; ++i) {
        unrestricted = unrestricted || trexio->mo_spin[i] != 0;
    }
    // An occupation above 1 is what makes it a TOTAL rather than a channel's own; a file with no
    // beta orbitals and none above 1 is already per channel and is left alone.
    double max_occupation = 0.0;
    for (size_t i = 0; trexio->mo_occupation && i < trexio->num_mo; ++i) {
        max_occupation = MAX(max_occupation, trexio->mo_occupation[i]);
    }
    const bool restricted = !unrestricted && max_occupation > 1.0;

    // The file's own numbers, unaltered, beside the per channel split below.
    md_qm_publish_series(sys, STR_LIT("trexio/mo/occupation"), STR_LIT("Occupation"), md_unit_none(), trexio->mo_occupation, trexio->mo_occupation ? trexio->num_mo : 0);
    if (trexio->mo_spin) {
        md_qm_publish_column(sys, STR_LIT("trexio/mo/spin"), STR_LIT("Spin"), md_unit_none(),
                             MD_ATTRIBUTE_TYPE_I64, trexio->mo_spin, sizeof(int64_t), trexio->num_mo);
    }

    md_gto_basis_t cart_basis = {0};
    if (!trexio->cartesian && !trexio_gto_basis_extract(&cart_basis, trexio, temp)) {
        return false;
    }

    md_attribute_id_t alpha_coeff_id = MD_ATTRIBUTE_INVALID;
    md_attribute_id_t alpha_occ_id   = MD_ATTRIBUTE_INVALID;
    md_attribute_id_t alpha_ener_id  = MD_ATTRIBUTE_INVALID;
    md_attribute_id_t alpha_sym_id   = MD_ATTRIBUTE_INVALID;

    for (int spin = 0; spin < 2; ++spin) {
        const bool beta = (spin == 1);
        if (beta && !unrestricted) {
            // One set of orbitals under two names: no copy, and the densities then come out as
            // twice alpha and zero without a special case anywhere.
            md_qm_alias(sys, alpha_coeff_id, STR_LIT("orbital/beta/coefficient"));
            md_qm_alias(sys, alpha_occ_id,   STR_LIT("orbital/beta/occupation"));
            md_qm_alias(sys, alpha_ener_id,  STR_LIT("orbital/beta/energy"));
            md_qm_alias(sys, alpha_sym_id,   STR_LIT("orbital/beta/symmetry"));
            break;
        }

        size_t count = 0;
        for (size_t i = 0; i < trexio->num_mo; ++i) {
            const bool is_beta = trexio->mo_spin && trexio->mo_spin[i] != 0;
            count += (is_beta == beta) ? 1 : 0;
        }
        if (count == 0) {
            continue;
        }

        double* energy = md_alloc(temp, sizeof(double) * count);
        double* occ    = md_alloc(temp, sizeof(double) * count);
        str_t*  sym    = md_alloc(temp, sizeof(str_t)  * count);
        double* coeff  = md_alloc(temp, sizeof(double) * count * num_cart);
        double* mdlib_ao = md_alloc(temp, sizeof(double) * num_ao);
        if (!energy || !occ || !sym || !coeff || !mdlib_ao) {
            MD_LOG_ERROR("TREXIO: failed to allocate scratch for %zu orbitals over %zu atomic orbitals", count, num_cart);
            return false;
        }

        bool have_symmetry = false;
        size_t m = 0;
        for (size_t i = 0; i < trexio->num_mo; ++i) {
            const bool is_beta = trexio->mo_spin && trexio->mo_spin[i] != 0;
            if (is_beta != beta) {
                continue;
            }
            energy[m] = trexio->mo_energy ? trexio->mo_energy[i] : 0.0;
            // Restricted: the occupation covers both spins, so a channel gets half of it.
            occ[m] = trexio->mo_occupation ? (restricted ? 0.5 * trexio->mo_occupation[i] : trexio->mo_occupation[i]) : 0.0;
            sym[m] = trexio->mo_symmetry ? trexio->mo_symmetry[i] : (str_t){0};
            have_symmetry = have_symmetry || !str_empty(sym[m]);

            for (size_t k = 0; k < num_ao; ++k) {
                mdlib_ao[perm[k]] = trexio->mo_coefficient[i * num_ao + k] * factor[k];
            }
            if (trexio->cartesian) {
                MEMCPY(coeff + m * num_cart, mdlib_ao, sizeof(double) * num_cart);
            } else if (md_qm_sph_to_cart_coefficients(coeff + m * num_cart, mdlib_ao, 1, &cart_basis) != 1) {
                MD_LOG_ERROR("TREXIO: could not expand orbital %zu into the Cartesian basis", i);
                return false;
            }
            m += 1;
        }

        const str_t energy_path = beta ? STR_LIT("orbital/beta/energy")      : STR_LIT("orbital/alpha/energy");
        const str_t occ_path    = beta ? STR_LIT("orbital/beta/occupation")  : STR_LIT("orbital/alpha/occupation");
        const str_t sym_path    = beta ? STR_LIT("orbital/beta/symmetry")    : STR_LIT("orbital/alpha/symmetry");
        const str_t coeff_path  = beta ? STR_LIT("orbital/beta/coefficient") : STR_LIT("orbital/alpha/coefficient");

        // Only when the file actually stated them: an absent path is how a consumer learns a block
        // is missing, and a column of zeros would read as a legitimate set of orbital energies.
        md_attribute_id_t ener_id = trexio->mo_energy ? md_qm_publish_series(sys, energy_path, STR_LIT("Energy"), md_unit_hartree(), energy, count) : MD_ATTRIBUTE_INVALID;
        md_attribute_id_t occ_id  = md_qm_publish_series(sys, occ_path,    STR_LIT("Occupation"), md_unit_none(),    occ,    count);
        md_attribute_id_t sym_id  = have_symmetry ? md_qm_publish_strings(sys, sym_path, STR_LIT("Symmetry"), sym, count) : MD_ATTRIBUTE_INVALID;

        md_attribute_id_t coeff_id = md_qm_publish_matrix(sys, coeff_path,
            beta ? STR_LIT("Beta Coefficient") : STR_LIT("Alpha Coefficient"),
            md_unit_none(), coeff, count, num_cart);

        if (!beta) {
            alpha_coeff_id = coeff_id;
            alpha_occ_id   = occ_id;
            alpha_ener_id  = ener_id;
            alpha_sym_id   = sym_id;
        }
    }

    return true;
}

// The AO overlap, in the SAME Cartesian order and convention as the coefficients above - so the two
// can be used together without further ceremony. Read the AO CONVENTION block in md_gto.h first:
// the Cartesian embedding of a spherical basis is rank deficient, so for a file that stored
// spherical data this matrix is SINGULAR. That is fine for Mulliken partitioning and tr(DS), and
// fatal for anything that inverts or factorises it.
static bool trexio_system_begin(md_system_t* sys, md_system_state_t* state, const trexio_t* trexio) {
    if (!sys->alloc) {
        MD_LOG_ERROR("System allocator not set");
        return false;
    }
    if (!state || !state->alloc) {
        MD_LOG_ERROR("State allocator not set");
        return false;
    }

    const size_t num_atoms = trexio->num_atoms;
    md_system_reset(sys);
    md_system_state_init(state, num_atoms);

    const size_t capacity = ROUND_UP(num_atoms, 16);
    md_array_resize(sys->atom.type_idx, capacity, sys->alloc);
    md_array_resize(sys->atom.flags,    capacity, sys->alloc);
    MEMSET(sys->atom.type_idx, 0, md_array_bytes(sys->atom.type_idx));
    MEMSET(sys->atom.flags,    0, md_array_bytes(sys->atom.flags));

    md_atom_type_find_or_add(&sys->atom.type, STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, sys->alloc);

    for (size_t i = 0; i < num_atoms; ++i) {
        state->x[i] = (float)trexio->coord[i].x;
        state->y[i] = (float)trexio->coord[i].y;
        state->z[i] = (float)trexio->coord[i].z;

        const md_atomic_number_t z = trexio->atomic_number[i];
        sys->atom.type_idx[i] = md_atom_type_find_or_add(&sys->atom.type, md_atomic_number_symbol(z), z,
                                                         md_atomic_number_mass(z), md_atomic_number_vdw_radius(z),
                                                         md_atomic_number_cpk_color(z), 0, sys->alloc);
    }

    sys->atom.count  = num_atoms;
    state->num_atoms = num_atoms;
    return true;
}

static bool trexio_publish(md_system_t* sys, const trexio_t* trexio, md_allocator_i* temp) {
    if (!sys->attributes.alloc) {
        MD_LOG_ERROR("Attribute table allocator not set; the system has not been initialised");
        return false;
    }

    md_qm_publish_str(sys, STR_LIT("trexio/metadata/code"),            STR_LIT("Code"),            trexio->code);
    md_qm_publish_str(sys, STR_LIT("trexio/metadata/package_version"), STR_LIT("Package Version"), trexio->package_version);
    md_qm_publish_str(sys, STR_LIT("trexio/basis_type"),               STR_LIT("Basis Type"),      trexio->basis_type);
    md_qm_publish_str(sys, STR_LIT("trexio/mo_type"),                  STR_LIT("MO Type"),         trexio->mo_type);

    if (trexio->has_electron_count) {
        md_qm_publish_scalar(sys, STR_LIT("trexio/electron_count/total"), STR_LIT("Electrons"),       md_unit_none(), (double)trexio->electron_num);
        md_qm_publish_scalar(sys, STR_LIT("trexio/electron_count/up"),    STR_LIT("Alpha Electrons"), md_unit_none(), (double)trexio->electron_up_num);
        md_qm_publish_scalar(sys, STR_LIT("trexio/electron_count/down"),  STR_LIT("Beta Electrons"),  md_unit_none(), (double)trexio->electron_dn_num);
    }
    if (trexio->has_nuclear_repulsion) {
        md_qm_publish_scalar(sys, STR_LIT("trexio/nuclear_repulsion_energy"), STR_LIT("Nuclear Repulsion Energy"), md_unit_hartree(), trexio->nuclear_repulsion);
    }

    md_qm_publish_atoms(sys, trexio->atomic_number, trexio->coord, trexio->num_atoms);

    if (trexio->num_shells > 0) {
        md_qm_publish_str(sys, STR_LIT("trexio/ao_convention"), STR_LIT("AO Convention"),
                           trexio->cartesian ? STR_LIT("cartesian") : STR_LIT("spherical"));

        md_gto_basis_t basis = {0};
        if (trexio_gto_basis_extract(&basis, trexio, temp)) {
            md_qm_publish_basis(sys, &basis);
        }
        if (!trexio_publish_orbitals(sys, trexio, temp)) {
            return false;
        }
        md_qm_publish_overlap(sys);
        md_qm_publish_orbital_densities(sys);
    }

    return true;
}

// ---------------------------------------------------------------------------
// Entry points
// ---------------------------------------------------------------------------

bool md_trexio_system_init_from_file(md_system_t* sys, md_system_state_t* state, str_t filename) {
    ASSERT(sys);

    char path[2048];
    str_copy_to_char_buf(path, sizeof(path), filename);

    h5_error_scope_t err = h5_error_scope_begin();
    hid_t file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file == H5I_INVALID_HID) {
        MD_LOG_ERROR("TREXIO: could not open '" STR_FMT "' as an HDF5 file", STR_ARG(filename));
        h5_error_scope_end(err);
        return false;
    }

    md_temp_scope_t temp_scope = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp_scope);

    trexio_t trexio = { .alloc = temp_arena };
    bool result = trexio_read_nucleus(&trexio, file);
    if (result) {
        trexio_read_metadata(&trexio, file);
        trexio_read_electron(&trexio, file);
        result = trexio_read_basis(&trexio, file);
    }
    if (result) {
        result = trexio_read_ao(&trexio, file);
    }
    if (result) {
        trexio_read_mo(&trexio, file);
    }

    H5Fclose(file);
    h5_error_scope_end(err);

    result = result
          && (trexio.num_shells == 0 || trexio_normalise_basis(&trexio))
          && trexio_system_begin(sys, state, &trexio)
          && trexio_publish(sys, &trexio, temp_arena);

    md_temp_end(temp_scope);
    return result;
}

bool md_trexio_file_is_trexio(str_t filename) {
    char path[2048];
    str_copy_to_char_buf(path, sizeof(path), filename);

    h5_error_scope_t err = h5_error_scope_begin();
    bool result = false;
    if (H5Fis_hdf5(path) > 0) {
        hid_t file = H5Fopen(path, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file != H5I_INVALID_HID) {
            result = trexio_has(file, "nucleus");
            if (result) {
                hid_t group = H5Gopen(file, "nucleus", H5P_DEFAULT);
                int64_t num = 0;
                result = (group != H5I_INVALID_HID) && trexio_read_i64_attr(&num, group, "nucleus_num") && num > 0;
                if (group != H5I_INVALID_HID) H5Gclose(group);
            }
            H5Fclose(file);
        }
    }
    h5_error_scope_end(err);
    return result;
}
