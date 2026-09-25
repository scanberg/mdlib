#include <md_tpr.h>
#include <md_xdr.h>

#include <md_system.h>
#include <md_util.h>
#include <md_nonbonded.h>

#include <core/md_common.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_array.h>
#include <core/md_hash.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <core/md_str.h>
#include <core/md_unit.h>

#include <stdlib.h>
#include <math.h>

// The layout follows GROMACS' own reader (src/gromacs/fileio/tpxio.cpp), which is the only
// specification the format has. The names of the version constants below are GROMACS' tpxv_*
// names, so a condition here can be checked against the corresponding one there.
//
// A file is an XDR encoded header followed by a body:
//
//   header   'VERSION ...' string, precision (4 or 8), tpx version and generation, atom count,
//            flags saying which of box / topology / coordinates / velocities / forces /
//            simulation parameters follow, and, since tpx 119, the byte size of the body.
//   body     box, topology, x, v, f, then the simulation parameters. Of those only the first part is
//            read: the periodic boundary type and the fields up to and including the non-bonded
//            interactions (cut-offs, modifiers, dielectric constants, Ewald tolerances).
//
// The simulation parameters are not meant to be read by anything but GROMACS itself: every tpx version
// may add or remove fields, and unlike the topology no generation number promises anything about them.
// So they are read only for versions whose layout is known here (up to TPX_VERSION_IR_MAX), following
// the conditions of GROMACS' do_inputrec field by field.
//
// Before tpx 119 the body is XDR encoded like the header. From 119 on it is written by GROMACS'
// in-memory serializer instead: still big endian, but every value takes its native size (a bool or
// uchar is 1 byte, a ushort 2) with no padding, and a string is a u64 length followed by its bytes.
// Everything else - ints, int64, floats, doubles - is identical in the two encodings.

// tpx versions, GROMACS' tpxv_* names
enum {
    tpxv_Pre96Version51 = 51,
    tpxv_Pre96Version53 = 53,
    tpxv_Pre96Version56 = 56,
    tpxv_Pre96Version57 = 57,                   // Oldest version GROMACS can not read
    tpxv_Pre96Version59 = 59,
    tpxv_Pre96Version60 = 60,
    tpxv_Pre96Version61 = 61,
    tpxv_Pre96Version62 = 62,
    tpxv_Pre96Version63 = 63,
    tpxv_Pre96Version65 = 65,
    tpxv_Pre96Version66 = 66,
    tpxv_Pre96Version67 = 67,
    tpxv_Pre96Version68 = 68,
    tpxv_Pre96Version69 = 69,
    tpxv_Pre96Version70 = 70,
    tpxv_Pre96Version72 = 72,
    tpxv_Pre96Version76 = 76,
    tpxv_Pre96Version77 = 77,
    tpxv_Pre96Version78 = 78,
    tpxv_Pre96Version79 = 79,
    tpxv_Pre96Version81 = 81,
    tpxv_Pre96Version82 = 82,
    tpxv_Pre96Version90 = 90,
    tpxv_Pre96Version93 = 93,
    tpxv_Pre96Version94 = 94,
    tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials = 98,
    tpxv_RemoveObsoleteParameters1 = 100,
    tpxv_IntermolecularBondeds = 103,
    tpxv_RemoveTwinRange = 108,
    tpxv_RemoveImplicitSolvation = 113,
    tpxv_GenericInternalParameters = 117,
    tpxv_VSite2FD = 118,
    tpxv_AddSizeField = 119,
    tpxv_StoreNonBondedInteractionExclusionGroup = 120,
    tpxv_VSite1 = 121,
    tpxv_MTS = 122,
    tpxv_RemoveTholeRfac = 127,
    tpxv_RemoveAtomtypes = 128,
    tpxv_EnsembleTemperature = 129,
    tpxv_MassRepartitioning = 131,
    tpxv_VerletBufferPressureTol = 133,
    tpxv_HandleMartiniBondedBStateParametersProperly = 134,
    tpxv_NNPotIFuncType = 137,
    tpxv_AwhHistogramTolerance = 138,           // GROMACS 2026
    tpxv_OutputControlInKeyValueTree = 139,
    tpxv_CmapBState = 140,
};

// The newest version whose simulation parameters can be read: the leading part of do_inputrec is known
// up to here. Checked against GROMACS 2026.1 and the development branch after it.
#define TPX_VERSION_IR_MAX tpxv_CmapBState

// tpx topology generations. The topology of a file is readable when its generation is known, even
// if the file version itself is newer than anything listed above.
enum {
    TPX_GENERATION_ADD_SIZE_FIELD = 27,
    TPX_GENERATION_MAX            = 29,
};

// GROMACS interaction function types (InteractionFunction in ifunc.h), in GROMACS' order. The
// order is the file format: the interaction lists of a molecule type are stored one per type, in
// exactly this order, and only types that existed when the file was written are present (see
// ftupd below).
enum {
    F_BONDS, F_G96BONDS, F_MORSE, F_CUBICBONDS, F_CONNBONDS, F_HARMONIC, F_FENEBONDS, F_TABBONDS,
    F_TABBONDSNC, F_RESTRBONDS, F_ANGLES, F_G96ANGLES, F_RESTRANGLES, F_LINEAR_ANGLES,
    F_CROSS_BOND_BONDS, F_CROSS_BOND_ANGLES, F_UREY_BRADLEY, F_QUARTIC_ANGLES, F_TABANGLES, F_PDIHS,
    F_RBDIHS, F_RESTRDIHS, F_CBTDIHS, F_FOURDIHS, F_IDIHS, F_PIDIHS, F_TABDIHS, F_CMAP, F_GB12, F_GB13,
    F_GB14, F_GBPOL, F_NPSOLVATION, F_LJ14, F_COUL14, F_LJC14_Q, F_LJC_PAIRS_NB, F_LJ, F_BHAM, F_LJ_LR,
    F_BHAM_LR, F_DISPCORR, F_COUL_SR, F_COUL_LR, F_RF_EXCL, F_COUL_RECIP, F_LJ_RECIP, F_DPD,
    F_POLARIZATION, F_WATER_POL, F_THOLE_POL, F_ANHARM_POL, F_POSRES, F_FBPOSRES, F_DISRES,
    F_DISRESVIOL, F_ORIRES, F_ORIRESDEV, F_ANGRES, F_ANGRESZ, F_DIHRES, F_DIHRESVIOL, F_CONSTR,
    F_CONSTRNC, F_SETTLE, F_VSITE1, F_VSITE2, F_VSITE2FD, F_VSITE3, F_VSITE3FD, F_VSITE3FAD,
    F_VSITE3OUT, F_VSITE4FD, F_VSITE4FDN, F_VSITEN, F_COM_PULL, F_DENSITYFITTING, F_EQM, F_ENNPOT,
    F_EPOT, F_EKIN, F_ETOT, F_ECONSERVED, F_TEMP, F_VTEMP, F_PDISPCORR, F_PRES, F_DVDL_CONSTR, F_DVDL,
    F_DKDL, F_DVDL_COUL, F_DVDL_VDW, F_DVDL_BONDED, F_DVDL_RESTRAINT, F_DVDL_TEMPERATURE,
    F_NRE
};
STATIC_ASSERT(F_NRE == 95, "The interaction function list must match GROMACS'");

// Function types added after the oldest supported version: a file older than 'version' does not
// have the type, so its interaction list is absent and every type number from 'ftype' up is one
// lower in the file than it is here. Ascending by type, then by version, as GROMACS' ftupd.
typedef struct ftupd_t {
    int version;
    int ftype;
} ftupd_t;

static const ftupd_t ftupd[] = {
    { tpxv_Pre96Version70, F_RESTRBONDS },
    { tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, F_RESTRANGLES },
    { tpxv_Pre96Version76, F_LINEAR_ANGLES },
    { tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, F_RESTRDIHS },
    { tpxv_RestrictedBendingAndCombinedAngleTorsionPotentials, F_CBTDIHS },
    { tpxv_Pre96Version65, F_CMAP },
    { tpxv_Pre96Version60, F_GB12 },
    { tpxv_Pre96Version61, F_GB13 },
    { tpxv_Pre96Version61, F_GB14 },
    { tpxv_Pre96Version72, F_GBPOL },
    { tpxv_Pre96Version72, F_NPSOLVATION },
    { tpxv_Pre96Version93, F_LJ_RECIP },
    { tpxv_Pre96Version76, F_ANHARM_POL },
    { tpxv_Pre96Version90, F_FBPOSRES },
    { tpxv_VSite1, F_VSITE1 },
    { tpxv_VSite2FD, F_VSITE2FD },
    { tpxv_GenericInternalParameters, F_DENSITYFITTING },
    { tpxv_NNPotIFuncType, F_ENNPOT },
    { tpxv_Pre96Version69, F_VTEMP },
    { tpxv_Pre96Version66, F_PDISPCORR },
    { tpxv_Pre96Version79, F_DVDL_COUL },
    { tpxv_Pre96Version79, F_DVDL_VDW },
    { tpxv_Pre96Version79, F_DVDL_BONDED },
    { tpxv_Pre96Version79, F_DVDL_RESTRAINT },
    { tpxv_Pre96Version79, F_DVDL_TEMPERATURE },
};

// Number of SimulationAtomGroupType entries (temperature coupling, energy output, ... QM/MM)
#define TPR_NUM_ATOM_GROUP_TYPES 10

// Sanity bound for any string in the file (names, the version string)
#define TPR_MAX_STRING_LENGTH (1 << 16)

// ### READER ###
// A cursor that knows which of the two encodings it is reading.

typedef struct tpr_reader_t {
    md_xdr_t xdr;
    bool     in_memory;     // Body written by the in-memory serializer (tpx >= 119)
    bool     dbl;           // Reals are doubles
    int      version;
} tpr_reader_t;

static inline bool rd_ok(const tpr_reader_t* r) {
    return md_xdr_ok(&r->xdr);
}

static inline void rd_fail(tpr_reader_t* r) {
    r->xdr.error = true;
}

static inline int32_t rd_int(tpr_reader_t* r) {
    int32_t v;
    md_xdr_read_i32(&r->xdr, &v);
    return v;
}

// An int that is a count or a size: negative means the file is corrupt
static inline size_t rd_count(tpr_reader_t* r) {
    const int32_t v = rd_int(r);
    if (v < 0) {
        rd_fail(r);
        return 0;
    }
    return (size_t)v;
}

static inline int64_t rd_int64(tpr_reader_t* r) {
    int64_t v;
    md_xdr_read_i64(&r->xdr, &v);
    return v;
}

static inline bool rd_bool(tpr_reader_t* r) {
    if (r->in_memory) {
        const uint8_t* p = md_xdr_take(&r->xdr, 1);
        return p && *p;
    }
    return rd_int(r) != 0;
}

static inline uint8_t rd_uchar(tpr_reader_t* r) {
    if (r->in_memory) {
        const uint8_t* p = md_xdr_take(&r->xdr, 1);
        return p ? *p : 0;
    }
    return (uint8_t)rd_int(r);
}

static inline uint16_t rd_ushort(tpr_reader_t* r) {
    if (r->in_memory) {
        const uint8_t* p = md_xdr_take(&r->xdr, 2);
        return p ? (uint16_t)((p[0] << 8) | p[1]) : 0;
    }
    return (uint16_t)rd_int(r);
}

static inline double rd_real(tpr_reader_t* r) {
    if (r->dbl) {
        double v;
        md_xdr_read_f64(&r->xdr, &v);
        return v;
    }
    float v;
    md_xdr_read_f32(&r->xdr, &v);
    return v;
}

static inline void skip_bytes(tpr_reader_t* r, size_t count, size_t size) {
    if (count > SIZE_MAX / size) {
        rd_fail(r);
        return;
    }
    md_xdr_skip(&r->xdr, count * size);
}

static inline void skip_ints (tpr_reader_t* r, size_t count) { skip_bytes(r, count, 4); }
static inline void skip_reals(tpr_reader_t* r, size_t count) { skip_bytes(r, count, r->dbl ? 8 : 4); }
static inline void skip_uchars(tpr_reader_t* r, size_t count) { skip_bytes(r, count, r->in_memory ? 1 : 4); }

// Reads count reals into dst as floats
static void rd_reals(tpr_reader_t* r, float* dst, size_t count) {
    const size_t size = r->dbl ? 8 : 4;
    const uint8_t* p = (count <= SIZE_MAX / size) ? md_xdr_take(&r->xdr, count * size) : NULL;
    if (!p) {
        rd_fail(r);
        MEMSET(dst, 0, count * sizeof(float));
        return;
    }
    if (r->dbl) {
        for (size_t i = 0; i < count; ++i) dst[i] = (float)md_xdr_load_f64(p + i * 8);
    } else {
        md_xdr_load_f32_array(dst, p, count);
    }
}

// The string is a view into the buffer
static str_t rd_string(tpr_reader_t* r) {
    str_t str = {0};
    if (r->in_memory) {
        uint64_t len;
        md_xdr_read_u64(&r->xdr, &len);
        const uint8_t* p = (len <= TPR_MAX_STRING_LENGTH) ? md_xdr_take(&r->xdr, (size_t)len) : NULL;
        if (!p) {
            rd_fail(r);
            return str;
        }
        str = (str_t){ (const char*)p, (size_t)len };
    } else {
        // GROMACS writes the length including the terminating zero, followed by an XDR string
        // (length without it, then the characters padded to a multiple of 4)
        rd_int(r);
        md_xdr_read_string(&r->xdr, &str, TPR_MAX_STRING_LENGTH);
    }
    // Writers have been known to include the terminating zero in the characters
    while (str.len > 0 && str.ptr[str.len - 1] == '\0') {
        str.len -= 1;
    }
    return str;
}

// ### TOPOLOGY ###

typedef struct tpr_symtab_t {
    size_t count;
    str_t* str;
} tpr_symtab_t;

static str_t rd_symstr(tpr_reader_t* r, const tpr_symtab_t* symtab) {
    const int32_t idx = rd_int(r);
    if (idx < 0 || (size_t)idx >= symtab->count) {
        rd_fail(r);
        return (str_t){0};
    }
    return symtab->str[idx];
}

static inline bool ftype_absent(int ftype, int version) {
    for (size_t k = 0; k < ARRAY_SIZE(ftupd); ++k) {
        if (version < ftupd[k].version && ftype == ftupd[k].ftype) {
            return true;
        }
    }
    return false;
}

// The parameters of one interaction type are not needed, only their size: ints and reals
static bool skip_iparams(tpr_reader_t* r, int ftype) {
    const int v = r->version;
    size_t ni = 0;
    size_t nr = 0;

    switch (ftype) {
    case F_ANGLES:
    case F_G96ANGLES:
    case F_BONDS:
    case F_G96BONDS:
    case F_HARMONIC:
    case F_IDIHS:
    case F_LINEAR_ANGLES:
        nr = 4;
        break;
    case F_RESTRANGLES:
    case F_RESTRDIHS:
        nr = (v < tpxv_HandleMartiniBondedBStateParametersProperly) ? 2 : 4;
        break;
    case F_FENEBONDS:
        nr = 2;
        break;
    case F_RESTRBONDS:
        nr = 8;
        break;
    case F_TABBONDS:
    case F_TABBONDSNC:
    case F_TABANGLES:
    case F_TABDIHS:
        nr = 2; ni = 1;
        break;
    case F_CROSS_BOND_BONDS:
        nr = 3;
        break;
    case F_CROSS_BOND_ANGLES:
        nr = 4;
        break;
    case F_UREY_BRADLEY:
        nr = (v >= tpxv_Pre96Version79) ? 8 : 4;
        break;
    case F_QUARTIC_ANGLES:
        nr = 6;
        break;
    case F_BHAM:
        nr = 3;
        break;
    case F_MORSE:
        nr = (v >= tpxv_Pre96Version79) ? 6 : 3;
        break;
    case F_CUBICBONDS:
        nr = 3;
        break;
    case F_CONNBONDS:
        break;
    case F_POLARIZATION:
        nr = 1;
        break;
    case F_ANHARM_POL:
        nr = 3;
        break;
    case F_WATER_POL:
        nr = 6;
        break;
    case F_THOLE_POL:
        nr = (v < tpxv_RemoveTholeRfac) ? 4 : 3;
        break;
    case F_LJ:
        nr = 2;
        break;
    case F_LJ14:
        nr = 4;
        break;
    case F_LJC14_Q:
        nr = 5;
        break;
    case F_LJC_PAIRS_NB:
        nr = 4;
        break;
    case F_PDIHS:
    case F_PIDIHS:
    case F_ANGRES:
    case F_ANGRESZ:
        nr = 4; ni = 1;
        break;
    case F_DISRES:
        ni = 2; nr = 4;
        break;
    case F_ORIRES:
        ni = 3; nr = 3;
        break;
    case F_DIHRES:
        if (v < tpxv_Pre96Version82) {
            ni = 2; nr = 3;
        } else {
            nr = 6;
        }
        break;
    case F_POSRES:
        nr = 12;
        break;
    case F_FBPOSRES:
        ni = 1; nr = 5;
        break;
    case F_CBTDIHS:
        nr = (v < tpxv_HandleMartiniBondedBStateParametersProperly) ? 6 : 12;
        break;
    case F_RBDIHS:
    case F_FOURDIHS:
        nr = 12;
        break;
    case F_CONSTR:
    case F_CONSTRNC:
    case F_SETTLE:
        nr = 2;
        break;
    case F_VSITE1:
        break;
    case F_VSITE2:
    case F_VSITE2FD:
        nr = 1;
        break;
    case F_VSITE3:
    case F_VSITE3FD:
    case F_VSITE3FAD:
        nr = 2;
        break;
    case F_VSITE3OUT:
    case F_VSITE4FD:
    case F_VSITE4FDN:
        nr = 3;
        break;
    case F_VSITEN:
        ni = 1; nr = 1;
        break;
    case F_GB12:
    case F_GB13:
    case F_GB14:
        if (v < tpxv_Pre96Version68)          nr += 4;
        if (v < tpxv_RemoveImplicitSolvation) nr += 5;
        break;
    case F_CMAP:
        ni = 2;
        break;
    default:
        MD_LOG_ERROR("TPR: Unknown interaction function type %d", ftype);
        rd_fail(r);
        return false;
    }

    skip_ints(r, ni);
    skip_reals(r, nr);
    return true;
}

// Number of atoms of the virtual site types, the first being the site itself
static int vsite_num_atoms(int ftype) {
    switch (ftype) {
    case F_VSITE1:    return 2;
    case F_VSITE2:
    case F_VSITE2FD:  return 3;
    case F_VSITE3:
    case F_VSITE3FD:
    case F_VSITE3FAD:
    case F_VSITE3OUT: return 4;
    case F_VSITE4FD:
    case F_VSITE4FDN: return 5;
    case F_VSITEN:    return 2;
    default:          return 0;
    }
}

// The two atom types GROMACS flags as chemical bonds (IF_CHEMBOND). SETTLE is the third kind.
static bool is_chemical_bond(int ftype) {
    switch (ftype) {
    case F_BONDS:
    case F_G96BONDS:
    case F_MORSE:
    case F_CUBICBONDS:
    case F_CONNBONDS:
    case F_FENEBONDS:
    case F_TABBONDS:
    case F_POLARIZATION:
    case F_ANHARM_POL:
    case F_CONSTR:
        return true;
    default:
        return false;
    }
}

static inline void push_pair(md_array(md_atom_pair_t)* pairs, int32_t a, int32_t b, size_t num_atoms, tpr_reader_t* r, md_allocator_i* alloc) {
    if (a < 0 || b < 0 || (size_t)a >= num_atoms || (size_t)b >= num_atoms) {
        rd_fail(r);
        return;
    }
    if (a == b) return;
    md_atom_pair_t pair = {{ MIN(a, b), MAX(a, b) }};
    md_array_push(*pairs, pair, alloc);
}

// Reads the interaction lists of a molecule type (or the intermolecular ones) and keeps the pairs
// that make up the connectivity. Atom indices are validated against num_atoms.
static void rd_ilists(tpr_reader_t* r, md_array(md_atom_pair_t)* pairs, size_t num_atoms, md_allocator_i* alloc) {
    for (int ftype = 0; ftype < F_NRE && rd_ok(r); ++ftype) {
        if (ftype_absent(ftype, r->version)) {
            continue;
        }

        const size_t n = rd_count(r);
        const uint8_t* p = (n <= SIZE_MAX / 4) ? md_xdr_take(&r->xdr, n * 4) : NULL;
        if (!p) {
            rd_fail(r);
            return;
        }
#define IATOM(i) md_xdr_load_i32(p + (size_t)(i) * 4)

        if (ftype == F_SETTLE) {
            // A SETTLE is one water: the oxygen and its two hydrogens. Files older than tpx 78
            // store only the oxygen, the hydrogens being the two atoms after it.
            const bool old = r->version < tpxv_Pre96Version78;
            const size_t stride = old ? 2 : 4;
            for (size_t i = 0; i + stride <= n; i += stride) {
                const int32_t o  = IATOM(i + 1);
                const int32_t h1 = old ? o + 1 : IATOM(i + 2);
                const int32_t h2 = old ? o + 2 : IATOM(i + 3);
                push_pair(pairs, o, h1, num_atoms, r, alloc);
                push_pair(pairs, o, h2, num_atoms, r, alloc);
            }
        } else if (is_chemical_bond(ftype)) {
            for (size_t i = 0; i + 3 <= n; i += 3) {
                push_pair(pairs, IATOM(i + 1), IATOM(i + 2), num_atoms, r, alloc);
            }
        } else if (vsite_num_atoms(ftype)) {
            // The site to its first constructing atom. An N-site is stored as one entry per
            // constructing atom, and only the first of those is wanted.
            const size_t stride = 1 + (size_t)vsite_num_atoms(ftype);
            int32_t prev_site = -1;
            for (size_t i = 0; i + stride <= n; i += stride) {
                const int32_t site = IATOM(i + 1);
                if (ftype == F_VSITEN && site == prev_site) continue;
                push_pair(pairs, site, IATOM(i + 2), num_atoms, r, alloc);
                prev_site = site;
            }
        }
#undef IATOM
    }
}

static int compare_pair(const void* a, const void* b) {
    const md_atom_pair_t* pa = (const md_atom_pair_t*)a;
    const md_atom_pair_t* pb = (const md_atom_pair_t*)b;
    if (pa->idx[0] != pb->idx[0]) return pa->idx[0] < pb->idx[0] ? -1 : 1;
    if (pa->idx[1] != pb->idx[1]) return pa->idx[1] < pb->idx[1] ? -1 : 1;
    return 0;
}

// Sorts and removes duplicates (a bond and a constraint between the same pair, say)
static size_t sort_unique_pairs(md_atom_pair_t* pairs, size_t count) {
    if (count < 2) return count;
    qsort(pairs, count, sizeof(md_atom_pair_t), compare_pair);
    size_t out = 1;
    for (size_t i = 1; i < count; ++i) {
        if (compare_pair(&pairs[i], &pairs[out - 1]) != 0) {
            pairs[out++] = pairs[i];
        }
    }
    return out;
}

static void rd_moltype(tpr_reader_t* r, md_tpr_moltype_t* mt, const tpr_symtab_t* symtab, md_allocator_i* alloc) {
    mt->name = rd_symstr(r, symtab);

    const size_t num_atoms = rd_count(r);
    const size_t num_res   = rd_count(r);
    // Every atom takes at least 36 bytes in the file, so anything claiming more than fits is corrupt
    if (!rd_ok(r) || num_atoms > md_xdr_remaining(&r->xdr) / 36 || num_res > num_atoms) {
        rd_fail(r);
        return;
    }

    mt->atoms = md_array_create(md_tpr_atom_t, num_atoms, alloc);
    mt->residues = md_array_create(md_tpr_residue_t, num_res, alloc);
    mt->num_atoms = num_atoms;
    mt->num_residues = num_res;
    if ((num_atoms && !mt->atoms) || (num_res && !mt->residues)) {
        rd_fail(r);
        return;
    }

    for (size_t i = 0; i < num_atoms; ++i) {
        md_tpr_atom_t* atom = &mt->atoms[i];
        atom->mass   = (float)rd_real(r);
        atom->charge = (float)rd_real(r);
        skip_reals(r, 2);                           // B state mass and charge
        atom->type_idx = rd_ushort(r);
        rd_ushort(r);                               // B state type
        atom->ptype = (uint8_t)rd_int(r);
        atom->residue = rd_int(r);
        atom->atomic_number = rd_int(r);
        if (atom->residue < 0 || (size_t)atom->residue >= num_res) {
            rd_fail(r);
        }
    }
    for (size_t i = 0; i < num_atoms; ++i) mt->atoms[i].name = rd_symstr(r, symtab);
    for (size_t i = 0; i < num_atoms; ++i) mt->atoms[i].type = rd_symstr(r, symtab);
    for (size_t i = 0; i < num_atoms; ++i) rd_symstr(r, symtab);    // B state type name

    for (size_t i = 0; i < num_res; ++i) {
        md_tpr_residue_t* res = &mt->residues[i];
        res->name = rd_symstr(r, symtab);
        if (r->version >= tpxv_Pre96Version63) {
            res->nr = rd_int(r);
            const uint8_t ic = rd_uchar(r);
            res->ic = (ic >= 32 && ic < 127) ? (char)ic : ' ';
        } else {
            res->nr = (int32_t)i + 1;
            res->ic = ' ';
        }
    }

    md_array(md_atom_pair_t) bonds = 0;
    rd_ilists(r, &bonds, num_atoms, alloc);
    mt->num_bonds = sort_unique_pairs(bonds, md_array_size(bonds));
    if (bonds) md_array_shrink(bonds, mt->num_bonds);
    mt->bonds = bonds;

    // Charge groups (obsolete but still written)
    const size_t num_cg = rd_count(r);
    skip_ints(r, num_cg + 1);

    // Exclusions: a list of lists, one list per atom, each including the atom itself
    const size_t num_lists = rd_count(r);
    const size_t num_elem  = rd_count(r);
    if (!rd_ok(r) || num_lists > num_atoms || num_elem > md_xdr_remaining(&r->xdr) / 4) {
        rd_fail(r);
        return;
    }
    const uint8_t* ranges = md_xdr_take(&r->xdr, (num_lists + 1) * 4);
    const uint8_t* elems  = md_xdr_take(&r->xdr, num_elem * 4);
    if (!ranges || (num_elem && !elems)) {
        rd_fail(r);
        return;
    }
    if (num_elem > 0 && num_atoms) {
        uint32_t* off  = md_alloc(alloc, sizeof(uint32_t) * (num_atoms + 1));
        uint32_t* excl = md_alloc(alloc, sizeof(uint32_t) * num_elem);
        uint32_t n = 0;
        off[0] = 0;
        for (size_t i = 0; i < num_atoms; ++i) {
            if (i < num_lists) {
                const int32_t beg = md_xdr_load_i32(ranges + i * 4);
                const int32_t end = md_xdr_load_i32(ranges + (i + 1) * 4);
                if (beg < 0 || end < beg || (size_t)end > num_elem) {
                    rd_fail(r);
                    break;
                }
                const uint32_t row = n;
                for (int32_t k = beg; k < end; ++k) {
                    const int32_t j = md_xdr_load_i32(elems + (size_t)k * 4);
                    if (j < 0 || (size_t)j >= num_atoms) {
                        rd_fail(r);
                        break;
                    }
                    if ((size_t)j == i) continue;
                    // Insertion into the sorted row, which is short
                    uint32_t y = n++;
                    while (y > row && excl[y - 1] > (uint32_t)j) {
                        excl[y] = excl[y - 1];
                        --y;
                    }
                    excl[y] = (uint32_t)j;
                }
            }
            off[i + 1] = n;
        }
        if (rd_ok(r) && n > 0) {
            mt->excl_offset = off;
            mt->excl = excl;
        } else {
            md_free(alloc, off, sizeof(uint32_t) * (num_atoms + 1));
            md_free(alloc, excl, sizeof(uint32_t) * num_elem);
        }
    }
}

static void rd_mtop(tpr_reader_t* r, md_tpr_data_t* data, str_t version_string, md_allocator_i* alloc, md_allocator_i* temp_alloc) {
    // ## Symbol table: every name in the topology, referred to by index from here on
    tpr_symtab_t symtab = {0};
    {
        const size_t count = rd_count(r);
        if (!rd_ok(r) || count > md_xdr_remaining(&r->xdr) / 4) {
            rd_fail(r);
            return;
        }
        str_t* views = md_alloc(temp_alloc, (count + 1) * sizeof(str_t));
        size_t total = version_string.len;
        for (size_t i = 0; i < count && rd_ok(r); ++i) {
            views[i] = rd_string(r);
            total += views[i].len;
        }
        if (!rd_ok(r)) return;

        // The strings are views into the file buffer: copy them into one block the data owns
        data->str_data = md_alloc(alloc, total + 1);
        data->str_size = total + 1;
        char* dst = data->str_data;
        MEMCPY(dst, version_string.ptr, version_string.len);
        data->version_string = (str_t){ dst, version_string.len };
        dst += version_string.len;
        for (size_t i = 0; i < count; ++i) {
            MEMCPY(dst, views[i].ptr, views[i].len);
            views[i].ptr = dst;
            dst += views[i].len;
        }
        *dst = '\0';
        symtab.count = count;
        symtab.str = views;
    }

    data->name = rd_symstr(r, &symtab);

    // ## Force field parameters
    // Skipped, except for the Lennard-Jones parameters of the non-bonded types. grompp writes the
    // non-bonded interactions first, as the full num_nb_types x num_nb_types table of type pairs.
    {
        const size_t num_nb_types = rd_count(r);
        const size_t num_types = rd_count(r);
        if (!rd_ok(r) || num_types > md_xdr_remaining(&r->xdr) / 4) {
            rd_fail(r);
            return;
        }
        int32_t* functype = md_alloc(temp_alloc, (num_types + 1) * sizeof(int32_t));
        for (size_t i = 0; i < num_types; ++i) {
            functype[i] = rd_int(r);
        }
        data->repulsion_power = 12.0;
        if (r->version >= tpxv_Pre96Version66) {
            md_xdr_read_f64(&r->xdr, &data->repulsion_power);   // Always a double
        }
        data->fudge_qq = (float)rd_real(r);
        if (!rd_ok(r) || num_nb_types > num_types) {
            rd_fail(r);
            return;
        }
        const size_t num_pairs = num_nb_types * num_nb_types;
        if (num_pairs > num_types) {
            rd_fail(r);
            return;
        }
        data->lj = md_array_create(md_tpr_lj_t, num_pairs, alloc);
        data->num_nb_types = num_nb_types;
        if (num_pairs) MEMSET(data->lj, 0, num_pairs * sizeof(md_tpr_lj_t));
        bool all_lj = num_pairs > 0;

        for (size_t i = 0; i < num_types && rd_ok(r); ++i) {
            // The type numbers in the file are those of the version that wrote it
            int ftype = functype[i];
            for (size_t k = 0; k < ARRAY_SIZE(ftupd); ++k) {
                if (r->version < ftupd[k].version && ftype >= ftupd[k].ftype) {
                    ftype += 1;
                }
            }
            if (i < num_pairs) {
                if (ftype == F_LJ) {
                    md_tpr_lj_t* lj = &data->lj[i];
                    lj->c6  = (float)rd_real(r);
                    lj->c12 = (float)rd_real(r);
                    continue;
                }
                all_lj = false;
            }
            skip_iparams(r, ftype);
        }
        data->nb_is_lj = all_lj;
        if (!all_lj && num_pairs) {
            // A table of mixed or other forms is not Lennard-Jones: leave nothing half filled in
            MEMSET(data->lj, 0, num_pairs * sizeof(md_tpr_lj_t));
        }
    }

    // ## Molecule types
    {
        const size_t count = rd_count(r);
        if (!rd_ok(r) || count > md_xdr_remaining(&r->xdr) / 8) {
            rd_fail(r);
            return;
        }
        data->moltypes = md_array_create(md_tpr_moltype_t, count, alloc);
        if (count) MEMSET(data->moltypes, 0, count * sizeof(md_tpr_moltype_t));
        data->num_moltypes = count;
        for (size_t i = 0; i < count && rd_ok(r); ++i) {
            rd_moltype(r, &data->moltypes[i], &symtab, alloc);
        }
    }

    // ## Molecule blocks
    size_t num_atoms = 0;
    {
        const size_t count = rd_count(r);
        if (!rd_ok(r) || count > md_xdr_remaining(&r->xdr) / 20) {
            rd_fail(r);
            return;
        }
        data->molblocks = md_array_create(md_tpr_molblock_t, count, alloc);
        data->num_molblocks = count;
        for (size_t i = 0; i < count && rd_ok(r); ++i) {
            md_tpr_molblock_t* mb = &data->molblocks[i];
            mb->moltype = rd_int(r);
            mb->nmol = rd_int(r);
            const int32_t atoms_per_mol = rd_int(r);
            if (mb->moltype < 0 || (size_t)mb->moltype >= data->num_moltypes || mb->nmol < 0 ||
                atoms_per_mol < 0 || (size_t)atoms_per_mol != data->moltypes[mb->moltype].num_atoms) {
                MD_LOG_ERROR("TPR: Molecule block %zu is inconsistent with its molecule type", i);
                rd_fail(r);
                return;
            }
            num_atoms += (size_t)mb->nmol * (size_t)atoms_per_mol;
            // Position restraint reference coordinates, A and B state
            skip_reals(r, rd_count(r) * 3);
            skip_reals(r, rd_count(r) * 3);
        }
    }

    const int32_t natoms = rd_int(r);
    if (!rd_ok(r) || natoms < 0 || (size_t)natoms != num_atoms) {
        MD_LOG_ERROR("TPR: The molecule blocks describe %zu atoms, the topology %d", num_atoms, natoms);
        rd_fail(r);
        return;
    }
    data->num_atoms = num_atoms;

    // ## Intermolecular interactions, in global atom indices
    if (r->version >= tpxv_IntermolecularBondeds) {
        if (rd_bool(r)) {
            md_array(md_atom_pair_t) bonds = 0;
            rd_ilists(r, &bonds, num_atoms, alloc);
            data->num_intermolecular_bonds = sort_unique_pairs(bonds, md_array_size(bonds));
            data->intermolecular_bonds = bonds;
        }
    }

    // ## Everything after this is skipped, but has to be read through to get to the coordinates

    // Atom types (removed in tpx 128)
    if (r->version < tpxv_RemoveAtomtypes) {
        const size_t count = rd_count(r);
        if (r->version < tpxv_RemoveImplicitSolvation) {
            skip_reals(r, 3 * count);
        }
        skip_ints(r, count);                        // Atomic numbers
        if (r->version >= tpxv_Pre96Version60 && r->version < tpxv_RemoveImplicitSolvation) {
            skip_reals(r, 2 * count);
        }
    }

    // CMAP grids, each extent x extent points of 4 reals
    if (r->version >= tpxv_Pre96Version65) {
        const size_t num_grids = rd_count(r);
        const size_t extent = rd_count(r);
        if (extent > 0 && num_grids > SIZE_MAX / 4 / extent / extent) {
            rd_fail(r);
            return;
        }
        skip_reals(r, num_grids * extent * extent * 4);
    }

    // Atom groups: the atom index lists of each group type, the group names, and the per atom
    // group numbers of each group type
    for (int i = 0; i < TPR_NUM_ATOM_GROUP_TYPES; ++i) {
        skip_ints(r, rd_count(r));
    }
    skip_ints(r, rd_count(r));
    for (int i = 0; i < TPR_NUM_ATOM_GROUP_TYPES; ++i) {
        skip_uchars(r, rd_count(r));
    }

    if (r->version >= tpxv_StoreNonBondedInteractionExclusionGroup) {
        const int64_t count = rd_int64(r);
        if (count < 0) {
            rd_fail(r);
            return;
        }
        skip_ints(r, (size_t)count);
    }
}

// ### SIMULATION PARAMETERS ###

// The leading part of GROMACS' do_inputrec, up to the non-bonded interactions and the Ewald parameters,
// which is all that is kept. Field by field in GROMACS' order, under GROMACS' version conditions.
static void rd_ir_nonbonded(tpr_reader_t* r, md_tpr_nonbonded_t* nb) {
    const int v = r->version;

    rd_int(r);                                      // Integrator
    if (v >= tpxv_Pre96Version62) {
        rd_int64(r);                                // nsteps
        rd_int64(r);                                // init_step
    } else {
        rd_int(r);
        rd_int(r);
    }
    rd_int(r);                                      // simulation_part

    if (v >= tpxv_MTS) {
        // Multiple time stepping: the levels are only counted when it is used, and none are then written
        // otherwise (a reading t_inputrec starts without levels)
        const bool use_mts = rd_bool(r);
        const size_t num_levels = use_mts ? rd_count(r) : 0;
        if (num_levels > 16) {
            rd_fail(r);
            return;
        }
        skip_ints(r, 2 * num_levels);               // Force groups and step factor per level
    }
    if (v >= tpxv_MassRepartitioning) {
        rd_real(r);                                 // Mass repartition factor
    }
    if (v >= tpxv_EnsembleTemperature) {
        rd_int(r);                                  // Ensemble temperature setting
        rd_real(r);                                 // Ensemble temperature
    }
    if (v >= tpxv_Pre96Version67 && v < tpxv_OutputControlInKeyValueTree) {
        rd_int(r);                                  // nstcalcenergy
    }

    int32_t scheme = MD_TPR_CUTOFF_SCHEME_GROUP;
    if (v >= tpxv_Pre96Version81) {
        scheme = rd_int(r);
        if (v < tpxv_Pre96Version94) {
            // The order of the two was inverted
            scheme = (scheme == 0) ? MD_TPR_CUTOFF_SCHEME_GROUP : MD_TPR_CUTOFF_SCHEME_VERLET;
        }
    }
    rd_int(r);                                      // Once ns_type
    rd_int(r);                                      // nstlist
    rd_int(r);                                      // Once ndelta
    rd_real(r);                                     // rtpi
    rd_int(r);                                      // nstcomm
    rd_int(r);                                      // comm_mode
    if (v < tpxv_RemoveObsoleteParameters1) {
        rd_int(r);                                  // nstcheckpoint
    }
    rd_int(r);                                      // nstcgsteep
    rd_int(r);                                      // nbfgscorr
    if (v < tpxv_OutputControlInKeyValueTree) {
        skip_ints(r, 6);                            // nstlog, nstxout, nstvout, nstfout, nstenergy, nstxout_compressed
    }
    if (v >= tpxv_Pre96Version59) {
        md_xdr_skip(&r->xdr, 16);                   // init_t and delta_t, doubles
    } else {
        skip_reals(r, 2);
    }
    if (v < tpxv_OutputControlInKeyValueTree) {
        rd_real(r);                                 // x_compression_precision
    }
    if (v >= tpxv_Pre96Version81) {
        rd_real(r);                                 // verletbuf_tol
    }
    if (v >= tpxv_VerletBufferPressureTol) {
        rd_real(r);                                 // Verlet buffer pressure tolerance
    }
    const float rlist = (float)rd_real(r);
    if (v >= tpxv_Pre96Version67 && v < tpxv_RemoveTwinRange) {
        rd_real(r);                                 // rlistlong
    }
    if (v >= tpxv_Pre96Version82 && v != tpxv_Pre96Version90) {
        rd_int(r);                                  // nstcalclr
    }

    const int32_t coulomb_type = rd_int(r);
    int32_t coulomb_modifier;
    if (v >= tpxv_Pre96Version81) {
        coulomb_modifier = rd_int(r);
    } else {
        coulomb_modifier = (scheme == MD_TPR_CUTOFF_SCHEME_VERLET) ? MD_TPR_MODIFIER_POT_SHIFT : MD_TPR_MODIFIER_NONE;
    }
    const float rcoulomb_switch = (float)rd_real(r);
    const float rcoulomb = (float)rd_real(r);

    const int32_t vdw_type = rd_int(r);
    int32_t vdw_modifier;
    if (v >= tpxv_Pre96Version81) {
        vdw_modifier = rd_int(r);
    } else {
        vdw_modifier = (scheme == MD_TPR_CUTOFF_SCHEME_VERLET) ? MD_TPR_MODIFIER_POT_SHIFT : MD_TPR_MODIFIER_NONE;
    }
    const float rvdw_switch = (float)rd_real(r);
    const float rvdw = (float)rd_real(r);
    const int32_t disp_corr = rd_int(r);
    const float epsilon_r = (float)rd_real(r);
    const float epsilon_rf = (float)rd_real(r);
    rd_real(r);                                     // Table extension

    if (v < tpxv_RemoveImplicitSolvation) {
        rd_int(r);
        rd_int(r);
        rd_real(r);
        rd_real(r);
        rd_int(r);
        skip_reals(r, 4);
        if (v >= tpxv_Pre96Version60) {
            rd_real(r);
            rd_int(r);
        }
        rd_real(r);
    }
    if (v >= tpxv_Pre96Version81) {
        rd_real(r);                                 // Fourier spacing
    }
    skip_ints(r, 4);                                // nkx, nky, nkz, pme_order
    const float ewald_rtol = (float)rd_real(r);
    float ewald_rtol_lj = ewald_rtol;
    if (v >= tpxv_Pre96Version93) {
        ewald_rtol_lj = (float)rd_real(r);
    }
    rd_int(r);                                      // Ewald geometry
    rd_real(r);                                     // Surface dielectric constant
    if (v < tpxv_RemoveObsoleteParameters1) {
        rd_bool(r);                                 // bOptFFT
    }
    int32_t ljpme_comb_rule = 0;
    if (v >= tpxv_Pre96Version93) {
        ljpme_comb_rule = rd_int(r);
    }

    // Whatever went wrong, what is read is only kept when all of it makes sense
    const bool sane = rd_ok(r) &&
        coulomb_type >= 0 && coulomb_type <= MD_TPR_COULOMB_FMM &&
        vdw_type >= 0 && vdw_type <= MD_TPR_VDW_PME &&
        coulomb_modifier >= 0 && coulomb_modifier <= MD_TPR_MODIFIER_FORCE_SWITCH &&
        vdw_modifier >= 0 && vdw_modifier <= MD_TPR_MODIFIER_FORCE_SWITCH &&
        disp_corr >= 0 && disp_corr <= MD_TPR_DISP_CORR_ALL_ENER &&
        (scheme == MD_TPR_CUTOFF_SCHEME_VERLET || scheme == MD_TPR_CUTOFF_SCHEME_GROUP) &&
        rvdw >= 0.0f && rvdw < 1000.0f && rcoulomb >= 0.0f && rcoulomb < 1000.0f &&
        rvdw_switch >= 0.0f && rvdw_switch <= rvdw + 1.0e-6f && rcoulomb_switch >= 0.0f &&
        epsilon_r >= 0.0f && epsilon_rf >= 0.0f && (ljpme_comb_rule == 0 || ljpme_comb_rule == 1);
    if (!sane) {
        MD_LOG_INFO("TPR: The non-bonded settings of the simulation parameters could not be read (tpx version %d)", v);
        return;
    }

    *nb = (md_tpr_nonbonded_t){
        .valid = true,
        .cutoff_scheme = scheme,
        .rlist = rlist,
        .vdw_type = vdw_type,
        .vdw_modifier = vdw_modifier,
        .rvdw_switch = rvdw_switch,
        .rvdw = rvdw,
        .coulomb_type = coulomb_type,
        .coulomb_modifier = coulomb_modifier,
        .rcoulomb_switch = rcoulomb_switch,
        .rcoulomb = rcoulomb,
        .epsilon_r = epsilon_r,
        .epsilon_rf = epsilon_rf,
        .disp_corr = disp_corr,
        .ewald_rtol = ewald_rtol,
        .ewald_rtol_lj = ewald_rtol_lj,
        .ljpme_comb_rule = ljpme_comb_rule,
    };
}

// ### FILE ###

bool md_tpr_data_parse_buffer(md_tpr_data_t* data, const void* buffer, size_t size, md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(alloc);
    MEMSET(data, 0, sizeof(md_tpr_data_t));
    data->pbc = MD_TPR_PBC_UNSET;

    if (!buffer || size == 0) {
        MD_LOG_ERROR("TPR: Empty input");
        return false;
    }

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);
    bool result = false;

    // ## Header (always XDR)
    tpr_reader_t hdr = { .xdr = md_xdr_init(buffer, size) };

    const str_t version_string = rd_string(&hdr);
    if (!rd_ok(&hdr) || !str_begins_with(version_string, STR_LIT("VERSION"))) {
        MD_LOG_ERROR("TPR: Not a GROMACS run input file");
        goto done;
    }

    const int32_t precision = rd_int(&hdr);
    if (precision != 4 && precision != 8) {
        MD_LOG_ERROR("TPR: Unknown precision, reals are %d bytes", precision);
        goto done;
    }
    hdr.dbl = precision == 8;

    const int32_t version = rd_int(&hdr);
    hdr.version = version;
    if (version >= tpxv_Pre96Version77 && version <= tpxv_Pre96Version79) {
        rd_string(&hdr);                            // File tag, misplaced in these versions
    }
    const int32_t generation = rd_int(&hdr);
    if (version >= tpxv_Pre96Version81) {
        rd_string(&hdr);                            // File tag
    }
    if (!rd_ok(&hdr)) {
        MD_LOG_ERROR("TPR: Truncated header");
        goto done;
    }
    if (version <= tpxv_Pre96Version57) {
        MD_LOG_ERROR("TPR: File version %d is older than GROMACS 4.0 and not supported", version);
        goto done;
    }
    if (generation > TPX_GENERATION_MAX) {
        MD_LOG_ERROR("TPR: Topology generation %d of this file (version %d) is newer than the supported %d", generation, version, TPX_GENERATION_MAX);
        goto done;
    }

    const int32_t natoms = rd_int(&hdr);
    const int32_t ngtc   = rd_int(&hdr);
    if (version < tpxv_Pre96Version62) {
        rd_int(&hdr);
        rd_real(&hdr);
    }
    if (version >= tpxv_Pre96Version79) {
        rd_int(&hdr);                               // Free energy state
    }
    rd_real(&hdr);                                  // Lambda
    const bool has_ir  = rd_bool(&hdr);
    const bool has_top = rd_bool(&hdr);
    const bool has_x   = rd_bool(&hdr);
    const bool has_v   = rd_bool(&hdr);
    const bool has_f   = rd_bool(&hdr);
    const bool has_box = rd_bool(&hdr);

    if (!rd_ok(&hdr) || natoms < 0 || ngtc < 0) {
        MD_LOG_ERROR("TPR: Truncated or corrupt header");
        goto done;
    }
    if (!has_top) {
        MD_LOG_ERROR("TPR: The file has no topology");
        goto done;
    }

    data->file_version = version;
    data->file_generation = generation;
    data->double_precision = hdr.dbl;

    // ## Body
    tpr_reader_t body = hdr;
    if (version >= tpxv_AddSizeField && generation >= TPX_GENERATION_ADD_SIZE_FIELD) {
        const int64_t body_size = rd_int64(&hdr);
        if (body_size > 0 && (uint64_t)body_size * 4 == md_xdr_remaining(&hdr.xdr)) {
            // GROMACS 2020 betas wrote the body with every byte padded to 4. Nothing reads those.
            MD_LOG_ERROR("TPR: The file was written by a beta version of GROMACS 2020, which is not supported");
            goto done;
        }
        const uint8_t* body_ptr = (body_size >= 0 && (uint64_t)body_size <= md_xdr_remaining(&hdr.xdr)) ? md_xdr_take(&hdr.xdr, (size_t)body_size) : NULL;
        if (!body_ptr) {
            MD_LOG_ERROR("TPR: The file is truncated");
            goto done;
        }
        body.xdr = md_xdr_init(body_ptr, (size_t)body_size);
        body.in_memory = true;
    }

    if (has_box) {
        float box[9];
        rd_reals(&body, box, 9);
        MEMCPY(data->box, box, sizeof(box));
        data->has_box = true;
        if (version >= tpxv_Pre96Version51) skip_reals(&body, 9);  // Relative box
        skip_reals(&body, 9);                                       // Box velocity
        if (version < tpxv_Pre96Version56) skip_reals(&body, 9);
    }
    if (ngtc > 0) {
        if (version < tpxv_Pre96Version69) skip_reals(&body, (size_t)ngtc);
        skip_reals(&body, (size_t)ngtc);            // Legacy temperature coupling state
    }

    rd_mtop(&body, data, version_string, alloc, temp_alloc);
    if (!rd_ok(&body)) {
        MD_LOG_ERROR("TPR: Failed to read the topology, the file is truncated or corrupt");
        goto done;
    }
    if (data->num_atoms != (size_t)natoms) {
        MD_LOG_ERROR("TPR: The header says %d atoms, the topology %zu", natoms, data->num_atoms);
        goto done;
    }

    const size_t count = data->num_atoms * 3;
    if (has_x) {
        data->x = md_array_create(float, count, alloc);
        rd_reals(&body, data->x, count);
    }
    if (has_v) {
        data->v = md_array_create(float, count, alloc);
        rd_reals(&body, data->v, count);
    }
    if (has_f) {
        skip_reals(&body, count);
    }
    if (!rd_ok(&body)) {
        MD_LOG_ERROR("TPR: Failed to read the coordinates, the file is truncated");
        goto done;
    }

    // The simulation parameters start with the periodic boundary type. They are nice to have, not
    // needed, so a file that ends here or cannot be read further is still fine: what is read of them
    // is read from a copy of the cursor, and only kept when complete.
    if (has_ir && version >= tpxv_Pre96Version53) {
        tpr_reader_t ir = body;
        const int32_t pbc = rd_int(&ir);
        if (rd_ok(&ir) && pbc >= MD_TPR_PBC_XYZ && pbc <= MD_TPR_PBC_UNSET) {
            data->pbc = pbc;
        }
        rd_bool(&ir);                               // Periodic molecules
        if (version <= TPX_VERSION_IR_MAX) {
            rd_ir_nonbonded(&ir, &data->nonbonded);
        } else {
            MD_LOG_INFO("TPR: The simulation parameters of tpx version %d are newer than this reader knows (%d), the non-bonded settings are not read", version, TPX_VERSION_IR_MAX);
        }
    }

    result = true;

done:
    md_temp_end(temp);
    if (!result) {
        md_tpr_data_free(data, alloc);
    }
    return result;
}

bool md_tpr_data_parse_file(md_tpr_data_t* data, str_t filename, md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(alloc);

    md_file_t file = {0};
    if (!md_file_open(&file, filename, MD_FILE_READ)) {
        MD_LOG_ERROR("TPR: Could not open file '" STR_FMT "'", STR_ARG(filename));
        return false;
    }

    bool result = false;
    const int64_t size = md_file_size(file);
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    void* buffer = size > 0 ? md_alloc(temp_alloc, (size_t)size) : NULL;
    if (!buffer) {
        MD_LOG_ERROR("TPR: File '" STR_FMT "' is empty or too large", STR_ARG(filename));
    } else if (md_file_read(file, buffer, (size_t)size) != (size_t)size) {
        MD_LOG_ERROR("TPR: Failed to read file '" STR_FMT "'", STR_ARG(filename));
    } else {
        result = md_tpr_data_parse_buffer(data, buffer, (size_t)size, alloc);
    }

    md_temp_end(temp);
    md_file_close(&file);
    return result;
}

void md_tpr_data_free(md_tpr_data_t* data, md_allocator_i* alloc) {
    ASSERT(data);
    ASSERT(alloc);
    for (size_t i = 0; i < md_array_size(data->moltypes); ++i) {
        md_array_free(data->moltypes[i].atoms, alloc);
        md_array_free(data->moltypes[i].residues, alloc);
        md_array_free(data->moltypes[i].bonds, alloc);
        const md_tpr_moltype_t* mt = &data->moltypes[i];
        if (mt->excl_offset) {
            md_free(alloc, mt->excl, sizeof(uint32_t) * MAX(mt->excl_offset[mt->num_atoms], 1));
            md_free(alloc, mt->excl_offset, sizeof(uint32_t) * (mt->num_atoms + 1));
        }
    }
    md_array_free(data->moltypes, alloc);
    md_array_free(data->molblocks, alloc);
    md_array_free(data->intermolecular_bonds, alloc);
    md_array_free(data->lj, alloc);
    md_array_free(data->x, alloc);
    md_array_free(data->v, alloc);
    if (data->str_data) {
        md_free(alloc, data->str_data, data->str_size);
    }
    MEMSET(data, 0, sizeof(md_tpr_data_t));
    data->pbc = MD_TPR_PBC_UNSET;
}

bool md_tpr_atoms_excluded(const md_tpr_data_t* data, size_t atom_a, size_t atom_b) {
    if (!data || atom_a >= data->num_atoms || atom_b >= data->num_atoms) return false;
    if (atom_a == atom_b) return true;
    // The molecule of atom_a
    size_t offset = 0;
    for (size_t b = 0; b < data->num_molblocks; ++b) {
        const md_tpr_molblock_t* mb = &data->molblocks[b];
        const md_tpr_moltype_t* mt = &data->moltypes[mb->moltype];
        const size_t block_atoms = (size_t)mb->nmol * mt->num_atoms;
        if (atom_a < offset + block_atoms) {
            if (!mt->excl_offset || mt->num_atoms == 0) return false;
            const size_t mol_beg = offset + (atom_a - offset) / mt->num_atoms * mt->num_atoms;
            if (atom_b < mol_beg || atom_b >= mol_beg + mt->num_atoms) return false;
            const uint32_t la = (uint32_t)(atom_a - mol_beg);
            const uint32_t lb = (uint32_t)(atom_b - mol_beg);
            uint32_t lo = mt->excl_offset[la];
            uint32_t hi = mt->excl_offset[la + 1];
            while (lo < hi) {
                const uint32_t mid = (lo + hi) / 2;
                if (mt->excl[mid] == lb) return true;
                if (mt->excl[mid] < lb) lo = mid + 1; else hi = mid;
            }
            return false;
        }
        offset += block_atoms;
    }
    return false;
}

// ### SYSTEM ###

float md_tpr_lj_vdw_radius(md_tpr_lj_t lj) {
    if (!(lj.c6 > 0.0f) || !(lj.c12 > 0.0f)) {
        return 0.0f;
    }
    // sigma = (c12/c6)^(1/6), the potential minimum is at r_min = 2^(1/6) sigma, and the van der
    // Waals radius is half of that: the distance two identical particles keep to each other
    const double sigma = pow((double)lj.c12 / (double)lj.c6, 1.0 / 6.0);
    return (float)(0.5 * pow(2.0, 1.0 / 6.0) * sigma * 10.0);   // nm -> Ångström
}

bool md_tpr_system_init_from_data(md_system_t* sys, md_system_state_t* state, const md_tpr_data_t* data) {
    ASSERT(sys);
    ASSERT(state);
    ASSERT(data);

    if (!sys->alloc) {
        MD_LOG_ERROR("System allocator not set");
        return false;
    }
    if (!state->alloc) {
        MD_LOG_ERROR("State allocator not set");
        return false;
    }

    md_system_reset(sys);
    const size_t num_atoms = data->num_atoms;
    md_system_state_init(state, num_atoms);
    if (num_atoms == 0) {
        MD_LOG_ERROR("TPR: The system has no atoms");
        return false;
    }

    md_allocator_i* alloc = sys->alloc;
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);

    // ## Residue numbering, as gmx does it when writing out the whole system: molecules of at most
    // max_renum residues are numbered consecutively after the highest number any other molecule uses
    // (see gmx_mtop_t::finalize). A single molecule keeps its own numbers.
    const int32_t max_renum = (data->num_molblocks == 1 && data->molblocks[0].nmol == 1) ? 0 : 1;
    int32_t next_renum = 0;
    for (size_t t = 0; t < data->num_moltypes; ++t) {
        const md_tpr_moltype_t* mt = &data->moltypes[t];
        if ((int64_t)mt->num_residues > max_renum) {
            for (size_t r = 0; r < mt->num_residues; ++r) {
                next_renum = MAX(next_renum, mt->residues[r].nr);
            }
        }
    }
    next_renum += 1;

    // ## Atom types, resolved once per molecule type atom
    // A type is everything the topology says about a particle that is not per atom: its name, element,
    // force field type, mass and particle type. Atoms sharing all of those share a type, so the type
    // carries the mass exactly and no per atom copy of it is needed.
    //   - Element and radius from the atomic number when the topology has one.
    //   - Otherwise the particle is a coarse grained bead (or a massless virtual site): no element, and
    //     a radius from the Lennard-Jones parameters of its non-bonded type when it has any.
    // The predefined bead tables then add what they know about the beads (backbone, side chain...).
    const size_t capacity = ROUND_UP(num_atoms, 16);
    md_array_resize(sys->atom.type_idx, capacity, alloc);
    md_array_resize(sys->atom.flags, capacity, alloc);
    MEMSET(sys->atom.type_idx, 0, capacity * sizeof(md_atom_type_idx_t));
    MEMSET(sys->atom.flags, 0, capacity * sizeof(md_flags_t));

    sys->atom.type.count = 0;
    md_atom_type_find_or_add(&sys->atom.type, STR_LIT("Unk"), 0, 0.0f, 0.0f, 0, 0, alloc);

    md_hashmap32_t type_map = { .allocator = md_temp_allocator(temp) };
    md_atom_type_idx_t** moltype_type = md_temp_alloc_array(temp, md_atom_type_idx_t*, data->num_moltypes + 1);
    for (size_t t = 0; t < data->num_moltypes; ++t) {
        const md_tpr_moltype_t* mt = &data->moltypes[t];
        moltype_type[t] = md_temp_alloc_array(temp, md_atom_type_idx_t, mt->num_atoms + 1);
        for (size_t i = 0; i < mt->num_atoms; ++i) {
            const md_tpr_atom_t* atom = &mt->atoms[i];
            const md_atomic_number_t z = (atom->atomic_number > 0 && atom->atomic_number < 119) ? (md_atomic_number_t)atom->atomic_number : 0;

            const struct {
                uint32_t nb_type;
                float    mass;
                uint32_t z;
                uint32_t ptype;
            } key_data = { atom->type_idx, atom->mass, z, atom->ptype };
            const uint64_t key = md_hash64_str(atom->name, md_hash64_str(atom->type, md_hash64(&key_data, sizeof(key_data), 0)));

            const uint32_t* cached = md_hashmap_get(&type_map, key);
            if (cached) {
                moltype_type[t][i] = (md_atom_type_idx_t)*cached;
                continue;
            }

            float radius = 0.0f;
            uint32_t color = 0;
            md_flags_t flags = 0;
            if (z) {
                radius = md_atomic_number_vdw_radius(z);
                color  = md_atomic_number_cpk_color(z);
            } else {
                const float lj_radius = md_tpr_lj_vdw_radius(md_tpr_lj_pair(data, atom->type_idx, atom->type_idx));
                radius = lj_radius > 0.0f ? lj_radius : md_atomic_number_vdw_radius(0);
                // A virtual site without Lennard-Jones is a charge site (TIP4P's M), not a bead
                if (atom->ptype != MD_TPR_PTYPE_VSITE || lj_radius > 0.0f) {
                    flags |= MD_FLAG_COARSE_GRAINED;
                }
            }

            const md_atom_type_idx_t type = md_atom_type_add(&sys->atom.type, atom->name, atom->type, z, atom->mass, radius, color, flags, alloc);

            md_hashmap_add(&type_map, key, (uint32_t)type);
            moltype_type[t][i] = type;
        }
    }

    // ## Atoms and residues
    float* charge = md_temp_alloc_array(temp, float, num_atoms);

    size_t ai = 0;
    for (size_t b = 0; b < data->num_molblocks; ++b) {
        const md_tpr_molblock_t* mb = &data->molblocks[b];
        const md_tpr_moltype_t* mt = &data->moltypes[mb->moltype];
        const bool renumber = (int64_t)mt->num_residues <= max_renum;

        for (int32_t m = 0; m < mb->nmol; ++m) {
            int32_t prev_res = -1;
            for (size_t i = 0; i < mt->num_atoms; ++i, ++ai) {
                const md_tpr_atom_t* atom = &mt->atoms[i];
                if (atom->residue != prev_res) {
                    const md_tpr_residue_t* res = &mt->residues[atom->residue];
                    const md_sequence_id_t seq_id = renumber ? next_renum + m * (int32_t)mt->num_residues + atom->residue : res->nr;
                    md_array_push(sys->component.atom_offset, (uint32_t)ai, alloc);
                    md_array_push(sys->component.name, make_label(res->name), alloc);
                    md_array_push(sys->component.seq_id, seq_id, alloc);
                    md_array_push(sys->component.flags, 0, alloc);
                    sys->component.count += 1;
                    prev_res = atom->residue;
                }

                const md_atom_type_idx_t type = moltype_type[mb->moltype][i];
                sys->atom.type_idx[ai] = type;
                sys->atom.flags[ai] = sys->atom.type.flags[type];
                charge[ai] = atom->charge;
            }
        }
        if (renumber) {
            next_renum += mb->nmol * (int32_t)mt->num_residues;
        }
    }
    ASSERT(ai == num_atoms);
    md_array_push(sys->component.atom_offset, (uint32_t)num_atoms, alloc);  // Final sentinel
    sys->atom.count = num_atoms;

    md_util_system_augment_atom_types(sys);

    // ## Coordinates, nm -> Ångström
    if (data->x) {
        for (size_t i = 0; i < num_atoms; ++i) {
            state->xyz[i] = vec3_set(data->x[i * 3 + 0] * 10.0f, data->x[i * 3 + 1] * 10.0f, data->x[i * 3 + 2] * 10.0f);
        }
    } else {
        MD_LOG_INFO("TPR: The file has no coordinates");
    }

    if (data->has_box && data->pbc != MD_TPR_PBC_NO) {
        float box[3][3];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                box[i][j] = data->box[i][j] * 10.0f;
            }
        }
        state->unitcell = md_unitcell_from_matrix_float(MD_AS_CONST_MAT3(box));
    }
    state->num_atoms = num_atoms;

    // ## Bonds, laid out per molecule
    {
        size_t count = data->num_intermolecular_bonds;
        for (size_t b = 0; b < data->num_molblocks; ++b) {
            count += (size_t)data->molblocks[b].nmol * data->moltypes[data->molblocks[b].moltype].num_bonds;
        }
        md_array_resize(sys->bond.pairs, count, alloc);
        md_array_resize(sys->bond.flags, count, alloc);

        const md_bond_flags_t flags = MD_BOND_FLAG_COVALENT | MD_BOND_FLAG_TOPOLOGY;
        size_t bi = 0;
        md_atom_idx_t offset = 0;
        for (size_t b = 0; b < data->num_molblocks; ++b) {
            const md_tpr_moltype_t* mt = &data->moltypes[data->molblocks[b].moltype];
            for (int32_t m = 0; m < data->molblocks[b].nmol; ++m) {
                for (size_t k = 0; k < mt->num_bonds; ++k, ++bi) {
                    sys->bond.pairs[bi].idx[0] = offset + mt->bonds[k].idx[0];
                    sys->bond.pairs[bi].idx[1] = offset + mt->bonds[k].idx[1];
                    sys->bond.flags[bi] = flags;
                }
                offset += (md_atom_idx_t)mt->num_atoms;
            }
        }
        for (size_t k = 0; k < data->num_intermolecular_bonds; ++k, ++bi) {
            sys->bond.pairs[bi] = data->intermolecular_bonds[k];
            sys->bond.flags[bi] = flags;
        }
        ASSERT(bi == count);
        sys->bond.count = count;
        md_bond_build_connectivity(&sys->bond, num_atoms, alloc);
    }

    md_util_system_infer_comp_flags(sys);

    // ## Attributes
    md_attributes_publish_atom_column(&sys->attributes, STR_LIT("atom/charge"), md_unit_elementary_charge(), 1, charge, num_atoms);
    if (data->v) {
        md_attributes_publish_atom_column(&sys->attributes, STR_LIT("atom/velocity"), md_unit_div(md_unit_nanometer(), md_unit_picosecond()), 3, data->v, num_atoms);
    }

    if (!str_empty(data->name)) {
        sys->description = str_copy(data->name, alloc);
    }

    // ## Non-bonded force field, when its interactions can be evaluated pair by pair
    {
        md_nb_forcefield_t* ff = md_alloc(alloc, sizeof(md_nb_forcefield_t));
        if (md_nb_forcefield_init_from_tpr(ff, data, alloc)) {
            sys->nonbonded = ff;
        } else {
            md_free(alloc, ff, sizeof(md_nb_forcefield_t));
        }
    }

    md_temp_end(temp);
    return true;
}

bool md_tpr_system_init_from_file(md_system_t* sys, md_system_state_t* state, str_t filename) {
    ASSERT(sys);
    md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
    md_allocator_i* temp_alloc = md_temp_allocator(temp);

    md_tpr_data_t data = {0};
    const bool result = md_tpr_data_parse_file(&data, filename, temp_alloc) && md_tpr_system_init_from_data(sys, state, &data);

    md_temp_end(temp);
    return result;
}
