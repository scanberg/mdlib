#include <md_infer.h>
#include <md_element.h>
#include "md_infer_internal.h"
#include "md_infer_profile.h"
#include <core/md_executor.h>
#include <core/md_arena_allocator.h>
#include <core/md_log.h>
#include <core/md_os.h>
#include <md_chem.h>
#include <md_system_fields.h>

// The inference of a system: md_system_infer runs the steps of md_bond.c, md_ring.c, md_atom.c, md_hierarchy.c,
// md_backbone.c, md_lipid.c, md_fragment.c and md_chem.c in order and gives the system what they infer.

// ### MOLECULES, RINGS AND BACKBONES, THEN CHEMISTRY ###
//
// Three branches of the inference which read the system and write only results of their own: the molecules, the rings,
// the backbones. None reads what another makes, so they run side by side, and their results are given to the system
// after, in the order of the steps. The chemistry follows them: it is perceived for the molecules which are not copies
// of another (md_molecule_data_t, STRUCTURES), and reads the rings.
typedef struct infer_branches_t {
    const md_system_t*       sys;
    bool molecules, rings, backbone;
    md_allocator_i*     alloc[3];           // Of the results of each branch
    md_molecule_data_t  molecule;
    md_ring_data_t      ring;
    md_backbone_data_t  backbone_data;
    double              ms[3];              // Of each branch, when profiled (mdi_infer_profile)
} infer_branches_t;

// The longest branch first: the molecules, then the rings, then the backbones
static void infer_branch(void* arg, uint32_t beg, uint32_t end, uint32_t thread_idx) {
    (void)thread_idx;
    infer_branches_t* br = (infer_branches_t*)arg;
    for (uint32_t b = beg; b < end; ++b) {
        const md_tick_t t0 = mdi_infer_profile ? md_tick_now() : 0;
        if (b == 0) {
            if (br->molecules) md_molecule_infer(&br->molecule, br->sys, br->alloc[0]);
        } else if (b == 1) {
            if (br->rings) md_ring_infer(&br->ring, br->sys, br->alloc[1]);
        } else {
            if (br->backbone) md_backbone_infer(&br->backbone_data, br->sys, br->alloc[2]);
        }
        if (mdi_infer_profile) br->ms[b] = md_tick_to_milliseconds(md_tick_now() - t0);
    }
}

// An array, anew in alloc
#define COPY_INTO(dst, src, alloc) do { (dst) = NULL; if (md_array_size(src)) md_array_push_array((dst), (src), md_array_size(src), (alloc)); } while (0)

static void molecule_data_copy(md_molecule_data_t* dst, const md_molecule_data_t* src, md_allocator_i* alloc) {
    *dst = (md_molecule_data_t){ .count = src->count };
    COPY_INTO(dst->offset,     src->offset,     alloc);
    COPY_INTO(dst->atom_idx,   src->atom_idx,   alloc);
    COPY_INTO(dst->parent_idx, src->parent_idx, alloc);
    COPY_INTO(dst->atom_slot,  src->atom_slot,  alloc);
    COPY_INTO(dst->type_idx,   src->type_idx,   alloc);
    for (size_t t = 0; t < src->type.count; ++t) md_molecule_type_data_push(&dst->type, src->type.name[t], src->type.flags[t], alloc);
    COPY_INTO(dst->structure_idx, src->structure_idx, alloc);
    COPY_INTO(dst->structure_ref, src->structure_ref, alloc);
    dst->structure_count  = src->structure_count;
    dst->structure_hash   = src->structure_hash;
    dst->topology_version = src->topology_version;
}

// ### PROFILING ###

mdi_infer_profile_t* mdi_infer_profile = NULL;

static void profile_record(const char* name, double ms) {
    mdi_infer_profile_t* p = mdi_infer_profile;
    if (p && p->count < MDI_INFER_PROFILE_MAX_STEPS) {
        p->name[p->count] = name;
        p->ms[p->count] = ms;
        p->count += 1;
    }
}

void mdi_infer_profile_step(int64_t* t, const char* name) {
    if (!mdi_infer_profile) return;
    const md_tick_t now = md_tick_now();
    profile_record(name, md_tick_to_milliseconds(now - *t));
    *t = now;
}

// The step which ended now, since the one before it
#define PROFILE_STEP(name) do { if (mdi_infer_profile) { const md_tick_t t_ = md_tick_now(); profile_record(name, md_tick_to_milliseconds(t_ - prof_t)); prof_t = t_; } } while (0)

// Try to fill in missing fields for molecule struct
// (Coordinates & Elements) -> Covalent Bonds
// (residues & Bonds)       -> Chains
// (Chains)                 -> Backbone
bool md_system_infer(md_system_t* sys, const md_system_state_t* state, md_infer_flags_t flags, const md_executor_t* executor) {
    MD_ASSERT(sys);

    if (sys->atom.count == 0) return false;
    if (!sys->alloc) {
        MD_LOG_ERROR("System allocator not set");
        return false;
    }
    md_tick_t prof_t = mdi_infer_profile ? md_tick_now() : 0;
    md_system_topology_changed(sys);

    // What is written is of the derived (md_system_t, ALLOCATORS): everything below but the colors of the types
    md_allocator_i* alloc = md_system_derived_alloc(sys);

    const bool cg = md_system_is_coarse_grained(sys);
    size_t num_particle_types = md_system_particle_type_count(sys);

    if (!sys->atom.flags) {
        md_array_resize(sys->atom.flags, sys->atom.count, alloc);
        md_memset(sys->atom.flags, 0, md_array_bytes(sys->atom.flags));
    }

    if (flags & MD_SYSTEM_INFER_COLOR_BIT) {
        if (sys->atom.type.color) {
            // @Note: The 0'th type is a sentinel type for representing invalid or unset atoms.
            // Hence it's skipped here
            for (size_t i = 1; i < num_particle_types; ++i) {
                uint32_t color = sys->atom.type.color[i];
                md_atomic_number_t z = sys->atom.type.z[i];
                if (color == 0) {
                    // Try to assign a color based on element if available
                    if (z != 0) {
                        color = md_atomic_number_cpk_color(z);
                        sys->atom.type.color[i] = color;
                    } else {
                        // Assign a random color for the type
                        float hue = (float)(i * 97 % 360) / 360.0f;
                        float chroma = 0.8f;
                        float light  = 1.0f;
                        float alpha  = 1.0f;
                        md_float3_t rgb = hcl_to_rgb(hue, chroma, light);
                        // Pack into 0xAARRGGBB
                        color = (uint32_t)(uint8_t)(alpha * 255) << 24 | (uint32_t)(uint8_t)(rgb.x * 255) << 16 | (uint32_t)(uint8_t)(rgb.y * 255) << 8 | (uint32_t)(uint8_t)(rgb.z * 255);
                        sys->atom.type.color[i] = color;
                    }
                }
            }
        }
    }
   
    PROFILE_STEP("setup");

    if (flags & MD_SYSTEM_INFER_BOND_BIT) {
        // Setting the bit means "derive this", unconditionally. A format which supplies its own
        // bonds (LAMMPS, PSF) is respected by the caller not setting the bit, not by a guard here.
        // Beads too: a coarse grained system gets the bonds of its model's topologies (mdi_infer_coarse_grained_bonds).
        // The bonds inferred, in place of those inferred before: those given (by a file, a topology, a user) stay. Without
        // any (most structure files) the inferred are the bonds, and nothing is merged.
        bool given = sys->bond.noncovalent_count > 0;
        for (size_t i = 0; i < sys->bond.count && !given; ++i) given = md_bond_origin(sys->bond.flags[i]) != MD_BOND_ORIGIN_INFERRED;
        md_bond_data_t bonds = {0};
        if (given) {
            md_temp_scope_t temp = md_temp_begin_avoid(alloc);
            md_bond_data_t inferred = {0};
            md_bond_infer(&inferred, sys, state, md_temp_allocator(temp), executor);
            md_bond_merge(&bonds, &sys->bond, &inferred, sys->atom.count, alloc);
            md_temp_end(temp);
        } else {
            md_bond_infer(&bonds, sys, state, alloc, executor);
        }
        md_bond_data_free(&sys->bond, alloc);
        sys->bond = bonds;
        md_system_topology_changed(sys);
        PROFILE_STEP("bonds");
    }

    // Whatever the bonds came from (inferred, a topology, a file, by hand): those between a metal and a non-metal are
    // coordination. Cheap, and every reader of the flags relies on it, so not behind a bit.
    md_system_infer_coordination(sys);
    PROFILE_STEP("coordination");

    // What the residues are (amino acids, nucleotides, waters, ions) and the roles of their atoms: by their names and by
    // the bonds, so after them. What a loader gives is the atoms and their names; everything below reads these.
    if (flags & MD_SYSTEM_INFER_RESIDUE_BIT) {
        md_system_infer_residue_flags(sys, state);
        PROFILE_STEP("residue flags");
    }

    // The parts before the molecules: a molecule holds its polymer strands whole (synthesize_molecule_links). They are
    // found from the instances the data gives, the bonds and the classified residues, which the molecules do not change.
    if (flags & MD_SYSTEM_INFER_PART_BIT) {
        md_chain_data_t chains = {0};
        md_chain_infer(&chains, sys, alloc);
        md_chain_data_free(&sys->chain, alloc);
        sys->chain = chains;
        md_part_data_t parts = {0};
        md_part_infer(&parts, sys, alloc);
        md_part_data_free(&sys->part, alloc);
        sys->part = parts;
        PROFILE_STEP("chains, parts");
    }

    bool molecules_inferred = false;
    // The molecules, the rings, the backbones: steps which read the system and make results of their own (MOLECULES,
    // RINGS AND BACKBONES), then the chemistry of the molecules and the rings, all given to the system in the order of
    // the steps
    {
        // Coarse grained systems get a hierarchy even without bonds (synthesize_hierarchy_links), and a residue is a
        // molecule whatever its bonds (synthesize_molecule_links): ions and waters given without bonds are molecules.
        // The chemistry finds the rings it needs where there are none.
        infer_branches_t br = {
            .sys       = sys,
            .molecules = (flags & MD_SYSTEM_INFER_MOLECULE_BIT) && (sys->bond.count || cg || sys->residue.count),
            .rings     = sys->bond.count && ((flags & MD_SYSTEM_INFER_MOLECULE_BIT) || ((flags & MD_SYSTEM_INFER_CHEMISTRY_BIT) && md_ring_data_count(&sys->ring) == 0)),
            .backbone  = (flags & MD_SYSTEM_INFER_BACKBONE_BIT) != 0,
        };
        const bool chemistry = (flags & MD_SYSTEM_INFER_CHEMISTRY_BIT) && sys->bond.count;
        if ((br.rings || br.molecules || chemistry) && sys->bond.count && !sys->bond.conn.offset) md_bond_build_connectivity(&sys->bond, sys->atom.count, alloc);
        PROFILE_STEP("connectivity");

        // Side by side on the executor the system's allocator, which is not for several threads, is the molecules' alone:
        // they are made in it, on whichever thread runs them, and the rings and the backbones each in an arena of its
        // own, copied into it after (they are small; the molecules are several arrays of the atoms, 19 ms to copy at 2M
        // atoms). But for an allocator this thread's temporary scopes may use meanwhile (md_temp_is_thread_arena),
        // where the molecules take an arena of their own as well. Otherwise one after the other, in the system's allocator.
        const bool side_by_side = md_executor_num_threads(executor) > 1;
        const bool molecules_own = side_by_side && md_temp_is_thread_arena(alloc);
        md_temp_scope_t temp = md_temp_begin_avoid(alloc);
        for (int i = 0; i < 3; ++i) {
            const bool own = side_by_side && (i > 0 || molecules_own);
            br.alloc[i] = own ? md_arena_allocator_create(md_get_heap_allocator(), MD_MEGABYTES(1)) : alloc;
        }
        md_parallel_for(side_by_side ? executor : NULL, 3, 1, infer_branch, &br);
        if (mdi_infer_profile) {
            if (br.molecules) profile_record("  molecules", br.ms[0]);
            if (br.rings)     profile_record("  rings", br.ms[1]);
            if (br.backbone)  profile_record("  backbones", br.ms[2]);
        }
        PROFILE_STEP("branches");

        // The chemistry of the molecules just inferred (the copies of one take its chemistry), in a temporary scope
        // until it is written. Protonation at pH 7 applies to residues without hydrogens only.
        mdi_chem_t* chem = NULL;
        if (chemistry) chem = mdi_chem_compute(sys, br.rings ? &br.ring : &sys->ring, state, MD_CHEM_FLAG_PROTONATE_PH7, br.molecules ? &br.molecule : NULL, md_temp_allocator(temp));
        if (chem) PROFILE_STEP("chemistry");

        if (br.molecules) {
            molecules_inferred = true;
            md_molecule_data_free(&sys->molecule, alloc);
            if (molecules_own) molecule_data_copy(&sys->molecule, &br.molecule, alloc);
            else              sys->molecule = br.molecule;
            // The fragments are of the molecules, looked up by their slots: gone with them until inferred again
            // (md_system_infer_fragments). Their types stay, for the colours to be kept.
            sys->fragment.count = 0;
            sys->fragment.template_count = 0;
            md_array_shrink(sys->fragment.mol_template, 0);
            md_array_shrink(sys->fragment.frag_offset, 0);
            md_array_shrink(sys->fragment.template_molecule, 0);
            md_array_shrink(sys->fragment.template_slot_offset, 0);
            md_array_shrink(sys->fragment.slot_frag, 0);
            md_array_shrink(sys->fragment.template_frag_offset, 0);
            md_array_shrink(sys->fragment.frag_type, 0);
            md_array_shrink(sys->fragment.frag_slot_offset, 0);
            md_array_shrink(sys->fragment.frag_slot, 0);
        }
        if (br.rings) {
            md_ring_data_free(&sys->ring, alloc);
            sys->ring = (md_ring_data_t){0};
            if (side_by_side) COPY_INTO(sys->ring.rings, br.ring.rings, alloc);
            else              sys->ring = br.ring;
        }
        // After the rings, before the hydrogen bond roles which read the orders
        PROFILE_STEP("write");
        if (chem) mdi_chem_write(sys, chem);
        if (chem) PROFILE_STEP("chemistry write");
        if (br.backbone) {
            md_backbone_data_t old = { .protein = sys->protein_backbone, .nucleic = sys->nucleic_backbone };
            md_backbone_data_free(&old, alloc);
            if (side_by_side) {
                sys->protein_backbone = (md_protein_backbone_data_t){ .range.count = br.backbone_data.protein.range.count, .segment.count = br.backbone_data.protein.segment.count };
                sys->nucleic_backbone = (md_nucleic_backbone_data_t){ .range.count = br.backbone_data.nucleic.range.count, .segment.count = br.backbone_data.nucleic.segment.count };
                COPY_INTO(sys->protein_backbone.range.offset,    br.backbone_data.protein.range.offset,    alloc);
                COPY_INTO(sys->protein_backbone.range.part_idx,  br.backbone_data.protein.range.part_idx,  alloc);
                COPY_INTO(sys->protein_backbone.segment.atoms,   br.backbone_data.protein.segment.atoms,   alloc);
                COPY_INTO(sys->protein_backbone.segment.res_idx, br.backbone_data.protein.segment.res_idx, alloc);
                COPY_INTO(sys->nucleic_backbone.range.offset,    br.backbone_data.nucleic.range.offset,    alloc);
                COPY_INTO(sys->nucleic_backbone.range.part_idx,  br.backbone_data.nucleic.range.part_idx,  alloc);
                COPY_INTO(sys->nucleic_backbone.segment.atoms,   br.backbone_data.nucleic.segment.atoms,   alloc);
                COPY_INTO(sys->nucleic_backbone.segment.res_idx, br.backbone_data.nucleic.segment.res_idx, alloc);
            } else {
                sys->protein_backbone = br.backbone_data.protein;
                sys->nucleic_backbone = br.backbone_data.nucleic;
            }
            md_system_set_chain_breaks(sys, br.backbone_data.breaks, br.backbone_data.num_breaks);
            if (!side_by_side) md_array_free(br.backbone_data.breaks, alloc);
        }
        for (int i = 0; i < 3; ++i) if (br.alloc[i] != alloc) md_arena_allocator_destroy(br.alloc[i]);
        md_temp_end(temp);
        PROFILE_STEP("write");
    }

    if (flags & MD_SYSTEM_INFER_LIPID_BIT) {
        // After the bonds, the rings (a sterol is told by its fused rings) and the classification of the residues
        md_lipid_data_t lipids = {0};
        md_lipid_infer(&lipids, sys, alloc);
        md_lipid_data_free(&sys->lipid, alloc);
        sys->lipid = lipids;
        md_system_infer_lipid_kinds(sys);
        PROFILE_STEP("lipids");
    }

    // Last of what depends on the classification of the residues, which the lipids are the last to change. Cheap,
    // and what a residue type is follows from the residues alone, so not behind a bit.
    md_system_infer_residue_types(sys);
    PROFILE_STEP("residue types");

    if (flags & MD_SYSTEM_INFER_FRAGMENT_BIT) {
        // Last: made of the molecules, the rings, the chemistry, the roles of the atoms and the lipid tails
        md_system_infer_fragments(sys);
        PROFILE_STEP("fragments");
    }

    // The structures of the molecules inferred are of the topology as it is now (md_system_molecule_structures_valid)
    if (molecules_inferred) sys->molecule.topology_version = sys->topology_version;

    // Record the state the topology above was inferred from. Written here, as part of performing
    // the inference, so it cannot disagree with what it produced.
    if (state != &sys->reference) {
        sys->reference.alloc = alloc;
        md_system_state_copy(&sys->reference, state);
    }
    PROFILE_STEP("reference");

    return true;
}

// ### INFERENCE WHILE THE SYSTEM IS IN USE ###

#define INFER_COPY_ARRAY(f) do { \
        work->sys.f = NULL; \
        if (md_array_size(sys->f)) md_array_push_array(work->sys.f, sys->f, md_array_size(sys->f), alloc); \
    } while (0);

#define INFER_COPY_STRINGS(f) do { \
        work->sys.f = NULL; \
        const size_t n = md_array_size(sys->f); \
        if (n) { \
            md_array_resize(work->sys.f, n, dalloc); \
            for (size_t i = 0; i < n; ++i) work->sys.f[i] = md_str_empty(sys->f[i]) ? (md_str_t){0} : md_alloc_str(dalloc, sys->f[i]); \
        } \
    } while (0);

#define INFER_COPY_DERIVED(f) do { \
        work->sys.f = NULL; \
        if (md_array_size(sys->f)) md_array_push_array(work->sys.f, sys->f, md_array_size(sys->f), dalloc); \
    } while (0);

bool md_system_infer_prepare(md_system_infer_work_t* work, const md_system_t* sys, const md_system_state_t* state, md_allocator_i* alloc) {
    MD_ASSERT(work);
    MD_ASSERT(sys);
    MD_ASSERT(alloc);
    md_memset(work, 0, sizeof(*work));

    // The core and the counts, shared. Not what is the system's alone: its description, its force field, its attributes
    // (which its owner may add to meanwhile).
    work->sys = *sys;
    work->sys.alloc = alloc;
    // What it infers, with the system's allocator of the derived, which the system frees it with once published: alloc
    // where the system has none (md_system_t, ALLOCATORS), and then alloc is kept with the system
    work->sys.derived_alloc = sys->derived_alloc;
    md_allocator_i* dalloc = md_system_derived_alloc(&work->sys);
    work->sys.nonbonded = NULL;
    work->sys.description = (md_str_t){0};
    md_memset(&work->sys.attributes, 0, sizeof(work->sys.attributes));
    work->sys.attributes.alloc = alloc;
    md_memset(&work->sys.reference, 0, sizeof(work->sys.reference));
    work->base_version = sys->topology_version;

    // The types of the atoms, which the owner may change meanwhile (a color, a radius): a few of them
    INFER_COPY_ARRAY(atom.type.name)
    INFER_COPY_ARRAY(atom.type.ff_type)
    INFER_COPY_ARRAY(atom.type.z)
    INFER_COPY_ARRAY(atom.type.mass)
    INFER_COPY_ARRAY(atom.type.radius)
    INFER_COPY_ARRAY(atom.type.color)
    INFER_COPY_ARRAY(atom.type.flags)

    // What inference writes, of its own: what the data gave of it is added to or replaced
    MD_SYSTEM_DERIVED_ARRAYS(INFER_COPY_DERIVED)
    MD_SYSTEM_DERIVED_STRINGS(INFER_COPY_STRINGS)

    // The one attribute inference reads: the formal charges the data gives, which the chemistry starts from
    // (md_chem_perceive)
    const md_attribute_t* q = md_attributes_find(&sys->attributes, MD_STR_LIT("atom/formal_charge"));
    if (q && md_attribute_slice_count(q, md_attribute_slice_all()) == sys->atom.count) {
        const float* values = (const float*)md_attribute_view(q, MD_ATTRIBUTE_TYPE_F32, 1, 1);
        if (values) {
            md_attributes_publish_atom_column(&work->sys.attributes, MD_STR_LIT("atom/formal_charge"), md_unit_none(), 1, values, sys->atom.count);
        }
    }

    // The coordinates, which the system's owner changes as it likes (a frame of a run): those inferred from, and kept
    // as the reference of what is inferred
    if (state && (state->xyz || state->num_atoms)) {
        work->sys.reference.alloc = dalloc;
        if (!md_system_state_copy(&work->sys.reference, state)) {
            md_system_infer_discard(work);
            return false;
        }
    }
    return true;
}

#undef INFER_COPY_ARRAY
#undef INFER_COPY_STRINGS
#undef INFER_COPY_DERIVED

bool md_system_infer_run(md_system_infer_work_t* work, md_infer_flags_t flags, const md_executor_t* executor) {
    MD_ASSERT(work);
    // The colors are of the types of the atoms, which stay the system's: given them before, which is what they are
    // shown with meanwhile
    flags &= ~MD_SYSTEM_INFER_COLOR_BIT;
    return md_system_infer(&work->sys, &work->sys.reference, flags, executor);
}

bool md_system_infer_publish(md_system_t* sys, md_system_infer_work_t* work) {
    MD_ASSERT(sys);
    MD_ASSERT(work);
    if (sys->topology_version != work->base_version) return false;
    // The core it shares is the system's still
    MD_ASSERT(sys->atom.count == work->sys.atom.count && sys->atom.type_idx == work->sys.atom.type_idx);
    MD_ASSERT(sys->residue.count == work->sys.residue.count && sys->residue.atom_offset == work->sys.residue.atom_offset);

    // Inferred with the system's allocator of the derived, which the system frees it with
    MD_ASSERT(work->sys.derived_alloc == sys->derived_alloc);
    // What the system had of the derived: given back (to an arena, which frees nothing, this is nothing)
    md_system_derived_free(sys);

#define MOVE(f) sys->f = work->sys.f;
    MD_SYSTEM_DERIVED_STRUCTS(MOVE)
    MD_SYSTEM_DERIVED_ARRAYS(MOVE)
    MD_SYSTEM_DERIVED_STRINGS(MOVE)
#undef MOVE
    sys->reference = work->sys.reference;

    // Of its own: the work's attributes go with it
    md_memset(work, 0, sizeof(*work));
    md_system_topology_changed(sys);
    return true;
}

void md_system_infer_discard(md_system_infer_work_t* work) {
    MD_ASSERT(work);
    // Its derived arrays and its reference, its own (the core is the system's)
    md_system_derived_free(&work->sys);
    md_memset(work, 0, sizeof(*work));
}
