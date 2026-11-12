#include <md_infer.h>
#include <md_element.h>
#include "md_infer_internal.h"

// The inference of a system: md_system_infer runs the steps of md_bond.c, md_ring.c, md_atom.c, md_hierarchy.c,
// md_backbone.c, md_lipid.c, md_fragment.c and md_chem.c in order and gives the system what they infer.

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
    md_system_topology_changed(sys);

    md_allocator_i* alloc = sys->alloc;

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
            md_temp_scope_t temp = md_temp_begin_avoid(sys->alloc);
            md_bond_data_t inferred = {0};
            md_bond_infer(&inferred, sys, state, md_temp_allocator(temp), executor);
            md_bond_merge(&bonds, &sys->bond, &inferred, sys->atom.count, sys->alloc);
            md_temp_end(temp);
        } else {
            md_bond_infer(&bonds, sys, state, sys->alloc, executor);
        }
        md_bond_data_free(&sys->bond, sys->alloc);
        sys->bond = bonds;
        md_system_topology_changed(sys);
    }

    // Whatever the bonds came from (inferred, a topology, a file, by hand): those between a metal and a non-metal are
    // coordination. Cheap, and every reader of the flags relies on it, so not behind a bit.
    md_system_infer_coordination(sys);

    // What the residues are (amino acids, nucleotides, waters, ions) and the roles of their atoms: by their names and by
    // the bonds, so after them. What a loader gives is the atoms and their names; everything below reads these.
    if (flags & MD_SYSTEM_INFER_RESIDUE_BIT) {
        md_system_infer_residue_flags(sys, state);
    }

    // The parts before the molecules: a molecule holds its polymer strands whole (synthesize_molecule_links). They are
    // found from the instances the data gives, the bonds and the classified residues, which the molecules do not change.
    if (flags & MD_SYSTEM_INFER_PART_BIT) {
        md_chain_data_t chains = {0};
        md_chain_infer(&chains, sys, sys->alloc);
        md_chain_data_free(&sys->chain, sys->alloc);
        sys->chain = chains;
        md_part_data_t parts = {0};
        md_part_infer(&parts, sys, sys->alloc);
        md_part_data_free(&sys->part, sys->alloc);
        sys->part = parts;
    }

    if (flags & MD_SYSTEM_INFER_MOLECULE_BIT) {
        // Coarse grained systems get a hierarchy even without bonds (synthesize_hierarchy_links), and a residue is a
        // molecule whatever its bonds (synthesize_molecule_links): ions and waters given without bonds are molecules
        if (sys->bond.count || cg || sys->residue.count) {
            md_molecule_data_t molecules = {0};
            md_molecule_infer(&molecules, sys, sys->alloc);
            md_molecule_data_free(&sys->molecule, sys->alloc);
            sys->molecule = molecules;
            // The fragments are of the molecules, looked up by their slots: gone with them until inferred again
            // (md_util_system_infer_fragments). Their types stay, for the colours to be kept.
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
        if (sys->bond.count) {
            md_ring_data_t rings = {0};
            md_ring_infer(&rings, sys, sys->alloc);
            md_ring_data_free(&sys->ring, sys->alloc);
            sys->ring = rings;
        }
    }
    
    if (flags & MD_SYSTEM_INFER_CHEMISTRY_BIT) {
        // After the rings, before the hydrogen bond roles which read the orders. Protonation at pH 7 applies to
        // residues without hydrogens only.
        if (sys->bond.count) {
            md_chem_perceive(sys, state, MD_CHEM_FLAG_PROTONATE_PH7);
        }
    }

    if (flags & MD_SYSTEM_INFER_BACKBONE_BIT) {
        md_backbone_data_t bb = {0};
        md_backbone_infer(&bb, sys, sys->alloc);
        md_backbone_data_t old = { .protein = sys->protein_backbone, .nucleic = sys->nucleic_backbone };
        md_backbone_data_free(&old, sys->alloc);
        sys->protein_backbone = bb.protein;
        sys->nucleic_backbone = bb.nucleic;
        md_system_set_chain_breaks(sys, bb.breaks, bb.num_breaks);
        md_array_free(bb.breaks, sys->alloc);
    }

    if (flags & MD_SYSTEM_INFER_LIPID_BIT) {
        // After the bonds, the rings (a sterol is told by its fused rings) and the classification of the residues
        md_lipid_data_t lipids = {0};
        md_lipid_infer(&lipids, sys, sys->alloc);
        md_lipid_data_free(&sys->lipid, sys->alloc);
        sys->lipid = lipids;
        md_system_infer_lipid_kinds(sys);
    }

    // Last of what depends on the classification of the residues, which the lipids are the last to change. Cheap,
    // and what a residue type is follows from the residues alone, so not behind a bit.
    md_system_infer_residue_types(sys);

    if (flags & MD_SYSTEM_INFER_FRAGMENT_BIT) {
        // Last: made of the molecules, the rings, the chemistry, the roles of the atoms and the lipid tails
        md_util_system_infer_fragments(sys);
    }

    // Record the state the topology above was inferred from. Written here, as part of performing
    // the inference, so it cannot disagree with what it produced.
    if (state != &sys->reference) {
        sys->reference.alloc = sys->alloc;
        md_system_state_copy(&sys->reference, state);
    }

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
            md_array_resize(work->sys.f, n, alloc); \
            for (size_t i = 0; i < n; ++i) work->sys.f[i] = md_str_empty(sys->f[i]) ? (md_str_t){0} : md_alloc_str(alloc, sys->f[i]); \
        } \
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
    MD_SYSTEM_DERIVED_ARRAYS(INFER_COPY_ARRAY)
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
        work->sys.reference.alloc = alloc;
        if (!md_system_state_copy(&work->sys.reference, state)) return false;
    }
    return true;
}

#undef INFER_COPY_ARRAY
#undef INFER_COPY_STRINGS

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

    md_allocator_i* alloc = sys->alloc;
    // What the system had of the derived: given back (to an arena, which frees nothing, this is nothing)
#define FREE_ARRAY(f) md_array_free(sys->f, alloc);
#define FREE_STRINGS(f) for (size_t i = 0; i < md_array_size(sys->f); ++i) { if (!md_str_empty(sys->f[i])) md_free_str(alloc, sys->f[i]); } md_array_free(sys->f, alloc);
    MD_SYSTEM_DERIVED_ARRAYS(FREE_ARRAY)
    MD_SYSTEM_DERIVED_STRINGS(FREE_STRINGS)
#undef FREE_ARRAY
#undef FREE_STRINGS
    md_system_state_free(&sys->reference);

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
