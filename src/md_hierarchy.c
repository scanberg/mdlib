#include "md_infer_internal.h"

// Inference of the levels of a system above its atoms: chains, instances, parts, molecules and residue types
// (md_util.h).

void md_chain_data_free(md_chain_data_t* chain, md_allocator_i* alloc) {
    MD_ASSERT(chain);
    MD_ASSERT(alloc);
    md_array_free(chain->id, alloc);
    md_array_free(chain->run_offset, alloc);
    md_array_free(chain->run, alloc);
    md_array_free(chain->res_chain, alloc);
    md_memset(chain, 0, sizeof(*chain));
}


static void clear_instances(md_system_t* sys, md_allocator_i* alloc) {
    md_array_free(sys->instance.id, alloc);
    md_array_free(sys->instance.auth_id, alloc);
    md_array_free(sys->instance.res_offset, alloc);
    md_array_free(sys->instance.part_type_idx, alloc);
    md_array_free(sys->instance.mol_type_idx, alloc);
    md_memset(&sys->instance, 0, sizeof(sys->instance));
}

bool md_system_instances_from_chain_ids(md_system_t* sys, const md_str_t res_chain_id[], const bool res_unit_begins[]) {
    if (!sys) {
        MD_LOG_ERROR("Missing system");
        return false;
    }
    MD_ASSERT(sys->alloc);
    md_allocator_i* alloc = sys->alloc;
    clear_instances(sys, alloc);
    md_system_topology_changed(sys);

    const size_t num_res = sys->residue.count;
    if (!res_chain_id || num_res == 0) return true;

    // No ids at all is no units: the data gives none
    bool any_id = false;
    for (size_t c = 0; c < num_res && !any_id; ++c) any_id = !md_str_empty(md_str_trim(res_chain_id[c]));
    if (!any_id) return true;

    // The runs of one id. Every residue is in one, as the instances are back to back; a run without an id is an
    // instance the data does not name (md_instance_data_t)
    size_t i = 0;
    while (i < num_res) {
        const md_str_t id = md_str_trim(res_chain_id[i]);
        size_t j = i + 1;
        while (j < num_res && md_str_eq(id, md_str_trim(res_chain_id[j])) && !(res_unit_begins && res_unit_begins[j])) ++j;
        md_array_push(sys->instance.id, (md_label_t){0}, alloc);
        md_array_push(sys->instance.auth_id, md_label_make(id), alloc);
        md_array_push(sys->instance.res_offset, (uint32_t)i, alloc);
        md_array_push(sys->instance.part_type_idx, -1, alloc);
        md_array_push(sys->instance.mol_type_idx, -1, alloc);
        sys->instance.count += 1;
        i = j;
    }
    md_array_push(sys->instance.res_offset, (uint32_t)num_res, alloc);
    return true;
}

bool md_chain_infer(md_chain_data_t* out, const md_system_t* sys, md_allocator_i* alloc) {
    MD_ASSERT(out);
    MD_ASSERT(sys);
    MD_ASSERT(alloc);
    MD_ASSERT(out != &sys->chain);
    md_chain_data_free(out, alloc);

    const size_t num_inst = sys->instance.count;
    const size_t num_res = sys->residue.count;
    if (num_inst == 0 || num_res == 0) return true;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    // The chain of each id, by the hash of the id
    md_hashmap32_t chain_of_id = { .allocator = temp_arena };
    md_array_resize(out->res_chain, num_res, alloc);
    for (size_t c = 0; c < num_res; ++c) out->res_chain[c] = -1;

    // The runs, in the order they are found, each with its chain. A chain is a label: only what carries an id is in one.
    md_array(md_urange_t) runs = 0;
    md_array(md_chain_idx_t) run_chain = 0;

    for (size_t i = 0; i < num_inst; ++i) {
        const md_str_t id = md_str_trim(md_system_instance_auth_id(sys, i));
        const md_urange_t residues = md_system_instance_residue_range(sys, i);
        if (md_str_empty(id) || residues.beg >= residues.end || residues.end > num_res) continue;
        const uint64_t key = md_hash64_str(id, 0);
        const uint32_t* found = md_hashmap_get(&chain_of_id, key);
        md_chain_idx_t chain;
        if (found) {
            chain = (md_chain_idx_t)*found;
        } else {
            chain = (md_chain_idx_t)out->count++;
            md_array_push(out->id, md_label_make(id), alloc);
            md_hashmap_add(&chain_of_id, key, (uint32_t)chain);
        }
        md_array_push(runs, residues, temp_arena);
        md_array_push(run_chain, chain, temp_arena);
    }

    // The runs of each chain, in the order of the chains, the consecutive ones of a chain joined
    if (out->count) {
        const size_t num_runs = md_array_size(runs);
        uint32_t* count = md_temp_alloc_zero_array(temp, uint32_t, out->count + 1);
        for (size_t r = 0; r < num_runs; ++r) count[run_chain[r] + 1] += 1;
        for (size_t c = 0; c < out->count; ++c) count[c + 1] += count[c];
        md_urange_t* sorted = md_temp_alloc_array(temp, md_urange_t, MD_MAX(1, num_runs));
        uint32_t* fill = md_temp_alloc_array(temp, uint32_t, out->count);
        md_memcpy(fill, count, sizeof(uint32_t) * out->count);
        for (size_t r = 0; r < num_runs; ++r) sorted[fill[run_chain[r]]++] = runs[r];

        md_array_resize(out->run_offset, out->count + 1, alloc);
        out->run_offset[0] = 0;
        for (size_t c = 0; c < out->count; ++c) {
            for (uint32_t r = count[c]; r < count[c + 1]; ++r) {
                const md_urange_t run = sorted[r];
                const size_t n = md_array_size(out->run);
                if (r > count[c] && out->run[n - 1].end == run.beg) {
                    out->run[n - 1].end = run.end;
                } else {
                    md_array_push(out->run, run, alloc);
                }
                for (uint32_t k = run.beg; k < run.end; ++k) out->res_chain[k] = (md_chain_idx_t)c;
            }
            out->run_offset[c + 1] = (uint32_t)md_array_size(out->run);
        }
    }

    md_temp_end(temp);
    return true;
}

static inline bool residue_kind_is_monomer(md_residue_kind_t kind) {
    return kind == MD_RESIDUE_KIND_AMINO_ACID || kind == MD_RESIDUE_KIND_NUCLEOTIDE;
}

// Whether a nucleotide is DNA (1) or RNA (2), 0 when it cannot be told. The name decides when it is one of the
// standard names (DA, DC, ... and A, C, ...). Otherwise (ADE, CYT, ... of CHARMM, which are both) the sugar does: a
// resolved nucleotide with a 2' oxygen is RNA, one without DNA.
static int nucleotide_type(const md_system_t* sys, size_t res_idx) {
    const md_str_t name = md_system_residue_name(sys, res_idx);
    if (md_util_resname_dna(name)) return 1;
    if (md_util_resname_rna(name)) return 2;
    md_nucleic_acid_atoms_t na;
    if ((md_system_residue_flags(sys, res_idx) & MD_RESIDUE_FLAG_RESOLVED) && mdi_residue_nucleotide_atoms(&na, sys, md_system_residue_atom_range(sys, res_idx))) {
        return mdi_bonded_atom_with(sys, na.c2, MD_Z_O, 0, 0, -1, -1) >= 0 ? 2 : 1;
    }
    return 0;
}

// What a part made of the residues in res_range is, from what its residues are. Two or more amino acids
// (nucleotides) which make up at least half of the residues are a peptide (nucleic acid), which lets the caps of a
// chain (ACE, NME) and modified residues the names do not cover come along.
static md_part_kind_t classify_part(const md_system_t* sys, md_urange_t res_range) {
    const size_t n = res_range.end - res_range.beg;
    if (n == 0) return MD_PART_KIND_UNKNOWN;
    if (n == 1) {
        return md_system_residue_kind(sys, res_range.beg) == MD_RESIDUE_KIND_WATER ? MD_PART_KIND_WATER : MD_PART_KIND_NON_POLYMER;
    }

    size_t num_amino = 0, num_nucl = 0, num_dna = 0, num_rna = 0;
    for (uint32_t c = res_range.beg; c < res_range.end; ++c) {
        const md_residue_kind_t kind = md_system_residue_kind(sys, c);
        if (kind == MD_RESIDUE_KIND_AMINO_ACID) {
            num_amino += 1;
        } else if (kind == MD_RESIDUE_KIND_NUCLEOTIDE) {
            num_nucl += 1;
            const int type = nucleotide_type(sys, c);
            num_dna += type == 1;
            num_rna += type == 2;
        }
    }

    if (num_amino >= 2 && 2 * num_amino >= n) {
        return MD_PART_KIND_PEPTIDE;
    }
    if (num_nucl >= 2 && 2 * num_nucl >= n) {
        if (num_dna == num_nucl) return MD_PART_KIND_DNA;
        if (num_rna == num_nucl) return MD_PART_KIND_RNA;
        return MD_PART_KIND_NUCLEIC;
    }
    return MD_PART_KIND_NON_POLYMER;
}

// The given part types are kept, the inferred ones (which come after them) dropped: the parts are about to infer them anew
// The types the data gave (those before the inferred ones), copied into out
static void copy_given_part_types(md_part_type_data_t* out, const md_part_type_data_t* type, md_allocator_i* alloc) {
    size_t n = type->count;
    while (n > 0 && (type->flags[n - 1] & MD_PART_TYPE_FLAG_INFERRED)) --n;
    for (size_t i = 0; i < n; ++i) {
        md_array_push(out->id, type->id[i], alloc);
        md_array_push(out->flags, type->flags[i], alloc);
        md_array_push(out->description, md_alloc_str(alloc, type->description[i]), alloc);
    }
    out->count = n;
}

void md_part_data_free(md_part_data_t* part, md_allocator_i* alloc) {
    MD_ASSERT(part);
    MD_ASSERT(alloc);
    md_array_free(part->res_offset, alloc);
    md_array_free(part->type_idx, alloc);
    md_part_type_data_free(&part->type, alloc);
    md_memset(part, 0, sizeof(*part));
}

// The type of a part which is not all of an instance of a given type: one per kind and sequence of residue names. That
// the type was inferred is carried by its flags, the description only says what it is: the kind of a polymer, the
// residue names of anything else.
static md_part_type_idx_t inferred_part_type(md_part_type_data_t* types, const md_system_t* sys, md_hashmap32_t* type_of_key, md_urange_t residues, md_allocator_i* alloc) {
    const md_part_kind_t kind = classify_part(sys, residues);
    const bool polymer = md_part_kind_is_polymer(kind);

    uint64_t key = md_hash64(&kind, sizeof(kind), 0);
    for (uint32_t c = residues.beg; c < residues.end; ++c) {
        key = md_hash64_str(md_system_residue_name(sys, c), key);
    }
    const uint32_t* found = md_hashmap_get(type_of_key, key);
    if (found) return (md_part_type_idx_t)*found;

    char buf[64];
    md_str_t desc = {0};
    if (polymer || kind == MD_PART_KIND_WATER) {
        desc = md_str_from_cstr(md_part_kind_name(kind));
    } else if (residues.end - residues.beg <= 4) {
        int len = 0;
        for (uint32_t c = residues.beg; c < residues.end; ++c) {
            const md_str_t name = md_system_residue_name(sys, c);
            len += md_snprintf(buf + len, sizeof(buf) - len, "%s" MD_STR_FMT, c > residues.beg ? "-" : "", MD_STR_ARG(name));
            len = MD_MIN(len, (int)sizeof(buf) - 1);
        }
        desc = (md_str_t){ buf, (size_t)len };
    } else {
        desc = md_system_residue_name(sys, residues.beg);
    }
    const md_part_type_idx_t type = md_part_type_data_push(types, MD_STR_LIT(""), md_part_type_flags_set_kind(MD_PART_TYPE_FLAG_INFERRED, kind), desc, alloc);
    md_hashmap_add(type_of_key, key, (uint32_t)type);
    return type;
}

bool md_part_infer(md_part_data_t* out, const md_system_t* sys, md_allocator_i* alloc) {
    MD_ASSERT(out);
    MD_ASSERT(sys);
    MD_ASSERT(alloc);
    MD_ASSERT(out != &sys->part);
    md_part_data_free(out, alloc);
    copy_given_part_types(&out->type, &sys->part.type, alloc);

    const size_t num_res = sys->residue.count;
    if (num_res == 0) return true;

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    // The bonds between residues, but those of a metal coordinated by residues: the residue bonded to the one before
    // it (by an atom of a backbone: a cap, a modified residue in a strand, rather than a glycan on a side chain), and
    // the last residue before each which it is bonded to
    uint64_t* bonded_to_prev = make_bitfield(num_res, temp_arena);
    uint64_t* backbone_to_prev = make_bitfield(num_res, temp_arena);
    int32_t* last_bonded_before = md_temp_alloc_array(temp, int32_t, num_res);
    for (size_t c = 0; c < num_res; ++c) last_bonded_before[c] = -1;
    if (sys->bond.count > 0) {
        md_residue_idx_t* atom_res = md_temp_alloc_array(temp, md_residue_idx_t, MD_MAX(1, sys->atom.count));
        md_memset(atom_res, -1, sizeof(md_residue_idx_t) * sys->atom.count);
        for (size_t c = 0; c < num_res; ++c) {
            const md_urange_t r = md_system_residue_atom_range(sys, c);
            for (uint32_t a = r.beg; a < r.end && a < sys->atom.count; ++a) atom_res[a] = (md_residue_idx_t)c;
        }
        for (size_t b = 0; b < sys->bond.count; ++b) {
            if (sys->bond.flags && (!md_bond_is_graph(sys->bond.flags[b]) || md_bond_kind(sys->bond.flags[b]) == MD_BOND_KIND_COORDINATE)) continue;
            const md_atom_pair_t pair = sys->bond.pairs[b];
            if ((size_t)pair.idx[0] >= sys->atom.count || (size_t)pair.idx[1] >= sys->atom.count) continue;
            const md_residue_idx_t x = atom_res[pair.idx[0]];
            const md_residue_idx_t y = atom_res[pair.idx[1]];
            if (x < 0 || y < 0 || x == y) continue;
            const md_residue_idx_t lo = MD_MIN(x, y), hi = MD_MAX(x, y);
            if (hi - lo == 1) {
                bitfield_set_bit(bonded_to_prev, hi);
                if ((md_atom_flags(&sys->atom, pair.idx[0]) | md_atom_flags(&sys->atom, pair.idx[1])) & MD_ATOM_FLAG_BACKBONE) {
                    bitfield_set_bit(backbone_to_prev, hi);
                }
            }
            last_bonded_before[hi] = MD_MAX(last_bonded_before[hi], lo);
        }
    }

    // The spans the parts are found in: the instances, and the residues before the first and after the last, which no
    // instance holds (all of them in a system without instances). A span which is an instance the data names (an id or
    // an author id) lets a strand go on over a gap while the numbering does.
    typedef struct {
        md_urange_t residues;
        md_instance_idx_t inst;     // -1 for residues in no instance
        bool named;
    } span_t;
    md_array(span_t) spans = 0;
    {
        uint32_t covered = 0;
        for (size_t i = 0; i < sys->instance.count; ++i) {
            md_urange_t r = md_system_instance_residue_range(sys, i);
            r.end = MD_MIN(r.end, (uint32_t)num_res);
            if (r.beg < covered) r.beg = covered;
            if (r.beg >= r.end) continue;
            if (covered < r.beg) {
                const span_t none = {{covered, r.beg}, -1, false};
                md_array_push(spans, none, temp_arena);
            }
            const bool named = !md_str_empty(md_system_instance_id(sys, i)) || !md_str_empty(md_system_instance_auth_id(sys, i));
            const span_t span = {r, (md_instance_idx_t)i, named};
            md_array_push(spans, span, temp_arena);
            covered = r.end;
        }
        if (covered < num_res) {
            const span_t none = {{covered, (uint32_t)num_res}, -1, false};
            md_array_push(spans, none, temp_arena);
        }
    }

    md_hashmap32_t type_of_key = { .allocator = temp_arena };

    md_array_push(out->res_offset, 0, alloc);
    for (size_t s = 0; s < md_array_size(spans); ++s) {
        const span_t span = spans[s];
        const md_urange_t residues = span.residues;
        const size_t first_part = md_array_size(out->res_offset) - 1;

        uint32_t beg = residues.beg;
        bool run_polymer = residue_kind_is_monomer(md_system_residue_kind(sys, beg));
        for (uint32_t j = residues.beg + 1; j <= residues.end; ++j) {
            bool goes_on = false;
            if (j < residues.end) {
                const md_residue_kind_t kj = md_system_residue_kind(sys, j);
                const md_residue_kind_t kp = md_system_residue_kind(sys, j - 1);
                const bool bonded = bitfield_test_bit(bonded_to_prev, j);
                if (kj == MD_RESIDUE_KIND_WATER || kj == MD_RESIDUE_KIND_ION || kp == MD_RESIDUE_KIND_WATER || kp == MD_RESIDUE_KIND_ION) {
                    goes_on = false;
                } else if (residue_kind_is_monomer(kj) && kj == kp) {
                    // A strand, over a gap unless the chain ends there: within an instance the data names while the
                    // numbering goes on, elsewhere only in a system without bonds while it runs on by one
                    if (bonded) {
                        goes_on = true;
                    } else if (!mdi_polymer_chain_ends_between(sys, j - 1, j)) {
                        const md_sequence_id_t seq_prev = md_system_residue_seq_id(sys, j - 1);
                        const md_sequence_id_t seq_j    = md_system_residue_seq_id(sys, j);
                        goes_on = span.named ? seq_j > seq_prev : (sys->bond.count == 0 && seq_j == seq_prev + 1);
                    }
                } else if (!residue_kind_is_monomer(kj) && !run_polymer) {
                    // Of no polymer: on while bonded to a residue of the run (the branches of a glycan)
                    goes_on = bonded || last_bonded_before[j] >= (int32_t)beg;
                } else {
                    // A cap or a modified residue bonded into a strand by its backbone, not a glycan on a side chain
                    goes_on = bitfield_test_bit(backbone_to_prev, j);
                }
            }
            if (goes_on) {
                run_polymer |= residue_kind_is_monomer(md_system_residue_kind(sys, j));
                continue;
            }
            md_array_push(out->res_offset, j, alloc);
            beg = j;
            if (j < residues.end) run_polymer = residue_kind_is_monomer(md_system_residue_kind(sys, j));
        }

        // The types: that of the instance where the part is all of it, else an inferred one
        const size_t end_part = md_array_size(out->res_offset) - 1;
        const md_part_type_idx_t given = span.inst >= 0 ? md_system_instance_part_type_idx(sys, (size_t)span.inst) : -1;
        for (size_t p = first_part; p < end_part; ++p) {
            const md_urange_t r = { out->res_offset[p], out->res_offset[p + 1] };
            md_part_type_idx_t type = -1;
            if (given >= 0 && (size_t)given < out->type.count && end_part - first_part == 1) {
                type = given;
                if (md_part_type_flags_kind(out->type.flags[type]) == MD_PART_KIND_UNKNOWN) {
                    // A given type which does not say what it is (a polymer entity without _entity_poly)
                    out->type.flags[type] = md_part_type_flags_set_kind(out->type.flags[type], classify_part(sys, r));
                }
            } else {
                type = inferred_part_type(&out->type, sys, &type_of_key, r, alloc);
            }
            md_array_push(out->type_idx, type, alloc);
        }
    }
    out->count = md_array_size(out->res_offset) - 1;
    MD_ASSERT(md_array_size(out->type_idx) == out->count);

    md_temp_end(temp);
    return true;
}

// The molecules (md_molecule_data_t) are the connected components of the bond graph, joined where the bonds do not
// hold a molecule together (see there).
// Breadth first traversal of the connected component containing 'seed'.
//
// Marks every atom it reaches in 'visited' and never clears it; the caller supplies a bitfield which
// is clear over the component. When out_atoms is supplied the atoms are written in BFS order and
// out_count receives how many there were. When out_pred is supplied it records, per global atom
// index, the atom each one was reached from, with the seed recording itself.
//
// Returns the last atom popped from the queue. BFS pops in non decreasing distance order, so that
// atom is at maximum edge distance from the seed - which is what the double sweep below needs.
static int bfs_traverse(int32_t* out_atoms, size_t* out_count, int32_t* out_pred, uint64_t* visited, md_fifo_t* queue, const md_bond_data_t* bond, int seed) {
    md_fifo_clear(queue);

    bitfield_set_bit(visited, seed);
    if (out_pred) out_pred[seed] = seed;
    md_fifo_push(queue, seed);

    size_t count = 0;
    int last = seed;

    while (!md_fifo_empty(queue)) {
        const int cur = md_fifo_pop(queue);
        last = cur;
        if (out_atoms) out_atoms[count] = cur;
        count += 1;

        md_bond_iter_t it = md_bond_iter(bond, cur);
        while (md_bond_iter_has_next(&it)) {
            const int next = md_bond_iter_atom_index(&it);
            md_bond_iter_next(&it);

            if (bitfield_test_bit(visited, next)) continue;

            bitfield_set_bit(visited, next);
            if (out_pred) out_pred[next] = cur;
            md_fifo_push(queue, next);
        }
    }

    if (out_count) *out_count = count;
    return last;
}

// Whether a bond holds a molecule together (md_molecule_data_t): a bond of the graph, but not a coordination bond
// between two residues, a metal and what surrounds it (the Mg of a GTP with its phosphates and the waters of its shell,
// the Zn of a zinc finger with its cysteines), which are molecules of their own. Unless a topology gives it: a bond
// within one of its molecules, as GROMACS has it. A coordination bond within a residue (the Fe of a heme) holds it.
static bool bond_holds_molecule(const md_system_t* sys, size_t b) {
    if (!sys->bond.flags) return true;
    const md_bond_flags_t flags = sys->bond.flags[b];
    if (!md_bond_is_graph(flags)) return false;
    if (md_bond_kind(flags) != MD_BOND_KIND_COORDINATE || md_bond_origin(flags) == MD_BOND_ORIGIN_TOPOLOGY) return true;
    if (!sys->residue.count || !sys->residue.atom_offset) return false;
    const md_atom_pair_t pair = sys->bond.pairs[b];
    return md_system_residue_find_by_atom_idx(sys, (size_t)pair.idx[0]) == md_system_residue_find_by_atom_idx(sys, (size_t)pair.idx[1]);
}

// The molecules of a topology are whole (md_instance_data_t::mol_type_idx): the atoms of each of its instances, joined
// to the first where nothing else joins them (the two strands of a duplex which martinize-dna writes as one molecule
// type, held together by an elastic network which is no bond of the graph)
static void link_topology_molecules(const md_system_t* sys, int32_t* set, md_array(md_atom_pair_t)* links, md_allocator_i* alloc) {
    if (!sys->instance.mol_type_idx) return;
    const size_t atom_count = md_system_atom_count(sys);
    for (size_t inst = 0; inst < sys->instance.count; ++inst) {
        if (md_system_instance_molecule_type_idx(sys, inst) < 0) continue;
        md_urange_t r = md_system_instance_atom_range(sys, inst);
        if (r.end > atom_count) r.end = (uint32_t)atom_count;
        for (uint32_t a = r.beg + 1; a < r.end; ++a) {
            if (uf_union(set, (int32_t)r.beg, (int32_t)a)) {
                md_atom_pair_t pair = {{ (int32_t)r.beg, (int32_t)a }};
                md_array_push(*links, pair, alloc);
            }
        }
    }
}

// Hierarchy links for coarse grained systems, used where bonds leave the topology disconnected.
//
// md_molecule_data_t is a spanning forest and unwrap only needs each atom's parent to stay within half a cell
// of it, not to be bonded to it. Coarse grained systems frequently have no bonds (md_util_system_infer skips
// bond inference for them) or have a distance guess that misses stretched ones. So the missing links are
// synthesized from the residue hierarchy instead, WITHOUT looking at coordinates:
//   - within a residue, every bead hangs off the residue's anchor (its backbone bead, else its first bead)
//   - the anchors of consecutive polymer residues (amino acid or nucleotide, same kind, consecutive sequence
//     ids, or of one part whatever their ids when the parts are known) are linked along the chain
// A link is only emitted between atoms that are not already connected, by bonds or by an earlier link, so a
// fully bonded system gets none and the result stays a forest over bonds + links.
//
// Deliberately NOT done: distance based linking. In a dense coarse grained membrane or solvent that percolates
// into one box spanning structure, and unwrapping that is meaningless. Non polymer residues (lipids, sterols)
// are not chained together either - a span of residues may hold hundreds of them in file order, not spatial order.
//
// The links are returned as pairs and never written to sys->bond: they are not bonds, must not be drawn and
// must not feed ring detection.
static md_array(md_atom_pair_t) synthesize_hierarchy_links(const md_system_t* sys, md_allocator_i* alloc) {
    md_array(md_atom_pair_t) links = 0;

    const size_t atom_count = md_system_atom_count(sys);
    const size_t res_count = sys->residue.count;
    if (atom_count == 0 || res_count == 0) {
        return links;
    }

    int32_t* set = md_alloc(alloc, sizeof(int32_t) * atom_count);
    for (size_t i = 0; i < atom_count; ++i) set[i] = (int32_t)i;

    for (size_t i = 0; i < sys->bond.count; ++i) {
        if (!bond_holds_molecule(sys, i)) continue;
        const md_atom_pair_t pair = sys->bond.pairs[i];
        if ((size_t)pair.idx[0] < atom_count && (size_t)pair.idx[1] < atom_count) {
            uf_union(set, pair.idx[0], pair.idx[1]);
        }
    }

    // Part of each residue, -1 when the parts are unknown (md_util_system_infer infers them before the molecules)
    int32_t* res_part = md_alloc(alloc, sizeof(int32_t) * res_count);
    for (size_t i = 0; i < res_count; ++i) res_part[i] = -1;
    for (size_t p = 0; p < sys->part.count; ++p) {
        md_urange_t range = md_system_part_residue_range(sys, p);
        for (uint32_t ci = range.beg; ci < range.end && ci < res_count; ++ci) {
            res_part[ci] = (int32_t)p;
        }
    }

    int32_t             prev_anchor = -1;
    size_t              prev_res    = 0;
    md_residue_kind_t prev_kind     = MD_RESIDUE_KIND_OTHER;

    for (size_t ci = 0; ci < res_count; ++ci) {
        md_urange_t range = md_system_residue_atom_range(sys, ci);
        if (range.end > atom_count) range.end = (uint32_t)atom_count;
        if (range.beg >= range.end) {
            prev_anchor = -1;
            continue;
        }

        int32_t anchor = (int32_t)range.beg;
        for (uint32_t i = range.beg; i < range.end; ++i) {
            if (md_atom_flags(&sys->atom, i) & MD_ATOM_FLAG_BACKBONE) {
                anchor = (int32_t)i;
                break;
            }
        }

        for (uint32_t i = range.beg; i < range.end; ++i) {
            if ((int32_t)i != anchor && uf_union(set, anchor, (int32_t)i)) {
                md_atom_pair_t pair = {{ anchor, (int32_t)i }};
                md_array_push(links, pair, alloc);
            }
        }

        const md_residue_kind_t ck = md_system_residue_kind(sys, ci);
        const md_residue_kind_t kind = residue_kind_is_monomer(ck) ? ck : MD_RESIDUE_KIND_OTHER;
        // Within a part the residues are of one strand whatever their numbering (a strand goes on over a gap,
        // md_part_infer), as a molecule holds its parts whole; without parts by consecutive numbers
        const bool same_part = res_part[prev_res] >= 0 && res_part[prev_res] == res_part[ci];
        const bool consecutive = !sys->residue.seq_id || md_system_residue_seq_id(sys, ci) == md_system_residue_seq_id(sys, prev_res) + 1;
        if (kind && prev_anchor >= 0 && prev_res + 1 == ci && prev_kind == kind &&
            res_part[prev_res] == res_part[ci] && (same_part || consecutive))
        {
            if (uf_union(set, prev_anchor, anchor)) {
                md_atom_pair_t pair = {{ prev_anchor, anchor }};
                md_array_push(links, pair, alloc);
            }
        }

        prev_anchor = kind ? anchor : -1;
        prev_res    = ci;
        prev_kind   = kind;
    }

    link_topology_molecules(sys, set, &links, alloc);
    return links;
}

// The atoms which join two residues of a polymer chain over a gap: the C and N of a peptide bond, the O3' and P of a
// phosphodiester bond (the O5' where the second has no phosphate: the 5' end of a strand, a phosphate not modelled),
// else the last atom of the first and the first of the second. Never an atom the residues do not have.
static md_atom_pair_t chain_gap_link(const md_system_t* sys, size_t res_idx) {
    const md_urange_t range_a = md_system_residue_atom_range(sys, res_idx);
    const md_urange_t range_b = md_system_residue_atom_range(sys, res_idx + 1);
    md_atom_pair_t pair = {{ (int32_t)range_a.end - 1, (int32_t)range_b.beg }};
    md_amino_acid_atoms_t aa_a, aa_b;
    md_nucleic_acid_atoms_t na_a, na_b;
    if (mdi_residue_amino_acid_atoms(&aa_a, sys, range_a) && mdi_residue_amino_acid_atoms(&aa_b, sys, range_b)) {
        if (aa_a.c >= 0) pair.idx[0] = aa_a.c;
        if (aa_b.n >= 0) pair.idx[1] = aa_b.n;
    } else if (mdi_residue_nucleotide_atoms(&na_a, sys, range_a) && mdi_residue_nucleotide_atoms(&na_b, sys, range_b)) {
        if (na_a.o3 >= 0) pair.idx[0] = na_a.o3;
        if (na_b.p >= 0) pair.idx[1] = na_b.p;
        else if (na_b.o5 >= 0) pair.idx[1] = na_b.o5;
    }
    return pair;
}

// The links which make a molecule whole where its bonds leave it in pieces (see md_molecule_data_t). The atoms of a
// residue are of one molecule: an atom which no bond reaches (the virtual site of a TIP4P water, an ion a file names
// with the residue it sits in) is linked to the first heavy atom of its residue. The residues of a polymer chain are of
// one molecule: over a gap (residues missing from a crystal structure, a chain of a topology whose bonds are not
// given) its two ends are linked. Between the pieces only: a residue or a chain its bonds join gets no link.
static md_array(md_atom_pair_t) synthesize_molecule_links(const md_system_t* sys, md_allocator_i* alloc) {
    md_array(md_atom_pair_t) links = 0;
    const size_t atom_count = md_system_atom_count(sys);
    if (atom_count == 0 || sys->residue.count == 0) return links;

    int32_t* set = md_alloc(alloc, sizeof(int32_t) * atom_count);
    for (size_t i = 0; i < atom_count; ++i) set[i] = (int32_t)i;
    for (size_t i = 0; i < sys->bond.count; ++i) {
        if (!bond_holds_molecule(sys, i)) continue;
        const md_atom_pair_t pair = sys->bond.pairs[i];
        if ((size_t)pair.idx[0] < atom_count && (size_t)pair.idx[1] < atom_count) {
            uf_union(set, pair.idx[0], pair.idx[1]);
        }
    }

    for (size_t ci = 0; ci < sys->residue.count; ++ci) {
        md_urange_t range = md_system_residue_atom_range(sys, ci);
        if (range.end > atom_count) range.end = (uint32_t)atom_count;
        if (range.beg >= range.end) continue;
        int32_t anchor = -1;
        for (uint32_t i = range.beg; i < range.end && anchor < 0; ++i) {
            if (md_atom_particle_kind(&sys->atom, i) != MD_PARTICLE_VIRTUAL_SITE && md_atom_atomic_number(&sys->atom, i) > MD_Z_H) anchor = (int32_t)i;
        }
        if (anchor < 0) anchor = (int32_t)range.beg;
        for (uint32_t i = range.beg; i < range.end; ++i) {
            if (uf_union(set, anchor, (int32_t)i)) {
                md_atom_pair_t pair = {{ anchor, (int32_t)i }};
                md_array_push(links, pair, alloc);
            }
        }
    }

    for (size_t p = 0; p < sys->part.count; ++p) {
        if (!md_part_kind_is_polymer(md_system_part_kind(sys, p))) continue;
        const md_urange_t residues = md_system_part_residue_range(sys, p);
        for (uint32_t ci = residues.beg; ci + 1 < residues.end; ++ci) {
            const md_urange_t ra = md_system_residue_atom_range(sys, ci);
            const md_urange_t rb = md_system_residue_atom_range(sys, ci + 1);
            if (ra.beg >= ra.end || rb.beg >= rb.end || rb.end > atom_count) continue;
            if (uf_find(set, (int32_t)ra.beg) == uf_find(set, (int32_t)rb.beg)) continue;
            const md_atom_pair_t pair = chain_gap_link(sys, ci);
            MD_ASSERT(pair.idx[0] >= 0 && pair.idx[1] >= 0);
            uf_union(set, pair.idx[0], pair.idx[1]);
            md_array_push(links, pair, alloc);
        }
    }
    link_topology_molecules(sys, set, &links, alloc);
    return links;
}

// The name of an inferred molecule type: what its parts are ("peptide + NAG-NAG"), or its formula where it has no parts
// (an XYZ file), C and H first and the other elements by symbol
static void molecule_type_name(char* buf, size_t cap, const md_system_t* sys, const md_part_idx_t* parts, size_t num_parts, const md_molecule_t* mol) {
    MD_ASSERT(cap > 0);
    buf[0] = '\0';
    int len = 0;
    if (num_parts > 0) {
        const size_t MAX_SHOWN = 3;
        for (size_t i = 0; i < num_parts && i < MAX_SHOWN; ++i) {
            const md_str_t desc = md_system_part_type_description(sys, (size_t)md_system_part_type_idx(sys, (size_t)parts[i]));
            len += md_snprintf(buf + len, cap - len, "%s" MD_STR_FMT, i > 0 ? " + " : "", MD_STR_ARG(desc));
            len = MD_MIN(len, (int)cap - 1);
        }
        if (num_parts > MAX_SHOWN) {
            len += md_snprintf(buf + len, cap - len, " + %zu more", num_parts - MAX_SHOWN);
        }
        return;
    }

    uint32_t count[MD_Z_Count] = {0};
    for (size_t i = 0; i < mol->count; ++i) {
        const md_atomic_number_t z = md_atom_atomic_number(&sys->atom, (size_t)mol->atom_idx[i]);
        if (z < MD_Z_Count) count[z] += 1;
    }
    const md_atomic_number_t first[2] = {MD_Z_C, MD_Z_H};
    for (int k = 0; k < 2; ++k) {
        if (count[MD_Z_C] == 0 && first[k] == MD_Z_H) break;    // Without carbon, H goes in order (Hill)
        const md_atomic_number_t z = first[k];
        if (!count[z]) continue;
        const md_str_t sym = md_atomic_number_symbol(z);
        len += count[z] > 1 ? md_snprintf(buf + len, cap - len, MD_STR_FMT "%u", MD_STR_ARG(sym), count[z]) : md_snprintf(buf + len, cap - len, MD_STR_FMT, MD_STR_ARG(sym));
        len = MD_MIN(len, (int)cap - 1);
        count[z] = 0;
    }
    // The rest by symbol: few elements, a selection each time is cheap
    for (;;) {
        md_atomic_number_t best = 0;
        for (int z = 1; z < MD_Z_Count; ++z) {
            if (count[z] && (!best || md_str_cmp_lex(md_atomic_number_symbol((md_atomic_number_t)z), md_atomic_number_symbol(best)) < 0)) best = (md_atomic_number_t)z;
        }
        if (!best) break;
        const md_str_t sym = md_atomic_number_symbol(best);
        len += count[best] > 1 ? md_snprintf(buf + len, cap - len, MD_STR_FMT "%u", MD_STR_ARG(sym), count[best]) : md_snprintf(buf + len, cap - len, MD_STR_FMT, MD_STR_ARG(sym));
        len = MD_MIN(len, (int)cap - 1);
        count[best] = 0;
    }
    if (len == 0) {
        md_snprintf(buf, cap, "%zu particles", mol->count);
    }
}

// The types of the molecules (md_molecule_type_data_t). The molecule which is all of an instance of a given type (a
// molecule of a topology) is of that type. The others are of inferred types, one for the molecules made alike: of the
// same part types in the same order, as many atoms of each particle type, and as many bonds between each pair of
// particle types. The inferred types of the last inference are dropped first.
static void infer_molecule_types(md_molecule_data_t* mol, const md_system_t* sys, md_allocator_i* alloc) {
    // The types the data gave (those before the inferred ones)
    {
        const md_molecule_type_data_t* type = &sys->molecule.type;
        size_t n = type->count;
        while (n > 0 && (type->flags[n - 1] & MD_MOLECULE_TYPE_FLAG_INFERRED)) --n;
        for (size_t i = 0; i < n; ++i) md_molecule_type_data_push(&mol->type, type->name[i], type->flags[i], alloc);
    }
    md_array_resize(mol->type_idx, mol->count, alloc);
    if (mol->count == 0) return;

    const size_t atom_count = sys->atom.count;
    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    int32_t* mol_of_atom = md_temp_alloc_array(temp, int32_t, MD_MAX(1, atom_count));
    for (size_t m = 0; m < mol->count; ++m) {
        for (uint32_t slot = mol->offset[m]; slot < mol->offset[m + 1]; ++slot) {
            mol_of_atom[mol->atom_idx[slot]] = (int32_t)m;
        }
    }

    // What each molecule is made of: its part types in order (hashed in sequence), and its particle types and bonds
    // (summed, as the order of its atoms and bonds is the file's)
    uint64_t* seq_key  = md_temp_alloc_zero_array(temp, uint64_t, mol->count);
    uint64_t* atom_sum = md_temp_alloc_zero_array(temp, uint64_t, mol->count);
    uint64_t* bond_sum = md_temp_alloc_zero_array(temp, uint64_t, mol->count);
    uint32_t* first_part = md_temp_alloc_array(temp, uint32_t, mol->count);
    uint32_t* num_parts  = md_temp_alloc_zero_array(temp, uint32_t, mol->count);
    for (size_t m = 0; m < mol->count; ++m) first_part[m] = UINT32_MAX;

    for (size_t p = 0; p < sys->part.count; ++p) {
        const md_urange_t atoms = md_system_part_atom_range(sys, p);
        if (atoms.beg >= atoms.end || atoms.beg >= atom_count) continue;
        const int32_t m = mol_of_atom[atoms.beg];
        const md_part_type_idx_t t = md_system_part_type_idx(sys, p);
        seq_key[m] = md_hash64(&t, sizeof(t), seq_key[m] + 1);
        if (first_part[m] == UINT32_MAX) first_part[m] = (uint32_t)p;
        num_parts[m] += 1;
    }
    for (size_t a = 0; a < atom_count; ++a) {
        const md_particle_type_idx_t t = sys->atom.type_idx ? sys->atom.type_idx[a] : 0;
        atom_sum[mol_of_atom[a]] += md_hash64(&t, sizeof(t), 0x9E3779B97F4A7C15ull);
    }
    for (size_t b = 0; b < sys->bond.count; ++b) {
        if (!bond_holds_molecule(sys, b)) continue;
        const md_atom_pair_t pair = sys->bond.pairs[b];
        if ((size_t)pair.idx[0] >= atom_count || (size_t)pair.idx[1] >= atom_count) continue;
        const md_particle_type_idx_t ta = sys->atom.type_idx ? sys->atom.type_idx[pair.idx[0]] : 0;
        const md_particle_type_idx_t tb = sys->atom.type_idx ? sys->atom.type_idx[pair.idx[1]] : 0;
        const uint32_t key = ((uint32_t)MD_MIN(ta, tb) << 16) | (uint32_t)MD_MAX(ta, tb);
        bond_sum[mol_of_atom[pair.idx[0]]] += md_hash64(&key, sizeof(key), 0xC2B2AE3D27D4EB4Full);
    }

    md_hashmap32_t type_of_key = { .allocator = temp_arena };
    md_array(md_part_idx_t) parts = 0;
    char name[128];

    for (size_t m = 0; m < mol->count; ++m) {
        md_molecule_t molecule = {0};
        molecule.count      = mol->offset[m + 1] - mol->offset[m];
        molecule.atom_idx   = mol->atom_idx + mol->offset[m];
        molecule.parent_idx = mol->parent_idx + mol->offset[m];

        // All of one instance of a given type
        md_molecule_type_idx_t type = -1;
        if (sys->instance.mol_type_idx && molecule.count > 0) {
            int32_t lo = INT32_MAX, hi = -1;
            for (size_t i = 0; i < molecule.count; ++i) {
                lo = MD_MIN(lo, molecule.atom_idx[i]);
                hi = MD_MAX(hi, molecule.atom_idx[i]);
            }
            const md_instance_idx_t inst = md_system_instance_find_by_atom_idx(sys, (size_t)lo);
            if (inst >= 0) {
                const md_urange_t r = md_system_instance_atom_range(sys, (size_t)inst);
                const md_molecule_type_idx_t given = md_system_instance_molecule_type_idx(sys, (size_t)inst);
                if (given >= 0 && (size_t)given < mol->type.count && r.beg == (uint32_t)lo && r.end == (uint32_t)hi + 1 && r.end - r.beg == molecule.count) {
                    type = given;
                }
            }
        }

        if (type < 0) {
            const size_t n = molecule.count;
            const uint64_t key = md_hash64_combine(md_hash64(&n, sizeof(n), seq_key[m]), md_hash64_combine(atom_sum[m], bond_sum[m]));
            const uint32_t* found = md_hashmap_get(&type_of_key, key);
            if (found) {
                type = (md_molecule_type_idx_t)*found;
            } else {
                md_array_shrink(parts, 0);
                for (uint32_t p = first_part[m]; p < sys->part.count && md_array_size(parts) < num_parts[m]; ++p) {
                    const md_urange_t atoms = md_system_part_atom_range(sys, p);
                    if (atoms.beg < atoms.end && atoms.beg < atom_count && mol_of_atom[atoms.beg] == (int32_t)m) {
                        md_array_push(parts, (md_part_idx_t)p, temp_arena);
                    }
                }
                molecule_type_name(name, sizeof(name), sys, parts, md_array_size(parts), &molecule);
                type = md_molecule_type_data_push(&mol->type, md_str_from_cstr(name), MD_MOLECULE_TYPE_FLAG_INFERRED, alloc);
                md_hashmap_add(&type_of_key, key, (uint32_t)type);
            }
        }
        mol->type_idx[m] = type;
    }

    md_temp_end(temp);
}

void md_molecule_data_free(md_molecule_data_t* mol, md_allocator_i* alloc) {
    MD_ASSERT(mol);
    MD_ASSERT(alloc);
    md_array_free(mol->offset, alloc);
    md_array_free(mol->atom_idx, alloc);
    md_array_free(mol->parent_idx, alloc);
    md_array_free(mol->atom_slot, alloc);
    md_array_free(mol->type_idx, alloc);
    md_molecule_type_data_free(&mol->type, alloc);
    md_memset(mol, 0, sizeof(*mol));
}

bool md_molecule_infer(md_molecule_data_t* out, const md_system_t* sys, md_allocator_i* alloc) {
    MD_ASSERT(out);
    MD_ASSERT(sys);
    MD_ASSERT(alloc);
    MD_ASSERT(out != &sys->molecule);
    md_molecule_data_free(out, alloc);

    md_temp_scope_t temp = md_temp_begin_avoid(alloc);
    md_allocator_i* temp_arena = md_temp_allocator(temp);

    const size_t atom_count = md_system_atom_count(sys);

    // Reverse map from atom to slot, one entry per atom in the system. Every atom is in one molecule: the slots are as
    // many as the atoms, and reserved at once rather than grown (the system's allocator is often an arena, which keeps
    // what a growing array leaves behind)
    md_array_resize(out->atom_slot, atom_count, alloc);
    md_array_ensure(out->atom_idx, atom_count, alloc);
    md_array_ensure(out->parent_idx, atom_count, alloc);

    // Which atoms have been assigned to a molecule. Persists across components.
    uint64_t* visited = make_bitfield(atom_count, temp_arena);

    // Reused by the sweeps below. Cleared over the atoms of the current component only, so the cost
    // stays proportional to the component rather than to the system.
    uint64_t* sweep_visited = make_bitfield(atom_count, temp_arena);

    // Atoms of the residue currently being processed, and the predecessor of each atom as recorded
    // by the second sweep. Only the entries belonging to the current component are meaningful.
    int32_t* component = md_temp_alloc_zero_array(temp, int32_t, atom_count);
    int32_t* pred      = md_temp_alloc_zero_array(temp, int32_t, atom_count);

    // Parallel queues:
    // - atom_queue   : global atom index to visit
    // - parent_queue : global index of the atom it was reached from (the root is its own parent)
    md_fifo_t atom_queue   = md_fifo_create(1024, temp_arena);
    md_fifo_t parent_queue = md_fifo_create(1024, temp_arena);

    // The traversal runs over the bonds which hold molecules together (bond_holds_molecule: not the coordination of a
    // metal by other residues), plus hierarchy links for coarse grained systems where bonds leave beads disconnected
    // (see synthesize_hierarchy_links), and otherwise the links which keep a residue and a polymer chain in one molecule
    // where no bond does (see synthesize_molecule_links). The links only shape the forest; sys->bond is untouched.
    const md_bond_data_t* graph = &sys->bond;
    md_bond_data_t link_graph = {0};
    {
        md_array(md_atom_pair_t) links = md_system_is_coarse_grained(sys) ? synthesize_hierarchy_links(sys, temp_arena) : synthesize_molecule_links(sys, temp_arena);
        const size_t num_links = md_array_size(links);
        // The bonds which hold molecules together (bond_holds_molecule): the coordination bonds between residues left out
        size_t num_held = 0;
        for (size_t b = 0; b < sys->bond.count; ++b) num_held += bond_holds_molecule(sys, b);
        if (num_links > 0 || num_held < sys->bond.count) {
            link_graph.count = num_held + num_links;
            md_array_resize(link_graph.pairs, link_graph.count, temp_arena);
            md_array_resize(link_graph.flags, link_graph.count, temp_arena);
            md_memset(link_graph.flags, 0, sizeof(md_bond_flags_t) * link_graph.count);
            size_t n = 0;
            for (size_t b = 0; b < sys->bond.count; ++b) {
                if (!bond_holds_molecule(sys, b)) continue;
                link_graph.pairs[n] = sys->bond.pairs[b];
                if (sys->bond.flags) link_graph.flags[n] = sys->bond.flags[b];
                n += 1;
            }
            if (num_links) md_memcpy(link_graph.pairs + n, links, sizeof(md_atom_pair_t) * num_links);
            md_bond_build_connectivity(&link_graph, atom_count, temp_arena);
            graph = &link_graph;
        }
    }

    // Sentinel start offset
    md_array_push(out->offset, 0, alloc);

    for (int i = 0; i < (int)atom_count; ++i) {
        if (bitfield_test_bit(visited, i)) continue;

        // Sweep 1: discover the component from an arbitrary seed and note an atom at maximum distance
        // from it. This is the sweep that claims the atoms in 'visited'; the atom list it produces is
        // what lets the following sweeps clear their bits without touching the rest of the system.
        size_t component_size = 0;
        const int far_a = bfs_traverse(component, &component_size, NULL, visited, &atom_queue, graph, i);

        // Locate the graph center: the atom whose greatest distance to any other atom in the
        // component is smallest. Standard double sweep - the midpoint of a longest path. Exact for
        // acyclic structures, a close approximation in the presence of rings.
        // Components of two atoms or fewer are their own center.
        int root = i;
        if (component_size > 2) {
            // Sweep 2: from far_a, find the opposite end of the path and the chain leading back
            for (size_t k = 0; k < component_size; ++k) bitfield_clear_bit(sweep_visited, component[k]);
            const int far_b = bfs_traverse(NULL, NULL, pred, sweep_visited, &atom_queue, graph, far_a);

            int length = 0;
            for (int a = far_b; a != pred[a]; a = pred[a]) length += 1;

            root = far_b;
            for (int k = 0; k < length / 2; ++k) root = pred[root];
        }

        // Sweep 3: emit the component in BFS order from the root, which establishes the
        // slot(parent) < slot(child) invariant md_molecule_data_t documents.
        for (size_t k = 0; k < component_size; ++k) bitfield_clear_bit(sweep_visited, component[k]);

        md_fifo_clear(&atom_queue);
        md_fifo_clear(&parent_queue);

        bitfield_set_bit(sweep_visited, root);
        md_fifo_push(&atom_queue, root);
        md_fifo_push(&parent_queue, root);   // the root is its own parent

        while (!md_fifo_empty(&atom_queue)) {
            const int cur        = md_fifo_pop(&atom_queue);
            const int parent_idx = md_fifo_pop(&parent_queue);

            out->atom_slot[cur] = (int32_t)md_array_size(out->atom_idx);
            md_array_push(out->atom_idx, cur, alloc);
            md_array_push(out->parent_idx, parent_idx, alloc);

            md_bond_iter_t it = md_bond_iter(graph, cur);
            while (md_bond_iter_has_next(&it)) {
                const int next = md_bond_iter_atom_index(&it);
                md_bond_iter_next(&it);

                if (bitfield_test_bit(sweep_visited, next)) continue;

                bitfield_set_bit(sweep_visited, next);
                md_fifo_push(&atom_queue, next);
                md_fifo_push(&parent_queue, cur);
            }
        }

        out->count += 1;
        md_array_push(out->offset, (uint32_t)md_array_size(out->atom_idx), alloc);
    }
    
    md_temp_end(temp);

    infer_molecule_types(out, sys, alloc);
    return true;
}

// ### RESIDUE TYPES ###

// The color a new component type starts out with: that of its standard monomer, that of the element of a monatomic
// ion, and otherwise one of its name, the same in every system
static uint32_t residue_type_default_color(const md_system_t* sys, size_t res_idx, md_str_t name, md_residue_kind_t kind, md_monomer_t monomer) {
    if (monomer != MD_MONOMER_NONE) {
        return md_monomer_color(monomer);
    }
    if (kind == MD_RESIDUE_KIND_ION) {
        const md_urange_t range = md_system_residue_atom_range(sys, res_idx);
        const md_atomic_number_t z = range.end > range.beg ? md_atom_atomic_number(&sys->atom, range.beg) : 0;
        if (z) return md_atomic_number_cpk_color(z);
    }
    // The golden ratio spreads the hues of names which hash alike
    const uint32_t h = (uint32_t)md_hash64_str(name, 0) * 0x9E3779B9u;
    const float hue = (float)(h >> 8) * (1.0f / 16777216.0f);
    const md_f32x3_t rgb = hcl_to_rgb(hue, 0.5f, 0.7f);    // Channels within [0.45, 0.95]
    const uint32_t r = (uint32_t)(MD_CLAMP(rgb.x, 0.0f, 1.0f) * 255.0f + 0.5f);
    const uint32_t g = (uint32_t)(MD_CLAMP(rgb.y, 0.0f, 1.0f) * 255.0f + 0.5f);
    const uint32_t b = (uint32_t)(MD_CLAMP(rgb.z, 0.0f, 1.0f) * 255.0f + 0.5f);
    return (0xFFu << 24) | (b << 16) | (g << 8) | r;
}

void md_system_set_residue_kind(md_system_t* sys, size_t res_idx, md_residue_kind_t kind) {
    MD_ASSERT(sys);
    md_residue_data_t* res = &sys->residue;
    if (res_idx >= res->count) return;
    MD_ASSERT(sys->alloc);
    if (md_array_size(res->type_idx) < res->count) {
        const size_t size = md_array_size(res->type_idx);
        md_array_resize(res->type_idx, res->count, sys->alloc);
        for (size_t i = size; i < res->count; ++i) res->type_idx[i] = -1;
    }
    const md_label_t name = res->name[res_idx];
    const md_residue_type_idx_t cur = res->type_idx[res_idx];
    if (cur >= 0 && (size_t)cur < res->type.count && res->type.kind[cur] == kind && md_label_eq(res->type.name[cur], name)) return;

    // The type of the residue before it is most often the one (a solvent, a membrane), else one of the types
    md_residue_type_idx_t type_idx = -1;
    const md_residue_type_idx_t prev = res_idx > 0 ? res->type_idx[res_idx - 1] : -1;
    if (prev >= 0 && (size_t)prev < res->type.count && res->type.kind[prev] == kind && md_label_eq(res->type.name[prev], name)) {
        type_idx = prev;
    }
    for (size_t t = 0; t < res->type.count && type_idx < 0; ++t) {
        if (res->type.kind[t] == kind && md_label_eq(res->type.name[t], name)) type_idx = (md_residue_type_idx_t)t;
    }
    if (type_idx < 0) {
        const md_str_t name_str = MD_LABEL_STR(res->name[res_idx]);
        const md_monomer_t monomer = md_monomer_from_residue_name(name_str, kind);
        const uint32_t color = residue_type_default_color(sys, res_idx, name_str, kind, monomer);
        type_idx = md_residue_type_add(&res->type, name_str, kind, monomer, color, sys->alloc);
    }
    res->type_idx[res_idx] = type_idx;
}

bool md_system_infer_residue_types(md_system_t* sys) {
    MD_ASSERT(sys);
    if (sys->residue.count == 0) return false;
    if (!sys->alloc) {
        MD_LOG_ERROR("System allocator not set");
        return false;
    }
    for (size_t c = 0; c < sys->residue.count; ++c) {
        if (md_system_residue_type_idx(sys, c) < 0) md_system_set_residue_kind(sys, c, MD_RESIDUE_KIND_OTHER);
    }
    return true;
}
