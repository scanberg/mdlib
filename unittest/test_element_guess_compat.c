#include "utest.h"
#include <md_util.h>
#include <md_molecule.h>
#include <core/md_allocator.h>
#include <core/md_arena_allocator.h>
#include <core/md_str.h>

// Test backward compatibility with md_util_element_guess
UTEST(element_guess_compat, basic_inference) {
    md_allocator_i* alloc = md_vm_arena_create(MEGABYTES(1));
    
    // Create a simple molecule structure
    md_molecule_t mol = {0};
    
    // Atom data
    const size_t atom_count = 5;
    mol.atom.count = atom_count;
    
    // Allocate arrays
    mol.atom.type = md_alloc(alloc, sizeof(md_label_t) * atom_count);
    mol.atom.resname = md_alloc(alloc, sizeof(md_label_t) * atom_count);
    mol.atom.element = md_alloc(alloc, sizeof(md_element_t) * atom_count);
    
    // Set up atom types and residue names
    // HOH water oxygen
    strncpy(mol.atom.type[0].buf, "O", sizeof(mol.atom.type[0].buf));
    mol.atom.type[0].len = 1;
    strncpy(mol.atom.resname[0].buf, "HOH", sizeof(mol.atom.resname[0].buf));
    mol.atom.resname[0].len = 3;
    mol.atom.element[0] = 0; // Start unknown
    
    // HOH water hydrogen
    strncpy(mol.atom.type[1].buf, "H1", sizeof(mol.atom.type[1].buf));
    mol.atom.type[1].len = 2;
    strncpy(mol.atom.resname[1].buf, "HOH", sizeof(mol.atom.resname[1].buf));
    mol.atom.resname[1].len = 3;
    mol.atom.element[1] = 0; // Start unknown
    
    // Alanine alpha carbon
    strncpy(mol.atom.type[2].buf, "CA", sizeof(mol.atom.type[2].buf));
    mol.atom.type[2].len = 2;
    strncpy(mol.atom.resname[2].buf, "ALA", sizeof(mol.atom.resname[2].buf));
    mol.atom.resname[2].len = 3;
    mol.atom.element[2] = 0; // Start unknown
    
    // Sodium ion
    strncpy(mol.atom.type[3].buf, "NA", sizeof(mol.atom.type[3].buf));
    mol.atom.type[3].len = 2;
    strncpy(mol.atom.resname[3].buf, "NA", sizeof(mol.atom.resname[3].buf));
    mol.atom.resname[3].len = 2;
    mol.atom.element[3] = 0; // Start unknown
    
    // Generic carbon
    strncpy(mol.atom.type[4].buf, "C1", sizeof(mol.atom.type[4].buf));
    mol.atom.type[4].len = 2;
    strncpy(mol.atom.resname[4].buf, "", sizeof(mol.atom.resname[4].buf));
    mol.atom.resname[4].len = 0;
    mol.atom.element[4] = 0; // Start unknown
    
    // Call the element guess function
    bool result = md_util_element_guess(mol.atom.element, atom_count, &mol);
    
    // Verify results
    EXPECT_TRUE(result);
    EXPECT_EQ(mol.atom.element[0], MD_Z_O);  // Water oxygen
    EXPECT_EQ(mol.atom.element[1], MD_Z_H);  // Water hydrogen
    EXPECT_EQ(mol.atom.element[2], MD_Z_C);  // Alanine alpha carbon (not calcium!)
    EXPECT_EQ(mol.atom.element[3], MD_Z_NA); // Sodium ion
    EXPECT_EQ(mol.atom.element[4], MD_Z_C);  // Generic carbon from C1
    
    md_vm_arena_destroy(alloc);
}