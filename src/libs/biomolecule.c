#ifndef BIOMOLECULE_H
#define BIOMOLECULE_H

#include "config.h"
#include "utility.h"

/*
 * Enum: molecule_class => classify type of biomolecule as either 
 * protein (PRO) or nucleic acid(DNA or RNA)
 */
enum molecule_class { DNA, PRO, RNA}; 

/* struct: Biomolecule => stores the necessary information related to a particular 
 * biomolecule.
 * data-members:
 *      filepath => Absolute of relative path to the pdb file of the biomolecule
 *      natoms => Number of atoms in the biomolecule
 *      type => Type of the biomolecule
 *      pdb_id => Four character PDB ID of the biomolecule
 *      name => Name of the biomolecule
 *      x, y, z => x, y and z coordinates of the atoms as mentioned in pdb file
 */
typedef struct biomolecule {
    char* filepath;
    int natoms;
    enum molecule_class type;
    char* pdb_id;
    char* name;
    float* x, y, z;
} biomolecule;

void create_biomolecule(biomolecule* bm, char* filepath);
void process_coordinates(biomolecule* bm);
config* rotate_coordinates(biomolecule* bm, int alpha, int beta, int gamma);

#endif
