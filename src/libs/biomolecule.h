#ifndef BIOMOLECULE_H
#define BIOMOLECULE_H

#include "utility.h"
/*
 * enum biomolecule_type => assign molecule types as either:
 *  Proteins => PRO     Deoxyribonucleic acid => DNA
 */
enum biomolecule_type { PRO, DNA, RNA};

// extern variables defined in main.c
extern float resolution;    
extern int size;

/*
 * struct biomolecule => store information regarding a biomol. being docked
 * data members :
 *      filepath    => absolute or relative path to coordinate file for the biomolecule
 *      pdb_id      => PDB Id of the biomolecule
 *      type        => type of biomolecule (PRO, DNA)
 *      x, y, z     => 3D coords of atoms of biomolecule as mentioned in pdb file
 *      no_atoms    => number of atoms in the biomolecule
 */
typedef struct biomolecule {
    char* filepath;
    char* pdb_id;
    enum biomolecule_type type;
    float *x, *y, *z;
    int no_atoms;
} biomolecule;

void create_biomolecule(biomolecule* bm, char* fp);

void process_coordinates(biomolecule* bm);

#endif
