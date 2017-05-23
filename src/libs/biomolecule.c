#include "biomolecule.h"
#include <stdlib.h>
#include <stdio.h>

void create_biomolecule(biomolecule* bm, char* fp) {
    bm = (biomolecule*)malloc(sizeof(biomolecule)); 
    bm->filepath = fp;
}
