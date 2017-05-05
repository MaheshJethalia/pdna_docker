#ifndef UTILITIES_H
#define UTILITIES_H
#include "molecule.h"
class Result {
    public:
        int xtrans, ytrans, ztrans;
        float alpha, beta, gamma;
        float GC, EI, VDWI;
        float TS;
        Result();
        Result(int x, int y, int z, float a, float b, float g, float gc, float ei, float vdwi, float ts);
        void PrintTabular();
};

float geometric_complementarity(const Molecule* p, const Molecule* d);
float electrostatic_interaction(const Molecule* p, const Molecule* d); 
float vanderwaal_interaction(const Molecule* p, const Molecule* d); 
float scoring_function(const Molecule* p, const Molecule* d); 


#endif
