#ifndef UTILITIES_H
#define UTILITIES_H
#include "molecule.h"
#include <string>

/*
 * class name: Result
 * data members:
 *      <x, y, z> trans: the translational coordinate offsets for the current result
 *      <alpha, beta, gamma>: the rotational coordinate offsets for the current result
 *      GC, ES, VDWS, TS: Geometric Complementarity, Electrostatic Score, Vander Waal Score and Total Score for the Result
 * member functions:
 *      GetString(): Return the result in a formatted string which can be pretty printed to stdout for written in a file
 */
class Result {
    public:
        int xtrans, ytrans, ztrans;
        int  alpha, beta, gamma;
        float GC, EI, VDWI;
        float TS;
        Result();
        Result(int x, int y, int z, int a, int b, int g, float gc, float ei, float vdwi, float ts);
        std::string GetString();
};

/*
 * Functions:
 *  geometric_complementarity(): obtain the GC score for the provided molecules
 *  electrostatic_score(): obtain the ES score for the provided molecules
 *  vanderwaal_score(): obtain the VDWS score for the provided molecules
 *  scoring_function(): obtain the TS of the provided molecules
 */

float geometric_complementarity(const Molecule* p, const Molecule* d);
float electrostatic_score(const Molecule* p, const Molecule* d); 
float vanderwaal_score(const Molecule* p, const Molecule* d); 
float scoring_function(const Molecule* p, const Molecule* d); 


#endif
