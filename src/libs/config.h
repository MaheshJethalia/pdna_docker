#ifndef CONFIG_H
#define CONFIG_H
#include <vector>
#include "molecule.h"

struct intCoordinate { int x; int y; int z;};

class RotationalAngle {
    public:
        int alpha, beta, gamma;
        RotationalAngle();
        RotationalAngle(int, int, int);
};

class Configuration {
    public:
        int noatoms;
        RotationalAngle angle;
        std::vector<intCoordinate> icoords;       
        std::vector<atom_type> atom_types;

        Configuration();
        void createConfiguration(Biomolecule* m, const RotationalAngle ang);
};

#endif
