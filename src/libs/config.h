#ifndef CONFIG_H
#define CONFIG_H
#include <vector>
#include "molecule.h"


class RotationalAngle {
    public:
        int alpha, beta, gamma;
        RotationalAngle();
        RotationalAngle(const int, const int, const int);
};

class Configuration {
    public:
        int noatoms;
        RotationalAngle angle;
        std::vector<Coordinate> coords;       
        std::vector<atom_type> atom_types;

        Configuration();
        void createConfiguration(Biomolecule* m, const RotationalAngle ang);
};

#endif
