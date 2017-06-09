#ifndef MOLECULE_H
#define MOLECULE_H

#include "config.h"
#include <string>
#include <vector>

enum biomol_type { NDEF, PRO, DNA };
enum molecule_type { NDEF, C, N};

class Coordinate {
    public:
        float x, y, z;

        Coordinate();
        Coordinate(float, float, float);
        Coordinate Rotate(RotationalAngle&, Coordinate center);
};

class Atom {
    public:
        Coordinate coord;
        molecule_type moltype;

        Atom();
        Atom(Coordinate&, molecule_type moltype);
};

class Biomolecule {
    public:
        int noatoms;
        biomol_type type;
        std::vector<Atom> atoms;

};
#endif
