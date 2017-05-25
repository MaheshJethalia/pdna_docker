#ifndef MOLECULE_H
#define MOLECULE_H

#include "config.h"
#include <string>
#include <vector>

enum biomol_type { PRO, DNA };

class Coordinate {
    public:
        int x, y, z;

        Coordinate();
        Coordinate(int, int, int);
        Coordinate Rotate(RotationalAngle&);
};

class Atom {
    public:
        Coordinate coord;
        biomol_type type;
        std::string name;

        Atom();
        Atom(Coordinate&, biomol_type, std::string&);
};

class Biomolecule {
    public:
        int noatoms;
        biomol_type type;
        Atom* atoms;
        std::string PDB;
};

#endif
