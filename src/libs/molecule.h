#ifndef MOLECULE_H
#define MOLECULE_H

#include "config.h"
#include <string>
#include <vector>

enum biomol_type { PRO, DNA };
enum atom_type { C, N};

class Coordinate {
    public:
        float x, y, z;

        Coordinate();
        Coordinate(const float, const float, const float);
        Coordinate Rotate(RotationalAngle&, const Coordinate center);
};

class Atom {
    public:
        Coordinate coord;
        atom_type type;

        Atom();
        Atom( const Coordinate c ,const  atom_type t);
};

class Biomolecule {
    public:
        int noatoms;
        biomol_type type;
        std::vector<Atom> atoms;
        Coordinate center;

        Biomolecule();
        Biomolecule(const std::string, const biomol_type t);
        void FilterCoordinates();
};
#endif
