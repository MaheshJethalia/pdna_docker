#include "config.h"

using namespace std;

// Default constructor for class RotationalAngle
RotationalAngle::RotationalAngle() {
    alpha = beta = gamma = 0;
}

// Parameterized constructor for class RotationalAngle
RotationalAngle::RotationalAngle( const int a, const int b, const int g){
    alpha = a; beta = b; gamma = g;
}

// Default constructor for class Configuration
Configuration::Configuration(){
    noatoms = 0;
}

// Function to generate a new configuration of Biomolecule m rotated by angles ang about the center
void Configuration::createConfiguration(Biomolecule* m, const RotationalAngle ang){
    noatoms = m->noatoms;       // set the no of atoms of config to that of parent biomolecule
    angle = ang;                // set the angle of rotation of the config to that supplied as param

    for(vector<Atom>::iterator it = m->atoms.begin(); it!= m->atoms.end(); it++) { // iterate through 
                                // all atoms in biomolecule m
        coords.push_back(it->coord.Rotate(ang, m->center));     // add the rotated coord of the atom to
                                // coords
    }
}
