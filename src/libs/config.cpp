#include "config.h"

// Define macros for easy traversal of the container

using namespace std;

RotationalAngle::RotationalAngle() {
    alpha = beta = gamma = 0;
}

RotationalAngle::RotationalAngle( int a, int b, int g){
    alpha = a; beta = b; gamma = g;
}


Configuration::Configuration(){
    noatoms = 0;
}

void Configuration::createConfiguration(Biomolecule* m, const RotationalAngle ang){
    noatoms = m->noatoms; 
    angle = ang;
    for(vector<Atom>::iterator it = m->atoms.begin(); it!= m->atoms.end(); it++) { 
    }
     
}
