#include "molecule.h"
#include <cmath>
#include <fstream>
using namespace std;

/*
* Default constructor for class Coordinate
 */
Coordinate::Coordinate() {
    x = 0; y = 0; z = 0;    //sets x, y, z coordinates to 0 by default
}

/*
 * Parameterised constructor for class Coordinate
 * Coordinate(xx, yy, zz) : sets the value of objects x, y and  coordinates to xx, yy, zz supplied
 * as parameters
 */
Coordinate::Coordinate(float tx, float ty, float tz){
    x = tx; y = ty; z = tz;
}

/*
 * Coordinate::Rotate(angle, center): Rotates the coordinates of the object about the [center] by the
 * [angles] about the specific axes x, y, z and return the value of the new coordinate. 
 * NOTE: This function does not change the original values of x,y,z of the object from which it was
 * called.
 */
Coordinate Coordinate::Rotate(RotationalAngle& angle, Coordinate center){
    Coordinate transformed;
    // performing rotation about z axis (gamma) first
    transformed.x = center.x + (float) ((x - center.x)*cos(angle.gamma) - (y - center.y)*sin(angle.gamma));
    transformed.y = center.y + (float)((x - center.x)*sin(angle.gamma) + (y - center.y)*cos(angle.gamma));
    transformed.z = z;

    //performing rotation about y axis (beta) second on the previously obtained coords
    transformed.y = center.y + (float)( (transformed.y - center.y)*cos(angle.alpha) - (transformed.z - center.z)*sin(angle.alpha) );
    transformed.z = center.z + (float)( (transformed.y - center.y)*sin(angle.alpha) + (transformed.z - center.z)*cos(angle.alpha) );

    //performing rotation about x axis (alpha) third on the previously obtained coords
    transformed.x = center.x + (float)( (transformed.z - center.z)*sin(angle.beta) + (transformed.x - center.x)*cos(angle.beta));
    transformed.z = center.z + (float)( (transformed.z - center.z)*cos(angle.beta) - (transformed.x - center.x)*sin(angle.beta));
    
    // returning the new rotated coordinate
    return transformed;
}


/*
 * Default constructor for class Atom
*/
Atom::Atom(){
    coord.x = coord.y = coord.z = 0;
    moltype = (molecule_type) 0;
}

/*
 * Parameterised constructor for class Atom:
 * Atom(Coordinate& c, molecule_type m) => sets the value of the coord of the atom to c and moltype to m
*/
Atom::Atom(Coordinate& c, molecule_type m){
    coord = c; moltype = m;
}

/*
 * Default constructor for class Biomolecule
 */
Biomolecule::Biomolecule(){
    noatoms = 0;
    type = NDEF;
}

/*
 * Parameterised constructor for class Biomolecule:
 * Biomolecule(string filename)=> reads the data from the filename and stores the coordinates from the same
 * in atoms[]
 */
Biomolecule::Biomolecule(string filename) {
    ifstream fin(filename);
    fin >> noatoms;
    float tx, ty, tz;
    Atom temp_atom;

    for(int i = 0; i < noatoms; i++){
        fin >> tx >> ty >> tz;
        temp_atom.coord.x = tx; temp_atom.coord.y = ty; temp_atom.coord.z = tz;
        atoms.push_back(temp_atom);
    }
    fin.close();
}