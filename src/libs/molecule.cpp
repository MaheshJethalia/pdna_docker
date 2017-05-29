#include "molecule.h"
#include <cmath>

using namespace std;

Coordinate::Coordinate() {
    x = 0; y = 0; z = 0;
}

Coordinate::Coordinate(int tx, int ty, int tz){
    x = tx; y = ty; z = tz;
}

Coordinate Coordinate::Rotate(RotationalAngle& angle, Coordinate center){
    Coordinate transformed;
    // performing rotation about z axis (gamma) first
    transformed.x = center.x + (x - center.x)*cos(angle.gamma) - (y - center.y)*sin(angle.gamma);
    transformed.y = center.y + (x - center.x)*sin(angle.gamma) + (y - center.y)*cos(angle.gamma);
    transformed.z = z;
    //performing rotation about y axis (beta) second on the previously obtained coords
    return transformed;
}
