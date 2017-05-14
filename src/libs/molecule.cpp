#include "molecule.h"
#include <vector>
#include <algorithm>

using namespace std;


/*
*	Implementation of class Molecule and other functions defined in molecule.h
*	NOTE: size, resolution are extern declared in molecule.h and defined in ../main.c
*/


/* 
 * Class Molecule Constructor: creates a 3d matrix as per the value of size set in main.cpp
 * and intializes it with the value 0
*/ 

Molecule::Molecule() {
	matrix.resize(size * size * size, 0);
	center_index = 0;
    xtop = xbot = ytop = ybot = ztop = zbot = 0;
    empty = 1;
}

/* 
 * Function to create a matrix of the macromolecule from the x, y and z coordiantes 
 * passed as floating vectors in the parameters
*/

int Molecule::CreateMatrix(vector<float>& X, vector<float>& Y, vector<float>& Z, int rho) {
	vector<int> vx, vy, vz;		// integer form of matrices X, Y, Z
	float xc, yc, zc;			// to hold index of geometric center of molecule
	xc = yc = zc = 0.0;

	for(int i = 0; i < (int)X.size(); i++) {			// creating vx, vy, vz from X, Y, Z
		xc += X[i]; yc += Y[i]; zc += Z[i];		// calculating sum of indices in each axis
		vx.push_back((int) X[i]/resolution);
		vy.push_back((int) Y[i]/resolution);
		vz.push_back((int) Z[i]/resolution);
	}	

	xc /= X.size(); yc /= Y.size(); zc /= Z.size();	// calculating avg of sum of indices = center
	center_index = index((int)xc, (int)yc, (int)zc);// setting value of center_index
    empty = 0;                                      // setting value of empty to false

	int xmin = *min_element(vx.begin(), vx.end());	// calculating the minimum index to 
	int ymin = *min_element(vy.begin(), vy.end());	// set all indices to positive by 
	int zmin = *min_element(vz.begin(), vz.end());	// adding the offset
	xmin = (xmin < 0) ? xmin : 0;
	ymin = (ymin < 0) ? ymin : 0;
	zmin = (zmin < 0) ? zmin : 0;

	xmin *= -1; ymin *= -1; zmin *= -1;

	for(int i = 0; i < (int) X.size(); i++) {		// adding the offset index to all values in vx, vy , vz
		vx[i] += xmin; 
		vy[i] += ymin;
		vz[i] += zmin;
	}
    xbot = ybot = zbot = xtop = ytop = ztop = 0;

	for(int i = 0; i < (int)X.size(); i++) {		// Creating the basic matrix
		matrix[index(vx[i], vy[i], vz[i])] = rho;	// set value to rho supplied as a parameter
        xtop = (xtop < vx[i]) ? vx[i] : xtop;       // checking and setting the values of the top and bottom coords
        ytop = (ytop < vy[i]) ? vy[i] : ytop;
        ztop = (ztop < vz[i]) ? vz[i] : ztop;
        xbot = (xbot > vx[i]) ? vx[i] : xbot;
        ybot = (ybot > vy[i]) ? vy[i] : ybot;
        zbot = (zbot > vz[i]) ? vz[i] : zbot;

	}

    return 0;
}

/* 
 * Function to shift the contents of the matrix of the molecule such that the geometric center
 * lies at the center of the matrix -> Index(size/2, size/2, size/2)
 */
int Molecule::CenterMatrix() {
    int xoff, yoff, zoff;                               // variables to store the offset for the geometric center

    int zc = center_index % size;                       // calculating the individual coords of geometric center from  center_index
    int yc = (center_index/size) % size;
    int xc = ((center_index/size)/size);

    xoff = (size/2) - xc; yoff = (size/2) - yc; zoff = (size/2) - zc; // setting values of offset variables

    for(int x = xtop; x >= xbot; x--) {                 // for loop to apply offset to each coordinate of the molecule
        for(int y = ytop; y >= ybot; y--) {
            for(int z = ztop; z >= zbot; z--) {
                matrix[index(x + xoff, y + yoff, z + zoff)] = matrix[index(x, y, z)];
                matrix[index(x, y, z)] = 0;
            }
        }
    }
    xbot += xoff; ybot += yoff; zbot += zoff;           // apply offset to coordinate data members which have shifted due to this
    xtop += xoff; ytop += yoff; ztop += zoff;

	return 0;
}

/*
 * Function to generate a surface for molecule 
 */
int Molecule::CreateSurface() {
    for(int x = xtop; x >= xbot; x--) {             // for loop to iterate through all occupied coordinates 
        for(int y = ytop; y >= ybot; y--) {
            for(int z = ztop; z >= zbot; z--) {
                if(is_surface_element(x, y, z, this))  // checking if current coordinate has a surface element
                    matrix[index(x, y, z)] = 1;     // setting value of surface element to 1
            }
        }
    }
	return 0;
}

/*
 * Function to check if the molecule is empty or not
 */
int Molecule::IsEmpty() {
	if(empty == 0)
		return 0;
	else
		return 1;
}

/* Function to return value of matrix(index(x, y, z)) in a smart manner, where missing coords return the value 0 instead of 
 * throwing an IndexOutOfBounds error; used in calculating the score for geometric complementarity.
 */
int Molecule::GetVal(int x, int y, int z) {
	if((x < 0 || x > size) || (y < 0 || y > size) || (z < 0 || z > size)) //check if indices are out of bounds
		return 0;
	else 
		return matrix[index(x, y, z)];
}

/*
 * Function to check if the element at the provided coords is an index element; used by molecule::CreateSurface()
 */
int is_surface_element(int x, int y, int z, const Molecule* m) {
    // Condn: top most coord elements will always be surface elements.
    if(x == m->xtop || x == m->xbot || y == m->ytop || y == m->ybot || z == m->ztop || z == m->zbot)
        return 1;
    // Condn: For all other elements, check if neighbouring coords are occupied or not. If any of them are not occupied, then
    // current element is on the surface of the molecule.
    for(int i = -1; i<= 1; i++) {
        for(int j = -1; j<=1; j++){
            for(int k = -1; k<= 1; k++){
                if(i == 0 && j ==0 && k == 0) continue;     // skip the coord of the element itself
                if(m->matrix[index(x + i, y + j, z + k)] == 0)         
                    return 1;
            }
        }
    }
    return 0;
}


