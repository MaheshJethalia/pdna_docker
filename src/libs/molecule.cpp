#include "molecule.h"
#include <vector>
#include <algorithm>
using namespace std;


/*
*	Implementation of class Molecule and other functions defined in molecule.h
*	NOTE: size is extern declared in molecule.h and defined in ../main.c
*/


/* Class Molecule Constructor: creates a 3d matrix as per the value of size set in main.cpp
 * and intializes it with the value 0*/

GCMatrix::GCMatrix() {
	matrix.resize(size * size * size, 0);
	center_index = 0;
    no_of_atoms = 0;
    xtop = xbot = ytop = ybot = ztop = zbot = 0;
}

/* Function to create a matrix of the macromolecule from the x, y and z xoordiantes 
 * passed as floating vectors in the parameters*/

int GCMatrix::CreateMatrix(vector<float>& X, vector<float>& Y, vector<float>& Z, int rho) {
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
	no_of_atoms = X.size();							// setting no_of_index

	int xmin = *min_element(vx.begin(), vx.end());	// calculating the minimum index to 
	int ymin = *min_element(vy.begin(), vy.end());	// set all indices to positive by 
	int zmin = *min_element(vz.begin(), vz.end());	// adding the offset
	xmin = (xmin < 0) ? xmin : 0;
	ymin = (ymin < 0) ? ymin : 0;
	zmin = (zmin < 0) ? zmin : 0;

	xmin *= -1; ymin *= -1; zmin *= -1;

	for(int i = 0; i < (int)X.size(); i++) {		// adding the offset index to all values in vx, vy , vz
		vx[i] += xmin; 
		vy[i] += ymin;
		vz[i] += zmin;
	}
    xbot = ybot = zbot = xtop = ytop = ztop = 0;

	for(int i = 0; i < (int)X.size(); i++) {		// Creating the basic matrix
		matrix[index(vx[i], vy[i], vz[i])] = rho;	// set value to rho supplied as a parameter
        xtop = (xtop < vx[i]) ? vx[i] : xtop;
        ytop = (ytop < vy[i]) ? vy[i] : ytop;
        ztop = (ztop < vz[i]) ? vz[i] : ztop;
        xbot = (xbot > vx[i]) ? vx[i] : xbot;
        ybot = (ybot > vy[i]) ? vy[i] : ybot;
        zbot = (zbot > vz[i]) ? vz[i] : zbot;

	}

    return 0;
}

/* Function to shift the contents of the matrix of the molecule such that the geometric center
 * lies at the center of the matrix -> Index(size/2, size/2, size/2)
 */
int GCMatrix::CenterMatrix() {
    int xoff, yoff, zoff;
    int zc = center_index % size;
    int yc = (center_index/size) % size;
    int xc = ((center_index/size)/size);
    xoff = (size/2) - xc; yoff = (size/2) - yc; zoff = (size/2) - zc;
    for(int x = xtop; x >= xbot; x--) {
        for(int y = ytop; y >= ybot; y--) {
            for(int z = ztop; z >= zbot; z--) {
                matrix[index(x + xoff, y + yoff, z + zoff)] = matrix[index(x, y, z)];
                matrix[index(x, y, z)] = 0;
            }
        }
    }
    xbot += xoff; ybot += yoff; zbot += zoff;
    xtop += xoff; ytop += yoff; ztop += zoff;

	return 0;
}

int GCMatrix::CreateSurface() {
    for(int x = xtop; x >= xbot; x--) {
        for(int y = ytop; y >= ybot; y--) {
            for(int z = ztop; z >= zbot; z--) {
                if(is_surface_element(x, y, z, this))
                    matrix[index(x, y, z)] = 1;
            }
        }
    }
	return 0;
}

int GCMatrix::IsEmpty() {
	if(no_of_atoms == 0)
		return 1;
	else
		return 0;
}

int GCMatrix::GetVal(int x, int y, int z) {
	if((x < 0 || x > size) || (y < 0 || y > size) || (z < 0 || z > size))
		return 0;
	else 
		return matrix[index(x, y, z)]; }

int index(int x, int y, int z) {
	return z + y * size + x * size * size;
}

int is_surface_element(int x, int y, int z, const GCMatrix* m) {
    if(x == m->xtop || x == m->xbot || y == m->ytop || y == m->ybot || z == m->ztop || z == m->zbot)
        return 1;
    for(int i = -1; i <= 1; i++){
        for(int j = -1; j<= 1; j++) {
            for(int k = -1; k <= 1; k++) {
                if(i==0 && j==0 && k==0)
                    continue;
                else if(m->matrix[index(x+i, y+j, z+k)] == 0)
                    return 1;
            }
        }
    }
    return 0;
}


