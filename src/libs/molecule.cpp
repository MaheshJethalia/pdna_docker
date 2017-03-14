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

Molecule::Molecule() {
	matrix.resize(size * size * size, 0);
	center_index = 0;
    no_of_atoms = 0;
}

/* Function to create a matrix of the macromolecule from the x, y and z xoordiantes 
 * passed as floating vectors in the parameters*/

int Molecule::CreateMatrix(vector<float>& X, vector<float>& Y, vector<float>& Z, int rho) {
	vector<int> vx, vy, vz;		// integer form of matrices X, Y, Z
	float xc, yc, zc;			// to hold index of geometric center of molecule
	xc = yc = zc = 0.0;

	for(int i = 0; i < X.size(); i++) {			// creating vx, vy, vz from X, Y, Z
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

	for(int i = 0; i < X.size(); i++) {		// adding the offset index to all values in vx, vy , vz
		vx[i] += xmin; 
		vy[i] += ymin;
		vz[i] += zmin;
	}

	for(int i = 0; i < X.size(); i++) {		// Creating the basic matrix
		matrix[index(vx[i], vy[i], vz[i])] = rho;	// -10 = value for core atom of protein
	}

    return 0;
}


int Molecule::CenterMatrix() {
	return 0;
}

int Molecule::CreateSurface() {
	return 0;
}

int Molecule::Rotate(int alpha, int beta, int gamma) {
	return 0;
}

int Molecule::IsEmpty() {
	if(no_of_atoms == 0)
		return 1;
	else
		return 0;
}

int Molecule::GetVal(int x, int y, int z) {
	if((x < 0 || x > size) || (y < 0 || y > size) || (z < 0 || z > size))
		return 0;
	else 
		return matrix[index(x, y, z)];
}

int index(int x, int y, int z) {
	return z + y * size + x * size * size;
}

