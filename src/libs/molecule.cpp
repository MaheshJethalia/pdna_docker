#include "molecule.h"
#include <vector>
#include <algorithm>
using namespace std;


/*
*	Implementation of class Molecule and other functions defined in molecule.h
*	NOTE: size is extern declared in molecule.h and defined in ../main.c
*/

Molecule::Molecule() {
	matrix.resize(size * size * size, 0);
	center_index = 0;
    no_of_atoms = 0;
}


int Molecule::CreateMatrix(vector<float>& X, vector<float>& Y, vector<float>& Z) {
	vector<int> vx, vy, vz;
	float xc, yc, zc;
	xc = yc = zc = 0.0;

	for(int i = 0; i < X.size(); i++) {
		xc += X[i]; yc += Y[i]; zc += Z[i];
		vx.push_back((int) X[i]/resolution);
		vy.push_back((int) Y[i]/resolution);
		vz.push_back((int) Z[i]/resolution);
	}	

	xc /= X.size(); yc /= Y.size(); zc /= Z.size();

	int xmin = min(0, min(vx.begin(), vx.end()));
	int ymin = min(0, min(vy.begin(), vy.end()));
	int zmin = min(0, min(vz.begin(), vz.end()));

	xmin *= -1; ymin *= -1; zmin *= -1;

	for(int i = 0; i < X.size(); i++) {
		vx[i] += xmin; 
		vy[i] += ymin;
		vz[i] += zmin;
	}

	for(int i = 0; i < X.size(); i++) {
		matrix[index(vx[i], vy[i], vz[i])] = 1;

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

