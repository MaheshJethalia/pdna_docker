#include "molecule.h"
#include <vector>

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

