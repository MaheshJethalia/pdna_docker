#include "molecule.h"

/*
*	Implementation of class Molecule and other functions defined in molecule.h
*	NOTE: size is extern declared in molecule.h and defined in ../main.c
*/
Molecule::Molecule() {
	matrix(vector<int>(size * size * size, 0));
	no_of_atoms = 0;
	center_index = 0;
}


int Molecule::CreateMatrix(vector<float>& X, vector<float>& Y, vector<float>& Z) {
}


int Molecule::IsEmpty() {
	if(no_of_atoms == 0)
		return 1
	else
		return 0;
}

int Molecule::GetVal(int index) {

}
int Indexv(int x, int y, int z) {
	return x + y * size + z * size * size;
}

int* Indexm(int n) {

}