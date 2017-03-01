#include "../header/molecule.h"

Molecule::Molecule() {
	matrix(vector<int>(size * size * size, 0));
	no_of_atoms = 0;
	center_of_mass = 0.0;
}


int Molecule::CreateMatrix(vector<float>& X, vector<float>& Y, vector<float>& Z) {

}


int RotateMatrix(float alpha, float beta, float gamma) {

}


int Molecule::IsEmpty() {
	if(no_of_atoms == 0)
		return 1
	else
		return 0;
}


int Indexv(int x, int y, int z) {
	return x + y * size + z * size * size;
}

int* Indexm(int n) {

}