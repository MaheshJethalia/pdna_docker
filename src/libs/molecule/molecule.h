#ifndef MOLECULE_H
#define MOLECULE_H

#include <vector>
#include <string>

class Molecule {
public:
	int vector< vector< vector<int> > > matrix;
	Molecule(int n);
	Rotate(float alpha, float beta, float gamma);
	Create(string s);
};

#endif