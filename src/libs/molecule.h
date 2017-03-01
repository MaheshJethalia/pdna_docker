/*
*	Info: Header file to create a class that abstracts a biomolecule 
*	Classes: Molecule
*	Global Variables: size
*	Functions: None
*/

#ifndef MOLECULE_H
#define MOLECULE_H

#include <vector>
#include <string>

/*
*	Global Variables :
*	1) int size : Stores the number of rows/columns for the 3D array used to represent the discretized 
*			  	  biomolecules vector val for the current input data of protein and dna biomolecules.
*				  Defined in 'src/libs/molecule.h'. Declared in 'src/main.cpp'. Used in 'src/molecule.cpp'
*
*	2) float resolution : 
*
*
*/
extern int size;
extern float resolution;

/*
*	Classes :
*	1) Molecule: Basic class to represent a biomolecule (protein or DNA) which contains discretized matrix
*			member variables: vector<int> matrix, int center_of_mass
*			member methods:	
*/
class Molecule{
public:
	std::vector<int> matrix;
	float center_of_mass;
	int no_of_atoms;

	Molecule();
	int CreateMatrix(std::vector<float>& X, std::vector<float>& Y, std::vector<float>& Z);
	int IsEmpty();
	int RotateMatrix(float alpha, float beta, float gamma);
};

/*
*	Functions :
*
*/
int Indexv(int x, int y, int z);
int* Indexm(int n);

#endif