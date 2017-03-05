/*
*	Info: Header file to create a class that abstracts a biomolecule 
*	Classes: Molecule
*	Global Variables: size, resolution
*	Global Functions: index(), get_center_index(), 
*/

#ifndef MOLECULE_H
#define MOLECULE_H

#include <vector>
#include <string>

/*
*	Global Variables :
*/
extern int size;
extern float resolution;

/*
*	 Class: Molecule
* 	 Structures: 
*	 member variables: vector<int> matrix, int center_index, int no_of_atoms
*	 member methods:	Molecule(), CreateMatrix(),	IsEmpty(), GetVal()
*/

class Molecule{
public:
	std::vector<int> matrix;
	int center_index;
    int no_of_atoms;

	Molecule();

	int CreateMatrix(std::vector<float>& X, std::vector<float>& Y, std::vector<float>& Z);
	int IsEmpty();
	int GetVal(int x, int y, int z);
	int CenterMatrix();
};

/*
*	Functions : index(x, y, z), 
*/
int index(int x, int y, int z);

#endif
