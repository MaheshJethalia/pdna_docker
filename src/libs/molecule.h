/*
*	Info: Header file to create a class that abstracts a biomolecule 
*	Classes: Molecule
*	Global Variables: size, resolution
*	Global Functions: index(), index_return() 
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
*	2) float resolution : A measure of the length in angstrom of each unit in the discretized structure of
*				  the molecule. For instance, if resolution = 1.9, then this means that the  distance between 
*				  self->matrix[Index(1, 0, 0)] and self->Matrix[Index(2, 0, 0)] is 1.9 Angstrom
*
*
*/
extern int size;
extern float resolution;

/*
*	Classes :
*	1) Molecule: Basic class to represent a biomolecule (protein or DNA) which contains discretized matrix
*
*			member variables: vector<int> matrix => store the discretized structure of the molecule, storing
*                      		  		1 in (i, j, k) if an atom of the molecule exists there and 0 if none exist.	
*				   			  int center_index => contains the index (i, j, k) of the geometric center of the 
*									molecule as it is created using the atomic coordinate matrices
*
*			member methods:	Molecule() => standard default constructor
*							CreateMatrix() => creates the discretized matrix for the biomolecule based on the 
*									atomic coordinate vectors it takes as parameters
*							IsEmpty() => check if the Molecule is empty, i.e. contains zero atoms
*							GetVal() => smartly get the value of self->matrix(x, y, z), handles out of bounds smartly 
*									checking with index limits and returning zero if it is out of bounds
*							
*/
class Molecule{
public:
	std::vector<int> matrix;
	int center_index;

	Molecule();
	int CreateMatrix(std::vector<float>& X, std::vector<float>& Y, std::vector<float>& Z);
	int IsEmpty();
	int GetVal(int x, int y, int z);

/*
*	Functions : index(x, y, z): returns the 1 D index of the 3 D matrix representation passed on as parameters
*				index_return(index): returns the coordinates {x, y, z} which corresponds to 1 D matrix as index
*
*/
int index(int x, int y, int z);
int* index_return(int n);

#endif