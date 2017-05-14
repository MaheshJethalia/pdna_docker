/*
*	Info: Header file to create a class that abstracts a biomolecule 
*/

#ifndef MOLECULE_H
#define MOLECULE_H

#include <vector>
#include <string>

// defining macros -> index(i, j, k)
#define index(x, y, z) ((z) + (y) * size + (x) * size * size)

/*
*	Global Variables (defined in main.cpp):
*/
extern int size; // Stores the value of the size of the matrices of the biomolecules to be docked.
extern float resolution; // Stores the value of the actual distance b/w 2 adjacent grid points in molecule::matrix

/*
*	 Class: Molecule
*	 Member properties:
*	    int empty: flag variable to check if molecule is properly initialized or not
*	    vector<int> matrix: represents the 3-d structure of the biomolecule
*	    int center_index: represents the index of the geometric center of the biomolecule
*	    int <x,y,z><top, bot>: represents the coords of the <top, bottom> atoms of the biomolecule in the <x, y, z> axes 
*	 Member functions:
*	    CreateMatrix(vector, vector, vector): creates molecule::matrix from x, y and z coords of atoms supplied to function as 3 vectors
*	    IsEmpty(): 
*	    GetVal(int, int, int):
*	    CenterMatrix():
*	    CreateSurface():
*
*/

class Molecule{
public:
    // defining data members:
    int empty;
	std::vector<int> matrix;
	int center_index;
    int xtop, ytop, ztop, xbot, ybot, zbot;

	Molecule(); // basic constructor

    // defining member functions:
	int CreateMatrix(std::vector<float>& X, std::vector<float>& Y, std::vector<float>& Z, int rho);
	int IsEmpty();
	int GetVal(int x, int y, int z);
	int CenterMatrix();
	int CreateSurface();
};

/*
*	Functions : 
*	    is_surface_element(int, int, int): Checks if element at the indices provided as an argument is a surface element or not
*/
int is_surface_element(int x, int y, int z, const Molecule* m);
#endif
