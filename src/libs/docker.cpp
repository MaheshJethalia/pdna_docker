#include "docker.h"

using namespace std;

/*
 * Default constructor for class Solution
 */
Solution::Solution(){ 
    translation = Coordinate(0, 0, 0);
    rotation = RotationalAngle(0, 0, 0);
    geometricScore = 0;
    RMS = 0;
}

/*
 * Parameterized constructor for class Solution
 */
Solution::Solution(const Coordinate c, const RotationalAngle r, const float gs){
    translation = c;
    rotation = r;
    geometricScore = gs;
}

/*
 * Default constructor for class Docker 
 */
Docker::Docker(){
    p = d = 0;
    nodecoys =0;
}

/*
 * Parameterised constructor for class Docker that creates the biomolecules *p and *d using the pdb files passed
 * as p_filename andd_filename via the parameters and sets the output file name to outfn
 */
Docker::Docker(const std::string p_filename,  const std::string d_filename, const int nd, std::string outfn){
    p = new Biomolecule(p_filename, PRO);
    d = new Biomolecule(d_filename, DNA);
    outfilename = outfn;
    nodecoys = nd;
}

/*
 * Main docking function
 */
void Docker::Dock(){
}
