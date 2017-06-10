#ifndef DOCKER_H
#define DOCKER_H
#include "molecule.h"
#include "config.h"
#include <vector>

class Solution {
    public:
        Coordinate translation;
        RotationalAngle rotation;
        float geometricScore;
        float RMS;
        
        Solution();
        Solution(const Coordinate trans, const RotationalAngle rot, const float gs);
};

class Docker {
    public:
        Biomolecule* p;
        Biomolecule* d;
        int nodecoys;
        std::vector<Solution> result_stack;
        std::string outfilename;

        Docker();
        Docker(const std::string p_filename,const std::string d_filename, const int nd, std::string outfn);
        void Dock();
};

#endif
