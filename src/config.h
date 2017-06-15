#ifndef CONFIG_H
#define CONFIG_H

/************** CUSTOMIZABLE VALUES **************/

// Sets the appropriate size of the grid for geometric scoring
#define GRID_SIZE 200

// Sets the steps in rotation around each axis of the DNA molecule
#define ROTATION_STEP 10

// Sets the distance covered per unit index of the space discretized matrix for geometric scoring
// NOTE: it is necessary that diameter(biomolecule) < GRID_SIZE * RESOLUTION else contraints aren't satisfied
#define RESOLUTION 1

// Defines the thickness of the protein surface for soft docking
#define SURFACE_THICKNESS 1.5

// Defines the penalty value per unit cell for geomteric scoring
#define RHO -10

// Sets flag for soft docking (1 => ON  | 0 => OFF)
#define SOFT_DOCK 1

// Sets the paths to the various input files
#define PROTEIN_PDB_FILEPATH ""
#define DNA_PDB_FILEPATH ""
#define OUTPUT_FILEPATH ""



/********** NOT TO BE CUSTOMIZED *************/
#define min(a, b) (((a)<(b))? (a) : (b))
#define max(a, b) (((a)>(b))? (a) : (b))
#define MEMORY_ERROR printf("\n******* ENCOUNTERED MEMORY ERROR: NOT ENOUGH MEMORY ********\n******* STOPPING *******\nERROR AT: ")
#define FILE_ERROR printf("\n******* ENCOUNTERED FILE ERROR: COULD NOT OPEN OR CREATE FILE *******\n********* STOPPING ********\nERROR AT: ")
#define length(x1, y1, z1, x2, y2, z2) sqrt(((x1)-(x2))*((x1)-(x2)) + ((y1)-(y2))*((y1)-(y2)) + ((z1)-(z2))*((z1)-(z2)))
#define index(x, y, z) ((z) + GRID_SIZE * ((y) + GRID_SIZE * (x)))
#define REAL 0
#define IM 1
#define PROTEIN_TMP_FILEPATH ""
#define DNA_TMP_FILEPATH ""
/********** NOT TO BE CUSTOMIZED *************/




#endif
