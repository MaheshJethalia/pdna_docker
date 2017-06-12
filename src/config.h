#ifndef CONFIG_H
#define CONFIG_H

/*** CUSTOMIZABLE VALUES ***/
#define GRID_SIZE 200
#define ROTATION_STEP 10
#define RESOLUTION 1
#define SURFACE_THICKNESS 1.5
#define NUMBER_OF_DECOYS 1000
#define RHO -10
#define SOFT_DOCK 1

/*** NOT TO BE CUSTOMIZED ***/
#define min(a, b) (((a)<(b))? (a) : (b))
#define max(a, b) (((a)>(b))? (a) : (b))
#define MEMORY_ERROR printf("\n******* ENCOUNTERED MEMORY ERROR: NOT ENOUGH MEMORY ********\n******* STOPPING *******\nERROR AT: ")
#define FILE_ERROR printf("\n******* ENCOUNTERED FILE ERROR: COULD NOT OPEN OR CREATE FILE *******\n********* STOPPING ********\nERROR AT: ")
#define length(x1, y1, z1, x2, y2, z2) sqrt(((x1)-(x2))*((x1)-(x2)) + ((y1)-(y2))*((y1)-(y2)) + ((z1)-(z2))*((z1)-(z2)))
#define index(x, y, z) ((z) + GRID_SIZE * ((y) + GRID_SIZE * (x)))
#define REAL 0
#define IM 1


#endif
