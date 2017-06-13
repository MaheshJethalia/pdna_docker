#include "docker.h"


float get_biomolecule_diameter(Biomolecule* m) {
    // Local variables
    float x1, x2, y1, y2, z1, z2; 
    int i;

    // Intitialising the local variables
    x1 = x2 = m->atoms[0].coordinate.x;
    y1 = y2 = m->atoms[0].coordinate.y;
    z1 = z2 = m->atoms[0].coordinate.z;
    
    // Calculating the coordinates of the extremes
    for(i=1; i < m->number_of_atoms; i++) {
        x1 = min(x1, m->atoms[i].coordinate.x);
        x2 = max(x2, m->atoms[i].coordinate.x);
        y1 = min(y1, m->atoms[i].coordinate.y);
        y2 = max(y2, m->atoms[i].coordinate.y);
        z1 = min(z1, m->atoms[i].coordinate.z);
        z2 = max(z2, m->atoms[i].coordinate.z);
    }

    // length from min coords to max coords ~ diameter of molecule
    return length(x1, y1, z1, x2, y2, z2);
}

int check_constraints(Biomolecule* A, Biomolecule* B) {

    // Checking that rotation step produces integer angles of rotation
    if(360 % ROTATION_STEP != 0) {
        printf("\n*** CONFIG ERROR: ROTATION_STEP must be an integer factor of 360 degrees ***");
        printf("\n*** ERROR IN=> Function: check_constraints | docker.c | 29 ***");
        return 0;
    }

    // Checking if rotation steps are sufficiently bounded or not
    if(!(ROTATION_STEP >= 5 && ROTATION_STEP <= 90)) {
        printf("\n*** CONFIG ERROR: ROTATION_STEP out of allowed range [5, 90] ***;");
        printf("\n*** ERROR IN=> Function: check_constraints | docker.c | 35 ***");
        return 0;
    }

    // Checking if biomolecule fits into specified grid or not
    if((get_biomolecule_diameter(A) > RESOLUTION * GRID_SIZE ) || (get_biomolecule_diameter(B) > RESOLUTION * GRID_SIZE)) {
        printf("\n*** CONFIG ERROR: Biomolecule not fitting in Grid size ***");
        printf("\n*** ERROR IN => Function: check_constraints | docker.c | 41 ***");
        return 0;
    }

    // All OK, return Good
    return 1;
}

int create_configuration(Configuration* config, const AngleOfRotation tangle, Biomolecule* tparent) {
    // Local variable declarations
    int flag = 1; 
    int i;

    // Allocating memory to config and handling memory errors
    if((config = (Configuration*)malloc(sizeof(Configuration))) == NULL) {
        MEMORY_ERROR;
        printf("\nFunction: create_configuration | docker.c | 60");

        // Cleaning up function heap allocations:
        free(config);

        // Return status flag
        flag = 0; return flag;
    }

    // If allocations are fine, continue creating the configuration
    config->parent = tparent;
    config->angle = tangle;
    config->number_of_atoms = tparent->number_of_atoms;
    config->center = tparent->center;

    // Allocating memory to config->atoms and handling memory errors
    if((config->atoms = (Atom*)malloc(sizeof(Atom) * config->number_of_atoms)) == NULL){
        MEMORY_ERROR;
        printf("\nFunction: create_configuration | docker.c | 72");

        // cleaning up function heap allocations
        free(config);
        
        return 0;
    }
    
    // All OK, copying atoms from config->parent to config
    for(i = 0; i < config->number_of_atoms; i++){
        config->atoms[i] = tparent->atoms[i];
        config->atoms[i].coordinate = rotate_coordinate(config->atoms[i].coordinate, tangle, config->center);
    }
    return flag;
}

void create_geometric_core(Configuration* config, Matrix* tmp, int init_value) {
    // Local variable declaration
    int i, j, k;

    tmp->parent_config = config;
    // Initialize tmp->val, then create the biomolecule core matrix
    for(i = 0; i < GRID_SIZE; i++) {
        for(j = 0; j < GRID_SIZE; j++) {
            for(k = 0; k < GRID_SIZE; k++) {
                tmp->value[index(i, j, k)] = 0.0;   // Initialising all values to 0.0
            }
        }
    }
    
    tmp->xmax = tmp->ymax = tmp->zmax = tmp->xmin = tmp->ymin = tmp->zmin = 0;
    for( i = 0; i < tmp->parent_config->number_of_atoms; i++) {
        int tx = (int) tmp->parent_config->atoms[i].coordinate.x / RESOLUTION;
        int ty = (int) tmp->parent_config->atoms[i].coordinate.y / RESOLUTION;
        int tz = (int) tmp->parent_config->atoms[i].coordinate.z / RESOLUTION;
        tmp->value[index(tx, ty, tz)] = init_value;        // Setting all core locations to RHO
        tmp->xmax = max(tmp->xmax, tx);             // Finding the bounds for core locations 
        tmp->ymax = max(tmp->ymax, ty);             // in the matrix value
        tmp->zmax = max(tmp->zmax, tz);             // for quicker lookup in the future
        tmp->xmin = min(tmp->xmin, tx);
        tmp->ymin = min(tmp->ymin, ty);
        tmp->zmin = min(tmp->zmin, tz);
    }

    
}

void create_geometric_surface(Matrix* matrix) {
    // Local variable declarations
    int i, j, k, p, q, r, flag;

    // Search through all core locations and set all normal surface locations to 1
    for( i = matrix->xmin; i <= matrix->xmax; i++) {
        for(j = matrix->ymin; j <= matrix->ymax; j++) {
            for(k = matrix->zmin; k <= matrix->zmax; k++) {
                if(matrix->value[index(i, j, k)] == RHO) {     // check if (i, j, k) stores a core == RHO
                    flag = 0;
                    for(p = -1; p <= +1; p++) {                 //check all neighbours of (i, j, k) for empty space
                        for(q = -1; q <= +1;q++) {
                            for(r = -1; r <= +1; r++){
                                if(i+p > 0 && j+q > 0 && k+r > 0 && i+p < GRID_SIZE && j+q < GRID_SIZE && k+r < GRID_SIZE && matrix->value[index(i+p, j+q, k+r)] == 0.0) 
                                    flag = 1;   // if empty neighbour exists, then current location is on surface 
                            }
                        }
                    }
                    if(flag == 1)
                        matrix->value[index(i, j, k)] = 1.0;    // set to surface value 1.0
                }
            }
        }
    } 
}

