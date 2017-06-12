#include "docker.h"

int generate_trig_tables(){
    int flag = 1;
    float* sin_table; float* cos_table;
    if(((sin_table = (float*)malloc(360 / ROTATION_STEP)) != NULL) && ((cos_table = (float*)malloc(360 / ROTATION_STEP)) != NULL)) {
        for(int i = 0; i < 360; i += ROTATION_STEP) {
            sin_table[i / ROTATION_STEP] = (float) sin((float)i);
            cos_table[i / ROTATION_STEP] = (float) cos((float)i);
        }
    }
    else flag = 0;
    return flag;
}

Coordinate create_coordinate(const float tx, const float ty, const float tz) {
    Coordinate tc; tc.x = tx, tc.y = ty; tc.z = tz;
    return tc;
}

Coordinate rotate_coordinate(Coordinate tc, const AngleOfRotation ta, const Coordinate tcenter){

    // Rotation about the X-axis (alpha)
    tc.y = tcenter.y + (tc.y - tcenter.y) * cos_table[ta.alpha / ROTATION_STEP] - (tc.z - tcenter.z) * sin_table[ta.alpha / ROTATION_STEP];
    tc.z = tcenter.z + (tc.y - tcenter.y) * sin_table[ta.alpha / ROTATION_STEP] + (tc.z - tcenter.z) * cos_table[ta.alpha / ROTATION_STEP];

    // Rotation about the Y-axis (beta)
    tc.x = tcenter.x + (tc.z - tcenter.z) * sin_table[ta.beta / ROTATION_STEP] + (tc.x - tcenter.x) * cos_table[ta.beta / ROTATION_STEP];
    tc.z = tcenter.z + (tc.z - tcenter.z) * cos_table[ta.beta / ROTATION_STEP] - (tc.x - tcenter.x) * sin_table[ta.beta / ROTATION_STEP];

    // Rotation about the Z-axis (gamma)
    tc.x = tcenter.x + (tc.x - tcenter.x) * cos_table[ta.gamma / ROTATION_STEP] - (tc.y - tcenter.y) * sin_table[ta.gamma / ROTATION_STEP];
    tc.y = tcenter.y + (tc.x - tcenter.x) * sin_table[ta.gamma / ROTATION_STEP] + (tc.y - tcenter.y) * cos_table[ta.gamma / ROTATION_STEP];
    
    return tc;
}

AngleOfRotation create_angle_of_rotation(const int talpha, const int tbeta, const int tgamma) {
    AngleOfRotation angle; 
    angle.alpha = talpha; angle.beta = tbeta; angle.gamma = tgamma;
    return angle;
}

Atom create_atom(const int tx, const int ty, const int tz, const AtomType ttype){
    Atom atom; 
    atom.coordinate = create_coordinate(tx, ty, tz);
    atom.type = ttype;
    return atom;
}

int read_pdb_to_biomolecule(const char* filename, Biomolecule* m){

    // Local variable declarations
    FILE* fptr;
    int flag = 1;
    float x, y, z;
    int t;
    int i;

    // Allocate memory to the biomolecule m and handling memory allocation error
    if((m = (Biomolecule*)malloc(sizeof(Biomolecule))) == NULL) {
        flag = 0;
        MEMORY_ERROR;
        printf("Function: read_pdb_to_biomolecule | preprocessing.c | 59\n");
        return 0;
    }
    
    // Open file <filename> and handle file error
    if((fptr = fopen(filename, "r")) == NULL) {
        flag = -1;
        FILE_ERROR;
        printf("Function: read_pdb_to_biomolecule | preprocessing.c | 67\n");
        return flag;
    }

    // Allocate memory to m->atoms and handle memory error
    fscanf(fptr, "%d", &(m->number_of_atoms));
    if((m->atoms = (Atom*)malloc(sizeof(Atom) * (m->number_of_atoms))) == NULL) {
        flag = 0;
        MEMORY_ERROR;
        printf("Function: read_pdb_to_biomolecule | preprocessing.c | 76\n");
        return flag;
    }

    // Finally if allocation is successfull, reading data into m
    for(i= 0; i < m->number_of_atoms; i++) {
        fscanf(fptr, "%f %f %f %d", &x, &y, &z, &t);
        m->atoms[i] = create_atom(x, y, z, (AtomType)t);
    }
    m->center = create_coordinate(0,0,0);
    m->max = create_coordinate(0,0,0);
    m->min = create_coordinate(0,0,0);
    fclose(fptr) ;
    return flag;
}

void filter_coordinates_of_biomolecule(Biomolecule* m) {

    // Local variable declarations
    int i;
    float xavg, yavg, zavg;
    xavg = yavg = zavg =  0;
    float xmin, ymin, zmin;
    xmin =ymin = zmin =0;
    
    // Finding the coordinates of center (geometric) and the most negative coordinates (min)
    for(i = 0; i < m->number_of_atoms; i++) {
        xavg += m->atoms[i].coordinate.x;
        yavg += m->atoms[i].coordinate.y;
        zavg += m->atoms[i].coordinate.z;
        xmin = min(xmin, m->atoms[i].coordinate.x);
        ymin = min(ymin, m->atoms[i].coordinate.y);
        zmin = min(zmin, m->atoms[i].coordinate.z);
    }
    m->center.x = xavg / m->number_of_atoms;
    m->center.y = yavg / m->number_of_atoms;
    m->center.z = zavg / m->number_of_atoms;

    // Shifting all coordinates by the necessary value to obtain all positive coordinates
    xmin = -1.0 * xmin; ymin = -1.0 * ymin; zmin = -1.0 * zmin;

    for(i = 0; i < m->number_of_atoms; i++) {
        m->atoms[i].coordinate.x += xmin;
        m->atoms[i].coordinate.y += ymin;
        m->atoms[i].coordinate.z += zmin;
    }
}


void center_coordinates_of_biomolecule(Biomolecule* m) {
    // Local varaible declarations (calculate the shift in coordinates needed)
    Coordinate grid_center = create_coordinate(GRID_SIZE * RESOLUTION / 2.0 , GRID_SIZE * RESOLUTION /2.0, GRID_SIZE *RESOLUTION / 2.0);
    Coordinate shift = create_coordinate(grid_center.x - m->center.x, grid_center.y - m->center.y, grid_center.z - m->center.z);
    int i;

    // Shift all atomic coordinates by the obtained shift-coordinates
    for(i = 0; i < m->number_of_atoms; i++) {
        m->atoms[i].coordinate.x += shift.x;
        m->atoms[i].coordinate.y += shift.y;
        m->atoms[i].coordinate.z += shift.z;
    }
}

