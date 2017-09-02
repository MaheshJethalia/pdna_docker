#ifndef DOCKER_H
#define DOCKER_H

#include <stdio.h>
#include <fftw3.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "config.h"

// Stores the table of precalculated sine and cosine values to boost computational speed
extern float* sin_table ;
extern float* cos_table ;

// Enum to store atom types
typedef enum AtomType { NDEF, C, N} AtomType;

typedef struct Coordinate {
    float x, y, z;
} Coordinate;

typedef struct AngleOfRotation {
    int alpha, beta, gamma;
} AngleOfRotation;

typedef struct Atom {
    Coordinate coordinate;    
    AtomType type;
} Atom;

typedef struct Biomolecule {
    Atom* atoms;
    Coordinate center;
    int number_of_atoms;
} Biomolecule;

typedef struct Cofiguration {
    Biomolecule* parent;
    AngleOfRotation angle;
    int number_of_atoms;
    Atom* atoms;
    Coordinate center;
} Configuration;

typedef struct Matrix {
    Configuration* parent_config;
    double* value;
    int xmax, ymax, zmax, xmin, ymin, zmin;

} Matrix;

typedef struct Decoy {
    Coordinate dna_translation;
    AngleOfRotation dna_rotation;
    float geometric_score;
    float refinement_score;
} Decoy;

typedef struct Result {
    Decoy* decoys;
    int number_of_results;
} Result;

extern int generate_trig_tables();

extern Coordinate create_coordinate(const float tx, const float ty, const float tz);

extern Coordinate rotate_coordinate(Coordinate tc, const AngleOfRotation ta, const Coordinate center);

extern AngleOfRotation create_angle_of_rotation(const int talpha, const int tbeta, const int tgamma);

extern Atom create_atom(const int tx, const int ty, const int tz, const AtomType ttype);

extern int read_pdb_to_biomolecule(const char* filename, Biomolecule* m);

extern void filter_coordinates_of_biomolecule(Biomolecule* m);

extern void center_coordinates_of_biomolecule(Biomolecule* m);

extern float get_biomolecule_diameter(Biomolecule* m);

extern int check_constraints(Biomolecule* A, Biomolecule* B);

extern int write_configuration_to_pdb(const char* filename, Configuration* m);

extern int create_configuration(Configuration* config, const AngleOfRotation tangle, Biomolecule* tparent);

extern void create_geometric_core(Configuration* config, Matrix* tmp, int init_value);

extern void create_geometric_surface(Matrix* matrix);

extern int dock_biomolecules(Biomolecule* P, Biomolecule* D, Result* result_stack);

extern Decoy insert_favourable_decoy_in_result_stack(double* geometric_correlation_function, AngleOfRotation current_angle);

extern int write_result_to_file(const char* filename, Result* result_stack);

extern void set_refinement_score(Configuration* protein, Configuration* dna, Decoy d);

#endif
