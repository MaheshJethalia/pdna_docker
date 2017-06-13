#include "docker.h"

 
int dock_biomolecules(Biomolecule* P, Biomolecule* D, Result* result_stack){
    // Declaring local variables here
    Configuration* protein, *dna;           // for storing current rotational configuration of protein and dna biomols
    Matrix *protein_space_matrix, *dna_space_matrix;      // for storing the space discretized matrix of protein and dna biomols
    double* geometric_correlation_function;          // for storing the geometric correlation of 
    fftw_complex* protein_space_matrix_transform, *dna_space_matrix_transform;  // FFT transforms
    fftw_plan fft_plan_protein, fft_plan_dna, ifft_plan;    // FFT plans
    AngleOfRotation current_angle;

    // Create matrices protein_space_matrix and dna_space_matrix for configurations of p and d respectively
    if((protein_space_matrix = (Matrix*) malloc (sizeof(Matrix))) == NULL) {
        MEMORY_ERROR;
        printf("Function: get_geometric_complementarity_score | scoring.h | 14");
        return 0;
    }

    if((dna_space_matrix = (Matrix*) malloc(sizeof(Matrix))) == NULL) {
        MEMORY_ERROR;
        printf("Function: get_geometric_complementarity_score | scoring.h | 20");
        return 0;
    }



    // Allocate memory to arrays to be transformed and handle errors
    if((protein_space_matrix->value = (double*) fftw_malloc (sizeof(double) * GRID_SIZE *GRID_SIZE *GRID_SIZE)) == NULL) {
        MEMORY_ERROR;
        printf("Function: get_geometric_complementarity_score | scoring.c | 26");
        return 0;
    }

    if((dna_space_matrix->value = (double*) fftw_malloc(sizeof(double) * GRID_SIZE * GRID_SIZE * GRID_SIZE)) == NULL) {
        MEMORY_ERROR;
        printf("Function: get_geometric_complementarity_score | scoring.c | 32");
        return 0;
    }

    if((geometric_correlation_function = (double*) fftw_malloc(sizeof(double) * GRID_SIZE * GRID_SIZE * GRID_SIZE)) == NULL) {
        MEMORY_ERROR;
        printf("Function: get_geometric_complementarity_score | scoring.c | 38");
        return 0;
    }

    if((protein_space_matrix_transform = (fftw_complex*)fftw_alloc_complex(sizeof(fftw_complex) * GRID_SIZE * GRID_SIZE * (GRID_SIZE / 2 + 1))) == NULL) {
        MEMORY_ERROR;
        printf("Function: get_geometric_complementarity_score | scoring.c | 44");
        return 0;
    }

    if((dna_space_matrix_transform = (fftw_complex*)fftw_alloc_complex(sizeof(fftw_complex) * GRID_SIZE * GRID_SIZE * (GRID_SIZE / 2 + 1))) == NULL) {
        MEMORY_ERROR;
        printf("Function: get_geometric_complementarity_score | scoring.c | 50");
        return 0;
    }

    // All allocations OK! Create fftw and ifftw plans
    fft_plan_protein = fftw_plan_dft_r2c_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE, protein_space_matrix->value, protein_space_matrix_transform, FFTW_MEASURE);
    fft_plan_dna = fftw_plan_dft_r2c_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE, dna_space_matrix->value, dna_space_matrix_transform , FFTW_MEASURE);
    ifft_plan = fftw_plan_dft_c2r_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE/2 + 1, protein_space_matrix_transform, geometric_correlation_function, FFTW_MEASURE);


    /******* All plans created. Now start the docking procedure. *******/

    // Create non rotational protein configuration:
    protein = dna = NULL;
    if(create_configuration(protein, create_angle_of_rotation(0, 0, 0), P) == 0) {
        // create_configuration returned an error. Clean up and exit. 
    }

    // Creating the space discretized matrix
    create_geometric_core(protein, protein_space_matrix, RHO);
    create_geometric_surface(protein_space_matrix);
    // Finished creating static protein space discretized matrix and it's surface 
    
    // Now ready to transform protein space matrix
    fftw_execute(fft_plan_protein);
    // Protein space matrix transformed and kept ready
    
    current_angle = create_angle_of_rotation(0, 0, 0);

    // Generate all possible rotational angles
    while(current_angle.alpha < 360) {
        while(current_angle.beta < 360) {
            while(current_angle.gamma < 360) {

                // Create dna configuration wrt new rotational angle current_angle
                if(create_configuration(dna, current_angle, D) == 0) {
                    // create_configuration returned error. Clean up and exit.
                    return 0;
                }

                // Memory allocation OK! Creating space matrix for dna configuration
                create_geometric_core(dna, dna_space_matrix, 1.0);
                // Finished creating the space discretized matrix of rotated dna config
                
                // Now ready to transform dna space matrix
                fftw_execute(fft_plan_dna);

                // finding the DFT of geometric correlation function now
                for(int i = 0; i<= (GRID_SIZE * GRID_SIZE * (GRID_SIZE/2 + 1)); i++) {

                    double real, im;
                    real = protein_space_matrix_transform[i][REAL];
                    im = protein_space_matrix_transform[i][IM];

                    protein_space_matrix_transform[i][REAL] = real * dna_space_matrix_transform[i][REAL] + im * dna_space_matrix_transform[i][IM];
                    protein_space_matrix_transform[i][IM] = real * dna_space_matrix_transform[i][IM] - im * dna_space_matrix_transform[i][REAL];
                }
                // protein_space_matrix_transform now stores the value of the transform of geometric correlation function
                    
                // Now ready to perform inverse transform
                fftw_execute(ifft_plan);

                // Need to extract highest correlation translations and store in Result stack
                insert_favourable_decoys_in_result_stack(result_stack, geometric_correlation_function, current_angle, protein->center);
                
                
                current_angle.gamma += ROTATION_STEP;
            }
            current_angle.beta += ROTATION_STEP;
        }
        current_angle.alpha += ROTATION_STEP;
    }


    /********** All iterations completed. Result Stack stores best geometric complementarity scores achieved ********/

    // Cleaning up the heap 

    return 1;
}

    
