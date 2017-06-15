#include "docker.h"


int main() {

    // Declaring local variables
    Biomolecule *protein, *dna;
    protein = dna = NULL;
    Result *result_stack = NULL;
    int dock_flag;

    // Reading biomolecule data into the pointers
    if(read_pdb_to_biomolecule(PROTEIN_TMP_FILEPATH, protein) == 0) {

        // Encountered error, print error message and clean up
        printf("\n****** STOPING DOCKER DUE TO ABOVE MENTIONED ERROR *****\n");
        printf("\n****** CLEANING UP ******\n");
        printf("\n****** NOW EXITING ******\n");
        return 0;
    }

    if(read_pdb_to_biomolecule(DNA_TMP_FILEPATH, dna) == 0) {
        // Encountered error, print error message and clean up
        printf("\n****** STOPING DOCKER DUE TO ABOVE MENTIONED ERROR *****\n");
        printf("\n****** CLEANING UP ******\n");
        free(protein);
        printf("\n****** NOW EXITING ******\n");
        return 0;
    }

    // All biomolecule data read, now generating trig tables
    if(generate_trig_tables() == 0) {
        // Encountered error, print error message and clean up
        printf("\n****** STOPING DOCKER DUE TO ABOVE MENTIONED ERROR *****\n");
        printf("\n****** CLEANING UP ******\n");
        free(protein);
        free(dna);
        printf("\n****** NOW EXITING ******\n");
        return 0;
    }

    // Trig tables generated, now filtering coordinates of biomolecules and centering the protein
    filter_coordinates_of_biomolecule(protein);
    filter_coordinates_of_biomolecule(dna);
    center_coordinates_of_biomolecule(protein);

    // All biomolecules have been processed, now checking constraints
    if(check_constraints(protein, dna) == 0) {
        // Encountered error, print error message and clean up
        printf("\n****** STOPING DOCKER DUE TO ABOVE MENTIONED ERROR *****\n");
        printf("\n****** CLEANING UP ******\n");
        free(protein);
        free(dna);
        free(sin_table);
        free(cos_table);
        printf("\n****** NOW EXITING ******\n");
        return 0;
    }

    // All costraints satisfied, now creating result stack
    result_stack = (Result*) malloc(sizeof(Result));

    // Check if result_stack allocation successfull
    if(result_stack == NULL) {
        // Encountered error, print error message and clean up
        printf("\n****** STOPING DOCKER DUE TO ABOVE MENTIONED ERROR *****\n");
        printf("\n****** CLEANING UP ******\n");
        free(protein);
        free(dna);
        free(sin_table);
        free(cos_table);
        printf("\n****** NOW EXITING ******\n");
        return 0;
    }

    // Initialize result stack
    result_stack->number_of_results = 0;
    result_stack->decoys = (Decoy*) malloc(sizeof(Decoy) * (360 / ROTATION_STEP) * (360 / ROTATION_STEP) * (360 / ROTATION_STEP));

    // Check if result_stack->decoys allocation successfull
    if(result_stack->decoys == NULL) {
        // Encountered error, print error message and clean up
        printf("\n****** STOPING DOCKER DUE TO ABOVE MENTIONED ERROR *****\n");
        printf("\n****** CLEANING UP ******\n");
        free(protein);
        free(dna);
        free(sin_table);
        free(cos_table);
        free(result_stack);
        printf("\n****** NOW EXITING ******\n");
        return 0;
    }

    // Result stack successfully created and initialized, now starting docking procedure
    dock_flag = dock_biomolecules(protein, dna, result_stack);

    // Check the status flag of dock_biomolecules and handle errors
    if(dock_flag == 0) {
        // Encountered error, print error message and clean up
        printf("\n****** STOPING DOCKER DUE TO ABOVE MENTIONED ERROR *****\n");
        printf("\n****** CLEANING UP ******\n");
        free(protein);
        free(dna);
        free(sin_table);
        free(cos_table);
        free(result_stack->decoys);
        free(result_stack);
        printf("\n****** NOW EXITING ******\n");
        return 0;
    }

    // docker procedure complete, now write output to file
    if(write_result_to_file(OUTPUT_FILEPATH, result_stack) == 0) {
        // Encountered error, print error message and clean up
        printf("\n****** STOPING DOCKER DUE TO ABOVE MENTIONED ERROR *****\n");
        printf("\n****** CLEANING UP ******\n");
        free(protein);
        free(dna);
        free(sin_table);
        free(cos_table);
        free(result_stack->decoys);
        free(result_stack);
        printf("\n****** NOW EXITING ******\n");
        return 0;
    }

    // Output has been written to file. All procedures finished. Cleaning up heap and exiting
    printf("\n******* DOCKING PROCEDURE COMPLETE ******\n");
    printf("\n******* OUTPUT WRITTEN TO CONFIGURED OUTPUT FILE *******\n");
    printf("\n******* CLEANING UP *******\n");
    free(protein);
    free(dna);
    free(sin_table);
    free(cos_table);
    free(result_stack->decoys);
    free(result_stack);
    printf("\n******* NOW EXITING *******\n");

    // Everything went fine, return standard status 0
    return 0;
}
