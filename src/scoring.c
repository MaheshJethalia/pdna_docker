#include "docker.h"

/****** NEED TO DEALLOCATE STUFF HERE ****************/
double* get_correlation_function(double* a, double* b) {
    // Declaring local variables:
    int i;    //counters
    double* c;
    double ar, ai;
    fftw_complex *out_a, *out_b;    // to store DFT of a and b;
    fftw_plan fft_plan, ifft_plan;            // plans for fourier transform and inverse fourier transform
    // NOTE: the DFT of c is stored by overwriting values in a to save memory

    /********  STEPS FOR PERFORMING THE "FORWARD" FOURIER TRANSFORMS ********/
    out_a = (fftw_complex*) malloc(sizeof(fftw_complex) * GRID_SIZE * GRID_SIZE * (GRID_SIZE/2 + 1)) ;
    if(out_a == NULL) {
        MEMORY_ERROR;
        printf("Function: get_correlation_function | scoring.c | 12");
        return NULL;
    }

    out_b = (fftw_complex*) malloc(sizeof(fftw_complex) * GRID_SIZE *GRID_SIZE * (GRID_SIZE/2 + 1));
    if(out_b == NULL) {
        MEMORY_ERROR;
        printf("Function: get_correlation_function | scoring.c | 19");
        return NULL;
    }

    // All okay. Proceeding with transform.
    fft_plan = fftw_plan_dft_r2c_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE, a, out_a, FFTW_MEASURE);
    fftw_execute(fft_plan);
    fft_plan = fftw_plan_dft_r2c_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE, b, out_b, FFTW_MEASURE);
    fftw_execute(fft_plan);

    /********* STEPS FOR CALCULATING TRANSFORMED CORRELATION *********/
    for(i = 0; i < GRID_SIZE * GRID_SIZE *(GRID_SIZE / 2 + 1); i++) {
        ar = out_a[i][REAL]; ai = out_a[i][IM];
        out_a[i][REAL] = ar * out_b[i][REAL] + ai * out_b[i][IM];
        out_a[i][IM] = -1.0 * ai * out_b[i][REAL] + ar * out_b[i][IM]; 
    }

    /********* STEPS FOR PERFORMING INVERSE TRANSFORM **********/
    c = (double*) malloc(sizeof(double) * GRID_SIZE * GRID_SIZE * GRID_SIZE);
    ifft_plan = fftw_plan_dft_c2r_3d(GRID_SIZE, GRID_SIZE, GRID_SIZE/2 + 1, out_a, c, FFTW_MEASURE);
    fftw_execute(ifft_plan);
    return c;
}


