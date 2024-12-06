// #######################################################
// Performs Anderson Mixing to solve for the W+(r) field
// #######################################################

#pragma once
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <memory>
#include "diblock.h"

class anderson {
    int nhMax_;                 // Maximum # histories (default: 10)

    // Mathematical arrays
    std::unique_ptr<std::unique_ptr<double[]>[]> DD_;               
    std::unique_ptr<double[]> U_;
    std::unique_ptr<double[]> V_;
    std::unique_ptr<double[]> C_;
    std::unique_ptr<std::unique_ptr<double[]>[]> Dh_;
    std::unique_ptr<std::unique_ptr<double[]>[]> wh_;
    std::unique_ptr<double[]> A_cpy_;
    std::unique_ptr<double[]> Y_cpy_;

    // Simulation constants derived from the input file (see lfts_params.h for details)
    int M_;

    public:
        // Constructor
        anderson(int M, int maxHist=10);

        // Perform Anderson Mixing
        int mix(diblock *dbc, int maxIter, double errTol, double *w);


    private:
        // Function to create a 2D array using unique_ptr
        std::unique_ptr<std::unique_ptr<double[]>[]> array2d(const int m, const int n);

        // Perform matrix LU decomposition
        void LUdecomp(double *A, double *Y, double *X, const int n);

};