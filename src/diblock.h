// ######################################################
// Provides public method: calc_concs(double *w), 
// to calculate concentrations (used in Anderson mixing)
// ######################################################

#pragma once
#include "step.h"
#include <iostream>

class diblock {

    // Diblock-specific variables              
    std::unique_ptr<double[]> qr_;  // Pointer to memory for propagators: q_{i}(r) and q^_{N+1-i}(r) are contigious in memory
    std::unique_ptr<double*[]> q1_; // Array of pointers to q_{j=i}(r), where j is the monomer index and i is array index
    std::unique_ptr<double*[]> q2_; // Array of pointers to q^_{j=N+1-i}(r), where j is the monomer index and i is array index
    std::unique_ptr<double[]> h_;   // Pointer to memory for hA(r) and hB(r)
    std::unique_ptr<step> Step_;    // Step object to get propagators for the next monomer

    // Simulation constants derived from the input file (see lfts_params.h for details)
    int M_;
    int NA_;
    int N_;

    public:
        // Constructor
        diblock(int NA, int NB, int *m, double *L, int M, int Mk);

        // Calculates phi-(r) and phi+(r): w+2*M -> phi-(0), w+3*M -> phi+(0).
        // Returns ln(Q)
        double calc_concs(double *w);

    private:
        // Calculate hA[r] and hB[r]: h -> hA[0], h+M -> hB[0] 
        void prepare_h(double *h, double *w, const int N, const int M);

        // Multiply and sum propagators for calculating either phiA[r] or phiB[r]
        void sum_phi(double *phi, double *q1, double *q2, const int M);

        // Normalise concentrations phi-(r) and phi+(r): w_+2*M -> phim[0],  w_+3*M -> phip[0]
        // Note: Destroys the contents of h[] to save memory.
        void normalize_phi(double *w, double *h, const double Q, const int N, const int M);
};