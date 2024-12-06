// #######################################################################################
// Provides the public method: void fwd(...), which takes the propagators of the previous
// monomer as input and returns the propagators of the next monomer as output
// #######################################################################################

#pragma once
#include <math.h>
#include <complex>
#include <fftw3.h>                                  // Fast Fourier transforms
#include <memory>

class step {

    struct deleter {
        // this alias is important for unique_ptr
        using pointer = fftw_plan;
        void operator()(pointer plan) const { fftw_destroy_plan(plan); }
    };

    // Step-specific variables
    std::unique_ptr<double[]> g_;                   // Bond potential Boltzmann weight, Fourier transformed and /M_
    std::unique_ptr<std::complex<double>[]> qk_;    // Fourier transforms of q_{j=i}(k) and q^_{j=N+1-i}(k) in contigious memory (for fftw_plan_many())
    std::unique_ptr<double[]> qr_;                  // q_{j=i}(r) and q^_{j=N+1-i}(r) in contigious memory (for fftw_plan_many())
    std::unique_ptr<void, deleter> qr_to_qk_;       // cufft plan to transform q1[r] and q2[r] to k-space
    std::unique_ptr<void, deleter> qk_to_qr_;       // cufft plan to transform q1[k] and q2[k] to real-space

    // Simulation constants derived from the input file (see lfts_params.h for details)
    std::unique_ptr<int[]> m_;
    int NA_;
    int NB_;
    int M_;
    int Mk_;


    public:
        // Constructor
        step(int NA, int NB, int *m, double *L, int M, int Mk);

        // Calculate propagators for the next monomer given propagators of previous monomer
        // q_in  = q{i}(r), q^{N+1-i}(r)
        // q_out = q{i+1}(r), q^{N-i}(r)
        // h_ = hA...hB
        void fwd(double* q_in, double* q_out, double *h, int i);

    private:
        // Calculate the Boltzmann weight of the bond potential in k-space, _g[k]
        void update_g_lookup(double *L);

};




