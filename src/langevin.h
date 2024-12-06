// #############################################################################
// Performs a langevin update of w-(r) and keeps track of symmetrised noise
// #############################################################################
#pragma once
#include <random>
#include <memory>

class langevin {

    std::normal_distribution<double> normRand_;     // Samples gaussian random noise
    std::unique_ptr<double[]> noise_cpu_;           // Array holding random noise for current step and previous step
    double *noise_cpu_new_;                         // Pointer to portion of memory for new noise in noise_cpu_[]
    double *noise_cpu_prev_;                        // Pointer to portion of memory for previous noise in noise_cpu_[]

    // Simulation constants derived from the input file (see lfts_params.h for details)
    int    M_;

    public:
        langevin(std::mt19937_64 &RNG, double sigma, int M);

        // Perform a Langevin update of the fields using symmetrised noise
        void step_wm(double* w, std::mt19937_64 &RNG, double XbN, double sigma, double dt);
};