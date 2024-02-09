// #############################################################################
// Performs a langevin update of w-(r) and keeps track of symmetrised noise
// #############################################################################
#pragma once
#include <random>

class langevin {
    std::normal_distribution<double> normRand_;     // Samples gaussian random noise
    double *noise_cpu_;                             // Array holding random noise for current step and previous step
    double *noise_cpu_new_;                         // Pointer to portion of memory for new noise in noise_cpu_[]
    double *noise_cpu_prev_;                        // Pointer to portion of memory for previous noise in noise_cpu_[]

    // Simulation constants derived from the input file (see lfts_params.h for details)
    int    M_;

    public:
        langevin(std::mt19937_64 &RNG, double sigma, int M) {
            M_ = M;

            // Allocate memory for Gaussian random noise
            noise_cpu_ = new double[2*M_];

            // Returns Gaussian random noise given a random number generator
            normRand_ = std::normal_distribution<double>(0.0, sigma);

            // Generate initial "previous" Gaussian random noise
            for (int r=0; r<M_; r++) noise_cpu_[r] = normRand_(RNG);
            noise_cpu_prev_ = noise_cpu_;
            noise_cpu_new_ = noise_cpu_ + M_;
        }

        ~langevin() {
            delete[] noise_cpu_;
        }

        // Perform a Langevin update of the fields using symmetrised noise
        void step_wm(double* w, std::mt19937_64 &RNG, double XbN, double sigma, double dt)
        {
            double *ptr_tmp;

            // Create new random noise
            for (int r=0; r<M_; r++) noise_cpu_new_[r] = normRand_(RNG);

            // Update the w-(r) field
            for (int r=0; r<M_; r++) w[r] += -(w[r+2*M_]+2*w[r]/XbN)*dt + 0.5*(noise_cpu_prev_[r]+noise_cpu_new_[r]);

            // Update the noise pointers to avoid copying data between steps
            ptr_tmp = noise_cpu_prev_;
            noise_cpu_prev_ = noise_cpu_new_;
            noise_cpu_new_ = ptr_tmp;
        }

};