// #############################################################################
// Performs a langevin update of w-(r) and keeps track of symmetrised noise
// #############################################################################

#include "langevin.h"

langevin::langevin(std::mt19937_64 &RNG, double sigma, int M)
    : noise_cpu_(std::make_unique<double[]>(2*M))               // Allocate memory for Gaussian random noise
{
    M_ = M;

    // Returns Gaussian random noise given a random number generator
    normRand_ = std::normal_distribution<double>(0.0, sigma);

    // Generate initial "previous" Gaussian random noise
    for (int r=0; r<M_; r++) noise_cpu_[r] = normRand_(RNG);
    noise_cpu_prev_ = noise_cpu_.get();
    noise_cpu_new_ = noise_cpu_.get() + M_;
}


// Perform a Langevin update of the fields using symmetrised noise
void langevin::step_wm(double* w, std::mt19937_64 &RNG, double XbN, double sigma, double dt)
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
