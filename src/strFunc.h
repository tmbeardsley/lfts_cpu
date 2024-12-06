// #####################################################################################
// Provides the public methods: void sample(...) and void save(...), 
// which take samples of the structure funtion, S(k), and save the spherically-averaged 
// S(k) to file.
// S(k) should only be calculated in simulations keeping L[] and XbN constant.
// #####################################################################################

#pragma once
#include <math.h>
#include <complex>
#include <fstream>
#include <fftw3.h>
#include <numeric>
#include <algorithm>
#include <memory>

class strFunc {

    struct deleter {
        // this alias is important for unique_ptr
        using pointer = fftw_plan;
        void operator()(pointer plan) const { fftw_destroy_plan(plan); }
    };
               
    std::unique_ptr<double[]> S_;                   // Collects sum of |wk_[k]| resulting from calls to: sample(double *w)                
    std::unique_ptr<double[]> K_;                   // Modulus of wavevector k
    double dK_;                                     // Maximum allowed difference used to define like wave vectors for spherical averaging
    double coeff_;                                  // A constant used in saving the structure function
    std::unique_ptr<int[]> wt_;                     // Weighting of contribution from wavevector k
    std::unique_ptr<int[]> P_;                      // Map transforming K_[] into ascending order
    int nsamples_;                                  // Number of structure function samples taken
    std::unique_ptr<std::complex<double>[]> wk_;    // w-(k)
    std::unique_ptr<void, deleter> wr_to_wk_;       // fftw plan transforming w-(r) to w-(k)

    // Simulation constants derived from the input file (see lfts_params.h for details)
    double chi_b_;
    int Mk_;

    public:
        // Constructor
        strFunc(double *w, int *m, double *L, int M, int Mk, double CV, double chi_b, double dK = 1E-5);

        // Sample norm(w-(k)) 
        void sample(double *w);

        // Output the spherically-averaged structure function to file
        void save(std::string fileName, int dp=8);

    private:
        // Calculate the wavevector moduli and store in K[]
        void calcK(double *K, int *_m, double *_L);

};
