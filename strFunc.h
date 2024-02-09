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
#include "sorts.h"

class strFunc {             
    double *S_;                     // Collects sum of |wk_[k]| resulting from calls to: sample(double *w)
    double *K_;                     // Modulus of wavevector k
    double dK_;                     // Maximum allowed difference used to define like wave vectors for spherical averaging
    double coeff_;                  // A constant used in saving the structure function
    int *wt_;                       // Weighting of contribution from wavevector k
    int *P_;                        // Map transforming K_[] into ascending order
    int nsamples_;                  // Number of structure function samples taken
    fftw_plan wr_to_wk_;            // fftw plan transforming w-(r) to w-(k)
    std::complex<double> *wk_;      // w-(k)

    // Simulation constants derived from the input file (see lfts_params.h for details)
    double chi_b_;
    int Mk_;

    public:
        // Constructor
        strFunc(double *w, int *m, double *L, int M, int Mk, double CV, double chi_b, double dK = 1E-5) {
            Mk_ = Mk;
            dK_ = dK;
            chi_b_ = chi_b;
            nsamples_ = 0;
            coeff_ = CV/(chi_b*chi_b*M*M);
            wk_ = new std::complex<double>[Mk];
            K_ = new double[Mk];
            wt_ = new int[Mk];
            P_ = new int[Mk];

            // Allocate memory for S(k) and zero the elements
            S_ = new double[Mk];
            for (int k=0; k<Mk; k++) S_[k] = 0.0;

            // Create an fftw plan for the Fourier transform: w-(r) -> w-(k)
            wr_to_wk_ = fftw_plan_dft_r2c(3, m, w, reinterpret_cast<fftw_complex*>(wk_), FFTW_PATIENT);

            // Populate the wavevector modulus array, K_
            calcK(K_, m, L);

            // Populate the map, P_, which puts the wavevector moduli, K_, into ascending order
            //calcSortedKMap(P_, K_);
            for (int k=0; k<Mk_; k++) P_[k] = k;
            sorts::quicksortMap(K_, P_, 0, Mk_-1);
        }

        // Sample norm(w-(k)) 
        void sample(double *w) {
            // Transform w-(r) to k-space to get w-(k)
            fftw_execute(wr_to_wk_);

            // Sample the norm of w-(k) for each wavevector and add to its sum
            for (int k=0; k<Mk_; k++) S_[k] += norm(wk_[k]);

            // Increment the number of samples
            nsamples_++;
        }

        // Output the spherically-averaged structure function to file
        void save(std::string fileName, int dp=8) {
            double S_sum = 0.0;
            int k, n_same = 0;
            std::ofstream out_stream;

            out_stream.open(fileName);
            out_stream.precision(dp);
            out_stream.setf(std::ios::fixed, std::ios::floatfield);

            // Spherical average of S(k)
            for (k=0; k<Mk_; k++) {
                // Take into account vector weighting from the FFT and sum S for repeated K-vectors
                S_sum += wt_[P_[k]] * ((coeff_/nsamples_)*S_[P_[k]] - 0.5/chi_b_);
                n_same += wt_[P_[k]];

                // Output value for current K-vector when difference in K exceeds tolerence dK_
                if ( (k==Mk_-1) || (fabs(K_[P_[k+1]]-K_[P_[k]]) > dK_) ) {
                    out_stream << K_[P_[k]] << "\t" << S_sum/n_same << std::endl;

                    // Reset summations for next K-vector
                    S_sum = 0.0;
                    n_same = 0;
                }
            } 
            out_stream.close();

        }

        // Destructor
        ~strFunc() {
            delete[] wk_;
            delete[] K_;
            delete[] wt_;
            delete[] P_;
            delete[] S_;
            fftw_destroy_plan(wr_to_wk_);
        }




    private:
        // Calculate the wavevector moduli and store in K[]
        void calcK(double *K, int *_m, double *_L) {

            int K0, K1, k;
            double kx_sq, ky_sq, kz_sq;

            for (k=0; k<Mk_; k++) wt_[k]=2;

            for (int k0=-(_m[0]-1)/2; k0<=_m[0]/2; k0++) {
                K0 = (k0<0)?(k0+_m[0]):k0;
                kx_sq = k0*k0/(_L[0]*_L[0]);

                for (int k1=-(_m[1]-1)/2; k1<=_m[1]/2; k1++) {
                    K1 = (k1<0)?(k1+_m[1]):k1;
                    ky_sq = k1*k1/(_L[1]*_L[1]);

                    for (int k2=0; k2<=_m[2]/2; k2++) {
                        kz_sq = k2*k2/(_L[2]*_L[2]);
                        k = k2 + (_m[2]/2+1)*(K1+_m[1]*K0);
                        K[k] = 2*M_PI*pow(kx_sq+ky_sq+kz_sq,0.5); 
                        if ((k2==0)||(k2==_m[2]/2)) wt_[k]=1;
                    }
                }
            }
        }

};
