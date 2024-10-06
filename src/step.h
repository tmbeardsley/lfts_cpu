// #######################################################################################
// Provides the public method: void fwd(...), which takes the propagators of the previous
// monomer as input and returns the propagators of the next monomer as output
// #######################################################################################

#pragma once
#include <math.h>
#include <complex>
#include <fftw3.h>                      // Fast Fourier transforms


class step {
    // Step-specific variables
    double *g_;                         // Bond potential Boltzmann weight, Fourier transformed and /M_
    std::complex<double> *qk_;          // Fourier transforms of q_{j=i}(k) and q^_{j=N+1-i}(k) in contigious memory (for fftw_plan_many())
    double *qr_;                        // q_{j=i}(r) and q^_{j=N+1-i}(r) in contigious memory (for fftw_plan_many())
    fftw_plan qr_to_qk_;                // cufft plan to transform q1[r] and q2[r] to k-space
    fftw_plan qk_to_qr_;                // cufft plan to transform q1[k] and q2[k] to real-space

    // Simulation constants derived from the input file (see lfts_params.h for details)
    int NA_;
    int NB_;
    int *m_;
    int M_;
    int Mk_;

    public:
        // Constructor
        step(int NA, int NB, int *m, double *L, int M, int Mk) {
            NA_ = NA;
            NB_ = NB;

            m_ = new int[3];
            for (int i=0; i<3; i++) m_[i] = m[i];

            M_ = M;
            Mk_ = Mk;

            // Allocate memory for g and calculate it
            g_ = new double[Mk_];
            update_g_lookup(L);

            // Allocate memory for q1[r] and q2[r], stored in contigious memory
            qr_ = new double[2*M_];

            // Allocate memory for q1[k] and q2[k], stored in contigious memory
            qk_ = new std::complex<double>[2*Mk_];

            // Configure cufft plans. cufftPlanMany used for batched processing
            qr_to_qk_ = fftw_plan_many_dft_r2c( 3,                                      //rank,
                                                m_,                                     //*n, 
                                                2,                                      //howmany,
                                                qr_,                                    //*in, 
                                                NULL,                                   //*inembed,
                                                1,                                      //istride, 
                                                M_,                                     //idist,
                                                reinterpret_cast<fftw_complex*>(qk_),   //*out, 
                                                NULL,                                   //*onembed,
                                                1,                                      //ostride, 
                                                Mk_,                                    //odist,
                                                FFTW_PATIENT                            //flags
                                            );

            qk_to_qr_ = fftw_plan_many_dft_c2r( 3,                                      //rank,
                                                m_,                                     //*n, 
                                                2,                                      //howmany,
                                                reinterpret_cast<fftw_complex*>(qk_),   //*in, 
                                                NULL,                                   //*inembed,
                                                1,                                      //istride, 
                                                Mk_,                                    //idist,
                                                qr_,                                    //*out, 
                                                NULL,                                   //*onembed,
                                                1,                                      //ostride, 
                                                M_,                                     //odist,
                                                FFTW_PATIENT                            //flags
                                            );
        }

        // Calculate propagators for the next monomer given propagators of previous monomer
        // q_in  = q{i}(r), q^{N+1-i}(r)
        // q_out = q{i+1}(r), q^{N-i}(r)
        // h_ = hA...hB
        void fwd(double* q_in, double* q_out, double *h, int i)
        {
            int r;

            // The qr_to_qk_ plan is linked to qr_[], so have to copy in data.
            for (r=0; r<2*M_; r++) qr_[r] = q_in[r];

            // Fourier transform q{i}(r),q^{N+1-i}(r) to k-space: qk_(k)
            fftw_execute(qr_to_qk_);

            // Multiply qk_(k) by g_(k)
            for (r=0; r<Mk_; r++) {
                qk_[r]      *= g_[r];
                qk_[r+Mk_]  *= g_[r];
            }

            // Fourier transform qk_(k) to real-space: qr_[r]
            fftw_execute(qk_to_qr_);

            // Multiply qr_[r] by hA[r] or hB[r] (depending on monomer index) to get q{i+1}(r)
            if (i < NA_) {
                for (r=0; r<M_; r++) q_out[r] = qr_[r]*h[r];
            } else {
                for (r=0; r<M_; r++) q_out[r] = qr_[r]*h[r+M_];
            }

            // Multiply qr_[r+M_] by hA[r] or hB[r] (depending on monomer index) to get q^{N-i}(r)
            if (i < NB_) {
                for (r=0; r<M_; r++) q_out[r+M_] = qr_[r+M_]*h[r+M_];
            } else {
                for (r=0; r<M_; r++) q_out[r+M_] = qr_[r+M_]*h[r];
            }
        }

        // Destructor
        ~step() {
            fftw_destroy_plan(qr_to_qk_);
            fftw_destroy_plan(qk_to_qr_);
            delete[] g_;
            delete[] qk_;
            delete[] m_;
            delete[] qr_;
        }


    private:
        // Calculate the Boltzmann weight of the bond potential in k-space, _g[k]
        void update_g_lookup(double *L) {
            int K0, K1, k, N=NA_+NB_;
            double K, kx_sq, ky_sq, kz_sq;

            for (int k0=-(m_[0]-1)/2; k0<=m_[0]/2; k0++) {
                K0 = (k0<0)?(k0+m_[0]):k0;
                kx_sq = k0*k0/(L[0]*L[0]);

                for (int k1=-(m_[1]-1)/2; k1<=m_[1]/2; k1++) {
                    K1 = (k1<0)?(k1+m_[1]):k1;
                    ky_sq = k1*k1/(L[1]*L[1]);

                    for (int k2=0; k2<=m_[2]/2; k2++) {
                        kz_sq = k2*k2/(L[2]*L[2]);
                        k = k2 + (m_[2]/2+1)*(K1+m_[1]*K0);
                        K = 2*M_PI*pow(kx_sq+ky_sq+kz_sq,0.5); 
                        g_[k] = exp(-K*K/(6.0*N))/M_; 
                    }
                }
            }
        }

};




