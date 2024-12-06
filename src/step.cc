// #######################################################################################
// Provides the public method: void fwd(...), which takes the propagators of the previous
// monomer as input and returns the propagators of the next monomer as output
// #######################################################################################

#include "step.h"


// Constructor
step::step(int NA, int NB, int *m, double *L, int M, int Mk) :
    g_(std::make_unique<double[]>(Mk)),                         // Allocate memory for g and calculate it
    qk_(std::make_unique<std::complex<double>[]>(2*Mk)),        // Allocate memory for q1[k] and q2[k], stored in contigious memory
    qr_(std::make_unique<double[]>(2*M)),                       // Allocate memory for q1[r] and q2[r], stored in contigious memory
                                                                // Configure cufft plans. cufftPlanMany used for batched processing
    qr_to_qk_(fftw_plan_many_dft_r2c(3, m, 2, qr_.get(), NULL, 1, M, reinterpret_cast<fftw_complex*>(qk_.get()), NULL, 1, Mk, FFTW_PATIENT), deleter()),
    qk_to_qr_(fftw_plan_many_dft_c2r(3, m, 2, reinterpret_cast<fftw_complex*>(qk_.get()), NULL, 1, Mk, qr_.get(), NULL, 1, M, FFTW_PATIENT), deleter()),
    m_(std::make_unique<int[]>(3)),                             // Lattice size
    NA_(NA),
    NB_(NB),
    M_(M),
    Mk_(Mk)
{

    for (int i=0; i<3; i++) m_[i] = m[i];

    update_g_lookup(L);

}


// Calculate propagators for the next monomer given propagators of previous monomer
// q_in  = q{i}(r), q^{N+1-i}(r)
// q_out = q{i+1}(r), q^{N-i}(r)
// h_ = hA...hB
void step::fwd(double* q_in, double* q_out, double *h, int i)
{
    // The qr_to_qk_ plan is linked to qr_[], so have to copy in data.
    for (int r=0; r<2*M_; r++) qr_[r] = q_in[r];

    // Fourier transform q{i}(r),q^{N+1-i}(r) to k-space: qk_(k)
    fftw_execute(qr_to_qk_.get());

    // Multiply qk_(k) by g_(k)
    for (int r=0; r<Mk_; r++) {
        qk_[r]      *= g_[r];
        qk_[r+Mk_]  *= g_[r];
    }

    // Fourier transform qk_(k) to real-space: qr_[r]
    fftw_execute(qk_to_qr_.get());

    // Multiply qr_[r] by hA[r] or hB[r] (depending on monomer index) to get q{i+1}(r)
    if (i < NA_) {
        for (int r=0; r<M_; r++) q_out[r] = qr_[r]*h[r];
    } else {
        for (int r=0; r<M_; r++) q_out[r] = qr_[r]*h[r+M_];
    }

    // Multiply qr_[r+M_] by hA[r] or hB[r] (depending on monomer index) to get q^{N-i}(r)
    if (i < NB_) {
        for (int r=0; r<M_; r++) q_out[r+M_] = qr_[r+M_]*h[r+M_];
    } else {
        for (int r=0; r<M_; r++) q_out[r+M_] = qr_[r+M_]*h[r];
    }
}


// Calculate the Boltzmann weight of the bond potential in k-space, _g[k]
void step::update_g_lookup(double *L) {
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





