// ######################################################
// Provides public method: calc_concs(double *w), 
// to calculate concentrations (used in Anderson mixing)
// ######################################################

#include "diblock.h"

// Constructor
diblock::diblock(int NA, int NB, int *m, double *L, int M, int Mk) :
    qr_(std::make_unique<double[]>(2*(NA+NB+1)*M)),                     // Allocate memory for qr_
    q1_(std::make_unique<double*[]>(NA+NB+1)),                          // Allocate arrays of pointers for q_{j=1...N}(r) and q^_{j=1...N}(r)
    q2_(std::make_unique<double*[]>(NA+NB+1)),
    h_(std::make_unique<double[]>(2*M)),                                // Allocate memory for h_
    Step_(std::make_unique<step>(NA, NB, m, L, M, Mk)),                 // New step object containing methods to get next monomer's propagators
    M_(M),
    NA_(NA),
    N_(NA+NB)
{
    // Assign pointers such that q_{1}(r) and q_{N}(r) are in contigious memory,
    // as are q_{2}(r) and q_{N-1}(r), q_{3}(r) and q_{N-2}(r)... etc. (required for fftw_plan_many_dft_r2c(...) in step() object)
    for (int i=1; i<=N_; i++) {
        q1_[i] = qr_.get() + 2*i*M_;
        q2_[N_+1-i] = qr_.get() + (2*i+1)*M_;
    }
}


// Calculates phi-(r) and phi+(r): w+2*M -> phi-(0), w+3*M -> phi+(0).
// Returns ln(Q)
double diblock::calc_concs(double *w) {
    int i, r;
    double Q;
    double *phiA=w+2*M_;
    double *phiB=w+3*M_;

    // Calculate hA[r] and hB[r]
    prepare_h(h_.get(), w, N_, M_);

    // Set initial conditions: q[1][r]=hA[r] and q^[N][r]=hB[r] for all r
    for (r=0; r<2*M_; r++) q1_[1][r] = h_[r];

    // Step the propagators q1 and q2 for each subsequent monomer (note q[i],q^[N+1-i]... contigious in memory)
    for (i=1; i<N_; i++) Step_->fwd(q1_[i], q1_[i+1], h_.get(), i);

    // Calculate single-chain partition function using a Thrust reduction sum
    Q = 0.0;
    for (r=0; r<M_; r++) Q += q1_[N_][r];
    Q /= M_;

    // Calculate concentrations
    for (r=0; r<2*M_; r++) phiA[r] = 0.0;
    for (i=1; i<=NA_; i++) sum_phi(phiA, q1_[i], q2_[i], M_);
    for (i=NA_+1; i<=N_; i++) sum_phi(phiB, q1_[i], q2_[i], M_);

    // Normalise concentrations.
    // NOTE: normalize_phi DESTROYS THE CONTENTS OF h_[] TO SAVE MEMORY
    normalize_phi(w, h_.get(), Q, N_, M_);

    return log(Q);
}


// Calculate hA[r] and hB[r]: h -> hA[0], h+M -> hB[0] 
void diblock::prepare_h(double *h, double *w, const int N, const int M) {
    double *wm=w, *wp=w+M;
    double *hA = h, *hB = h+M;

    for (int r=0; r<M; r++) {
        hA[r] = exp(-(wm[r]+wp[r])/N);
        hB[r] = exp(-(-wm[r]+wp[r])/N);
    }
}


// Multiply and sum propagators for calculating either phiA[r] or phiB[r]
void diblock::sum_phi(double *phi, double *q1, double *q2, const int M) {
    for (int r=0; r<M; r++) {
        phi[r] += q1[r]*q2[r];
    }
}


// Normalise concentrations phi-(r) and phi+(r): w_+2*M -> phim[0],  w_+3*M -> phip[0]
// Note: Destroys the contents of h[] to save memory.
void diblock::normalize_phi(double *w, double *h, const double Q, const int N, const int M) {
    double *phiA=h, *phiB=h+M, *hA=h, *hB=h+M;
    double *phim=w+2*M, *phip=w+3*M;

    for (int r=0; r<M; r++) {
        phiA[r] = phim[r]/(N*Q*hA[r]);
        phiB[r] = phip[r]/(N*Q*hB[r]);
        phim[r] = phiA[r] - phiB[r];
        phip[r] = phiA[r] + phiB[r];
    }
}
