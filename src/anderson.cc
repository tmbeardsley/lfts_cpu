// #######################################################
// Performs Anderson Mixing to solve for the W+(r) field
// #######################################################

#include "anderson.h"

// Constructor
anderson::anderson(int M, int maxHist) :
    nhMax_(maxHist),
    DD_(array2d(maxHist+1, maxHist+1)),
    U_(std::make_unique<double[]>((maxHist*maxHist))),
    V_(std::make_unique<double[]>((maxHist))),
    C_(std::make_unique<double[]>((maxHist))),
    Dh_(array2d(maxHist+1, M)),
    wh_(array2d(maxHist+1, M)),
    A_cpy_(std::make_unique<double[]>((maxHist*maxHist))),
    Y_cpy_(std::make_unique<double[]>((maxHist))),
    M_(M)
{
    const int DIM = nhMax_+1;
}


int anderson::mix(diblock *dbc, int maxIter, double errTol, double *w) {
    double lambda, err=1.0, S1;
    int    k, m, n, nh, r;
    double omp_sum = 0.0;

    for (k=1; k<maxIter && err>errTol; k++) {

        // Calculate concentrations. In: w-(r) and w+(r). Out: phi-(r) and phi+(r).
        dbc->calc_concs(w);

        // Copy arrays into working memory for current iteration
        for (r=0; r<M_; r++) { 
            Dh_[0][r] = w[r+3*M_] - 1.0;
            wh_[0][r] = w[r+M_];
        }

        // Sum of (phi-(r)-1.0)^2
        S1 = 0.0;
        for (r=0; r<M_; r++) S1 += pow(Dh_[0][r], 2.0);

        // Update mixing error
        err = pow(S1/M_,0.5);
        lambda = 1.0-pow(0.9,double(k));

        // Update the number of histories
        nh = (k<nhMax_+1)?k-1:nhMax_;

        // Perform summations for each history
        for (m=0; m<=nh; m++) {
            DD_[0][m] = 0.0;
            for (r=0; r<M_; r++) DD_[0][m] += Dh_[0][r]*Dh_[m][r];
        }

        if (k<2) {
            // Simple mixing
            for (r=0; r<M_; r++) w[r+M_] += lambda*Dh_[0][r];
        } else {   
            // Anderson mixing
            for (m=1; m<=nh; m++) {
                V_[m-1] = DD_[0][0]-DD_[0][m];
                for (n=1; n<=m; n++) {
                    U_[(m-1)*nh+n-1] = U_[(n-1)*nh+m-1] = DD_[0][0]-DD_[0][m]-DD_[0][n]+DD_[n][m];
                }
            }

            // Solve for small matrix C_[] on the host using LU decomposition (U_[] and V_[] unchanged)
            LUdecomp(U_.get(),V_.get(),C_.get(),nh);

            // Initial simple mixing step: updates w+(r)
            for (r=0; r<M_; r++) w[r+M_] += lambda*Dh_[0][r];

            // Perform Anderson Mixing for each history: updates w+(r)
            for (n=1; n<=nh; n++) {
                for (r=0; r<M_; r++) w[r+M_] += C_[n-1]*( (wh_[n][r]+lambda*Dh_[n][r]) - (wh_[0][r]+lambda*Dh_[0][r]) );
            }
        }

        // Field and deviation of current step become history n
        n=1+(k-1)%(nhMax_);
        for (r=0; r<M_; r++) {
            Dh_[n][r] = Dh_[0][r];
            wh_[n][r] = wh_[0][r];
        }
        DD_[n][n] = DD_[0][0];
        for (m=1; m<n; m++) DD_[m][n] = DD_[0][m];
        for (m=n+1; m<=nh; m++) DD_[n][m] = DD_[0][m];
    }

    return k;
}


// Function to create a 2D array using unique_ptr
std::unique_ptr<std::unique_ptr<double[]>[]> anderson::array2d(const int m, const int n) {
    auto a = std::make_unique<std::unique_ptr<double[]>[]>(m);  // Array of unique_ptr
    for (int i = 0; i < m; ++i) {
        a[i] = std::make_unique<double[]>(n);  // Allocate array for each row
    }
    return a;
}


void anderson::LUdecomp(double *A, double *Y, double *X, const int n) {
    int s;

    // Make copies of A and Y since LU will changes matrix/vector contents.
    for (int i=0; i<n; i++) {
        Y_cpy_[i] = Y[i];
        for (int j=0; j<n; j++) {
            A_cpy_[i*n+j] = U_[i*n+j];
        }
    }

    // Create matrix and vector views to use gsl library functions
    gsl_matrix_view a = gsl_matrix_view_array(A_cpy_.get(), n, n);
    gsl_vector_view y = gsl_vector_view_array(Y_cpy_.get(), n);
    gsl_vector_view x = gsl_vector_view_array(X, n);
    gsl_permutation *p = gsl_permutation_alloc(n);

    // Solve for x using LU decomposition via gsl library
    gsl_linalg_LU_decomp(&a.matrix, p, &s);
    gsl_linalg_LU_solve(&a.matrix, p, &y.vector, &x.vector);
    gsl_permutation_free(p);
}

