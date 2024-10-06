// #######################################################
// Performs Anderson Mixing to solve for the W+(r) field
// #######################################################

#pragma once
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

class anderson {
    int nhMax_;                 // Maximum # histories (default: 10)

    // Mathematical arrays
    double **DD_;               
    double *U_;
    double *V_;
    double *C_;
    double *Dh_mem_;
    double **Dh_;
    double *wh_mem_;
    double **wh_;

    // Simulation constants derived from the input file (see lfts_params.h for details)
    int M_;

    public:
        // Constructor
        anderson(int M, int maxHist=10) {
            nhMax_ = maxHist;
            const int DIM = nhMax_+1;
            M_ = M;

            // Mathematical array memory
            DD_ = array2d(DIM, DIM);
            U_ = new double[nhMax_*nhMax_];
            V_ = new double[nhMax_];
            C_ = new double[nhMax_];
            Dh_mem_ = new double[DIM*M_];
            Dh_ = new double*[DIM];
            for (int i=0; i<DIM; i++) Dh_[i] = Dh_mem_ + i*M_;
            wh_mem_ = new double[DIM*M_];
            wh_ = new double*[DIM];
            for (int i=0; i<DIM; i++) wh_[i] = wh_mem_ + i*M_;
        }
        

        int mix(diblock *dbc, int maxIter, double errTol, double *w) {
            double lambda, err=1.0, S1;
            int    k, m, n, nh, r;

            for (k=1; k<maxIter && err>errTol; k++) {

                // Calculate concentrations. In: w-(r) and w+(r). Out: phi-(r) and phi+(r).
                dbc->calc_concs(w);

                // Copy arrays into working memory for current iteration
                for (r=0; r<M_; r++) Dh_[0][r] = w[r+3*M_] - 1.0;
                for (r=0; r<M_; r++) wh_[0][r] = w[r+M_];

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
                    LUdecomp(U_,V_,C_,nh);

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

        // Destructor
        ~anderson() {
            deleteArray2d(DD_, nhMax_+1);
            delete[] U_;
            delete[] V_;
            delete[] C_;
            delete[] Dh_mem_;
            delete[] wh_mem_;
            delete[] Dh_;
            delete[] wh_;
        }

        private:
            // Return a 2d array of dimensions (m,n)
            double** array2d(const int m, const int n) {
	            double** a = new double*[m];
	            for (int i=0; i<m; i++) a[i] = new double[n];
                return a;
            }

            // Deallocate memory of a 2d array
            void deleteArray2d(double **a, const int m) {
                for (int i=0; i<m; i++) delete[] a[i];
                delete [] a;
            }

            void LUdecomp(double *A, double *Y, double *X, const int n) {
                double *A_cpy;
                double *Y_cpy;
                int s;

                // Make copies of A and Y since LU will changes matrix/vector contents.
                A_cpy = new double[n*n];
                Y_cpy = new double[n];
                for (int i=0; i<n; i++) {
                    Y_cpy[i] = Y[i];
                    for (int j=0; j<n; j++) {
                        A_cpy[i*n+j] = U_[i*n+j];
                    }
                }

                // Create matrix and vector views to use gsl library functions
                gsl_matrix_view a = gsl_matrix_view_array(A_cpy, n, n);
                gsl_vector_view y = gsl_vector_view_array(Y_cpy, n);
                gsl_vector_view x = gsl_vector_view_array(X, n);
                gsl_permutation *p = gsl_permutation_alloc(n);

                // Solve for x using LU decomposition via gsl library
                gsl_linalg_LU_decomp(&a.matrix, p, &s);
                gsl_linalg_LU_solve(&a.matrix, p, &y.vector, &x.vector);
                gsl_permutation_free(p);
                delete[] A_cpy, Y_cpy;
            }

};