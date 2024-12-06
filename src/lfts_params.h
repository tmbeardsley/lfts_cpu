// #############################################################################
// A class for dealing with the L-FTS simulation input parameters.
// Reads the input parameters from file. Not fully encapsulated as *m and 
// *L are exposed for convenience. Will act as a helper class that  
// automatically updates the values of derived parameters in future code 
// iterations (e.g., box-altering move and adaptive time step).
// #############################################################################

#pragma once
#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>
#include <iostream>


class lfts_params {

    // Input file parameters
    int    N_;                  // Total polymer length
    int    NA_;                 // Length of polymer A-block
    int    NB_;
    double chi_b_;              // Bare chi of the simulation
    double C_;                  // Dimensionless concentration, Nbar^0.5
    double dt_;                 // Langevin time step multiplied by N
    int    isXeN_;              // Whether XN in the input file is bare or effective
    int    m_[3];               // Number of mesh points [mx,my,mz]
    double L_[3];               // Dimensions of simulation box [Lx,Ly,Lz] in units of aN^0.5
    int    equil_its_;          // Number of Langevin steps for equilibration
    int    sim_its_;            // Number of Langevin steps for statistics
    int    sample_freq_;        // Number of Langevin steps between samples
    int    save_freq_;          // Number of steps between saving statistics to file
    int    loadType_;           // Whether to load from file or create a new field

    // Derived parameters
    int    M_;                  // Total number of field mesh points
    int    Mk_;                 // Total number of field mesh points in k-space
    double V_;                  // Volume of the simulation box
    double sigma_;              // Standard deviation of random noise multiplied by N
    double XeN_;                // Effective chi multiplied by N

    public: 
        lfts_params(std::string inputFile);

        ~lfts_params();

        // Print parameters that were read from the input file to standard output
        void outputParameters();

        // Getters
        inline int N() { return N_; }
        inline int NA() { return NA_; }
        inline int NB() { return NB_; }
        inline double XbN() { return chi_b_; }
        inline double XeN() { return XeN_; }
        inline double C() { return C_; }
        inline double dt() { return dt_; }
        inline int isXeN() { return isXeN_; }
        inline int mx() { return m_[0]; }
        inline int my() { return m_[1]; }
        inline int mz() { return m_[2]; }
        inline int m(int dim) { return m_[dim]; }
        inline int* m() { return m_; }
        inline double Lx() { return L_[0]; }
        inline double Ly() { return L_[1]; }
        inline double Lz() { return L_[2]; }
        inline double L(int dim) { return L_[dim]; }
        inline double* L() { return L_; }
        inline int equil_its() { return equil_its_; }
        inline int sim_its() { return sim_its_; }
        inline int sample_freq() { return sample_freq_; }
        inline int save_freq() { return save_freq_; }
        inline int loadType() { return loadType_; }
        inline int M() { return M_; }
        inline int Mk() { return Mk_; }
        inline double V() { return V_; }
        inline double sigma() { return sigma_; }
        inline double n() { return C_*V_; }                // Total number of polymers in the system

        void saveOutputParams(std::string fileName, bool append=false);


    private:

        // Read the simulation input parameters from file (first line)
        void read_Input_Params(std::string fileName);

        // Calculate derived parameters (should only be called after read_Input_Params())
        void calculate_Derived_Params();

        // Calculate z_infinity (discrete chain)
        double z_inf_discrete(double *L, int *m, int N, double Nbar, int tmax=100);

};