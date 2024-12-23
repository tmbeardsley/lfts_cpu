// #############################################################################
// A class for dealing with the L-FTS simulation input parameters.
// Reads the input parameters from file. Not fully encapsulated as *m and 
// *L are exposed for convenience. Will act as a helper class that  
// automatically updates the values of derived parameters in future code 
// iterations (e.g., box-altering move and adaptive time step).
// #############################################################################

#include "lfts_params.h"


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


lfts_params::lfts_params(std::string inputFile) {
    read_Input_Params(inputFile);
}


lfts_params::~lfts_params() {
}


// Print parameters that were read from the input file to standard output
void lfts_params::outputParameters() {
    std::cout << "N = "             << N_           << std::endl;
    std::cout << "NA = "            << NA_          << std::endl;
    std::cout << "NB = "            << NB_          << std::endl;
    std::cout << "chi_b = "         << chi_b_       << std::endl;
    std::cout << "C = "             << C_           << std::endl;
    std::cout << "dt = "            << dt_          << std::endl;
    std::cout << "XeN = "           << XeN_         << std::endl;
    std::cout << "isXeN = "         << isXeN_       << std::endl;
    std::cout << "m[0] = "          << m_[0]        << std::endl;
    std::cout << "m[1] = "          << m_[1]        << std::endl;
    std::cout << "m[2] = "          << m_[2]        << std::endl;
    std::cout << "L[0] = "          << L_[0]        << std::endl;
    std::cout << "L[1] = "          << L_[1]        << std::endl;
    std::cout << "L[2] = "          << L_[2]        << std::endl;
    std::cout << "equil_its = "     << equil_its_   << std::endl;
    std::cout << "sim_its = "       << sim_its_     << std::endl;
    std::cout << "sample_freq = "   << sample_freq_ << std::endl;
    std::cout << "save_freq_ = "    << save_freq_   << std::endl;
    std::cout << "loadType_ = "     << loadType_    << std::endl;
    std::cout << "M_ = "            << M_           << std::endl;
    std::cout << "Mk_ = "           << Mk_          << std::endl;
    std::cout << "V_ = "            << V_           << std::endl;
    std::cout << "sigma_ = "        << sigma_       << std::endl;
}


void lfts_params::saveOutputParams(std::string fileName, bool append) {
    double XN_out = chi_b_;
    std::ofstream outstream;
    if (append) outstream.open(fileName,std::ios_base::app);
    else outstream.open(fileName);
    if (isXeN_ == 1) XN_out = XeN_;
    outstream << N_ << " " << NA_ << " " << XN_out << " " << C_ << " " << dt_ << " " << isXeN_ << std::endl;
    outstream << m_[0] << " " << m_[1] << " " << m_[2] << " " << L_[0] << " " << L_[1] << " " << L_[2] << std::endl;
    outstream << equil_its_ << " " << sim_its_ << " " << sample_freq_ << " " << save_freq_ << " " << 1 << std::endl;
    outstream.close();
}


// Read the simulation input parameters from file (first line)
void lfts_params::read_Input_Params(std::string fileName) {
    std::ifstream instream;
    instream.open(fileName);
    
    // Read the simulation parameters
    instream >> N_ >> NA_ >> chi_b_ >> C_ >> dt_ >> isXeN_;
    instream >> m_[0] >> m_[1] >> m_[2] >> L_[0] >> L_[1] >> L_[2];
    instream >> equil_its_ >> sim_its_ >> sample_freq_ >> save_freq_ >> loadType_;

    // Redefine variables to contain and integer number of sample_freq_ periods
    equil_its_ = (equil_its_/sample_freq_)*sample_freq_;
    sim_its_ = (sim_its_/sample_freq_)*sample_freq_;  

    // Transform from Xe to Xb if necessary
    double z_inf = z_inf_discrete(L_, m_, N_, C_*C_);
    if (isXeN_ == 1) {
        XeN_ = chi_b_;
        chi_b_ = chi_b_/z_inf;
    } else {
        XeN_ = chi_b_*z_inf;
    }

    calculate_Derived_Params();
}

// Calculate derived parameters (should only be called after read_Input_Params())
void lfts_params::calculate_Derived_Params() {
    M_ = m_[0]*m_[1]*m_[2];
    Mk_ = m_[0]*m_[1]*(m_[2]/2+1);
    V_ = L_[0]*L_[1]*L_[2];
    sigma_ = sqrt(2.0*M_*dt_/(C_*V_));
    NB_ = N_-NA_;
}

// Calculate z_infinity (discrete chain)
double lfts_params::z_inf_discrete(double *L, int *m, int N, double Nbar, int tmax) {
    double sum, prod=1.0, X, ell, R0, R0dL[3];
    sum = 0.5;
    R0 = sqrt((double) N);
    for (int i=0; i<3; i++) { 
        R0dL[i] = R0*L[i]/m[i];
        prod *= R0dL[i];
    }
    ell = pow(prod, 1.0/3.0);

    for (int t = 1; t <= tmax; t++) {
        X = (M_PI/ell)*sqrt(t/6.0);
        prod = 1.0;
        for (int i=0; i<3; i++) prod *= erf(X*ell/R0dL[i]);
        sum += pow(sqrt(M_PI)/(2*X), 3.0) * prod;
    }
    sum *= 2*R0/(pow(ell, 3.0)*sqrt(Nbar));
    X = (M_PI/ell)*sqrt((0.5+tmax)/6.0);

    return 1.0 - sum - 3*R0/(ell*sqrt(M_PI*Nbar)*X);
}

