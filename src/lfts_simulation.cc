// ######################################################################################
// Exposes public methods to perform an L-FTS simulation: equilibrate() and statistics()
// ######################################################################################

#include "lfts_simulation.h"


lfts_simulation::lfts_simulation(std::string inputFile) {

    // Check that input file exists before proceeding
    if (!file_IO::isValidFile(inputFile)) {
        std::cout << "ERROR => Cannot open one of the input files." << std::endl;
        exit(1);
    }

    // Read simulation parameters from the input file and allocate temporary host memory for fields
    std::cout << "Creating lfts_params object..." << std::endl;
    P_ = std::make_unique<lfts_params>(inputFile);
    P_->outputParameters();
    M_=P_->M();

    // Set up random number generator
    RNG_.seed(123456789);

    // Allocate memory for field array
    w_ = std::make_unique<double[]>(4*M_);

    // Create a new diblock object
    std::cout << "Creating diblock object..." << std::endl;
    dbc_ = std::make_unique<diblock>(P_->NA(), P_->NB(), P_->m(), P_->L(), M_, P_->Mk());

    // Create a new anderson mixing object
    std::cout << "Creating anderson object..." << std::endl;
    AM_ = std::make_unique<anderson>(M_, 10);

    // Set up a langevin object to upstate w-(r) at each step
    std::cout << "Creating langevin object..." << std::endl;
    Langevin_ = std::make_unique<langevin>(RNG_, P_->sigma(), M_);

    // Create new structure function object.
    // Note: Sk_ must be constructed before read_Input_Fields(...) as w-[r] contents get destroyed during fftw plan setup
    Sk_ = std::make_unique<strFunc>(w_.get(), P_->m(), P_->L(), M_, P_->Mk(), P_->n(), P_->XbN());

    // Read w-[r] and w+[r] from the input file
    if (P_->loadType() == 1) { 
        std::cout << "Loading input field..." << std::endl;
        file_IO::readArray(w_.get(), inputFile, 2*M_, 3);
    }
    else generate_field(w_.get(), P_->loadType());

    // Perform an initial mix to get phi-(r) and phi+(r) from the input fields
    AM_->mix(dbc_.get(),200,1e-4,w_.get());

    // Output initial fields
    saveStdOutputFile("w_0");
    file_IO::saveArray(w_.get()+2*M_, "phi_0", 2*M_);
}


// Destructor
lfts_simulation::~lfts_simulation() {

}


// Equilibration loop, during which statistics are NOT sampled
void lfts_simulation::equilibrate() {
    int it;
    for (it=1; it<=P_->equil_its(); it++) {

        // Perform a Langevin step with symmetrised noise to update w-(r)
        Langevin_->step_wm(w_.get(), RNG_, P_->XbN(), P_->sigma(), P_->dt());

        // Calculate saddle point value of w+(r), phi-(r) and phi+(r)
        AM_->mix(dbc_.get(),200,1e-4,w_.get());
        std::cout << "lnQ = " << dbc_->calc_concs(w_.get()) << std::endl;

        // Save to file every save_freq_ steps
        if (it%P_->save_freq()==0) { 
            saveStdOutputFile("w_eq_" + std::to_string(it));
            file_IO::saveArray(w_.get()+2*M_, "phi_eq_"+std::to_string(it), 2*M_);
        }
    }
    // Final save to file at end of equilibration period
    saveStdOutputFile("w_eq_" + std::to_string(it-1));
    file_IO::saveArray(w_.get()+2*M_, "phi_eq_"+std::to_string(it-1), 2*M_);
}


// Statistics loop, during which statistics are sampled
void lfts_simulation::statistics() {
    int it;
    for (it=1; it<=P_->sim_its(); it++) {

        // Perform a Langevin step with symmetrised noise to update w-(r)
        Langevin_->step_wm(w_.get(), RNG_, P_->XbN(), P_->sigma(), P_->dt());

        // Calculate saddle point value of w+(r), phi-(r) and phi+(r)
        AM_->mix(dbc_.get(),200,1e-4,w_.get());
        std::cout << "lnQ = " << dbc_->calc_concs(w_.get()) << std::endl;

        // Sample statistics every sample_freq_ steps
        if (it%P_->sample_freq()==0) {
            Sk_->sample(w_.get());
        }

        // Save fields to file every save_freq_ steps
        if (it%P_->save_freq()==0) { 
            saveStdOutputFile("w_st_" + std::to_string(it));
            file_IO::saveArray(w_.get()+2*M_, "phi_st_"+std::to_string(it), 2*M_);
            Sk_->save("struct_st_"+ std::to_string(it));
        }
    }
    // Final save to file at end of equilibration period
    saveStdOutputFile("w_st_" + std::to_string(it-1));
    file_IO::saveArray(w_.get()+2*M_, "phi_st_"+std::to_string(it-1), 2*M_);
    Sk_->save("struct_st_"+ std::to_string(it-1));
}


// Calculate the diblock copolymer Hamiltonian
double lfts_simulation::getH() {
    int r;
    double wp_sum=0.0, wm2_sum=0.0;

    // Calculate the natural log of the partition function
    double lnQ = dbc_->calc_concs(w_.get());

    // Calculate the sums of w+(r) and w-(r)^2
    for (r=0; r<M_; r++) wp_sum += w_[r+M_];
    for (r=0; r<M_; r++) wm2_sum += pow(w_[r],2.0);

    // Return the Hamiltonian
    return -lnQ + (wm2_sum/P_->XbN() - wp_sum)/M_;
}


// Save data in a standard format to be used as in input file
void lfts_simulation::saveStdOutputFile(std::string fileName) {
    P_->saveOutputParams(fileName);
    file_IO::saveArray(w_.get(), fileName, 2*M_, true);
}


void lfts_simulation::generate_field(double *w, int loadType) {
    switch (loadType) {
        case 2:
            field_generator::create_lamellar(w, P_->XbN(), P_->m()); break;
        default:
            // Create a random field with noise of amplitude XN/2
            std::uniform_real_distribution<double> uDist(0.0, 1.0);
            for (int r=0; r<M_; r++) {
                w[r] = P_->XbN()*(uDist(RNG_)-0.5);
                w[r+M_] = 0.0;
            }
            break;
    }
}
