// ######################################################################################
// Exposes public methods to perform an L-FTS simulation: equilibrate() and statistics()
// ######################################################################################

#pragma once
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include "diblock.h"
#include "anderson.h"
#include "strFunc.h"
#include <random>
#include "field_generator.h"
#include "lfts_params.h"
#include "file_IO.h"
#include "langevin.h"
#include <memory>

class lfts_simulation {

    std::unique_ptr<double[]> w_;           // Array containing: N*w-(r), N*w+(r), phi-(r), phi+(r)            
    std::unique_ptr<diblock> dbc_;          // Diblock object for calculating phi-(r) and phi+(r)
    std::unique_ptr<anderson> AM_;          // Anderson mixing object to solve for w+(r)
    std::unique_ptr<langevin> Langevin_;    // Langevin object to update w-(r) at each step
    std::unique_ptr<strFunc> Sk_;           // StrFunc object for dealing with sampling and calculating the structure function
    std::mt19937_64 RNG_;                   // Random number generator
    std::unique_ptr<lfts_params> P_;        // Object to hold the simulation parameters - automatically updates derived parameters

    int M_;                                 // Total number of field mesh points (constant - contained in lfts_params object but copied for tidier code)

    public:

        // Constructor
        lfts_simulation(std::string inputFile);

        // Destructor
        ~lfts_simulation();

        // Equilibration loop, during which statistics are NOT sampled
        void equilibrate();

        // Statistics loop, during which statistics are sampled
        void statistics();

        // Calculate the diblock copolymer Hamiltonian
        double getH();

    private:

        // Save data in a standard format to be used as in input file
        void saveStdOutputFile(std::string fileName);

        void generate_field(double *w, int loadType);

};