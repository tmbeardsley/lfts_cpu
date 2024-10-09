// #######################################################################
// Creates a new lfts_simulation(...) object, passing in input file to 
// its contructor and subsequently accessing its public methods to 
// equilibrate the system, gather statistics from the equilibrated system
// and finally output the system's energy
// #######################################################################

#include <string>
#include <iostream>
#include "lfts_simulation.h"
#include <chrono>

using namespace std;

int main (int argc, char *argv[])
{
    // Get input file name from command-line argument
    if (argc != 2) {
        cout << "Please supply a single input file name as a command-line argument." << endl << endl;
        exit(1);
    }
    string inputFile(argv[1]);


    cout << "Test message" << endl;

    // New lfts_simulation object with input file name specified
    lfts_simulation *lfts_sim = new lfts_simulation(inputFile);
    
    // Time the equilibration period
    cout << "Starting Equilibration..." << endl;
    auto start = std::chrono::steady_clock::now();
    lfts_sim -> equilibrate();
    auto end = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "Equlibration time = " << duration.count() << "secs" << endl << endl;


    // Time the statistics gathering period
    cout << "Starting Statistics..." << endl;
    start = std::chrono::steady_clock::now();
    lfts_sim -> statistics();
    end = std::chrono::steady_clock::now();
    duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "Statistics time = " << duration.count() << "secs" << endl << endl;

    // Output the final energy of the system
    cout.precision(6);
    cout << lfts_sim->getH() << endl;

    delete lfts_sim;
    return 0;
}
