#!/bin/bash
#$ -cwd

mkdir build

if [ "$1" == "USE_OMP" ]; then
    echo "Using OpenMP."
    echo "Set the desired number of threads (n) in your environment before running the compiled program via: export OMP_NUM_THREADS=n."
    g++ -O3 -std=c++14 ./src/fts_cpu.cc ./src/strFunc.cc ./src/langevin.cc ./src/step.cc ./src/diblock.cc ./src/anderson.cc ./src/lfts_simulation.cc ./src/lfts_params.cc -o ./build/lfts-cpu -fopenmp -lfftw3 -lfftw3_omp -lgsl -lgslcblas -lm -DUSE_OMP
else
    g++ -O3 -std=c++14 ./src/fts_cpu.cc ./src/strFunc.cc ./src/langevin.cc ./src/step.cc ./src/diblock.cc ./src/anderson.cc ./src/lfts_simulation.cc ./src/lfts_params.cc -o ./build/lfts-cpu -lfftw3 -lgsl -lgslcblas -lm
fi



