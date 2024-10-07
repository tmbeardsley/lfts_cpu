#!/bin/bash
#$ -cwd

mkdir build
g++ -O2 -std=c++14 ./src/fts_cpu.cc -o ./build/lfts-cpu -lfftw3 -lgsl -lgslcblas -lm