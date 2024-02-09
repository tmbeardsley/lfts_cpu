#!/bin/bash
#$ -cwd

g++ -O2 -std=c++14 -lgsl -lgslcblas -lfftw3 fts_cpu.cc -o fts_cpu
