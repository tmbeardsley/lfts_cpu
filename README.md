![Build Status](https://github.com/tmbeardsley/lfts_cpu/actions/workflows/build_and_test.yml/badge.svg)
# Langevin Field-Theoretic Simulation of Diblock Copolymers on CPUs

## Description
See https://www.tbeardsley.com/projects/lfts/fts_gpu for a detailed discussion of this project.<br>

## Required Dependencies
GSL - GNU Scientific Library (https://www.gnu.org/software/gsl/)<br>
FFTW3 Fast Fourier Transform Library (https://www.fftw.org/)<br>

## Compiling
Two methods of compiling the program are available:<br>
<ol>
  <li><b>comp.sh</b>
    <br>
    A simple bash script to create a 'build' directory containing the compiled program code: lfts-cpu.<br><br>
    On a Linux system, run the bash script from the top directory via:<br>
    <b>sh comp.sh</b>
    <br><br>
  </li>
  <li><b>CMake</b>
    <br>
    CMakeLists.txt specifies the required commands for CMake to create (and run) Makefiles, which create a 'build' directory and compile the program code as: lfts-cpu.<br><br>
    From the top directory, run: <br>
    <b>cmake -B build</b><br>
    <b>cmake --build build</b>
  </li>
</ol>


## Running the program
After compilation the executable file, lfts-cpu, resides in the 'build' directory. An input file must be supplied to the executable at the command line, examples of which are contained in the 'input_files' folder. 
For example, from the top level of the directory tree, the program could be run via: <br><br>
<b>./build/lfts-cpu ./input_files/input</b>


## Input Files
The input_files directory contains example input files that can be supplied to the program from the command line.

### Input file format
Line 1: <em>N NA XN C Ndt isXeN</em><br>
Line 2: <em>mx my mz Lx Ly Lz</em><br>
Line 3: <em>n_eq n_st n_smpl save_freq loadType</em><br>
Lines 4->(M+3): W-(r)<br>
Lines (M+4)->(2M+3): w+(r)<br>

Note: A real-space position r = (x,y,z) corresponds to a mesh point position r_m = (i,j,k), where i=0->mx-1, j=0->my-1 and k=0->mz-1 are integers. The elements of the fields, W-(r) and w+(r), are then written in ascending order of the row-major index: p = mx\*(i\*my+j)+k.

#### Parameter descriptions
| Parameter | Description |
| :---: | --- |
| <em>N</em> | number of monomers in a single polymer chain (integer) |
| <em>NA</em> | number of monomers in the A-block of a polymer chain (integer) |
| <em>XN</em> | interaction strength between A and B-type monomers (double) |
| <em>C</em> | square root of the invariant polymerisation index, Nbar (double) |
| <em>Ndt</em> | size of the time step in the Langevin update of W-(r) (double) |
| <em>isXeN</em> | whether the parameter XN is in terms of bare (isXeN=0) or effective (isXeN=1) chi (integer) |
| <em>mx, my, mz</em> | number of mesh points in the x, y, and z dimensions of the simulation box (integers) |
| <em>Lx, Ly, Lz</em> | dimensions of the simulation box (in units of the polymer end-to-end length, R0) in the x, y, and z dimensions (doubles) |
| <em>n_eq</em> | number of langevin steps performed to equilibrate the system (integer) |
| <em>n_st</em> | number of langevin steps performed after equilibration has ended, during which statistics are sampled (integer) |
| <em>n_smpl</em> | number of steps between samples being taken in the statistics period (integer) |
| <em>save_freq</em> | number of steps between saving outputs to file |
| <em>loadType</em> | whether to load the W-(r) and w+(r) fields from the proceeding file lines (loadType=1), start from a disordered state (loadType=0) or start from a (300) lamellar phase (loadType=2) |
| M | Total number of mesh points (M= mx\*my\*mz), such that the proceeding 2*M lines of the file can hold the W-(r) and w+(r) fields that are to be loaded |

## Output files
#### w_eq_<step_number>
The state of the W-(r) and w+(r) fields at simulation step number <step_number> during the equilibration period. First three lines are simulation parameters so it can be used as an input file.<br>

#### w_st_<step_number>
The state of the W-(r) and w+(r) fields at simulation step number <step_number> during the statistics gathering period. First three lines are simulation parameters so it can be used as an input file.<br>

#### phi_eq_<step_number>
The state of the phi-(r) and phi+(r) fields at simulation step number <step_number> during the equilibration period.<br>

#### phi_eq_<step_number>
The state of the phi-(r) and phi+(r) fields at simulation step number <step_number> during the statistics gathering period.<br>

