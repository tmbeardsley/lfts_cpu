#!/bin/bash
#$ -cwd

########################################################################
### This script puts output files in a format that can be visualised ###
### by a paraview volume plot (https://www.paraview.org/)            ###
########################################################################

# Check for the correct number of arguments.
if [ "$#" -ne 4 ]
then
  echo "Error -> The following arguments are required..."
  echo "(1) Input file name"
  echo "(2) mx - number of lattice points in x-direction."
  echo "(3) my - number of lattice points in y-direction."
  echo "(4) mz - number of lattice points in z-direction."
  exit 1
fi

# Set number of mesh points in the x, y, z directions from script arguments.
mx=${2}
my=${3}
mz=${4}

# Arbitrarily set the distance between mesh points for paraview.
dx=0.1
dy=0.1
dz=0.1

# Calculate the total number of mesh points.
((M=mx*my*mz))

# Write out parameters to console.
echo 'mx='$mx
echo 'my='$my
echo 'mz='$mz
echo 'M='$M
echo 'dx='$dx
echo 'dy='$dy
echo 'dz='$dz

# Get the name of the input file.
in_file=${1}

# Create an output file name with .vtk extension for paraview.
out_file=${1}.vtk

# Write header information into the output .vtk file.
echo '# vtk DataFile Version 2.0' > ./${out_file}
echo "CT Density" >> ./${out_file}
echo "ASCII" >> ./${out_file}
echo >> ./${out_file}
echo "DATASET STRUCTURED_POINTS" >> ./${out_file}
echo "DIMENSIONS ${mx} ${my} ${mz}" >> ./${out_file}
echo "ORIGIN 0.000000 0.000000 0.000000" >> ./${out_file}
echo "SPACING ${dx} ${dy} ${dz}" >> ./${out_file}
echo >> ./${out_file}
echo "POINT_DATA ${M}" >> ./${out_file}
echo "SCALARS scalars float" >> ./${out_file}
echo "LOOKUP_TABLE default" >> ./${out_file}
echo >> ./${out_file}

# Write the first M lines of the input file to the output file (after the header).
awk -F '\t' '{print $1} NR=='"$M"'{exit}' ${in_file} >> ./${out_file}