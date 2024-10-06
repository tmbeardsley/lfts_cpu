# Specify the minimum version of CMake required to build this project.
cmake_minimum_required(VERSION 3.5)

# Define the project name.
project(LFTS-CPU)

# Set the C++ standard to C++14.
set(CMAKE_CXX_STANDARD 14)

# search path for CMake modules to be loaded by the include() or find_package() commands before checking the default modules that come with CMake.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Use find_package to find and check for the gsl library.
# The REQUIRED keyword makes sure the build process stops if the library is not found.
# Variables GSL_INCLUDE_DIRS and GSL_LDFLAGS are created here.
find_package(GSL REQUIRED)
message(STATUS "GSL_INCLUDE_DIRS = ${GSL_INCLUDE_DIRS}")
message(STATUS "GSL_LDFLAGS = ${GSL_LDFLAGS}")
include_directories(${GSL_INCLUDE_DIRS})

# Use find_package to find and check for the gsl library.
# The REQUIRED keyword makes sure the build process stops if the library is not found.
# Variables GSL_INCLUDE_DIRS and GSL_LDFLAGS are created here.
find_package(FFTW3 REQUIRED)
message(STATUS "FFTW3_INCLUDE_DIRS = ${FFTW3_INCLUDE_DIRS}")
message(STATUS "FFTW3_LIBRARIES = ${FFTW3_LIBRARIES}")
include_directories(${FFTW3_INCLUDE_DIRS})

# Define the executable target "lfts-cpu" and specify "fts_cpu.cc" as the source file.
add_executable(lfts-cpu ./src/fts_cpu.cc)

# Link the GSL libraries to the executable. This includes linking necessary runtime libraries found by pkg-config.
target_link_libraries(lfts-cpu ${GSL_LDFLAGS})
target_link_libraries(lfts-cpu ${FFTW3_LIBRARIES})