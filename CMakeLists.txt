# Specify the minimum version of CMake required to build this project.
cmake_minimum_required(VERSION 3.14)

# Define the project name.
project(LFTS-CPU)

# Set the C++ standard to C++17.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)





# FetchContent is a modern cmake module that allows you to automatically download
# and incorporate third-party code into the project
include(FetchContent)
# Declare an external content source (named googletest)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
# After this line the GoogleTest library will be downloaded, configured, and ready to use
FetchContent_MakeAvailable(googletest)




# search path for CMake modules to be loaded by the include() or find_package() commands before checking the default modules that come with CMake.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Use find_package to find and check for the gsl library.
# The REQUIRED keyword makes sure the build process stops if the library is not found.
# Variables GSL_INCLUDE_DIRS and GSL_LDFLAGS are created here.
find_package(GSL REQUIRED)
message(STATUS "GSL_INCLUDE_DIRS = ${GSL_INCLUDE_DIRS}")
message(STATUS "GSL_LDFLAGS = ${GSL_LDFLAGS}")
include_directories(${GSL_INCLUDE_DIRS})

# Use find_package to find and check for the OpenMP library.
find_package(OpenMP REQUIRED)

# Use find_package to find and check for the FFTW3 library.
# The REQUIRED keyword makes sure the build process stops if the library is not found.
# Variables FFTW3_INCLUDE_DIRS and FFTW3_LIBRARIES are created here.
find_package(FFTW3 REQUIRED)
message(STATUS "FFTW3_INCLUDE_DIRS = ${FFTW3_INCLUDE_DIRS}")
message(STATUS "FFTW3_LIBRARIES = ${FFTW3_LIBRARIES}")
include_directories(${FFTW3_INCLUDE_DIRS})

# Define the executable target "lfts-cpu" and specify "fts_cpu.cc" as the source file.
add_executable(lfts-cpu ./src/fts_cpu.cc ./src/strFunc.cc ./src/langevin.cc ./src/step.cc ./src/diblock.cc ./src/anderson.cc ./src/lfts_simulation.cc ./src/lfts_params.cc)

# Link the GSL libraries to the executable. This includes linking necessary runtime libraries found by pkg-config.
target_link_libraries(lfts-cpu ${GSL_LDFLAGS})

# Link the FFTW3 libraries to the executable. This includes linking necessary runtime libraries found by pkg-config.
target_link_libraries(lfts-cpu ${FFTW3_LIBRARIES})


######## OPTION TO COMPILE WITH OPENMP FOR FFTW3 ########
option(USE_OMP "Include OpenMP for threads in fftw3" OFF)
if (USE_OMP)
  message(STATUS "USING OPENMP!")
  # Link OpenMP flags and libraries to the executable.
  target_compile_definitions(lfts-cpu PRIVATE USE_OMP)
  target_link_libraries(lfts-cpu 
    OpenMP::OpenMP_CXX
    -lfftw3_omp
  )
endif()




# Allows CMake to recognize and manage test cases.
# Without this, CMake would not generate the necessary test infrastructure (such as adding tests to the CTest suite).
enable_testing()

# CMake will compile ./tests/file_IO_test.cc and produce an executable file that can be run during the test phase.
add_executable(
  fileIO_tests
  ./tests/file_IO_test.cc
)

# Link the hello_test executable should be linked with GoogleTest's main library.
# The GTest::gtest_main library contains the main function required to run GoogleTest test cases.
target_link_libraries(
  fileIO_tests
  GTest::gtest_main
)






# Include the GoogleTest CMake module, which provides utility functions to make working with GoogleTest easier. 
# Module provides functions like gtest_discover_tests to help discover and organize GoogleTest cases in the project.
include(GoogleTest)

# Helper function from the GoogleTest CMake module.
# Searches for all test cases using the GoogleTest framework within the trial_test executable and adds them to the CMake test suite.
# Can then run CTest, a testing framework integrated with CMake, without having to manually define each test.
gtest_discover_tests(fileIO_tests)