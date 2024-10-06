# FIND PATH AND LIBRARIES OF FFTW3 LIBRARY

find_path(FFTW3_INCLUDE_DIR fftw3.h)
find_library(FFTW3_LIBRARY fftw3)

if(FFTW3_INCLUDE_DIR AND FFTW3_LIBRARY)
    set(FFTW3_FOUND TRUE)
else()
    set(FFTW3_FOUND FALSE)
endif()

if(FFTW3_FOUND)
    # The below variables can be accessed in CMakeLists.txt
    set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
    set(FFTW3_LIBRARIES ${FFTW3_LIBRARY})
endif()