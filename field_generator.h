// #############################################################################
// Provides functions for creating input fields
// #############################################################################
#pragma once
#include <math.h>

namespace field_generator {

    // Create a lamellar w-(r) field with a cosine wave of size Amplitude.
    // (kx, ky, kz) is the wavevector defining the lamellar orientation. m[] is the mesh size.
    void create_lamellar(double *w, double Amplitude, int *m, int kx=3, int ky=0, int kz=0) {
        int r, M;

        M = m[0]*m[1]*m[2];
        for (int mx=0; mx<m[0]; mx++) {
            for (int my=0; my<m[1]; my++) {
                for (int mz=0; mz<m[2]; mz++) {
                    r = m[2] * (mx*m[1] + my) + mz;	    // Row Major Indexing
                    w[r] = Amplitude * cos(2.0*M_PI*( kx*((double)mx)/m[0] + ky*((double)my)/m[1] + kz*((double)mz)/m[2]));
                    w[r+M] = 0.0;
                }
            }
        }
    }
}