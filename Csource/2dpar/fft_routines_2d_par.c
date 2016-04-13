/**
 * @file fft_routines_2d_par.c
 * @brief Provides FFT routines for Euler equations
 *
 * @author Claudio Viotti 
 * @date 3/7/13.
 * @note Copyright (c) 2013 Claudio Viotti. All rights reserved.
 */

#include "fft_routines_2d_par.h"

/**
 * @brief FFT routine in two dimensions
 *
 * @param u a double array pointer that contains the function values
 * @param hu a fftw_complex array pointer that contains the Fourier coefficients (afterwards)
 * @param plan a fftw_plan data structure
 */

void fft_2d(double* u, fftw_complex* hu, fftw_plan plan){
    
    
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            /* Copy real data inside f skipping padded elements 
            (note the Ny+2 striding over the first dimension) */
            f[(fNy+2)*i + j]=u[(fNy+2)*i + j];
            
        }
        
    }
    
    
    fftw_execute(plan);
    
    
//    for (i=0; i<local_Nx; i++) {
//    
//        for (j=0; j<Ny/2+1; j++) {
//            
//            hu[(fNy/2+1)*i + j]=hf[(fNy/2+1)*i + j]/Nx/Ny;
//            
//        }
//        
//    }
    
    
    
    for (i=0; i<local_N; i++) {
        
        hu[i] = hf[i]/Nx/Ny;
        
    }
    
}

/**
 * @brief IFFT routine in two dimensions
 *
 * @param hu a fftw_complex array pointer that contains the Fourier coefficients (afterwards)
 * @param u a double array pointer that contains the function values
 * @param plan a fftw_plan data structure
 */

void ifft_2d(fftw_complex* hu, double* u, fftw_plan plan){
    
    
//    for (i=0; i<local_Nx; i++) {
//        
//        for (j=0; j<Ny/2+1; j++) {
//            
//            hf[(Ny/2+1)*i + j]=hu[(Ny/2+1)*i + j];
//            
//        }
//        
//    }

    for (i=0; i<local_N; i++) {
        
        hf[i] = hu[i];
            
    }
    
    
    fftw_execute(plan);
    
    
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            /* Copy real data inside f skipping padded elements 
            (note the Ny+2 striding over the first dimension) */
            u[(Ny+2)*i + j]=f[(Ny+2)*i + j];
            
        }
        
    }
    
}


