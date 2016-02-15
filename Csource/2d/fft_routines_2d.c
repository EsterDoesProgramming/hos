/**
 * @file fft_routines.c
 * @brief Provides FFT routines for Euler equations
 *
 * @author Claudio Viotti 
 * @date 3/7/13.
 * @note Copyright (c) 2013 Claudio Viotti. All rights reserved.
 */

#include "fft_routines_2d.h"


/**
 * @brief FFT routine in two dimensions
 *
 * @param u a double array pointer that contains the function values
 * @param hu a fftw_complex array pointer that contains the Fourier coefficients (afterwards)
 * @param plan a fftw_plan data structure
 */

void fft_2d(double* u, fftw_complex* hu, fftw_plan plan){
    
    int index;
    
    for (int i=0; i<Nx; i++) {
        
        for (int j=0; j<Ny; j++) {
            
            index=Ny*i + j;
    
            f[index]=u[index];
            
        }
        
    }
    
    
    fftw_execute(plan);
    
    
    for (int i=0; i<Nx; i++) {
    
        for (int j=0; j<Ny/2+1; j++) {
            
            index=(Ny/2+1)*i + j;
            
            /* FFTW3 does not include the normalization of the Fourier coefficients */
            hu[index]=hf[index]/Nx/Ny;
            
        }
        
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
    
    int index;
    
    for (int i=0; i<Nx; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
            
            index=(Ny/2+1)*i + j;
            
            hf[index]=hu[index];
            
        }
        
    }

    
    fftw_execute(plan);
    
    
    for (int i=0; i<Nx; i++) {
        
        for (int j=0; j<Ny; j++) {
            
            index=Ny*i + j;
            
            u[index]=f[index];
            u[index]=u[index];
            
        }
        
    }
    
}



/**
 * @brief FFT routine in two dimensions for large arrays
 *
 * @param u a double array pointer that contains the function values
 * @param hu a fftw_complex array pointer that contains the Fourier coefficients (afterwards)
 * @param plan a fftw_plan data structure
 */

void fft_2d_large(double* u, fftw_complex* hu, fftw_plan plan){
    
    int index;
    
    
    for (int i=0; i<Nxl; i++) {
        
        for (int j=0; j<Nyl; j++) {
            
            index=i*Nyl + j;
    
            ff[index]=u[index];
            
        }
        
    }
    
    
    fftw_execute(plan);
    
    
    for (int i=0; i<Nxl; i++) {
    
        for (int j=0; j<Nyl/2+1; j++) {
            
            index=(Nyl/2+1)*i + j;
            
            hu[index]=hff[index]/Nxl/Nyl;
            
        }
        
    }
    
}

/**
 * @brief IFFT routine in two dimensions for large arrays
 *
 * @param hu a fftw_complex array pointer that contains the Fourier coefficients (afterwards)
 * @param u a double array pointer that contains the function values
 * @param plan a fftw_plan data structure
 */

void ifft_2d_large(fftw_complex* hu, double* u, fftw_plan plan){
    
    int index;
    
    
    for (int i=0; i<Nxl; i++) {
        
        for (int j=0; j<Nyl/2+1; j++) {
            
            index=(Nyl/2+1)*i + j;
            
            hff[index]=hu[index];
            
        }
        
    }

    
    fftw_execute(plan);
    
    
    for (int i=0; i<Nxl; i++) {
        
        for (int j=0; j<Nyl; j++) {
            
            index=Nyl*i + j;
            
            u[index]=ff[index];
            u[index]=u[index];
            
        }
        
    }
    
}




