/**
 * @file operators.c
 * @brief provides basic routines to use in the HOS euler model
 *
 * @author Claudio Viotti 
 * @date 7/10/13
 * @note Copyright (c) 2013 Claudio Viotti. All rights reserved.
 */

#include "operators.h"

/**
 * @brief Computes the discrete spectral derivative of hu in x-direction
 * @param hu a fftw_complex pointer
 * @param hu_x a fftw_complex pointer
 */

/* Spectral x-derivative */
void Dx(fftw_complex* hu, fftw_complex* hu_x){


    int index;
    double kx;
    

    for (int i=0; i<Nx/2; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
            
            kx = i*Kx0;
            index = (Ny/2+1)*i + j;
            hu_x[index] = I*kx*hu[index];

        }
        
    }


    for (int i=Nx/2; i<Nx; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
            
            kx = -(Nx-i)*Kx0;
            index = (Ny/2+1)*i + j;
            hu_x[index] = I*kx*hu[index];

        }
        
    }

}

/**
 * @brief Computes the discrete spectral derivative of hu in y-direction
 * @param hu a fftw_complex pointer
 * @param hu_y a fftw_complex pointer
 */

/* Spectral y-derivative */
void Dy(fftw_complex* hu, fftw_complex* hu_y){


    int index;
    double ky;
    
    
    for (int i=0; i<Nx; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
            
            ky = j*Ky0;
            index = (Ny/2+1)*i + j;
            hu_y[index] = I*ky*hu[index];

        }
        
    }

}

/**
 * @brief Computes the discrete derivative of hu in z-direction
 * @param hu a fftw_complex pointer
 * @param hu_z a fftw_complex pointer
 */


/* Spectral z-derivative for deep water */
void Dz(fftw_complex* hu, fftw_complex* hu_z){


    int index;
    double kx, ky, kz;


    for (int i=0; i<Nx/2; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
            
            kx = i*Kx0;
            ky = j*Ky0;
            kz = kx*kx + ky*ky;
            kz = sqrt(kz);
    
            index = (Ny/2+1)*i + j;
            hu_z[index] = kz*hu[index];

        }
        
    }


    for (int i=Nx/2; i<Nx; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
            
            kx = -(Nx-i)*Kx0;
            ky = j*Ky0;
            kz = kx*kx + ky*ky;
            kz = sqrt(kz);
    
            index = (Ny/2+1)*i + j;
            hu_z[index] = kz*hu[index];

        }
        
    }

}

/**
 * @brief writes a fftw_complex pointer array into a larger array with zero-padding
 *
 * @param hu a fftw_complex array pointer
 * @param hu_large a fftw_complex array pointer of larger dimension
 */

/* Copy a complex Nx*(Ny/2+1) array into a 3/2*Nx*(3*Ny/4+1) array with zero-padding */
void Extend(fftw_complex* hu, fftw_complex* hu_large){


    int index, index_large;
    
    
    for (int i=0; i<Nx/2; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
            
            index = (Ny/2+1)*i + j;
            index_large = (Nyl/2+1)*i + j;
            hu_large[index_large] = hu[index];

        }
        
        for (int j=Ny/2+1; j<Nyl/2+1; j++) {
            
            index_large = (Nyl/2+1)*i + j;
            hu_large[index_large] = 0.0;

        }
        
    }

    for (int i=Nx/2; i<(Nxl-Nx/2); i++) {
    //for (int i=Nx/2; i<Nx; i++) {
        
        for (int j=0; j<Nyl/2+1; j++) {
            
            index_large = (Nyl/2+1)*i + j;
            hu_large[index_large] = 0.0;

        }
        
    }

    for (int i=Nx/2; i<Nx; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
            
            index = (Ny/2+1)*i + j;
            //index_large = (Nyl/2+1)*(i+Nx/2) + j;
            index_large = (Nyl/2+1)*(i+Nxl-Nx) + j;
            hu_large[index_large] = hu[index];
            
        }
        
        for (int j=Ny/2+1; j<Nyl/2+1; j++) {
            
            //index_large = (Nyl/2+1)*(i+Nx/2) + j;
            index_large = (Nyl/2+1)*(i+Nxl-Nx) + j;
            hu_large[index_large] = 0.0;
            
        }
        
    }


}

/**
 * @brief Chops a complex 3/2*Nx-by-(3*Ny/4+1) array into a Nx-by-(Ny/2+1) array
 *
 * @param hu_large a fftw_complex array pointer
 * @param hu a fftw_complex array pointer of lesser dimension
 *
 * @note This basically reverses "Extend"
 */

/* Chop a complex 3/2*Nx-by-(3*Ny/4+1) array into a Nx-by-(Ny/2+1) array */
void Shrink(fftw_complex* hu_large, fftw_complex* hu){


    int index, index_large;
    
    
    for (int i=0; i<Nx/2; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
            
            index = (Ny/2+1)*i + j;
            index_large = (Nyl/2+1)*i + j;
            hu[index] = hu_large[index_large];

        }
        
    }
    
    for (int i=Nx/2; i<Nx; i++) {
        
        for (int j=0; j<Ny/2+1; j++) {
            
            index = (Ny/2+1)*i + j;
            //index_large = (Nyl/2+1)*(i+Nx/2) + j;
            index_large = (Nyl/2+1)*(i+Nxl-Nx) + j;
            hu[index] = hu_large[index_large];
            
        }
        
    }

    for (int j=0; j<Ny/2+1; j++) {
        
        index = (Ny/2+1)*(Nx/2) + j;
        hu[index] = 0.0;
        //index = (Ny/2+1)*(Nx/2-1) + j;
        //hu[index] = 0.0;
        //index = (Ny/2+1)*(Nx/2+1) + j;
        //hu[index] = 0.0;
    }

    
}


/**
 * @brief computes a multiplication of two vectors using zero padding for dealiasing.
 *
 * @param hu1 a fftw_complex pointer
 * @param hu2 a fftw_complex pointer
 * @param hprod a fftw_complex pointer - the product of hu1 and hu2
 */


/* Compute dealiased product between two arrays. */
/* The operation can be executed fully in-place hu1=hu2=hprod. */
void MultDea(fftw_complex* hu1, fftw_complex* hu2, fftw_complex* hprod){
    
    
    int index;
    
    
    Extend(hu1, htemp1_large);
    ifft_2d_large(htemp1_large, temp1_large, ifftp_large);
    
    Extend(hu2, htemp1_large);
    ifft_2d_large(htemp1_large, temp2_large, ifftp_large);

    for (int i=0; i<Nxl; i++) {
        
        for (int j=0; j<Nyl; j++) {
    
            index = Nyl*i + j;
            temp1_large[index] = temp1_large[index]*temp2_large[index];
    
        }
    
    }
    
    fft_2d_large(temp1_large, htemp1_large, fftp_large);
    Shrink(htemp1_large, hprod);
    
}

/**
 * @brief compute product between two arrays
 *
 * @param hu1 a fftw_complex array pointer
 * @param hu2 a fftw_complex array pointer
 * @param hprod a fftw_complex array pointer - the product of hu1 and hu2
 */


/* Compute product between two arrays. */
/* The operation can be executed fully in-place hu1=hu2=hprod. */
void Mult(fftw_complex* hu1, fftw_complex* hu2, fftw_complex* hprod){
    
    
    int index;
    
    
    ifft_2d(hu1, temp1, ifftp);
    ifft_2d(hu2, temp2, ifftp);


    for (int i=0; i<Nx; i++) {
        
        for (int j=0; j<Ny; j++) {
    
            index = Ny*i + j;
            temp1[index] = temp1[index]*temp2[index];
    
        }
    
    }
    
    fft_2d(temp1, hprod, fftp);
    
}


/**
 * @brief computes the sum of two arrays
 *
 * @param hu1 a fftw_complex array pointer
 * @param hu2 a fftw_complex array pointer
 * @param hsum a fftw_complex array pointer - the sum of hu1 and hu2
 */

void Sum(fftw_complex* hu1, fftw_complex* hu2, fftw_complex* hsum){
    
    
    int index, Nyhal;
    
    Nyhal = Ny/2+1;
    
    for (int i=0; i<Nx; i++) {
        
        for (int j=0; j<Nyhal; j++) {
    
            index = Nyhal*i + j;
            hsum[index] = hu1[index] + hu2[index];
    
        }
        
    }    
    
}
