//
//  fft_routines.h
//  euler
//
//  Created by Claudio Viotti on 3/7/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//

#include "euler_2d_par.h"


#ifndef euler_fft_routines_h
#define euler_fft_routines_h
#endif

/* FFT-related global variables */
TYPE_REAL          *f;
#ifdef USE_DOUBLES
fftw_complex    *hf;
fftw_plan       fftp;
fftw_plan       ifftp;
#else
fftwf_complex    *hf;
fftwf_plan       fftp;
fftwf_plan       ifftp;
#endif

ptrdiff_t alloc_local;
ptrdiff_t local_Nx;
ptrdiff_t local_Nyhpo;
ptrdiff_t local_N;
ptrdiff_t local_0_start;
ptrdiff_t local_1_start;
ptrdiff_t i;
ptrdiff_t j;

ptrdiff_t fNx;
ptrdiff_t fNy;

#ifdef USE_DOUBLES
void fft_2d(TYPE_REAL* u, fftw_complex* hu, fftw_plan plan);
void ifft_2d(fftw_complex* hu, TYPE_REAL* u, fftw_plan plan);
void fft_2d_large(TYPE_REAL* u, fftw_complex* hu, fftw_plan plan);
void ifft_2d_large(fftw_complex* hu, TYPE_REAL* u, fftw_plan plan);
#else
void fft_2d(TYPE_REAL* u, fftwf_complex* hu, fftwf_plan plan);
void ifft_2d(fftwf_complex* hu, TYPE_REAL* u, fftwf_plan plan);
void fft_2d_large(TYPE_REAL* u, fftwf_complex* hu, fftwf_plan plan);
void ifft_2d_large(fftwf_complex* hu, TYPE_REAL* u, fftwf_plan plan);
#endif
