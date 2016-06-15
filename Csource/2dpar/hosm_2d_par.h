/**
 * @file  hosm_2d_par.h
 * @brief header file for hosm euler code
 *
 * @author Nicole Beisiegel
 * @date   2/17/16
 * @note   Copyright (c) 2016 Nicole Beisiegel. All rights reserved.
 */


#include "euler_2d_par.h"
#include "fft_routines_2d_par.h"
#include "operators_par.h"

#ifndef euler_model_h
#define euler_model_h
#endif

void rhs_hos_setup();
void rhs_hos_clean();
TYPE_REAL RampFun(const TYPE_REAL t);

#ifdef USE_DOUBLES
void rhs_test(fftw_complex* rhs, fftw_complex* u);
void rhs_hos(fftw_complex* rhs, fftw_complex* u, TYPE_REAL t);

void Zvel(fftw_complex* hu, fftw_complex* hZvelM, fftw_complex* hZvelM2, fftw_complex* hZvel2M, fftw_complex* hZvel2M2, TYPE_REAL t);
void ZvelLinear(const fftw_complex* hu, fftw_complex* hZvelLinear);

TYPE_REAL Hamiltonian(const fftw_complex* heta, const fftw_complex* heta_t, const fftw_complex* hphi);
#else
void rhs_test(fftwf_complex* rhs, fftwf_complex* u);
void rhs_hos(fftwf_complex* rhs, fftwf_complex* u, TYPE_REAL t);

void Zvel(fftwf_complex* hu, fftwf_complex* hZvelM, fftwf_complex* hZvelM2, fftwf_complex* hZvel2M, fftwf_complex* hZvel2M2, TYPE_REAL t);
void ZvelLinear(const fftwf_complex* hu, fftwf_complex* hZvelLinear);

TYPE_REAL Hamiltonian(const fftwf_complex* heta, const fftwf_complex* heta_t, const fftwf_complex* hphi);

#endif
