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


void rhs_test(fftw_complex* rhs, fftw_complex* u);
void rhs_hos_setup();
void rhs_hos_clean();
void rhs_hos(fftw_complex* rhs, fftw_complex* u, double t);


void Zvel(fftw_complex* hu, fftw_complex* hZvelM, fftw_complex* hZvelM2, fftw_complex* hZvel2M, fftw_complex* hZvel2M2, double t);
void ZvelLinear(const fftw_complex* hu, fftw_complex* hZvelLinear);

double Hamiltonian(const fftw_complex* heta, const fftw_complex* heta_t, const fftw_complex* hphi);
double RampFun(const double t);
