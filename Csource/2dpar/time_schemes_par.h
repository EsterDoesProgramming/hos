//
//  time_schemes.h
//  euler
//
//  Created by Claudio Viotti on 3/9/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//


#include "euler_2d_par.h"
#include "fft_routines_2d_par.h"
#include "model_2d_par.h"



#ifndef euler_time_schemes_h
    #define euler_time_schemes_h
#endif

#ifdef USE_DOUBLES
void sol_update_RK(fftw_complex* u_old,TYPE_REAL* t,TYPE_REAL dt,char* dtflag);
#else
void sol_update_RK(fftwf_complex* u_old,TYPE_REAL* t,TYPE_REAL dt,char* dtflag);
#endif

void Setup_TimeScheme(int scheme_flg);
