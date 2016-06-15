//
//  euler_2d.h
//  euler
//
//  Created by Claudio Viotti on 3/7/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//


#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<tgmath.h>
#include<fftw3-mpi.h>
#include "mpi.h"
#include "hdf5.h"

#define H5_DATATYPE "H5T_IEEE_F64LE"
#define SAVE_FILE_BUFSIZE 128
#define RUNTIME_DATA_BUFSIZE 32

#define FFT_TRANSPOSE 1

#define NLevs 2 // Number of levels used in the HOS expansion.
                // Particular cases are:
                // Nlev=0 --> linear waves;
                // Nlev=2 --> Zakharov equations.


///////////////////////
/* Global variables. */

#ifdef USE_DOUBLES
typedef double TYPE_REAL;
#else
typedef float TYPE_REAL;
#endif

/* Main parameters */
int             Nx;
int             Ny;
//ptrdiff_t       N;
int             runsubid;
TYPE_REAL          Lx;
TYPE_REAL          Ly;
TYPE_REAL          Kx0;
TYPE_REAL          Ky0;
TYPE_REAL          g;

TYPE_REAL          T;
TYPE_REAL          dtsave;
int             saveflg;    // =1 -> basic output, >1 -> extended output

#ifdef USE_DOUBLES
fftw_complex*   hetan;
fftw_complex*   hphin;
#else
fftwf_complex*   hetan;
fftwf_complex*   hphin;
#endif

/* MPI related variables */
int         mpi_size;
int         mpi_rank;
MPI_Comm    comm;
MPI_Info    info;


/* Dommermuth ramping */
int         rampflg;
TYPE_REAL      Tramp;


/* Definitions global for debugging pourpose */
hid_t       savefileid;
hid_t       savefileid2;
hid_t       FileID;
char        savefile_buff[SAVE_FILE_BUFSIZE];
char        savefile2_buff[SAVE_FILE_BUFSIZE];
herr_t      status;


/* Global arrays for temporary storage */
TYPE_REAL*         temp1;
TYPE_REAL*         temp2;
#ifdef USE_DOUBLES
fftw_complex*   htemp1;
fftw_complex*   htemp2;
#else
fftwf_complex*   htemp1;
fftwf_complex*   htemp2;
#endif

/* Runtime datafile */
char        runtime_data_buff[RUNTIME_DATA_BUFSIZE];
FILE*       runtime_fid;
TYPE_REAL      Ham;
TYPE_REAL      Ham_glob;

