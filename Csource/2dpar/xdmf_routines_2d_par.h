//
//  xdmf_routines_2d_par.h
//  euler
//
//  Created by Nicole Beisiegel on 06/02/16.
//

#include"euler_2d_par.h"

#ifndef euler_xdmf_routines_2d_h
#define euler_xdmf_routines_2d_h

#define DATASETNAME "ExtendibleArray"
#define RANK         2

#endif


void init_xml(FILE *xmlfile, TYPE_REAL time, int Nx, int Ny, TYPE_REAL Lx, TYPE_REAL Ly);
void write_xml(FILE *xmlfile, TYPE_REAL time, char *datafile);
void close_xml(FILE *xmlfile);
