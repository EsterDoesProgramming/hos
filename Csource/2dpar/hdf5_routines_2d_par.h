//
//  hdf5_routines_2d.h
//  euler
//
//  Created by Claudio Viotti on 6/24/13.
//  Copyright (c) 2013 Claudio Viotti. All rights reserved.
//

#include"euler_2d_par.h"

#ifndef euler_hdf5_routines_2d_h
#define euler_hdf5_routines_2d_h

#define DATASETNAME "ExtendibleArray"
#define RANK         2


#endif


void get_params(char* filename);
void get_ic_2d(char* filename, TYPE_REAL* u1, TYPE_REAL* u2, TYPE_REAL* u3);
hid_t create_file_2d(char* filename);
hid_t open_file_2d(char* filename);
herr_t close_file_2d(hid_t file);
void write_header_2d(hid_t file, TYPE_REAL time);
void write_field_2d(hid_t file, TYPE_REAL* eta, TYPE_REAL* phi, TYPE_REAL* bat);
void write_extra_2d(hid_t file, TYPE_REAL* array1, TYPE_REAL* array2);
