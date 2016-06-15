/**
 * @file  hdf5_routines_2d_par.c
 * @brief Parallelized HDF5 routines for the HOS Euler model
 *
 * @author Claudio Viotti, Nicole Beisiegel
 * @date   06/24/2013.
 */

#include "hdf5_routines_2d_par.h"
#include "fft_routines_2d_par.h"


void get_params(char* filename){
    
    hid_t       h5_file, h5_data;
    herr_t      status;
    double      Nd[1];
    
    
    h5_file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    
    /* Nx */
    h5_data = H5Dopen(h5_file, "/Nx", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);
    
    Nx = (int) Nd[0];
    
    /* Ny */
    h5_data = H5Dopen(h5_file, "/Ny", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);
    
    Ny = (int) Nd[0];
    
    /* Lx */
    h5_data = H5Dopen(h5_file, "/Lx", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);

    Lx = (TYPE_REAL) Nd[0];
    
    /* Ly */
    h5_data = H5Dopen(h5_file, "/Ly", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);

    Ly = (TYPE_REAL) Nd[0];
    
    /* g  */
    h5_data = H5Dopen(h5_file, "/g", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);

    g = (TYPE_REAL) Nd[0];
    
    /* runsubid */
    h5_data = H5Dopen(h5_file, "/runsubid", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);

    runsubid = (int) Nd[0];

    /* T */
    h5_data = H5Dopen(h5_file, "/T", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);

    T = (TYPE_REAL) Nd[0];
    
    /* dtsave */
    h5_data = H5Dopen(h5_file, "/dtsave", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);

    dtsave = (TYPE_REAL) Nd[0];
    
    /* saveflg */
    h5_data = H5Dopen(h5_file, "/saveflg", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);
    
    saveflg = (int) Nd[0];
    
    /* rampflg */
    h5_data = H5Dopen(h5_file, "/rampflg", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);
    
    rampflg = (int) Nd[0];
    
    if ( rampflg!=0 && rampflg!=1 ) {
    
        printf("Illegal value assigned to rampflag, setting to default value (0).\n");
        rampflg = 0;
    
    }
    
    /* Tramp */
    h5_data = H5Dopen(h5_file, "/Tramp", H5P_DEFAULT);
    status  = H5Dread(h5_data, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, Nd);

    Tramp = (TYPE_REAL) Nd[0];
   
    status = H5Dclose(h5_data);
    status = H5Fclose(h5_file);
    
    Kx0 = 2*M_PI/Lx;
    Ky0 = 2*M_PI/Ly;
    
}


void get_ic_2d(char* filename, TYPE_REAL* u1, TYPE_REAL* u2, TYPE_REAL* u3){
    
    hid_t       h5_file, h5_dataset, h5_memspace, h5_filespace, FORMAT_H5T_IEEE;
    herr_t      status;
    hid_t       plist_id;
    hsize_t     dimsf[2], count[2], offset[2];    
    
    #ifdef USE_DOUBLES
    FORMAT_H5T_IEEE = H5Tcopy(H5T_IEEE_F64LE);
    #else
    FORMAT_H5T_IEEE = H5Tcopy(H5T_IEEE_F32LE);
    #endif


    plist_id = H5Pcreate(H5P_DATASET_XFER);
    h5_file  = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    
    dimsf[0] = Nx;
    dimsf[1] = Ny;
    count[0] = dimsf[0]/mpi_size;
    count[1] = dimsf[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    
    
    h5_dataset = H5Dopen(h5_file, "/eta0", H5P_DEFAULT);
    
    h5_memspace = H5Screate_simple(RANK, count, NULL);
    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    /* Read data and write it inside a temporary array, then copy the data inside the padded target array */
    status = H5Dread(h5_dataset, FORMAT_H5T_IEEE, h5_memspace, h5_filespace, H5P_DEFAULT, temp1);
    
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
	  u1[(Ny+2)*i + j] = (TYPE_REAL) temp1[Ny*i + j];
        
        }
        
    }
    
    h5_dataset = H5Dopen(h5_file, "/phi0", H5P_DEFAULT);
    
    h5_memspace = H5Screate_simple(RANK, count, NULL);
    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    /* Read data and write it inside a temporary array, then copy the data inside the padded target array */
    status = H5Dread(h5_dataset, FORMAT_H5T_IEEE, h5_memspace, h5_filespace, H5P_DEFAULT, temp1);
    
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
	  u2[(Ny+2)*i + j]= (TYPE_REAL) temp1[Ny*i + j];
            
        }
        
    }
    h5_dataset = H5Dopen(h5_file, "/bat0", H5P_DEFAULT);
    
    h5_memspace = H5Screate_simple(RANK, count, NULL);
    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    /* Read bathymetry  data and write it inside a temporary array, then copy the data inside the padded target array */
    status = H5Dread(h5_dataset, FORMAT_H5T_IEEE, h5_memspace, h5_filespace, H5P_DEFAULT, temp1);
    
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
	  u3[(Ny+2)*i + j] = (TYPE_REAL) temp1[Ny*i + j];
            
        }
        
    }
    
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    H5Fclose(h5_file);
    
}


hid_t create_file_2d(char* filename){
    
    hid_t       h5_file;
    hid_t       plist_id;
    
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);
    
    h5_file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
        
    return h5_file;
    
}


hid_t open_file_2d(char* filename){
    
    hid_t       h5_file;
    
    h5_file = H5Fopen( filename, H5F_ACC_RDWR, H5P_DEFAULT );
    
    return h5_file;
    
}


herr_t close_file_2d(hid_t file){
    
    herr_t      status;
     
    status = H5Fclose(file);
    
    return status;
    
}


void write_header_2d(hid_t file, TYPE_REAL time){
    
    hid_t       data, fid, FORMAT_H5T_IEEE;
    herr_t      status;
    hsize_t     fdim [2];
    TYPE_REAL      var;
    
    #ifdef USE_DOUBLES
    FORMAT_H5T_IEEE = H5Tcopy(H5T_IEEE_F64LE);
    #else
    FORMAT_H5T_IEEE = H5Tcopy(H5T_IEEE_F32LE);
    #endif

    /* For some reason we get a deadlock of only one process executes this part */
    //if (mpi_rank==0) {
  
        fdim[0]=1;
        fdim[1]=1;
        
        fid = H5Screate_simple(2, fdim, NULL);
    
        var = time;

        data = H5Dcreate(file, "time", FORMAT_H5T_IEEE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(data, FORMAT_H5T_IEEE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&var);

        status = H5Dclose(data);

        var = g;
        data = H5Dcreate(file, "g", FORMAT_H5T_IEEE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(data, FORMAT_H5T_IEEE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&g);
        status = H5Dclose(data);

        var = (TYPE_REAL) Nx;
        data = H5Dcreate(file, "Nx", FORMAT_H5T_IEEE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(data, FORMAT_H5T_IEEE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&var);
        status = H5Dclose(data);
        
        var = (TYPE_REAL) Ny;
        data = H5Dcreate(file, "Ny", FORMAT_H5T_IEEE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(data, FORMAT_H5T_IEEE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&var);
        status = H5Dclose(data);
        
        data = H5Dcreate(file, "Lx", FORMAT_H5T_IEEE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(data, FORMAT_H5T_IEEE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Lx);
        status = H5Dclose(data);
        
        data = H5Dcreate(file, "Ly", FORMAT_H5T_IEEE, fid, H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(data, FORMAT_H5T_IEEE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Ly);
        status = H5Dclose(data);
        
        status = H5Sclose(fid);


    //}
    
}


void write_field_2d(hid_t h5_file, TYPE_REAL* eta, TYPE_REAL* phi, TYPE_REAL* bat){
    

    hid_t       h5_dataset, h5_memspace, h5_filespace, FORMAT_H5T_IEEE;
    herr_t      status;
    hid_t       plist_id;
    hsize_t     dimsf[2], count[2], offset[2];
    
    #ifdef USE_DOUBLES
    FORMAT_H5T_IEEE = H5Tcopy(H5T_IEEE_F64LE);
    #else
    FORMAT_H5T_IEEE = H5Tcopy(H5T_IEEE_F32LE);
    #endif

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    dimsf[0] = Nx;
    dimsf[1] = Ny;
    count[0] = dimsf[0]/mpi_size;
    count[1] = dimsf[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    
    
    h5_filespace = H5Screate_simple(RANK, dimsf, NULL);
    h5_dataset = H5Dcreate(h5_file, "eta", FORMAT_H5T_IEEE, h5_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(h5_filespace);

    h5_memspace = H5Screate_simple(RANK, count, NULL);

    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    /* Copy data inside a temporary array to remove padding, then write inside the datafile */
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            temp1[Ny*i + j] = eta[(Ny+2)*i + j];
            
        }
        
    }
    
    status = H5Dwrite(h5_dataset, FORMAT_H5T_IEEE, h5_memspace, h5_filespace, plist_id, temp1);
   
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    
    
    h5_filespace = H5Screate_simple(RANK, dimsf, NULL);
    h5_dataset = H5Dcreate(h5_file, "phi", FORMAT_H5T_IEEE, h5_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(h5_filespace);
    
    h5_memspace = H5Screate_simple(RANK, count, NULL);

    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    /* Copy data inside a temporary array to remove padding, then write inside the datafile */
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            temp1[Ny*i + j] = phi[(Ny+2)*i + j];
            
        }
        
    }
    
    status = H5Dwrite(h5_dataset, FORMAT_H5T_IEEE, h5_memspace, h5_filespace, plist_id, temp1);
       
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    

    h5_filespace = H5Screate_simple(RANK, dimsf, NULL);
    h5_dataset = H5Dcreate(h5_file, "bat", FORMAT_H5T_IEEE, h5_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(h5_filespace);

    h5_memspace = H5Screate_simple(RANK, count, NULL);

    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    /* Copy data inside a temporary array to remove padding, then write inside the datafile */
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            temp1[Ny*i + j] = bat[(Ny+2)*i + j];
            
        }
        
    }
    
    status = H5Dwrite(h5_dataset, FORMAT_H5T_IEEE, h5_memspace, h5_filespace, plist_id, temp1);
   
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    

    
}

void write_coordinates_2d(hid_t h5_file, TYPE_REAL* coords){

    hid_t       h5_dataset, h5_memspace, h5_filespace, FORMAT_H5T_IEEE;
    herr_t      status;
    hid_t       plist_id;
    hsize_t     dimsf[2], count[2], offset[2];
  
    #ifdef USE_DOUBLES
    FORMAT_H5T_IEEE = H5Tcopy(H5T_IEEE_F64LE);
    #else
    FORMAT_H5T_IEEE = H5Tcopy(H5T_IEEE_F32LE);
    #endif  
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    dimsf[0] = Nx*Ny;
    dimsf[1] = 2;

    count[0] = dimsf[0]/mpi_size;
    count[1] = dimsf[1];

    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    
    /* Create HDF5 file */
    h5_filespace = H5Screate_simple(RANK, dimsf, NULL);
    h5_dataset = H5Dcreate(h5_file, "xy", FORMAT_H5T_IEEE, h5_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(h5_filespace);

    h5_memspace = H5Screate_simple(RANK, count, NULL);

    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    plist_id = H5Pcreate(H5P_DATASET_XFER); 
    status = H5Dwrite(h5_dataset, FORMAT_H5T_IEEE, h5_memspace, h5_filespace, plist_id, coords);

    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id); 

}

void write_extra_2d(hid_t h5_file, TYPE_REAL* array1, TYPE_REAL* array2){
    

    hid_t       h5_dataset, h5_memspace, h5_filespace, FORMAT_H5T_IEEE;
    herr_t      status;
    hid_t       plist_id;
    hsize_t     dimsf[2], count[2], offset[2];
  
    #ifdef USE_DOUBLES
    FORMAT_H5T_IEEE = H5Tcopy(H5T_IEEE_F64LE);
    #else
    FORMAT_H5T_IEEE = H5Tcopy(H5T_IEEE_F32LE);
    #endif
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    dimsf[0] = Nx;
    dimsf[1] = Ny;
    count[0] = dimsf[0]/mpi_size;
    count[1] = dimsf[1];
    offset[0] = mpi_rank * count[0];
    offset[1] = 0;
    
    
    h5_filespace = H5Screate_simple(RANK, dimsf, NULL);
    h5_dataset = H5Dcreate(h5_file, "Array1", FORMAT_H5T_IEEE, h5_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(h5_filespace);

    h5_memspace = H5Screate_simple(RANK, count, NULL);

    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    /* Copy data inside a temporary array to remove padding, then write inside the datafile */
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            temp1[Ny*i + j] = array1[(Ny+2)*i + j];
            
        }
        
    }
    
    status = H5Dwrite(h5_dataset, FORMAT_H5T_IEEE, h5_memspace, h5_filespace, plist_id, temp1);
   
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    
    
    h5_filespace = H5Screate_simple(RANK, dimsf, NULL);
    h5_dataset = H5Dcreate(h5_file, "Array2", FORMAT_H5T_IEEE, h5_filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(h5_filespace);
    
    h5_memspace = H5Screate_simple(RANK, count, NULL);

    h5_filespace = H5Dget_space(h5_dataset);
    H5Sselect_hyperslab(h5_filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    
    /* Copy data inside a temporary array to remove padding, then write inside the datafile */
    for (i=0; i<local_Nx; i++) {
        
        for (j=0; j<Ny; j++) {
            
            temp1[Ny*i + j] = array2[(Ny+2)*i + j];
            
        }
        
    }
    
    status = H5Dwrite(h5_dataset, FORMAT_H5T_IEEE, h5_memspace, h5_filespace, plist_id, temp1);
       
    H5Dclose(h5_dataset);
    H5Sclose(h5_filespace);
    H5Sclose(h5_memspace);
    H5Pclose(plist_id);
    
    
}

