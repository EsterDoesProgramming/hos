/**         
 * @file  xdmf_routines_2d_par.c                           
 * @brief Internal routines to generate XML schema for reading hdf5 files with paraview                         
 *                           
 * @author Nicole Beisiegel               
 * @date   06/02/2016.                                     
 */

#include "xdmf_routines_2d_par.h"
#include "hdf5_routines_2d_par.h"

#ifdef USE_DOUBLES
#define TYPE_DATA Float
#else
#define TYPE_DATA Single
#endif


/**                                                                             
 * @brief Initizing XML scheme for h5 output files (needed for processing files with paraview)              
 *         
 * @param
 * @return 
 */

void init_xml(FILE *xmlfile, TYPE_REAL time, int Nx, int Ny, TYPE_REAL Lx, TYPE_REAL Ly){

    int     j, k, indx=0;
    TYPE_REAL  WaveNumbers[Nx*Ny][2]; 
    TYPE_REAL  pi;
    
    pi = (TYPE_REAL) M_PI;

    /* Compute Vectors with wave numbers in X and Y direction and write to HDF5  */
    for (j=0; j<Nx; j++){
      for (k=0; k<Ny; k++){
	WaveNumbers[indx][0] = 2 * k * pi / Ly;
	WaveNumbers[indx][1] = 2 * j * pi / Lx;
     
	indx += 1;
	  }
    }

    FileID = create_file_2d("data_coordinates.h5");                             
                      
    write_coordinates_2d(FileID, WaveNumbers);    
    status = close_file_2d(FileID);  

    /* Create XMF file */
    fprintf( xmlfile,"<Xdmf>\n");  
    fprintf( xmlfile,"<Domain>\n");
    fprintf( xmlfile,"<Grid Name=\"CellTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");


    /* Print all information for zeroth time step */
    fprintf( xmlfile,"<Grid Name=\"Cells\" GridType=\"Uniform\">\n");
    
    fprintf( xmlfile, "<Time Value=\"%f\"/>\n", time);
    fprintf( xmlfile, "<Topology TopologyType=\"Polyvertex\" NodesPerElement=\"%d\">\n",Nx*Ny);
    fprintf( xmlfile, "</Topology>\n");
    
    /* Coordinates of the grid */
    fprintf( xmlfile, "<Geometry GeometryType=\"XY\">\n");

    fprintf( xmlfile, "<DataItem DataType=\"TYPE_DATA\" Dimensions=\"%d 2\" Format=\"HDF\">\n", Nx*Ny);
    fprintf( xmlfile, "data_coordinates.h5:/xy");
    fprintf( xmlfile, "</DataItem>\n");

    fprintf( xmlfile, "</Geometry>\n");

    /* Values on the grid points --> Initial Conditions */
    fprintf( xmlfile, "<Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"Phi\">\n");
    fprintf( xmlfile, "<DataItem DataType=\"TYPE_DATA\" Dimensions=\"%d 1\" Format=\"HDF\">data0.1.h5:/phi</DataItem>\n",Nx*Ny);
    fprintf( xmlfile, "</Attribute>\n");

    fprintf( xmlfile, "<Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"Eta\">\n");
    fprintf( xmlfile, "<DataItem DataType=\"TYPE_DATA\" Dimensions=\"%d 1\" Format=\"HDF\">data0.1.h5:/eta</DataItem>\n",Nx*Ny);
    fprintf( xmlfile, "</Attribute>\n");

    fprintf( xmlfile, "<Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"Bathy\">\n");
    fprintf( xmlfile, "<DataItem DataType=\"TYPE_DATA\" Dimensions=\"%d 1\" Format=\"HDF\">data0.1.h5:/bat</DataItem>\n",Nx*Ny);
    fprintf( xmlfile, "</Attribute>\n");
    
    fprintf( xmlfile, "</Grid>\n");
} 


/**                                                                             
 * @brief Writes temporal
 *         
 * @param
 * @return 
 */

void write_xml(FILE *xmlfile, TYPE_REAL time, char *datafile){

   /* Print all information for zeroth time step */
    fprintf( xmlfile,"<Grid Name=\"Cells\" GridType=\"Uniform\">\n");
    
    fprintf( xmlfile, "<Time Value=\"%f\"/>\n", time);
    fprintf( xmlfile, "<Topology TopologyType=\"Polyvertex\" NodesPerElement=\"%d\">\n",Nx*Ny);
    fprintf( xmlfile, "</Topology>\n");
    
    /* Coordinates of the grid */
    fprintf( xmlfile, "<Geometry GeometryType=\"XY\">\n");

    fprintf( xmlfile, "<DataItem DataType=\"TYPE_DATA\" Dimensions=\"%d 2\" Format=\"HDF\">\n", Nx*Ny);
    fprintf( xmlfile, "data_coordinates.h5:/xy");
    fprintf( xmlfile, "</DataItem>\n");

    fprintf( xmlfile, "</Geometry>\n");

    /* Values on the grid points --> Initial Conditions */
    fprintf( xmlfile, "<Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"Phi\">\n");
    fprintf( xmlfile, "<DataItem DataType=\"TYPE_DATA\" Dimensions=\"%d 1\" Format=\"HDF\">%s:/phi</DataItem>\n",Nx*Ny, datafile);
    fprintf( xmlfile, "</Attribute>\n");

    fprintf( xmlfile, "<Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"Eta\">\n");
    fprintf( xmlfile, "<DataItem DataType=\"TYPE_DATA\" Dimensions=\"%d 1\" Format=\"HDF\">%s:/eta</DataItem>\n",Nx*Ny, datafile);
    fprintf( xmlfile, "</Attribute>\n");

    fprintf( xmlfile, "<Attribute AttributeType=\"Scalar\" Center=\"Node\" Name=\"Bathy\">\n");
    fprintf( xmlfile, "<DataItem DataType=\"TYPE_DATA\" Dimensions=\"%d 1\" Format=\"HDF\">%s:/bat</DataItem>\n",Nx*Ny, datafile);
    fprintf( xmlfile, "</Attribute>\n");
    
    fprintf( xmlfile, "</Grid>\n");

}


/**                                                                             
 * @brief Finalizes and closes the XML scheme
 *         
 * @param
 * @return 
 */

void close_xml(FILE *xmlfile){
  
    fprintf( xmlfile,"</Grid>\n"); 
    fprintf( xmlfile,"</Domain>\n");
    fprintf( xmlfile,"</Xdmf>\n");
    
    fclose(xmlfile);

}

