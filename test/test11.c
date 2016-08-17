/** Test to write poh5 format with poh5.[ch].
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "poh5.h"

#define FILE "poh5_test11.h5"

#define GLEV 4
#define RLEV 0

/*========================================================================*/

int main(){

  phid_t    file_id;   /* identifiers */
  phid_t    g_gid; /* group id for /Grd/Hgrd */
  int   status;

  int      i, j, k, l, ijkl;

  /* data preparation */

  int32_t glevel = GLEV;
  int32_t rlevel = RLEV;
  int32_t grid_topology = 0;
  int is_complete = 1;
  int32_t num_tot_rgn = 10*pow(4,RLEV);
  int32_t gall1d = pow(2,GLEV-RLEV)+2;

  int32_t my_pe = 0;
  int32_t num_of_pe = 1;
  

  int32_t num_of_var = 6; /* tentative; \todo how to set this ? */

  int32_t num_of_rgn=num_tot_rgn;
  int32_t *rgnid_data;
  double  *var_data1; /* var_data[l,ij] */
  float   *var_data2; /* var_data[l,ij] */

  char description[64];
  char note[256];

  char vname[16];
  char vdscr[64]  ;
  char vnote[256] ;
  char vunit[16];
  char layername[16];
  int32_t num_of_layer ;
  int64_t datasize;
  int32_t step;
  int64_t time_start;
  int64_t time_end;

  status=poh5_setup();

  printf("num_tot_rgn=%d\n",num_tot_rgn);
  printf("gall1d=%d\n",gall1d);

  rgnid_data = (int32_t *)malloc( num_tot_rgn * sizeof(int32_t));

  for (i=0;i<num_tot_rgn;i++){
    rgnid_data[i] = i;
  } 

  /* 
   * Start Creating poh5 file.
   */

  /* file_id = H5Fopen(FILE, H5F_ACC_RDWR, H5P_DEFAULT); */
  /* file_id = H5Fcreate(FILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); */
  file_id = poh5_open_file(FILE, POH5_FWRITE);


  /* write global attributes */
  sprintf(description,"%s","description of this file.");
  sprintf(note,"%s","long long note of this file.");
  status = poh5_write_global_attr(file_id, glevel, rlevel, grid_topology, is_complete,
                                  my_pe, num_of_pe,
                                  num_of_rgn, rgnid_data, description, note, num_of_var );



  /*
   * create Grd
   *
   */
  g_gid = poh5_create_hgrid( file_id, num_of_rgn, gall1d );


  /*
   * write one variable 
   *
   */


  {
    /* setup variable */
    sprintf(vname,"%s","var1");
    sprintf(vdscr,"%s","description of var1.");
    sprintf(vnote,"%s","long long note of var1.");
    sprintf(vunit,"%s","unit1");
    sprintf(layername,"%s","ZSALL12");
    num_of_layer = 12;

    datasize = gall1d * gall1d * num_of_layer * num_of_rgn;
    var_data1 = (double *)malloc( datasize * sizeof(double) );
    ijkl = 0;
    for (l=0;l<num_of_rgn;l++){
      for ( k=0;k<num_of_layer;k++){
        for ( j=0;j<gall1d;j++){
          for (i=0;i<gall1d;i++){
            var_data1[ijkl++] = i + j*100 + k*10000 + l*1000000;
          }
        }
      }
    }
  }

  phid_t var_gid = poh5_create_variable(file_id,
                                       vname, vdscr, vnote, vunit, layername, num_of_layer, 
                                       gall1d, POH5_REAL8, num_of_rgn); 


  step = 1;  time_start = 1;  time_end = 10;
  poh5_write_variable_data(var_gid, step, time_start, time_end, POH5_REAL8, var_data1);

  step = 2;  time_start = 11;  time_end = 20;
  poh5_write_variable_data(var_gid, step, time_start, time_end, POH5_REAL8, var_data1);

  step = 3;  time_start = 21;  time_end = 30;
  poh5_write_variable_data(var_gid, step, time_start, time_end, POH5_REAL8, var_data1);


  poh5_close_variable(var_gid);



  /*
   * another variable
   *
   */


  {
    /* setup variable */
    sprintf(vname,"%s","variable2");
    sprintf(vdscr,"%s","another variable2.");
    sprintf(vnote,"%s","another notes for variable2.");
    sprintf(layername,"%s","ZSSFC1");
    num_of_layer = 1;
    datasize = gall1d * gall1d * num_of_layer * num_of_rgn;
    var_data2 = (float *)malloc( datasize * sizeof(float) );
    ijkl = 0;
    for (l=0;l<num_of_rgn;l++){
      for ( k=0;k<num_of_layer;k++){
        for ( j=0;j<gall1d;j++){
          for (i=0;i<gall1d;i++){
            var_data2[ijkl++] = i + j*100 + k*10000 + l*1000000;
          }
        }
      }
    }
  }

  var_gid = poh5_create_variable(file_id,vname,vdscr, vnote, vunit, layername, num_of_layer,
                                 gall1d, POH5_REAL4, num_of_rgn); 

  step = 1;  time_start = 1;  time_end = 30;
  poh5_write_variable_data(var_gid, step, time_start, time_end, POH5_REAL4, var_data2);

  step = 2;  time_start = 31;  time_end = 60;
  poh5_write_variable_data(var_gid, step, time_start, time_end, POH5_REAL4, var_data2);

  step = 3;  time_start = 61;  time_end = 90;
  poh5_write_variable_data(var_gid, step, time_start, time_end, POH5_REAL4, var_data2);

  step = 4;  time_start = 91;  time_end = 120;
  poh5_write_variable_data(var_gid, step, time_start, time_end, POH5_REAL4, var_data2);

  poh5_close_variable( var_gid );


  /* yet another variable but no data */
  /* num_of_layer = 1; */

  var_gid = poh5_create_variable(file_id,
                                 "v3",
                                 "vacant variable", 
                                 "this variable is to check about update_num_of_var()",
                                 "ND",
                                 "ND",
                                 1,
                                 gall1d, POH5_INTEGER4, num_of_rgn); 

  poh5_close_variable( var_gid );


  /* Done */
       
  
  poh5_close_file(file_id);

  printf("...done.\n");


    
  return 0;


}


/**************************************************************************/

