/** Test to read poh5 format with poh5.[ch].
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "poh5.h"

#define INFILE  "poh5_test11.h5"
#define OUTFILE "poh5_test12.h5"

/* #define GLEV 4 */
/* #define RLEV 0 */

void show_usage(char *argv0){
  printf("Usage:\n");
  printf("%% %s varname step\n",argv0);
  exit(1);
}



/*========================================================================*/

int main(int argc, char *argv[]){

  phid_t fid, gid;  /* hdf5 identifiers */



  /* global attributes */
  int32_t glevel;
  int32_t rlevel;
  int32_t grid_topology ; /* 0 for normal ico grid */
  int is_complete; /* 0 or 1 */
  char description[POH5_HMID];
  char note[POH5_HLNG];

  int32_t my_pe;
  int32_t num_of_pe;
  
  int32_t num_of_var;
  int32_t num_of_rgn;
  int32_t *rgnid;

  /* variable attributes */
  char vname[POH5_HSHT]="hoge"; /* specified by command line */
  char vname0[POH5_HSHT]; /* read from file */
  char vdscr[POH5_HMID];
  char vnote[POH5_HLNG];
  char vunit[POH5_HSHT];
  char lname[POH5_HSHT];
  int32_t nlayer ;
  int32_t nsteps = 0;
  int64_t  *vdatai8 = NULL;
  int32_t  *vdatai4 = NULL;
  real64_t *vdatar8 = NULL;
  real32_t *vdatar4 = NULL;
  int64_t *ts = NULL;
  int64_t *te = NULL;

  real64_t *grddata = NULL;

  int32_t num_tot_rgn;
  int32_t gall1d;

  int64_t datasize;
  int32_t step = 0;
  int64_t time_start;
  int64_t time_end;


  
  int ret;

  if ( argc < 2 ) show_usage(argv[0]);
  if ( argc > 2 ) { sscanf(argv[2],"%d",&step); };
  if ( argc > 1 ) { sscanf(argv[1],"%s",vname);}

  printf("vname,step:%s,%d\n",vname,step);


  ret = poh5_setup();
  /* Open existing poh5 file */
  fid = poh5_open_file(INFILE,POH5_FREAD);


  /* read global attributes */
  ret = poh5_read_global_attr(fid,
                                 &glevel, &rlevel, &grid_topology, &is_complete,
                                 &my_pe, &num_of_pe,
                                 &num_of_rgn, &rgnid, description, note, &num_of_var );

  gall1d = pow(2,glevel-rlevel)+2;
  num_tot_rgn = 10*pow(4,rlevel);
  {
    printf("*** Global Attributes:***\n");
    printf("%15s: %d\n","glevel"       , glevel);
    printf("%15s: %d\n","rlevel"       , rlevel);
    printf("%15s: %d\n","grid_topology", grid_topology);
    printf("%15s: %d\n","is_complete"  , is_complete);
    printf("%15s: %d\n","my_pe"        , my_pe);
    printf("%15s: %d\n","num_of_pe"    , num_of_pe);
    printf("%15s: %s\n","description"  , description);
    printf("%15s: %s\n","note"         , note);
    printf("%15s: %d\n","gall1d",gall1d);
    printf("%15s: %d\n","num_of_rgn"   , num_of_rgn);
    printf("%15s: %d\n","num_of_var"  , num_of_var);

    printf("%15s: ","rgnid");
    for(int i=0;i<num_of_rgn;i++){printf("%d ",rgnid[i]);};printf("\n");

  }

  /* \todo below should be independent function() */
  if (    ( (num_tot_rgn == num_of_rgn) &&   is_complete ) 
       || ( (num_tot_rgn != num_of_rgn) && ! is_complete ) ){
      printf("Number of regions and complete flag seems consistent.\n");
  }else{
      printf("Number of regions and complete flag seems NOT consistent.\n");
  }
    

  /* Open and read out hgrid */
  gid = poh5_open_hgrid(fid);
  if ( gid < 0 ) {
          fprintf(stderr,"ERR: Failed to open /Grd/Hgrid.\n");
          exit(1);
  }
  

  datasize = 3 * gall1d * gall1d * num_of_rgn;
  grddata = (real64_t *)malloc(datasize * sizeof(real64_t) );
  ret = poh5_read_hgrid_data(gid,"grd_x",grddata);

  printf("grd_x:\n");
  for(int i=0;i<datasize;i++) printf("%10d,%14.4f\n",i,grddata[i]);



  /* Open variable */
  gid = poh5_open_variable(fid, vname);
  if ( gid < 0 ) {
    printf("ERR: No such variable: %s\n",vname);
    exit(1);
  }

  int dtype; /* HIO_{INTEGER|REAL}{4|8} */
  ret = poh5_read_variable_attr(gid, 
                                   vname0, vdscr, vnote, vunit, lname,
                                   &nlayer, &nsteps,&dtype);

  {
    printf("*** Variable \"%s\": ***\n",vname);
    printf(" %15s: %s\n","vdscr",vdscr);
    printf(" %15s: %s\n","vnote",vnote);
    printf(" %15s: %s\n","vunit",vunit);
    printf(" %15s: %s\n","layername",lname);
    printf(" %15s: %d\n","num_of_layer",nlayer);
    printf(" %15s: %d\n","num_of_steps",nsteps);
    /* printf("%15s: ","time_start"); */
    /* for(int i=0;i<nsteps;i++){printf("%ld ",(ts[i]));};printf("\n"); */
    /* printf("%15s: ","time_end"); */
    /* for(int i=0;i<nsteps;i++){printf("%ld ",(te[i]));};printf("\n"); */
  }

    /* note that step is 1-base, not 0-base */
  if ( (step > nsteps) || (step < 1) ){
    printf("ERR: given step(%d) is out of bounds(%d).\n",(int)step,(int)nsteps);
    exit(1);
  }



  datasize = gall1d * gall1d * nlayer * num_of_rgn;
  switch(dtype){
  case POH5_INTEGER4 : 
    printf(" %15s: %s\n","datatype","INTEGER 4byte\n");
    vdatai4 = (int32_t *)malloc(datasize * sizeof(int32_t) );
    ret = poh5_read_variable_data(gid,  step, &time_start, &time_end, dtype , vdatai4);
    break;
  case POH5_INTEGER8 : 
    printf(" %15s: %s\n","datatype","INTEGER 8byte\n");
    vdatai8 = (int64_t *)malloc(datasize * sizeof(int64_t) );
    ret = poh5_read_variable_data(gid,  step, &time_start, &time_end, dtype, vdatai8);
    break;
  case POH5_REAL4:
    printf(" %15s: %s\n","datatype","REAL 4byte\n");
    vdatar4 = (real32_t *)malloc(datasize * sizeof(real32_t) );
    ret = poh5_read_variable_data(gid,  step, &time_start, &time_end, dtype, vdatar4);
    break;
  case POH5_REAL8:
    printf(" %15s: %s\n","datatype","REAL 8byte\n");
    vdatar8 = (real64_t *)malloc(datasize * sizeof(real64_t) );
    ret = poh5_read_variable_data(gid, step, &time_start, &time_end, dtype, vdatar8);
    break;
  };
    

  printf(" %15s: %d\n","given step",step);
  printf(" %15s: %ld\n","time_start",time_start);
  printf(" %15s: %ld\n","time_end",time_end);

  switch(dtype){
  case POH5_INTEGER4 : 
    for(int i=0;i<datasize;i++) printf("%10d,%14d\n",i,(int)vdatai4[i]);break;
    break;
  case POH5_INTEGER8 : 
    for(int i=0;i<datasize;i++) printf("%10d,%14d\n",i,(int)vdatai8[i]);break;
    break;
  case POH5_REAL4:
    for(int i=0;i<datasize;i++) printf("%10d,%14.4f\n",i,vdatar4[i]);break;
    break;
  case POH5_REAL8:
    for(int i=0;i<datasize;i++) printf("%10d,%14.4f\n",i,vdatar8[i]);break;
    break;

  };


  printf("...done.\n");

  return 0;


}
