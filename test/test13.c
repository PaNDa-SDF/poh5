/** Test to read poh5 format with poh5.[ch].
 * Read out ALL variables in the file
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "poh5.h"


void show_usage(char *argv0){
  printf("Usage:\n");
  printf("%% %s poh5_file\n",argv0);
  exit(1);
}



/*========================================================================*/

int main(int argc, char *argv[]){


  phid_t    fid,gid;   /* poh5 identifiers */

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
  char vname[POH5_HSHT]="hoge";
  char vdscr[POH5_HMID]  ;
  char vnote[POH5_HLNG] ;
  char vunit[POH5_HSHT];
  char lname[POH5_HSHT];
  int32_t nlayer ;
  int32_t nsteps;
  int64_t  *vdatai8 = NULL;
  int32_t  *vdatai4 = NULL;
  real64_t  *vdatar8 = NULL;
  real32_t  *vdatar4 = NULL;
  int64_t *ts = NULL;
  int64_t *te = NULL;

  int64_t time_start;
  int64_t time_end;

  real64_t *grddata = NULL;

  int32_t num_tot_rgn;
  int32_t gall1d;
  int64_t datasize;

  int dtype; /* POH5_{INTEGER|REAL}{4|8} */

  int ret;
  

  if ( argc != 2 ) show_usage(argv[0]);

  printf("filename:%s\n",argv[1]);


  ret = poh5_setup();
  /* Open existing poh5 file */
  fid = poh5_open_file(argv[1], POH5_FREAD);

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
    printf("%15s: %d\n","num_of_var"   , num_of_var);

    printf(" %15s: ","rgnid");
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
  for ( int n = 0; n < num_of_var; n++ ){
          gid = poh5_open_variable_by_idx(fid, n);
          if ( gid < 0 ) {
                  printf("ERR: No such variable: %s\n",vname);
                  exit(1);
          }
          poh5_read_variable_attr( gid, 
                                   vname, vdscr, vnote, vunit,
                                   lname, &nlayer, &nsteps, &dtype);
    
          {
                  printf("*** Variable no.%d: ***\n",(int)n);
                  printf(" %15s: %s\n","vname",vname);
                  printf(" %15s: %s\n","vdscr",vdscr);
                  printf(" %15s: %s\n","vnote",vnote);
                  printf(" %15s: %s\n","vunit",vunit);
                  printf(" %15s: %s\n","layername",lname);
                  printf(" %15s: %d\n","num_of_layer",nlayer);
                  printf(" %15s: %d\n","num_of_steps",nsteps);
          }
  }
  printf("...done.\n");
  
  return 0;


}
