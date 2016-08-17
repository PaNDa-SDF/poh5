/**
 * @file poh5.h
 * @brief Header file for poh5 API.
 * @author Takahiro INOUE @ RIST
 * @date 2016/08/12
 */


#ifndef __POH5_H__
#define __POH5_H__


#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* #include "hdf5.h" */

#define POH5_VERSION 94

#ifndef POH5_DO_SHUF
#define POH5_DO_SHUF
#endif

#ifndef POH5_DO_GZIP
#define POH5_DO_GZIP 6
#endif

/*
 * \note below must be consistent with the one onf hio.h
 *
 */

/* character length */
#define POH5_HLNG 256
#define POH5_HMID  64
#define POH5_HSHT  16

/* data type */
#define POH5_REAL4     0
#define POH5_REAL8     1
#define POH5_INTEGER4  2
#define POH5_INTEGER8  3

/* flie open mode */
#define POH5_FREAD   0
#define POH5_FWRITE  1
#define POH5_FAPPEND 2

/* must be consistent with the ones defined in NICAM sources. */
typedef float (real32_t);
typedef double(real64_t);

/* must be consistent with hid_t defined in H5Ipoublic.h */
typedef int phid_t;


/*========================================================================*/
extern int poh5_setup();


/*========================================================================*/
extern phid_t poh5_open_file( 
                             const char *fname, /**< [in] file name */
                             const int mode);   /**< [in] file open mode */


/*========================================================================*/
extern int poh5_close_file(
                           const phid_t fid); /**< [in] file id */


/*========================================================================*/
extern int poh5_write_global_attr(
                    const phid_t    file_id,      /**< [in] poh5 file_id */
                    const int32_t  glevel,        /**< [in] glevel */
                    const int32_t  rlevel,        /**< [in] rlevel */
                    const int32_t  grid_topology, /**< [in] grid_topology */
                    const int      is_complete,   /**< [in] 1 if complete file */
                    const int32_t  my_pe,         /**< [in] my pe. */
                    const int32_t  num_of_pe,     /**< [in] num of pe. */
                    const int32_t  num_of_rgn,    /**< [in] num of region this file contains. */
                    const int32_t *rgnid,         /**< [in] array of region ids */
                    const char    *description,   /**< [in] description of this file */
                    const char    *note,          /**< [in] longer note of this file */
                    const int32_t  num_of_var);   /**< [in] number of data in this file. */


/*========================================================================*/
extern int poh5_read_global_attr(
                    const phid_t file_id,    /**< [in] poh5 file_id */
                    int32_t *glevel,         /**< [out] glevel */
                    int32_t *rlevel,         /**< [out] rlevel */
                    int32_t *grid_topology,  /**< [out] grid_topology */
                    int     *is_complete,    /**< [out] 1 if complete file */
                    int32_t *my_pe,          /**< [in] my pe. */
                    int32_t *num_of_pe,      /**< [in] num of pe. */
                    int32_t *num_of_rgn,     /**< [out] num of region this file contains. */
                    int32_t *rgnid[],        /**< [out] array of region id's */
                    char    description[64], /**< [out] description of this file */
                    char    note[256] ,      /**< [out] longer note of this file */ 
                    int32_t *num_of_var);    /**< [out] number of data in this file. */


/*========================================================================*/
extern phid_t poh5_create_variable(
                    const phid_t file_id, /**< [in] poh5 file id*/
                    const char *name,     /**< [in] variable name */
                    const char *dscr,     /**< [in] variable description */
                    const char *note,     /**< [in] variable long note */
                    const char *unit,     /**< [in] variable unit */
                    const char *lname,    /**< [in] name of vertical layer */
                    const int32_t nlayer, /**< [in] number of layer */
                    const int gall1d,     /**< [in] gall1d */
                    const int dtype,      /**< [in] POH5_{REAL|INTEGER}{4|8} */
                    const int num_of_rgn);/**< [in] number of regions */


/*========================================================================*/
extern void poh5_write_variable_data(
                    const phid_t var_gid,     /**< gid for variable opened by open_variable().*/
                    const int32_t step,       /**< current step???? */
                    const int64_t time_start, /**< start time this data represents */
                    const int64_t time_end,   /**< end time this data represents */
                    const int dtype,          /**< [in] POH5_{REAL|INTEGER}{4|8} */
                    const void *var_data);    /**< variable data at current step */


/*========================================================================*/
extern void poh5_write_variable_data_1rgn(
                    const phid_t v_gid, /**<[in] poh5 file id .*/
                    const int32_t step,      /**< step count */
                    const int32_t ll, /**< rgn number, 0-based */
                    const int64_t time_start, /**< start time this data represents */
                    const int64_t time_end,   /**< end time this data represents */
                    const int dtype,      /**< [in] POH5_{REAL|INTEGER}{4|8} */
                    const void *var_data);   /**< variable data at current step */


/*========================================================================*/
extern phid_t poh5_open_variable(
                    const phid_t file_id, /**<[in] poh5 file id .*/
                    const char *vname);   /**<[in] variable name */


/*========================================================================*/
extern phid_t poh5_open_variable_by_idx(
                    const phid_t file_id, /**<[in] poh5 file id .*/
                    const int idx);      /**<[in] variable index */


/*========================================================================*/
extern int poh5_close_variable(
                    const phid_t v_gid); /**<[in] group id .*/
  

/*========================================================================*/
extern int poh5_read_variable_attr(
                    const phid_t var_gid, /**< [in] group id of variable */
                    char *name,       /**< [out] variable name  */
                    char *dscr,       /**< [out] variable description */
                    char *note,          /**< [out] variable long note */
                    char *unit,          /**< [out] variable unit */
                    char *lname,         /**< [out] name of vertical layer */
                    int32_t *nlayer,     /**< [out] number of layers */
                    int32_t *nsteps,     /**< [out] number of steps */
                    int *dtype);          /**< [out] POH5_{REAL4,REAL8,INTEGER4,INTEGER8} */


/*========================================================================*/
extern int poh5_read_variable_time( 
                    const phid_t v_gid, /**< [in] group id of variable */
                    int64_t ts[],         /**< [out] start time of this step */
                    int64_t te[]);         /**< [out] end time of this step */



/*========================================================================*/
extern int poh5_read_variable_data(
                    const phid_t var_gid, /**< [in] group id of variable */
                    const int32_t step,        /**<  [in] step counter */
                    int64_t *time_start, /**< [out] start time of this step */
                    int64_t *time_end,   /**< [out] end time of this step */
                    const int dtype,     /**< [in] POH5_{REAL4,REAL8,INTEGER4,INTEGER8} */
                    void *vdata);          /**< [out] data */


/*========================================================================*/
extern int poh5_create_hgrid(
                    const phid_t file_id, /**< [in] poh5 file id */
                    const int num_of_rgn,
                    const int gall1d);


/*========================================================================*/
extern phid_t poh5_open_hgrid(
                    const int file_id); /**< [in] poh5 file id */


/*========================================================================*/
extern int poh5_write_hgrid_data(
                    const int v_gid, /**[in] group_id for /Grd/Hgrid */
                    const char* dname,  /**< [in] "grd_x" or "grd_xt" */
                    const void *var_data);   /**< hgrid data at current step */


/*========================================================================*/
extern int poh5_read_hgrid_data(
                    const int v_gid, /**[in] group_id for /Grd/Hgrid */
                    const char* dname,  /**< [in] "grd_x" or "grd_xt" */
                    void *gdata);   /**< hgrid data at current step */


#endif /* __POH5_H__ */
