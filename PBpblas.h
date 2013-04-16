/* ---------------------------------------------------------------------
*
*  -- PBLAS routine (version 2.0) --
*     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
*     and University of California, Berkeley.
*     April 1, 1998
*
*  ---------------------------------------------------------------------
*/
/*
*  This file includes PBLAS definitions. All PBLAS routines include this
*  file.
*
*  ---------------------------------------------------------------------
*  #define macro constants
*  ---------------------------------------------------------------------
*/
/*
*  ---------------------------------------------------------------------
*  Function prototypes
*  ---------------------------------------------------------------------
*/
#ifdef __STDC__

void           PB_freebuf_     ( void );

void           picopy_         ( int *,     int *,     int *,
                                 int *,     int *,     int *,
                                 int *,     int *,     int *,
                                 int *,     int * );
void           pscopy_         ( int *,     float *,   int *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     int * );
void           pdcopy_         ( int *,     double *,  int *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     int * );
void           pccopy_         ( int *,     float *,   int *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     int * );
void           pzcopy_         ( int *,     double *,  int *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     int * );

void           psswap_         ( int *,     float *,   int *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     int * );
void           pdswap_         ( int *,     double *,  int *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     int * );
void           pcswap_         ( int *,     float *,   int *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     int * );
void           pzswap_         ( int *,     double *,  int *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     int * );

void           psaxpy_         ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     int * );
void           pdaxpy_         ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     int * );
void           pcaxpy_         ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     int * );
void           pzaxpy_         ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     int * );

void           psscal_         ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pdscal_         ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );
void           pcscal_         ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pcsscal_        ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pzscal_         ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );
void           pzdscal_        ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );

void           psasum_         ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pdasum_         ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );
void           pscasum_        ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pdzasum_        ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );

void           psnrm2_         ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pdnrm2_         ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );
void           pscnrm2_        ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pdznrm2_        ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );

void           psdot_          ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     int * );
void           pddot_          ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     int * );
void           pcdotc_         ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     int * );
void           pcdotu_         ( int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     int * );
void           pzdotc_         ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     int * );
void           pzdotu_         ( int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     int * );

void           psamax_         ( int *,     float *,   int *,
                                 float *,   int *,     int *,
                                 int *,     int * );
void           pdamax_         ( int *,     double *,  int *,
                                 double *,  int *,     int *,
                                 int *,     int * );
void           pcamax_         ( int *,     float *,   int *,
                                 float *,   int *,     int *,
                                 int *,     int * );
void           pzamax_         ( int *,     double *,  int *,
                                 double *,  int *,     int *,
                                 int *,     int * );

void           psgemv_         ( char *,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pdgemv_         ( char *,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );
void           pcgemv_         ( char *,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pzgemv_         ( char *,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );

void           psagemv_        ( char *,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pdagemv_        ( char *,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );
void           pcagemv_        ( char *,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 int * );
void           pzagemv_        ( char *,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 int * );

void           psger_          ( int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int * );
void           pdger_          ( int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int * );
void           pcgerc_         ( int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int * );
void           pcgeru_         ( int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int * );
void           pzgerc_         ( int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int * );
void           pzgeru_         ( int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int * );

void           pssymv_         ( char *,  int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     int * );
void           pdsymv_         ( char *,  int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     int * );
void           pchemv_         ( char *,  int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     int * );
void           pzhemv_         ( char *,  int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     int * );

void           psasymv_        ( char *,  int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     int * );
void           pdasymv_        ( char *,  int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     int * );
void           pcahemv_        ( char *,  int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     int * );
void           pzahemv_        ( char *,  int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     int * );

void           pssyr_          ( char *,  int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int * );
void           pdsyr_          ( char *,  int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int * );
void           pcher_          ( char *,  int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int * );
void           pzher_          ( char *,  int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int * );

void           pssyr2_         ( char *,  int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int * );
void           pdsyr2_         ( char *,  int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int * );
void           pcher2_         ( char *,  int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int * );
void           pzher2_         ( char *,  int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int * );

void           pstrmv_         ( char *,  char *,  char *,
                                 int *,     float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int * );
void           pdtrmv_         ( char *,  char *,  char *,
                                 int *,     double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int * );
void           pctrmv_         ( char *,  char *,  char *,
                                 int *,     float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int * );
void           pztrmv_         ( char *,  char *,  char *,
                                 int *,     double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int * );

void           psatrmv_        ( char *,  char *,  char *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     int * );
void           pdatrmv_        ( char *,  char *,  char *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     int * );
void           pcatrmv_        ( char *,  char *,  char *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     int * );
void           pzatrmv_        ( char *,  char *,  char *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     int * );

void           pstrsv_         ( char *,  char *,  char *,
                                 int *,     float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int * );
void           pdtrsv_         ( char *,  char *,  char *,
                                 int *,     double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int * );
void           pctrsv_         ( char *,  char *,  char *,
                                 int *,     float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int *,
                                 int * );
void           pztrsv_         ( char *,  char *,  char *,
                                 int *,     double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int *,
                                 int * );

void           psgeadd_        ( char *,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int * );
void           pdgeadd_        ( char *,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int * );
void           pcgeadd_        ( char *,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int * );
void           pzgeadd_        ( char *,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int * );

void           psgemm_         ( char *,  char *,  int *,
                                 int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int * );
void           pdgemm_         ( char *,  char *,  int *,
                                 int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int * );
void           pcgemm_         ( char *,  char *,  int *,
                                 int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     float *,   int *,
                                 int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int * );
void           pzgemm_         ( char *,  char *,  int *,
                                 int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     double *,  int *,
                                 int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int * );

void           pssymm_         ( char *,  char *,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int * );
void           pdsymm_         ( char *,  char *,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int * );
void           pcsymm_         ( char *,  char *,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int * );
void           pzsymm_         ( char *,  char *,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int * );
void           pchemm_         ( char *,  char *,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int * );
void           pzhemm_         ( char *,  char *,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int * );

void           pssyr2k_        ( char *,  char *,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int * );
void           pdsyr2k_        ( char *,  char *,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int * );
void           pcsyr2k_        ( char *,  char *,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int * );
void           pzsyr2k_        ( char *,  char *,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int * );
void           pcher2k_        ( char *,  char *,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int * );
void           pzher2k_        ( char *,  char *,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int * );

void           pssyrk_         ( char *,  char *,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int * );
void           pdsyrk_         ( char *,  char *,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int * );
void           pcsyrk_         ( char *,  char *,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int * );
void           pzsyrk_         ( char *,  char *,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int * );
void           pcherk_         ( char *,  char *,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int * );
void           pzherk_         ( char *,  char *,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int * );

void           pstradd_        ( char *,  char *,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int * );
void           pdtradd_        ( char *,  char *,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int * );
void           pctradd_        ( char *,  char *,  int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int * );
void           pztradd_        ( char *,  char *,  int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int * );

void           pstran_         ( int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int * );
void           pdtran_         ( int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int * );
void           pctranc_        ( int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int * );
void           pztranc_        ( int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int * );
void           pctranu_        ( int *,     int *,     float *,
                                 float *,   int *,     int *,
                                 int *,     float *,   float *,
                                 int *,     int *,     int * );
void           pztranu_        ( int *,     int *,     double *,
                                 double *,  int *,     int *,
                                 int *,     double *,  double *,
                                 int *,     int *,     int * );

void           pstrmm_         ( char *,  char *,  char *,
                                 char *,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int * );
void           pdtrmm_         ( char *,  char *,  char *,
                                 char *,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int * );
void           pctrmm_         ( char *,  char *,  char *,
                                 char *,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int * );
void           pztrmm_         ( char *,  char *,  char *,
                                 char *,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int * );

void           pstrsm_         ( char *,  char *,  char *,
                                 char *,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int * );
void           pdtrsm_         ( char *,  char *,  char *,
                                 char *,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int * );
void           pctrsm_         ( char *,  char *,  char *,
                                 char *,  int *,     int *,
                                 float *,   float *,   int *,
                                 int *,     int *,     float *,
                                 int *,     int *,     int * );
void           pztrsm_         ( char *,  char *,  char *,
                                 char *,  int *,     int *,
                                 double *,  double *,  int *,
                                 int *,     int *,     double *,
                                 int *,     int *,     int * );
#else

void           PB_freebuf_     ();
void           PB_topget_      ();
void           PB_topset_      ();

void           picopy_         ();
void           pscopy_         ();
void           pdcopy_         ();
void           pccopy_         ();
void           pzcopy_         ();

void           psswap_         ();
void           pdswap_         ();
void           pcswap_         ();
void           pzswap_         ();

void           psaxpy_         ();
void           pdaxpy_         ();
void           pcaxpy_         ();
void           pzaxpy_         ();

void           psscal_         ();
void           pdscal_         ();
void           pcscal_         ();
void           pcsscal_        ();
void           pzscal_         ();
void           pzdscal_        ();

void           psasum_         ();
void           pdasum_         ();
void           pscasum_        ();
void           pdzasum_        ();

void           psnrm2_         ();
void           pdnrm2_         ();
void           pscnrm2_        ();
void           pdznrm2_        ();

void           psdot_          ();
void           pddot_          ();
void           pcdotc_         ();
void           pcdotu_         ();
void           pzdotc_         ();
void           pzdotu_         ();

void           psamax_         ();
void           pdamax_         ();
void           pcamax_         ();
void           pzamax_         ();

void           psgemv_         ();
void           pdgemv_         ();
void           pcgemv_         ();
void           pzgemv_         ();

void           psagemv_        ();
void           pdagemv_        ();
void           pcagemv_        ();
void           pzagemv_        ();

void           psger_          ();
void           pdger_          ();
void           pcgerc_         ();
void           pcgeru_         ();
void           pzgerc_         ();
void           pzgeru_         ();

void           pssymv_         ();
void           pdsymv_         ();
void           pchemv_         ();
void           pzhemv_         ();

void           psasymv_        ();
void           pdasymv_        ();
void           pcahemv_        ();
void           pzahemv_        ();

void           pssyr_          ();
void           pdsyr_          ();
void           pcher_          ();
void           pzher_          ();

void           pssyr2_         ();
void           pdsyr2_         ();
void           pcher2_         ();
void           pzher2_         ();

void           pstrmv_         ();
void           pdtrmv_         ();
void           pctrmv_         ();
void           pztrmv_         ();

void           psatrmv_        ();
void           pdatrmv_        ();
void           pcatrmv_        ();
void           pzatrmv_        ();

void           pstrsv_         ();
void           pdtrsv_         ();
void           pctrsv_         ();
void           pztrsv_         ();

void           psgeadd_        ();
void           pdgeadd_        ();
void           pcgeadd_        ();
void           pzgeadd_        ();

void           psgemm_         ();
void           pdgemm_         ();
void           pcgemm_         ();
void           pzgemm_         ();

void           pssymm_         ();
void           pdsymm_         ();
void           pcsymm_         ();
void           pchemm_         ();
void           pzsymm_         ();
void           pzhemm_         ();

void           pssyr2k_        ();
void           pdsyr2k_        ();
void           pcsyr2k_        ();
void           pcher2k_        ();
void           pzsyr2k_        ();
void           pzher2k_        ();

void           pssyrk_         ();
void           pdsyrk_         ();
void           pcsyrk_         ();
void           pcherk_         ();
void           pzsyrk_         ();
void           pzherk_         ();

void           pstradd_        ();
void           pdtradd_        ();
void           pctradd_        ();
void           pztradd_        ();

void           pstran_         ();
void           pdtran_         ();
void           pctranc_        ();
void           pctranu_        ();
void           pztranc_        ();
void           pztranu_        ();

void           pstrmm_         ();
void           pdtrmm_         ();
void           pctrmm_         ();
void           pztrmm_         ();

void           pstrsm_         ();
void           pdtrsm_         ();
void           pctrsm_         ();
void           pztrsm_         ();

#endif
