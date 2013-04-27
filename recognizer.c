/*******************************************************************************
* FILENAME:    recognizer.c
* DESCRIPTION: main program for PBlacs-enabled character recognizer
* AUTHOR:      James Matthew Welch [JMW]
* SCHOOL:      Arizona State University
* CLASS:       CSE598: High Performance Computing
* INSTRUCTOR:  Dr. Gil Speyer
* SECTION:     20520
* TERM:        Spring 2013
*******************************************************************************/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"
#include "mkl_scalapack.h"
#include "PBblacs.h"
#include "PBpblas.h"

#define MPI_IO 1
#define DEBUG 1
//#define DESCRIPTORS
//#define DOTHEMATH

#define AA(i,j) AA[(i)*M+(j)]

//void Cblacs_get (int context, int request, int* value);
//void Cblacs_gridinfo (int context, int* np_row, int* np_col, int* my_row, int* my_col);
//void Cblacs_pinfo(int*, int*);
//void Cblacs_get(int, int, int*);
//void Cblacs_gridinit(int*, char*, int, int);
//void Cblacs_gridinfo(int, int*, int*, int*, int*);
//void Cblacs_barrier(int , char*);
//void Cblacs_gridexit(int);
int numroc_(int *, int *, int *, int *, int *);

void descinit_(int *, int *, int *, int *, int *,
		int *, int *, int *, int *, int *);

int main(int argc, char **argv) {
	int i, j, k;
	/* initialize timing harness */
	/************  MPI ***************************/
	int myrank_mpi, nprocs_mpi;
	int errorcode=0;
	MPI_Init( &argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);

	/************  BLACS ***************************/
	int context=0, nprow, npcol, myrow, mycol, nb, ndims=2;
#ifdef ORIGINAL	
	int ictxt, nprow, npcol, myrow, mycol, nb;
#endif /* ORIGINAL */
	int info,itemp;
	int ZERO=0,ONE=1; /* constants for passing into pblacs routines */

#ifdef MPI_IO
	/* read in data with MPI-IO */
	char* testFilename = "/scratch/jwelch4/testcharacter.bin";
	char* refFilename  = "/scratch/jwelch4/referenceset.bin";
	MPI_File  refFile;
	MPI_File testFile;
	int fileStatus;

	Cblacs_pinfo( &myrank_mpi, &nprocs_mpi ) ;
	/* processor array dims = sqrt (nprocs)  - assume square integer */
	nprow = npcol = (int)sqrt(nprocs_mpi);

	/* matrix size is already known: 72*(128*128)= 72*16384 */
	int M=72,N=16384;
	nb = 2;	/* number of blocks per side */

	/* create the darray to describe the block distribution */
	MPI_Datatype darray;
	/* number of elements of oldtype in each dimension of global array */
	int gsizes[2] = {M, N};
	/* distribution of array in each dimension */
	int distribs[2] = {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC};
	/* distribution argument in each dimension */
	int dargs[2] = {nb, nb};
	/* size of process grid in each dimension */
	int psizes[2] = {nprow, npcol};

	/* create the darray corresponding to the referenceFile */
	MPI_Type_create_darray(nprocs_mpi, myrank_mpi, ndims, 
		gsizes, distribs, dargs,  psizes,
		MPI_ORDER_FORTRAN, MPI_DOUBLE, &darray);

	MPI_Type_commit(&darray);

#ifdef DEBUG
	printf("P(%d): Begin File read\n", myrank_mpi);
#endif
	fileStatus = MPI_File_open(MPI_COMM_WORLD, refFilename,  MPI_MODE_CREATE, MPI_INFO_NULL, &refFile);
		
		
		
		
		
		
//	fileStatus = MPI_File_open(MPI_COMM_WORLD, testFilename, MPI_MODE_RDONLY, MPI_INFO_NULL, testFile);
#ifdef DEBUG
	printf("P(%d): End File read\n", myrank_mpi);
#endif
#endif /* MPI_IO */
#ifdef DESCRIPTORS
	/* initialize descriptors for each array that is to be shared among the
	 * global process grid
	 * A mxn, U mxm, S mxn, V nxn, T 1xn, X nxn, x 1xn, D nxn  
	 *		m=samples (letters)=72
	 *		n=linear binary file length  = 128x128 = 16384
	*/
	/* determine number of processor rows/columns based on nprocs_mpi */
	Cblacs_get( 0, 0, &context );
	Cblacs_gridinit( &context, "Row", nprow, npcol );
	Cblacs_gridinfo( context, &nprow, &npcol, &myrow, &mycol );

	/* determine the size of the matrix */

	int descA[9],descT[9],descU[9],descS[9],descV[9],descX[9],descx[9],descD[9];

	/* get num rows, columns for each local matrix */
	/* the following is based on netlib.org/scalapack/tools/numroc.f
	 * but it disagrees with the examples in Dr. Speyer's code namely: 
	 *		the third argument to numroc_ should be the coordinate of the processor
	 *		the fifth argument to numroc_ should be nprocs (total num processors) */
	int mA = numroc_( &M, &nb, &myrank_mpi, &ZERO, &nprocs_mpi);
	int nA = numroc_( &N, &nb, &myrank_mpi, &ZERO, &nprocs_mpi);
	int mU = numroc_( &M, &nb, &myrank_mpi, &ZERO, &nprocs_mpi);
	int nU = numroc_( &N, &nb, &myrank_mpi, &ZERO, &nprocs_mpi);
	int mS = numroc_( &M, &nb, &myrank_mpi, &ZERO, &nprocs_mpi);
	int nS = numroc_( &N, &nb, &myrank_mpi, &ZERO, &nprocs_mpi);
	int mV = numroc_( &N, &nb, &myrank_mpi, &ZERO, &nprocs_mpi);
	int nV = numroc_( &N, &nb, &myrank_mpi, &ZERO, &nprocs_mpi);
	int mX = numroc_( &N, &nb, &myrank_mpi, &ZERO, &nprocs_mpi);
	int nX = numroc_( &N, &nb, &myrank_mpi, &ZERO, &nprocs_mpi);
	int mD = numroc_( &N, &nb, &myrank_mpi, &ZERO, &nprocs_mpi);
	int nD = numroc_( &N, &nb, &myrank_mpi, &ZERO, &nprocs_mpi);
	int mx = 1;
	int nx = numroc_( &N, &nb, &myrank_mpi, &ZERO, &nprocs_mpi);
	int mT = 1;
	int nT = numroc_( &N, &nb, &myrank_mpi, &ZERO, &nprocs_mpi);

	/* initialize the descriptors for each global array */
	descinit_(descA, &M,   &N, &nb, &nb, &ZERO, &ZERO, &context, &mA, &info);
	descinit_(descU, &M,   &N, &nb, &nb, &ZERO, &ZERO, &context, &mU, &info);
	descinit_(descS, &M,   &N, &nb, &nb, &ZERO, &ZERO, &context, &mS, &info);
	descinit_(descV, &N,   &N, &nb, &nb, &ZERO, &ZERO, &context, &mV, &info);
	descinit_(descX, &N,   &N, &nb, &nb, &ZERO, &ZERO, &context, &mX, &info);
	descinit_(descD, &N,   &N, &nb, &nb, &ZERO, &ZERO, &context, &mD, &info);
	descinit_(descx, &ONE, &N, &nb, &nb, &ZERO, &ZERO, &context, &mx, &info);
	descinit_(descT, &ONE, &N, &nb, &nb, &ZERO, &ZERO, &context, &mT, &info);

	/* initialize each array with the appropriate contents */
	double *gA = (double*) malloc(M*N*sizeof(double));	/* array to hold the global reference matrix, A */
	double *gT = (double*) malloc(1*N*sizeof(double));	/* array to hold the global test matrix, T */
	double *A = (double*) malloc(mA*nA*sizeof(double));	/* array to hold the local reference matrix, A */
	double *U = (double*) malloc(mU*nU*sizeof(double));	/* array to hold the orthonormal matrix, U */
	double *S = (double*) malloc(mS*nS*sizeof(double));	/* array to hold the singular values matrix, S */
	double *V = (double*) malloc(mV*nV*sizeof(double));	/* array to hold the orthonormal matrix, V */
	double *X = (double*) malloc(mX*nX*sizeof(double));	/* array to hold the matrix, X */
	double *D = (double*) malloc(mD*nD*sizeof(double));	/* array to hold the matrix, D */
	double *x = (double*) malloc(mx*nx*sizeof(double));	/* array to hold the matrix, x */
	double *T = (double*) malloc(mT*nT*sizeof(double));	/* array to hold the test matrix, T */
	
	/* fill the local arrays (reference, test) from the global arrays */
	/* TODO: linker error preventing MPI-IO */
#endif /* DESCRIPTORS */
#ifdef ORIGINAL
	nprow = 2; npcol = 2; nb =2;
	Cblacs_pinfo( &myrank_mpi, &nprocs_mpi ) ;
	Cblacs_get( -1, 0, &ictxt );
	Cblacs_gridinit( &ictxt, "Row", nprow, npcol );
	Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

	/* determine the size of the matrix */
	int M=5;
	double *AA = (double*) malloc(M*M*sizeof(double));
	/* fill the matrix */
	for(i=0;i<M;i++ )
		for(j=0;j<M;j++)
			AA[i*M+j]=(2*i+3*j+1);

	double *X = (double*) malloc(M*sizeof(double));
	X[0]=1;X[1]=1;X[2]=0;X[3]=0;X[4]=1;

	int descA[9],descx[9],descy[9];
	int mA = numroc_( &M, &nb, &myrow, &ZERO, &nprow );
	int nA = numroc_( &M, &nb, &mycol, &ZERO, &npcol );
	int nx = numroc_( &M, &nb, &myrow, &ZERO, &nprow );
	int my = numroc_( &M, &nb, &myrow, &ZERO, &nprow );
	descinit_(descA, &M,   &M,  &nb,  &nb,  &ZERO, &ZERO, &ictxt, &mA, &info);
	descinit_(descx, &M, &ONE,  &nb, &ONE,  &ZERO, &ZERO, &ictxt, &nx, &info);
	descinit_(descy, &M, &ONE,  &nb, &ONE,  &ZERO, &ZERO, &ictxt, &my, &info);
	double *x = (double*) malloc(nx*sizeof(double));
	double *y = (double*) calloc(my,sizeof(double));
	double *A = (double*) malloc(mA*nA*sizeof(double));
	int sat,sut;
	for(i=0;i<mA;i++){
		for(j=0;j<nA;j++){
			sat= (myrow*nb)+i+(i/nb)*nb;
			sut= (mycol*nb)+j+(j/nb)*nb;
			A[j*mA+i]=AA(sat,sut);
		}
	}


	for(i=0;i<nx;i++){
		sut= (myrow*nb)+i+(i/nb)*nb;
		x[i]=X[sut];
	}
#endif
#ifdef DOTHEMATH
	/* calculate mean vector of reference set */
	
	/* subtract the mean from the test image and the reference set (mean center) */

	/* perform svd on the normalized A to get basis matrices, U & V	*/
	/* pdgesvd_() */

	/* Multiply X=U'A  and  x=U'T (corrected from x=UU'T) */

	/* Find mimimum:
	 * Subtract x from each column of X to make D */

	/* caculate sqrt(D'D) */

	/* find the minimum diagonal element, this is the matching character */

	/* print number of matching character */

#endif /* DOTHEMATH */
#ifdef ORIGINAL
	double alpha = 1.0; double beta = 0.0;
	pdgemv_("N",&M,&M,&alpha,A,&ONE,&ONE,descA,x,&ONE,&ONE,descx,&ONE,&beta,y,&ONE,&ONE,descy,&ONE);

	Cblacs_barrier(ictxt,"A");
	for(i=0;i<my;i++)
		printf("rank=%d %.2f \n", myrank_mpi,y[i]);
#endif /* ORIGINAL */
	Cblacs_barrier(context, "A");
	Cblacs_gridexit( 0 );
#ifdef DEBUG
	printf("P(%d) Waiting on barrier\n", myrank_mpi);
#endif
	MPI_Barrier(MPI_COMM_WORLD);
#ifdef MPI_IO
	/* close file with MPI-IO */
	MPI_File_close( &refFile);
//	MPI_File_close(&testFile);
#endif /* MPI_IO */

	MPI_Finalize();
	/* finalize timing harness */
	return 0;
}

/* algorithm for Character-Recognizer-MPI
 * (0) Use MPI-IO to read in the test and reference sets
 * (1) Normalization: Calculate a "mean" vector of all values in the character
 *		set.  Subtract mean from the test image (T) & from each character in 
 *		the set (A).
 * (2) Perform SVD on the normalized A to get basis matrices, U & V
 * (3) Projections: Multiply X=U'A  and  x=U'T (corrected from UU'T)
 * (4) Find Minimum: Subtract x from each column of X to make D,
 *		the "difference matrix"*
 * (5) Calculate DD', find the minimum diagonal element. The column of this
 *		element is the matching character
 */

/* QUESTIONS:  
 *	(3) is the leading dimension the number of rows?? - no. LLD_A=number of rows of the local array that stores the blocks of A
 *	*/
