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
#include <time.h>
#include <sys/time.h>

#define ALERT
//#define DESCRIPTORS
#define DOTHEMATH
#define READTESTFILE

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

double rowMax(double * u, int rowLength);
double rowMin(double * u, int rowLength);

int main(int argc, char **argv) {
    /* timing variables */
    struct timeval startTime, endTime, ioStart, ioEnd;
    long seconds, useconds;
    double preciseTime;
    /* Get start time - executed by all nodes since no rank assigned yet*/
	gettimeofday(&startTime, NULL);

	/************  MPI ***************************/
	int myrank_mpi, nprocs_mpi;
	int errorcode=0;
	MPI_Init( &argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);

	/************  BLACS ***************************/
	int context=0, nprow, npcol, myrow, mycol, nb, ndims=2;
	int i, j, k;
	int info,itemp;
	int ZERO=0,ONE=1; /* constants for passing into pblacs routines */
	double *A, *U, *S, *V, *X, *D, *x, *T; /* local arrays */

	/* timing harness for IOPS */
	gettimeofday(&ioStart, NULL);

	/* read in data with MPI-IO */
	char* testFilename = "/scratch/jwelch4/testcharacter.bin";
	char* refFilename  = "/scratch/jwelch4/referenceset.bin";
	MPI_File  refFile = NULL;
	MPI_File testFile = NULL;
	int fileStatus;
	MPI_Status * mpistatus=0;

	/* get rank and processor count for the BLACS grid */
	Cblacs_pinfo( &myrank_mpi, &nprocs_mpi ) ;
	/* idefine processor grid topology: 
	 * processor array dims = sqrt (nprocs)  - assume square integer */
	nprow = npcol = (int)sqrt(nprocs_mpi);

	/* matrix size is already known: 72*(128*128)= 72*16384 */
	int M=72,N=16384, nElements;
	nb = 2;	/* number of blocks per side */

	/* set variables for darray */
	MPI_Datatype darray;/* create the darray to describe the block distribution */
	int gsizes[2] = {M, N};/* number of elements of oldtype in each dimension of global array */
	int locsize = 0;
	int distribs[2] = {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC};/* distribution of array in each dimension */
	int dargs[2] = {nb, nb};/* distribution argument in each dimension */
	int psizes[2] = {nprow, npcol};/* size of process grid in each dimension */

	/* create the darray corresponding to the referenceFile */
	MPI_Type_create_darray(nprocs_mpi, myrank_mpi, ndims, 
		gsizes, distribs, dargs,  psizes,
		MPI_ORDER_FORTRAN, MPI_DOUBLE, &darray);
	MPI_Type_commit(&darray);

	fileStatus = MPI_File_open(MPI_COMM_WORLD, refFilename,  MPI_MODE_RDONLY, 
			MPI_INFO_NULL, &refFile);
	char* datarep = "native";
	int mA = numroc_( &M, &nb, &myrank_mpi, &ZERO, &nprocs_mpi);
	int nA = numroc_( &N, &nb, &myrank_mpi, &ZERO, &nprocs_mpi);
	MPI_Type_size(darray, &locsize);
	nElements = locsize / 8; /* doubles are 8 bytes per element */
#ifdef DEBUG
	printf("P(%d): N=%d\n", myrank_mpi, nElements); fflush(stdout);
#endif
	A = malloc(nElements*sizeof(double));  /* array to hold the local reference matrix, A */
#ifdef DEBUG
	printf("P(%d): %d doubles allocated for A.\n", myrank_mpi, nElements);
#endif
	MPI_File_set_view(refFile, 0, MPI_DOUBLE, darray, datarep, MPI_INFO_NULL);
	MPI_File_read_all(refFile, A, nElements, MPI_DOUBLE, mpistatus); 
#ifdef DEBUG
	printf("P(%d): REF read successful \n", myrank_mpi);fflush(stdout);
	printf("P(%d): max(A[0][])=%3.2f\n", myrank_mpi, rowMax(A, nA));
	printf("P(%d): min(A[0][])=%3.2f\n", myrank_mpi, rowMin(A, nA));
#endif

	//MPI_Type_free(&darray); /* clear darray just in case */
	int mT = 1;
	int nT = numroc_( &N, &nb, &myrank_mpi, &ZERO, &nprocs_mpi);

	/* read in test file */
	gsizes[0]=1;	/* only one row in test matrix */

	MPI_Type_create_darray(nprocs_mpi, myrank_mpi, ndims, 
		gsizes, distribs, dargs,  psizes,
		MPI_ORDER_FORTRAN, MPI_DOUBLE, &darray);
	MPI_Type_commit(&darray);
#ifdef DEBUG
	printf("P(%d): darray created successfully\n", myrank_mpi);
#endif
	fileStatus = MPI_File_open(MPI_COMM_WORLD, testFilename,  MPI_MODE_RDONLY, 
			MPI_INFO_NULL, &testFile);
	MPI_Type_size(darray, &locsize);
	nElements = locsize / 8; /* doubles are 8 bytes per element */
#ifdef DEBUG
	printf("P(%d): N=%d\n", myrank_mpi, nElements); fflush(stdout);
#endif
	T = calloc(nElements,sizeof(double));	/* array to hold the test matrix, T */
#ifdef DEBUG
	printf("P(%d): %d doubles allocated for T.\n", myrank_mpi, nElements);
#endif
	MPI_File_set_view(testFile, 0, MPI_DOUBLE, darray, datarep, MPI_INFO_NULL);
#ifdef DEBUG
	printf("P(%d): begin TEST File read\n", myrank_mpi); fflush(stdout);
#endif
	MPI_File_read_all(testFile, T, nElements, MPI_DOUBLE, mpistatus); 
#ifdef DEBUG
	printf("P(%d): End TEST File read\n", myrank_mpi);
	printf("P(%d): nElements=%d, nT=%d\n", myrank_mpi, nElements,nT);
	printf("P(%d): max(T)=%3.2f\n", myrank_mpi, rowMax(T, nElements));
	printf("P(%d): min(T)=%3.2f\n", myrank_mpi, rowMin(T, nElements));
	fflush(stdout);
#endif

	MPI_File_close(&refFile);
	MPI_File_close(&testFile);
#ifdef DEBUG
	printf("P(%d): REF& TEST files closed\n", myrank_mpi);
#endif
	/* close timing harness for IOPS */
    /* timekeeping - only needs to be executed by root*/ 
    if(myrank_mpi == 0) {
		gettimeofday(&ioEnd, NULL);
		seconds  = ioEnd.tv_sec  - startTime.tv_sec;
		useconds = ioEnd.tv_usec - startTime.tv_usec;
		preciseTime = seconds + useconds/1000000.0;
		printf(" IO Time = %3.4f\n", preciseTime );  
		fflush(stdout);
	}

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

#ifdef DESCRIPTORS
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
	U = malloc(mU*nU*sizeof(double));	/* array to hold the orthonormal matrix, U */
	S = malloc(mS*nS*sizeof(double));	/* array to hold the singular values matrix, S */
	V = malloc(mV*nV*sizeof(double));	/* array to hold the orthonormal matrix, V */
	X = malloc(mX*nX*sizeof(double));	/* array to hold the matrix, X */
	D = malloc(mD*nD*sizeof(double));	/* array to hold the matrix, D */
	x = malloc(mx*nx*sizeof(double));	/* array to hold the matrix, x */
	
	/* fill the local arrays (reference, test) from the global arrays */
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

	MPI_Finalize();
	/* finalize timing harness */
    /* timekeeping - only needs to be executed by root*/ 
    if(myrank_mpi == 0) {
		gettimeofday(&endTime, NULL);
		seconds  = endTime.tv_sec  - startTime.tv_sec;
		useconds = endTime.tv_usec - startTime.tv_usec;
		preciseTime = seconds + useconds/1000000.0;
		printf(" Total Time = %3.4f\n", preciseTime );  
		fflush(stdout);
	}
	return 0;
}
/* Local Function Declarations */
double findMaxMag(double** u, int length){
	int i,j;
	double gMax=u[0][0];
	for(j=0; j < length; ++j){
		for(i=0; i<length; ++i){
			if(u[i][j] > gMax){
				gMax = u[i][j];
			}
		}	
	}
	return gMax;
}

double rowMin(double * u, int rowLength){
	int i;
	double minimum=0;
	for (i = 0; i < rowLength; i++) {
		if(u[i] < minimum)
			minimum = u[i];
	}
	return minimum;
}

double rowMax(double * u, int rowLength){
	int i;
	double maximum=0;
	for (i = 0; i < rowLength; i++) {
		if(u[i] > maximum)
			maximum = u[i];
	}
	return maximum;
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
 *	(4) why is T elements getting distributed assymetrically? [8192, 8192, 0, 0] 
 *
 *	*/
