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
#define DEBUG
#define DESCRIPTORS
//#define DOTHEMATH
#define READTESTFILE

#define A(i,j) A[(i)*mA+(j)]
#define Cp(i,j) Cp[(i)*mA+(j)]
#define Z(i,j) Z[(i)*mZ+(j)]
#define x(i) x[i];

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
	int context=0, nprow, npcol, myrow, mycol, ndims=2;
	int nb = 4;	/* number of blocks per side */
	int i, j, k, count=0;
	int info=0,itemp;
	int ZERO=0,ONE=1; /* constants for passing into pblacs routines */
	double *A, *Cp, *T, *U, *UT, *S, *V, *VT, *X, *D, *DT, *x, *y; /* local arrays */
	/* matrix size is already known: 72*(128*128)= 72*16384 */
	int M=16384, N=72, nElements;
#ifdef VERBOSE
	printf("M=%d, N=%d\n", M, N);
#endif

	/* array descriptors for matrices */
	int descA[9],descCp[9],descT[9],descU[9],descUT[9],descS[9],descV[9],descVT[9],descX[9],
		descD[9],descDT[9],descx[9],descy[9];

	/* define processor grid topology: 
	 * Square grid, nprocs x nprocs 
	 * processor array dims = sqrt (nprocs)  - assume square integer */
	nprow = npcol = (int)sqrt(nprocs_mpi);
	/* get rank and processor count for the BLACS grid */
	Cblacs_pinfo( &myrank_mpi, &nprocs_mpi ) ;
	/* establish blacs grid & topology */
	Cblacs_get( 0, 0, &context );
	Cblacs_gridinit( &context, "R", nprow, npcol );
	Cblacs_gridinfo( context, &nprow, &npcol, &myrow, &mycol );
	
#ifdef VERBOSE
	printf("nprow=%d, npcol=%d\n", nprow, npcol);
#endif
	
	/* determine number of processor rows/columns based on nprocs_mpi */
	/* use numroc to get local num rows/columns for each local matrix 
	 * Argument analysis: numroc (n, nb, iproc, isrcproc, nprocs);
	 * 1) n		= num rows/cols of global array
	 * 2) nb	= row/col block size
	 * 3) iproc	= the row/col coord of the processor
	 * 4) isrcproc = the process row/col over which the first row/col of the global matrix is distributed
	 * 5) nprocs = number of rows (nprow) or columns (npcol) in the process grid */
	int mA  = numroc_( &M, &nb, &myrow, &ZERO, &nprow);
	int nA  = numroc_( &N, &nb, &mycol, &ZERO, &npcol);
	int mT  = numroc_( &M, &nb, &myrow, &ZERO, &nprow);
	int nT  = numroc_(&ONE,&nb, &mycol, &ZERO, &npcol);
	int mU  = numroc_( &M, &nb, &myrow, &ZERO, &nprow);
	int nU  = numroc_( &M, &nb, &mycol, &ZERO, &nprow);
	int mUT = numroc_( &M, &nb, &myrow, &ZERO, &nprow);
	int nUT = numroc_( &M, &nb, &mycol, &ZERO, &npcol);
	int mS  = numroc_( &M, &nb, &myrow, &ZERO, &nprow);
	int nS  = numroc_( &N, &nb, &mycol, &ZERO, &npcol);
	int mV  = numroc_( &N, &nb, &myrow, &ZERO, &nprow);
	int nV  = numroc_( &N, &nb, &mycol, &ZERO, &npcol);
	int mVT = numroc_( &N, &nb, &myrow, &ZERO, &nprow);
	int nVT = numroc_( &N, &nb, &mycol, &ZERO, &npcol);
	int mX  = numroc_( &M, &nb, &myrow, &ZERO, &nprow);
	int nX  = numroc_( &N, &nb, &mycol, &ZERO, &npcol);
	int mD  = numroc_( &M, &nb, &myrow, &ZERO, &nprow);
	int nD  = numroc_( &M, &nb, &mycol, &ZERO, &npcol);
	int mx  = numroc_( &M, &nb, &myrow, &ZERO, &nprow);
	int nx  = numroc_(&ONE,&nb, &mycol, &ZERO, &npcol); 
	int my  = numroc_( &N, &nb, &myrow, &ZERO, &nprow);
	int ny  = numroc_(&ONE,&nb, &mycol, &ZERO, &npcol);   

#ifdef VERBOSE
	printf("P(%d):mA=%d,nA=%d  mU=%d,nU=%d  mS=%d,nS=%d  mV=%d,nV=%d  mX=%d,nX=%d  mD=%d,nD=%d  mx=%d,nx=%d  mT=%d,nT=%d\n", 
			myrank_mpi,
		   	mA, nA,
		   	mU, nU,
		   	mS, nS,
		   	mV, nV,
		   	mX, nX,
		   	mD, nD,
		   	mx, nx,
			mT, nT);
	fflush(stdout);
#endif

	/* initialize each array with the appropriate contents */
	U  = malloc(mU * nU * sizeof(double));	/* array to hold the orthonormal matrix, U */
	UT = malloc(mUT* nUT* sizeof(double));	/* array to hold the orthonormal matrix, U */
	S  = malloc(mS * nS * sizeof(double));	/* array to hold the singular values matrix, S, returned as a vector from scalapack */
	V  = malloc(mV * nV * sizeof(double));	/* array to hold the orthonormal matrix, V */
	VT = malloc(nV * mV * sizeof(double));	/* array to hold the orthonormal matrix, VT */
	X  = malloc(mX * nX * sizeof(double));	/* array to hold the matrix, X */
	D  = malloc(mD * nD * sizeof(double));	/* array to hold the matrix, D */
	DT = malloc(nD * mD * sizeof(double));	/* array to hold the matrix, D */
	x  = malloc(mx * nx * sizeof(double));	/* array to hold the matrix, x */
	y  = malloc(my * ny * sizeof(double));	/* array to hold the ones matrix, y */

	/* IO with reference matrix and test matrix */
	/* timing harness for IOPS */
	gettimeofday(&ioStart, NULL);

	/* read in data with MPI-IO */
	char* testFilename = "/scratch/jwelch4/testcharacter.bin";
	char* refFilename  = "/scratch/jwelch4/referenceset.bin";
	MPI_File  refFile = NULL;
	MPI_File testFile = NULL;
	int fileStatus;
	MPI_Status * mpistatus=0;

	/* set variables for darray */
	MPI_Datatype darrayA, darrayT;/* create the darray to describe the block distribution */
	int gsizes[2] = {M, N};/* number of elements of oldtype in each dimension of global array */
	int locsize = 0;
	/* distribution of array in each dimension */
	int distribs[2] = {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC};
	int dargs[2] = {nb, nb};/* distribution argument in each dimension --> square block grid */
	int psizes[2] = {nprow, npcol};/* size of process grid in each dimension */

	/* create the darrayA corresponding to the referenceFile */
	MPI_Type_create_darray(nprocs_mpi, myrank_mpi, ndims, 
		gsizes, distribs, dargs,  psizes,
		MPI_ORDER_FORTRAN, MPI_DOUBLE, &darrayA);
	MPI_Type_commit(&darrayA);

	fileStatus = MPI_File_open(MPI_COMM_WORLD, refFilename,  MPI_MODE_RDONLY, 
			MPI_INFO_NULL, &refFile);
	char* datarep = "native";
	MPI_Type_size(darrayA, &locsize);
	nElements = locsize / 8; /* doubles are 8 bytes per element */
#ifdef VERBOSE
	printf("P(%d): nElements=%d\n", myrank_mpi, nElements); fflush(stdout);
#endif
	A		= malloc(nElements*sizeof(double));  /* array to hold the local reference matrix, A */
	Cp	= malloc(nElements*sizeof(double));  /* array to hold the local copy of ref matrix, Cp */
	MPI_File_set_view(refFile, 0, MPI_DOUBLE, darrayA, datarep, MPI_INFO_NULL);
	MPI_File_read_all(refFile, A, nElements, MPI_DOUBLE, mpistatus); 
#ifdef VERBOSE
	printf("P(%d): REF read successful \n", myrank_mpi);fflush(stdout);
	printf("P(%d): min(A[0][])=%3.2f, max(A[0][])=%3.2f\n",
		   	myrank_mpi, rowMin(A, nA), rowMax(A, nA));
	fflush(stdout);
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	/* read in test file */
	gsizes[0]=1;	/* only one row in test matrix */
	gsizes[1]=M;	/* TODO: this is an odd contradiction - thought this should be a col.mjr */
	/* TODO: above seems inconsistent, but seg-faults if gsizes[1]=1 instead
	 * this seems inconsistent, because A is column-major (samples==columns)
	 * so T should be column major as well.  We'll see when the math is
	 * working */
	MPI_Type_create_darray(nprocs_mpi, myrank_mpi, ndims, 
		gsizes, distribs, dargs,  psizes,
		MPI_ORDER_FORTRAN, MPI_DOUBLE, &darrayT);
	MPI_Type_commit(&darrayT);
	fileStatus = MPI_File_open(MPI_COMM_WORLD, testFilename,  MPI_MODE_RDONLY, 
			MPI_INFO_NULL, &testFile);
	MPI_Type_size(darrayT, &locsize);
#ifdef VERBOSE
	printf("P(%d): Testfile opened successfully\n", myrank_mpi);fflush(stdout);
#endif

	MPI_File_set_view(testFile, 0, MPI_DOUBLE, darrayT, datarep, MPI_INFO_NULL);
	nElements = locsize / 8; /* doubles are 8 bytes per element */
#ifdef VERBOSE
	printf("P(%d): nElements=%d\n", myrank_mpi, nElements); fflush(stdout);
#endif
	/* TODO: see if can switch back to malloc : calloc is slightly slower than malloc */
	T = calloc(nElements,sizeof(double));	/* array to hold the test matrix, T */
#ifdef VERBOSE
	printf("P(%d): %d doubles allocated for T.\n", myrank_mpi, nElements);
	printf("P(%d): begin TEST File read\n", myrank_mpi); fflush(stdout);
#endif

	MPI_File_read_all(testFile, T, nElements, MPI_DOUBLE, mpistatus); 
#ifdef VERBOSE
	printf("P(%d): End TEST File read\n", myrank_mpi);
	printf("P(%d): nElements=%d, mT=%d, nT=%d,  min(T)=%3.2f, max(T)=%3.2f \n", 
			myrank_mpi, nElements, mT, nT, rowMin(T, nElements), rowMax(T, nElements));
	fflush(stdout);
#endif

	MPI_File_close(&refFile);
	MPI_File_close(&testFile);
#ifdef VERBOSE
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

#ifdef DESCRIPTORS
	/* initialize descriptors for each array that is to be shared among the
	 * global process grid
	 * A mxn, U mxm, S mxn, V nxn, T nx1, X mxn, x nx1, D mxn  
	 *		m=samples (letters)=72
	 *		n=linear binary file length  = 128x128 = 16384
	*/
	/* initialize the descriptors for each global array */
	descinit_(descA,  &M,   &N, &nb, &nb, &ZERO, &ZERO, &context, &mA, &info);
	descinit_(descCp,  &M,   &N, &nb, &nb, &ZERO, &ZERO, &context, &mA, &info);
	descinit_(descU,  &M,   &M, &nb, &nb, &ZERO, &ZERO, &context, &mU, &info);
	descinit_(descUT, &M,   &M, &nb, &nb, &ZERO, &ZERO, &context, &mUT, &info);
	descinit_(descS,  &M,   &N, &nb, &nb, &ZERO, &ZERO, &context, &mS, &info);
	descinit_(descV,  &N,   &N, &nb, &nb, &ZERO, &ZERO, &context, &mV, &info);
	descinit_(descVT, &N,   &N, &nb, &nb, &ZERO, &ZERO, &context, &mV, &info);
	descinit_(descX,  &M,   &N, &nb, &nb, &ZERO, &ZERO, &context, &mX, &info);
	descinit_(descD,  &M,   &N, &nb, &nb, &ZERO, &ZERO, &context, &mD, &info);
	descinit_(descDT, &N,   &M, &nb, &nb, &ZERO, &ZERO, &context, &nD, &info);
	descinit_(descx,  &M, &ONE, &nb, &nb, &ZERO, &ZERO, &context, &mx, &info);
	descinit_(descy,  &N, &ONE, &nb, &nb, &ZERO, &ZERO, &context, &my, &info);
	descinit_(descT,  &M, &ONE, &nb, &nb, &ZERO, &ZERO, &context, &mT, &info);
	
#endif /* DESCRIPTORS */
	double alpha = 1.0; double beta = 0.0; double n_reciprocal = 1.0 / N;
	/* Function prototypes: 
	 * Y = aAX  + bY	PvAXPY( N, ALPHA, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )  
	 * C = aAB  + bC	PvGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA, DESCA, B, IB, JB, DESCB, BETA, C, IC, JC, DESCC ) 
	 * C = bC   + aA'	PvTRAN( M, N, ALPHA, A, IA, JA, DESCA, BETA, C, IC, JC, DESCC ) 
	 * Y = aAX  + bY	PvGEMV( TRANS="N", M, N, ALPHA, A, IA, JA, DESCA, X, IX, JX, DESCX, INCX, BETA, Y, IY, JY, DESCY, INCY ) 
	 * Y = aA'X + bY	PvGEMV( TRANS="T", M, N, ALPHA, A, IA, JA, DESCA, X, IX, JX, DESCX, INCX, BETA, Y, IY, JY, DESCY, INCY ) 
	 * */
	
#ifdef DOTHEMATH
	/* calculate mean vector of reference set */
	/* make y a ones matrix for calculating mean column of A */
	for (i = 0; i < my; i++) {
		y[i]=1;
	}
	// sum vectors with: X = AY + 0*X
	// x = 1/72 * A * y + 0 * x
	alpha = 1.0/72; // 1/N for mean
	beta = 0.0; // do not add the additoinal vector 
	pdgemv_("N",&M,&N,&alpha,
			A,&ONE,&ONE,descA,
			y,&ONE,&ONE,descy,&ONE,
			&beta,x,&ONE,&ONE,descx,&ONE);

#ifdef VERYVERBOSE
	if (myrank_mpi==0) {
		for (i = 0; i < mx; i++) {
			printf("%3.3f ", x[i]);
		}
	}
#endif
	
	Cblacs_barrier(context, "A");
	/* (1) subtract the mean from the test image and the reference set (mean center) 
	 * T is the test image
	 * x is the mean of the reference characters */
	for(i=0;i<mT;++i){ 
		//printf("T[%d]=%3.2f  x[%d]=%3.2f\n",i,T[i],i,x[i]);
		T[i] = T[i] - x[i];
	}
	/* Subtract the mean vector from each column of A
	 * inefficient because can't find a scalapack routine for it */
	/* TODO: block the loops for more efficient cache utilization */
	for (j = 0; i < nA; j++) {
		for(i = 0; i < mA; i++)
			A(i,j) = A(i,j) - x(i);
	}

#ifdef DEBUG
	printf("P(%d):mA=%d,nA=%d  mU=%d,nU=%d  mS=%d,nS=%d  mVT=%d,nVT=%d\n", 
			myrank_mpi,
		   	mA, nA,
		   	mU, nU,
		   	mS, nS,
		   	mV, nV);
	fflush(stdout);
#endif
	/* (2) perform svd on the normalized A to get basis matrices, U & V	*/
	/* make a copy of A because it will be overwritten by SVD */
	int lwork = mA*nA;
	double * work = NULL;
#ifdef DEBUG		
	Cblacs_barrier(context, "A");
	printf("P(%d): done with mean centering, begin copy A\n", myrank_mpi);
	fflush(stdout);
#endif
	for (i = 0; i < mA; i++) {
		Cp[i] = A[i];
	}
#ifdef DEBUG		
	Cblacs_barrier(context, "A");
	printf("P(%d): done with copy A, begin SVD\n", myrank_mpi);
	fflush(stdout);
#endif
	/* pdgesvd_() 
	 * PDGESVD(JOBU,JOBVT,M,N,    JOBU, JOBVT = 'V' if want values returned
	 *		A,IA,JA,DESCA,S,
	 *		U,IU,JU,DESCU,
	 *		VT,IVT,JVT,DESCVT,
	 *		WORK,LWORK,INFO)*/
	work =  malloc(lwork * sizeof(double));
	pdgesvd_("V", "N", &M, &N, 
			A, &mA, &nA, descA, S,
			U, &mU, &nU, descU,
			VT,&mV, &nV, descVT,
			work, &lwork, &info); /* when lwork=0, work is ignored */
	printf("P(%d): lwork = %d\n",myrank_mpi,lwork);
#ifdef DEBUG
	printf("P(%d): SVD complete\n", myrank_mpi);
#endif
	/* (3) Projections: Multiply X=U'A  and  x=U'T (corrected from UU'T) */
	// C = alphaA'B+betaC
	// pdgemm (transa, transb, m, n, k, 
	//		alpha, a, ia, ja, desc_a, 
	//		b, ib, jb, desc_b, 
	//		beta, c, ic, jc, desc_c);
	//
	// X=U'A by pdgemm(): X = alphaU'A+betaX (transa="T", transb="N")
	// since transa="T", k=numrows of U
	beta = 0.0;
	pdgemm_("T","N", &M, &N, &mU,
			&alpha, U, &ONE, &ONE, descU,
			A, &ONE, &ONE, descA,
			&beta, X, &ONE, &ONE, descX);
			
	// x=U'T+0*x by pdgemv() 
	pdgemv_("T",&M,&M,
			&alpha,U,&ONE,&ONE,descU,
			T,&ONE,&ONE,descT,&ONE,
			&beta,x,&ONE,&ONE,descx,&ONE);

	/* (4) Find mimimum:
	 * Subtract x from each column of X to make D */
	// D(:,j) = X(:,j) + -1 * x
	/* Subtract x vector from each column of X to form D
	 * inefficient because can't find a scalapack routine for it */
	/* TODO: block the loops for more efficient cache utilization */
	for (j = 0; i < nX; j++) {
		for(i = 0; i < mX; i++)
			D[i*mX+j] = X[i*mX +j] - x(i);
	}

	/* (5) caculate D'D	(corrected from sqrt(D'D) which is unnecessary) 
	 * copy D to another matrix then run pdgemv with transpose the copy of D */
	/* use VT as D'D since it's NxN and don't need it anymore */
	pdgemm_("T","N", &N, &N, &mD,
			&alpha, DT, &ONE, &ONE, descDT,
			D, &ONE, &ONE, descD, 
			&beta, VT, &ONE, &ONE, descVT);

	/* find the minimum diagonal element, this is the matching character */
	int min  = VT[0];
	for (i = 0; i < mV; i++) {
		if( VT[i*mV + i] < min ) {
			min = VT[i*mV + i];
			count = i;
		}
	}

#endif
	/* print number of matching character */
	if(myrank_mpi==0) printf("The Character is number %d  !!!!!    ;-)\n", count);
#ifdef DEBUG
	printf("P(%d) Waiting on barrier\n", myrank_mpi);
#endif
	Cblacs_barrier(context, "A");
	Cblacs_gridexit( 0 );
	/* free allocated arrays */
	free(A);
	free(Cp);
	free(T);
	free(U);
	free(UT);
	free(S); 
	free(V);
	free(VT);
	free(X);
	free(D);
	free(x);
	free(y);

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
	MPI_Finalize();
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
 *	*/
