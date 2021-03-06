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
/* Function prototypes: 
	 * Y = aAX  + bY	PvAXPY( N, ALPHA, X, IX, JX, DESCX, INCX, Y, IY, JY, DESCY, INCY )  
	 * C = aAB  + bC	PvGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA, DESCA, B, IB, JB, DESCB, BETA, C, IC, JC, DESCC ) 
	 * C = bC   + aA'	PvTRAN( M, N, ALPHA, A, IA, JA, DESCA, BETA, C, IC, JC, DESCC ) 
	 * Y = aAX  + bY	PvGEMV( TRANS="N", M, N, ALPHA, A, IA, JA, DESCA, X, IX, JX, DESCX, INCX, BETA, Y, IY, JY, DESCY, INCY ) 
	 * Y = aA'X + bY	PvGEMV( TRANS="T", M, N, ALPHA, A, IA, JA, DESCA, X, IX, JX, DESCX, INCX, BETA, Y, IY, JY, DESCY, INCY ) */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "mpi.h"
#include "mkl_scalapack.h"
#include "PBblacs.h"
#include "PBpblas.h"
#include <time.h>
#include <sys/time.h>

//#define DEBUG

#define A(i,j) A[(i)*mA+(j)]
#define Cp(i,j) Cp[(i)*mA+(j)]
#define Z(i,j) Z[(i)*mZ+(j)]
#define x(i) x[i];

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
	int nb = 2;	/* number of blocks per side */
	int i, j, k, idx, count=0;
	int info=0,itemp;
	int ZERO=0,ONE=1,NEG=-1; /* constants for passing into pblacs routines */
	double *A, *Cp, *T, *U, *UT, *S, *V, *VT, *X, *D, *DT, *x, *y; /* local arrays */
	double ScalarOne[1] = {1};
	/* matrix size is already known: 72*(128*128)= 72*16384 */
	int M=16384, N=72, nElements;
#ifdef VERBOSE
	printf("M=%d, N=%d\n", M, N);
#endif

	/* array descriptors for matrices */
	int descA[9],descCp[9],descT[9],descU[9],descUT[9],descS[9],descV[9],descVT[9],descX[9],
		descD[9],descDT[9],descx[9],descy[9], descScalarOne[9];

	/* define processor grid topology: 
	 * Square grid, nprocs x nprocs 
	 * processor array dims = sqrt (nprocs)  - assume square integer */
	if(nprocs_mpi == 4)
		nprow = npcol = 2;
	else if(nprocs_mpi == 16)
		nprow = npcol = 4;
	else if(nprocs_mpi == 64)
		nprow = npcol = 16;
	else 
		nprow = npcol = 1;

	if (myrank_mpi==0) {
		printf("Character Recognizer with SVD!\n");
		printf("Using %d processors in %d x %d grid\n", nprocs_mpi,nprow,npcol);
	}
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
	int mO  = numroc_(&ONE,&nb, &myrow, &ZERO, &nprow);
	int nO  = numroc_(&ONE,&nb, &mycol, &ZERO, &npcol);   


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
	U  = calloc(mU * nU , sizeof(double));	/* array to hold the orthonormal matrix, U */
	UT = calloc(mUT* nUT, sizeof(double));	/* array to hold the orthonormal matrix, U */
	S  = calloc(N , sizeof(double));	/* array to hold the singular values matrix, S, returned as a vector from scalapack */
	V  = calloc(mV * nV , sizeof(double));	/* array to hold the orthonormal matrix, V */
	VT = calloc(nV * mV , sizeof(double));	/* array to hold the orthonormal matrix, VT */
	X  = calloc(mX * nX , sizeof(double));	/* array to hold the matrix, X */
	D  = calloc(mD * nD , sizeof(double));	/* array to hold the matrix, D */
	DT = calloc(nD * mD , sizeof(double));	/* array to hold the matrix, D */
	x  = calloc(mx * nx , sizeof(double));	/* array to hold the matrix, x */
	y  = calloc(my * ny , sizeof(double));	/* array to hold the ones matrix, y */

	/* initialize work array for SVD later */
	int lwork = -1;//mA*nA;
	double * work = calloc(mA*nA,sizeof(double));

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
	A	= calloc(nElements,sizeof(double));  /* array to hold the local reference matrix, A */
	Cp	= calloc(nElements,sizeof(double));  /* array to hold the local copy of ref matrix, Cp */
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
#ifdef DEBUG
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

#ifdef VERBOSE
	/* attempt to print "ASCII art" of the character */
	int mZ=128;
	double * Z = calloc(mZ*mZ,sizeof(double));
	for (i = 0; i < mZ*mZ; i++) {
		Z[i]=1.0;
	}

	pdlacp3_( &mZ, &ZERO, T, descT, Z, &mZ, &NEG, &NEG, &ZERO);

	Cblacs_barrier(context, "A");
	printf("P(%d): End pdlacp3 \n", myrank_mpi);
	printf("P(%d): mT=%d, nT=%d,  min(T)=%3.2f, max(T)=%3.2f \n", 
			myrank_mpi, mT, nT, rowMin(T, mT*nT), rowMax(T, mT*nT));
	fflush(stdout);
	if(myrank_mpi == 0){	
		printf("\n\n");
		for (i = 0; i < mZ*mZ; i++) {
			if( Z[i] > 0 ) {
				printf("X");
			}else{
				printf(".");
			}
			if((i+1)%128==0)
				printf("\n");
		}
		printf("\n\n");
		fflush(stdout);
	}

	Cblacs_barrier(context, "A");
	MPI_Abort(MPI_COMM_WORLD, errorcode);
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
	 * A mxn, U mxm, S mxn, V nxn, T nx1, X mxn, x nx1, D mxn  
	 *		m=samples (letters)=72
	 *		n=linear binary file length  = 128x128 = 16384
	*/
	/* initialize the descriptors for each global array */
	descinit_(descA,  &M,   &N, &nb, &nb, &ZERO, &ZERO, &context, &mA, &info);
	descinit_(descCp, &M,   &N, &nb, &nb, &ZERO, &ZERO, &context, &mA, &info);
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
	descinit_(descScalarOne, &ONE, &ONE, &nb, &nb, &ZERO, &ZERO, &context, &ONE, &info);
	//if(info){printf("P(%d): Descriptor Init fail, code(%d)\n", myrank_mpi, info);}
	
	double alpha = 1.0; double beta = 0.0; double n_reciprocal = 1.0 / N;
	/* calculate mean vector of reference set */
	/* make y a ones matrix for calculating mean column of A */
	for (i = 0; i < my; i++) {
		y[i]=1;
	}
	// sum vectors with: X = AY + 0*X
	// x = 1/72 * A * y + 0 * x
	alpha = 1.0/N; // 1/N for mean
	beta = 0.0; // do not add the additional vector 
	/* should &mA,&nA be &M,&N instead?? */
	pdgemv_("N",&mA,&nA,
			&alpha,A,&ONE,&ONE,descA,
			       y,&ONE,&ONE,descy,&ONE,
			&beta, x,&ONE,&ONE,descx,&ONE);

	/* (1) subtract the mean from the test image and the reference set (mean center) */
	// sum vectors with: C = AY + 0*X
	//* C = aAB  + bC	PvGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA, DESCA, B, IB, JB, DESCB, BETA, C, IC, JC, DESCC ) 
	// A = -1 * x * y + 1.0 * A
	// 
	alpha=-1.0;
	beta = 1.0;
	/* should &mA,&nA be &M,&N instead?? */
	pdgemm_("N","T", &M , &N , &ONE,
			&alpha, x, &ONE, &ONE, descx,
			y, &ONE, &ONE, descy,
			&beta, A, &ONE, &ONE, descA);
	/* T is the test image
	 * x is the mean of the reference characters */
	// mean center T with: T = -1*Ay + 1.0*T
	// x = -1.0/72 * A * y + 1.0 * T
	alpha = -1.0/72;
	beta = 1.0;
	pdgemv_("N",&M,&ONE,
			&alpha,A,&ONE,&ONE,descA,
			       y,&ONE,&ONE,descy,&ONE,
			&beta, T,&ONE,&ONE,descT,&ONE);
#ifdef DEBUG
	printf("p(%d): Mean Centering Complete\n", myrank_mpi);
	Cblacs_barrier(context, "A");
#endif

	Cblacs_barrier(context, "A");
	/* (2) perform svd on the normalized A to get basis matrices, U & V	*/
	/*PDGESVD(JOBU,JOBVT,M,N,    JOBU, JOBVT = 'V' if want values returned
	 *		A,IA,JA,DESCA,S,		all of IX,IX should be &ONE typically
	 *		U,IU,JU,DESCU,
	 *		VT,IVT,JVT,DESCVT,
	 *		WORK,LWORK,INFO)*/
	 /*workopt is a temp variable to hold size of the work array */
	double wkopt;
	/* SVD */
	pdgesvd_("V", "V", &M, &N, 
			A, &ONE, &ONE, descA, 
			S,
			U, &ONE, &ONE, descU,
			VT,&ONE, &ONE, descVT,
			&wkopt, &lwork, &info); /* init lwork=-1, to get size of work*/
	lwork = (int)wkopt;
#ifdef DEBUG
	printf("P(%d): First SVD complete: lwork = %d, wkopt = %3.2f, work-len=%d\n", myrank_mpi, lwork, wkopt, mA*nA);
	fflush(stdout);
#endif
	work = calloc(lwork , sizeof(double));

#ifdef VERBOSE
	printf("P(%d): done with malloc work, begin SVD\n", myrank_mpi);
	printf("P(%d): lwork = %d\n",myrank_mpi,lwork); 
	fflush(stdout);
	Cblacs_barrier(context, "A");
#endif

	/* copy A back after SVD overwrites it */
	for ( i = 0; i < nElements; i++) {
		A[i] = Cp[i];
	}
	/* SVD */
	pdgesvd_("V", "V", &M, &N, 
			A, &ONE, &ONE, descA, 
			S,
			U, &ONE, &ONE, descU,
			VT,&ONE, &ONE, descVT,
			work, &lwork, &info);
#ifdef DEBUG
	printf("P(%d): Second SVD complete\n", myrank_mpi);
#ifdef VERBOSE
	if(myrank_mpi == 0){
		for (i = 0; i < N; i++) {
			printf("P(%d): S[%d]=%2.4f\n",myrank_mpi,i, S[i]);
		}
	}
#endif
	fflush(stdout);
	Cblacs_barrier(context, "A");
#endif

	/* copy A back after SVD overwrites it */
	for ( i = 0; i < nElements; i++) {
		A[i] = Cp[i];
	}
	Cblacs_barrier(context, "A");
	/* (3) Projections: Multiply X=U'A  and  x=U'T (corrected from UU'T) */
	// C = alphaA'B+betaC
	// pdgemm (transa, transb, m, n, k, 
	//		alpha, a, ia, ja, desc_a, 
	//		b, ib, jb, desc_b, 
	//		beta, c, ic, jc, desc_c);
	//
	// X=U'A by pdgemm(): X = alphaU'A+betaX (transa="T", transb="N")
	// since transa="T", k=numrows of U
	alpha = 1.0;
	beta = 0.0;
	pdgemm_("T","N", &M, &N, &mU,
			&alpha, U, &ONE, &ONE, descU,
			A, &ONE, &ONE, descA,
			&beta, X, &ONE, &ONE, descX);
			
	Cblacs_barrier(context, "A");
	// x=U'T+0*x by pdgemv() 
	pdgemv_("T",&M,&M,
			&alpha,U,&ONE,&ONE,descU,
			T,&ONE,&ONE,descT,&ONE,
			&beta,x,&ONE,&ONE,descx,&ONE);

#ifdef DEBUG
	printf("P(%d): projections complete\n", myrank_mpi);
	fflush(stdout);
#endif

	Cblacs_barrier(context, "A");
	/* (4) Find mimimum:
	 * Subtract x from each column of X to make D */
	// D(:,j) = X(:,j) - x
	/* Subtract x vector from each column of X to form D
	 * inefficient because can't find a scalapack routine for it */
	/* TODO: block the loops for more efficient cache utilization */
	// sum vectors with: C = AY + 0*X
	//* C = aAB  + bC	PvGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, IA, JA, DESCA, B, IB, JB, DESCB, BETA, C, IC, JC, DESCC ) 
	// X = -1 * x * y + 1.0 * X
	alpha=-1.0;
	beta = 1.0;
	/* should &mA,&nA be &M,&N instead?? */
	pdgemm_("N","T", &M, &N, &ONE,
			&alpha, x, &ONE, &ONE, descx,
					y, &ONE, &ONE, descy,
			&beta,	X, &ONE, &ONE, descX);
//	/* copy X into D	*/
//	for (i = 0; i < mD*nD; i++) {
//		D[i] = X[i];
//	}
#ifdef DEBUG
	printf("P(%d): Difference Matrix calculated\n", myrank_mpi);
	fflush(stdout);
#endif

	Cblacs_barrier(context, "A");
	/* (5) caculate D'D	(corrected from sqrt(D'D) which is unnecessary) 
	 * copy D to another matrix then run pdgemv with transpose the copy of D */
	/* use VT as D'D since it's NxN and don't need it anymore */
	alpha=1.0; beta=0.0;
	/* TODO: originally &M, &N, &mD  was seg-faulting, try alternate*/
	pdgemm_("T","N", &N, &N, &N,
			&alpha, X, &ONE, &ONE, descX,
			        X, &ONE, &ONE, descX, 
			&beta, VT, &ONE, &ONE, descVT);

#ifdef DEBUG
	printf("P(%d): D'D calculated\n", myrank_mpi);
	fflush(stdout);
#endif

	Cblacs_barrier(context, "A");
	/* find the minimum diagonal element, this is the matching character */
	/* TODO: how do we do this across all nodes? 
	 * PDLACP3( M, I, A, DESCA, B, LDB, II, JJ, REV )*/
	int mB=N+1;
	 NEG=-1;
	double * B = calloc(mB*mB,sizeof(double));

	pdlacp3_( &mB, &ZERO, VT, descVT, B, &mB, &NEG, &NEG, &ZERO);

	if(myrank_mpi == 0){	
		//for (i = 0; i < mB; i++) {
		//	printf("B[%d,%d]=%3.2f\n", i,i,B[i*mB+i]);
		//}
		int minVal  = VT[0];
		for (i = 0; i < mB; i++) {
			idx = i*mB + i;
			if( B[idx] < minVal ) {
				minVal = B[idx];
				count = i;	/* count is used here as the matching sample of the reference vector */
			}
		}

		/* print number of matching character */
		printf("P(%d): Character is number %d, value=%3.2f  !!!!!  ;-)\n", 
				myrank_mpi, count+1, VT[count*mV+count]);
	}
#ifdef DEBUG
	printf("P(%d) Waiting on barrier\n", myrank_mpi);
	Cblacs_barrier(context, "A");
#endif

	/* free allocated arrays */
	free(A); free(Cp); free(T); free(U); free(UT); free(S); free(V); free(VT); 
	free(X); free(D); free(x); free(y);

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
	Cblacs_barrier(context, "A");
	Cblacs_gridexit( 0 );
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

