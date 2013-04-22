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
	MPI_Init( &argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);
	/************  BLACS ***************************/
	int ictxt, nprow, npcol, myrow, mycol,nb;
	int info,itemp;
	int ZERO=0,ONE=1; /* constants for passing into pblacs routines */
	char* testFilename = "testcharacter.bin";
	char* refFilename = "referenceset.bin";

	/* read in data with MPI-IO */

	nprow = 2; npcol = 2; nb =2;
	Cblacs_pinfo( &myrank_mpi, &nprocs_mpi ) ;
	Cblacs_get( -1, 0, &ictxt );
	Cblacs_gridinit( &ictxt, "Row", nprow, npcol );
	Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );

	/* determine the size of the matrix */
	/* matrix size is already known at 72*(128*128)= 72*16384 */
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
	descinit_(descA, &M,   &M,   &nb,  &nb,  &ZERO, &ZERO, &ictxt, &mA,  &info);
	descinit_(descx, &M, &ONE,   &nb, &ONE,  &ZERO, &ZERO, &ictxt, &nx, &info);

	descinit_(descy, &M, &ONE,   &nb, &ONE,  &ZERO, &ZERO, &ictxt, &my, &info);
	double *x = (double*) malloc(nx*sizeof(double));
	double *y = (double*) calloc(my,sizeof(double));
	double *A = (double*) malloc(mA*nA*sizeof(double));
	int sat,sut;
	for(i=0;i<mA;i++)

		for(j=0;j<nA;j++){
			sat= (myrow*nb)+i+(i/nb)*nb;
			sut= (mycol*nb)+j+(j/nb)*nb;
			A[j*mA+i]=AA(sat,sut);
		}


	for(i=0;i<nx;i++){
		sut= (myrow*nb)+i+(i/nb)*nb;
		x[i]=X[sut];
	}

	double alpha = 1.0; double beta = 0.0;
	pdgemv_("N",&M,&M,&alpha,A,&ONE,&ONE,descA,x,&ONE,&ONE,descx,&ONE,&beta,y,&ONE,&ONE,descy,&ONE);

	Cblacs_barrier(ictxt,"A");
	for(i=0;i<my;i++)
		printf("rank=%d %.2f \n", myrank_mpi,y[i]);
	Cblacs_gridexit( 0 );
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
 * (3) Projections: Multiply X=U'A  and  x=UU'T 
 * (4) Find Minimum: Subtract x from each column of X to make D,
 *		the "difference matrix"*
 * (5) Calculate DD', find the minimum diagonal element. The column of this
 *		element is the matching character
 */

/* QUESTIONS:  
 * (1) does this program need to run with various numbe of processors?
 *		1,4,16,64?
 *	*/
