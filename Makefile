#/*******************************************************************************
# * FILENAME:    Makefile
# * DESCRIPTION: Makefile for pdgemv example PBlacs program 
# * AUTHOR:      James Matthew Welch [JMW]
# * SCHOOL:      Arizona State University
# * CLASS:       CSE598: High Performance Computing
# * INSTRUCTOR:  Dr. Gil Speyer
# * SECTION:     20520
# * TERM:        Spring 2013
# *******************************************************************************/
#
# prepare_env.sh loads the openmp module and creates the environment variables $IMKLPATH and $MKLPATH

recognizer: recognizer.c
	mpicc -o recognizer recognizer.c -I$(IMKLPATH) -L$(MKLPATH) -lmkl_scalapack -lmkl_blacs_openmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -openmp -limf

all: recognizer demo tagit

debug: recognizer.c
	mpicc -o recognizer recognizer.c -DDEBUG -I$(IMKLPATH) -L$(MKLPATH) -lmkl_scalapack -lmkl_blacs_openmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -openmp -limf

verbose: recognizer.c tagit
	mpicc -o recognizer recognizer.c -DDEBUG -DVERBOSE  -I$(IMKLPATH) -L$(MKLPATH) -lmkl_scalapack -lmkl_blacs_openmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -openmp -limf

demo:  pdgemv.c prepare_env.sh
	mpicc -o pdgemv pdgemv.c -I$(IMKLPATH) -L$(MKLPATH) -lmkl_scalapack -lmkl_blacs_openmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -openmp -limf

clean: 
	rm -rf pdgemv *.o recognizer

tidy: tagit clean
	rm -rf *~ .*~ *.swp .*.swp

rundemo: pdgemv
	mpiexec -np 4 ./pdgemv

run: recognizer
	mpiexec -np 4 ./recognizer

tagit:  
	ctags --c-kinds=+defglmpstux -R *.c *.h
# location of the mkl_scalapack library: 
#/packages/intel/cmkl/10.0.2.018/lib/em64t/libmkl_scalapack.a 
