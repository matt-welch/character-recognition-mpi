# Makefile for pdgemv example PBlacs program
# assumes existence of $MKLPATH environment variable??
#
pdgemv:  pdgemv.c prepare_env.sh
	. prepare_env.sh
	mpicc -o pdgemv pdgemv.c -I$IMKLPATH -L$MKLPATH /packages/intel/cmkl/10.0.2.018/lib/em64t/libmkl_scalapack.a -lmkl_blacs_openmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -openmp
