#!/bin/bash
module load openmpi/1.4.3-intel
module list
export MKLPATH=/packages/intel/cmkl/10.0.2.018/lib/em64t
export IMKLPATH=/packages/intel/cmkl/10.0.2.018/include
echo $MKLPATH
echo $IMKLPATH

