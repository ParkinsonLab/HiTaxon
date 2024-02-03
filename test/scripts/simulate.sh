#!/usr/bin/env bash
TESTDIR=$1 

module load intel/2019u4 gcc/8.3.0 r/4.1.2
#${SCINET_R_ROOT}/bin/Rscript --no-restore "$@"

#Use Polyester to simulate test environments
Rscript simulate.R $TESTDIR
