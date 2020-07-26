#!/bin/bash

# this is an old file which I don't use anymore

source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup uboonecode v08_00_00_06 -q e17:prof


# Not explicitly using any ifdh commands in g4numi_job.sh right now, but
# leave these in, just in case.
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh
setup ifdhc #v2_2_3
export IFDH_GRIDFTP_EXTRA="-st 10" #set ifdh cp stall timeout to 10 sec
export IFDH_CP_MAXRETRIES=2


# We might actually need ifdhc v
