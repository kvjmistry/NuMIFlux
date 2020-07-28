#!/bin/bash

run(){
  COLOR="\033[1;33m"
  DEFAULT="\033[0m"
  echo -e "${COLOR}-> ${1}${DEFAULT}";
  eval ${1};
}

# New instructions

# jobsub_submit + client
run 'source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh'
run 'setup jobsub_client'
run 'setup ifdhc' #v2_2_3
run 'export IFDH_GRIDFTP_EXTRA="-st 10"'
run 'export IFDH_CP_MAXRETRIES=2'

echo -e "\033[1;31m Now re-setting some environmental variables \033[0m"
run 'EXPERIMENT=uboone'
run 'JOBSUB_GROUP=uboone'
run 'SAM_STATION=uboone'
run 'GROUP=uboone'
run 'SAM_EXPERIMENT=uboone'
run 'SAM_GROUP=uboone'

echo "remember to setup the proxy if you havent already done so"
echo "voms-proxy-destroy"
echo "kx509"
echo "voms-proxy-init -noregen -voms fermilab:/fermilab/uboone/Role=Analysis"
