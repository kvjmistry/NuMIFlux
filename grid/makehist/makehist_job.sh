#!/bin/bash

# Cluster details
echo Start  `date`
echo Site:${GLIDEIN_ResourceName}
echo "the worker node is " `hostname` "OS: "  `uname -a`
whoami
id

# Spit out the environmental variables sent via -e option in jobsub
# Process number is jobnumber in jobid url
echo
echo "Run: ${RUN}"
echo "Files per job: ${FILES_PER_JOB}"
echo "IN PROCESS: ${PROCESS}"
echo "SHIFT PROCESS: ${PROCESS_SHIFT}"
echo
PROCESS=$((PROCESS + PROCESS_SHIFT)) # shift process numbers
echo "Process ID is ${PROCESS}"

echo
echo "======== cd to CONDOR_DIR_INPUT ========"
cd $CONDOR_DIR_INPUT

echo
echo "======== ls ========"
ls

echo
echo "======== SETUP THE ENV ========"
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh
echo "setup uboonecode v07_08_00 -q e17:prof"
setup uboonecode v08_00_00_14 -q e17:prof
setup ifdhc #v2_2_3
export IFDH_GRIDFTP_EXTRA="-st 10" #set ifdh cp stall timeout to 10 sec
export IFDH_CP_MAXRETRIES=5

# Try exporting all the env variables in interactive library path as not all of them seem to get set...
ifdh cp -D /pnfs/uboone/resilient/users/kmistry/files/Environmental_variables.txt .
for e in $(cat Environmental_variables.txt)
do
	export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$e
done

# copy a tarball and untar it if you need it
#ifdh cp /pnfs/uboone/resilient/users/kmistry/tars/PPFX_uboone_bugfix_slim_withtilt.tar intar.tar
#tar -vxf intar.tar

echo
echo "======== COPY THE EXECUTABLES ========"
echo "ifdh cp /pnfs/uboone/resilient/users/kmistry/exe/makehist makehist"
ifdh cp /pnfs/uboone/resilient/users/kmistry/exe/makehist makehist
chmod u+x makehist

echo
echo "======== COPY THE FILELIST ========"
ifdh cp /pnfs/uboone/persistent/users/kmistry/PPFX/makehist/files_run${RUN}.list files.list

# Now get cut the file list up into the relavent chunk
LINE1=$((PROCESS*${FILES_PER_JOB} + 1))
LINE2=$((PROCESS*${FILES_PER_JOB} + (${FILES_PER_JOB})))

echo
echo "LINE1: ${LINE1}, LINE2: ${LINE2}"

# Create the subfiles list
sed -n "${LINE1},${LINE2}p" < files.list | tee subfiles.list

echo
echo "======== Now copying over the files: ========"
# Now copy the files over to the local directory
for f in $(cat subfiles.list)
do
	echo $f
	ifdh cp -D $f .
done;

echo
echo "======== DOING AN LS  ========"
ls -ltrh

echo
echo "======== EXECUTING 2D makehist ========"
echo "./makehist g4*.root"
./makehist uboone g4*.root

# Rename the output file
echo
echo "======== Renaming the output file ========"
echo "mv output.root output_uboone_${PROCESS}.root"
mv output.root output_uboone_${PROCESS}.root

echo
echo "Moving output to CONDOR_DIR_MAKEHIST"
echo "ifdh cp -D output_*.root $CONDOR_DIR_MAKEHIST"
ifdh cp -D output_*.root $CONDOR_DIR_MAKEHIST

# Now exit gracefully
echo END:  `date`
exit 0
