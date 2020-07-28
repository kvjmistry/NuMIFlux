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
echo "IN PROCESS:    ${PROCESS}"
echo "SHIFT PROCESS: ${PROCESS_SHIFT}"
echo "HPSET:         ${HPSET}"
echo "Input list:    ${FLIST}" 
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
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup
source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setup 
export PRODUCTS=${PRODUCTS}:/cvmfs/fermilab.opensciencegrid.org/products/common/db/
echo "setup uboonecode v08_08_00_19 -q e17:prof"
setup uboonecode v08_00_00_19 -q e17:prof
echo "setup ifdhc"
setup ifdhc #v2_2_3
export IFDH_GRIDFTP_EXTRA="-st 10" #set ifdh cp stall timeout to 10 sec
export IFDH_CP_MAXRETRIES=5

# For some reason these ups products dont get auto setup
echo "setup sqlite v3_20_01_00"
echo "gcc v7_3_0"
setup sqlite v3_20_01_00
setup gcc v7_3_0

echo
echo "======== COPY THE EXECUTABLES ========"
if [ ${HPSET} -eq 0 ]; then 
	# copy over the default file and also the filelist
	echo "ifdh cp /pnfs/uboone/resilient/users/kmistry/exe/makehist_v2 makehist"
	ifdh cp /pnfs/uboone/resilient/users/kmistry/exe/makehist_v2 makehist
	echo "ifdh cp ${FLIST} files.list"
	ifdh cp ${FLIST} files.list
elif [ ${HPSET} -eq 1 ]; then 
	# copy over the default file and also the filelist
	echo "ifdh cp /pnfs/uboone/resilient/users/kmistry/exe/makehist_v2 makehist"
	ifdh cp /pnfs/uboone/resilient/users/kmistry/exe/makehist_v2 makehist
	echo "ifdh cp ${FLIST} files.list"
	ifdh cp ${FLIST} files.list
elif [ ${HPSET} -eq 2 ]; then
	# copy over the set 2 file and also the filelist
	echo "ifdh cp /pnfs/uboone/resilient/users/kmistry/exe/makehist_v2_set2 makehist"
	ifdh cp /pnfs/uboone/resilient/users/kmistry/exe/makehist_v2_set2 makehist
	echo "ifdh cp ${FLIST} files.list"
	ifdh cp ${FLIST} files.list
elif [ ${HPSET} -eq 3 ]; then
	# copy over the set 3 file and also the filelist
	echo "ifdh cp /pnfs/uboone/resilient/users/kmistry/exe/makehist_v2_set3 makehist"
	ifdh cp /pnfs/uboone/resilient/users/kmistry/exe/makehist_v2_set3 makehist
	echo "ifdh cp ${FLIST} files.list"
	ifdh cp ${FLIST} files.list
fi

chmod u+x makehist

echo
echo "======== Split the FILELIST ========"

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
if [ ${HPSET} -eq 0 ]; then 
	echo "./makehist g4*.root"
	./makehist g4*.root
else
	# This is to include the multisims
	echo "./makehist g4*.root CV"
	./makehist g4*.root CV
fi  

# Rename the output file
echo
echo "======== Renaming the output file ========"
echo "mv output.root output_uboone_${PROCESS}_Set${HPSET}.root"
mv output.root output_uboone_${PROCESS}_Set${HPSET}.root

echo
echo "Moving output to CONDOR_DIR_MAKEHIST"
echo "ifdh cp -D output_*.root $CONDOR_DIR_MAKEHIST"
ifdh cp -D output_*.root $CONDOR_DIR_MAKEHIST

# Now exit gracefully
echo END:  `date`
exit 0
