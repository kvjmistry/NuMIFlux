#!/usr/bin/env python

import argparse
import os.path
import sys
import tarfile
import re
import subprocess

parser = argparse.ArgumentParser(description='Submit beam data jobs.')

parser.add_argument("-f", "--first-run", dest="firstrun",
                    required=True,
                    help="First run.")
parser.add_argument("-l", "--last-run", dest="lastrun",
                    required=True,
                    help="Last run.")
parser.add_argument("-i", "--input-filelist", dest="inputlist",
                    default="bnbcos.list",
                    help="Path to input file list. (default=bnbcos.list)")
parser.add_argument("-o", "--output-path", dest="outputpath",
                    default="/pnfs/uboone/scratch/users/%s/bnb_redecay_to_gsimple"%os.environ['USER'],
                    help="Path where to copy final output. (default=/pnfs/uboone/scratch/users/%s/bnb_redecay_to_gsimple)"%os.environ['USER'])
parser.add_argument("-d", "--debug",action='store_true',
                    help="Will not delete submission files in the end. Useful for debugging and will only print the submission command on screen.")

args = parser.parse_args()

#now create jobfiles_*.tar that is shipped with the job
#this includes the executable
#tar -cf jobfiles.tar --transform '!^[^/]*/!!' file1 file2 file3
tarfilename="jobfiles_%i.tar.bz2"%os.getpid()
outtar = tarfile.open(tarfilename, mode='w:bz2')
outtar.add("make_tree",arcname="make_tree")
outtar.add(args.inputlist,arcname=args.inputlist)
outtar.close()

ofstr='''
#!/bin/bash

export _RUN_NUMBER=$((PROCESS+%(firstrun)s))

echo "Running $0 on "$HOSTNAME >>${_RUN_NUMBER}.out 2>&1
echo "Cluster: " ${CLUSTER} >>${_RUN_NUMBER}.out 2>&1
echo "Process: " ${PROCESS} >>${_RUN_NUMBER}.out 2>&1

echo " Sourcing everything...." >>${_RUN_NUMBER}.out 2>&1

source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setup >>${_RUN_NUMBER}.out 2>&1
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh >>${_RUN_NUMBER}.out 2>&1
source /cvmfs/fermilab.opensciencegrid.org/products/larsoft/setup >>${_RUN_NUMBER}.out 2>&1

setup larsoftobj v1_13_00 -q e10:prof  >>${_RUN_NUMBER}.out 2>&1
setup uboonecode v06_26_01_07 -q e10:prof  >>${_RUN_NUMBER}.out 2>&1

mv ${_RUN_NUMBER}.out ${_CONDOR_SCRATCH_DIR}/.
cd ${_CONDOR_SCRATCH_DIR}

cp $INPUT_TAR_FILE . 
tar -jvxf `basename ${INPUT_TAR_FILE}` >>${_RUN_NUMBER}.out 2>&1

fullpath=`head -n ${_RUN_NUMBER} %(inputfile)s | tail -n 1`
echo $fullpath >>${_RUN_NUMBER}.out 2>&1

echo " Copying File...." >>${_RUN_NUMBER}.out 2>&1
ifdh cp $fullpath ${PWD}/test.root >>${_RUN_NUMBER}.out 2>&1

echo "What's in here? "
ls >>${_RUN_NUMBER}.out 2>&1

echo " Run time! " >>${_RUN_NUMBER}.out 2>&1

./make_tree test.root  >>${_RUN_NUMBER}.out 2>&1

mkdir ${_RUN_NUMBER}

cp ${_RUN_NUMBER}.out ${_RUN_NUMBER}/.
cp out*.root ${_RUN_NUMBER}/.

ifdh mkdir %(outputdir)s/
ifdh mkdir %(outputdir)s/${_RUN_NUMBER}
ifdh cp -r ${_RUN_NUMBER} %(outputdir)s/${_RUN_NUMBER}/

'''%{'inputfile':args.inputlist,'firstrun':args.firstrun,'outputdir':args.outputpath}

runjobfname="runjob_%i.sh"%os.getpid()
of=open(runjobfname,'w')
of.write(ofstr)
of.close()

cmd="jobsub_submit --memory=2000MB --group=uboone -N %i --tar_file_name=dropbox://%s file://%s"%(int(args.lastrun)-int(args.firstrun)+1,os.path.abspath(tarfilename),os.path.abspath(runjobfname))

if (not args.debug):
    print "Running submit cmd:"
    print cmd
    os.system(cmd)
else:
    print "Would have ran:"
    print cmd

#Delete temp files unless debugging
if (not args.debug):
    os.remove(tarfilename)
    os.remove(runjobfname)



