# bash script that takes a directory and puts the full path of the files in that directory
# into a file list so we can submit grid jobs
# USAGE:
# source make_list.sh <input path> <beamline run number>
# e.g. source make_list.sh /pnfs/uboone/scratch/users/kmistry/g4numi/me000z200i/run1/files/ 1

for file in $(ls $1 ); do
  echo "${1}/${file}" >> filelist_fullpath_run${2}.list
done
