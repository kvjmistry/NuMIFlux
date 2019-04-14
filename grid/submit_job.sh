# Notes, it took about 15 mins to copy over 10 files
# Files sizes for 500k POT are around 1GB per file, so refrain from submitting > 25 files_per job otherwise more disk space will need to be requested.
# This can be avoided if you manage to get xrootd streaming to work! I just cant seem to get it working for some unknown reason.

python ProcessMakehist.py --n_jobs=1 --files_per_job=10 --run=0
