#!/usr/bin/env python
##################################################
# Makehist script to run on the grid
##################################################
import os, optparse, random, shutil, tarfile, sys
import subprocess, string

PWD = os.getenv("PWD")

##################################################
# Job Defaults
##################################################
N_JOBS                = 1
RUN            	      = 0
FILES_PER_JOB         = 1

def main():
  options = get_options()

  # Create a beam config directory
  MAKEHISTDIR = "/pnfs/uboone/scratch/users/{USER}/makehist/".format( USER = os.getenv("USER"))

  if os.path.isdir(MAKEHISTDIR) == False:
      print MAKEHISTDIR, " directory doen't exist, so creating...\n"
      os.mkdir(MAKEHISTDIR)

  # Create a run number directory
  RUNDIR = "/pnfs/uboone/scratch/users/{USER}/makehist/run{RUN}/".format( USER = os.getenv("USER"), RUN = options.run)

  if os.path.isdir(RUNDIR) == False:
      print RUNDIR, " directory doen't exist, so creating...\n"
      os.mkdir(RUNDIR)

  # Create a output file directory  
  OUTDIR = "/pnfs/uboone/scratch/users/{USER}/makehist/run{RUN}/files/".format( USER = os.getenv("USER"), RUN = options.run)

  if os.path.isdir(OUTDIR) == False:
      print OUTDIR, " directory doen't exist, so creating...\n"
      os.mkdir(OUTDIR)

  # Create a log file directory  
  LOGDIR = "/pnfs/uboone/scratch/users/{USER}/makehist/run{RUN}/log/".format( USER = os.getenv("USER"), RUN = options.run)

  if os.path.isdir(LOGDIR) == False:
      print LOGDIR, " directory doen't exist, so creating...\n"
      os.mkdir(LOGDIR)

  logfile = LOGDIR + "/makehist_{RUN}_\$PROCESS.log".format(RUN = options.run)

  COPYPATH = "{LOGDIR}/makehist_job.sh".format(LOGDIR = LOGDIR)
  shutil.copy('./makehist_job.sh', COPYPATH)

  print "\nOutput logfile(s):",logfile

  submit_command = ("jobsub_submit {GRID} {MEMORY} -N {NJOBS} -dMAKEHIST {OUTDIR} "
      "-G uboone "
      "-e FILES_PER_JOB={FILES_PER_JOB} " 
      "-e RUN={RUN} "
      "-L {LOGFILE} "
      "file://{LOGDIR}makehist_job.sh".format(
      GRID       = ("--OS=SL6 -g "
                    "--resource-provides=usage_model=DEDICATED,OPPORTUNISTIC,OFFSITE "
                    "--role=Analysis "
		    "--expected-lifetime=1h "),
      MEMORY     =  "--memory 800MB ",
      NJOBS      = options.n_jobs,
      OUTDIR     = OUTDIR,
      FILES_PER_JOB = options.files_per_job,
      RUN        = options.run,
      LOGFILE    = logfile,
      LOGDIR     = LOGDIR)
  )

  #Ship it
  print "\nSubmitting to grid:\n"+submit_command+"\n"
  status = subprocess.call(submit_command, shell=True)

def get_options():
  parser       = optparse.OptionParser(usage="usage: %prog [options]")
  grid_group   = optparse.OptionGroup(parser, "Grid Options")

  grid_group.add_option("--n_jobs",
        default = N_JOBS, type=int,
        help = "Number of jobs. Default = %default.")

  grid_group.add_option("--run",
        default = RUN, type=int,
        help = "Tag on the end of outfiles. Doubles as random # seed. Default = %default.")

  grid_group.add_option("--files_per_job",
        default = FILES_PER_JOB, type=int,
        help = "The number of files to process in one job. Default = %default.")

  parser.add_option_group(grid_group)

  options, remainder = parser.parse_args()

  options = finalize_options(options)

  return options

def finalize_options(options):
  print 'run',                        options.run
  print 'n_jobs',                     options.n_jobs
  print 'files_per_job',              options.files_per_job

  return options


if __name__ == "__main__":
  sys.exit(main())
