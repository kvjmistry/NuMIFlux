#!/usr/bin/env python
##################################################
# Makehist script to run on the grid
##################################################
import os, optparse, random, shutil, tarfile, sys
import subprocess, string

PWD = os.getenv("PWD")

##################################################
# Job Defaults -- if you dont set these, then it will use these values
##################################################
N_JOBS                = 1
RUN                   = 0
FILES_PER_JOB         = 1
MEMORY                = 700 # MB
LIFETIME              = 3   # Hrs
PROCESS_SHIFT         = 0

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

  # scratch area from which to send tarfile/config files to grid
  CACHEDIR = "/pnfs/uboone/scratch/users/{USER}/makehist/run{RUN}/CACHE/".format( USER = os.getenv("USER"), RUN = options.run)
  
  if os.path.isdir(CACHEDIR) == False:
      print CACHEDIR, " directory doen't exist, so creating...\n"
      os.mkdir(CACHEDIR)

  
  CACHE = CACHEDIR + str(random.randint(10000,99999)) + "/"
  os.mkdir(CACHE)

  COPYPATH = "{CACHE}/makehist_job.sh".format(CACHE = CACHE)
  shutil.copy('./makehist_job.sh', COPYPATH)

  print "\nOutput logfile(s):",logfile

  submit_command = ("jobsub_submit {GRID} {MEMORY_} -N {NJOBS} -dMAKEHIST {OUTDIR} "
      "-G uboone "
      "-e FILES_PER_JOB={FILES_PER_JOB} " 
      "-e RUN={RUN} "
      "-e PROCESS_SHIFT={PROCESS_SHIFT} "
      "-L {LOGFILE} "
      "file://{CACHE}makehist_job.sh".format(
      GRID       = ("--OS=SL6 -g "
                    "--resource-provides=usage_model=DEDICATED,OPPERTUNISTIC,OFFSITE "
                    "--role=Analysis "
                    "--expected-lifetime={LIFETIME}h ".format(LIFETIME = options.lifetime)),
      MEMORY_     =  "--memory {MEMORY}MB ".format(MEMORY = options.memory),
      NJOBS      = options.n_jobs,
      OUTDIR     = OUTDIR,
      FILES_PER_JOB = options.files_per_job,
      RUN        = options.run,
      PROCESS_SHIFT = options.process_shift,
      LOGFILE    = logfile,
      CACHE      = CACHE)
  )

  # Ship it on user input of yes
  print "\nSubmitting to grid:\n"+submit_command+"\n\n"
  submit = query_yes_no("Would you like to proceed to submit?")
  if submit == False:
      sys.exit()
  
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
  
  grid_group.add_option("--memory",
        default = MEMORY, type=int,
        help = "The memory for the job in MB. Default = %default.")

  grid_group.add_option("--lifetime",
        default = LIFETIME, type=int,
        help = "The expected lifetime for the job in hours. Default = %default.")
  
  grid_group.add_option("--process_shift",
        default = PROCESS_SHIFT, type=int,
        help = "Shift the process number, to re-run over a job that failed. Default = %default.")

  parser.add_option_group(grid_group)

  options, remainder = parser.parse_args()

  options = finalize_options(options)

  return options

def finalize_options(options):
  print 'run',                        options.run
  print 'n_jobs',                     options.n_jobs
  print 'files_per_job',              options.files_per_job
  print 'memory',                     options.memory
  print 'lifetime',                   options.lifetime
  print 'process_shift',              options.process_shift

  return options

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")

if __name__ == "__main__":
  sys.exit(main())
