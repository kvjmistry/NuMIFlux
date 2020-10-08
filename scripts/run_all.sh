# This bash script will execute all updates scripts and make all the plots for the technote
# To run do: source run_all.sh
source run_scripts.sh all fhc uboone
source run_scripts.sh all rhc uboone
source run_scripts.sh parent fhc
source run_scripts.sh parent rhc
source run_scripts.sh gsimple fhc
source run_scripts.sh gsimple rhc
source run_scripts.sh beamline fhc
source run_scripts.sh beamline rhc
source run_scripts.sh eventrate fhc
source run_scripts.sh eventrate rhc
