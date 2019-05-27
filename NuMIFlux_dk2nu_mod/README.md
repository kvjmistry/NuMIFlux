# NuMIFlux
Code to calculate the NuMI Flux at MicroBooNE.

Wiki page: https://cdcvs.fnal.gov/redmine/projects/ubooneoffline/wiki/NuMI_Flux_Histograms.


## Updates  
Krishan: Updated scripts to use the dk2nu file format as an input as opposed to flugg.


## First Usage

From the MicroBooNE gpvm machines:

```
source SetupNuMIFlux.sh

root -l compileNuMIFlux.C

python RunNuMIFlux.py

```

## General Instructions

Setup the environment with:

```
source SetupNuMIFlux.sh

```

Then to run   
```
python RunNuMIFlux.py

```

If you make any edits to the NuMIFlux.(cc/hh) scripts re-compile with:

```
root -l compileNuMIFlux.C

```

