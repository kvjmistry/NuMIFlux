# NuMIFlux
Code to calculate the NuMI Flux at MicroBooNE.

Wiki page: https://cdcvs.fnal.gov/redmine/projects/ubooneoffline/wiki/NuMI_Flux_Histograms.


##Updates  
Krishan: Updated scripts to use the dk2nu file format as an input as opposed to flugg.


## First Usage

From the MicroBooNE gpvm machines:

```
git clone git@github.com:marcodeltutto/NuMIFlux.git
cd NuMIFlux/
source SetupNuMIFlux.sh

root -l
> .L NuMIFlux.cc+
> .q

```

## General Instructions

Run with:

```
source SetupNuMIFlux.sh

root -l
> .L NuMIFlux.cc+
> .q
```

or use the python script:

```
source SetupNuMIFlux.sh
python RunNuMIFlux.py
```

Compile with:

```
root -l
> .L NuMIFlux.cc+
```


Generate new FluxNtuple_C.so (couldn't get to working with dk2nu so removed dependancy):

```
cd FluggNtuple
root -l
.L FluxNtuple.C+
```
