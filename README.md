NuMI Flux Validation
====================
Scripts for NuMI flux testing and validation.

This is a set of scripts initially written bt Andy Mastbaum, I am just developing
with it. 

Flux Histograms
---------------
To plot flux histograms, run create an art ROOT file with PPFX weights:

    lar -c fcl/fluxreader_source_nova.fcl <input files>
    lar -c fcl/numippfx_test_one.fcl <output of previous step>

Then build this code and run:

    make
    bin/makehist <art ROOT files>

This produces a ROOT file `output.root` with histograms for the central value
and multisim universes.


## Validation Scripts:  
Most of the `plot_comparison` scripts have not been kept up to date at this present time. 
This means that they may not work out of the box and would need some modifications.

To produce the covariance matrix, correlation matrix and flux diagrams for microboone,
run this script. This script has mostly been kept up to date, but may need some minor changes
in the future.
`root -l 'plot_uboone_flux.C("output.root","numu")'`

I plan on moving the `Fluxystematics.C` and `FluxSyst_functions.h` files to another repository since these
are what were used for the nue flux sytematics and so dont really have a place anymore in this repo.

`checkfiles.C` was a script that was made to open root files and see if they are zombie. This was made due to the nova beamline files being broken and so I am not sure this is needed anymore.

`flugg_comparison.C`is a simple script that counts the number of kaons in a file. This was to try an explain the differences in decay at rest species in flugg and dk2nu. Since I believe there is actually a difference between the two, this script serves no purpose, but I will keep it here for now anyway.




Credits
-------
Geometry tools in src/geo are from LArLite: https://github.com/larlight/larlite.

