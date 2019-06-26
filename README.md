NuMI Flux Validation
====================
Scripts for NuMI flux and validation.

This is a set of scripts initially written by Andy Mastbaum, Joseph Zennamo and Zarko Pavolvic that has built upon. 

Flux Histograms
---------------
To plot flux histograms, run create an art ROOT file with PPFX weights:

    lar -c fcl/FluxReader/fluxreader_uboone.fcl <input files>
    lar -c fcl/FluxReader/ppfx_uboone.fcl <output of previous step>

Then build this code, in the main directory, run:

    make
    
This will compile the src/makehist.cpp into an executable file which can then be run with the command:
    
    bin/makehist <art ROOT file> <detector>

where `<detector>` can be either `uboone` or `nova`.

If you want to run over many files then specify the path to the files and use a `*`
wildcard. 

After running, a ROOT file `output.root` with histograms for the central value
and multisim universes is produced.


## Validation Scripts:  
This development area was initially designed for ppfx validation with the nova config. 
The nova config has not been kept up-to-date with the histogram format in `output.root`.
To get the scripts in the scripts/nova directory to work, path and file names will need
to be updated.

A bash script has been setup `scripts/run_all.sh` to run all the plotting scripts in one go. 
These are run in root batch mode and so will be faster than just running each individual script
which takes a lot of time opening the root canvases. All the resulting plots will be creating into
a plots folder that will be created. 

To run the script, simply do `source run_all.sh` when in the scripts directory.

The input file paths have been hardcoded. If a user wants to update these paths, they will, as it stands,
need to edit the paths in each of the scripts, shouldnt be too much effort I hope! 



Credits
-------
Geometry tools in src/geo are from LArLite: https://github.com/larlight/larlite.

