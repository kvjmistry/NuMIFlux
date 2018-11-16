NuMI Flux Validation
====================
Scripts for NuMI flux testing and validation.

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

Credits
-------
Geometry tools in src/geo are from LArLite: https://github.com/larlight/larlite.

