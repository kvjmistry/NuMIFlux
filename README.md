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


Latest Updates:
Can now run new plot comparison all function with arguments:
// mippon/mipoff, input file name, name of product of weights, mode 
	root -l '/uboone/app/users/kmistry/PPFX/numi-validation/scripts/plot_comparison_all.C("mippoff","output_all.root","noThinKaon", "numu")'


Credits
-------
Geometry tools in src/geo are from LArLite: https://github.com/larlight/larlite.

