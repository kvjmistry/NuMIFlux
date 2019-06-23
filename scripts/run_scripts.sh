# To run this script do source run_scripts.sh <mode> <horn_mode>
# where <horn_mode> = fhc/rhc and <mode> = all/parent/gsimple/beamline
mode=$1
horn=$2

if [ $mode == "all" ]
then
	root -l -q -b 'plot_uboone_flux.C("'$horn'","nue")'
	root -l -q -b 'plot_uboone_flux.C("'$horn'","nuebar")'
	root -l -q -b 'plot_uboone_flux.C("'$horn'","numu")'
	root -l -q -b 'plot_uboone_flux.C("'$horn'","numubar")'
elif [ $mode == "parent" ]
then 
	root -l -q -b 'plot_parent_flux.C("'$horn'","nue")'
	root -l -q -b 'plot_parent_flux.C("'$horn'","nuebar")'
	root -l -q -b 'plot_parent_flux.C("'$horn'","numu")'
	root -l -q -b 'plot_parent_flux.C("'$horn'","numubar")'
elif [ $mode == "gsimple" ]
then 
	root -l -q -b 'plot_gsimple_flux.C("'$horn'","nue")'
	root -l -q -b 'plot_gsimple_flux.C("'$horn'","nuebar")'
	root -l -q -b 'plot_gsimple_flux.C("'$horn'","numu")'
	root -l -q -b 'plot_gsimple_flux.C("'$horn'","numubar")'
elif [ $mode == "beamline" ]
then 
	root -l -q -b 'plot_beamline_flux.C("nue")'
	root -l -q -b 'plot_beamline_flux.C("nuebar")'
	root -l -q -b 'plot_beamline_flux.C("numu")'
	root -l -q -b 'plot_beamline_flux.C("numubar")'
else 
	echo "Error in running the script, run by un_scripts.sh <mode> <horn_mode> where <horn_mode> = fhc/rhc and <mode> = all/parent/gsimple/beamline"
fi
