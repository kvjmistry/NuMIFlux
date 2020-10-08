# To run this script do source run_scripts.sh <mode> <horn_mode>
# where <horn_mode> = fhc/rhc and <mode> = all/parent/gsimple/beamline/eventrate
mode=$1
horn=$2
detector=$3

if [ $mode == "all" ]
then
	root -l -q -b 'plot_uboone_flux.C("'$horn'","nue", "'$detector'")'
	root -l -q -b 'plot_uboone_flux.C("'$horn'","nuebar", "'$detector'")'
	root -l -q -b 'plot_uboone_flux.C("'$horn'","numu", "'$detector'")'
	root -l -q -b 'plot_uboone_flux.C("'$horn'","numubar", "'$detector'")'
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
	root -l -q -b 'plot_beamline_flux.C("nue", "'$horn'")'
	root -l -q -b 'plot_beamline_flux.C("nuebar", "'$horn'")'
	root -l -q -b 'plot_beamline_flux.C("numu", "'$horn'")'
	root -l -q -b 'plot_beamline_flux.C("numubar", "'$horn'")'
elif [ $mode == "eventrate" ]
then 
	root -l -q -b 'plot_event_rates.C("'$horn'")'
else 
	echo "Error in running the script, run by run_scripts.sh <mode> <horn_mode> where <horn_mode> = fhc/rhc and <mode> = all/parent/gsimple/beamline/eventrate"
fi
