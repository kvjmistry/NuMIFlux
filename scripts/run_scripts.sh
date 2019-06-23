# To run this script do source run_scripts.sh <mode> <horn_mode>
# where <horn_mode> = fhc/rhc and <mode> = all/parent/gsimple
mode=$1
horn=$2

if [ $mode == "all" ]
then
	root -l -q 'plot_uboone_flux.C("'$horn'","nue")'
	root -l -q 'plot_uboone_flux.C("'$horn'","nuebar")'
	root -l -q 'plot_uboone_flux.C("'$horn'","numu")'
	root -l -q 'plot_uboone_flux.C("'$horn'","numubar")'
elif [ $mode == "parent" ]
then 
	root -l -q 'plot_parent_flux.C("'$horn'","nue")'
	root -l -q 'plot_parent_flux.C("'$horn'","nuebar")'
	root -l -q 'plot_parent_flux.C("'$horn'","numu")'
	root -l -q 'plot_parent_flux.C("'$horn'","numubar")'
elif [ $mode == "gsimple" ]
then 
	root -l -q 'plot_gsimple_flux.C("'$horn'","nue")'
	root -l -q 'plot_gsimple_flux.C("'$horn'","nuebar")'
	root -l -q 'plot_gsimple_flux.C("'$horn'","numu")'
	root -l -q 'plot_gsimple_flux.C("'$horn'","numubar")'
else 
	echo "Enter Mode from: all, parent, gsimple"
fi
