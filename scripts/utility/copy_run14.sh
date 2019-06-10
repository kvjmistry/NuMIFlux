# This script will copy files via ifdh to the current directory. It will first check if that file exists before copying. 
# It works by looping over the displayed entities from running the ls command. 

# -8 is the posion for the file e.g. g4numiv6_minervame_me000z200i_409_0006.root in the full string /pnfs/numix/flux/ME_FHC_G4NuMI_v6r1_0th/me000z200i/g4numi/ +2

for i in $(ifdh ls /pnfs/numix/flux/g4numi_syst_3rd_ana/test_g4numi_v6/me000z200i/minervame/run0014);do
	#echo "$i"
	if [ $i == "" -o  $i == "/pnfs/numix/flux/g4numi_syst_3rd_ana/test_g4numi_v6/me000z200i/minervame/run0014/" ]; then echo "skipping copy of $i";continue; fi
	#if [ $i == "/pnfs/nova/data/flux/g4numi/v6r1/me000z200i/g4numiv6_minervame_me000z200i_409_0006.root" ]; then echo "skipping copy of $i";continue; fi
	#if [ !  -e $(echo "$i" | cut -d'/' -f9) ]; then echo " copying file $(echo "$i" | cut -d'/' -f9)";ifdh cp $i ./; fi
	if [ !  -e $(echo "$i" | cut -d'/' -f10) ]; then echo " copying file $(echo "$i" | cut -d'/' -f10)";ifdh cp $i ./; fi
done
