# This script will copy files via ifdh to the current directory. It will first check if that file exists before copying. 
# It works by looping over the displayed entities from running the ls command. 

# -8 is the posion for the file e.g. g4numiv6_minervame_me000z200i_409_0006.root in the full string /pnfs/numix/flux/ME_FHC_G4NuMI_v6r1_0th/me000z200i/g4numi/ +2

for i in $(ifdh ls /pnfs/minerva/mc_generation_flux/mc-flux/mc/g4numiv6/00/00/00/06);do
	#echo "$i"
	#echo "$i" | cut -c66-100;
	#echo "g4numiv6_dk2nu_minerva13_le010z185i";
	if [ $i == "" -o  $i == "/pnfs/minerva/mc_generation_flux/mc-flux/mc/g4numiv6/00/00/00/06/" ]; then echo "skipping copy of $i";continue; fi
	if  [ "$(echo "$i" | cut -c66-100)" != "g4numiv6_dk2nu_minerva13_le010z185i" ] && [ "$(echo "$i" | cut -c66-99)" != "g4numiv6_dk2nu_minerva1_le010z185i" ] && [ "$(echo "$i" | cut -c66-99)" != "g4numiv6_dk2nu_minerva3_le010z185i" ]; then echo "skipping copy of $i";continue; fi
	#if [ !  -e $(echo "$i" | cut -d'/' -f12) ]; then echo " copying file $(echo "$i" | cut -d'/' -f12)"; fi
	if [ !  -e $(echo "$i" | cut -d'/' -f12) ]; then echo " copying file $(echo "$i" | cut -d'/' -f12)";ifdh cp $i ./; fi
done
