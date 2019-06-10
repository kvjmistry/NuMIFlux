# This script will copy files via ifdh to the current directory. It will first check if that file exists before copying. 
# It works by looping over the displayed entities from running the ls command. 

for i in $(ifdh ls /pnfs/nova/data/flux/g4numi/v6r1/me000z200i);do
	#echo "$i"
	if [ $i == "/pnfs/nova/data/flux/g4numi/v6r1/me000z200i/log/" -o $i == "/pnfs/nova/data/flux/g4numi/v6r1/me000z200i/" ]; then echo "skipping copy of $i";continue; fi
	#if [ $i == "/pnfs/nova/data/flux/g4numi/v6r1/me000z200i/g4numiv6_minervame_me000z200i_409_0006.root" ]; then echo "skipping copy of $i";continue; fi
	#if [ $i == "/pnfs/nova/data/flux/g4numi/v6r1/me000z200i/g4numiv6_minervame_me000z200i_414_0006.root" ]; then echo "skipping copy of $i";continue; fi
	#if [ $i == "/pnfs/nova/data/flux/g4numi/v6r1/me000z200i/g4numiv6_minervame_me000z200i_459_0006.root" ]; then echo "skipping copy of $i";continue; fi
	#if [ $i == "/pnfs/nova/data/flux/g4numi/v6r1/me000z200i/g4numiv6_minervame_me000z200i_473_0006.root" ]; then echo "skipping copy of $i";continue; fi
	#if [ $i == "/pnfs/nova/data/flux/g4numi/v6r1/me000z200i/g4numiv6_minervame_me000z200i_70_0006.root" ]; then echo "skipping copy of $i";continue; fi
	if [ !  -e $(echo "$i" | cut -d'/' -f9) ]; then echo " copying file $(echo "$i" | cut -d'/' -f9)";ifdh cp $i ./; fi
done
