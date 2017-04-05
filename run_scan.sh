#!/bin/bash



#for f in `ls ../../data/garfield_data/*.root` 
#do 
#	echo "File: $f"; 
#	python extract_waveform.py $f
#done
#exit 0






for f in `ls ../../data/garfield_data/*_Y.txt` 
do 
	echo "File: $f"; 

	#python diffusion_effect.py $f
	#continue

	for t in 1.0
	do 
		echo "tau: $t";
		for d in 100.0
		do 
			echo "distance: $d";
			for freq in 0.5 1.0 2.0
			do 
				echo "freq: $freq";			
				python waveform_ana.py --tau $t --diffusion --filter gaus $f --distance $d --freq $freq --batch; 
			done			
		done
	done
done
