#!/bin/bash

for f in `ls ../../data/garfield_data/*.txt` 
do 
	echo "File: $f"; 
	for t in 0.5 1 2 5 10
	do 
		echo "tau: $t";
		for freq in 0.5 1.0 2.0 5.0 10.0
		do 
			echo "freq: $t";			
			python waveform_ana.py -t $t $f --freq $freq --batch; 
		done

	done
done
