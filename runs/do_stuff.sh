#!/bin/bash

#while [ $i -lt 18 ] ; do
for i in {0..17}; do
#for i in {0..05}; do
#for i in {6..11}; do
#for i in {12..17}; do
	let angle=5*i;
#echo angle is "$angle"
#echo '%04d' "$angle"
	b="$( printf '%02d' "$angle" )"
	echo "$b"
	#mkdir b"$b"/
	cd b"$b"/
	#cp ../ss.par .
	#python3.6 /Users/hagrusa/planetary/scripts/init.py ../ceres ../imp "$angle"
	
	#cp init.ss ss.0000000000
	
#	let angle=5*i;
#	b="$( printf '%02d' "$angle" )"
#	echo "$b"
#	cd b"$b"/
#	../pkdgrav_pthread_DEM_linux -sz 4 ss.par > output 2>&1 &
#	cd ../

	#do analysis:
	python /home/hagrusa/planetary/scripts/mass_v_time.py ./ mass_time_b"$b".txt
	cp mass_time_b"$b".txt /home/hagrusa/planetary/runs/results/
	cd ../

done


