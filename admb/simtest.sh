#!/bin/sh
 
#echo "Hello, world!"



counter=1

while [ $counter -le 100 ]
do
	  cd OM/simple 
	./lagrangian_OM
	
	cd ../../mov_est/simple/
	./lagrangian_estlog
	
	cd ../../../simeval/SimResult_3areas_tau04/
	
	#source(simRun.r, chdir = TRUE) | R --vanilla --slave
	
	Rscript simRun.r

	
	cd ../../admb
	
	counter=$(( $counter + 1 ))

	

 done



