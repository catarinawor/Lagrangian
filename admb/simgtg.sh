#!/bin/sh
 
#echo "Hello, world!"



counter=1

while [ $counter -le 200 ]
do
	  cd OM/gtg 
	./lagrangian_OM_gtg
	
	cd ../../mov_est/gtg/
	./lagrangian_est_gtg
	
	cd ../../../simeval/SimResult_gtg_3areas_tau1_delta2
	
	#source(simRun.r, chdir = TRUE) | R --vanilla --slave
	
	Rscript simRun.r

	
	cd ../../admb
	
	counter=$(( $counter + 1 ))

 done



