#!/bin/sh
 
#echo "Hello, world!"



counter=1

while [ $counter -le 100 ]
do
	  cd OM/gtg 
	./lagrangian_OM_gtg
	
	cd ../../mov_est/gtg/
	./lagrangian_est_gtg 1>"$counter-est.log"  2>"$counter-est-error.log"
	
	cd ../../../simeval/SimResult_gtg_5areas_tau04_delta2
	
	#source(simRun.r, chdir = TRUE) | R --vanilla --slave
	
	Rscript simRun.r

	
	cd ../../admb
	
	counter=$(( $counter + 1 ))

 done



