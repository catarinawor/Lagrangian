#!/bin/sh
 
#echo "Hello, world!"



counter=1
while [ $counter -le 3 ]
do
	  cd OM/gtg 
	./lagrangian_OM_gtg
	
	cd ../../mov_est/gtg/
	./lagrangian_est_gtg
	
	cd ../../../simeval/SimResult_gtg_5areas_tau04
	
	#source(simRun.r) | R --vanilla --slave
	
	Rscript simRun.r
	
	cd ../../admb
	((counter++))

 done


