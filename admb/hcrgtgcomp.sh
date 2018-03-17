#!/bin/bash
 
#echo "Hello, world!"


tp=0.5

hrtp=0.05

for i in {1..10}; 
do 

 	sp=`echo "$i * $hrtp"| bc -l`
	it=`echo "$tp"| bc -l`
 	

	echo "$sp"
	
	cd ../R/OM_dat/ 
	Rscript  runHCRsearch.R  $sp  $it

	cd ../../admb/ 

	counter=1

	while [ $counter -le 100 ]
	do

	 cd ../R/OM_dat/ 
		Rscript runWTerr.R

	cd ../../admb/OM/gtg 

		./lagrangian_OM_gtg
	
	 cd ../../../R/OM_dat/res_search/ 

	Rscript read_result.R

	
	cd ../../../admb
	
	counter=$(( $counter + 1 ))	

	done

done

	

