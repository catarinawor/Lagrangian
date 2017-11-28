#!/bin/bash
 
#echo "Hello, world!"


	cd ../R/OM_dat/ 
	Rscript  runHCRstquo.R  

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


	



