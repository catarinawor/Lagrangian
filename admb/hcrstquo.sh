#!/bin/sh
 
#echo "Hello, world!"


	cd ../R/OM_dat/ 
	Rscript  runHCRstquo.r  

	cd ../../admb/ 

	counter=1

	while [ $counter -le 10 ]
	do

	 cd ../R/OM_dat/ 
		Rscript runWTerr.r

	cd ../../admb/OM/gtg 

		./lagrangian_OM_gtg
	
	 cd ../../../R/OM_dat/res_search/ 

	Rscript read_result.R

	
	cd ../../../admb
	
	counter=$(( $counter + 1 ))	

	done


	



