#!/usr/local/bin/Rscript --slave
args <- commandArgs(trailingOnly=TRUE)


source("/Users/catarinawor/Documents/Lagrangian/R/OM_dat/OM_dat.R")
source("/Users/catarinawor/Documents/Lagrangian/R/OM_dat/reset_seed.R")





if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else{
	writeHCRdat(slope=args[1],intercept=args[2])
	reset_seed()
}





#nsim=100
#slopes<-seq(0.05,1,by=0.05)
#interceptos<-seq(0.05,1,by=0.05)



#for(sp in 1:length(slopes)){
#	for(it in 1:length(interceptos)){
#		writeHCRdat(DIR= "/Users/catarinawor/Documents/Lagrangian/admb/OM/gtg/",slope=slopes[sp],intercept=interceptos[it])
		
#		for(i in 1:nsim){
				
			#run lagrangian

			#read data
#			setwd("/Users/catarinawor/Documents/Lagrangian/R")3
#			source("read.admb.R")

#			sim_gtg <- read.rep("../admb/OM/gtg/CL_gtg.rep")


#		}
#	}
#}




	











