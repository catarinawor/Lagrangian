sssource("../../read.admb.R")



read.pin <- function(fn)
{
	# The following reads is catarina's effort to read in the pin file
	# Created By Steven Martell
	options(warn=-1)  #Suppress the NA message in the coercion to double
	lines<- readLines(fn)
	
	myline<-seq(2, length(lines), by=2)

	A=list()
	for(i in 1:length(myline)){
		
		A[[i]] <- as.vector(read.table(text=lines[myline[i]], header=FALSE))
	}
 
	return(A)
}




readOutput <- function(dir)
{
	setwd(dir)
	sim <- read.rep("OM/gtg/CL_gtg.rep")	
	HCR <- read.pin("OM/gtg/HCR.dat")
	seed<- scan("OM/gtg/seed.txt")
	C <- list(sim,HCR,seed)
	return( C );
}



	#seed<-scan("../../../admb/OM/gtg/seed.txt")	
	#HCR <- read.pin("../../../admb/OM/gtg/HCR.dat")
	sims<-readOutput("../../../admb/")
	slo<-round(sims[[2]][[4]]*100)
	int<-round(sims[[2]][[5]]*10)
	file.name <- paste("HCR",sims[[3]],"_0",slo,"_0",int,".Rdata",sep="")
	#setwd("../R/OM_dat/res_search")
	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_BesT_nostRec")
	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_stquo_nostRec")
	
	#put a 600 thousand ton limit on TAC.
	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_stquo_nostRec_limit")	
	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_stquo_onestRec_limit")	
	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_stquo_twostRec_limit")	
		

	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_Best_nostRec_limit")
	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_Best_onestRec_limit")
	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_Best_twostRec_limit")
	
	#no cap	
	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_Best_nostRec")
	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_Best_onestRec")
	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_Best_twostRec")

	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_stquo_nostRec")	
	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_stquo_onestRec")	
	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_stquo_twostRec")	
		
	
	#A50=3.0 , sd=1	
	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_Best_nostRec_m2")
	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_Best_onestRec_m2")
	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_Best_twostRec_m2")

	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_stquo_nostRec_m2")	
	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_stquo_onestRec_m2")	
	#setwd("/Volumes/3T_dom_media/Catarina/HCR_sims_stquo_twostRec_m2")	
		
	#on angrycrab
	#setwd("/home/wisdom/result/HCR_sims_Best_nostRec_m2")
	#setwd("/home/wisdom/result/HCR_sims_Best_onestRec_m2")
	#setwd("/home/wisdom/result/HCR_sims_Best_twostRec_m2")

	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_sims_Best_twostRec_m2")


	#setwd("/home/wisdom/result/HCR_sims_stquo_nostRec_m2")	
	#setwd("/home/wisdom/result/HCR_sims_stquo_onestRec_m2")	
	#setwd("/home/wisdom/result/HCR_sims_stquo_twostRec_m2")	
		
	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_sims_stquo_nostRec_m2")
	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_sims_stquo_onestRec_m2")
	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_sims_stquo_twostRec_m2")

	#zoology cluster stuff
	#setwd("~/results/HCR_linear_nostrec")
	#setwd("~/results/HCR_linear_onestrec")
	#setwd("~/results/HCR_linear_twostrec")
	#setwd("~/results/HCR_100_twostrec/")


	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_linear_nostRec")
	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_linear_onestRec")
	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_linear_twostRec")

	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_stquo_nostRec_100")
	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_stquo_onestRec_100")
	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_stquo_twostRec_100")

	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_linear_nostRec_100")
	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_linear_onestRec_100")
	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_linear_twostRec_100")

	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_stquo_nostRec_100_m2")
	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_stquo_onestRec_100_m2")
	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_stquo_twostRec_100_m2")

	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_linear_nostRec_100_m2")
	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_linear_onestRec_100_m2")
	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_linear_twostRec_100_m2")

	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_linear_nostRec_100_m3")
	
	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_stquo_onestRec_100_m2")
	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_stquo_twostRec_100_m2")
	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_stquo_onestRec_100_m3")
	setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_stquo_twostRec_100_m3")

	#setwd("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_stquo_nostRec_100_m3")


	save(sims,file=file.name)




	











