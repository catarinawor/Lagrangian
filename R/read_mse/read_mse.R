#==============================================
#Title: read_mse
#Author: Catarina Wor
#date: Oct 7th 2016
#Code to read admb assessment outputs and store them for each assessment run
#
#==============================================
setwd("/Users/catarinawor/Documents/Lagrangian/R")
source("read.admb.R")


readOutput <- function(dir)
{
	setwd(dir)
	assmt <- read.rep("lagrangian_SA.rep")
	par <- read.fit("lagrangian_SA")	
	C <- list(assmt,par)
	return( C );
}


	setwd("/Users/catarinawor/Documents/Lagrangian/admb/")
	seed<- scan("seed.txt")

	setwd("/Users/catarinawor/Documents/Lagrangian/admb/stock_assessment")
	tmp <- read.rep("lagrangian_SA.rep")
	assmt_nyr<- tmp$nyr	
	
	sims<-readOutput("/Users/catarinawor/Documents/Lagrangian/admb/stock_assessment")
	setwd("/Users/catarinawor/Documents/Lagrangian/R/read_mse/result/")

	file.name <- paste("assess_yr",assmt_nyr,".Rdata",sep="")
	
	system(paste("mkdir",paste("run","SA",seed,sep="")))
	setwd(paste("/Users/catarinawor/Documents/Lagrangian/R/read_mse/result/","run","SA",seed,sep=""))


	save(sims,file=file.name)



##argsim = paste("./lagrangian_OM") 
##argest = paste("./lagrangian_est") 
###
##for(i in 1:5){
##	setwd("/Users/catarinawor/Documents/Lagrangian")
##	system(argsim)
##	system(argest)
##	seed<-scan("seed.txt")	
##	file.name <- paste("simest",seed,".Rdata",sep="")
##	sims<-readOutput("/Users/catarinawor/Documents/Lagrangian/")
##	setwd("/Users/catarinawor/Documents/Lagrangian/SimResult")
##	save(sims,file=file.name)
##}





#argsim = paste("./lagrangian_OM") 
#argest = paste("./lagrangian_est") 
#
#for(i in 1:2){
#
#	system(argsim)
# 	system(argest)
#	seed<-scan("seed.txt")	
#	file.name <- paste("simest",seed,".Rdata",sep="")
#	sims<-readOutput("/Users/catarinawor/Documents/Lagrangian/")
#	setwd("/Users/catarinawor/Documents/Lagrangian/SimResult")
#	save(sims,file=file.name)
#}
#
#	





#run.Simulation=function(N=1){
#	true_pars<-matrix(NA, nrow=N)
#	est_pars <-matrix(NA, nrow=N)
#	
#	argsim = paste("./lagrangian_OM") 
#	argest = paste("./lagrangian_est") 
#	
#	for(i in 1:N){
#		system(argsim)
#		system(argest)
#		sim = read.rep("lagrangian_OM.rep")
#		est = read.rep("lagrangian_est.rep")
#
#		true_pars[i,] <- c(sim$"mo", exp(sim$"log_tau_c"),sim$"maxPos50",sim$"maxPossd",sim$"cvPos")  
#		est_pars[i,] <- c(est$"mo",exp(est$"log_tau_c"),est$"maxPos50",est$"maxPossd",est$"cvPos")
#    
#    }
#}



#sim = read.rep("lagrangian_OM.rep")
#est = read.rep("lagrangian_est.rep")
#
#nomes <- names(sim)
