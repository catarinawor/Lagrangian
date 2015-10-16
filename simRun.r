#==============================================
#Title: SimRun
#Author: Catarina Wor
#date: Sep 26 2015
#shameful code to run simulations from R because 
#it's staurday and I'm too lazy to write a makefile
#
#==============================================
setwd("/Users/catarinawor/Documents/Lagrangian/")
source("read.admb.R")



##readOutput <- function(dir)
##{
##	setwd(dir)
##	sim = read.rep("lagrangian_OM.rep")
##	est = read.rep("lagrangian_est.rep")
##	C <- c(sim,est)
##	return( C );
##}
##
##
##	seed<-scan("seed.txt")	
##	file.name <- paste("simest",seed,".Rdata",sep="")
##	sims<-readOutput("/Users/catarinawor/Documents/Lagrangian/")
##	setwd("/Users/catarinawor/Documents/Lagrangian/SimResult")
##	save(sims,file=file.name)




readOutput <- function(dir)
{
	setwd(dir)
	sim = read.rep("lagrangian_OM.rep")
	est = read.rep("lagrangian_est.rep")
	par = read.fit("lagrangian_est")
	C <- list(sim,est,par)
	return( C );
}



	seed<-scan("seed.txt")	
	file.name <- paste("simest",seed,".Rdata",sep="")
	sims<-readOutput("/Users/catarinawor/Documents/Lagrangian/")
	setwd("/Users/catarinawor/Documents/Lagrangian/SimResult")
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






