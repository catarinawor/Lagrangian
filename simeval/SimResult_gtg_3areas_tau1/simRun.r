#==============================================
#Title: SimRun
#Author: Catarina Wor
#date: Sep 26 2015
#shameful code to run simulations from R because 
#it's staurday and I'm too lazy to write a makefile
#
#==============================================
source("../../R/read.admb.R")


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
	sim <- read.rep("OM/gtg/lagrangian_OM_gtg.rep")
	est <- read.rep("mov_est/gtg/lagrangian_est_gtg.rep")
	par <- read.fit("mov_est/gtg/lagrangian_est_gtg")
	guess <- read.pin("mov_est/gtg/lagrangian_est_gtg.pin")
	seed<- scan("OM/gtg/seed.txt")
	C <- list(sim,est,par,seed,guess)
	return( C );
}


	seed<-scan("../../admb/OM/gtg/seed.txt")	
	file.name <- paste("simest",seed,".Rdata",sep="")
	sims<-readOutput("../../admb/")
	setwd("../simeval/SimResult_gtg_3areas_tau1/")
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






