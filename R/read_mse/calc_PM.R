#==============================================
#Title: calc_PM
#Author: Catarina Wor
#date: Oct 21th 2016
#Read results form MSE runs and calculate performanace measures
#
#==============================================

#TODO
#
#OMplots
#
#iSCAM plots
#

 
DIR<-"/Users/catarinawor/Documents/Lagrangian/R/read_mse/result"
setwd(DIR)
iSCAMdir<-list.files(DIR,pattern="runiSCAM",full.name=TRUE)

plotlib<-"/Users/catarinawor/Documents/Lagrangian/R/read_mse/Rplots"
.RFILES     <- list.files(.LIB,pattern="\\.[Rr]$")

OMdat<-list()
SAdat<-list()
SApar<-list()

for(i in 1:length(iSCAMdir)){

	Rfiles <- list.files(iSCAMdir[i],pattern="\\.Rdata",full.name=TRUE)
	rfil<-basename(Rfiles)

	

	load(grep("OM",Rfiles,value=TRUE))
	
	seed<-paste("seed",OM[[1]]$seed,sep="")
	OMdat[[seed]]<-OM[[1]]


	



	SAfiles<-grep("assess",Rfiles,value=TRUE)
	
	indSAdat<-list(length=length(SAfiles))


	for(j in 1:length(SAfiles)){
		
		load(SAfiles[j])
		
		names(sims$assmt)
		yr=paste("yr",sims$assmt$nyr,sep="")

		SAdat[[seed]][[yr]]<-sims$assmt
		SApar[[seed]][[yr]]<-sims$par

	}

}











