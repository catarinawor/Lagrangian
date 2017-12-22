#==============================================================
# Read in the results from the HCR search - no assessment error
#Author: Catarina Wor
#Date: Sep 14 2017

#==============================================================


library(cowplot)

#DIR<-"/Volumes/3T_dom_media/Catarina/HCR_sims_coarse_onestRec"
#DIR<-"/Volumes/3T_dom_media/Catarina/HCR_sims_coarse_nostRec"

#DIR<-"/Volumes/3T_dom_media/Catarina/HCR_sims_Best_nostRec"
DIRlim<-"/Volumes/3T_dom_media/Catarina/HCR_sims_Best_nostRec_limit"

#Rfiles <- list.files(DIR,pattern="\\.Rdata",full.name=TRUE)
Rfileslim <- list.files(DIRlim,pattern="\\.Rdata",full.name=TRUE)


plotlib<-"/Users/catarinawor/Documents/Lagrangian/R/plot_routines/HCR_noerror"
plotfiles <- list.files(plotlib,pattern="plot",full.name=TRUE)


	

#SIMSdat<-list()
#for(i in 1:length(Rfiles)){
#	
#	load(Rfiles[i])
#	
#	SIMSdat[[i]]<-sims
#
#}
#
#length(SIMSdat)


SIMSdatlim<-list()


for(i in 1:length(Rfileslim)){
	

	load(Rfileslim[i])
	
	SIMSdatlim[[i]]<-sims


}

length(SIMSdatlim)



#DIRsq<-"/Volumes/3T_dom_media/Catarina/HCR_sims_stquo_nostRec"
DIRsqlim<-"/Volumes/3T_dom_media/Catarina/HCR_sims_stquo_nostRec_limit"


#SIMSsq<-list()
SIMSsqlim<-list()

#Rfilessq <- list.files(DIRsq,pattern="\\.Rdata",full.name=TRUE)
Rfilessqlim <- list.files(DIRsqlim,pattern="\\.Rdata",full.name=TRUE)

#for(y in 1:length(Rfilessq)){
#	
#
#	load(Rfilessq[y])
#	
#	SIMSsq[[y]]<-sims
#}

length(SIMSsq)



for(y in 1:length(Rfilessqlim)){
	

	load(Rfilessqlim[y])
	
	SIMSsqlim[[y]]<-sims

}

length(SIMSsqlim)


#load graphing routinesfor
for(p in 1:length(plotfiles)){
	source(plotfiles[p])
}

#normal<-plot_logUtility( SIMSdat,SIMSsq , sv=TRUE, nome="nostrongrec",nations=TRUE)
lim<-plot_logUtility( SIMSdatlim,SIMSsqlim , sv=TRUE, nome="nostrongrec",nations=TRUE)

#logu<-plot_grid(normal,lim,nrow=2)
#		print(logu)

Ynormal<-plot_Yield( SIMSdat, SIMSsq , sv=TRUE, nome="nostrongrec",nations=TRUE)
Ylim<-plot_Yield( SIMSdatlim, SIMSsqlim , sv=TRUE, nome="nostrongrec",nations=TRUE)

yld<-plot_grid(Ynormal,Ylim,nrow=2)
print(yld)

plot_AAV( SIMSdat, SIMSsq , sv=TRUE, nome="nostrongrec",nations=TRUE)
plot_AAV( SIMSdatlim, SIMSsqlim , sv=TRUE, nome="nostrongrec",nations=TRUE)


plot_closures( SIMSdat, SIMSsq , sv=TRUE, nome="nostrongrec",nations=TRUE)
plot_closures( SIMSdatlim, SIMSsqlim , sv=TRUE, nome="nostrongrec",nations=TRUE)


plot_pTAC( SIMSdat, SIMSsq , SIMSsq , sv=TRUE, nome="nostrongrec",nations=TRUE)


plot_B40( SIMSdat , SIMSsq , sv=TRUE, nome="nostrongrec",nations=TRUE)
plot_B40( SIMSdatlim , SIMSsqlim , sv=TRUE, nome="nostrongrec",nations=TRUE)





#Msq<-SIMSsq
