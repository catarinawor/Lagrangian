#==============================================
#Title: test the gtg simulation stuff
#Author: Catarina Wor
#date:June 2017
#
#==============================================

rm(list=ls()); 

library(plyr)
library(data.table)
library(ggplot2)
library(reshape2)
#if (Sys.info()["nodename"] =="sager")  setwd("~/Dropbox/LSRA/length_SRA/sim_est_lsra")
setwd("/Users/catarinawor/Documents/Lagrangian")
source("R/read.admb.R")





#parameter estimates
.SIMDIRS   <- c("/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_3areas_tau04")

.SIMNAME<-list(length(.SIMDIRS))

estn<-list(length(.SIMDIRS))
pbias<-list(length(.SIMDIRS))
maxgrad<-list(length(.SIMDIRS))
initvals<-list(length(.SIMDIRS))
initvals_bad<-list(length(.SIMDIRS))



for( i in 1:length(.SIMDIRS)){
  .SIMNAME[[i]]   <- list.files(.SIMDIRS[i],pattern="\\.Rdata", full.name=TRUE)
  
  tmp_estn<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=5)
  tmp_pbias<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=5)
  tmp_maxgrad<-vector(length=length(.SIMNAME[[i]]))
  tmp_initvals<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=5)
  

  for( j in 1:length(.SIMNAME[[i]])){
    
  	load(.SIMNAME[[i]][j])
	true_pars <- c(sims[[1]]$"mo",sims[[1]]$"cvPos",sims[[1]]$"maxPos50",sims[[1]]$"maxPossd",mean(sims[[1]]$"Fmult")) 
    #parameters
    tmp_estn[j,]<-exp(sims[[3]]$est[1:5])
    tmp_pbias[j,]<-((tmp_estn[j,]-true_pars)/true_pars)*100
    tmp_maxgrad[j]<-sims[[3]]$maxgrad
    tmp_initvals[j,]<-exp(unlist(sims[[5]][1:5]))
   }

  tmp_estn<- tmp_estn[tmp_maxgrad<=1.0000e-01,]
  tmp_pbias<- tmp_pbias[tmp_maxgrad<=1.0000e-01,]
  
  estn[[i]]<-tmp_estn
  pbias[[i]]<-tmp_pbias
  maxgrad[[i]]<-tmp_maxgrad
  initvals[[i]]<-tmp_initvals[tmp_maxgrad<=1.0000e-01,]
  initvals_bad[[i]]<-tmp_initvals[tmp_maxgrad>1.0000e-01,]


}

titulos<-c("3 areas, tau=0.4, B = 1.0")
  #"5 areas, tau=1.0, B = 1.0"
#,"5 areas, tau= 1.0, B = 2.0",
#"3 areas, tau=1.0, B = 1.0","3 areas, tau=1.0, B = 2.0",
#"5 areas, tau=0.4, B = 1.0","5 areas, tau= 0.4, B = 2.0",
#,"3 areas, tau=0.4, B = 2.0"

#setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")
#setwd("/Users/catarinawor/Documents/hake/Lag_Model_paper/")
#pdf("Figure3.pdf", width=14, height=7)
#par(mfcol=c(2,4))
for( i in 1:length(.SIMDIRS)){
  boxplot(pbias[[i]],names=c(expression("t"[0]),expression("CV"),expression("a"[50]),
    expression(sigma["X"["max"]]),expression("q")),ylim=c(-50,50),main=titulos[i],cex.axis=1.5,
    cex.lab=2,cex.main=2,cex=1.6)
  abline(h=0)
  text(4, y = 8, labels = nrow(pbias[[i]]), cex=2)
}
mtext("Relative Error", 2, line = -2, outer = TRUE, font=2)
#dev.off()


#============================================================
#gtg
#parameter estimates
.SIMDIRS   <- c("/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_gtg_5areas_tau04")

.SIMNAME<-list(length(.SIMDIRS))

estn<-list(length(.SIMDIRS))
pbias<-list(length(.SIMDIRS))
maxgrad<-list(length(.SIMDIRS))
initvals<-list(length(.SIMDIRS))
initvals_bad<-list(length(.SIMDIRS))



for( i in 1:length(.SIMDIRS)){
  .SIMNAME[[i]]   <- list.files(.SIMDIRS[i],pattern="\\.Rdata", full.name=TRUE)
  
  tmp_estn<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=5)
  tmp_pbias<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=5)
  tmp_maxgrad<-vector(length=length(.SIMNAME[[i]]))
  tmp_initvals<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=5)
  

  for( j in 1:length(.SIMNAME[[i]])){
    
    load(.SIMNAME[[i]][j])
    
  true_pars <- c(sims[[1]]$"mo",sims[[1]]$"cvPos",sims[[1]]$"maxPos50",sims[[1]]$"maxPossd",3) 
    #parameters
    tmp_estn[j,]<-exp(sims[[3]]$est[1:5])
    tmp_pbias[j,]<-((tmp_estn[j,]-true_pars)/true_pars)*100
    tmp_maxgrad[j]<-sims[[3]]$maxgrad
    tmp_initvals[j,]<-exp(unlist(sims[[5]][1:5]))
   }

  tmp_estn<- tmp_estn[tmp_maxgrad<=1.0000e-01,]
  tmp_pbias<- tmp_pbias[tmp_maxgrad<=1.0000e-01,]
  
  estn[[i]]<-tmp_estn
  pbias[[i]]<-tmp_pbias
  maxgrad[[i]]<-tmp_maxgrad
  initvals[[i]]<-tmp_initvals[tmp_maxgrad<=1.0000e-01,]
  initvals_bad[[i]]<-tmp_initvals[tmp_maxgrad>1.0000e-01,]


}

titulos<-c("5 areas, tau=0.4, B = 1.0")
  #"5 areas, tau=1.0, B = 1.0"
#,"5 areas, tau= 1.0, B = 2.0",
#"3 areas, tau=1.0, B = 1.0","3 areas, tau=1.0, B = 2.0",
#"5 areas, tau=0.4, B = 1.0","5 areas, tau= 0.4, B = 2.0",
#,"3 areas, tau=0.4, B = 2.0"

#setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")
#setwd("/Users/catarinawor/Documents/hake/Lag_Model_paper/")
#pdf("Figure3.pdf", width=14, height=7)
#par(mfcol=c(2,4))
for( i in 1:length(.SIMDIRS)){
  boxplot(pbias[[i]],names=c(expression("t"[0]),expression("CV"),expression("a"[50]),
    expression(sigma["X"["max"]]),expression("q")),ylim=c(-50,50),main=titulos[i],cex.axis=1.5,
    cex.lab=2,cex.main=2,cex=1.6)
  abline(h=0)
  text(4, y = 8, labels = nrow(pbias[[i]]), cex=2)
}
mtext("Relative Error", 2, line = -2, outer = TRUE, font=2)
#dev.off()




#=========================================================
#try to come up with a way of plotting movement for all simulations




