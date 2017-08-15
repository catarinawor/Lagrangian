#==============================================
#Title: Read in simulation outcomes and plot thes
#Author: Catarina Wor
#date: Oct 6 2015
#
#==============================================


library(ggplot2)
library(reshape2)



.PWD        <- "/Users/catarinawor/Documents/Lagrangian/"
.THEME      <- theme_bw(16)

.SIMDIRS   <- list.dirs(paste(.PWD,"simeval",sep=""), full.name=TRUE)







########################################################################
#simeval one scenario


setwd("/Users/catarinawor/Documents/Lagrangian/admb/OM/simple")
sim <- read.rep("lagrangian_OM.rep")

sim$tau_c
nomes <- names(sim)

true_pars <- c(sim$"mo",sim$"cvPos",sim$"maxPos50",sim$"maxPossd")  



#parameter estimates
.SIMDIRS   <- c("/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_5areas_tau1",
  "/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_5areas_tau1_delta2",
  "/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_3areas_tau1",
  "/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_3areas_tau1_delta2",
  "/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_5areas_tau04",
  "/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_5areas_tau04_delta2",
  "/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_3areas_tau04",
  "/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_3areas_tau04_delta2")

.SIMNAME<-list(length(.SIMDIRS))

estn<-list(length(.SIMDIRS))
pbias<-list(length(.SIMDIRS))
maxgrad<-list(length(.SIMDIRS))
initvals<-list(length(.SIMDIRS))
initvals_bad<-list(length(.SIMDIRS))



for( i in 1:length(.SIMDIRS)){
  .SIMNAME[[i]]   <- list.files(.SIMDIRS[i],pattern="\\.Rdata", full.name=TRUE)
  
  tmp_estn<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=4)
  tmp_pbias<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=4)
  tmp_maxgrad<-vector(length=length(.SIMNAME[[i]]))
  tmp_initvals<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=4)
  

  for( j in 1:length(.SIMNAME[[i]])){
    load(.SIMNAME[[i]][j])

    #parameters
    tmp_estn[j,]<-exp(sims[[3]]$est[1:4])
    tmp_pbias[j,]<-((tmp_estn[j,]-true_pars)/true_pars)*100
    tmp_maxgrad[j]<-sims[[3]]$maxgrad
    tmp_initvals[j,]<-exp(unlist(sims[[5]][1:4]))
   }

  tmp_estn<- tmp_estn[tmp_maxgrad<=1.0000e-01,]
  tmp_pbias<- tmp_pbias[tmp_maxgrad<=1.0000e-01,]
  
  estn[[i]]<-tmp_estn
  pbias[[i]]<-tmp_pbias
  maxgrad[[i]]<-tmp_maxgrad
  initvals[[i]]<-tmp_initvals[tmp_maxgrad<=1.0000e-01,]
  initvals_bad[[i]]<-tmp_initvals[tmp_maxgrad>1.0000e-01,]


}


#========================================================================
# Description of the list sims
# [[1]] -> sim- read.rep("lagrangian_OM.rep") 
# [[2]] -> est -read.rep("lagrangian_est.rep")
# [[3]] -> par - read.fit("lagrangian_est")
# [[4]] -> seed
# [[5]] -> pin - in a list
#========================================================================


titulos<-c("5 areas, tau=1.0, B = 1.0","5 areas, tau= 1.0, B = 2.0",
"3 areas, tau=1.0, B = 1.0","3 areas, tau=1.0, B = 2.0",
"5 areas, tau=0.4, B = 1.0","5 areas, tau= 0.4, B = 2.0",
"3 areas, tau=0.4, B = 1.0","3 areas, tau=0.4, B = 2.0")

setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")
#setwd("/Users/catarinawor/Documents/hake/JTC_talk")
pdf("single_version_simeval.pdf", width=14, height=7)
par(mfcol=c(2,4))
for( i in 1:length(.SIMDIRS)){
  boxplot(pbias[[i]],names=c(expression("t"[0]),expression("CV"),expression("a"[50]),
    expression(sigma["X"["max"]])),ylim=c(-10,10),main=titulos[i],cex.axis=1.5,
    cex.lab=2,cex.main=2,cex=1.6)
  abline(h=0)
  text(4, y = 8, labels = nrow(pbias[[i]]), cex=2)
}
mtext("% Bias", 2, line = -2, outer = TRUE, font=2)
dev.off()


#=========================================================================
#gtg
setwd("/Users/catarinawor/Documents/Lagrangian/admb/OM/gtg")
gtg <- read.rep("lagrangian_OM_gtg.rep")


#parameter estimates
.SIMDIRS   <- c("/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_gtg_5areas_tau1")

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

    true_pars <- c(sims[[1]]$"mo",sims[[1]]$"cvPos",sims[[1]]$"maxPos50",sims[[1]]$"maxPossd",mean(sims[[1]]$Fmult-(.1*.1/2)))  


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


titulos<-c("5 areas, tau=1.0, B = 1.0")

#setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")
#setwd("/Users/catarinawor/Documents/hake/JTC_talk")
#pdf("single_version_simeval.pdf", width=14, height=7)
par(mfcol=c(1,1))
for( i in 1:length(.SIMDIRS)){
  boxplot(pbias[[i]],names=c(expression("t"[0]),expression("CV"),expression("a"[50]),
    expression(sigma["X"["max"]]),expression("q")),ylim=c(-50,50),main=titulos[i],cex.axis=1.5,
    cex.lab=2,cex.main=2,cex=1.6)
  abline(h=0)
  text(4, y = 8, labels = nrow(pbias[[i]]), cex=2)
}
mtext("% Bias", 2, line = -2, outer = TRUE, font=2)
#dev.off()



