rm(list=ls()); 

library(plyr)
library(data.table)
library(ggplot2)
library(reshape2)
#if (Sys.info()["nodename"] =="sager")  setwd("~/Dropbox/LSRA/length_SRA/sim_est_lsra")
setwd("/Users/catarinawor/Documents/Lagrangian")
source("R/read.admb.R")


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

titulos<-c("5 areas, tau=1.0, B = 1.0","5 areas, tau= 1.0, B = 2.0",
"3 areas, tau=1.0, B = 1.0","3 areas, tau=1.0, B = 2.0",
"5 areas, tau=0.4, B = 1.0","5 areas, tau= 0.4, B = 2.0",
"3 areas, tau=0.4, B = 1.0","3 areas, tau=0.4, B = 2.0")

setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")
#setwd("/Users/catarinawor/Documents/hake/Lag_Model_paper/")
pdf("Figure3.pdf", width=14, height=7)
par(mfcol=c(2,4))
for( i in 1:length(.SIMDIRS)){
  boxplot(pbias[[i]],names=c(expression("t"[0]),expression("CV"),expression("a"[50]),
    expression(sigma["X"["max"]])),ylim=c(-10,10),main=titulos[i],cex.axis=1.5,
    cex.lab=2,cex.main=2,cex=1.6)
  abline(h=0)
  text(4, y = 8, labels = nrow(pbias[[i]]), cex=2)
}
mtext("Relative Error", 2, line = -2, outer = TRUE, font=2)
dev.off()

#=======================================================================
#gtg
#=======================================================================

rm(list=ls()); 

#if (Sys.info()["nodename"] =="sager")  setwd("~/Dropbox/LSRA/length_SRA/sim_est_lsra")
setwd("/Users/catarinawor/Documents/Lagrangian/")
source("R/read.admb.R")

sim_gtg <- read.rep("admb/OM/gtg/lagrangian_OM_gtg.rep")
#est_gtg <- read.rep("lagrangian_est_gtg.rep")




true_pars_gtg <- c(sim_gtg$"mo",sim_gtg$"cvPos",sim_gtg$"maxPos50",sim_gtg$"maxPossd")

indpar<-c(1,2,3,4)



.SIMDIRSGTG   <- c("/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_gtg_5areas_tau1",
                   "/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_gtg_5areas_tau1_delta2",
                   "/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_gtg_3areas_tau1",
                   "/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_gtg_3areas_tau1_delta2",
                   "/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_gtg_5areas_tau04",
                   "/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_gtg_5areas_tau04_delta2",
                   "/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_gtg_3areas_tau04",
                   "/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_gtg_3areas_tau04_delta2")



.SIMNAME<-list(length(.SIMDIRSGTG))

estn_gtg<-list(length(.SIMDIRSGTG))
pbias_gtg<-list(length(.SIMDIRSGTG))
maxgrad_gtg<-list(length(.SIMDIRSGTG))
initvals_gtg<-list(length(.SIMDIRSGTG))
initvals_bad_gtg<-list(length(.SIMDIRSGTG))

for( i in 1:length(.SIMDIRSGTG)){
  .SIMNAME[[i]]   <- list.files(.SIMDIRSGTG[i],pattern="\\.Rdata", full.name=TRUE)
  
  tmp_estn<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=4)
  tmp_pbias<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=length(true_pars_gtg))
  tmp_maxgrad<-vector(length=length(.SIMNAME[[i]]))
  tmp_initvals<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=4)
  

  for( j in 1:length(.SIMNAME[[i]])){
    load(.SIMNAME[[i]][j])

    #parameters
    tmp_estn[j,]<-exp(sims[[3]]$est)
    #tmp_pbias[j,]<-((round(tmp_estn[j,],2)-true_pars_gtg)/true_pars_gtg)*100
    
    for(a in 1:(length(true_pars_gtg))){
        tmp_pbias[j,a]<-((tmp_estn[j,indpar[a]]-true_pars_gtg[a])/true_pars_gtg[a])*100
    }

    

    tmp_maxgrad[j]<-sims[[3]]$maxgrad
    tmp_initvals[j,]<-exp(unlist(sims[[5]][1:4]))
   }

  tmp_estn<- tmp_estn[tmp_maxgrad<=1.0000e-03,]
  tmp_pbias<- tmp_pbias[tmp_maxgrad<=1.0000e-03,]
  
  estn_gtg[[i]]<-tmp_estn
  pbias_gtg[[i]]<-tmp_pbias
  maxgrad_gtg[[i]]<-tmp_maxgrad
  initvals_gtg[[i]]<-tmp_initvals[tmp_maxgrad<=1.0000e-04,]
  initvals_bad_gtg[[i]]<-tmp_initvals[tmp_maxgrad>1.0000e-04,]


}

indAyr<-rep(1:3,1200)[2161:3600]
titulos<-c("5 areas, tau=1.0, B = 1.0","5 areas, tau= 1.0, B = 2.0",
"3 areas, tau=1.0, B = 1.0","3 areas, tau=1.0, B = 2.0",
"5 areas, tau=0.4, B = 1.0","5 areas, tau= 0.4, B = 2.0",
"3 areas, tau=0.4, B = 1.0","3 areas, tau=0.4, B = 2.0")


#setwd("/Users/catarinawor/Documents/hake/Thesis/figs")
#setwd("/Users/catarinawor/Documents/hake/JTC_talk")
#pdf("quatroscn.pdf")

setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")
pdf("GTG_version_simeval.pdf", width=14, height=7)
par(mfcol=c(2,4))
for( i in 1:length(.SIMDIRSGTG)){
  boxplot(pbias_gtg[[i]],names= c(expression("t"[0]),expression("CV"),expression("a"[50]),
    expression(sigma["X"["max"]])),ylim=c(-50,60),main=titulos[i],
  cex.axis=1.5,cex.lab=2,cex.main=2,cex=1.6)
  abline(h=0)
  text(4, y = -28, labels = nrow(pbias_gtg[[i]]),cex=2)
}
mtext("Relative Error", 2, line = -2, outer = TRUE, font=2)
dev.off()






#======================================


