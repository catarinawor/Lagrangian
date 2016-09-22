#=======================================================
#comparisons of gtg and non-gtg-estimation models
#=======================================================

rm(list=ls()); 
#if (Sys.info()["nodename"] =="sager")  setwd("~/Dropbox/LSRA/length_SRA/sim_est_lsra")
setwd("/Users/catarinawor/Documents/Lagrangian/")
source("read.admb.R")


simple<-read.rep("lagrangian_est.rep")
#sim_gtg <- read.rep("lagrangian_OM_gtg.rep")
sim <- read.rep("lagrangian_OM.rep")
gtg_sim <- read.rep("lagrangian_OM_gtg.rep")

gtg <- read.rep("lagrangian_est_gtg.rep")

gtg_sim$rep_yr

plot(gtg_sim$SB[(gtg_sim$rep_yr*12+1):length(gtg_sim$SB)], type="l") 
lines(gtg$SB[(20*12+1):length(gtg$SB)], col="red")

cbind(
gtg$SB[1:10],
gtg_sim$SB[1:10])

head(gtg$Nage)


plot(apply(gtg_sim$Nage[(gtg_sim$rep_yr*12+1):length(gtg_sim$SB),-1],1,sum), type="l") 
lines(apply(gtg$Nage[,-1],1,sum), col="red")

names(sim)
names(gtg_sim)
names(gtg)


head(sim$Nage)
head(gtg_sim$Nage)

plot(apply(sim$Nage,1,sum), type="l") 
lines(apply(gtg_sim$Nage[-1,],1,sum), col="red")

plot(sim$SB, type="l") 
lines(gtg_sim$SB, col="red")

