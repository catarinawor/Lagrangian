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

#================================================================
rm(list=ls()); 
#if (Sys.info()["nodename"] =="sager")  setwd("~/Dropbox/LSRA/length_SRA/sim_est_lsra")
setwd("/Users/catarinawor/Documents/Lagrangian/R")
source("read.admb.R")

setwd("/Users/catarinawor/Documents/Lagrangian/admb")
simple_est<-read.rep("mov_est/simple/firstrun.rep")
#sim_gtg <- read.rep("lagrangian_OM_gtg.rep")
simple_sim <- read.rep("OM/simple/lagrangian_OM.rep")

names(simple_est)
names(simple_sim)

par(mfrow=c(1,2))
plot(simple_est$SB[(20*12+1):length(simple_est$SB)], type="l", main="est") 
lines(simple_sim$SB[(simple_sim$rep_yr*12+1):length(simple_sim$SB)], col="red", main="sim",type="l")
plot(simple_sim$SB[(simple_sim$rep_yr*12+1):length(simple_sim$SB)], col="red", main="sim",type="l")
lines(simple_est$SB[(20*12+1):length(simple_est$SB)], type="l", main="est") 

simple_est$Nage[1:15,]
simple_sim$Nage[1:15,]

simple_est$Effarea[4:12,-c(1:9)]
simple_sim$Effarea[4:12,-c(1:9)]


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


