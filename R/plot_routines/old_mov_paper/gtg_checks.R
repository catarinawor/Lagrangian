#========================================================================================
# Biomass comparison between gtg and non gtg
#========================================================================================

#need to run with no effort
rm(list=ls()); 

library(plyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(animation)
library(ggmap)
#if (Sys.info()["nodename"] =="sager")  setwd("~/Dropbox/LSRA/length_SRA/sim_est_lsra")
setwd("/Users/catarinawor/Documents/Lagrangian/R")
source("read.admb.R")

#sim_gtg <- read.rep("../admb/OM/gtg/lagrangian_OM_gtg.rep")
sim_gtg <- read.rep("../admb/OM/gtg/CL_gtg.rep")


names(sim_gtg)
apply(sim_gtg$histTAC,1,sum)
apply(sim_gtg$yYieldNat[100:130,],1,sum)/apply(sim_gtg$histTAC,1,sum)

sim_gtg$yYieldNat[100:130,]/sim_gtg$histTAC



plot(apply(sim_gtg$yNage,1,sum),type="l")

plot(sim_gtg$ytB,type="l")
abline(h=sim_gtg$Bo)

plot(sim_gtg$fmsy,type="l")


matplot(t(sim_gtg$comm_obsCatLen), type="b")

dim(sim_gtg$yNage)



dim(sim_gtg$comm_obsCatLen)

sim_gtg$comm_obsCatLen[1:49,]
plot(sim_gtg$SB[sim_gtg$indyr>51], type="l")

sim_gtg$wa
sim_gtg$ywa


plot(apply(sim_gtg$yNage[51:102,],1,sum),type="l")
lines(apply(sim_gtg$comm_obsCatLen,1,sum),type="l",lwd=2, col="blue")
lines(apply(sim_gtg$comm_obsCatage,1,sum),type="l",lwd=2, col="red")

plot(apply(sim_gtg$yNage[52:102,],1,sum),type="l")
lines(apply(sim_gtg$yCatchtotalage[52:102,],1,sum),type="l", col="blue")
lines(apply(sim_gtg$comm_obsCatLen,1,sum),type="l",lwd=2, col="darkgreen")
lines(apply(sim_gtg$comm_obsCatage,1,sum),type="l",lwd=2, col="red")



plot(sim_gtg$ywa)
lines(sim_gtg$wa,type="l",lwd=2, col="red")


sim_gtg$Wage

sim_gtg$Lage


dim(est_gtg$VulB)
head((est_gtg$VulB))
summary(est_gtg$VulB)
