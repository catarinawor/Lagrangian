#=======================================================
#comparisons of gtg and non-gtg-estimation models
#=======================================================

rm(list=ls()); 
#if (Sys.info()["nodename"] =="sager")  setwd("~/Dropbox/LSRA/length_SRA/sim_est_lsra")
setwd("/Users/catarinawor/Documents/Lagrangian/R")
source("read.admb.R")


#sim_gtg <- read.rep("lagrangian_OM_gtg.rep")
sim <- read.rep("../admb/OM/gtg/CL_gtg.rep")



names(sim)

apply(sim$yNage,1,sum)


sim$yNage[,1]


length(sim$ytB)

plot(sim$ytB, type="b", lwd=2, ylim=c(sim$Bo*0.05,sim$Bo*1.5))
abline(h=sim$Bo)
abline(h=sim$Bo*0.4,lty=3)
abline(h=sim$Bo*0.1,lty=3)



#================================================================
