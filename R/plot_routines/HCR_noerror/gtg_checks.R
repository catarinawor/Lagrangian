#========================================================================================
# Check output from gtg 
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

sim_gtg <- read.rep("../admb/OM/gtg/lagrangian_OM_gtg.rep")
hcr <- read.rep("../admb/OM/gtg/CL_gtg.rep")

names(sim_gtg)

sim_gtg$yYieldNat






dim(sim_gtg$Effarea)

sim_gtg$indyr

agg_eff<-aggregate(sim_gtg$Effarea,list(sim_gtg$indyr),sum)
plot(apply(agg_eff[,-1],1,sum))





dim(sim_gtg$VulB)

vbsimc<-(sim_gtg$VulB[which(sim_gtg$VulB[,1]>(12*(sim_gtg$rep_yr-30))),])





vbsimc<-rename(as.data.frame(vbsimc), c("V1"="time", "V2"= "group","V3"="1", "V4"="2", "V5"= "3",
  "V6"="4", "V7"="5", "V8"= "6","V9"="7", "V10"="8", "V11"= "9","V12"="10", "V13"="11", "V14"= "12",
  "V15"="13", "V16"="14", "V17"= "15","V18"="16", "V19"="17", "V20"= "18","V21"="19", "V22"="20"))


vbsimm<-melt(vbsimc,id=c("time","group"))

head(vbsimm)

totvbsimm<-aggregate(vbsimm$value, by=list( vbsimm$time ),sum)


totvbsimm$type<-"simulated"



dfa2<-rename(as.data.frame(totvbsimm),c("Group.1"= "time","x"="vb","type"="type"))
summary(dfa2)

dim(dfa2)

names(sim_gtg)
plot(sim_gtg$ytB, type="b")
abline(h=hcr$Bo)
abline(h=hcr$Bo*.40, lty=2)
abline(h=hcr$Bo*.10, lty=2)


sim_gtg$ytB/(hcr$Bo)

hcr$yNage
sim_gtg$

dfb<-aggregate(dfa2$vb,by=list(dfa2$time,dfa2$type), sum)

pb<-ggplot(dfb, aes(x=Group.1,y=x,color=Group.2))
pb<-pb+geom_line()
pb





#m=.2
#prop_ng=c( 0.05,0.1,0.2,0.3,0.2,0.1,0.05)
#
#lxo<-NULL
#lxo[1]<-1
#for(a in 2:20){
#	lxo[a]<-lxo[a-1]*exp(-m)
#}
#lxo[20]=lxo[20]/(1-exp(-m))


#lxg<-matrix(NA,ncol=length(prop_ng),nrow=20)
#for(g in 1:length(prop_ng)){

#lxg[1,g]<-1*prop_ng[g]
#for(a in 2:20){
#	lxg[a,g]<-lxg[a-1,g]*exp(-m)
#}
#lxg[20,g]=lxg[20,g]/(1-exp(-m))
#}
#
#apply(lxg,1,sum)