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

sim_gtg <- read.rep("../admb/OM/gtg/lagrangian_OM_gtg.rep")
est_gtg <- read.rep("../admb/mov_est/gtg/firstrun.rep")

names(sim_gtg)
names(est_gtg)

dim(est_gtg$VulB)
head((est_gtg$VulB))
summary(est_gtg$VulB)

dim(sim_gtg$VulB)

vbestc<-(est_gtg$VulB)
vbsimc<-(sim_gtg$VulB[which(sim_gtg$VulB[,1]>(12*(sim_gtg$rep_yr-30))),])


dim(vbestc)
dim(vbsimc)

head(vbestc)
head(vbsimc)


vbestc<-rename(as.data.frame(vbestc), c("V1"="time", "V2"= "group","V3"="1", "V4"="2", "V5"= "3",
  "V6"="4", "V7"="5", "V8"= "6","V9"="7", "V10"="8", "V11"= "9","V12"="10", "V13"="11", "V14"= "12",
  "V15"="13", "V16"="14", "V17"= "15","V18"="16", "V19"="17", "V20"= "18","V21"="19", "V22"="20"))
vbsimc<-rename(as.data.frame(vbsimc), c("V1"="time", "V2"= "group","V3"="1", "V4"="2", "V5"= "3",
  "V6"="4", "V7"="5", "V8"= "6","V9"="7", "V10"="8", "V11"= "9","V12"="10", "V13"="11", "V14"= "12",
  "V15"="13", "V16"="14", "V17"= "15","V18"="16", "V19"="17", "V20"= "18","V21"="19", "V22"="20"))


vbestm<-melt(vbestc,id=c("time","group"))
head(vbestm)
vbestm$time=vbestm$time+480
vbsimm<-melt(vbsimc,id=c("time","group"))

head(vbsimm)

totvbestm<-aggregate(vbestm$value, by=list(vbestm$group, vbestm$time ),sum)
totvbsimm<-aggregate(vbsimm$value, by=list(vbsimm$group, vbsimm$time ),sum)


totvbsimm$type<-"simulated"

totvbestm$type<-"estimated"


dfa<-rbind(totvbestm,totvbsimm)
summary(dfa2)
dfa[1:100,]

dfa2<-rename(as.data.frame(dfa),c("Group.1"= "group", "Group.2"="time","x"="vb","type"="type"))
summary(dfa2)

pa<-ggplot(dfa2, aes(x=time,y=vb,color=type))
pa<-pa+geom_line()
pa<-pa+facet_wrap(~group,scale="free")
pa<-pa+geom_vline(aes(xintercept=12*70)pa
pa


dfb<-aggregate(dfa2$vb,by=list(dfa2$time,dfa2$type), sum)

pb<-ggplot(dfb, aes(x=Group.1,y=x,color=Group.2))
pb<-pb+geom_line()
pb


#==================================
#Nage
dim(sim_gtg$Nage)
dim(est_gtg$Nage)

nestc<-(est_gtg$Nage[which(est_gtg$VulB[,1]>(12*30)),])
nsimc<-(sim_gtg$Nage[which(sim_gtg$VulB[,1]>(12*sim_gtg$rep_yr)),])


dim(nestc)
dim(nsimc)

head(nestc)
head(nsimc)


vbestc<-rename(as.data.frame(vbestc), c("V1"="time", "V2"= "group","V3"="1", "V4"="2", "V5"= "3",
  "V6"="4", "V7"="5", "V8"= "6","V9"="7", "V10"="8", "V11"= "9","V12"="10", "V13"="11", "V14"= "12",
  "V15"="13", "V16"="14", "V17"= "15","V18"="16", "V19"="17", "V20"= "18","V21"="19", "V22"="20"))
vbsimc<-rename(as.data.frame(vbsimc), c("V1"="time", "V2"= "group","V3"="1", "V4"="2", "V5"= "3",
  "V6"="4", "V7"="5", "V8"= "6","V9"="7", "V10"="8", "V11"= "9","V12"="10", "V13"="11", "V14"= "12",
  "V15"="13", "V16"="14", "V17"= "15","V18"="16", "V19"="17", "V20"= "18","V21"="19", "V22"="20"))


vbestm<-melt(vbestc,id=c("time","group"))
head(vbestm)
vbestm$time=vbestm$time+480
vbsimm<-melt(vbsimc,id=c("time","group"))

head(vbsimm)

totvbestm<-aggregate(vbestm$value, by=list(vbestm$group, vbestm$time ),sum)
totvbsimm<-aggregate(vbsimm$value, by=list(vbsimm$group, vbsimm$time ),sum)


totvbsimm$type<-"simulated"

totvbestm$type<-"estimated"


dfa<-rbind(totvbestm,totvbsimm)
summary(dfa2)
dfa[1:100,]

dfa2<-rename(as.data.frame(dfa),c("Group.1"= "group", "Group.2"="time","x"="vb","type"="type"))
summary(dfa2)

pa<-ggplot(dfa2, aes(x=time,y=vb,color=type))
pa<-pa+geom_line()
pa<-pa+facet_wrap(~group,scale="free")
pa

#===========================================================
#compare simple model - sim and est

rm(list=ls()); 
#if (Sys.info()["nodename"] =="sager")  setwd("~/Dropbox/LSRA/length_SRA/sim_est_lsra")
setwd("/Users/catarinawor/Documents/Lagrangian/R")
source("read.admb.R")

setwd("/Users/catarinawor/Documents/Lagrangian/admb")
simple_est<-read.rep("mov_est/simple/firstrun.rep")
#sim_gtg <- read.rep("lagrangian_OM_gtg.rep")
simple_sim <- read.rep("OM/simple/lagrangian_OM.rep")
simple_est$Nage[1:15,]
simple_sim$Nage[1:15,]
simple_est$Effarea[4:12,-c(1:9)]
simple_sim$Effarea[4:12,-c(1:9)]


par(mfrow=c(1,2))
plot(simple_est$SB[(70*12+1):length(simple_est$SB)], type="l", main="est") 
lines(simple_sim$SB[(simple_sim$rep_yr*12+1):length(simple_sim$SB)], col="red", main="sim",type="l")
plot(simple_sim$SB[(simple_sim$rep_yr*12+1):length(simple_sim$SB)], col="red", main="sim",type="l")
lines(simple_est$SB[(70*12+1):length(simple_est$SB)], type="l", main="est")
plot(simple_est$SB[(70*12+1):length(simple_est$SB)], type="l", main="est")
lines(simple_sim$SB[(simple_sim$rep_yr*12+1):length(simple_sim$SB)], col="red", main="sim",type="l")






#=============================================================
#gtg model
rm(list=ls()); 
#if (Sys.info()["nodename"] =="sager")  setwd("~/Dropbox/LSRA/length_SRA/sim_est_lsra")
setwd("/Users/catarinawor/Documents/Lagrangian/R")
source("read.admb.R")

setwd("/Users/catarinawor/Documents/Lagrangian/admb")
gtg_est<-read.rep("mov_est/gtg/firstrun.rep")
#sim_gtg <- read.rep("lagrangian_OM_gtg.rep")
gtg_sim <- read.rep("OM/gtg/lagrangian_OM_gtg.rep")
gtg_est$Nage[1:15,]
gtg_sim$Nage[1:15,]
gtg_est$Effarea[4:12,-c(1:9)]
gtg_sim$Effarea[4:12,-c(1:9)]

names(gtg_sim)
names(gtg_est)
length(gtg_est$SB)
length(gtg_sim$SB)

plot(gtg_est$SB[(70*12+1):length(gtg_est$SB)], type="l", main="est",lwd=2)
lines(gtg_sim$SB[(70*12+1):length(gtg_sim$SB)], col="red", main="sim",type="l")




