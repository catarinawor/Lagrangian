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

sim <- read.rep("../admb/OM/simple/lagrangian_OM.rep")
est<- read.rep("../admb/mov_est/simple/firstrun.rep")

names(sim)
names(est)

dim(est$VulB)
head((est$VulB))
summary(est$VulB)

dim(sim$VulB)
head((sim$VulB))


vbestc<-apply(est$VulB,1,sum)
vbsimc<-apply(sim$VulB,1,sum)

plot(vbestc, ylim=c(0,3))
lines(vbsimc, lwd=2)


plot(est$yCatchtotalobs )
lines(est$yCatchtotal, lwd=2)

dim(vbestc)
dim(vbsimc)

head(vbestc)
head(vbsimc)







totvbestm<-agvbestm<-melt(vbestc,id=c("time","group"))
head(vbestm)
vbestm$time=vbestm$time+480
vbsimm<-melt(vbsimc,id=c("time","group"))

head(vbsimm)
gregate(vbestm$value, by=list(vbestm$group, vbestm$time ),sum)
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
plot(simple_est$SB[(20*12+1):length(simple_est$SB)], type="l", main="est") 
lines(simple_sim$SB[(simple_sim$rep_yr*12+1):length(simple_sim$SB)], col="red", main="sim",type="l")
plot(simple_sim$SB[(simple_sim$rep_yr*12+1):length(simple_sim$SB)], col="red", main="sim",type="l")
lines(simple_est$SB[(20*12+1):length(simple_est$SB)], type="l", main="est")
plot(simple_est$SB[(20*12+1):length(simple_est$SB)], type="l", main="est")
lines(simple_sim$SB[(simple_sim$rep_yr*12+1):length(simple_sim$SB)], col="red", main="sim",type="l")




########################################################################
#simeval one scenario




setwd("/Users/catarinawor/Documents/Lagrangian/admb/OM/simple")
sim <- read.rep("lagrangian_OM.rep")

sim$tau_c
nomes <- names(sim)


#parameter estimates
.SIMDIRS   <- c("/Users/catarinawor/Documents/Lagrangian/simeval/SimResult_5areas_tau04")

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

    true_pars <- c(sim$"mo",sim$"cvPos",sim$"maxPos50",sim$"maxPossd",mean(sim$Fmult))  


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


titulos<-c("5 areas, tau=1.0, B = 0.4")

#setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")
#setwd("/Users/catarinawor/Documents/hake/JTC_talk")
#pdf("single_version_simeval.pdf", width=14, height=7)
par(mfcol=c(1,1))
for( i in 1:length(.SIMDIRS)){
  boxplot(pbias[[i]],names=c(expression("t"[0]),expression("CV"),expression("a"[50]),
    expression(sigma["X"["max"]]),expression("q")),ylim=c(-10,10),main=titulos[i],cex.axis=1.5,
    cex.lab=2,cex.main=2,cex=1.6)
  abline(h=0)
  text(4, y = 8, labels = nrow(pbias[[i]]), cex=2)
}
mtext("% Bias", 2, line = -2, outer = TRUE, font=2)
#dev.off()


#===============================================================
#come up with an indicator of the impact of biased parameter estimates on movement
 library(ggplot2)
dim(sim$VBarea)
head(sim$VBarea)

indmonth<-rep(1:12,100)
indyr<-rep(1:100,each=12)
length(indmonth)
barplot(sim$VBarea[1200-5,])


head(sims[[2]]$VBarea)

ALLDF<-list(length(.SIMDIRS))

for( i in 1:length(.SIMDIRS)){
  .SIMNAME[[i]]   <- list.files(.SIMDIRS[i],pattern="\\.Rdata", full.name=TRUE)
  
  dfall<-NULL
  

  for( j in 1:length(.SIMNAME[[i]])){

    load(.SIMNAME[[i]][j])
    names(sims[[2]])

    dfsi<-cbind(indmonth,indyr,melt(sims[[1]]$VBarea, value.name="vb"))
    names(dfsi)<-c("indmonth", "indyr","tstp","area","vb")
    dfsi1<-dfsi[dfsi$indyr>71,]
    dfsi1$area<-as.factor(as.numeric(dfsi1$area)+29)
   
    dfsi2<-aggregate(dfsi1$vb,by=list(dfsi1$area,dfsi1$indmonth), median)
    names(dfsi2)<-c("area", "month","medianvb")
    dfsi2$type<-"simulated"
    

    dfes<-cbind(indmonth,indyr,melt(sims[[2]]$VBarea, value.name="vb"))
    names(dfes)<-c("indmonth", "indyr","tstp","area","vb")
    dfes1<-dfes[dfes$indyr>71,]
    dfes1$area<-as.factor(as.numeric(dfes1$area)+29)
   
    dfes2<-aggregate(dfes1$vb,by=list(dfes1$area,dfes1$indmonth), median)
    names(dfes2)<-c("area", "month","medianvb")
    dfes2$type<-"estimated"

    dfall<-rbind(dfall,dfsi2,dfes2)
  }
ALLDF[[i]]<-dfall
}


summary(dfall)


p<-ggplot(dfall, aes(x=as.factor(area), y=medianvb, color=type))
p<-p+geom_boxplot()
p<-p+geom_vline(xintercept = 49-29)
p<-p+facet_wrap(~month)
p<-p+theme_bw(16)
p




dfe<-cbind(indmonth,indyr,melt(sim$VBarea, value.name="vb"))
names(dfe)<-c("indmonth", "indyr","tstp","area","vb")
dfe1<-dfe[dfe$indyr>71,]
dfe1$area<-as.factor(as.numeric(dfe1$area)+29)
head(dfe1)


aggregate(dfe1$vb,by=list(dfe1$area,dfe1$indmonth), median)

p<-ggplot(dfe1, aes(x=(area), y=vb))
p<-p+geom_boxplot()
p<-p+geom_vline(xintercept = 49-29)
p<-p+facet_wrap(~indmonth)
p

boxplot(sim$VBarea[which(indmonth==7)[71:100],])

plot(apply(sim$VBarea[which(indmonth==7)[71:100],],2,median),type="l")



dfe<-melt(sim$VBarea[which(indmonth==7)[71:100],])
head(dfe)

p<-ggplot(dfe, aes(x=as.numeric(Var2), y=value))
p<-p+geom_line(aes(color=as.factor(Var1)))
p

#doing it for a series of simulations

.SIMNAME[[1]]   <- list.files(.SIMDIRS[i],pattern="\\.Rdata", full.name=TRUE)
  
  length(sims)
    length(sims[[1]])
    names(sims[[1]])
    length(sims[[2]])
    names(sims[[2]])

    sims[[1]]$phiE

finaldf<-NULL
for( j in 1:length(.SIMNAME[[i]])){
    load(.SIMNAME[[1]][j])

    
    indmonth<-rep(1:12,100)

    

    for(m in 1:12){
    	tmpsim<-apply(sims[[1]]$"VBarea"[which(indmonth==m)[71:100],1:26],2,median)
    	tmpest<-apply(sims[[2]]$"VBarea"[which(indmonth==m)[71:100],1:26],2,median)
    	
    	tsdf<-data.frame(value=tmpsim, month=m, src="simulated", area=30:55)
    	tedf<-data.frame(value=tmpest, month=m, src="estimated", area=30:55)

    	finaldf<-rbind(finaldf,tsdf,tedf)
    }
   }

dim(finaldf)
head(finaldf)
summary


p<-ggplot(finaldf, aes(x=as.factor(area), y=value))
p<-p+geom_boxplot( aes(color=src, fill=src))
p<-p+ facet_wrap(~month)
p<-p+ geom_vline(aes(xintercept=which(30:55==49)))
p


