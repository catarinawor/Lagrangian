#=========================================================================
# graphing lagrangian model structure and sim-est
# Author: Catarina Wor
# Jun 1st 2015 - updated in November 26th 2015
#=========================================================================

rm(list=ls()); 

library(plyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(animation)
library(ggmap)
#if (Sys.info()["nodename"] =="sager")  setwd("~/Dropbox/LSRA/length_SRA/sim_est_lsra")
setwd("/Users/catarinawor/Documents/Lagrangian/")
source("read.admb.R")



sim <- read.rep("lagrangian_OM.rep")
sim_gtg <- read.rep("lag_OM_gtg_new.rep")
est <- read.rep("lagrangian_est.rep")


nomes <- names(sim)

true_pars <- c(sim$"mo",sim$"cvPos",sim$"maxPos50",sim$"maxPossd")  
est_pars  <- c(est$"mo",est$"cvPos",est$"maxPos50",est$"maxPossd")

#parameter plot
par(mfrow=c(1,1))
barplot(t(matrix(c(true_pars,est_pars),ncol=2)),names.arg = c("mo","cvPos","maxPos50","maxPossd"),beside=T)

#Check simulations plots
names(sim)
names(sim_gtg)


dim(sim$VulB) 

vbgtg<-matrix(0,ncol=length(sim_gtg$sage:sim_gtg$nage),nrow=length(sim_gtg$indyr))

for(i in 1:length(sim_gtg$indyr)){
    vbgtg[i,]<-apply(sim_gtg$VulB[sim_gtg$VulB[,1]==i,-c(1,2)],2,sum)
  }

#par(mfrow=c(1,2))
plot(apply(sim$VulB,1,sum), ylim=c(0,max(apply(sim$VulB,1,sum))),type="l")
lines(apply(vbgtg,1,sum), col="red")
# plot(apply(sim_gtg$VulB[,-1],1,sum),ylim=c(0,max(apply(sim_gtg$VulB[,-1],1,sum))))
head(sim$VulB)
head(sim$Nage)


plot(apply(sim$VulB,1,sum)-apply(vbgtg,1,sum),type="l")

#============================================================================
#Plot VBarea
#============================================================================


head(sim_gtg$tVBarea)
head(sim$VBarea)

#catch plot

summary(sim$"CatchNatAge")
dim(sim$"CatchNatAge")
head(sim$"CatchNatAge")

meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

#============================================================================
#Plot difference catches in the last 30 years
#============================================================================


CatageFRall<-data.frame(sim$"CatchNatAge"[sim$"CatchNatAge"[,1]>(max(sim$"CatchNatAge"[,1])-12*30),])
CatageFRestall<-data.frame(est$"CatchNatAge"[est$"CatchNatAge"[,1]>(max(est$"CatchNatAge"[,1])-12*30),])
CatageFRgtgall<-data.frame(sim_gtg$"CatchNatAge"[sim_gtg$"CatchNatAge"[,1]>(max(sim_gtg$"CatchNatAge"[,1])-12*30),])

CatageFRall[,1]<-rep(meses,30)[CatageFRall[,1]-(max(sim$"CatchNatAge"[,1])-12*30)]
CatageFRestall[,1]<-rep(meses,30)[CatageFRestall[,1]-(max(est$"CatchNatAge"[,1])-12*30)]
CatageFRgtgall[,1]<-rep(meses,30)[CatageFRgtgall[,1]-(max(sim_gtg$"CatchNatAge"[,1])-12*30)]

CatageFRall<-rename(CatageFRall, c("V1"="month", "V2"="area", "V3"="1", "V4"="2", "V5"= "3",
  "V6"="4", "V7"="5", "V8"= "6","V9"="7", "V10"="8", "V11"= "9","V12"="10", "V13"="11", "V14"= "12",
  "V15"="13", "V16"="14", "V17"= "15","V18"="16", "V19"="17", "V20"= "18","V21"="19", "V22"="20"))

CatageFRestall<-rename(CatageFRestall, c("V1"="month", "V2"="area", "V3"="1", "V4"="2", "V5"= "3",
  "V6"="4", "V7"="5", "V8"= "6","V9"="7", "V10"="8", "V11"= "9","V12"="10", "V13"="11", "V14"= "12",
  "V15"="13", "V16"="14", "V17"= "15","V18"="16", "V19"="17", "V20"= "18","V21"="19", "V22"="20"))

CatageFRgtgall<-rename(CatageFRgtgall, c("V1"="month", "V2"="area", "V3"="1", "V4"="2", "V5"= "3",
  "V6"="4", "V7"="5", "V8"= "6","V9"="7", "V10"="8", "V11"= "9","V12"="10", "V13"="11", "V14"= "12",
  "V15"="13", "V16"="14", "V17"= "15","V18"="16", "V19"="17", "V20"= "18","V21"="19", "V22"="20"))



# bias in gtg estimation


biasCtgtg<-cbind(CatageFRgtgall[,c(1,2)],((CatageFRestall[,-c(1,2)]-CatageFRgtgall[,-c(1,2)])/CatageFRgtgall[,-c(1,2)])*100)

biasCNplot<-melt(biasCtgtg, id=c("month","area"),variable.name="age")
biasCNplot<-arrange(transform(biasCNplot,month=factor(month,levels=meses)),month)


for(i in 1:(ncol(biasCNplot))){
  biasCNplot[is.nan(biasCNplot[,i]),i] <- NA
}


head(biasCNplot)

p <- ggplot(biasCNplot) 
p <- p + geom_boxplot(aes(x=(area), y=value,fill=age), alpha=0.5)
#p <- p + geom_line(aes(x=as.numeric(area), y=value, colour=age, alpha=0.5))
#p <- p + geom_bar(aes(x=(area), y=value,fill=age),stat="identity", alpha=0.5)
p <- p + facet_wrap(~month,ncol=4)
#p <- p + geom_vline(xintercept=48.5, linetype=3,alpha=0.3)
p <- p +  scale_y_continuous(limits=c(-.5,.5))
p <- p + theme_bw()+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title.x= element_blank())
p

#============================================================================
#Plot difference catches in the last year
#============================================================================

names(sim_gtg)
names(sim)

head(sim_gtg$"CatchNatAge")

tail(sim_gtg$"CatchNatAge")


CatageFR<-data.frame(sim$"CatchNatAge"[sim$"CatchNatAge"[,1]>(max(sim$"CatchNatAge"[,1])-12),])
CatageFRest<-data.frame(est$"CatchNatAge"[est$"CatchNatAge"[,1]>(max(est$"CatchNatAge"[,1])-12),])
CatageFRgtg<-data.frame(sim_gtg$"CatchNatAge"[sim_gtg$"CatchNatAge"[,1]>(max(sim_gtg$"CatchNatAge"[,1])-12),])




CatageFR[,1]<-meses[CatageFR[,1]-(max(sim$"CatchNatAge"[,1])-12)]
CatageFRest[,1]<-meses[CatageFRest[,1]-(max(est$"CatchNatAge"[,1])-12)]
CatageFRgtg[,1]<-meses[CatageFRgtg[,1]-(max(sim_gtg$"CatchNatAge"[,1])-12)]

CatageFR<-rename(CatageFR, c("V1"="month", "V2"="area", "V3"="1", "V4"="2", "V5"= "3",
  "V6"="4", "V7"="5", "V8"= "6","V9"="7", "V10"="8", "V11"= "9","V12"="10", "V13"="11", "V14"= "12",
  "V15"="13", "V16"="14", "V17"= "15","V18"="16", "V19"="17", "V20"= "18","V21"="19", "V22"="20"))

CatageFRest<-rename(CatageFRest, c("V1"="month", "V2"="area", "V3"="1", "V4"="2", "V5"= "3",
  "V6"="4", "V7"="5", "V8"= "6","V9"="7", "V10"="8", "V11"= "9","V12"="10", "V13"="11", "V14"= "12",
  "V15"="13", "V16"="14", "V17"= "15","V18"="16", "V19"="17", "V20"= "18","V21"="19", "V22"="20"))

CatageFRgtg<-rename(CatageFRgtg, c("V1"="month", "V2"="area", "V3"="1", "V4"="2", "V5"= "3",
  "V6"="4", "V7"="5", "V8"= "6","V9"="7", "V10"="8", "V11"= "9","V12"="10", "V13"="11", "V14"= "12",
  "V15"="13", "V16"="14", "V17"= "15","V18"="16", "V19"="17", "V20"= "18","V21"="19", "V22"="20"))

CNplot<-melt(CatageFR, id=c("month","area"),variable.name="age")
CNplot<-arrange(transform(CNplot,month=factor(month,levels=meses)),month)

CNplotest<-melt(CatageFRest, id=c("month","area"),variable.name="age")
CNplotest<-arrange(transform(CNplotest,month=factor(month,levels=meses)),month)

CNplotgtg<-melt(CatageFRgtg, id=c("month","area"),variable.name="age")
CNplotgtg<-arrange(transform(CNplotgtg,month=factor(month,levels=meses)),month)


p <- ggplot(CNplot) 
#p <- p + geom_line(aes(x=as.numeric(area), y=value, colour=age, alpha=0.5))
p <- p + geom_bar(aes(x=(area), y=value,fill=age),stat="identity", alpha=0.5)
p <- p + facet_wrap(~month,ncol=4)


#p <- p + geom_vline(xintercept=48.5, linetype=3,alpha=0.3)
p <- p + theme_bw()+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title.x= element_blank())
p


method<-c(rep("sim",nrow(CNplot)),rep("est",nrow(CNplotest)),rep("gtg",nrow(CNplotest)))

CFR<-cbind(rbind(CNplot,CNplotest,CNplotgtg),method)


for(a in 1:length(sim$sage:sim$nage)){

  ag=(sim$sage:sim$nage)[a]
  Cnplotmp<-CFR[CFR$age==ag,]

  setwd("/Users/catarinawor/Documents/Lagrangian/catageplot")
  p <- ggplot(Cnplotmp) 
  p <- p + geom_bar(aes(x=(area), y=value,fill=method),stat="identity",position="dodge", alpha=0.5)
  p <- p + facet_wrap(~month,ncol=4)
  p <- p + theme_bw()+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title.x= element_blank())
  p<-p+labs(title=paste("age",ag, sep=" "))
  p
  ggsave(file=paste("ctage",ag,".pdf", sep=""))
}


method<-c(rep("sim",nrow(CNplot)),rep("gtg",nrow(CNplotest)))

CFR<-cbind(rbind(CNplot,CNplotgtg),method)


for(a in 1:length(sim$sage:sim$nage)){

  ag=(sim$sage:sim$nage)[a]
  Cnplotmp<-CFR[CFR$age==ag,]

  setwd("/Users/catarinawor/Documents/Lagrangian/catageplot")
  p <- ggplot(Cnplotmp) 
  p <- p + geom_bar(aes(x=(area), y=value,fill=method),stat="identity",position="dodge", alpha=0.5)
  p <- p + facet_wrap(~month,ncol=4)
  p <- p + theme_bw()+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title.x= element_blank())
  p<-p+labs(title=paste("age",ag, sep=" "))
  p
  ggsave(file=paste("ctage",ag,".pdf", sep=""))
}


#============================================================================
#Compare OM- simple and gtg 
#============================================================================

names(sim)











indAyr<-rep(1:sim$fisharea,length(sim$syr:sim$nyr))






#Catch of each nation 
CN1sim<-apply(sim$"CatchNatAge"[2161:3600,],1,sum)[indAyr==1]
CN1est<-apply(est$"CatchNatAge",1,sum)[indAyr==1]
CN2sim<-apply(sim$"CatchNatAge"[2161:3600,],1,sum)[indAyr==2]
CN2est<-apply(est$"CatchNatAge",1,sum)[indAyr==2]
CN3sim<-apply(sim$"CatchNatAge"[2161:3600,],1,sum)[indAyr==3]
CN3est<-apply(est$"CatchNatAge",1,sum)[indAyr==3]


plot(rep(1:12,40),CN1sim, type="l")
lines(rep(1:12,40),CN1est,col="red")

CN1pbias<-((CN1est-CN1sim)/CN1sim)*100
CN2pbias<-((CN2est-CN2sim)/CN2sim)*100
CN3pbias<-((CN3est-CN3sim)/CN3sim)*100


par(mfrow=c(3,1))
boxplot(matrix(CN1pbias, ncol=12,byrow=T), ylim=c(-100,400))
abline(h=0, col = "red")
boxplot(matrix(CN2pbias, ncol=12,byrow=T),ylim=c(-100,400))
abline(h=0, col = "red")
boxplot(matrix(CN3pbias, ncol=12,byrow=T),ylim=c(-100,400))
abline(h=0, col = "red")

par(mfrow=c(3,1))
boxplot(matrix(CN1pbias, ncol=12,byrow=T), ylim=c(-20,20))
abline(h=0, col = "red")
boxplot(matrix(CN2pbias, ncol=12,byrow=T),ylim=c(-20,20))
abline(h=0, col = "red")
boxplot(matrix(CN3pbias, ncol=12,byrow=T),ylim=c(-20,20))
abline(h=0, col = "red")

plot(CN1pbias,CN1sim)
plot(CN1pbias,CN1est)
plot(CN1est,CN1sim,ylim=range(CN1sim),xlim=range(CN1sim))

CN1est[which.max(CN1pbias)]
CN1sim[which.max(CN1pbias)]



#==============================================================================
#==============================================================================
#migration at age - JTC talk
#==============================================================================
#==============================================================================

ages <- sim$nage:sim$sage

maxPos50<-sim$maxPos50
maxPossd<-sim$maxPossd

ys<-round(approx(ages,(1/(1+exp(-(ages-maxPos50)/maxPossd)))*(sim$narea-sim$sarea)+sim$sarea, n=1000)$y,1)
xs<-approx(ages,(1/(1+exp(-(ages-maxPos50)/maxPossd)))*(sim$narea-sim$sarea)+sim$sarea, n=1000)$x

xs[which(ys==42)]

xs[which(ys==49)]

#setwd("/Users/catarinawor/Documents/hake/JTC_talk")
pdf("movAtAge.pdf",width = 7, height = 7)
par(mfrow=c(1,1),xpd=FALSE)
plot(ages,(1/(1+exp(-(ages-maxPos50)/maxPossd)))*(sim$narea-sim$sarea)+sim$sarea, ylab="AREA",xlab="AGE",
  lwd=4, type="l", col="black", font=2, font.lab=2,cex=1.5,cex.axis=1.5,cex.lab=1.5)
segments(xs[which(ys==42)], 42, x1 =xs[which(ys==42)] , y1 = 33, col="darkred", lwd=2)
segments(xs[which(ys==42)], 42, x1 =0 , y1 = 42, col="darkred", lwd=2)
segments(xs[which(ys==49)], 49, x1 =xs[which(ys==49)] , y1 = 33, col="blue", lwd=2)
segments(xs[which(ys==49)], 49, x1 =0 , y1 = 49, col="blue", lwd=2)
text(2, 41.5,"CA",col="darkred", font=2)
text(2, 42.5,"OR",col="darkred", font=2)
text(2, 48.5,"U.S.A.", col="blue", font=2)
text(2, 49.5,"Canada", col="blue", font=2)

#dev.off()


#==============================================================================
#==============================================================================
# simulation estimation
#==============================================================================
#==============================================================================



#=======================================================================
#Simulation evaluation graphs

#parameter estimates
.SIMDIRS   <- c("/Users/catarinawor/Documents/Lagrangian/SimResult_3areas_tau1",
  "/Users/catarinawor/Documents/Lagrangian/SimResult_3areas_tau04",
  "/Users/catarinawor/Documents/Lagrangian/SimResult_5areas_tau1",
  "/Users/catarinawor/Documents/Lagrangian/SimResult_5areas_tau04")

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
    tmp_estn[j,]<-exp(sims[[3]]$est)
    tmp_pbias[j,]<-((tmp_estn[j,]-true_pars)/true_pars)*100
    tmp_maxgrad[j]<-sims[[3]]$maxgrad
    tmp_initvals[j,]<-exp(unlist(sims[[5]][1:4]))
   }

  tmp_estn<- tmp_estn[tmp_maxgrad<=1.0000e-04,]
  tmp_pbias<- tmp_pbias[tmp_maxgrad<=1.0000e-04,]
  
  estn[[i]]<-tmp_estn
  pbias[[i]]<-tmp_pbias
  maxgrad[[i]]<-tmp_maxgrad
  initvals[[i]]<-tmp_initvals[tmp_maxgrad<=1.0000e-04,]
  initvals_bad[[i]]<-tmp_initvals[tmp_maxgrad>1.0000e-04,]


}


#========================================================================
# Description of the list sims
# [[1]] -> sim- read.rep("lagrangian_OM.rep") 
# [[2]] -> est -read.rep("lagrangian_est.rep")
# [[3]] -> par - read.fit("lagrangian_est")
# [[4]] -> seed
# [[5]] -> pin - in a list
#========================================================================


titulos<-c("3 areas, tau=1.0","3 areas, tau=0.4","5 areas, tau=1.0","5 areas, tau=0.4")


setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")
#setwd("/Users/catarinawor/Documents/hake/JTC_talk")
pdf("quatroscn.pdf", width=12, height=10)
par(mfcol=c(2,2))
for( i in 1:length(.SIMDIRS)){
  boxplot(pbias[[i]],names= c(expression("t"[0]),expression("cv"),expression("a"[50]),expression("sd"["X"["max"])),ylim=c(-10,10),main=titulos[i],cex.axis=1.5,cex.lab=1.5,cex.main=1.5)
  abline(h=0)
  text(4, y = 8, labels = nrow(pbias[[i]]))
}
dev.off()




#=======================================================================

library(ggplot2)
library(reshape2)
library(animation)
library(ggmap)
library("plyr")

#======================================================================== 
# Graphs for base case scenario
#======================================================================== 
#Has the biomass stabilized - non error only

rm(list=ls()); 
#if (Sys.info()["nodename"] =="sager")  setwd("~/Dropbox/LSRA/length_SRA/sim_est_lsra")
setwd("/Users/catarinawor/Documents/Lagrangian/")
source("read.admb.R")

sim <- read.rep("lagrangian_OM.rep")
est <- read.rep("lagrangian_est.rep")

names(sim)
VBplot<-sim$VBarea[(nrow(sim$PosX)-11):nrow(sim$PosX),]

PosXplot<-sim$PosX[(nrow(sim$PosX)-11):nrow(sim$PosX),5]

x<-sim$sarea:sim$narea
propXplot<-matrix(NA, nrow=12,ncol=length(x))


for(mth in 1:12){
  propXplot[mth,]<-dnorm(x,PosXplot[mth],sim$varPos[mth])
}
Fish<-c(propXplot)


EffAplot<-sim$Effarea[(nrow(sim$PosX)-11):nrow(sim$PosX),]
Effort<-c(EffAplot)/max(c(EffAplot))*0.12
meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")


Month<-rep(meses,length(x))
Latitude<-rep(x, each=length(meses))


df<- data.frame(Month,Latitude,Fish,Effort)
df2<-arrange(transform(df,Month=factor(Month,levels=meses)),Month)

#dftry<-df[df$Month=="Jun",]


setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")
p <- ggplot(df2, aes(x=Latitude, y=-Effort)) 
p <- p + geom_bar(stat="identity", alpha=0.5, aes(fill="Effort") )
p <- p + geom_vline(xintercept=48.5, linetype=3,alpha=0.3)
p <- p + geom_line(data=df2, aes(x=Latitude, y=-Fish, colour="Abundance"))
p <- p + facet_wrap(~Month,ncol=4)
p <- p + coord_flip()
p <- p + theme_bw()+theme(legend.title=element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title.x= element_blank())
p <- p + scale_colour_grey() +scale_fill_grey() 
p
ggsave(file="model_example.pdf")

#========================================================================================
# GTG version
#========================================================================================

rm(list=ls()); 


library(plyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(animation)
library(ggmap)


#if (Sys.info()["nodename"] =="sager")  setwd("~/Dropbox/LSRA/length_SRA/sim_est_lsra")
setwd("/Users/catarinawor/Documents/Lagrangian/")
source("read.admb.R")

sim_gtg <- read.rep("lagrangian_OM_gtg.rep")
#est_gtg <- read.rep("lagrangian_est_gtg.rep")



head(sim_gtg$propVBarea)
head(sim_gtg$Nage)
head(sim_gtg$Effarea)

cores=rainbow(sim_gtg$ngroup)
plot(apply(sim_gtg$Nage[sim_gtg$Nage[,2]==1,-c(1,2)],1,sum), ylim=c(0,10000))
for(i in 1:sim_gtg$ngroup){
  lines(apply(sim_gtg$Nage[sim_gtg$Nage[,2]==i,-c(1,2)],1,sum), col=cores[i], lwd=2)
}

#propBarea indices are year, group and area
#sim_gtg$totB
#sim_gtg$Nage
#PosXplot<-matrix(NA, nrow=36,ncol=length(sim_gtg$sage:sim_gtg$nage)+1)
meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

tmpXplot<-data.frame(sim_gtg$propVBarea[sim_gtg$propVBarea[,1]>(max(sim_gtg$propVBarea[,1])-12),])
tmpXplot[,1]<-meses[tmpXplot[,1]-(max(sim_gtg$propVBarea[,1])-12)]
tmpXplot<-rename(tmpXplot, c("V1"="month", "V2"="group", "V3"= "Latitude","V4"="1", "V5"="2", "V6"= "3",
  "V7"="4", "V8"="5", "V9"= "6","V10"="7", "V11"="8", "V12"= "9","V13"="10", "V14"="11", "V15"= "12",
  "V16"="13", "V17"="14", "V18"= "15","V19"="16", "V20"="17", "V21"= "18","V22"="19", "V23"="20"))

Xplot<-melt(tmpXplot, id=c("month","group", "Latitude"),variable.name="age")

Xplot<-Xplot[Xplot$age==4,]
#Xplot<-Xplot[Xplot$group==5,]
Xplot$group<-as.factor(Xplot$group)

Xplot<-arrange(transform(Xplot,month=factor(month,levels=meses)),month)

tmpEffplot<-data.frame(month= 1:12,matrix(sim_gtg$Effarea[(nrow(sim_gtg$Effarea)-11):nrow(sim_gtg$Effarea),],ncol=31,
  dimnames=list(1:12,(sim_gtg$sarea:sim_gtg$narea))))
tmpEffplot<-setnames(tmpEffplot, old = paste("X",30:60, sep=""), new = as.character(30:60))

Effplot<-melt(tmpEffplot, id="month",variable.name="Latitude",value.name="effort")
Effplot$month<-meses[Effplot$month]
Effplot$effort<-Effplot$effort/max(Effplot$effort)*max(Xplot$value)
Effplot$Latitude<-as.numeric(Effplot$Latitude)+sim_gtg$sarea-1
df<-merge(x = Xplot, y = Effplot, by = c("Latitude","month"), all = TRUE)
head(df)

#setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")
p <- ggplot(Xplot) 
p <- p + geom_line(aes(x=as.numeric(Latitude), y=as.numeric(value), colour=group))
p <- p + geom_vline(xintercept=48.5, linetype=3,alpha=0.3)
p <- p + facet_wrap(~month,ncol=4)
#p <- p + geom_bar(data=Effplot,aes(x=as.numeric(Latitude), y=effort,fill="effort"),stat="identity", alpha=0.5 )
#p <- p + coord_flip()
#p <- p + theme_bw()+theme(panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank(),
#        axis.title.x= element_blank())
p <- p + scale_colour_grey() +scale_fill_grey(name=" ") 
p
#ggsave(file="gtg_model_example.pdf")


#groups combined


meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

tmpXplot<-data.frame(sim_gtg$propVBarea[sim_gtg$propVBarea[,1]>(max(sim_gtg$propVBarea[,1])-12),])
tmpXplot[,1]<-meses[tmpXplot[,1]-(max(sim_gtg$propVBarea[,1])-12)]
tmpXplot<-rename(tmpXplot, c("V1"="month", "V2"="group", "V3"= "Latitude","V4"="1", "V5"="2", "V6"= "3",
  "V7"="4", "V8"="5", "V9"= "6","V10"="7", "V11"="8", "V12"= "9","V13"="10", "V14"="11", "V15"= "12",
  "V16"="13", "V17"="14", "V18"= "15","V19"="16", "V20"="17", "V21"= "18","V22"="19", "V23"="20"))

Xplot<-melt(tmpXplot, id=c("month","group", "Latitude"),variable.name="age")
#Xplot<-Xplot[as.numeric(Xplot$age)>3,]

Xplot$group<-as.factor(Xplot$group)
Xplot<-arrange(transform(Xplot,month=factor(month,levels=meses)),month)


Xplot$month<-as.numeric(Xplot$month)
head(Xplot)
Xplot.v<-subset(Xplot,select= value)
Xplot.f<-subset(Xplot,select= -c(group,value))

newXplot<- aggregate(x=Xplot.v,by=Xplot.f,FUN=sum)



Effplot$effort<-Effplot$effort/max(Effplot$effort)*max(Xplot$value)
Effplot<-arrange(transform(Effplot,month=factor(month,levels=meses)),month)

#newXplot<-newXplot[,newXplot$value>0.0]
head(newXplot)
head(Effplot)

meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
newXplot$month<-meses[newXplot$month]


p1 <- ggplot(Effplot) 
p1 <- p1+   geom_bar(aes(x=Latitude, y=effort),stat="identity")
p1 <- p1 + facet_wrap(~month,ncol=4)
p1 <- p1 + geom_line(data=newXplot, aes(x=as.numeric(Latitude), y=as.numeric(value), col=age))
p1 <- p1 + geom_vline(xintercept=48.5, linetype=3,alpha=0.3)
p1 <- p1 + theme_bw()+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title.x= element_blank())
p1

#===============================================================
#===============================================================
# GTG parameter estimates scenarios
#===============================================================
#===============================================================


true_pars_gtg <- c(sim_gtg$"mo",sim_gtg$"cvPos",sim_gtg$"maxPos50",sim_gtg$"maxPossd")

indpar<-c(1,2,3,3,3,4)

#parameter estimates
.SIMDIRSGTG   <- c("/Users/catarinawor/Documents/Lagrangian/SimResult_gtg_3areas_tau1",
  "/Users/catarinawor/Documents/Lagrangian/SimResult_gtg_3areas_tau04",
  "/Users/catarinawor/Documents/Lagrangian/SimResult_gtg_5areas_tau1",
  "/Users/catarinawor/Documents/Lagrangian/SimResult_gtg_5areas_tau04")

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

  tmp_estn<- tmp_estn[tmp_maxgrad<=1.0000e-04,]
  tmp_pbias<- tmp_pbias[tmp_maxgrad<=1.0000e-04,]
  
  estn_gtg[[i]]<-tmp_estn
  pbias_gtg[[i]]<-tmp_pbias
  maxgrad_gtg[[i]]<-tmp_maxgrad
  initvals_gtg[[i]]<-tmp_initvals[tmp_maxgrad<=1.0000e-04,]
  initvals_bad_gtg[[i]]<-tmp_initvals[tmp_maxgrad>1.0000e-04,]


}

indAyr<-rep(1:3,1200)[2161:3600]
titulos<-c("3 areas, tau=1.0","3 areas, tau=0.4","5 areas, tau=1.0","5 areas, tau=0.4")


#setwd("/Users/catarinawor/Documents/hake/Thesis/figs")
#setwd("/Users/catarinawor/Documents/hake/JTC_talk")
#pdf("quatroscn.pdf")



par(mfcol=c(1,1))
for( i in 1:length(.SIMDIRSGTG)){
  boxplot(pbias_gtg[[i]],names= c(expression("t"[0]),expression("cv"),
    expression("a"[50,s]),expression("a"[50,m]),expression("a"[50,f]),expression("sd"["X"["max"]))
    ,ylim=c(-10,10),main=" ")
  abline(h=0)
  text(4, y = 8, labels = nrow(pbias_gtg[[i]]))
}









for( i in 1:length(.SIMDIRS)){
  .SIMNAME[[i]]   <- list.files(.SIMDIRS[i],pattern="\\.Rdata", full.name=TRUE)
  
  tmp_estn<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=4)
  tmp_pbias<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=4)
  tmp_maxgrad<-vector(length=length(.SIMNAME[[i]]))
  tmp_initvals<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=4)
  

  for( j in 1:length(.SIMNAME[[i]])){
    load(.SIMNAME[[i]][j])

    #parameters
    tmp_estn[j,]<-exp(sims[[3]]$est)
    tmp_pbias[j,]<-((round(tmp_estn[j,],2)-true_pars)/true_pars)*100
    tmp_maxgrad[j]<-sims[[3]]$maxgrad
    tmp_initvals[j,]<-exp(unlist(sims[[5]][1:4]))
   }

  tmp_estn<- tmp_estn[tmp_maxgrad<=1.0000e-04,]
  tmp_pbias<- tmp_pbias[tmp_maxgrad<=1.0000e-04,]
  
  estn[[i]]<-tmp_estn
  pbias[[i]]<-tmp_pbias
  maxgrad[[i]]<-tmp_maxgrad
  initvals[[i]]<-tmp_initvals[tmp_maxgrad<=1.0000e-04,]
  initvals_bad[[i]]<-tmp_initvals[tmp_maxgrad>1.0000e-04,]


}







#End of gtg plotting
#===============================================================
#Weight at age plotting


WA<-c(0.13,0.41,0.47,0.48,0.54,0.57,0.62,0.66,0.72,0.70,1.16,1.02,0.95,0.97,1.06,1.06,1.06,1.06,1.06,1.06)
plot(WA)
lines(WA+0.1)
lines(WA-0.1, col="blue")


#VBplot<-sim$VBarea[(nrow(sim$PosX)-11):nrow(sim$PosX),]
#for(g in 1:sim_gtg$ngroup){
#  PosXplot[1:12+(12*(g-1))]<-(sim_gtg$PosX[sim_gtg$PosX[,1]==g,5])[sim_gtg$indyr==max(sim_gtg$indyr)]
#
#}

x<-sim_gtg$sarea:sim_gtg$narea
propXplot<-matrix(NA, nrow=12*sim_gtg$ngroup,ncol=length(x)+2)

for(g in 1:sim_gtg$ngroup){
  for(mth in 1:12){
    propXplot[mth+(12*(g-1)),1]<-g
    propXplot[mth+(12*(g-1)),2]<-mth
    propXplot[mth+(12*(g-1)),-c(1,2)]<-dnorm(x,PosXplot[mth+(12*(g-1))],sim_gtg$varPos[g,mth])
  }
}

Fish<-c(propXplot[,-c(1,2)])
Group<-as.factor(rep(c(propXplot[,1]),length(x)))
Month<-meses[rep(c(propXplot[,2]),length(x))]
Latitude<-rep(x, each=length(meses)*sim_gtg$ngroup)


EffAplot<-sim_gtg$Effarea[(nrow(sim_gtg$Effarea)-11):nrow(sim_gtg$Effarea),]


Effort<-c(rbind(EffAplot,EffAplot,EffAplot))
meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")



df<- data.frame(Month,Group,Latitude,Fish,Effort)
df2<-arrange(transform(df,Month=factor(Month,levels=meses)),Month)


p <- ggplot(df2) 
p <- p + geom_line(data=df2, aes(x=Latitude, y=-Fish, colour=Group))
p <- p + geom_vline(xintercept=48.5, linetype=3,alpha=0.3)
p <- p + facet_wrap(~Month,ncol=4)
p <- p + geom_bar(aes(x=Latitude, y=-Effort,fill="Effort"),stat="identity", alpha=0.5 )
p <- p + coord_flip()
p <- p + theme_bw()+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title.x= element_blank())
p <- p + scale_colour_grey() +scale_fill_grey() 
p


p <- ggplot(df2, aes(x=Latitude, y=-Effort)) 
p <- p + geom_bar(stat="identity", alpha=0.5, aes(fill="Effort") )
p <- p + geom_vline(xintercept=48.5, linetype=3,alpha=0.3)
p <- p + geom_line(data=df2, aes(x=Latitude, y=-Fish, colour="Abundance"))
p <- p + facet_wrap(~Month,ncol=4)
p <- p + coord_flip()
p <- p + theme_bw()+theme(legend.title=element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title.x= element_blank())
p <- p + scale_colour_grey() +scale_fill_grey() 
p




p <- p + geom_bar(stat="identity", alpha=0.5, aes(fill="Effort") )
p <- p + geom_vline(xintercept=48.5, linetype=3,alpha=0.3)
p
p <- ggplot(df2) 
p <- p + geom_line(data=df2, aes(x=Latitude, y=-Fish, colour=Group))
p <- p + facet_wrap(~Month,ncol=4)
p <- p + coord_flip()
p <- p + theme_bw()+theme(legend.title=element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title.x= element_blank())
p <- p + scale_colour_grey() +scale_fill_grey() 
p


#========================================================================================
# Old version
#========================================================================================


meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
cores<-gray.colors(length(agep)+1)
setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")
pdf("fish_mov.pdf", width=8, height=6)
par(mfrow=c(3,4), oma = c(5,4,3,0) + 0.1, mar = c(2,2,3,1) + 0.1 )

for(mth in 1:12){
  for( i in 1:(length(PosXplot[mth,]))){
    plot(x,dnorm(x,PosXplot[mth,i],sim$varPos[i]),type="l", lwd=2, col=cores[i],main=meses[mth],xlab="",ylim=c(0,.2), ylab=" ", cex.main=2,cex.lab=2)
      abline(v=48.9)
    if(mth==1){
      legend("topright", legend=agep,  col = cores, border = "n", lwd=2, bty="n")
    }
    ##if(m==5){
    ##  polygon(c(x[30:100],x[100]), c(dnorm(x[30:100],PosX[mth,4],varPos[4]),0),col="blue")
    ##}
    par(new=T)
  }
par(new=F)
}
title(xlab = expression("Latitude "(degree)),
       outer = TRUE, line = 3,cex.main=3,cex.lab=3,font.main=2,)
 dev.off()

#setwd("/Users/catarinawor/Documents/hake/Proposal/Proposal_rev_mtng")
#pdf("fish_mov.pdf", width=6, height=4)






#========================================================================
#maps and animation plots
#========================================================================
library(ggplot2)
library(reshape2)


ntsp=1:((sim$nyr-sim$syr+1)*(sim$nmon-est$smon+1))
ages=sim$nage-sim$sage+1
indmonth= rep(sim$smon:sim$nmon,(sim$nyr-sim$syr+1))


EffNatAgeSim<-matrix(sim$EffNatAge,ncol=(sim$nage-sim$sage+3),dimnames=list(1:(((sim$nyr-sim$syr+1)*(sim$nmon-est$smon+1))*sim$nations),c("tstp","nation",sim$sage:sim$nage)))
EffNatAgeEst<-matrix(est$EffNatAge,ncol=(est$nage-est$sage+3),dimnames=list(1:(((est$nyr-sim$syr+1)*(est$nmon-est$smon+1))*est$nations),c("tstp","nation",est$sage:est$nage)))

EffNatAgeSim<-as.data.frame(EffNatAgeSim)
 
 ENAsim<-NULL

for(i in 1:sim$nations)
{
	tmp <- melt(EffNatAgeSim[EffNatAgeSim$nation==i,sim$sage:sim$nage+2])
	tmp2 <- cbind(tim=rep(ntsp,ages),tmp, nat=rep(i,nrow(tmp)),met=rep("simulated",nrow(tmp)))
	ENAsim <- rbind(ENAsim,tmp2)
}

names(ENAsim) <- c("time","age", "effort","nations","method") 


#============####============##
EffNatAgeEst<-matrix(est$EffNatAge,ncol=(est$nage-est$sage+3),dimnames=list(1:(((est$nyr-sim$syr+1)*(est$nmon-est$smon+1))*est$nations),c("tstp","nation",est$sage:est$nage)))
EffNatAgeEst<-as.data.frame(EffNatAgeEst)
 
 ENAest<-NULL

for(i in 1:est$nations)
{
	tmp <- melt(EffNatAgeEst[EffNatAgeEst$nation==i,est$sage:est$nage+2])
	tmp2 <- cbind(tim=rep(ntsp,ages),tmp, nat=rep(i,nrow(tmp)),met=rep("estimated",nrow(tmp)))
	ENAest <- rbind(ENAest,tmp2)
}

names(ENAest) <- c("time","age", "effort","nations","method") 

ENA<- rbind(ENAsim,ENAest)

ENA <- ENA[order(ENA$time),]

maxeff<-max(ENA$effort)
meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
mets<-c("simulated","estimated")

setwd("/Users/catarinawor/Documents/Lagrangian/anime/effort")

library("animation")

ani.options(ani.dev = "pdf", ani.type = "pdf",ani.width=8, ani.height=4)
saveLatex(
#if you don't have latex installed you need to install it or use other function such as saveGIF

for(i in 1:length(ntsp))
{
  i=6
  df <-ENA[ENA$time==i,]
  for(n in 1:est$nations )
  {
    for(m in 1:(length(mets)) )
    {
      if(min(df[df$nations==n&df$method==mets[m],]$effort)>0)
      {
        tmp3 <- df[df$nations==n&df$method==mets[m],]$effort/max(df[df$nations==n&df$method==mets[m],]$effort)
        df[df$nations==n&df$method==mets[m],]$effort = tmp3
      }
    }
  }
    
  y<- ggplot(df, aes(x=age,y=effort, group=method)) + theme_bw()
  y <- y+ facet_grid(.~nations)
  y <- y+ geom_line(size=2, aes(colour=method))
  y <- y + scale_y_continuous(limits=c(0,1))
  y <- y + labs(title=meses[indmonth[i]])
  

print(y)
}
,
pdflatex = "/usr/texbin/pdflatex",latex.filename = "selnation.tex")


sumNatEff<-apply(EffNatAgeEst[,-(1:2)],1,sum)
mon<- rep(1:12,est$nyr)
yr<- rep(1:est$nyr,each=12)
totNatEff<- cbind(EffNatAgeEst[,(1:2)],mon,yr,sumNatEff)
head(totNatEff)




p <- ggplot(totNatEff, aes(x=mon,y=sumNatEff,fill=as.factor(nation))) + theme_bw()
p <- p + geom_bar(stat="identity", position="dodge")
p <- p+facet_wrap(~yr,ncol=4)
p

EffNatAgeEst<-as.data.frame(EffNatAgeEst)
 


for(i in 1:length(ntsp))
{
  df <-ENA[ENA$time==i,]
  for(n in 1:est$nations )
  {
    for(m in 1:(length(mets)) )
    {
      if(min(df[df$nations==n&df$method==mets[m],]$effort)>0)
      {
        tmp3 <- df[df$nations==n&df$method==mets[m],]$effort/max(df[df$nations==n&df$method==mets[m],]$effort)
        df[df$nations==n&df$method==mets[m],]$effort = tmp3
      }
    }
  }
    
  y<- ggplot(df, aes(x=age,y=effort, group=method)) + theme_bw()
  y <- y + facet_grid(.~nations)
  y <- y + geom_line(size=2, aes(colour=method))
  y <- y + scale_y_continuous(limits=c(0,1))
  y <- y + labs(title=meses[indmonth[i]])
  

print(y)
}



#====================================================================================================================
setwd("/Users/catarinawor/Documents/Lagrangian/anime/VBarea")

VBareaSim<-matrix(sim$VBarea, ncol=(sim$narea-sim$sarea+1),dimnames=list(ntsp,sim$sarea:sim$narea))
VBareaEst<-matrix(est$VBarea, ncol=(est$narea-est$sarea+1),dimnames=list(ntsp,est$sarea:est$narea))

VBplotSim<-cbind(melt(VBareaSim),rep("simulated",nrow(melt(VBareaSim))))
VBplotEst<-cbind(melt(VBareaEst),rep("estimated",nrow(melt(VBareaEst))))


names(VBplotSim)<- c("time","area", "VB","method")
names(VBplotEst)<- c("time","area", "VB","method")

lat<-VBareaplot$area
lon<-rep(-131,length(lat))
VBareaplot<-rbind(VBplotSim,VBplotEst)

lat<-VBareaplot$area
lon<-rep(-131,length(lat))

VBareaplot<-cbind(VBareaplot,lat,lon)

library(ggmap)


basemap<-get_map(location = c(lon = -125, lat = 45),
    zoom = 5, maptype = "terrain")

saveLatex(
#if you don't have latex installed you need to install it or use other function such as saveGIF


for(i in 1:ntsp)
{
ex1<-VBareaplot[VBareaplot$time==i ,]
		minVB<-min(VBareaplot$VB)
		maxVB<-max(VBareaplot$VB)

p2<-  ggmap(basemap,
    extent = "panel",
    ylab = "Latitude",
    xlab = "Longitude")
p2 <- p2 + geom_line(y=48.5, linetype=2, colour="grey60")
p2 <- p2 + geom_point(alpha=0.8,aes(size=VB, group=method,colour=method),data=ex1) 
p2 <- p2 + labs(title=meses[indmonth[i]])
p2 <- p2 + scale_size_area(limits=c(minVB,maxVB),max_size = 10,breaks = c(50,100,200,400,600), labels = c(50,100,200,400,600), name = "vulnerable biomass")
print(p2)   


}
,
pdflatex = "/usr/texbin/pdflatex")


#==================================================================================
#Random graphs
#==================================================================================
 #=====================================================================
#just one simulation scenario

#with Tau =40
.SIMDIRS   <- "/Users/catarinawor/Documents/Lagrangian/SimResult_tau40"
.SIMNAME   <- list.files(.SIMDIRS,pattern="\\.Rdata", full.name=TRUE)

estn<-matrix(NA,nrow=length(.SIMNAME),ncol=4)
pbias<-matrix(NA,nrow=length(.SIMNAME),ncol=4)



for( i in 1:length(.SIMNAME )){
  
    load(.SIMNAME[i])

    #parameters
    estn[i,]<-exp(sims[[3]]$est)
    pbias[i,]<-((estn[i,]-true_pars)/true_pars)*100

}


boxplot(pbias,names= c("mo","cvPos","maxPos50","maxPossd"),ylim=c(-100,100),main="5 areas, tau=40")
  abline(h=0)




#=====================================================================
#


CN1pb_median<-matrix(NA,nrow=length(.SIMNAME),ncol=12)
CN1pb_sd<-matrix(NA,nrow=length(.SIMNAME),ncol=12)
CN2pb_median<-matrix(NA,nrow=length(.SIMNAME),ncol=12)
CN2pb_sd<-matrix(NA,nrow=length(.SIMNAME),ncol=12)
CN3pb_median<-matrix(NA,nrow=length(.SIMNAME),ncol=12)
CN3pb_sd<-matrix(NA,nrow=length(.SIMNAME),ncol=12)






for( i in 1:length(.SIMNAME)){
  load(.SIMNAME[i])

  #parameters
  estn[i,]<-exp(sims[[3]]$est)
  pbias[i,]<-((estn[i,]-true_pars)/true_pars)*100


  #catch
  #Catch of each nation 
  CN1sim<-apply(sims[[1]]$"CatchNatAge"[2161:3600,],1,sum)[indAyr==1]
  CN1est<-apply(sims[[2]]$"CatchNatAge",1,sum)[indAyr==1]
  CN2sim<-apply(sims[[1]]$"CatchNatAge"[2161:3600,],1,sum)[indAyr==2]
  CN2est<-apply(sims[[2]]$"CatchNatAge",1,sum)[indAyr==2]
  CN3sim<-apply(sims[[1]]$"CatchNatAge"[2161:3600,],1,sum)[indAyr==3]
  CN3est<-apply(sims[[2]]$"CatchNatAge",1,sum)[indAyr==3]

  CN1pbias<-((CN1est-CN1sim)/CN1sim)*100
  CN2pbias<-((CN2est-CN2sim)/CN2sim)*100
  CN3pbias<-((CN3est-CN3sim)/CN3sim)*100

  CN1pb_median[i,]<-apply(matrix(CN1pbias, ncol=12,byrow=T),2,median)
  CN2pb_median[i,]<-apply(matrix(CN2pbias, ncol=12,byrow=T),2,median)
  CN3pb_median[i,]<-apply(matrix(CN3pbias, ncol=12,byrow=T),2,median)

  CN1pb_sd[i,]<-sqrt(apply(matrix(CN1pbias, ncol=12,byrow=T),2,var))
  CN2pb_sd[i,]<-sqrt(apply(matrix(CN2pbias, ncol=12,byrow=T),2,var))
  CN3pb_sd[i,]<-sqrt(apply(matrix(CN3pbias, ncol=12,byrow=T),2,var))


}

#parameter estimate plot
#setwd("/Users/catarinawor/Documents/hake/Thesis/figs")
setwd("/Users/catarinawor/Documents/hake/JTC_talk")
pdf("3terr_scn_tau04.pdf")
boxplot(pbias,names= c("mo","cvPos","maxPos50","maxPossd"),ylim=c(-10,10))
abline(h=0)
dev.off()


#catch plots
par(mfcol=c(3,2))
boxplot(CN1pb_median)
abline(h=0, col="red")
boxplot(CN2pb_median)
abline(h=0, col="red")
boxplot(CN3pb_median)
abline(h=0, col="red")

boxplot(CN1pb_sd,ylim=c(0,1000))
boxplot(CN2pb_sd,ylim=c(0,1000))
boxplot(CN3pb_sd,ylim=c(0,1000))

#========================================================================
# Description of the list sims
# [[1]] -> sim- read.rep("lagrangian_OM.rep") 
# [[2]] -> est -read.rep("lagrangian_est.rep")
# [[3]] -> par - read.fit("lagrangian_est")
# [[4]] -> seed
#========================================================================


ls() 
length(sims)
names(sims[[2]])

sims[[1]]$"nyr"
dim(sims[[1]]$"CatchNatAge")


sims[[1]]$"CatchNatAge"[,1]

ind<- rep(1:3,1200)
par(mfrow=c(3,1))
plot(rep(1:12,100),apply(sims[[1]]$"CatchNatAge",1,sum)[ind==1], type="l")
lines(rep(1:12,100),apply(sims[[2]]$"CatchNatAge",1,sum)[ind==1],col="red")
plot(rep(1:12,100),apply(sims[[1]]$"CatchNatAge",1,sum)[ind==2], type="l")
lines(rep(1:12,100),apply(sims[[2]]$"CatchNatAge",1,sum)[ind==2],col="red")
plot(rep(1:12,100),apply(sims[[1]]$"CatchNatAge",1,sum)[ind==3], type="l")
lines(rep(1:12,100),apply(sims[[2]]$"CatchNatAge",1,sum)[ind==3],col="red")

catvar<-(apply(sims[[2]]$"CatchNatAge",1,sum)[ind==1]-apply(sims[[1]]$"CatchNatAge",1,sum)[ind==1])/apply(sims[[1]]$"CatchNatAge",1,sum)[ind==1]

boxplot(catvar, ylim=c(-1,50))

names(sims[[2]])

plot(apply(sims[[1]]$"CatchNatAge",1,sum)[1201:2400])
plot(apply(sims[[1]]$"CatchNatAge",1,sum)[2401:3600])


plot(estn)

mydat<-data.frame(parameter=rep(c("mo","cvPos","maxPos50","maxPossd"),each=100),estimate=c(estn))

plot(mydat)
points(true_pars[c(2,3,4,1)],col="red",pch=16)

