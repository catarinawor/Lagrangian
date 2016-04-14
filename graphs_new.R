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

sim_gtg <- read.rep("lag_OM_gtg_new.rep")
sim_gtg <- read.rep("lag_OM_gtg_new.rep")


#======================================================================================
# it seems that vbarea is messed up. is Nage messed up as well?
#======================================================================================
plot(sim_gtg$prop_ng)

head(sim_gtg$Areano)

arno<-data.frame(sim_gtg$Areano[sim_gtg$Areano[,1]>(max(sim_gtg$Areano[,1])-12),])
arno1<-arno[arno[,1]==1189

Nag<-data.frame(sim_gtg$Nage[sim_gtg$Nage[,1]>(max(sim_gtg$Nage[,1])-12),])
Nag1<-Nag[Nag[,1]==1189,]

cbind(floor(arno1$V3),Nag1$V3)

ar<-unique(floor(arno1$V3))

nums<-NULL

for(i in 1:(length(ar))){
  nums[i]<-sum(Nag1$V3[floor(arno1$V3)==ar[i]])
}

plot(ar,nums, type="l")
plot(floor(arno1$V3),Nag1$V3,type="l")

head(sim_gtg$Nage)
meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")


#======================================================================================
# plotting Bareaage
#======================================================================================

names(sim_gtg)
dim(sim_gtg$BAreaAge)
head(sim_gtg$BAreaAge)

meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

tmpXplot<-data.frame(sim_gtg$BAreaAge[sim_gtg$BAreaAge[,1]>(max(sim_gtg$BAreaAge[,1])-12),])
tmpXplot[,1]<-meses[tmpXplot[,1]-(max(sim_gtg$tVBarea[,1])-12)]

tmpXplot<-rename(tmpXplot, c("V1"="month", "V2"="area", "V3"="1", "V4"="2", "V5"= "3",
  "V6"="4", "V7"="5", "V8"= "6","V9"="7", "V10"="8", "V11"= "9","V12"="10", "V13"="11", "V14"= "12",
  "V15"="13", "V16"="14", "V17"= "15","V18"="16", "V19"="17", "V20"= "18","V21"="19", "V22"="20"))



Xplot<-melt(tmpXplot, id=c("month","area"),variable.name="age")

Xplot<-Xplot[Xplot$area>30&Xplot$area<60,]
Xplot


Xplot<-arrange(transform(Xplot,month=factor(month,levels=meses)),month)

tmpEffplot<-data.frame(month= 1:12,matrix(sim_gtg$Effarea[(nrow(sim_gtg$Effarea)-11):nrow(sim_gtg$Effarea),],ncol=31,
  dimnames=list(1:12,(sim_gtg$sarea:sim_gtg$narea))))
tmpEffplot<-setnames(tmpEffplot, old = paste("X",30:60, sep=""), new = as.character(30:60))

Effplot<-melt(tmpEffplot, id="month",variable.name="area",value.name="effort")
Effplot$area<-as.numeric(Effplot$area)+29
Effplot<-Effplot[Effplot$area>30&Effplot$area<60,]

Effplot$month<-meses[Effplot$month]
Effplot$effort<-Effplot$effort/max(Effplot$effort)*max(Xplot$value)
df<-merge(x = Xplot, y = Effplot, by = c("area","month"), all = TRUE)
head(df)

Xplot2<-Xplot[Xplot$age==6,]

p <- ggplot(Xplot2) 
p <- p + geom_line(aes(x=as.numeric(area), y=value, colour=age, alpha=0.5))
#p <- p + geom_bar(aes(x=(area), y=value,fill=age),stat="identity", alpha=0.5)
p <- p + facet_wrap(~month,ncol=4)
p <- p + geom_vline(xintercept=48.5, linetype=3,alpha=0.3)
p <- p + theme_bw()+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title.x= element_blank())
p


setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")
p <- ggplot(Xplot) 
p <- p + geom_line(aes(x=as.numeric(area), y=value, colour=age, alpha=0.5))
#p <- p + geom_bar(aes(x=(area), y=value,fill=age),stat="identity", alpha=0.5)
p <- p + facet_wrap(~month,ncol=4)
p <- p + geom_vline(xintercept=48.5, linetype=3,alpha=0.3)
p <- p + theme_bw()+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title.x= element_blank())
p
ggsave(file="gtg_bars_expl_graph.pdf")

setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")

p <- ggplot(Xplot) 
p <- p + geom_line(aes(x=as.numeric(area), y=value, colour=age, alpha=0.5))
#p <- p + geom_bar(aes(x=(area), y=value,fill=age),stat="identity", alpha=0.5)
p <- p + facet_wrap(~month,ncol=4)
p <- p + geom_vline(xintercept=48.5, linetype=3,alpha=0.3)
p <- p + geom_bar(data=Effplot,aes(x=area, y=effort),colour="black",stat="identity" )
p <- p + theme_bw()+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title.x= element_blank())
p

ggsave(file="gtg_bars_expl_graph.pdf")


p <- ggplot(Xplot) 
p <- p + geom_bar(data=Effplot,aes(x=area, y=effort),colour="black",stat="identity" )
p <- p + geom_line(aes(x=as.numeric(area), y=value, colour=age),alpha=0.5, size=1.5)
#p <- p + geom_bar(aes(x=(area), y=value,fill=age),stat="identity", alpha=0.5)
p <- p + facet_wrap(~month,ncol=4)
p <- p + geom_vline(xintercept=48.5, linetype=3,alpha=0.3)
p <- p + theme_bw()+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title.x= element_blank())
p
ggsave(file="gtg_lines_expl_graph.pdf")







#==========================================================================
#Simulation evaluation
#
#==========================================================================
#"/Users/catarinawor/Documents/Lagrangian/SimResult_3areas_tau1",
#  "/Users/catarinawor/Documents/Lagrangian/SimResult_3areas_tau04",
#  "/Users/catarinawor/Documents/Lagrangian/SimResult_5areas_tau1",

setwd("/Users/catarinawor/Documents/Lagrangian/")
source("read.admb.R")

sim_gtg <- read.rep("lag_OM_gtg_new.rep")

true_pars <- c(sim_gtg$"mo",sim_gtg$"cvPos",sim_gtg$"maxPos50",sim_gtg$"maxPossd")  

.SIMDIRS   <- c("/Users/catarinawor/Documents/Lagrangian/SimResult_gtg_5areas_tau1",
  "/Users/catarinawor/Documents/Lagrangian/SimResult_gtg_5areas_tau04")

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
#pdf("quatroscn.pdf", width=12, height=10)
par(mfcol=c(1,1))
for( i in 1:length(.SIMDIRS)){
  
  boxplot(pbias[[i]],names= c(expression("t"[0]),expression("cv"),expression("a"[50]),
    expression("sd"["X"["max"]])),ylim=c(-50,50),main=titulos[i],cex.axis=1.5,cex.lab=1.5,
  cex.main=1.5)
  abline(h=0)
  text(4, y = 8, labels = nrow(pbias[[i]]))

}
#dev.off()

#===========================================================================
#catch plot


#catch plot

summary(sims[[1]]$"CatchNatAge")
dim(sim$"CatchNatAge")
head(sim$"CatchNatAge")
meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

CatageFRall<-data.frame(sim$"CatchNatAge"[sim$"CatchNatAge"[,1]>(max(sim$"CatchNatAge"[,1])-12*30),])
CatageFRest<-data.frame(est$"CatchNatAge"[est$"CatchNatAge"[,1]>(max(est$"CatchNatAge"[,1])-12*30),])


CatageFRgtg<-data.frame(sim_gtg$"CatchNatAge"[sim_gtg$"CatchNatAge"[,1]>(max(sim_gtg$"CatchNatAge"[,1])-12*30),])





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

CatageFRgtg<-rename(CatageFRest, c("V1"="month", "V2"="area", "V3"="1", "V4"="2", "V5"= "3",
  "V6"="4", "V7"="5", "V8"= "6","V9"="7", "V10"="8", "V11"= "9","V12"="10", "V13"="11", "V14"= "12",
  "V15"="13", "V16"="14", "V17"= "15","V18"="16", "V19"="17", "V20"= "18","V21"="19", "V22"="20"))


CNplot<-melt(CatageFR, id=c("month","area"),variable.name="age")
CNplot<-arrange(transform(CNplot,month=factor(month,levels=meses)),month)

p <- ggplot(CNplot) 
#p <- p + geom_line(aes(x=as.numeric(area), y=value, colour=age, alpha=0.5))
p <- p + geom_bar(aes(x=(area), y=value,fill=age),stat="identity", alpha=0.5)
p <- p + facet_wrap(~month,ncol=4)
#p <- p + geom_vline(xintercept=48.5, linetype=3,alpha=0.3)

p <- p + theme_bw()+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title.x= element_blank())
p

CNplotest<-melt(CatageFRest, id=c("month","area"),variable.name="age")
CNplotest<-arrange(transform(CNplotest,month=factor(month,levels=meses)),month)

CNplotgtg<-melt(CatageFRgtg, id=c("month","area"),variable.name="age")
CNplotgtg<-arrange(transform(CNplotgtg,month=factor(month,levels=meses)),month)


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

