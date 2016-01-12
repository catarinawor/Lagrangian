#=========================================================================
# graphing lagrangian model structure and sim-est
# Author: Catarina Wor
# Jun 1st 2015 - updated in November 26th 2015
#=========================================================================

rm(list=ls()); 
#if (Sys.info()["nodename"] =="sager")  setwd("~/Dropbox/LSRA/length_SRA/sim_est_lsra")
setwd("/Users/catarinawor/Documents/Lagrangian/")
source("read.admb.R")

sim <- read.rep("lagrangian_OM.rep")
est <- read.rep("lagrangian_est.rep")

nomes <- names(sim)

true_pars <- c(sim$"mo",sim$"cvPos",sim$"maxPos50",sim$"maxPossd")  
est_pars  <- c(est$"mo",est$"cvPos",est$"maxPos50",est$"maxPossd")

#parameter plot
par(mfrow=c(1,1))
barplot(t(matrix(c(true_pars,est_pars),ncol=2)),names.arg = c("mo","cvPos","maxPos50","maxPossd"),beside=T)

#catch plot
indAyr<-rep(1:3,1200)[2161:3600]

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

#migration at age



ages <- sim$nage:sim$sage

maxPos50<-sim$maxPos50
maxPossd<-sim$maxPossd

ys<-round(approx(ages,(1/(1+exp(-(ages-maxPos50)/maxPossd)))*(sim$narea-sim$sarea)+sim$sarea, n=1000)$y,1)
xs<-approx(ages,(1/(1+exp(-(ages-maxPos50)/maxPossd)))*(sim$narea-sim$sarea)+sim$sarea, n=1000)$x

xs[which(ys==42)]

xs[which(ys==49)]


setwd("/Users/catarinawor/Documents/hake/JTC_talk")

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

dev.off()



#=======================================================================
#Simulation evaluation graphs

#parameter estimates
.SIMDIRS   <- "/Users/catarinawor/Documents/Lagrangian/SimResult"
.SIMNAME   <- list.files(.SIMDIRS,pattern="\\.Rdata", full.name=TRUE)

estn<-matrix(NA,nrow=length(.SIMNAME),ncol=4)
pbias<-matrix(NA,nrow=length(.SIMNAME),ncol=4)

#========================================================================
# Description of the list sims
# [[1]] -> sim- read.rep("lagrangian_OM.rep") 
# [[2]] -> est -read.rep("lagrangian_est.rep")
# [[3]] -> par - read.fit("lagrangian_est")
# [[4]] -> seed
#========================================================================

indAyr<-rep(1:3,1200)[2161:3600]

for( i in 1:length(.SIMNAME)){
  load(.SIMNAME[i])

  #parameters
  estn[i,]<-exp(sims[[3]]$est)
  pbias[i,]<-((estn[i,]-true_pars)/true_pars)*100

  

}

setwd("/Users/catarinawor/Documents/hake/Thesis/figs")
pdf("3terr_scn_tau0.2.pdf")
boxplot(pbias,names= c("mo","cvPos","maxPos50","maxPossd"),)
abline(h=0)
dev.off()

#with Tau =0.4
.SIMDIRS   <- "/Users/catarinawor/Documents/Lagrangian/SimResult_tau04"
.SIMNAME   <- list.files(.SIMDIRS,pattern="\\.Rdata", full.name=TRUE)

estn<-matrix(NA,nrow=length(.SIMNAME),ncol=4)
pbias<-matrix(NA,nrow=length(.SIMNAME),ncol=4)

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

#=======================================================================

library(ggplot2)
library(reshape2)
library(animation)
library(ggmap)

#======================================================================== 
# Graphs for base case scenario
#======================================================================== 
#Has the biomass stabilized -non error only




PosXplot<-sim$PosX[(nrow(sim$PosX)-11):nrow(sim$PosX),c(1,3,5,7,9,11,13,15,17,19)]

x<-baseF$sarea:baseF$narea
agep<-(baseF$sage:baseF$nage)[c(1,3,5,7,9,11,13,15,17,19)]


meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
cores<-rainbow(length(agep))
#setwd("/Users/catarinawor/Documents/hake/Proposal/Proposal_rev_mtng")
#pdf("fish_mov.pdf", width=6, height=4)
par(mfrow=c(3,4), mar=c(4.1,4.1,4.1,0.3),c(4, 1, 1, 1) )
for(mth in 1:12){
  for( i in 1:(length(PosXplot[mth,]))){
    plot(x,dnorm(x,PosXplot[mth,i],baseF$varPos[i]),type="l", lwd=2, col=cores[i],main=meses[mth],xlab="",ylim=c(0,.2), ylab=" ", cex.main=3,cex.lab=2)
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
mtext(expression("Latitude "(degree)), side = 1, line = -1, outer=T, cex=1.5, font=2)



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

