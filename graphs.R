#=========================================================================
#graphing simulation -estiation routine for lagrangian model
# Author: Catarina Wor
# Jun 1st 2015
#=========================================================================

rm(list=ls()); 
#if (Sys.info()["nodename"] =="sager")  setwd("~/Dropbox/LSRA/length_SRA/sim_est_lsra")
setwd("/Users/catarinawor/Documents/Lagrangian/")
source("read.admb.R")

sim = read.rep("lagrangian_OM.rep")
est = read.rep("lagrangian_est.rep")

nomes <- names(sim)

true_pars <- c(sim$"mo", exp(sim$"log_tau_c"),sim$"maxPos50",sim$"maxPossd",sim$"cvPos")  
est_pars <- c(est$"mo",exp(est$"log_tau_c"),est$"maxPos50",est$"maxPossd",est$"cvPos")

#parameter plot
par(mfrow=c(1,1))
barplot(t(matrix(c(true_pars,est_pars),ncol=2)),names.arg = c("mo","tau_c","maxPos50","maxPossd","cvPos"),beside=T)



#========================================================================
#multiplot function from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#======================================================================================================

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

