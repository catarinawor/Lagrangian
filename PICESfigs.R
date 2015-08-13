#=========================================================================
# Routine for generating figures for PICES presentations
# Author: Catarina Wor
# Aug 5th 2015
#=========================================================================

rm(list=ls()); 

setwd("/Users/catarinawor/Documents/Lagrangian/")
source("read.admb.R")

#read data from lowF scenario
setwd("/Users/catarinawor/Documents/Lagrangian/lowF")
loF<-read.rep("lagrangian_OM.rep")

#read data from highF scenario
setwd("/Users/catarinawor/Documents/Lagrangian/highF")
hiF<-read.rep("lagrangian_OM.rep")

#========================================================================
library(ggplot2)
library(reshape2)
library(animation)
library(ggmap)


ntsp <- 1:((loF$nyr-loF$syr+1)*(loF$nmon-loF$smon+1))
nage <- loF$nage-loF$sage+1
ages <- loF$nage:loF$sage
indmonth <- rep(loF$smon:loF$nmon,(loF$nyr-loF$syr+1))
indyr <-  rep(loF$syr:loF$nyr,each=(loF$nmon-loF$smon+1))

maxPos501<-loF$maxPos501
maxPos502<-loF$maxPos502
maxPossd1<-loF$maxPossd1
maxPossd2<-loF$maxPossd2

plot(ages,(1/(1+exp(-(ages-maxPos501)/maxPossd1)))*(loF$narea-loF$sarea)+loF$sarea, ylab="area", lwd=3, type="l", col="deeppink2")
lines(ages,1/(1+exp(-(ages-maxPos502)/maxPossd2))*(loF$narea-loF$sarea)+loF$sarea, col="olivedrab4",lwd=3)
legend("topleft",c("high abundance", "low abundance"), bty="n", col=c("deeppink2","olivedrab4"),
  lwd=3,pch=NULL,text.col=c("deeppink2","olivedrab4"), text.font=2)



par(mfrow=c(1,1))
plot(ntsp,loF$SB, type="l", ylim=c(0,300))
lines(ntsp,hiF$SB,  col="blue")



VBarealoF<-matrix(loF$VBarea, ncol=(loF$narea-loF$sarea+1),dimnames=list(ntsp,loF$sarea:loF$narea))
VBareahiF<-matrix(hiF$VBarea, ncol=(hiF$narea-hiF$sarea+1),dimnames=list(ntsp,hiF$sarea:hiF$narea))

VBplotloF<-cbind(melt(VBarealoF),rep("low F",nrow(melt(VBarealoF))))
VBplothiF<-cbind(melt(VBareahiF),rep("high F",nrow(melt(VBareahiF))))

names(VBplotloF)<- c("time","area", "VB","Scenario")
names(VBplothiF)<- c("time","area", "VB","Scenario")

VBareaplot<-rbind(VBplotloF,VBplothiF)
lat<-VBareaplot$area
lon<-rep(-131,length(lat))
nation<-rep(1,nrow(VBareaplot))
nation[VBareaplot$lat>48.1]<-2

VBareaplot<-cbind(VBareaplot,lat,lon,month=indmonth,yr=indyr)

summary(VBareaplot)

#this figure is probably not useful for this scenario

meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

basemap<-get_map(location = c(lon = -125, lat = 45),
    zoom = 5, maptype = "terrain")
setwd("/Users/catarinawor/Documents/hake/PICES_conference/presentation/Fscn_anime")

##saveLatex(
###if you don't have latex installed you need to install it or use other function such as saveGIF
##  for(i in 121:max(ntsp))
##  {
##      ex1<-VBareaplot[VBareaplot$time==i ,]
##      minVB<-min(VBareaplot$VB)
##      maxVB<-max(VBareaplot$VB)
##      
##      p2<-  ggmap(basemap,
##          extent = "panel",
##          ylab = "Latitude",
##          xlab = "Longitude")
##      p2 <- p2 + geom_line(y=48.5, linetype=2, colour="grey60")
##      p2 <- p2 + geom_point(alpha=0.8,aes(size=VB, group=Scenario,colour=Scenario),data=ex1) 
##      p2 <- p2 + facet_wrap(~Scenario)
##      p2 <- p2 + labs(title=paste(meses[indmonth[i]],", year: ",yr[indmonth[i]]))
##      p2 <- p2 + scale_size_area(limits=c(minVB,maxVB),max_size = 10,breaks = c(50,100,200,400,600), labels = c(50,100,200,400,600), name = "vulnerable biomass")
##      print(p2)   
##  }
##  ,
##  pdflatex = "/usr/texbin/pdflatex")



meanVBloFnat1<-NULL
meanVBloFnat2<-NULL

meanVBhiFnat1<-NULL
meanVBhiFnat2<-NULL

for(i in 1:12){

  meanVBloFnat1[i]<-mean(VBareaplot$VB[VBareaplot$Scenario=="low F"&VBareaplot$area<48.1&VBareaplot$month==i&VBareaplot$yr>10])
  meanVBloFnat2[i]<-mean(VBareaplot$VB[VBareaplot$Scenario=="low F"&VBareaplot$area>48.1&VBareaplot$month==i&VBareaplot$yr>10])
  meanVBhiFnat1[i]<-mean(VBareaplot$VB[VBareaplot$Scenario=="high F"&VBareaplot$area<48.1&VBareaplot$month==i&VBareaplot$yr>10])
  meanVBhiFnat2[i]<-mean(VBareaplot$VB[VBareaplot$Scenario=="high F"&VBareaplot$area>48.1&VBareaplot$month==i&VBareaplot$yr>10])

}

meanVBloFnat1p<-meanVBloFnat1/(meanVBloFnat1+meanVBloFnat2)
meanVBloFnat2p<-meanVBloFnat2/(meanVBloFnat1+meanVBloFnat2)
meanVBhiFnat1p<-meanVBhiFnat1/(meanVBhiFnat1+meanVBhiFnat2)
meanVBhiFnat2p<-meanVBhiFnat2/(meanVBhiFnat1+meanVBhiFnat2)


VulB<-c(meanVBloFnat1,meanVBloFnat2,meanVBhiFnat1,meanVBhiFnat2,meanVBloFnat1p,meanVBloFnat2p,meanVBhiFnat1p,meanVBhiFnat2p)
month<-rep(1:12,8)
nation<-rep(c(rep(1,12),rep(2,12),rep(1,12),rep(2,12)),2)
scenario<-rep(c(rep("low F",24),rep("high F",24)),2)
measure<-c(rep("nominal",48),rep("proportion",48))


df<- data.frame(VulB,month,nation,scenario,measure)

p <- ggplot(df, aes(x=as.factor(month), y=VulB, fill=as.factor(nation)))
p <- p + geom_bar(stat = "identity")
p <- p + facet_grid(measure~scenario,scales = "free_y")
p <- p + labs(x="month", y="vulnerable biomass",fill="nation")
p <- p + theme_bw()
p 



names(loF)

VBarealoF<-matrix(loF$VBarea, ncol=(loF$narea-loF$sarea+1),dimnames=list(ntsp,loF$sarea:loF$narea))
VBareahiF<-matrix(hiF$VBarea, ncol=(hiF$narea-hiF$sarea+1),dimnames=list(ntsp,hiF$sarea:hiF$narea))

EffarealoF<-matrix(loF$Effarea, ncol=(loF$narea-loF$sarea+1),dimnames=list(ntsp,loF$sarea:loF$narea))
EffareahiF<-matrix(hiF$Effarea, ncol=(hiF$narea-hiF$sarea+1),dimnames=list(ntsp,hiF$sarea:hiF$narea))

EffplotloF<-cbind(melt(EffarealoF),rep("low F",nrow(melt(EffarealoF))))
EffplothiF<-cbind(melt(EffareahiF),rep("high F",nrow(melt(EffareahiF))))

VBplotloF<-cbind(melt(VBarealoF),rep("low F",nrow(melt(VBarealoF))))
VBplothiF<-cbind(melt(VBareahiF),rep("high F",nrow(melt(VBareahiF))))

names(VBplotloF)<- c("time","area", "value","Scenario")
names(VBplothiF)<- c("time","area", "value","Scenario")

names(EffplotloF)<- c("time","area", "value","Scenario")
names(EffplothiF)<- c("time","area", "value","Scenario")

VBEffareaplot<-rbind(VBplotloF,VBplothiF,EffplotloF,EffplothiF)


lat<-VBEffareaplot$area
lon<-rep(-131,length(lat))
nation<-rep(1,nrow(VBEffareaplot))
nation[VBEffareaplot$lat>48.1]<-2
variable<-c(rep("Biomass",length(lat)/2),rep("Effort",length(lat)/2))

VBEffareaplot<-cbind(VBEffareaplot,lat,lon,month=indmonth,yr=indyr,variable)
#Effareaplot<-cbind(Effareaplot,lat,lon,month=indmonth,yr=indyr,variable)

sdvalue<-c(
      VBEffareaplot$value[VBEffareaplot$variable=="Biomass"],
      (VBEffareaplot$value[VBEffareaplot$variable=="Effort"]/mean(VBEffareaplot$value[VBEffareaplot$variable=="Effort"]))
      *mean(VBEffareaplot$value[VBEffareaplot$variable=="Biomass"]))
      
VBEffareaplot<-cbind(VBEffareaplot,sdvalue)

summary(VBEffareaplot)
anos<- VBEffareaplot$yr
saveLatex(
#if you don't have latex installed you need to install it or use other function such as saveGIF
  for(i in 121:max(ntsp))
  {
    

      ex1<-VBEffareaplot[VBEffareaplot$time==i ,]

      minVB<-min(VBEffareaplot$value[VBEffareaplot$variable=="Biomass"])
      maxVB<-max(VBEffareaplot$value[VBEffareaplot$variable=="Biomass"])
      minEff<-min(VBEffareaplot$value[VBEffareaplot$variable=="Effort"])
      maxEff<-max(VBEffareaplot$value[VBEffareaplot$variable=="Effort"])
      


      p3<-  ggmap(basemap,
          extent = "panel",
          ylab = "Latitude",
          xlab = "Longitude")
      p3 <- p3 + geom_line(y=48.5, linetype=2, colour="grey60")
      p3 <- p3 + geom_point(alpha=0.8,aes(size=sdvalue, shape=variable, color=variable),data=ex1) 
      p3 <- p3 + scale_shape_manual(values=c(16,21)) + scale_fill_discrete(na.value=NA, guide="none")
      p3 <- p3 + scale_color_manual(values=c("red", "black")) + scale_size_area(guide = "none")
      p3 <- p3 + labs(title=paste(meses[indmonth[i]],", year: ",anos[indmonth[i]]))     
      p3 <- p3 + facet_wrap(~Scenario)
      
      print(p3)     

  }
  ,
  pdflatex = "/usr/texbin/pdflatex")













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

meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
anos<- VBEffareaplot$yr

basemap<-get_map(location = c(lon = -125, lat = 45),
    zoom = 5, maptype = "terrain")
setwd("/Users/catarinawor/Documents/hake/PICES_conference/presentation/Fscn_anime")

#saveLatex(
#if you don't have latex installed you need to install it or use other function such as saveGIF
#  for(i in 121:max(ntsp))
#  {
    i=6

      ex1<-VBareaplot[VBareaplot$time==i ,]
      ex2<-Effareaplot[Effareaplot$time==i ,]
      minVB<-min(VBareaplot$VB)
      maxVB<-max(VBareaplot$VB)
      minEff<-min(Effareaplot$Eff)
      maxEff<-max(Effareaplot$Eff)
      
      p3<-  ggmap(basemap,
          extent = "panel",
          ylab = "Latitude",
          xlab = "Longitude")
      p3 <- p3 + geom_line(y=48.5, linetype=2, colour="grey60")
      p3 <- p3 + geom_point(alpha=0.8,aes(size=VB, group=Scenario),data=ex1) 
      p3 <- p3 + facet_wrap(~Scenario)
      p3 <- p3 + labs(title=paste(meses[indmonth[i]],", year: ",anos[indmonth[i]]))
      p3 <- p3 + scale_size_area(limits=c(minVB,maxVB),max_size = 10,breaks = c(50,100,200,400,600), labels = c(50,100,200,400,600), name = "vulnerable biomass")
       

      p2<-  ggmap(basemap,
          extent = "panel",
          ylab = "Latitude",
          xlab = "Longitude")
      p2 <- p2 + geom_line(y=48.5, linetype=2, colour="grey60")
      p2 <- p2 + geom_point(alpha=0.8,aes(size=Eff, group=Scenario),data=ex2,shape=1)
      p2 <- p2 + facet_wrap(~Scenario)     
      p2 <- p2 + scale_size_area(limits=c(minEff,maxEff),max_size = 10,  name = "Effort")
      

      multiplot(p3,p2) 


      

  #}
  #,
  #pdflatex = "/usr/texbin/pdflatex")

















EffNatAgeloF <-matrix(loF$EffNatAge,ncol=(loF$nage-loF$sage+3),dimnames=list(1:(((loF$nyr-loF$syr+1)*(loF$nmon-loF$smon+1))*loF$nations),c("tstp","nation",loF$sage:loF$nage)))
EffNatAgeloF<-as.data.frame(EffNatAgeloF)
 
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
