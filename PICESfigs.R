#=========================================================================
# Routine for generating figures for PICES presentations
# Author: Catarina Wor
# Aug 5th 2015
#=========================================================================

rm(list=ls()); 

setwd("/Users/catarinawor/Documents/Lagrangian/")
source("read.admb.R")
baseF<-read.rep("lagrangian_OM.rep")

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

#======================================================================== 
# Graphs for base case scenario
#======================================================================== 

#remodel data base for plotting
ntsp <- 1:((baseF$nyr-baseF$syr+1)*(baseF$nmon-baseF$smon+1))
nage <- baseF$nage-baseF$sage+1
ages <- baseF$nage:baseF$sage
indmonth <- rep(baseF$smon:baseF$nmon,(baseF$nyr-baseF$syr+1))
indyr <-  rep(baseF$syr:baseF$nyr,each=(baseF$nmon-baseF$smon+1))
meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")



VBareabaseF<-matrix(baseF$VBarea, ncol=(baseF$narea-baseF$sarea+1),dimnames=list(ntsp,baseF$sarea:baseF$narea))
EffareabaseF<-matrix(baseF$Effarea, ncol=(baseF$narea-loF$sarea+1),dimnames=list(ntsp,baseF$sarea:baseF$narea))

VBplotbaseF<-melt(VBareabaseF)
EffplotbaseF<-melt(EffareabaseF)

names(VBplotbaseF)<- c("time","area", "value")
names(EffplotbaseF)<- c("time","area", "value")

VBEffareaplotbase<-rbind(VBplotbaseF,EffplotbaseF)

lat<-VBEffareaplotbase$area
lon<-rep(-131,length(lat))
nation<-rep(1,nrow(VBEffareaplotbase))
nation[VBEffareaplotbase$lat>48.1]<-2
variable<-c(rep("Biomass",length(lat)/2),rep("Effort",length(lat)/2))

VBEffareaplotbase<-cbind(VBEffareaplotbase,lat,lon,month=indmonth,yr=indyr,variable)

#rescale variables to be plotted in the same graph
sdvalue<-c(
      VBEffareaplotbase$value[VBEffareaplotbase$variable=="Biomass"],
      (VBEffareaplotbase$value[VBEffareaplotbase$variable=="Effort"]/mean(VBEffareaplotbase$value[VBEffareaplotbase$variable=="Effort"]))
      *mean(VBEffareaplotbase$value[VBEffareaplotbase$variable=="Biomass"]))
      
#dataframe to use in plotting
VBEffareaplotbase<-cbind(VBEffareaplotbase,sdvalue)


yr<- VBEffareaplotbase$yr
basemap<-get_map(location = c(lon = -125, lat = 45),
    zoom = 5, maptype = "terrain")
setwd("/Users/catarinawor/Documents/hake/PICES_conference/presentation/Fbase_anime")

#saveLatex( #not using savelatex anymore due to poor resolution
#if you don't have latex installed you need to install it or use other function such as saveGIF
  for(i in 121:max(ntsp)){
     
      ex1<-VBEffareaplotbase[VBEffareaplotbase$time==i ,]
      graphics.off()
      
      p2<-  ggmap(basemap,
          extent = "panel",
          ylab = "Latitude",
          xlab = "Longitude")
      p2 <- p2 + geom_line(y=48.5, linetype=2, colour="grey60")
      p2 <- p2 + geom_point(alpha=0.8,aes(size=sdvalue*10, shape=variable, color=variable),data=ex1) 
      p2 <- p2 + scale_shape_manual(values=c(16,21)) + scale_fill_discrete(na.value=NA, guide="none")
      p2 <- p2 + scale_color_manual(values=c("red", "black")) + scale_size_area(guide = "none")
      p2 <- p2 + labs(title=paste(meses[indmonth[i]],", year: ",yr[indmonth[i]]))
      
      print(p2) 

      ggsave(filename =paste0("Rplot",i-120,".png"))

      #png(filename = paste0("maracuja",i,".png"),width = 960, height = 960, units = "px", pointsize = 12)

  }
  

  #,
  #pdflatex = "/usr/texbin/pdflatex",)


#======================================================================== 
# Graphs for high and low abundance scenarios case scenario
#======================================================================== 

ntsp <- 1:((loF$nyr-loF$syr+1)*(loF$nmon-loF$smon+1))
nage <- loF$nage-loF$sage+1
ages <- loF$nage:loF$sage
indmonth <- rep(loF$smon:loF$nmon,(loF$nyr-loF$syr+1))
indyr <-  rep(loF$syr:loF$nyr,each=(loF$nmon-loF$smon+1))

maxPos501<-loF$maxPos501
maxPos502<-loF$maxPos502
maxPossd1<-loF$maxPossd1
maxPossd2<-loF$maxPossd2

setwd("/Users/catarinawor/Documents/hake/PICES_conference/presentation")

pdf("movAtAge.pdf",width = 7, height = 7)

plot(ages,(1/(1+exp(-(ages-maxPos501)/maxPossd1)))*(loF$narea-loF$sarea)+loF$sarea, ylab="AREA",xlab="AGE",
  lwd=4, type="l", col="black", font=2, font.lab=2,cex=1.5,cex.axis=1.5,cex.lab=1.5)
lines(ages,1/(1+exp(-(ages-maxPos502)/maxPossd2))*(loF$narea-loF$sarea)+loF$sarea, col="chartreuse3",lwd=4)
legend("topleft",c("high abundance", "low abundance"), bty="n", col=c("black","chartreuse3"),
  lwd=3,pch=NULL,text.col=c("black","chartreuse3"), text.font=2, cex=1.3)

dev.off()




