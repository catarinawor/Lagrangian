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
#Has the biomass stabilized -non error only
names(baseF)
plot(baseF$SB)


#remodel data base for plotting
ntsp <- 1:((baseF$nyr-baseF$syr+1)*(baseF$nmon-baseF$smon+1))
nage <- baseF$nage-baseF$sage+1
ages <- baseF$nage:baseF$sage
indmonth <- rep(baseF$smon:baseF$nmon,(baseF$nyr-baseF$syr+1))
indyr <-  rep(baseF$syr:baseF$nyr,each=(baseF$nmon-baseF$smon+1))
meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")


VBareabaseF<-matrix(baseF$VBarea, ncol=(baseF$narea-baseF$sarea+1),dimnames=list(ntsp,baseF$sarea:baseF$narea))
EffareabaseF<-matrix(baseF$Effarea, ncol=(baseF$narea-baseF$sarea+1),dimnames=list(ntsp,baseF$sarea:baseF$narea))

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
  for(i in 601:624){
     
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
      #p2 <- p2 + continuous_scale(,scale_name="size")
      p2 <- p2 + labs(title=meses[indmonth[i]],x="Longitude",y="Latitude") 
      
      print(p2) 

      ggsave(filename =paste0("Rplot",i-600,".png"))

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
print(p)
ggsave(filename ="availEffect.png") 



names(loF)

VBarealoF<-matrix(loF$VBarea, ncol=(loF$narea-loF$sarea+1),dimnames=list(ntsp,loF$sarea:loF$narea))
VBareahiF<-matrix(hiF$VBarea, ncol=(hiF$narea-hiF$sarea+1),dimnames=list(ntsp,hiF$sarea:hiF$narea))

EffarealoF<-matrix(loF$Effarea, ncol=(loF$narea-loF$sarea+1),dimnames=list(ntsp,loF$sarea:loF$narea))
EffareahiF<-matrix(hiF$Effarea, ncol=(hiF$narea-hiF$sarea+1),dimnames=list(ntsp,hiF$sarea:hiF$narea))

EffplotloF<-cbind(melt(EffarealoF),rep("High Biomass",nrow(melt(EffarealoF))))
EffplothiF<-cbind(melt(EffareahiF),rep("Low Biomass",nrow(melt(EffareahiF))))

VBplotloF<-cbind(melt(VBarealoF),rep("High Biomass",nrow(melt(VBarealoF))))
VBplothiF<-cbind(melt(VBareahiF),rep("High Biomass",nrow(melt(VBareahiF))))

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

setwd("/Users/catarinawor/Documents/hake/PICES_conference/presentation/Fscn_anime")

anos<- VBEffareaplot$yr

meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

basemap<-get_map(location = c(lon = -125, lat = 45),
    zoom = 5, maptype = "terrain")


#saveLatex( # not using save latex anymore due to poor resolution
#if you don't have latex installed you need to install it or use other function such as saveGIF
  for(i in 601:624)
  {
    

      ex1<-VBEffareaplot[VBEffareaplot$time==i ,]

      

      p3<-  ggmap(basemap,
          extent = "panel",
          ylab = "Latitude",
          xlab = "Longitude")
      p3 <- p3 + geom_line(y=48.5, linetype=2, colour="grey60")
      p3 <- p3 + geom_point(alpha=0.8,aes(size=sdvalue, shape=variable, color=variable),data=ex1) 
      p3 <- p3 + scale_shape_manual(values=c(16,21)) + scale_fill_discrete(na.value=NA, guide="none")
      p3 <- p3 + scale_color_manual(values=c("red", "black")) + scale_size_area(guide = "none")
      p3 <- p3 + labs(title=meses[indmonth[i]],x="Longitude",y="Latitude")     
      p3 <- p3 + facet_wrap(~Scenario)
      
      print(p3)

      ggsave(filename =paste0("Rplot",i-600,".png"))     

  }
  #,
  #pdflatex = "/usr/texbin/pdflatex")





#======================================================================== 
# Graphs for environmental influence scenario 
#======================================================================== 
#read data from environmental scenario
setwd("/Users/catarinawor/Documents/Lagrangian/environ")
evr<-read.rep("lagrangian_OM.rep")

names(evr)

ntsp <- 1:((evr$nyr-evr$syr+1)*(evr$nmon-evr$smon+1))
ages <- evr$nage:evr$sage
indmonth <- rep(evr$smon:evr$nmon,(evr$nyr-evr$syr+1))
indyr <-  rep(evr$syr:evr$nyr,each=(evr$nmon-evr$smon+1))


Emaxpos<-evr[["maxpos50E"]]

mxMAt<-matrix(NA, ncol=length(ages), nrow=length(Emaxpos))

for(i in 1:length(Emaxpos)){
  mxMAt[i,]<- (1/(1+exp(-(ages-Emaxpos[i])/evr$maxPossd)))*(evr$narea-evr$sarea)+evr$sarea
}

mxMAt<-mxMAt[order(Emaxpos),]
cores<-rainbow(length(Emaxpos), alpha = .7,end=4/6)

setwd("/Users/catarinawor/Documents/hake/PICES_conference/presentation")

pdf("movAtAge_clim.pdf",width = 7, height = 6)

plot(ages,mxMAt[1,], col= cores[1], lwd=3, type="l", ylab="AREA", xlab="AGE",font=2, font.lab=2,cex=1.5,cex.axis=1.5,cex.lab=1.5)
for(i in 1:length(Emaxpos)){
  lines(ages,mxMAt[i,], col= cores[i], lwd=3, type="l")
}
legend("topleft",c("warm","cold"), bty="n", col=c(cores[1],cores[length(cores)]),
  lwd=3,pch=NULL, text.font=2, cex=1.3)

dev.off()

#==============================================

#Animation with warmest and coldest year on simulation record. 

VBareaevr0<-matrix(evr$VBarea, ncol=length(evr$sarea:evr$narea),dimnames=list(ntsp,evr$sarea:evr$narea))
VBareaevr<-cbind(VBareaevr0[601:nrow(VBareaevr0),])
head(VBareaevr)

VBareaMax<-VBareaevr[which(indyr[601:length(indyr)]==which.max(Emaxpos[50:length(Emaxpos)])+50),]
VBareaMin<-VBareaevr[which(indyr[601:length(indyr)]==which.min(Emaxpos[50:length(Emaxpos)])+50),]


VBplotMax<-cbind(melt(VBareaMax),rep("Cold",nrow(melt(VBareaMax))))
VBplotMin<-cbind(melt(VBareaMin),rep("Warm",nrow(melt(VBareaMin))))

Effareaevr<-matrix(evr$Effarea, ncol=length(evr$sarea:evr$narea),dimnames=list(ntsp,evr$sarea:evr$narea))

EffareaMax<-Effareaevr[which(indyr[601:length(indyr)]==which.max(Emaxpos[50:length(Emaxpos)])+50),]
EffareaMin<-Effareaevr[which(indyr[601:length(indyr)]==which.min(Emaxpos[50:length(Emaxpos)])+50),]

EffplotMax<-cbind(melt(EffareaMax),rep("Cold",nrow(melt(EffareaMax))))
EffplotMin<-cbind(melt(EffareaMin),rep("Warm",nrow(melt(EffareaMin))))

names(VBplotMax)<- c("time","area", "value","Scenario")
names(VBplotMin)<- c("time","area", "value","Scenario")

names(EffplotMax)<- c("time","area", "value","Scenario")
names(EffplotMin)<- c("time","area", "value","Scenario")

VBEffareaplotEVR<-rbind(VBplotMax,VBplotMin,EffplotMax,EffplotMin)

lat<-VBEffareaplotEVR$area
lon<-rep(-131,length(lat))
nation<-rep(1,nrow(VBEffareaplotEVR))
nation[VBEffareaplotEVR$lat>48.1]<-2
variable<-c(rep("Biomass",length(lat)/2),rep("Effort",length(lat)/2))

VBEffareaplotEVR<-cbind(VBEffareaplotEVR,lat,lon,month=1:12,variable)

sdvalue<-c(
      VBEffareaplotEVR$value[VBEffareaplotEVR$variable=="Biomass"],
      (VBEffareaplotEVR$value[VBEffareaplotEVR$variable=="Effort"]/mean(VBEffareaplotEVR$value[VBEffareaplotEVR$variable=="Effort"]))
      *mean(VBEffareaplotEVR$value[VBEffareaplotEVR$variable=="Biomass"]))
      
VBEffareaplotEVR<-cbind(VBEffareaplotEVR,sdvalue)

setwd("/Users/catarinawor/Documents/hake/PICES_conference/presentation/evr_anime")

meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

basemap<-get_map(location = c(lon = -125, lat = 45),
    zoom = 5, maptype = "terrain")


#saveLatex( # not using save latex anymore due to poor resolution
#if you don't have latex installed you need to install it or use other function such as saveGIF
  for(i in 1:12){
    
      ex1<-VBEffareaplotEVR[VBEffareaplotEVR$month==i ,]

      

      p3<-  ggmap(basemap,
          extent = "panel",
          ylab = "Latitude",
          xlab = "Longitude")
      p3 <- p3 + geom_line(y=48.5, linetype=2, colour="grey60")
      p3 <- p3 + geom_point(alpha=0.8,aes(size=sdvalue*20, shape=variable, color=variable),data=ex1) 
      p3 <- p3 + scale_shape_manual(values=c(16,21)) + scale_fill_discrete(na.value=NA, guide="none")
      p3 <- p3 + scale_color_manual(values=c("red", "black")) + scale_size_area(guide = "none")
      p3 <- p3 + labs(title=meses[indmonth[i]],x="Longitude",y="Latitude")     
      p3 <- p3 + facet_wrap(~Scenario)
      
      print(p3)

      ggsave(filename =paste0("Rplot",i,".png"))     

  }
  #,
  #pdflatex = "/usr/texbin/pdflatex")

plot(evr$SB)
VBEffareaplotEVR[VBEffareaplotEVR$Scenario=="low",]

envir<-evr$envir[50:100]
cols <- c("blue", "red")[(envir > 0) + 1]  

## Pass the colors in to barplot()

setwd("/Users/catarinawor/Documents/hake/PICES_conference/presentation")

pdf("envIndex.pdf",width = 7, height = 7)   
barplot(envir, col = cols,border=NA,ylim=c(-2,2))
dev.off()

#====================================================================
#availability plot
#done with exaggerated temperature gradients

system("./lagrangian_OM -ind warm.dat")
setwd("/Users/catarinawor/Documents/Lagrangian")
warm<-read.rep("lagrangian_OM.rep")

system("./lagrangian_OM -ind cold.dat")
setwd("/Users/catarinawor/Documents/Lagrangian")
cold<-read.rep("lagrangian_OM.rep")

names(warm)
dim(warm$VBarea)

warmVBarea<-warm$VBarea[1189:1200,]
coldVBarea<-cold$VBarea[1189:1200,]

bdr<-48
Nat1<-warm$sarea:bdr
Nat2<-(bdr+1):warm$narea


VBareawarmNat1<-apply(warmVBarea[,Nat1-30],1,sum)
VBareawarmNat2<-apply(warmVBarea[,Nat2-30],1,sum)

VBareawarmNat1p<-VBareawarmNat1/(VBareawarmNat1+VBareawarmNat2)
VBareawarmNat2p<-VBareawarmNat2/(VBareawarmNat1+VBareawarmNat2)


VBareacoldNat1<-apply(coldVBarea[,Nat1-30],1,sum)
VBareacoldNat2<-apply(coldVBarea[,Nat2-30],1,sum)

VBareacoldNat1p<-VBareacoldNat1/(VBareacoldNat1+VBareacoldNat2)
VBareacoldNat2p<-VBareacoldNat2/(VBareacoldNat1+VBareacoldNat2)

VBNation<-c(VBareawarmNat1,VBareawarmNat2,VBareacoldNat1,VBareacoldNat2,VBareawarmNat1p,VBareawarmNat2p,VBareacoldNat1p,VBareacoldNat2p)
nation<-rep(c(rep(1,warm$nmon),rep(2,warm$nmon)),4)
scenario<-rep(c(rep("Warm",nrow(warmVBarea)*2), rep("Cold",nrow(coldVBarea)*2)),2)
month<-rep(1:12,8)
measure<-c(rep("Nominal",length(month)/2),rep("Proportion",length(month)/2))

## need to calculate measure still

df<- data.frame(VBNation,month,nation,scenario,measure)

p <- ggplot(df, aes(x=as.factor(month), y=VBNation, fill=as.factor(nation)))
p <- p + geom_bar(stat = "identity")
p <- p + facet_grid(measure~scenario,scales = "free_y")
p <- p + labs(x="month", y="vulnerable biomass",fill="nation")
p <- p + theme_bw()
print(p)
setwd("/Users/catarinawor/Documents/hake/PICES_conference/presentation")
ggsave(filename ="availEffectEvr.png") 
