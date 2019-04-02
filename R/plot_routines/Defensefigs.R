#=========================================================================
# Routine for generating figures for Defense presentation
# Author: Catarina Wor
# Feb 2018
# In Nov 2018 added more figures for the PBS seminar
#=========================================================================

source("/Users/catarinawor/Documents/Lagrangian/R/read.admb.R")
setwd("/Users/catarinawor/Documents/Lagrangian/admb/OM/simple/")
OM<-read.rep("lagrangian_OM.rep")


#========================================================================
library(ggplot2)
library(reshape2)
library(animation)
library(ggmap)
library(plyr)

#======================================================================== 
# Graphs for base case scenario
#======================================================================== 

#======================================================
#remodel data base for plotting maps
names(OM)
ntsp <- 1:length(OM$indyr)
ages <- OM$sage:OM$nage
nage <- OM$nage
indmonth <- OM$indyr
indyr <-  OM$indyr
meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")


VBarea<-matrix(OM$VBarea, ncol=(OM$narea-OM$sarea+1),dimnames=list(ntsp,OM$sarea:OM$narea))
Effarea<-matrix(OM$Effarea, ncol=(OM$narea-OM$sarea+1),dimnames=list(ntsp,OM$sarea:OM$narea))

VBplot<-melt(VBarea)
Effplot<-melt(Effarea)

names(VBplot)<- c("time","area", "value")
names(Effplot)<- c("time","area", "value")

VBEffareaplot<-rbind(VBplot,Effplot)

lat<-VBEffareaplot$area
lon<-rep(-131,length(lat))
nation<-rep(1,nrow(VBEffareaplot))
nation[VBEffareaplot$lat>48.1]<-2
variable<-c(rep("Biomass",length(lat)/2),rep("Effort",length(lat)/2))

VBEffareaplot<-cbind(VBEffareaplot,lat,lon,month=indmonth,yr=indyr,variable)

#rescale variables to be plotted in the same graph
sdvalue<-c(
      VBEffareaplot$value[VBEffareaplot$variable=="Biomass"],
      (VBEffareaplot$value[VBEffareaplot$variable=="Effort"]/mean(VBEffareaplot$value[VBEffareaplot$variable=="Effort"&VBEffareaplot$value>0]))*mean(VBEffareaplot$value[VBEffareaplot$variable=="Biomass"]))
      
#dataframe to use in plotting
VBEffareaplot<-cbind(VBEffareaplot,sdvalue)


yr<- VBEffareaplot$yr
basemap<-get_stamenmap(location = c(lon = -125, lat = 45),
    zoom = 5, maptype = "terrain")
#setwd("/Users/catarinawor/Documents/hake/PICES_conference/presentation/Fbase_anime")
summary(VBEffareaplot)
#saveLatex( #not using savelatex anymore due to poor resolution
#if you don't have latex installed you need to install it or use other function such as saveGIF
  for(i in 601:624){
     
      ex1<-VBEffareaplot[VBEffareaplot$time==i ,]
      graphics.off()
      
      p2<-  ggmap(basemap,
          extent = "panel",
          ylab = "Latitude",
          xlab = "Longitude")
      p2 <- p2 + geom_line(y=48.5, linetype=2, colour="grey60")
      p2 <- p2 + geom_point(alpha=0.8,aes(size=sdvalue, shape=variable, color=variable),data=ex1) 
      p2 <- p2 + scale_shape_manual(values=c(16,21)) + scale_fill_discrete(na.value=NA, guide="none")
      p2 <- p2 + scale_color_manual(values=c("red", "black")) + scale_size_area(guide = "none") 
      #p2 <- p2 + continuous_scale(,scale_name="size")
      p2 <- p2 + labs(title=meses[indmonth[i]],x="Longitude",y="Latitude") 
      
      print(p2) 

      #ggsave(filename =paste0("Rplot",i-600,".png"))

      #png(filename = paste0("maracuja",i,".png"),width = 960, height = 960, units = "px", pointsize = 12)

  }
  

  #,
  #pdflatex = "/usr/texbin/pdflatex",)

#====================================================================
#Still figures
#====================================================================

library(ggplot2)
library(reshape2)
library(animation)
library(ggmap)
library(plyr)


setwd("/Users/catarinawor/Documents/Lagrangian/R")
source("read.admb.R")

sim <- read.rep("../admb/OM/simple/lagrangian_OM.rep")
sim_gtg <- read.rep("../admb/OM/gtg/lagrangian_OM_gtg.rep")

meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
names(sim)
dim(sim$propVBarea)
head(sim$propVBarea)
 
indmonth<-sim$indmonth

simvb<-(cbind(data.frame(time=sim$propVBarea[,1],
  area=sim$propVBarea[,2],
  type="single"),
  sim$propVBarea[,3:ncol(sim$propVBarea)]))


names(simvb)[4:ncol(simvb)]<-paste("A",sim$sage:sim$nage, sep="")

head(simvb)
molsimvb<-melt(simvb, id.vars=list("area", "time", "type"))
head(molsimvb)


names(sim_gtg)

gtgvb<-cbind(data.frame(time=sim_gtg$propVBarea[,1],
  group=sim_gtg$propVBarea[,2],
  area=sim_gtg$propVBarea[,3],  
  type="multiple"),
  sim_gtg$propVBarea[,4:ncol(sim_gtg$propVBarea)])

names(gtgvb)[5:ncol(gtgvb)]<-paste("A",sim_gtg$sage:sim_gtg$nage, sep="")
head(gtgvb)

ngtgvb<-aggregate(gtgvb[,5:ncol(gtgvb)], by=list(gtgvb$time,gtgvb$area,gtgvb$type), sum)

head(ngtgvb)

names(ngtgvb)[1]<-"time"
names(ngtgvb)[2]<-"area"
names(ngtgvb)[3]<-"type"

molgtgvb<-melt(ngtgvb, id.vars=list("area", "time", "type"))
head(molgtgvb)

df1<-rbind(molgtgvb,molsimvb)

df2<-df1[df1$variable=="A1"|df1$variable=="A6",]
#df2<-df1




df2$time2<-meses[indmonth[df2$time]]

df2<-df2[df2$time>1188,]
summary(df2)
df2<-arrange(transform(df2,Month=factor(time2,levels=meses)),Month)

df2$age<-as.factor(as.numeric(df2$variable))
#df2$values<-df2$value#/max(df2$value)
df2$values<-df2$value/max(df2$value)

summary(df2)




p <- ggplot(df2, aes(x=area,y=values))
p <- p+geom_line(aes(color=type, lty=age),alpha=0.8, size=1.)
p <- p+facet_wrap(~Month, ncol =3)
p <- p + geom_vline(aes(xintercept=49), linetype=3,alpha=0.8)
p <- p+ theme_bw(16)+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank())
#p <- p + scale_linetype_manual(breaks=c("1","6"), values=c(3,1))
p <- p + ylab("Relative Biomass")
p <- p + xlab("Latitude (areas)") + scale_color_brewer(palette="Dark2")
p
setwd("/Users/catarinawor/Documents/Thesis_presentation/Figures")
ggsave(file="mov_age6.pdf", plot=p)


dim(sim$Effarea)
head(sim$Effarea)
head(sim_gtg$Effarea)

effsim<-melt(sim$Effarea)

summary(effsim)

names(effsim)[1]<-"time"
names(effsim)[2]<-"area"
names(effsim)[3]<-"values"

summary(effsim)

effsim$area<-as.numeric(effsim$area)+29
summary(effsim$area)
effsim$dat<-"effort"
effsim$type<-"single"


effgtg<-melt(sim_gtg$Effarea)


names(effgtg)[1]<-"time"
names(effgtg)[2]<-"area"
names(effgtg)[3]<-"values"



effgtg$area<-as.numeric(effgtg$area)+29
effgtg$dat<-"effort"
effgtg$type<-"multiple"

summary(effgtg)
 
Eff<-rbind(effsim,effgtg)
Eff<-Eff[Eff$time>1188,]

Eff$time2<-meses[indmonth[Eff$time]]
Eff$values<-Eff$values/max(Eff$values)

Eff<-arrange(transform(Eff,Month=factor(time2,levels=meses)),Month)
summary(Eff)






p <- ggplot(df2, aes(x=area,y=values))
p <- p+geom_line(aes(color=type, lty=age),alpha=0.8, size=1.)
p <- p+facet_wrap(~Month, ncol =3)
p <- p + geom_vline(aes(xintercept=49), linetype=3,alpha=0.8)
p <- p+ theme_bw(16)+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank())
#p <- p + scale_linetype_manual(breaks=c("1","6"), values=c(3,1))
p <- p+ geom_bar(data=Eff,aes(x=area,y=values,fill=type),alpha=0.8,stat = "identity", position=position_dodge())
p <- p + ylab("Relative Biomass/Effort") + scale_fill_brewer(palette="Dark2")
p <- p + xlab("Latitude (areas)") + scale_color_brewer(palette="Dark2")
p
setwd("/Users/catarinawor/Documents/Thesis_presentation/Figures")
ggsave(file="mov_age6_effort.pdf", plot=p)



df2<-df1




df2$time2<-meses[indmonth[df2$time]]

df2<-df2[df2$time>1188,]
summary(df2)
df2<-arrange(transform(df2,Month=factor(time2,levels=meses)),Month)

df2$age<-as.factor(as.numeric(df2$variable))
#df2$values<-df2$value#/max(df2$value)
df2$values<-df2$value/max(df2$value)

summary(df2)





df3<- aggregate(df2$values, by=list( Month=df2$Month,area=df2$area, type=df2$type),sum)

summary(df3)
df3$values<-df3$x/max(df3$x)



p <- ggplot(df3, aes(x=area,y=values))
p <- p+geom_line(aes(color=type),alpha=0.8, size=1.)
p <- p+facet_wrap(~Month, ncol =3)
p <- p + geom_vline(aes(xintercept=48.5), linetype=3,alpha=0.8)
p <- p+ theme_bw(16)+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank())
p <- p + geom_bar(data=Eff,aes(x=area,y=values,fill=type),alpha=0.8,stat = "identity", position=position_dodge())
p <- p + ylab("Relative Biomass/Effort") + scale_fill_brewer(palette="Dark2")
p <- p + xlab("Latitude (areas)") + scale_color_brewer(palette="Dark2")
p
setwd("/Users/catarinawor/Documents/Thesis_presentation/Figures")
ggsave(file="VB_effort.pdf", plot=p)


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
sim_gtg <- read.rep("../admb/OM/gtg/lagrangian_OM_gtg.rep")

names(sim)


meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

head(tmpXplot)
tmpXplot<-data.frame(sim$propVBarea[sim$propVBarea[,1]>(max(sim$propVBarea[,1])-12),])
tmpXplot[,1]<-meses[tmpXplot[,1]-(max(sim$propVBarea[,1])-12)]
tmpXplot<-rename(tmpXplot, c("V1"="month", "V2"= "Latitude","V3"="1", "V4"="2", "V5"= "3",
  "V6"="4", "V7"="5", "V8"= "6","V9"="7", "V10"="8", "V11"= "9","V12"="10", "V13"="11", "V14"= "12",
  "V15"="13", "V16"="14", "V17"= "15","V18"="16", "V19"="17", "V20"= "18","V21"="19", "V22"="20"))


Xplot<-melt(tmpXplot, id=c("month", "Latitude"),variable.name="age")
Xplot$group<-as.factor(rep(0,length(Xplot$month)))
Xplot<-arrange(transform(Xplot,month=factor(month,levels=meses)),month)
head(Xplot)



summary(sim_gtg$propVBarea)


tmpXplotgtg<-data.frame(sim_gtg$propVBarea[sim_gtg$propVBarea[,1]>1188,])

tmpXplotgtg[,1]<-meses[tmpXplotgtg[,1]-(max(sim_gtg$propVBarea[,1])-12)]
tmpXplotgtg<-rename(tmpXplotgtg, c("V1"="month", "V2"="group", "V3"= "Latitude","V4"="1", "V5"="2", "V6"= "3",
  "V7"="4", "V8"="5", "V9"= "6","V10"="7", "V11"="8", "V12"= "9","V13"="10", "V14"="11", "V15"= "12",
  "V16"="13", "V17"="14", "V18"= "15","V19"="16", "V20"="17", "V21"= "18","V22"="19", "V23"="20"))
head(tmpXplotgtg)
dim(tmpXplotgtg)

Xplotgtg<-melt(tmpXplotgtg, id=c("month","group", "Latitude"),variable.name="age")
Xplotgtg[1:30,]

#Xplotgtg<-arrange(transform(Xplotgtg,month=factor(month,levels=meses)),month)


d1=Xplot[Xplot$age==5,]
d2=Xplotgtg[Xplotgtg$age==5,]

summary(d1)
d2[d2$month=="Jul",]



Xplot$model<-rep("single",length(Xplot$age))
Xplotgtg$model<-rep("multiple",length(Xplotgtg$age))

Xplot.v<-subset(Xplotgtg,select= value)
Xplot.f<-subset(Xplotgtg,select= -c(group,value))

newXplotgtg<- aggregate(x=Xplot.v,by=Xplot.f,FUN=sum)
head(newXplotgtg)


combXplotgtg<-data.frame(month=newXplotgtg$month,
                        Latitude=newXplotgtg$Latitude,
                        age=newXplotgtg$age,
                        value=newXplotgtg$value,
                        group=rep("-1",length(newXplotgtg$value)),
                        model=rep("sum_multiple",length(newXplotgtg$value)))
head(combXplotgtg)
head(Xplotgtg)
summary(combXplotgtg)
summary(Xplotgtg)
head(Xplot)

dfall<-rbind(Xplotgtg,Xplot,combXplotgtg)

#df<-rbind(d1,d2)
head(dfall)
summary(dfall)

dfalljul<-dfall[dfall$age==5&dfall$month=="Jul",]
dfalljul$Latitude<-as.numeric(dfalljul$Latitude)
#dfalljul$value<-dfalljul$value/max(dfalljul$value)

unique(dfalljul$group)

dfalljul[dfalljul$group==10,]


p <- ggplot(dfalljul) 
p <- p + geom_line(aes(x=Latitude, y=value, lty=group, colour=model),size=1.2)
p <- p + geom_vline(xintercept=48.5, linetype=3)
#p <- p + facet_wrap(~month,ncol=4)
p <- p + theme_bw()+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks.y = element_blank(),
        axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold")) 
p <- p + ylab("Biomass")
p <- p + scale_linetype_manual(breaks=c(as.factor(1:20),"0","-1"), values=c(5,1,rep(1,20)))
p <- p + guides(lty=FALSE) + scale_color_brewer(palette="Dark2")
p

setwd("/Users/catarinawor/Documents/Thesis_presentation/Figures")
ggsave(file="gtg_simple_comp.pdf", plot=p)


#====================================================================
#month by month - (almost all cohorts)

source("/Users/catarinawor/Documents/Lagrangian/R/read.admb.R")
setwd("/Users/catarinawor/Documents/Lagrangian/admb/OM/simple/")
sim<-read.rep("lagrangian_OM.rep")


meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")


tmpXplot<-data.frame(sim$propVBarea[sim$propVBarea[,1]>(max(sim$propVBarea[,1])-12),])
tmpXplot[,1]<-meses[tmpXplot[,1]-(max(sim$propVBarea[,1])-12)]
tmpXplot<-rename(tmpXplot, c("V1"="month", "V2"= "Latitude","V3"="1", "V4"="2", "V5"= "3",
  "V6"="4", "V7"="5", "V8"= "6","V9"="7", "V10"="8", "V11"= "9","V12"="10", "V13"="11", "V14"= "12",
  "V15"="13", "V16"="14", "V17"= "15","V18"="16", "V19"="17", "V20"= "18","V21"="19", "V22"="20"))
dim(tmpXplot)
head(tmpXplot)


Xplot<-melt(tmpXplot, id=c("month", "Latitude"),variable.name="age")
Xplot$group<-as.factor(rep(0,length(Xplot$month)))
Xplot<-arrange(transform(Xplot,month=factor(month,levels=meses)),month)
aggregate(Xplot$value,list(Xplot$age),sum) 


head(Xplot)
dim(Xplot)



Xplot10<-Xplot[as.numeric(Xplot$age)%%2==0&as.numeric(Xplot$age)<20,] 

aggregate(Xplot10$value,list(Xplot10$age),sum) 


p2 <- ggplot(Xplot10) 
p2 <- p2 + geom_line(aes(x=Latitude, y=value, colour=age),size=1.8)
p2 <- p2 + geom_vline(xintercept=48.5, linetype=3,size=1.8)
p2 <- p2 + facet_wrap(~month,ncol=4)
p2 <- p2 + theme_bw(18)+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks.y = element_blank(),
        axis.text.y=element_blank(),axis.text.x=element_text(face="bold"),
        axis.title=element_text(face="bold"), strip.text=element_text(face="bold")) 
p2 <- p2 + ylab("Biomass")
p2
ggsave(filename ="age_move.pdf")

#===============================
#animation version


setwd("/Users/catarinawor/Documents/Lagrangian/anime")
library("animation")

ani.options(ani.dev = "pdf", ani.type = "pdf",ani.width=8, ani.height=6)
saveLatex(
#if you don't have latex installed you need to install it or use other function such as saveGIF

for(i in 1:length(meses)){

    Xplotmes<-Xplot10[Xplot10$month==meses[i],]


    p2 <- ggplot(Xplotmes) 
    p2 <- p2 + geom_line(aes(x=Latitude, y=value, colour=age),size=1.8)
    p2 <- p2 + geom_vline(xintercept=48.5, linetype=3,size=1.8)
    p2 <- p2 + theme_bw(18)+theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),axis.ticks.y = element_blank(),
            axis.text.y=element_blank(),
            axis.text.x=element_text(face="bold"),axis.title=element_text(face="bold"),
            plot.title=element_text(face="bold",hjust=.5) )
    p2 <- p2 + ylab("Biomass") + labs(title=meses[i])
    print(p2)
    #ggsave(filename =paste0("Rplot",i,".pdf")) 


}
,
pdflatex = "/usr/texbin/pdflatex",latex.filename = "movecurve.tex")



