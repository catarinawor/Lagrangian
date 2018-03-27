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

head(sim$propVBarea)
dim(sim$propVBarea)

sim$CatchNatAge

sim$yNage

plot(sim$SB)

meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

head()
head(tmpXplot)
tmpXplot<-data.frame(sim$propVBarea[sim$propVBarea[,1]>(max(sim$propVBarea[,1])-12),])
tmpXplot[,1]<-meses[tmpXplot[,1]-(max(sim$propVBarea[,1])-12)]
tmpXplot<-rename(tmpXplot, c("V1"="month", "V2"= "Latitude","V3"="1", "V4"="2", "V5"= "3",
  "V6"="4", "V7"="5", "V8"= "6","V9"="7", "V10"="8", "V11"= "9","V12"="10", "V13"="11", "V14"= "12",
  "V15"="13", "V16"="14", "V17"= "15","V18"="16", "V19"="17", "V20"= "18","V21"="19", "V22"="20"))

apply(tmpXplot[,-c(1:3)],2,sum)

Xplot<-melt(tmpXplot, id=c("month", "Latitude"),variable.name="age")
Xplot$group<-as.factor(rep(0,length(Xplot$month)))
Xplot<-arrange(transform(Xplot,month=factor(month,levels=meses)),month)
head(Xplot)


head(sim_gtg$propVBarea)

(sim_gtg$SB)
tail(sim_gtg$propVBarea)

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


head(Xplotgtg)

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

setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")
p <- ggplot(dfalljul) 
p <- p + geom_line(aes(x=Latitude, y=value, lty=group, colour=model),size=1.2)
p <- p + geom_vline(xintercept=48.5, linetype=3)
#p <- p + facet_wrap(~month,ncol=4)
p <- p + theme_bw()+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks.y = element_blank(),
        axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold")) 
p <- p + ylab("Biomass")
p <- p + scale_linetype_manual(breaks=c(as.factor(1:20),"0","-1"), values=c(5,1,rep(1,20)))
p <- p + scale_colour_grey(start = 0.6, end = 0.1) +guides(lty=FALSE)
p
#setwd("/Users/catarinawor/Documents/hake/Lag_Model_paper")
#ggsave(file="Figure1.pdf")

#presentation version

p <- ggplot(dfalljul) 
p <- p + geom_line(aes(x=Latitude, y=value, lty=group, colour=model),size=1.2)
p <- p + geom_vline(xintercept=48.5, linetype=3)
#p <- p + facet_wrap(~month,ncol=4)
p <- p + theme_bw(18)+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks.y = element_blank(),
        axis.text=element_text(size=18,face="bold"),axis.title=element_text(size=18,face="bold"),
         legend.text=element_text(size=18,face="bold"), legend.title=element_text(size=18,face="bold")) 
p <- p + ylab("Biomass")
p <- p + scale_linetype_manual(breaks=c(as.factor(1:20),"0","-1"), values=c(5,1,rep(1,20)))
p <- p + scale_colour_brewer(palette="Dark2") +guides(lty=FALSE)
p <- p + annotate("text", x = 54, y = 0.045, label = "Canada", fontface =2, size= theme_get()$text[["size"]]/2)
p <- p + annotate("text", x = 44, y = 0.045, label = "U.S.A.", fontface =2, size= theme_get()$text[["size"]]/2)
p



############need to run with  effort ############


setwd("/Users/catarinawor/Documents/Lagrangian/R")
source("read.admb.R")

sim <- read.rep("../admb/OM/simple/lagrangian_OM.rep")
sim_gtg <- read.rep("../admb/OM/gtg/lagrangian_OM_gtg.rep")

meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

indmonth<-sim$indmonth

names(sim)
dim(sim$propVBarea)
head(sim$propVBarea)
 

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

df2<-df1[df1$variable=="A1"|df1$variable=="A5",]



df2$time2<-meses[indmonth[df2$time]]

df2<-df2[df2$time>1188,]
summary(df2)
df2<-arrange(transform(df2,Month=factor(time2,levels=meses)),Month)

df2$age<-as.factor(as.numeric(df2$variable))
#df2$values<-df2$value#/max(df2$value)
df2$values<-df2$value/max(df2$value)

#df2$values[df2$variable=="A1"]<-df2$value[df2$variable=="A1"]/max(df2$value[df2$variable=="A1"])
#df2$values[df2$variable=="A5"]<-df2$value[df2$variable=="A5"]/max(df2$value[df2$variable=="A5"])


#df2$values[df2$variable=="A2"&df2$type=="gtg"]<-df2$value[df2$variable=="A2"&df2$type=="gtg"]/max(df2$value[df2$variable=="A2"&df2$type=="gtg"])
#df2$values[df2$variable=="A2"&df2$type=="simple"]<-df2$value[df2$variable=="A2"&df2$type=="simple"]/max(df2$value[df2$variable=="A2"&df2$type=="simple"])

#df2$values[df2$variable=="A5"&df2$type=="gtg"]<-df2$value[df2$variable=="A5"&df2$type=="gtg"]/max(df2$value[df2$variable=="A5"&df2$type=="gtg"])
#df2$values[df2$variable=="A5"&df2$type=="simple"]<-df2$value[df2$variable=="A5"&df2$type=="simple"]/max(df2$value[df2$variable=="A5"&df2$type=="simple"])



names(sim)
sim$Effarea
dim(sim$Effarea)
head(sim$Effarea)

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


p<-ggplot(df2, aes(x=area,y=values))
p<-p+geom_line(aes(color=type, lty=age),alpha=0.8, size=1.)
p<-p+facet_wrap(~Month, ncol =3)
p <- p + geom_vline(aes(xintercept=49), linetype=3,alpha=0.8)
p<-p+ geom_bar(data=Eff,aes(x=area,y=values,fill=type),alpha=0.8,stat = "identity", position=position_dodge())
p<-p+ theme_bw(16)+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank())
p <- p + scale_color_grey("model") + scale_fill_grey("model")
p <- p + scale_linetype_manual(breaks=c("1","5"), values=c(3,1))
p <- p + ylab("Relative Biomass/Effort")
p <- p + xlab("Latitude (areas)")
p
#setwd("/Users/catarinawor/Documents/hake/Lag_Model_paper")
#ggsave(file="Figure2.pdf")

?geom_bar

#presentation version


#no effort

df3<-df2[df2$Month=="Jan"|df2$Month=="Apr"|df2$Month=="Jul"|df2$Month=="Oct",]


p<-ggplot(df3, aes(x=area,y=values))
p<-p+geom_line(aes(color=type, lty=age),alpha=0.8, size=2.)
p <- p + annotate("text", x = 54, y = 0.85, label = "Canada", fontface =2, size= theme_get()$text[["size"]]/2)
p <- p + annotate("text", x = 44, y = 0.85, label = "U.S.A.", fontface =2, size= theme_get()$text[["size"]]/2)
p<-p+facet_wrap(~Month, ncol =2)
p <- p + geom_vline(aes(xintercept=49), linetype=3,alpha=0.8)
#p<-p+ geom_bar(data=Eff,aes(x=area,y=values,fill=type),alpha=0.8,stat = "identity", position=position_dodge())
p<-p+ theme_bw(18)+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(),
        axis.text=element_text(size=18,face="bold"),axis.title=element_text(size=18,face="bold"),
         legend.text=element_text(size=18,face="bold"), legend.title=element_text(size=18,face="bold"))
p <- p + scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2")
p <- p + scale_linetype_manual(breaks=c("1","5"), values=c(3,1))
p <- p + ylab("Relative Biomass/Effort")
p <- p + xlab("Latitude")
p


head(Eff)
Eff3<-Eff[Eff$time2=="Jan"|Eff$time2=="Apr"|Eff$time2=="Jul"|Eff$time2=="Oct",]


p<-ggplot(df3, aes(x=area,y=values))
p<-p+geom_line(aes(color=type, lty=age),alpha=0.8, size=2.)
p <- p + annotate("text", x = 54, y = 0.85, label = "Canada", fontface =2, size= theme_get()$text[["size"]]/2)
p <- p + annotate("text", x = 44, y = 0.85, label = "U.S.A.", fontface =2, size= theme_get()$text[["size"]]/2)
p<-p+facet_wrap(~Month, ncol =2)
p <- p + geom_vline(aes(xintercept=49), linetype=3,alpha=0.8)
p<-p+ geom_bar(data=Eff3,aes(x=area,y=values,fill=type),alpha=0.8,stat = "identity", position=position_dodge())
p<-p+ theme_bw(16)+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(),
        axis.text=element_text(size=18,face="bold"),axis.title=element_text(size=18,face="bold"),
         legend.text=element_text(size=18,face="bold"), legend.title=element_text(size=18,face="bold"))
p <- p + scale_color_brewer(palette="Dark2") + scale_fill_brewer(palette="Dark2")
p <- p + scale_linetype_manual(breaks=c("1","5"), values=c(3,1))
p <- p + ylab("Relative Biomass/Effort")
p <- p + xlab("Latitude")
p