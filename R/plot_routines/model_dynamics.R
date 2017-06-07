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
sim_gtg$propVBarea[1:2400,1:3]

tmpXplotgtg<-data.frame(sim_gtg$propVBarea[sim_gtg$propVBarea[,1]>( calc_numbers_at_age),])

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
setwd("/Users/catarinawor/Documents/hake/Lag_Model_paper")
ggsave(file="Figure1.pdf")

############need to run with  effort ############