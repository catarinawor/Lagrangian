#==============================================
#Title:lat_ct
#Author: Catarina Wor
#date: Oct 25th 2016
#Plots of real data! Catch vs Latitude for Pacific hake fleet
# code adapted from iSCAM stuff (Martell et al)
#==============================================
library(ggplot2)
library(reshape2)


dat<-read.csv("/Users/catarinawor/Documents/Lagrangian/R/real_data/HakeDataForCat.csv")

summary(dat)
names(dat)


dat$LatBin

dat$LatBinw<-floor(dat$LatBin)

p <- ggplot(dat,aes(x=as.factor(LatBinw), y=HakeWt))
p <- p + geom_boxplot()
p <- p + facet_wrap(~Year,scales="free")
p

p <- ggplot(dat,aes(x=as.factor(LatBinw), y=HakeWt))
p <- p + geom_boxplot() 
p <- p + facet_wrap(~Month)
p <- p + ylim(0,5000)
p

p <- ggplot(dat,aes(x=as.factor(LatBinw), y=HakeWt))
p <- p + geom_bar(stat="identity")
p <- p + facet_wrap(~Month)
p <- p + ylim(0,130000)
p


p <- ggplot(dat,aes(x=as.factor(LatBin), y=HakeWt))
p <- p + geom_bar(stat="identity")
p <- p + facet_wrap(~Month)
p <- p + ylim(0,30000)
p


p <- ggplot(dat,aes(x=as.factor(LatBinw), y=Pop))
p <- p + geom_bar(stat="identity")
p <- p + facet_wrap(~Month)
#p <- p + ylim(0,30000)
p

p <- ggplot(dat,aes(x=as.factor(LatBinw), y=Pop))
p <- p + geom_boxplot()
p <- p + facet_wrap(~Month)
p <- p + ylim(0,1)
p


p <- ggplot(dat,aes(x=LatBin, y=HakeWt))
p <- p + geom_bar(stat="identity")
p <- p + facet_wrap(~Year,scales="free")
p

datb<-melt(dat, id.vars = c("X", "Year", "Month","VesselType","LatBin","LatBinw","HakeWt"))

pp <- ggplot(datb,aes(x=LatBin, y=value, fill=variable))
pp <- pp + geom_bar(stat="identity")
pp <- pp + facet_wrap(~Year,scales="free")
pp

pp <- ggplot(datb,aes(x=LatBinw, y=value, fill=variable))
pp <- pp + geom_bar(stat="identity")
pp <- pp + facet_wrap(~Year,scales="free")
pp

pp <- ggplot(datb,aes(x=as.factor(LatBinw), y=value, fill=variable))
pp <- pp + geom_boxplot()
pp <- pp + facet_wrap(~Month)
pp




hakeprop<-dat$HakeWt/max(dat$HakeWt) 
widowprop<-dat$Widow/max(dat$Widow) 
Popprop<-dat$Pop/max(dat$Pop)  
DKBprop<-dat$DKB/max(dat$DKB) 

datc<-cbind(dat,hakeprop,widowprop,Popprop,DKBprop)
datall<-melt(datc, id.vars = c("X", "Year", "Month","VesselType","LatBin","LatBinw","HakeWt","Pop","DKB","Widow", "hakeprop"))

summary(datall)

pp <- ggplot(datall,aes(x=as.factor(LatBinw), y=value, fill=variable))
pp <- pp + geom_boxplot()
pp <- pp + facet_wrap(~Month)
pp



pp <- ggplot(datall,aes(x=LatBin, y=value, fill=variable))
pp <- pp + geom_bar(stat="identity")
pp <- pp + facet_wrap(~Year,scales="free")
pp



pp1 <- ggplot(datall,aes(x=LatBin, y=value, fill=variable))
pp1 <- pp1 + geom_bar(stat="identity")
pp1 <- pp1 + facet_wrap(~Month)
pp1





pp1 <- ggplot(datall[datall$Year==2008,],aes(x=LatBin, y=value, fill=variable))
pp1 <- pp1 + geom_bar(stat="identity")
pp1 <- pp1 + facet_wrap(~Month)
pp1 <- pp1 + geom_vline(xintercept = 42)
pp1 <- pp1 + geom_vline(xintercept = 46.2)
pp1

pp1 <- ggplot(datall[datall$Year==2009,],aes(x=LatBin, y=value, fill=variable))
pp1 <- pp1 + geom_bar(stat="identity")
pp1 <- pp1 + facet_wrap(~Month)
pp1 <- pp1 + geom_vline(xintercept = 42)
pp1 <- pp1 + geom_vline(xintercept = 46.2)
pp1

pp1 <- ggplot(datall[datall$Year==2010,],aes(x=LatBin, y=value, fill=variable))
pp1 <- pp1 + geom_bar(stat="identity")
pp1 <- pp1 + facet_wrap(~Month)
pp1 <- pp1 + geom_vline(xintercept = 42)
pp1 <- pp1 + geom_vline(xintercept = 46.2)
pp1

pp1 <- ggplot(datall[datall$Year==2011,],aes(x=LatBin, y=value, fill=variable))
pp1 <- pp1 + geom_bar(stat="identity")
pp1 <- pp1 + facet_wrap(~Month)
pp1 <- pp1 + geom_vline(xintercept = 42)
pp1 <- pp1 + geom_vline(xintercept = 46.2)
pp1

pp1 <- ggplot(datall[datall$Year==2012,],aes(x=LatBin, y=value, fill=variable))
pp1 <- pp1 + geom_bar(stat="identity")
pp1 <- pp1 + facet_wrap(~Month)
pp1 <- pp1 + geom_vline(xintercept = 42)
pp1 <- pp1 + geom_vline(xintercept = 46.2)
pp1


pp1 <- ggplot(datall[datall$Year==2013,],aes(x=LatBin, y=value, fill=variable))
pp1 <- pp1 + geom_bar(stat="identity")
pp1 <- pp1 + facet_wrap(~Month)
pp1 <- pp1 + geom_vline(xintercept = 42)
pp1 <- pp1 + geom_vline(xintercept = 46.2)
pp1

pp1 <- ggplot(datall[datall$Year==2014,],aes(x=LatBin, y=value, fill=variable))
pp1 <- pp1 + geom_bar(stat="identity")
pp1 <- pp1 + facet_wrap(~Month)
pp1 <- pp1 + geom_vline(xintercept = 42)
pp1 <- pp1 + geom_vline(xintercept = 46.2)
pp1

pp1 <- ggplot(datall[datall$Year==2015,],aes(x=LatBin, y=value, fill=variable))
pp1 <- pp1 + geom_bar(stat="identity")
pp1 <- pp1 + facet_wrap(~Month)
pp1 <- pp1 + geom_vline(xintercept = 42)
pp1 <- pp1 + geom_vline(xintercept = 46.2)
pp1

pp1 <- ggplot(datall[datall$Year==2016,],aes(x=LatBin, y=value, fill=variable))
pp1 <- pp1 + geom_bar(stat="identity")
pp1 <- pp1 + facet_wrap(~Month)
pp1 <- pp1 + geom_vline(xintercept = 42)
pp1 <- pp1 + geom_vline(xintercept = 46.2)
pp1


pp1 <- ggplot(datall,aes(x=LatBin, y=value, fill=variable))
pp1 <- pp1 + geom_bar(stat="identity")
pp1 <- pp1 + facet_wrap(~Month)
pp1








p1 <- ggplot(dat,aes(x=LatBin, y=Pop ))
p1 <- p1 + geom_bar(stat="identity")
p1 <- p1 + facet_wrap(~Year,scales="free")
p1

p2 <- ggplot(dat,aes(x=LatBin, y=Pop ))
p2 <- p2 + geom_bar(stat="identity")
p2 <- p2 + facet_wrap(~Month,scales="free")
p2



p3 <- ggplot(dat,aes(x=LatBin, y=DKB ))
p3 <- p3 + geom_bar(stat="identity")
p3 <- p3 + facet_wrap(~Year,scales="free")
p3

p4 <- ggplot(dat,aes(x=LatBin, y=DKB ))
p4 <- p4 + geom_bar(stat="identity")
p4 <- p4 + facet_wrap(~Month,scales="free")
p4

names(dat)
datall<-melt(dat, id.vars = c("X", "Year", "Month","VesselType","LatBin"))
summary(datall)


p <- ggplot(datall,aes(x=LatBin, y=value, fill=variable))
p <- p + geom_bar(stat="identity")
p <- p + facet_wrap(~Year,scales="free")
p


