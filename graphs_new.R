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
Xplot<-Xplot


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

setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")
p <- ggplot(Xplot) 
#p <- p + geom_line(aes(x=as.numeric(area), y=value, colour=age, alpha=0.5))
p <- p + geom_bar(aes(x=(area), y=value,fill=age),stat="identity", alpha=0.5)
p <- p + facet_wrap(~month,ncol=4)
p <- p + geom_vline(xintercept=48.5, linetype=3,alpha=0.3)
p <- p + theme_bw()+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title.x= element_blank())
ggsave(file="gtg_bars_expl_graph.pdf")

setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")

p <- ggplot(Xplot) 
#p <- p + geom_line(aes(x=as.numeric(area), y=value, colour=age, alpha=0.5))
p <- p + geom_bar(aes(x=(area), y=value,fill=age),stat="identity", alpha=0.5)
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













dim(sim_gtg$Nage)

sim_gtg$Nage[,1]
sum(sim_gtg$Effarea)
NallG<-matrix(NA, nrow=1200,ncol=length(sim_gtg$sage:sim_gtg$nage))

for(i in 1:1200){
	NallG[i,]<-apply(sim_gtg$Nage[sim_gtg$Nage[,1]==i,],2,sum)[-c(1,2)]
}

NallG[,1]

plot(apply(sim_gtg$tVBarea[,-1],1,sum), type="l")
plot(apply(NallG,1,sum), type="l")
abline(v=c(13,25,37,49))



tmpXplot<-data.frame(sim_gtg$tVBarea[sim_gtg$tVBarea[,1]>(max(sim_gtg$tVBarea[,1])-12),])
tmpXplot[,1]<-meses[tmpXplot[,1]-(max(sim_gtg$tVBarea[,1])-12)]

tmpXplot<-rename(tmpXplot, c("V1"="month","V2"="30", "V3"="31", "V4"= "32",
  "V5"="33", "V6"="34", "V7"= "35","V8"="36", "V9"="37", "V10"= "38","V11"="39", "V12"="40", "V13"= "41",
  "V14"="42", "V15"="43", "V16"= "44","V17"="45", "V18"="46", "V19"= "47","V20"="48", "V21"="49",
  "V22"="50", "V23"="51", "V24"= "52","V25"="53", "V26"="54", "V27"= "55","V28"="56", "V29"="57",
   "V30"= "58","V31"="59", "V32"="60"))
tmpXplot<-tmpXplot[,-c(2,32)]

Xplot<-melt(tmpXplot, id=c("month"),variable.name="area")

Xplot<-arrange(transform(Xplot,month=factor(month,levels=meses)),month)

p <- ggplot(Xplot) 
p <- p + geom_bar(aes(x=as.numeric(area), y=value,fill="value"),stat="identity", alpha=0.5)
p <- p + facet_wrap(~month,ncol=4)
p


tmpXplote<-data.frame(sim_gtg$tVBarea[205:216,])
tmpXplote[,1]<-meses[tmpXplot[,1]-(max(sim_gtg$tVBarea[,1])-12)]

tmpXplot<-rename(tmpXplot, c("V1"="month","V2"="30", "V3"="31", "V4"= "32",
  "V5"="33", "V6"="34", "V7"= "35","V8"="36", "V9"="37", "V10"= "38","V11"="39", "V12"="40", "V13"= "41",
  "V14"="42", "V15"="43", "V16"= "44","V17"="45", "V18"="46", "V19"= "47","V20"="48", "V21"="49",
  "V22"="50", "V23"="51", "V24"= "52","V25"="53", "V26"="54", "V27"= "55","V28"="56", "V29"="57",
   "V30"= "58","V31"="50", "V32"="60"))

Xplot<-melt(tmpXplot, id=c("month"),variable.name="area")

Xplot<-arrange(transform(Xplot,month=factor(month,levels=meses)),month)

p <- ggplot(Xplot) 
p <- p + geom_bar(aes(x=(area), y=value,fill="value"),stat="identity", alpha=0.5)
p <- p + facet_wrap(~month,ncol=4)
p <- p + geom_vline(xintercept=18.5, linetype=3,alpha=0.3)
p <- p + theme_bw()+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title.x= element_blank())
p


tmpEffplot<-data.frame(month= 1:12,matrix(sim_gtg$Effarea[(nrow(sim_gtg$Effarea)-11):nrow(sim_gtg$Effarea),],ncol=31,
  dimnames=list(1:12,(sim_gtg$sarea:sim_gtg$narea))))
tmpEffplot<-setnames(tmpEffplot, old = paste("X",30:60, sep=""), new = as.character(30:60))

Effplot<-melt(tmpEffplot, id="month",variable.name="Latitude",value.name="effort")
Effplot$month<-meses[Effplot$month]
Effplot$effort<-Effplot$effort/max(Effplot$effort)*max(Xplot$value)
Effplot$Latitude<-as.numeric(Effplot$Latitude)+sim_gtg$sarea-1
df<-merge(x = Xplot, y = Effplot, by = c("Latitude","month"), all = TRUE)
head(df)

setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")
p <- ggplot(Xplot) 
p <- p + geom_line(aes(x=Latitude, y=value, colour=group))
p <- p + geom_vline(xintercept=48.5, linetype=3,alpha=0.3)
p <- p + facet_wrap(~month,ncol=4)
p <- p + geom_bar(data=Effplot,aes(x=as.numeric(Latitude), y=effort,fill="effort"),stat="identity", alpha=0.5 )
#p <- p + coord_flip()
p <- p + theme_bw()+theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.title.x= element_blank())
p <- p + scale_colour_grey() +scale_fill_grey(name=" ") 
p
ggsave(file="gtg_model_example.pdf")



tmp<-sim_gtg$PosX[sim_gtg$PosX[,1]>1188,]

meses<-c("Jan", "Feb", "Mar","Apr", "May", "Jun","Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

tmp<-data.frame(tmp)
tmp[,1]<-meses[tmp[,1]-(max(tmp[,1])-12)]

JAN<-tmp[tmp[,1]=="Jan",]
MAY<-tmp[tmp[,1]=="May",]


summary(sim_gtg$PosX[1,])