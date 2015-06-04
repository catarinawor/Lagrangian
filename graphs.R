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

setwd("/Users/catarinawor/Documents/Lagrangian/anime")

library("animation")

#saveLatex(
#if you don't have latex installed you need to install it or use other function such as saveGIF

#for(i in 1:length(ntsp))
#{
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
	y<- y+ geom_text(data = NULL, x = 7 , y = .9, label = meses[indmonth[i]] )
print(y)
#}
#,
#pdflatex = "/usr/texbin/pdflatex")




#====================================================================================================================

#plot effort and abundance
effagesim<-melt(matrix(sim$Effage,ncol=(sim$nage-sim$sage+1),dimnames=list(a=1:((sim$nyr-sim$syr+1)*(sim$nmon-est$smon+1)),b=sim$sage:sim$nage)))
effageest<-melt(matrix(est$Effage,ncol=(est$nage-est$sage+1),dimnames = list(a=1:((est$nyr-est$syr+1)*(est$nmon-est$smon+1)),b=est$sage:est$nage)))

names(effagesim) <- c("time", "age", "effort") 
names(effageest) <- c("time", "age", "effort") 

    
effagesim=effagesim[order(effagesim$time),]

vo<- ggplot(effagesim[29:42,], aes(age, effort)) + theme_bw()
vo <-vo+ geom_line(size=2)
vo

v <- ggplot(effagesim[22:28,], aes(age, z = effort)) + theme_bw()
v <- v +  stat_contour(aes(colour=..level..,fill=..level..)) 
v <- v + stat_contour(geom="polygon", aes(fill=..level..))

v2 <- ggplot(effageest, aes(time, age, z = effort)) + theme_bw()
v2 <- v2 +  stat_contour(aes(colour=..level..,fill=..level..)) 
v2 <- v2 + stat_contour(geom="polygon", aes(fill=..level..))

multiplot(v,v2,cols=2)



df<-data.frame(time=rep(1:length(indyr[indmonth==12]),3),Abundance=c(SB[indmonth==12],apply(VulB,1,sum)[indmonth==12],apply(Nage,1,sum)[indmonth==12]),type=rep(c("SB","VB","N"),each=ntstp/12))

p <- ggplot(df, aes(x=time,y=Abundance,group=type)) + theme_bw()
p <- p+ geom_line(aes(colour=type),size=2)
p



#plot biomass by area

VBareaplot<-melt(VBarea) 
names(VBareaplot) <- c("time", "area", "VB")
p1 <- ggplot(VBareaplot, aes(time, area, z = VB)) + theme_bw()
p1 <- p1 +  stat_contour(aes(colour=..level..,fill=..level..))+ scale_fill_gradient(low = "olivedrab3", high = "magenta3")
p1 <- p1 + stat_contour(geom="polygon", aes(fill=..level..,colour=..level..))+ scale_colour_gradient(low = "olivedrab3", high = "magenta3")
#p1 <- p1 + scale_colour_gradient(low = "olivedrab3", high = "magenta3")
p1

multiplot(p1,p,cols=1)


VBareaplot<-melt(VBarea[325:432,]) 
names(VBareaplot) <- c("time", "area", "VB")
p1 <- ggplot(VBareaplot, aes(time, area, z = VB)) + theme_bw()
p1 <- p1 +  stat_contour(aes(colour=..level..,fill=..level..))+ scale_fill_gradient(low = "olivedrab3", high = "magenta3")
p1 <- p1 + stat_contour(geom="polygon", aes(fill=..level..,colour=..level..))+ scale_colour_gradient(low = "olivedrab3", high = "magenta3")

#p1 <- p1 + scale_colour_gradient(low = "olivedrab3", high = "magenta3")
p1

#======================================================================================================
#animation stuff

# How would each frame of the animations look like?

VBareaplot<-melt(VBarea) 
names(VBareaplot) <- c("time", "area", "VB")
x<-as.factor(rep(1,ntstp*nareas))
nation<- rep("Canada",ntstp*nareas)
nation[VBareaplot$area<nationareas[1]]<-"U.S.A"
Month<-(as.factor(indmonth))
lat<-VBareaplot$area+29
lon<-rep(-133,length(lat))
VBareaplot<-cbind(indyr,indmonth,lat,lon,VBareaplot,x,nation,Month)

