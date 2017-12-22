#==============================================================
# Plot the ratio between the biomass and TAc for each nations
#Author: Catarina Wor
#Date: Sep 14 2017

#==============================================================




#M<-SIMSdat

require(reshape2)
require(tidyr)
require(ggplot2)
library(scales)
library(RColorBrewer)
library(cowplot)


source("/Users/catarinawor/Documents/Lagrangian/R/plot_routines/HCR_noerror/calc_quantile.R")

plot_DynVB <- function( M,  Msq, sv=F, nome="",nations=F, limite=.4)
{
	cat("plot_Bnat_TAC")

	gg_color_hue <- function(n) {
  		hues = seq(15, 375, length = n + 1)
  		hcl(h = hues, l = 65, c = 88)[1:n]
	}


	n <- length( M )

	names(M[[1]][[1]])

	dim(M[[1]][[1]]$tVBarea)
	length(M[[1]][[1]]$indmonth)

	slope_hcr<-NULL
	intercept_hcr<-NULL

	melt(M[[1]][[1]]$tVBarea[M[[1]][[1]]$indmonth==7,])


	VB_total<-NULL
	


	for(i in 1:n){

		slope_hcr[i]<-M[[i]][[1]]$slope_hcr
		intercept_hcr[i]<-M[[i]][[1]]$intercept_hcr

		tmp<-melt(M[[i]][[1]]$tVBarea[M[[i]][[1]]$indmonth==7,])
		colnames(tmp)<-c("year","area","VB")
		tmp2<-cbind(tmp,slope=slope_hcr[i],intercept=intercept_hcr[i])		
		
		VB_total<-rbind(VB_total,tmp2)

	}

	VB_total_selec<-VB_total[VB_total$intercept<.2&
	VB_total$slope!=0&VB_total$slope!=.2&VB_total$slope!=.3&
	VB_total$slope!=.4&VB_total$slope<.6&VB_total$year>49,]
	summary(VB_total_selec)

	nsq<- length( Msq )

	slope_hcr_sq<-NULL
	intercept_hcr_sq<-NULL
	VB_total_sq<-NULL
	


	for(y in 1:nsq){

		slope_hcr_sq[y]<-999
		intercept_hcr_sq[y]<-999

		tmp<-melt(Msq[[y]][[1]]$tVBarea[Msq[[y]][[1]]$indmonth==7,])
		colnames(tmp)<-c("year","area","VB")
		tmp2<-cbind(tmp,slope=slope_hcr_sq[y],intercept=intercept_hcr_sq[y])		
		
		VB_total_sq<-rbind(VB_total_sq,tmp2)
	}



	VB<-rbind(VB_total_selec,VB_total_sq)
	names(VB)

	p<- ggplot(VB)
	p<- p+geom_boxplot(aes(x=area,y=VB, fill=as.factor(slope)),alpha=.3, outlier.shape=NA)
	p<- p +scale_color_manual(values = c(gg_color_hue(2),"black") ,name="exploitation rate",labels=c("0.1", "0.5", "40:10 rule"))
	p<- p +scale_fill_manual(values = c(gg_color_hue(2),"black"),name="exploitation rate",labels=c("0.1", "0.5", "40:10 rule"))
	p <- p + labs(x= "Area", y= " Biomass (1000 t) ")
	p <- p + theme_bw(16) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	p <- p + scale_x_discrete(labels=as.character(seq(30,60)))
	p <- p + geom_vline(xintercept=19.5)
	p <- p + coord_cartesian(ylim=c(0,220))
	p <- p + annotate("text",x=21,y=200,label="Canada", fontface =2,size=7)
	p <- p + annotate("text",x=2,y=200,label="U.S.A.", fontface =2,size=7)
			
			
	p




	if(sv==TRUE){
		setwd("/Users/catarinawor/Documents/Lagrangian/report/HCR_coarse")
		ggsave(paste(nome,"catch_trajectory.pdf",sep=""), plot=p, width = 12, height = 6)
	}

	return(plun)
	
	
}

