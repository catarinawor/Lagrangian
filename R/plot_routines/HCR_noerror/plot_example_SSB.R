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

plot_DynCatch <- function( M,  Msq, sv=F, nome="",nations=F, limite=.4)
{
	cat("plot_Bnat_TAC")

	gg_color_hue <- function(n) {
  		hues = seq(15, 375, length = n + 1)
  		hcl(h = hues, l = 65, c = 88)[1:n]
	}

	names(M[[1]][[1]])
	length(M[[1]][[1]]$ySB)

	n <- length( M )


	slope_hcr<-NULL
	intercept_hcr<-NULL
	SSB_total<-matrix(ncol=n, nrow=111)
	



	for(i in 1:n){

		slope_hcr[i]<-M[[i]][[1]]$slope_hcr
		intercept_hcr[i]<-M[[i]][[1]]$intercept_hcr
		SSB_total[,i]<-M[[i]][[1]]$ySB[50:160]
		
	}



	nsq<- length( Msq )

	slope_hcr_sq<-NULL
	intercept_hcr_sq<-NULL
	SSB_total_sq<-matrix(ncol=nsq, nrow=111)
	


	for(y in 1:nsq){

		slope_hcr_sq[y]<-999
		intercept_hcr_sq[y]<-999
		SSB_total_sq[,y]<-Msq[[y]][[1]]$ySB[50:160]

	}



	sp<- unique(slope_hcr)
	int<-unique(intercept_hcr)
	SSB_lo<-matrix(nrow=length(sp)*length(int)+1,ncol=111)
	SSB_hi<-matrix(nrow=length(sp)*length(int)+1,ncol=111)
	SSB_avg<-matrix(nrow=length(sp)*length(int)+1,ncol=111)

	
	intercept_dois<- NULL
	slope_dois<-NULL

	a <-1
	for(si in 1:length(sp)){
		for(ii in 1:length(int)){

			intercept_dois[a]<-int[ii]
			slope_dois[a]<-sp[si]

			

			ncol(Yield_nat2)

			tmp2<-SSB_total[,intersect(which(intercept_hcr==int[ii]),which(slope_hcr==sp[si]))]

			

			SSB_lo[a,]<-apply(tmp2,1,calc_quantile)[1,]
			SSB_hi[a,]<-apply(tmp2,1,calc_quantile)[5,]
			SSB_avg[a,]<-apply(tmp2,1,calc_quantile)[3,]


			a <- a+1
		}
	}

	intercept_dois[a]<-4010
	slope_dois[a]<-4010

	tmp2<-SSB_total_sq

	SSB_lo[a,]<-apply(tmp2,1,calc_quantile)[1,]
	SSB_hi[a,]<-apply(tmp2,1,calc_quantile)[5,]
	SSB_avg[a,]<-apply(tmp2,1,calc_quantile)[3,]
	
	

	lossb<-melt(SSB_lo)
	hissb<-melt(SSB_hi)
	avgssb<-melt(SSB_avg)
	

	ssb<- data.frame(year=(1966:2076)[avgssb$Var2], hcr=avgssb$Var1,
		 intercept=intercept_dois[avgssb$Var1],slope=slope_dois[avgssb$Var1],
		mean=avgssb$value,high=hissb$value,low=lossb$value)

	ssb$high[ssb$year<2016]=NA
	ssb$low[ssb$year<2016]=NA

	

	ssb_o<-ssb[(ssb$intercept==.1|ssb$intercept==4010)&(ssb$slope==.1|ssb$slope==.3|ssb$slope==4010),]
	
	pb<- ggplot(ssb_o)
	pb<- pb+geom_ribbon(aes(x=year,ymin=low, ymax=high, fill=as.factor(slope)),alpha=.3)
	pb<- pb+geom_line(aes(x=year,y=mean, color=as.factor(slope)))	
	pb<- pb +scale_color_manual(values = c(gg_color_hue(2),"black") ,name="exploitation rate",labels=c("0.1", "0.3", "40:10 rule"))
	pb<- pb +scale_fill_manual(values = c(gg_color_hue(2),"black"),name="exploitation rate",labels=c("0.1", "0.3", "40:10 rule"))
	pb <- pb + labs(x= "", y= "SSB (1000t) ")
	pb <- pb + theme_bw(16) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	pb

	
	

	if(sv==TRUE){
		setwd("/Users/catarinawor/Documents/Lagrangian/figures/HCR")
		ggsave(paste(nome,"ssb_trajectory.pdf",sep=""), plot=pb, width = 10, height = 7)
	}

		return(plun)

	
}

