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


	n <- length( M )


	slope_hcr<-NULL
	intercept_hcr<-NULL
	Yield_total<-matrix(ncol=n, nrow=111)
	Yield_nat1<-matrix(ncol=n, nrow=111)
	Yield_nat2<-matrix(ncol=n, nrow=111)



	for(i in 1:n){

		slope_hcr[i]<-M[[i]][[1]]$slope_hcr
		intercept_hcr[i]<-M[[i]][[1]]$intercept_hcr
		Yield_total[,i]<-(apply(M[[i]][[1]]$yYieldNat[50:160,],1,sum))
		Yield_nat1[,i]<-(M[[i]][[1]]$yYieldNat[50:160,1])
		Yield_nat2[,i]<-(M[[i]][[1]]$yYieldNat[50:160,2])
	}



	nsq<- length( Msq )

	slope_hcr_sq<-NULL
	intercept_hcr_sq<-NULL
	Yield_total_sq<-matrix(ncol=nsq, nrow=111)
	Yield_nat1_sq<-matrix(ncol=nsq, nrow=111)
	Yield_nat2_sq<-matrix(ncol=nsq, nrow=111)


	for(y in 1:nsq){

		slope_hcr_sq[y]<-999
		intercept_hcr_sq[y]<-999
		Yield_total_sq[,y]<-apply(Msq[[y]][[1]]$yYieldNat[50:160,],1,sum)
		Yield_nat1_sq[,y]<-Msq[[y]][[1]]$yYieldNat[50:160,1]
		Yield_nat2_sq[,y]<-Msq[[y]][[1]]$yYieldNat[50:160,2]
	}




	sp<- unique(slope_hcr)
	int<-unique(intercept_hcr)
	Yield_lo_nat2<-matrix(nrow=length(sp)*length(int)+1,ncol=111)
	Yield_hi_nat2<-matrix(nrow=length(sp)*length(int)+1,ncol=111)
	Yield_avg_nat2<-matrix(nrow=length(sp)*length(int)+1,ncol=111)

	Yield_lo_nat1<-matrix(nrow=length(sp)*length(int)+1,ncol=111)
	Yield_hi_nat1<-matrix(nrow=length(sp)*length(int)+1,ncol=111)
	Yield_avg_nat1<-matrix(nrow=length(sp)*length(int)+1,ncol=111)

	Yield_lo_tot<-matrix(nrow=length(sp)*length(int)+1,ncol=111)
	Yield_hi_tot<-matrix(nrow=length(sp)*length(int)+1,ncol=111)
	Yield_avg_tot<-matrix(nrow=length(sp)*length(int)+1,ncol=111)
	Yield_m_tot<-matrix(nrow=length(sp)*length(int)+1,ncol=111)

	intercept_dois<- NULL
	slope_dois<-NULL

	a <-1
	for(si in 1:length(sp)){
		for(ii in 1:length(int)){

			intercept_dois[a]<-int[ii]
			slope_dois[a]<-sp[si]

			

			ncol(Yield_nat2)

			tmp2<-Yield_nat2[,intersect(which(intercept_hcr==int[ii]),which(slope_hcr==sp[si]))]

			

			Yield_lo_nat2[a,]<-apply(tmp2,1,calc_quantile)[1,]
			Yield_hi_nat2[a,]<-apply(tmp2,1,calc_quantile)[5,]
			Yield_avg_nat2[a,]<-apply(tmp2,1,calc_quantile)[3,]

			tmp1<-Yield_nat1[,intersect(which(intercept_hcr==int[ii]),which(slope_hcr==sp[si]))]
			Yield_lo_nat1[a,]<-apply(tmp1,1,calc_quantile)[1,]
			Yield_hi_nat1[a,]<-apply(tmp1,1,calc_quantile)[5,]
			Yield_avg_nat1[a,]<-apply(tmp1,1,calc_quantile)[3,]
		
			tmp0<-Yield_total[,intersect(which(intercept_hcr==int[ii]),which(slope_hcr==sp[si]))]
			Yield_lo_tot[a,]<-apply(tmp0,1,calc_quantile)[1,]
			Yield_hi_tot[a,]<-apply(tmp0,1,calc_quantile)[5,]
			Yield_avg_tot[a,]<-apply(tmp0,1,calc_quantile)[3,]
			Yield_m_tot[a,]<-apply(tmp0,1,mean)

			a <- a+1
		}
	}

	intercept_dois[a]<-4010
	slope_dois[a]<-4010

	tmp2<-Yield_nat2_sq

	Yield_lo_nat2[a,]<-apply(tmp2,1,calc_quantile)[1,]
	Yield_hi_nat2[a,]<-apply(tmp2,1,calc_quantile)[5,]
	Yield_avg_nat2[a,]<-apply(tmp2,1,calc_quantile)[3,]
	
	tmp1<-Yield_nat1_sq
	Yield_lo_nat1[a,]<-apply(tmp1,1,calc_quantile)[1,]
	Yield_hi_nat1[a,]<-apply(tmp1,1,calc_quantile)[5,]
	Yield_avg_nat1[a,]<-apply(tmp1,1,calc_quantile)[3,]
		
	tmp0<-Yield_total_sq
	Yield_lo_tot[a,]<-apply(tmp0,1,calc_quantile)[1,]
	Yield_hi_tot[a,]<-apply(tmp0,1,calc_quantile)[5,]
	Yield_avg_tot[a,]<-apply(tmp0,1,calc_quantile)[3,]
	Yield_m_tot[a,]<-apply(tmp0,1,mean)

	lonat2<-melt(Yield_lo_nat2)
	hinat2<-melt(Yield_hi_nat2)
	avgnat2<-melt(Yield_avg_nat2)
	

	nat2<- data.frame(year=(1966:2076)[avgnat2$Var2], hcr=avgnat2$Var1,
		 intercept=intercept_dois[avgnat2$Var1],slope=slope_dois[avgnat2$Var1],
		mean=avgnat2$value,high=hinat2$value,low=lonat2$value)

	nat2$high[nat2$year<2016]=NA
	nat2$low[nat2$year<2016]=NA

	nat2_o<-nat2[(nat2$intercept==.1|nat2$intercept==4010)&(nat2$slope==.1|nat2$slope==.5|nat2$intercept==4010),]
	
	p<- ggplot(nat2_o)
	p<- p+geom_ribbon(aes(x=year,ymin=low, ymax=high, fill=as.factor(slope)),alpha=.3)
	p<- p+geom_line(aes(x=year,y=mean, color=as.factor(slope)))	
	p<- p +scale_color_manual(values = c(gg_color_hue(2),"black") ,name="exploitation rate",labels=c("0.1", "0.5", "40:10 rule"))
	p<- p +scale_fill_manual(values = c(gg_color_hue(2),"black"),name="exploitation rate",labels=c("0.1", "0.5", "40:10 rule"))
	p <- p + labs(x= "", y= " Yield (1000t) ")
	p <- p + theme_bw(16) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	p

	
	lonat1<-melt(Yield_lo_nat1)
	hinat1<-melt(Yield_hi_nat1)
	avgnat1<-melt(Yield_avg_nat1)
	

	nat1<- data.frame(year=(1966:2076)[avgnat1$Var2], hcr=avgnat1$Var1,
		 intercept=intercept_dois[avgnat1$Var1],slope=slope_dois[avgnat1$Var1],
		mean=avgnat1$value,high=hinat1$value,low=lonat1$value)

	nat1$high[nat1$year<2016]=NA
	nat1$low[nat1$year<2016]=NA

	nat1_o<-nat1[(nat1$intercept==.1|nat1$intercept==4010)&(nat1$slope==.1|nat1$slope==.5|nat1$intercept==4010),]
	
	p1<- ggplot(nat1_o)
	p1<- p1+geom_ribbon(aes(x=year,ymin=low, ymax=high, fill=as.factor(slope)),alpha=.3)
	p1<- p1+geom_line(aes(x=year,y=mean, color=as.factor(slope)))	
	p1<- p1 +scale_color_manual(values = c(gg_color_hue(2),"black") ,name="exploitation rate",labels=c("0.1", "0.5", "40:10 rule"))
	p1<- p1 +scale_fill_manual(values = c(gg_color_hue(2),"black"),name="exploitation rate",labels=c("0.1", "0.5", "40:10 rule"))
	p1 <- p1 + labs(x= "Year", y= " Yield (1000t) ")
	p1 <- p1 + theme_bw(16) +theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

	p1


	plun<-plot_grid(p1,p,nrow=2,labels=c("Nation 1","Nation 2"),hjust = -1.5,vjust = 3)
		print(plun)



	if(sv==TRUE){
		setwd("/Users/catarinawor/Documents/Lagrangian/figures/HCR")
		ggsave(paste(nome,"catch_trajectory.pdf",sep=""), plot=plun, width = 10, height = 7)
	}

		return(plun)

	
}

