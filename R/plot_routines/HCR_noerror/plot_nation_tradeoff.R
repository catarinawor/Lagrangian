#==============================================================
# Read in the results from the HCR search - no assessment error
#Author: Catarina Wor
#Date: Sep 14 2017

#==============================================================





require(reshape2)
require(tidyr)
require(ggplot2)
library(scales)
library(RColorBrewer)
library(cowplot)


plot_tradeoff <- function( M1 ,M2,M3, Msq1,Msq2,Msq3, sv=F, nome="",nations=F)
{
	cat("plot_tradeoff")

	n <- c(length( M1 ),length( M2 ),length( M3 ))

	mm<-list(M1 ,M2,M3)
	mmsq<-list(Msq1,Msq2,Msq3)

	slope_hcr<-NULL
	intercept_hcr<-NULL
	utility_total<-NULL
	utility_total2<-NULL
	utility_nat1<-NULL
	utility_nat2<-NULL
	yield_total<-NULL
	yield_nat1<-NULL
	yield_nat2<-NULL
	B40_total<-NULL
	scenario<-NULL

	i=1
	for(d in 1:3){

		M<-mm[[d]]

		for(ii in 1:n[d]){

			iniyr<-M[[ii]][[1]]$"nyr"+1
			fimyr<-M[[ii]][[1]]$"proj_yr"
		
			B40_total[i]<-sum((M[[ii]][[1]]$ytB/M[[ii]][[1]]$Bo)[iniyr:fimyr]<0.4)/length(iniyr:fimyr)

			slope_hcr[i]<-M[[ii]][[1]]$slope_hcr
			intercept_hcr[i]<-M[[ii]][[1]]$intercept_hcr
			utility_total[i]<-mean(log(apply(M[[ii]][[1]]$yYieldNat[101:160,],1,sum)+1))
			utility_total2[i]<-mean(log(M[[ii]][[1]]$yYieldNat[101:160,1]+1))*sum(log(M[[ii]][[1]]$yYieldNat[101:160,2]+1))
			utility_nat1[i]<-mean(log(M[[ii]][[1]]$yYieldNat[101:160,1]+1))
			utility_nat2[i]<-mean(log(M[[ii]][[1]]$yYieldNat[101:160,2]+1))
			yield_total[i]<-mean(apply(M[[ii]][[1]]$yYieldNat[101:160,],1,sum))
			yield_nat1[i]<-mean((M[[ii]][[1]]$yYieldNat[101:160,1]))
			yield_nat2[i]<-mean((M[[ii]][[1]]$yYieldNat[101:160,2]))

			scenario[i]<-d
			i=i+1
		}

	}

	length(slope_hcr)


	slope_hcr_sq<-NULL
	intercept_hcr_sq<-NULL
	utility_total_sq<-NULL
	utility_total2_sq<-NULL
	utility_nat1_sq<-NULL
	utility_nat2_sq<-NULL
	yield_total_sq<-NULL
	yield_nat1_sq<-NULL
	yield_nat2_sq<-NULL
	B40_total_sq<-NULL
	scenario_sq<-NULL

	nsq<- c(length( Msq1 ),length( Msq2 ),length( Msq3 ))
	y<-1

	for(dd in 1:3){

		Msq<-mmsq[[dd]]

		for(yy in 1:nsq[[dd]]){

			slope_hcr_sq[y]<-999
			intercept_hcr_sq[y]<-999
			utility_total_sq[y]<-sum(log(apply(Msq[[yy]][[1]]$yYieldNat[101:160,],1,sum)+1))
			utility_total2_sq[y]<-mean(log(Msq[[yy]][[1]]$yYieldNat[101:160,1]+1))*sum(log(Msq[[yy]][[1]]$yYieldNat[101:160,2]+1))
			utility_nat1_sq[y]<-mean(log(Msq[[yy]][[1]]$yYieldNat[101:160,1]+1))
			utility_nat2_sq[y]<-mean(log(Msq[[yy]][[1]]$yYieldNat[101:160,2]+1))
			yield_total_sq[y]<-mean(apply(Msq[[yy]][[1]]$yYieldNat[101:160,],1,sum))
			yield_nat1_sq[y]<-mean((Msq[[yy]][[1]]$yYieldNat[101:160,1]))
			yield_nat2_sq[y]<-mean((Msq[[yy]][[1]]$yYieldNat[101:160,2]))
			B40_total_sq[y]<-sum((Msq[[yy]][[1]]$ytB/Msq[[yy]][[1]]$Bo)[iniyr:fimyr]<0.4)/length(iniyr:fimyr)
			scenario_sq[y]<-dd
			y=y+1
		}

	}

	dflu<-data.frame(slope_hcr=slope_hcr,intercept_hcr=intercept_hcr,
		utility_total2=utility_total2,utility_total=utility_total,utility_nat1=utility_nat1,
		utility_nat2=utility_nat2,yield_total=yield_total,yield_nat1=yield_nat1,
		yield_nat2=yield_nat2,	b40=B40_total, scenario=scenario)

	dflusq<-data.frame(slope_hcr=slope_hcr_sq,intercept_hcr=intercept_hcr_sq,
		utility_total2=utility_total2_sq,utility_total=utility_total_sq,utility_nat1=utility_nat1_sq,
		utility_nat2=utility_nat2_sq,yield_total=yield_total_sq,yield_nat1=yield_nat1_sq,
		yield_nat2=yield_nat2_sq,b40=B40_total_sq,scenario=scenario_sq)

	dflu_plot<-aggregate(dflu,list(dflu$slope_hcr,dflu$intercept_hcr, dflu$scenario), mean)

	dflu_plotsq<-aggregate(dflusq,list(dflusq$slope_hcr,dflusq$intercept_hcr, dflusq$scenario), mean)




	summary(dflu_plot)
	head(dflu_plotsq)

	dflu_selec<-dflu_plot

	dfa<-rbind(dflu_selec,dflu_plotsq)

	dfa$slope<-as.character(dfa$slope_hcr)
	dfa$slope[dfa$slope=="999"]<-"40:10"

	dfa$intercept<-as.character(dfa$intercept_hcr)
	dfa$intercept[dfa$intercept=="999"]<-"40:10"

	dfa$ut2tf<- -999
	dfa$ut2tf[which(dfa$slope=="40:10")]<-dfa$utility_total2[which(dfa$slope=="40:10")]

	dfa$nt1tf<- -999
	dfa$nt1tf[which(dfa$slope=="40:10")]<-dfa$utility_nat1[which(dfa$slope=="40:10")]

	dfa$nt2tf<- -999
	dfa$nt2tf[which(dfa$slope=="40:10")]<-dfa$utility_nat2[which(dfa$slope=="40:10")]

	#dfa$scn<-"Base case movement"
	#dfa$scn[dfa$scenario==2]<-"early movement"
	#dfa$scn[dfa$scenario==3]<-"late movement"

	dfa$scn<-"no strong recruitment"
	dfa$scn[dfa$scenario==2]<-"one strong recruitment"
	dfa$scn[dfa$scenario==3]<-"two strong recruitments"


	dfa$b40tf<- -999
	dfa$b40tf[dfa$slope=="40:10"]<-dfa$b40[which(dfa$slope=="40:10")]

summary(dfa)

	dfc<-dfa[dfa$intercept_hcr<0.4|dfa$intercept_hcr>2,] 
	dfcr<-dfc[dfc$slope_hcr<0.5|dfc$slope_hcr>2,] 
	summary(dfcr)

	dfcr$slope<-as.factor(dfcr$slope)
	dfcr$intercept<-as.factor(dfcr$intercept)
	levels(dfcr$intercept)

	dfcr$utility_nat2_sc[dfcr$scenario==1]<-dfcr$utility_nat2[dfcr$scenario==1]/max(dfcr$utility_nat2[dfcr$scenario==1])
	dfcr$utility_nat2_sc[dfcr$scenario==2]<-dfcr$utility_nat2[dfcr$scenario==2]/max(dfcr$utility_nat2[dfcr$scenario==2])
	dfcr$utility_nat2_sc[dfcr$scenario==3]<-dfcr$utility_nat2[dfcr$scenario==3]/max(dfcr$utility_nat2[dfcr$scenario==3])
	
	dfcr$utility_nat1_sc[dfcr$scenario==1]<-dfcr$utility_nat1[dfcr$scenario==1]/max(dfcr$utility_nat1[dfcr$scenario==1])
	dfcr$utility_nat1_sc[dfcr$scenario==2]<-dfcr$utility_nat1[dfcr$scenario==2]/max(dfcr$utility_nat1[dfcr$scenario==2])
	dfcr$utility_nat1_sc[dfcr$scenario==3]<-dfcr$utility_nat1[dfcr$scenario==3]/max(dfcr$utility_nat1[dfcr$scenario==3])
	

	dfcr$yield_nat2_sc[dfcr$scenario==1]<-dfcr$yield_nat2[dfcr$scenario==1]/max(dfcr$yield_nat2[dfcr$scenario==1])
	dfcr$yield_nat2_sc[dfcr$scenario==2]<-dfcr$yield_nat2[dfcr$scenario==2]/max(dfcr$yield_nat2[dfcr$scenario==2])
	dfcr$yield_nat2_sc[dfcr$scenario==3]<-dfcr$yield_nat2[dfcr$scenario==3]/max(dfcr$yield_nat2[dfcr$scenario==3])
	
	dfcr$yield_nat1_sc[dfcr$scenario==1]<-dfcr$yield_nat1[dfcr$scenario==1]/max(dfcr$yield_nat1[dfcr$scenario==1])
	dfcr$yield_nat1_sc[dfcr$scenario==2]<-dfcr$yield_nat1[dfcr$scenario==2]/max(dfcr$yield_nat1[dfcr$scenario==2])
	dfcr$yield_nat1_sc[dfcr$scenario==3]<-dfcr$yield_nat1[dfcr$scenario==3]/max(dfcr$yield_nat1[dfcr$scenario==3])
	

	dfcr<-dfcr[order(dfcr$intercept),]

	dft<-dfcr[dfcr$slope_hcr==999,]

	summary(dfcr)

	#presentation version
	p<-ggplot(dfcr)
	p<-p+geom_path(aes(y=utility_nat2_sc,x=utility_nat1_sc,color=slope),size=1.5)
	p<-p+geom_point(aes(y=utility_nat2_sc,x=utility_nat1_sc,color=slope, shape=intercept),size=8,alpha=.6)
	p<-p+ geom_text(data=dft,aes(y=dft$utility_nat2_sc,x=dft$utility_nat1_sc,label="40:10"),fontface = "bold",colour="black")				
	p<-p + facet_wrap(~scn, scales="free")
	p<-p+ ylab("average log yield - Canada") + xlab("average log yield - U.S.A.")
	p <- p  + theme_bw(18)+ theme(legend.position = "none")
	p <- p + guides(color=guide_legend("harvest rate"))
	p <- p + theme(axis.text = element_text(face="bold", size=18),
  			axis.title = element_text(face="bold", size=18),
  			strip.text = element_text(face="bold", size=18))
	print(p)

	

	p1 <- ggplot(dfcr)
	p1 <- p1+geom_path(aes(y=yield_nat2_sc,x=yield_nat1_sc,color=slope),size=1.5)
	p1 <- p1+geom_point(aes(y=yield_nat2_sc,x=yield_nat1_sc,color=slope, shape=intercept),size=8,alpha=.6)
	p1 <- p1+ geom_text(data=dft,aes(y=dft$yield_nat2_sc,x=dft$yield_nat1_sc,label="40:10"),fontface = "bold",colour="black")				
	
	p1 <- p1+facet_wrap(~scn, scales="free")
	p1 <- p1+ ylab(" average yield - Canada ") + xlab(" average yield - U.S.A.")
	p1 <- p1  + theme_bw(18)+ theme(legend.position = "bottom")
	p1 <- p1 + guides(color=guide_legend("harvest rate"), shape=guide_legend("biomass threshold"),size=guide_legend("none"))
	p1 <- p1 +  theme(axis.text = element_text(face="bold", size=18),
  			axis.title = element_text(face="bold", size=18),
  			strip.text = element_text(face="bold", size=18))
	#p1 <-p1+ geom_text(aes(y=nt2tf,x=nt1tf,label="10:40"),fontface = "bold",colour="black", position = position_nudge(y = -0.02))				
	#, size=guide_legend("biomass threshold")
	p1

	pto<-plot_grid(p,p1,nrow=2,rel_heights  = c(1, 1.3))
	pto



	p<-ggplot(dfcr)
	p<-p+geom_path(aes(y=utility_nat2_sc,x=utility_nat1_sc,color=slope))
	p<-p+geom_point(aes(y=utility_nat2_sc,x=utility_nat1_sc,color=slope, shape=intercept),size=4,alpha=.6)
	p<-p+ geom_text(data=dft,aes(y=dft$utility_nat2_sc,x=dft$utility_nat1_sc,label="40:10"),fontface = "bold",colour="black")				
	p<-p + facet_wrap(~scn, scales="free")
	p<-p+ ylab("average log utility - Nation 2") + xlab("average log utility - Nation 1")
	p <- p  + theme_bw(16)+ theme(legend.position = "none")
	p <- p + guides(color=guide_legend("harvest rate"))
	p
	

	p1 <- ggplot(dfcr)
	p1 <- p1+geom_path(aes(y=yield_nat2_sc,x=yield_nat1_sc,color=slope))
	p1 <- p1+geom_point(aes(y=yield_nat2_sc,x=yield_nat1_sc,color=slope, shape=intercept),size=4,alpha=.6)
	#p1 <- p1+ geom_text(data=dft,aes(y=dft$yield_nat2_sc,x=dft$yield_nat1_sc,label="40:10"),fontface = "bold",colour="black")				
	
	p1 <- p1+facet_wrap(~scn, scales="free")
	p1 <- p1+ ylab(" average yield - Nation 2") + xlab(" average yield - Nation 1")
	p1 <- p1  + theme_bw(16)+ theme(legend.position = "bottom")
	p1 <- p1 + guides(color=guide_legend("harvest rate"), shape=guide_legend("biomass threshold"),size=guide_legend("none"))
	p1
	#p<-p+ geom_text(aes(y=nt2tf,x=nt1tf,label="10:40"),fontface = "bold",colour="black", position = position_nudge(y = -0.02))				
	#, size=guide_legend("biomass threshold")

	pto<-plot_grid(p,p1,nrow=2,rel_heights  = c(1, 1.3))
	pto
	

	if(sv==TRUE){
			#setwd("/Users/catarinawor/Documents/Lagrangian/report/manuscript/ICES")
			#ggsave(paste(nome,"tradeoff_nations_mov_scn.pdf",sep=""), plot=pto, width = 11, height = 9)
			setwd("/Users/catarinawor/Documents/")
			ggsave(paste("tradeoff.pdf",sep=""), plot=pto, width = 13, height = 10)
		
		}

		return(pto)

#dflu<-dflu[dflu$intercept_hcr<0.6,]

	#summary(dflu)


	
	
}


plot_tradeoff_later <- function( M1 ,M2,M3, Msq1,Msq2,Msq3, sv=F, nome="",nations=F)
{
	cat("plot_tradeoff")

	n <- c(length( M1 ),length( M2 ),length( M3 ))

	mm<-list(M1 ,M2,M3)
	mmsq<-list(Msq1,Msq2,Msq3)

	slope_hcr<-NULL
	intercept_hcr<-NULL
	utility_total<-NULL
	utility_total2<-NULL
	utility_nat1<-NULL
	utility_nat2<-NULL
	yield_total<-NULL
	yield_nat1<-NULL
	yield_nat2<-NULL
	B40_total<-NULL
	scenario<-NULL

	i=1
	for(d in 1:3){

		M<-mm[[d]]

		for(ii in 1:n[d]){

			iniyr<-M[[ii]][[1]]$"nyr"+41
			fimyr<-M[[ii]][[1]]$"proj_yr"
		
			B40_total[i]<-sum((M[[ii]][[1]]$ytB/M[[ii]][[1]]$Bo)[iniyr:fimyr]<0.4)/length(iniyr:fimyr)

			slope_hcr[i]<-M[[ii]][[1]]$slope_hcr
			intercept_hcr[i]<-M[[ii]][[1]]$intercept_hcr
			utility_total[i]<-mean(log(apply(M[[ii]][[1]]$yYieldNat[iniyr:160,],1,sum)+1))
			utility_total2[i]<-mean(log(M[[ii]][[1]]$yYieldNat[iniyr:160,1]+1))*sum(log(M[[ii]][[1]]$yYieldNat[101:160,2]+1))
			utility_nat1[i]<-mean(log(M[[ii]][[1]]$yYieldNat[iniyr:160,1]+1))
			utility_nat2[i]<-mean(log(M[[ii]][[1]]$yYieldNat[iniyr:160,2]+1))
			yield_total[i]<-mean(apply(M[[ii]][[1]]$yYieldNat[iniyr:160,],1,sum))
			yield_nat1[i]<-mean((M[[ii]][[1]]$yYieldNat[iniyr:160,1]))
			yield_nat2[i]<-mean((M[[ii]][[1]]$yYieldNat[iniyr:160,2]))

			scenario[i]<-d
			i=i+1
		}

	}

	length(slope_hcr)


	slope_hcr_sq<-NULL
	intercept_hcr_sq<-NULL
	utility_total_sq<-NULL
	utility_total2_sq<-NULL
	utility_nat1_sq<-NULL
	utility_nat2_sq<-NULL
	yield_total_sq<-NULL
	yield_nat1_sq<-NULL
	yield_nat2_sq<-NULL
	B40_total_sq<-NULL
	scenario_sq<-NULL

	nsq<- c(length( Msq1 ),length( Msq2 ),length( Msq3 ))
	y<-1

	for(dd in 1:3){

		Msq<-mmsq[[dd]]

		for(yy in 1:nsq[[dd]]){

			slope_hcr_sq[y]<-999
			intercept_hcr_sq[y]<-999
			utility_total_sq[y]<-sum(log(apply(Msq[[yy]][[1]]$yYieldNat[101:160,],1,sum)+1))
			utility_total2_sq[y]<-mean(log(Msq[[yy]][[1]]$yYieldNat[101:160,1]+1))*sum(log(Msq[[yy]][[1]]$yYieldNat[101:160,2]+1))
			utility_nat1_sq[y]<-mean(log(Msq[[yy]][[1]]$yYieldNat[101:160,1]+1))
			utility_nat2_sq[y]<-mean(log(Msq[[yy]][[1]]$yYieldNat[101:160,2]+1))
			yield_total_sq[y]<-mean(apply(Msq[[yy]][[1]]$yYieldNat[101:160,],1,sum))
			yield_nat1_sq[y]<-mean((Msq[[yy]][[1]]$yYieldNat[101:160,1]))
			yield_nat2_sq[y]<-mean((Msq[[yy]][[1]]$yYieldNat[101:160,2]))
			B40_total_sq[y]<-sum((Msq[[yy]][[1]]$ytB/Msq[[yy]][[1]]$Bo)[iniyr:fimyr]<0.4)/length(iniyr:fimyr)
			scenario_sq[y]<-dd
			y=y+1
		}

	}

	dflu<-data.frame(slope_hcr=slope_hcr,intercept_hcr=intercept_hcr,
		utility_total2=utility_total2,utility_total=utility_total,utility_nat1=utility_nat1,
		utility_nat2=utility_nat2,yield_total=yield_total,yield_nat1=yield_nat1,
		yield_nat2=yield_nat2,	b40=B40_total, scenario=scenario)

	dflusq<-data.frame(slope_hcr=slope_hcr_sq,intercept_hcr=intercept_hcr_sq,
		utility_total2=utility_total2_sq,utility_total=utility_total_sq,utility_nat1=utility_nat1_sq,
		utility_nat2=utility_nat2_sq,yield_total=yield_total_sq,yield_nat1=yield_nat1_sq,
		yield_nat2=yield_nat2_sq,b40=B40_total_sq,scenario=scenario_sq)

	dflu_plot<-aggregate(dflu,list(dflu$slope_hcr,dflu$intercept_hcr, dflu$scenario), mean)

	dflu_plotsq<-aggregate(dflusq,list(dflusq$slope_hcr,dflusq$intercept_hcr, dflusq$scenario), mean)




	summary(dflu_plot)
	head(dflu_plotsq)

	dflu_selec<-dflu_plot

	dfa<-rbind(dflu_selec,dflu_plotsq)

	dfa$slope<-as.character(dfa$slope_hcr)
	dfa$slope[dfa$slope=="999"]<-"40:10"

	dfa$intercept<-as.character(dfa$intercept_hcr)
	dfa$intercept[dfa$intercept=="999"]<-"40:10"

	dfa$ut2tf<- -999
	dfa$ut2tf[which(dfa$slope=="40:10")]<-dfa$utility_total2[which(dfa$slope=="40:10")]

	dfa$nt1tf<- -999
	dfa$nt1tf[which(dfa$slope=="40:10")]<-dfa$utility_nat1[which(dfa$slope=="40:10")]

	dfa$nt2tf<- -999
	dfa$nt2tf[which(dfa$slope=="40:10")]<-dfa$utility_nat2[which(dfa$slope=="40:10")]

	dfa$scn<-"Base case movement"
	dfa$scn[dfa$scenario==2]<-"early movement"
	dfa$scn[dfa$scenario==3]<-"late movement"

	dfa$b40tf<- -999
	dfa$b40tf[dfa$slope=="40:10"]<-dfa$b40[which(dfa$slope=="40:10")]

summary(dfa)

	dfc<-dfa[dfa$intercept_hcr<0.4|dfa$intercept_hcr>2,] 
	dfcr<-dfc
	summary(dfcr)

	dfcr$slope<-as.factor(dfcr$slope)
	dfcr$intercept<-as.factor(dfcr$intercept)
	levels(dfcr$intercept)

	dfcr$utility_nat2_sc[dfcr$scenario==1]<-dfcr$utility_nat2[dfcr$scenario==1]/max(dfcr$utility_nat2[dfcr$scenario==1])
	dfcr$utility_nat2_sc[dfcr$scenario==2]<-dfcr$utility_nat2[dfcr$scenario==2]/max(dfcr$utility_nat2[dfcr$scenario==2])
	dfcr$utility_nat2_sc[dfcr$scenario==3]<-dfcr$utility_nat2[dfcr$scenario==3]/max(dfcr$utility_nat2[dfcr$scenario==3])
	
	dfcr$utility_nat1_sc[dfcr$scenario==1]<-dfcr$utility_nat1[dfcr$scenario==1]/max(dfcr$utility_nat1[dfcr$scenario==1])
	dfcr$utility_nat1_sc[dfcr$scenario==2]<-dfcr$utility_nat1[dfcr$scenario==2]/max(dfcr$utility_nat1[dfcr$scenario==2])
	dfcr$utility_nat1_sc[dfcr$scenario==3]<-dfcr$utility_nat1[dfcr$scenario==3]/max(dfcr$utility_nat1[dfcr$scenario==3])
	

	dfcr$yield_nat2_sc[dfcr$scenario==1]<-dfcr$yield_nat2[dfcr$scenario==1]/max(dfcr$yield_nat2[dfcr$scenario==1])
	dfcr$yield_nat2_sc[dfcr$scenario==2]<-dfcr$yield_nat2[dfcr$scenario==2]/max(dfcr$yield_nat2[dfcr$scenario==2])
	dfcr$yield_nat2_sc[dfcr$scenario==3]<-dfcr$yield_nat2[dfcr$scenario==3]/max(dfcr$yield_nat2[dfcr$scenario==3])
	
	dfcr$yield_nat1_sc[dfcr$scenario==1]<-dfcr$yield_nat1[dfcr$scenario==1]/max(dfcr$yield_nat1[dfcr$scenario==1])
	dfcr$yield_nat1_sc[dfcr$scenario==2]<-dfcr$yield_nat1[dfcr$scenario==2]/max(dfcr$yield_nat1[dfcr$scenario==2])
	dfcr$yield_nat1_sc[dfcr$scenario==3]<-dfcr$yield_nat1[dfcr$scenario==3]/max(dfcr$yield_nat1[dfcr$scenario==3])
	

	dfcr<-dfcr[order(dfcr$intercept),]

	dft<-dfcr[dfcr$slope_hcr==999,]

	p<-ggplot(dfcr)
	p<-p+geom_path(aes(y=utility_nat2_sc,x=utility_nat1_sc,color=slope))
	p<-p+geom_point(aes(y=utility_nat2_sc,x=utility_nat1_sc,color=slope, shape=intercept),size=4,alpha=.6)
	p<-p+ geom_text(data=dft,aes(y=dft$utility_nat2_sc,x=dft$utility_nat1_sc,label="40:10"),fontface = "bold",colour="black")				
	p<-p + facet_wrap(~scn, scales="free")
	p<-p+ ylab("average log utility - Nation 2") + xlab("average log utility - Nation 1")
	p <- p  + theme_bw(16)+ theme(legend.position = "none")
	p <- p + guides(color=guide_legend("harvest rate"))
	p
	

	p1 <- ggplot(dfcr)
	p1 <- p1+geom_path(aes(y=yield_nat2_sc,x=yield_nat1_sc,color=slope))
	p1 <- p1+geom_point(aes(y=yield_nat2_sc,x=yield_nat1_sc,color=slope, shape=intercept),size=4,alpha=.6)
	#p1 <- p1+ geom_text(data=dft,aes(y=dft$yield_nat2_sc,x=dft$yield_nat1_sc,label="40:10"),fontface = "bold",colour="black")				
	
	p1 <- p1+facet_wrap(~scn, scales="free")
	p1 <- p1+ ylab(" average yield - Nation 2") + xlab(" average yield - Nation 1")
	p1 <- p1  + theme_bw(16)+ theme(legend.position = "bottom")
	p1 <- p1 + guides(color=guide_legend("harvest rate"), shape=guide_legend("biomass threshold"),size=guide_legend("none"))
	p1
	#p<-p+ geom_text(aes(y=nt2tf,x=nt1tf,label="10:40"),fontface = "bold",colour="black", position = position_nudge(y = -0.02))				
	#, size=guide_legend("biomass threshold")

	pto<-plot_grid(p,p1,nrow=2,rel_heights  = c(1, 1.3))
	pto
	

	if(sv==TRUE){
			setwd("/Users/catarinawor/Documents/Lagrangian/report/manuscript/ICES")
			ggsave(paste(nome,"tradeoff_nations_mov_scn.pdf",sep=""), plot=pto, width = 11, height = 9)
		}

		return(pto)

#dflu<-dflu[dflu$intercept_hcr<0.6,]

	#summary(dflu)


	
	
}




