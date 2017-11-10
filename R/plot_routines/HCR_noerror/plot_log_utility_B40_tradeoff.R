#==============================================================
# Read in the results from the HCR search - no assessment error
#Author: Catarina Wor
#Date: Sep 14 2017

#==============================================================




#M<-SIMSdatlim
#Msq<-SIMSsqlim

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
			utility_total[i]<-sum(log(apply(M[[ii]][[1]]$yYieldNat[101:160,],1,sum)+1))
			utility_total2[i]<-sum(log(M[[ii]][[1]]$yYieldNat[101:160,1]+1))*sum(log(M[[ii]][[1]]$yYieldNat[101:160,2]+1))
			utility_nat1[i]<-sum(log(M[[ii]][[1]]$yYieldNat[101:160,1]+1))
			utility_nat2[i]<-sum(log(M[[ii]][[1]]$yYieldNat[101:160,2]+1))
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
			utility_total2_sq[y]<-sum(log(Msq[[yy]][[1]]$yYieldNat[101:160,1]+1))*sum(log(Msq[[yy]][[1]]$yYieldNat[101:160,2]+1))
			utility_nat1_sq[y]<-sum(log(Msq[[yy]][[1]]$yYieldNat[101:160,1]+1))
			utility_nat2_sq[y]<-sum(log(Msq[[yy]][[1]]$yYieldNat[101:160,2]+1))

			B40_total_sq[y]<-sum((Msq[[yy]][[1]]$ytB/Msq[[yy]][[1]]$Bo)[iniyr:fimyr]<0.4)/length(iniyr:fimyr)
			scenario_sq[y]<-dd
			y=y+1
		}

	}

	dflu<-data.frame(slope_hcr=slope_hcr,intercept_hcr=intercept_hcr,
		utility_total2=utility_total2,utility_total=utility_total,utility_nat1=utility_nat1,
		utility_nat2=utility_nat2,b40=B40_total, scenario=scenario)

	dflusq<-data.frame(slope_hcr=slope_hcr_sq,intercept_hcr=intercept_hcr_sq,
		utility_total2=utility_total2_sq,utility_total=utility_total_sq,utility_nat1=utility_nat1_sq,
		utility_nat2=utility_nat2_sq,b40=B40_total_sq,scenario=scenario_sq)

	dflu_plot<-aggregate(dflu,list(dflu$slope_hcr,dflu$intercept_hcr, dflu$scenario), mean)

	dflu_plotsq<-aggregate(dflusq,list(dflusq$slope_hcr,dflusq$intercept_hcr, dflusq$scenario), mean)




	summary(dflu_plot)
	head(dflu_plotsq)

	dflu_selec<-dflu_plot

	dfa<-rbind(dflu_selec,dflu_plotsq)

	dfa$slope<-as.character(dfa$slope_hcr)
	dfa$slope[dfa$slope=="999"]<-"10:40"

	dfa$intercept<-as.character(dfa$intercept_hcr)
	dfa$intercept[dfa$intercept=="999"]<-"10:40"

	dfa$ut2tf<- -999
	dfa$ut2tf[which(dfa$slope=="10:40")]<-dfa$utility_total2[which(dfa$slope=="10:40")]

	dfa$nt1tf<- -999
	dfa$nt1tf[which(dfa$slope=="10:40")]<-dfa$utility_nat1[which(dfa$slope=="10:40")]

	dfa$nt2tf<- -999
	dfa$nt2tf[which(dfa$slope=="10:40")]<-dfa$utility_nat2[which(dfa$slope=="10:40")]

	dfa$scn<-"no strong recruitment"
	dfa$scn[dfa$scenario==2]<-"one strong recruitment"
	dfa$scn[dfa$scenario==3]<-"two strong recruitments"

	dfa$b40tf<- -999
	dfa$b40tf[dfa$slope=="10:40"]<-dfa$b40[which(dfa$slope=="10:40")]

	head(dfa)



	p<-ggplot(dfa)
	p<-p+geom_line(aes(y=utility_total2,x=b40,color=slope))
	p<-p+geom_point(aes(y=utility_total2,x=b40,color=slope, size=intercept),alpha=.6)
	#p<-p+ geom_text(aes(y=dfa$utility_total2[dfa$slope=="10:40"],x=dfa$b40[dfa$slope=="10:40"],label="10:40"),fontface = "bold",colour="black", position = position_nudge(y = -0.02))				
	p<-p+facet_wrap(~scn)
	p<-p+ xlim(c(0,.9))
	p<-p+ ylab("log utility - Total") + xlab("% time Bt < 0.4 Bo")
	p<-p+ geom_text(aes(y=ut2tf,x=b40tf,label="10:40"),fontface = "bold",colour="black", position = position_nudge(y = -0.02))				
	p <- p  + theme_bw(16)+ theme(legend.position = "bottom")
	p

	p1<-ggplot(dfa)
	p1<-p1+geom_line(aes(y=utility_nat1,x=b40,color=slope))
	p1<-p1+geom_point(aes(y=utility_nat1,x=b40,color=slope, size=intercept),alpha=.6)
	p1<-p1+facet_wrap(~scn)
	p1<-p1+ ylab("log utility - Nation 1") + xlab("")
	p1<-p1+ geom_text(aes(y=nt1tf,x=b40tf,label="10:40"),fontface = "bold",colour="black", position = position_nudge(y = -0.02))				
	p1<-p1+ xlim(c(0,.9))+ylim(c(0,max(dfa$utility_nat1)))
	p1 <- p1  + theme_bw(16)+ theme(legend.position = "none") 
	p1


	p2<-ggplot(dfa)
	p2<-p2+ geom_line(aes(y=utility_nat2,x=b40,color=slope))
	p2<-p2+ geom_point(aes(y=utility_nat2,x=b40,color=slope, size=intercept),alpha=.6)
	p2<-p2+ facet_wrap(~scn)
	p2<-p2+ ylab("log utility - Nation 2") + xlab("")
	p2<-p2+ geom_text(aes(y=nt2tf,x=b40tf,label="10:40"),fontface = "bold",colour="black", position = position_nudge(y = -0.02))				
	p2<-p2+ xlim(c(0,.9))+ylim(c(0,max(dfa$utility_nat2)))
	p2 <- p2  + theme_bw(16)+ theme(legend.position = "none") 
	p2

	pto<-plot_grid(p1,p2,p,nrow=3,rel_heights  = c(1,1, 1.3))


	if(sv==TRUE){
			setwd("/Users/catarinawor/Documents/Lagrangian/figures/HCR")
			ggsave(paste(nome,"tradeoff.pdf",sep=""), plot=pto, width = 9, height = 11)
		}

		return(pto)

#dflu<-dflu[dflu$intercept_hcr<0.6,]

	#summary(dflu)


	
	
}





plot_logUtility_all <- function( M , sv=F, nome="",nations=F)
{
	cat("plot_logUtility")

	#M<-SIMSdat

	length(M)

	n <- length( M )


	slope_hcr<-NULL
	intercept_hcr<-NULL
	utility_total<-NULL
	utility_total2<-NULL
	utility_nat1<-NULL
	utility_nat2<-NULL





	for(i in 1:n){

		slope_hcr[i]<-M[[i]][[1]]$slope_hcr
		intercept_hcr[i]<-M[[i]][[1]]$intercept_hcr
		utility_total[i]<-sum(log(apply(M[[i]][[1]]$yYieldNat[101:160,],1,sum)+0.01))
		utility_total2[i]<-sum(log(M[[i]][[1]]$yYieldNat[101:160,1]+1))*sum(log(M[[i]][[1]]$yYieldNat[101:160,2]+1))
		utility_nat1[i]<-sum(log(M[[i]][[1]]$yYieldNat[101:160,1]+0.01))
		utility_nat2[i]<-sum(log(M[[i]][[1]]$yYieldNat[101:160,2]+0.01))


	}


	slope_hcr_sq<-NULL
	intercept_hcr_sq<-NULL
	utility_total_sq<-NULL
	utility_total2_sq<-NULL
	utility_nat1_sq<-NULL
	utility_nat2_sq<-NULL

	nsq<- length( Msq )

	for(y in 1:nsq){

		slope_hcr_sq[y]<-999
		intercept_hcr_sq[y]<-999
		utility_total_sq[y]<-sum(log(apply(Msq[[y]][[1]]$yYieldNat[101:160,],1,sum)+0.01))
		utility_total2_sq[y]<-sum(log(Msq[[y]][[1]]$yYieldNat[101:160,1]+1))*sum(log(Msq[[y]][[1]]$yYieldNat[101:160,2]+1))
		utility_nat1_sq[y]<-sum(log(Msq[[y]][[1]]$yYieldNat[101:160,1]+0.01))
		utility_nat2_sq[y]<-sum(log(Msq[[y]][[1]]$yYieldNat[101:160,2]+0.01))


	}



	dflu<-data.frame(slope_hcr=slope_hcr,intercept_hcr=intercept_hcr,
		utility_total2=utility_total2,utility_total=utility_total,utility_nat1=utility_nat1,utility_nat2=utility_nat2)

	dflusq<-data.frame(slope_hcr=slope_hcr_sq,intercept_hcr=intercept_hcr_sq,
		utility_total2=utility_total2_sq,utility_total=utility_total_sq,utility_nat1=utility_nat1_sq,utility_nat2=utility_nat2_sq)


	#dflu<-dflu[dflu$intercept_hcr<0.6,]

	#summary(dflu)

	dflu_plot<-aggregate(dflu,list(dflu$slope_hcr,dflu$intercept_hcr), mean)

	dflu_plotsq<-aggregate(dflusq,list(dflusq$slope_hcr,dflusq$intercept_hcr), mean)
	
	utility_total2_comparison<-which.min(abs(dflu_plot$utility_total2-dflu_plotsq$utility_total2))
	utility_total_comparison<-which.min(abs(dflu_plot$utility_total-dflu_plotsq$utility_total))
	utility_nat1_comparison<-which.min(abs(dflu_plot$utility_nat1-dflu_plotsq$utility_nat1))
	utility_nat2_comparison<-which.min(abs(dflu_plot$utility_nat2-dflu_plotsq$utility_nat2))


	
	dflu_plotn<-data.frame(slope_hcr=rep(dflu_plot$slope_hcr,3),intercept_hcr=rep(dflu_plot$intercept_hcr,3),
		utility=c(dflu_plot$utility_total/max(dflu_plot$utility_total),dflu_plot$utility_nat1/max(dflu_plot$utility_nat1),
			dflu_plot$utility_nat2/max(dflu_plot$utility_nat2)),
		nation=rep(c("Total", "Nation 1", "Nation 2"),each=length(dflu_plot$intercept_hcr)))

	
	myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

	

	if(nations==F){
		plu <- ggplot(dflu_plot, aes(x=intercept_hcr,y=slope_hcr,z=utility_total))
		plu <- plu+geom_raster(aes(fill=utility_total))
		plu <- plu + scale_fill_gradientn(colours = myPalette(4))
		plu <- plu + theme_bw(16)
		print(plu)
	}else{
		#plun <- ggplot(dflu_plotn, aes(x=intercept_hcr,y=slope_hcr,z=utility))
		#plun <- plun +geom_raster(aes(fill=utility))
		#plun <- plun+facet_wrap(~nation)
		#plun  <- plun  + scale_fill_gradientn(colours = myPalette(6),name  ="mean log utility")
		#plun  <- plun  + theme_bw(16)
		#plun  <- plun  + ylab("Slope") + xlab("Intercept")

		#plun <- plun+geom_contour()
		#print(plun)
	
		tlabs<-round(range(dflu_plot$utility_total))+c(1,-1)
		plut <- ggplot(dflu_plot, aes(x=intercept_hcr,y=slope_hcr,z=utility_total))
		plut <- plut +geom_raster(aes(fill=utility_total))
		plut  <- plut  + scale_fill_gradientn(colours = myPalette(6),name  ="log utility",breaks=tlabs,labels=round(tlabs))
		plut  <- plut  + theme_bw(16) + theme(legend.position="top")
		plut  <- plut  + ylab("") + xlab("Intercept")
		plut  <- plut  + ggtitle("Total sum")
		plut  <- plut  +  theme(plot.title = element_text(hjust = 0.5))
		#plut  <- plut  + geom_text(aes(x=intercept_hcr[which.max(utility_total)],y=slope_hcr[which.max(utility_total)],label=round(max(dflu_plot$utility_total)/1000,2)),fontface = "bold",colour="gray90")
		plut  <- plut  + geom_text(aes(x=intercept_hcr[which.max(utility_total)],y=slope_hcr[which.max(utility_total)],label="max"),fontface = "bold",colour="gray90")
		plut  <- plut  + geom_point(aes(x=intercept_hcr[utility_total_comparison],y=slope_hcr[utility_total_comparison]),colour="gray90", size=10)
		#plut  <- plut  + geom_text(aes(x=intercept_hcr[utility_total_comparison],y=slope_hcr[utility_total_comparison],label=round(dflu_plotsq$utility_total/1000,2)),fontface = "bold",colour="black")		
		plut  <- plut  + geom_text(aes(x=intercept_hcr[utility_total_comparison],y=slope_hcr[utility_total_comparison],label="10:40"),fontface = "bold",colour="black")		
		plut  <- plut  + theme(legend.position = "none") 
		#plut

		tlabs<-round(range(dflu_plot$utility_total2))+c(1,-1)
		plut2 <- ggplot(dflu_plot, aes(x=intercept_hcr,y=slope_hcr,z=utility_total2))
		plut2 <- plut2 +geom_raster(aes(fill=utility_total2))
		plut2  <- plut2  + scale_fill_gradientn(colours = myPalette(6),name  ="log utility",breaks=tlabs,labels=c("low","high"))
		plut2  <- plut2  + theme_bw(16) + theme(legend.position="top")
		plut2  <- plut2  + ylab("") + xlab("Intercept")
		plut2  <- plut2  + ggtitle("Total")
		plut2  <- plut2  + guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                label.position = "bottom"))
		#plut2  <- plut2  + geom_text(aes(x=intercept_hcr[which.max(utility_total2)],y=slope_hcr[which.max(utility_total2)],label=round(max(dflu_plot$utility_total2)/1000,2)),fontface = "bold",colour="gray90")
		plut2  <- plut2  + geom_text(aes(x=intercept_hcr[which.max(utility_total2)],y=slope_hcr[which.max(utility_total2)],label="max"),fontface = "bold",colour="gray90")
		plut2  <- plut2  + geom_point(aes(x=intercept_hcr[utility_total2_comparison],y=slope_hcr[utility_total2_comparison]),colour="gray90", size=10)
		#plut2  <- plut2  + geom_text(aes(x=intercept_hcr[utility_total2_comparison],y=slope_hcr[utility_total2_comparison],label=round(dflu_plotsq$utility_total2/1000,2)),fontface = "bold",colour="black")
		plut2  <- plut2  + geom_text(aes(x=intercept_hcr[utility_total2_comparison],y=slope_hcr[utility_total2_comparison],label="10:40"),fontface = "bold",colour="black")
		plut2  <- plut2  +  theme(plot.title = element_text(hjust = 0.5))
		#plut2  <- plut2  + theme(legend.position = "none") 
		plut2

		n1labs<-round(range(dflu_plot$utility_nat1))+c(1,-1)
		plun1 <- ggplot(dflu_plot, aes(x=intercept_hcr,y=slope_hcr,z=utility_nat1))
		plun1 <- plun1 + geom_raster(aes(fill=utility_nat1))
		plun1 <- plun1 + scale_fill_gradientn(colours = myPalette(6),name  ="log utility",breaks=n1labs,labels=c("low","high"))
		plun1 <- plun1 + theme_bw(16) + theme(legend.position="top")
		plun1 <- plun1 + ylab("Slope") + xlab("Intercept")
		plun1 <- plun1 + ggtitle("Nation 1")
		plun1 <- plun1 + guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                label.position = "bottom"))
		plun1 <- plun1 +  theme(plot.title = element_text(hjust = 0.5))
		#plun1 <- plun1  + geom_text(aes(x=intercept_hcr[which.max(utility_nat1)],y=slope_hcr[which.max(utility_nat1)],label=round(max(dflu_plot$utility_nat1)/1000,2)),fontface = "bold",colour="gray90")
		plun1 <- plun1  + geom_text(aes(x=intercept_hcr[which.max(utility_nat1)],y=slope_hcr[which.max(utility_nat1)],label="max"),fontface = "bold",colour="gray90")
		plun1  <- plun1  + geom_point(aes(x=intercept_hcr[utility_nat1_comparison],y=slope_hcr[utility_nat1_comparison]),colour="gray90", size=10)
		#plun1  <- plun1  + geom_text(aes(x=intercept_hcr[utility_nat1_comparison],y=slope_hcr[utility_nat1_comparison],label=round(dflu_plotsq$utility_nat1/1000,2)),fontface = "bold",colour="black")
		plun1  <- plun1  + geom_text(aes(x=intercept_hcr[utility_nat1_comparison],y=slope_hcr[utility_nat1_comparison],label="10:40"),fontface = "bold",colour="black")	
		#plun1 <- plun1  + theme(legend.position = "none") 
		plun1

		n2labs<-round(range(dflu_plot$utility_nat2))+c(1,-1)
		plun2 <- ggplot(dflu_plot, aes(x=intercept_hcr,y=slope_hcr,z=utility_nat2))
		plun2 <- plun2 + geom_raster(aes(fill=utility_nat2))
		plun2 <- plun2  + scale_fill_gradientn(colours = myPalette(6),name ="log utility",breaks=n2labs,labels=c("low","high"))
		plun2 <- plun2  + theme_bw(16) + theme(legend.position="top")
		plun2 <- plun2  + ylab("") + xlab("Intercept")
		plun2 <- plun2  + ggtitle("Nation 2")
		plun2 <- plun2  + guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                label.position = "bottom"))
		plun2 <- plun2  +  theme(plot.title = element_text(hjust = 0.5))
		#plun2  <- plun2  + geom_text(aes(x=intercept_hcr[which.max(utility_nat2)],y=slope_hcr[which.max(utility_nat2)],label=round(max(dflu_plot$utility_nat2)/1000,2)),fontface = "bold",colour="gray90")
		plun2  <- plun2  + geom_text(aes(x=intercept_hcr[which.max(utility_nat2)],y=slope_hcr[which.max(utility_nat2)],label="max"),fontface = "bold",colour="gray90")
		plun2  <- plun2  + geom_point(aes(x=intercept_hcr[utility_nat2_comparison],y=slope_hcr[utility_nat2_comparison]),colour="gray90", size=10)
		plun2  <- plun2  + geom_text(aes(x=intercept_hcr[utility_nat2_comparison],y=slope_hcr[utility_nat2_comparison],label="10:40"),fontface = "bold",colour="black")
		#plun2  <- plun2  + geom_text(aes(x=intercept_hcr[utility_nat2_comparison],y=slope_hcr[utility_nat2_comparison],label=round(dflu_plotsq$utility_nat2/1000,2)),fontface = "bold",colour="black")
		#plun2  <- plun2  + theme(legend.position = "none") 
		plun2

		plun<-plot_grid(plun1,plun2,plut2,ncol=3)
		print(plun)

		if(sv==TRUE){
			setwd("/Users/catarinawor/Documents/Lagrangian/report/HCR_coarse")
			ggsave(paste(nome,"logUtility.pdf",sep=""), plot=plun, width = 15, height = 5)
		}

		return(plun)
	}
	
}











	

