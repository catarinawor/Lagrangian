#==============================================================
# Read in the results from the HCR search - no assessment error
#Author: Catarina Wor
#Date: Sep 14 2017

#==============================================================




#M<-SIMSdat
#Msq<-SIMSsq

require(reshape2)
require(tidyr)
require(ggplot2)
library(scales)
library(RColorBrewer)


plot_Yield <- function( M, Msq, sv=F, nome="",nations=F)
{
	cat("plot_Yield")


	n <- length( M )


	slope_hcr<-NULL
	intercept_hcr<-NULL
	Yield_total<-NULL
	Yield_nat1<-NULL
	Yield_nat2<-NULL


	for(i in 1:n){

		slope_hcr[i]<-M[[i]][[1]]$slope_hcr
		intercept_hcr[i]<-M[[i]][[1]]$intercept_hcr
		Yield_total[i]<-mean(apply(M[[i]][[1]]$yYieldNat[101:160,],1,sum))
		Yield_nat1[i]<-mean((M[[i]][[1]]$yYieldNat[101:160,1]))
		Yield_nat2[i]<-mean((M[[i]][[1]]$yYieldNat[101:160,2]))
	}


	slope_hcr_sq<-NULL
	intercept_hcr_sq<-NULL
	Yield_total_sq<-NULL
	Yield_nat1_sq<-NULL
	Yield_nat2_sq<-NULL

	nsq<- length( Msq )

	for(y in 1:nsq){

		slope_hcr_sq[y]<-999
		intercept_hcr_sq[y]<-999
		Yield_total_sq[y]<-mean(apply(Msq[[y]][[1]]$yYieldNat[101:160,],1,sum))
		Yield_nat1_sq[y]<-mean((Msq[[y]][[1]]$yYieldNat[101:160,1]))
		Yield_nat2_sq[y]<-mean((Msq[[y]][[1]]$yYieldNat[101:160,2]))
	}


	dfY<-data.frame(slope_hcr=slope_hcr,intercept_hcr=intercept_hcr,
		Yield_total=Yield_total,Yield_nat1=Yield_nat1,Yield_nat2=Yield_nat2)

	dfYsq<-data.frame(slope_hcr=slope_hcr_sq,intercept_hcr=intercept_hcr_sq,
		Yield_total=Yield_total_sq,Yield_nat1=Yield_nat1_sq,Yield_nat2=Yield_nat2_sq)

	#summary(dfY)

	dfY_plot<-aggregate(dfY,list(dfY$slope_hcr,dfY$intercept_hcr), median)

	dfY_plotsq<-aggregate(dfYsq,list(dfYsq$slope_hcr,dfYsq$intercept_hcr), median)
	
	Yield_total_comparison<-which.min(abs(dfY_plot$Yield_total-dfY_plotsq$Yield_total))
	Yield_nat1_comparison<-which.min(abs(dfY_plot$Yield_nat1-dfY_plotsq$Yield_nat1))
	Yield_nat2_comparison<-which.min(abs(dfY_plot$Yield_nat2-dfY_plotsq$Yield_nat2))



	#plu <- ggplot(dflu_plot, aes(x=intercept_hcr,y=slope_hcr))
	#plu <- plu+geom_point(aes(color=utility_total,size=utility_total))
	#plu


	dfY_plotn<-data.frame(slope_hcr=rep(dfY_plot$slope_hcr,3),intercept_hcr=rep(dfY_plot$intercept_hcr,3),
		Yield=c(dfY_plot$Yield_total/max(dfY_plot$Yield_total),dfY_plot$Yield_nat1/max(dfY_plot$Yield_nat1),
			dfY_plot$Yield_nat2/max(dfY_plot$Yield_nat2)),nation=rep(c("Total", "Nation 1", "Nation 2"),each=length(dfY_plot$intercept_hcr)))

	#dfY_plotn<-data.frame(slope_hcr=rep(dfY_plot$slope_hcr,3),intercept_hcr=rep(dfY_plot$intercept_hcr,3),
	#	Yield=c(dfY_plot$Yield_total,dfY_plot$Yield_nat1,
	#		dfY_plot$Yield_nat2),nation=rep(c("Total", "Nation 1", "Nation 2"),each=length(dfY_plot$intercept_hcr)))

	
	myPalette <- colorRampPalette(rev(brewer.pal(8, "Spectral")))
	#myPalette <- colorRampPalette((brewer.pal(8, "Greys")))
	

	if(nations==F){
		pY <- ggplot(dfY_plot, aes(x=intercept_hcr,y=slope_hcr,z=Yield_total))
		pY <- pY+geom_raster(aes(fill=Yield_total))
		pY <- pY + scale_fill_gradientn(colours = myPalette(4))
		pY <- pY + theme_bw(16)
		print(pY)
	}else{
		#pYn <- ggplot(dfY_plotn, aes(x=intercept_hcr,y=slope_hcr,z=Yield))
		#pYn <- pYn +geom_raster(aes(fill=Yield))
		#pYn <- pYn+facet_wrap(~nation)
		#pYn <- pYn + scale_fill_gradientn(colours = myPalette(4),name  ="scaled mean Yield")
		#pYn <- pYn + theme_bw(16)
		#pYn <- pYn + ylab("Slope") + xlab("Intercept")

		tlabs<-round(range(dfY_plot$Yield_total))+c(1,-1)
		pYt <- ggplot(dfY_plot, aes(x=intercept_hcr,y=slope_hcr,z=Yield_total))
		pYt <- pYt +geom_raster(aes(fill=Yield_total))
		pYt <- pYt + scale_fill_gradientn(colours = myPalette(11),name  ="Yield",breaks=tlabs,labels=c("low","high"))
		pYt <- pYt + theme_bw(16) + theme(legend.position="right")
		pYt <- pYt + ylab("") + xlab("Intercept")
		pYt <- pYt  + ggtitle("Total")
		pYt <- pYt  +  theme(plot.title = element_text(hjust = 0.5))
		pYt  <- pYt  + guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5))
		pYt  <- pYt  + geom_text(aes(x=dfY_plot$intercept_hcr[which.max(dfY_plot$Yield_total)],y=dfY_plot$slope_hcr[which.max(dfY_plot$Yield_total)],label="max"),fontface = "bold",colour="gray90")
		pYt  <- pYt  + geom_point(aes(x=dfY_plot$intercept_hcr[Yield_total_comparison],y=dfY_plot$slope_hcr[Yield_total_comparison]),colour="gray90", size=10)
		pYt  <- pYt  + geom_text(aes(x=dfY_plot$intercept_hcr[Yield_total_comparison],y=dfY_plot$slope_hcr[Yield_total_comparison],label="40:10"),fontface = "bold",colour="black")
		pYt

		n1labs<-round(range(dfY_plot$Yield_nat1))+c(1,-1)
		pYn1 <- ggplot(dfY_plot, aes(x=intercept_hcr,y=slope_hcr,z=Yield_nat1))
		pYn1 <- pYn1 +geom_raster(aes(fill=Yield_nat1))
		pYn1 <- pYn1 + scale_fill_gradientn(colours = myPalette(11),name  ="Yield",breaks=n1labs,labels=c("low","high"))
		pYn1 <- pYn1 + theme_bw(16) + theme(legend.position="none")
		pYn1 <- pYn1 + ylab("Slope") + xlab("Intercept")
		pYn1 <- pYn1  + ggtitle("Nation 1")
		pYn1 <- pYn1  +  theme(plot.title = element_text(hjust = 0.5))
		pYn1 <- pYn1  + guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                label.position = "bottom"))
		pYn1  <- pYn1  + geom_text(aes(x=dfY_plot$intercept_hcr[which.max(dfY_plot$Yield_nat1)],y=dfY_plot$slope_hcr[which.max(dfY_plot$Yield_nat1)],label="max"),fontface = "bold",colour="gray90")
		pYn1  <- pYn1  + geom_point(aes(x=dfY_plot$intercept_hcr[Yield_nat1_comparison],y=dfY_plot$slope_hcr[Yield_nat1_comparison]),colour="gray90", size=10)
		pYn1  <- pYn1  + geom_text(aes(x=dfY_plot$intercept_hcr[Yield_nat1_comparison],y=dfY_plot$slope_hcr[Yield_nat1_comparison],label="40:10"),fontface = "bold",colour="black")
		pYn1

		n2labs<-round(range(dfY_plot$Yield_nat2))+c(1,-1)
		pYn2 <- ggplot(dfY_plot, aes(x=intercept_hcr,y=slope_hcr,z=Yield_nat2))
		pYn2 <- pYn2 +geom_raster(aes(fill=Yield_nat2))
		pYn2 <- pYn2 + scale_fill_gradientn(colours = myPalette(11),name  ="Yield",breaks=n2labs,labels=c("low","high"))
		pYn2 <- pYn2 + theme_bw(16) + theme(legend.position="none")
		pYn2 <- pYn2 + ylab("") + xlab("Intercept")
		pYn2 <- pYn2  + ggtitle("Nation 2")
		pYn2 <- pYn2  +  theme(plot.title = element_text(hjust = 0.5))
		pYn2 <- pYn2  + guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                label.position = "bottom"))
		pYn2  <- pYn2  + geom_text(aes(x=dfY_plot$intercept_hcr[which.max(dfY_plot$Yield_nat2)],y=dfY_plot$slope_hcr[which.max(dfY_plot$Yield_nat2)],label="max"),fontface = "bold",colour="gray90")
		pYn2  <- pYn2  + geom_point(aes(x=dfY_plot$intercept_hcr[Yield_nat2_comparison],y=dfY_plot$slope_hcr[Yield_nat2_comparison]),colour="gray90", size=10)
		pYn2  <- pYn2  + geom_text(aes(x=dfY_plot$intercept_hcr[Yield_nat2_comparison],y=dfY_plot$slope_hcr[Yield_nat2_comparison],label="40:10"),fontface = "bold",colour="black")

		pYn2

		pYn<-plot_grid(pYn1,pYn2,pYt,ncol=3,rel_widths = c(1,1, 1.2))
		#pYn <- pYn+geom_contour()
		print(pYn)

		if(sv==TRUE){
			setwd("/Users/catarinawor/Documents/Lagrangian/report/HCR_coarse")
			ggsave(paste(nome,"_smYield.pdf",sep=""), plot=pYn, width = 15, height = 5)
		}

		return(pYn)

	}
	
}


plot_Yield_later <- function( M, Msq, sv=F, nome="",nations=F)
{
	cat("plot_Yield")


	n <- length( M )


	slope_hcr<-NULL
	intercept_hcr<-NULL
	Yield_total<-NULL
	Yield_nat1<-NULL
	Yield_nat2<-NULL


	for(i in 1:n){

		slope_hcr[i]<-M[[i]][[1]]$slope_hcr
		intercept_hcr[i]<-M[[i]][[1]]$intercept_hcr
		Yield_total[i]<-mean(apply(M[[i]][[1]]$yYieldNat[141:160,],1,sum))
		Yield_nat1[i]<-mean((M[[i]][[1]]$yYieldNat[141:160,1]))
		Yield_nat2[i]<-mean((M[[i]][[1]]$yYieldNat[141:160,2]))
	}


	slope_hcr_sq<-NULL
	intercept_hcr_sq<-NULL
	Yield_total_sq<-NULL
	Yield_nat1_sq<-NULL
	Yield_nat2_sq<-NULL

	nsq<- length( Msq )

	for(y in 1:nsq){

		slope_hcr_sq[y]<-999
		intercept_hcr_sq[y]<-999
		Yield_total_sq[y]<-mean(apply(Msq[[y]][[1]]$yYieldNat[141:160,],1,sum))
		Yield_nat1_sq[y]<-mean((Msq[[y]][[1]]$yYieldNat[141:160,1]))
		Yield_nat2_sq[y]<-mean((Msq[[y]][[1]]$yYieldNat[141:160,2]))
	}


	dfY<-data.frame(slope_hcr=slope_hcr,intercept_hcr=intercept_hcr,
		Yield_total=Yield_total,Yield_nat1=Yield_nat1,Yield_nat2=Yield_nat2)

	dfYsq<-data.frame(slope_hcr=slope_hcr_sq,intercept_hcr=intercept_hcr_sq,
		Yield_total=Yield_total_sq,Yield_nat1=Yield_nat1_sq,Yield_nat2=Yield_nat2_sq)

	#summary(dfY)

	dfY_plot<-aggregate(dfY,list(dfY$slope_hcr,dfY$intercept_hcr), median)

	dfY_plotsq<-aggregate(dfYsq,list(dfYsq$slope_hcr,dfYsq$intercept_hcr), median)
	
	Yield_total_comparison<-which.min(abs(dfY_plot$Yield_total-dfY_plotsq$Yield_total))
	Yield_nat1_comparison<-which.min(abs(dfY_plot$Yield_nat1-dfY_plotsq$Yield_nat1))
	Yield_nat2_comparison<-which.min(abs(dfY_plot$Yield_nat2-dfY_plotsq$Yield_nat2))



	#plu <- ggplot(dflu_plot, aes(x=intercept_hcr,y=slope_hcr))
	#plu <- plu+geom_point(aes(color=utility_total,size=utility_total))
	#plu


	dfY_plotn<-data.frame(slope_hcr=rep(dfY_plot$slope_hcr,3),intercept_hcr=rep(dfY_plot$intercept_hcr,3),
		Yield=c(dfY_plot$Yield_total/max(dfY_plot$Yield_total),dfY_plot$Yield_nat1/max(dfY_plot$Yield_nat1),
			dfY_plot$Yield_nat2/max(dfY_plot$Yield_nat2)),nation=rep(c("Total", "Nation 1", "Nation 2"),each=length(dfY_plot$intercept_hcr)))

	#dfY_plotn<-data.frame(slope_hcr=rep(dfY_plot$slope_hcr,3),intercept_hcr=rep(dfY_plot$intercept_hcr,3),
	#	Yield=c(dfY_plot$Yield_total,dfY_plot$Yield_nat1,
	#		dfY_plot$Yield_nat2),nation=rep(c("Total", "Nation 1", "Nation 2"),each=length(dfY_plot$intercept_hcr)))

	
	myPalette <- colorRampPalette(rev(brewer.pal(8, "Spectral")))
	#myPalette <- colorRampPalette((brewer.pal(8, "Greys")))
	

	if(nations==F){
		pY <- ggplot(dfY_plot, aes(x=intercept_hcr,y=slope_hcr,z=Yield_total))
		pY <- pY+geom_raster(aes(fill=Yield_total))
		pY <- pY + scale_fill_gradientn(colours = myPalette(4))
		pY <- pY + theme_bw(16)
		print(pY)
	}else{
		#pYn <- ggplot(dfY_plotn, aes(x=intercept_hcr,y=slope_hcr,z=Yield))
		#pYn <- pYn +geom_raster(aes(fill=Yield))
		#pYn <- pYn+facet_wrap(~nation)
		#pYn <- pYn + scale_fill_gradientn(colours = myPalette(4),name  ="scaled mean Yield")
		#pYn <- pYn + theme_bw(16)
		#pYn <- pYn + ylab("Slope") + xlab("Intercept")

		tlabs<-round(range(dfY_plot$Yield_total))+c(1,-1)
		pYt <- ggplot(dfY_plot, aes(x=intercept_hcr,y=slope_hcr,z=Yield_total))
		pYt <- pYt +geom_raster(aes(fill=Yield_total))
		pYt <- pYt + scale_fill_gradientn(colours = myPalette(11),name  ="Yield",breaks=tlabs,labels=c("low","high"))
		pYt <- pYt + theme_bw(16) + theme(legend.position="right")
		pYt <- pYt + ylab("") + xlab("Intercept")
		pYt <- pYt  + ggtitle("Total")
		pYt <- pYt  +  theme(plot.title = element_text(hjust = 0.5))
		pYt  <- pYt  + guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5))
		pYt  <- pYt  + geom_text(aes(x=dfY_plot$intercept_hcr[which.max(dfY_plot$Yield_total)],y=dfY_plot$slope_hcr[which.max(dfY_plot$Yield_total)],label="max"),fontface = "bold",colour="gray90")
		pYt  <- pYt  + geom_point(aes(x=dfY_plot$intercept_hcr[Yield_total_comparison],y=dfY_plot$slope_hcr[Yield_total_comparison]),colour="gray90", size=10)
		pYt  <- pYt  + geom_text(aes(x=dfY_plot$intercept_hcr[Yield_total_comparison],y=dfY_plot$slope_hcr[Yield_total_comparison],label="40:10"),fontface = "bold",colour="black")
		pYt

		n1labs<-round(range(dfY_plot$Yield_nat1))+c(1,-1)
		pYn1 <- ggplot(dfY_plot, aes(x=intercept_hcr,y=slope_hcr,z=Yield_nat1))
		pYn1 <- pYn1 +geom_raster(aes(fill=Yield_nat1))
		pYn1 <- pYn1 + scale_fill_gradientn(colours = myPalette(11),name  ="Yield",breaks=n1labs,labels=c("low","high"))
		pYn1 <- pYn1 + theme_bw(16) + theme(legend.position="none")
		pYn1 <- pYn1 + ylab("Slope") + xlab("Intercept")
		pYn1 <- pYn1  + ggtitle("Nation 1")
		pYn1 <- pYn1  +  theme(plot.title = element_text(hjust = 0.5))
		pYn1 <- pYn1  + guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                label.position = "bottom"))
		pYn1  <- pYn1  + geom_text(aes(x=dfY_plot$intercept_hcr[which.max(dfY_plot$Yield_nat1)],y=dfY_plot$slope_hcr[which.max(dfY_plot$Yield_nat1)],label="max"),fontface = "bold",colour="gray90")
		pYn1  <- pYn1  + geom_point(aes(x=dfY_plot$intercept_hcr[Yield_nat1_comparison],y=dfY_plot$slope_hcr[Yield_nat1_comparison]),colour="gray90", size=10)
		pYn1  <- pYn1  + geom_text(aes(x=dfY_plot$intercept_hcr[Yield_nat1_comparison],y=dfY_plot$slope_hcr[Yield_nat1_comparison],label="40:10"),fontface = "bold",colour="black")
		pYn1

		n2labs<-round(range(dfY_plot$Yield_nat2))+c(1,-1)
		pYn2 <- ggplot(dfY_plot, aes(x=intercept_hcr,y=slope_hcr,z=Yield_nat2))
		pYn2 <- pYn2 +geom_raster(aes(fill=Yield_nat2))
		pYn2 <- pYn2 + scale_fill_gradientn(colours = myPalette(11),name  ="Yield",breaks=n2labs,labels=c("low","high"))
		pYn2 <- pYn2 + theme_bw(16) + theme(legend.position="none")
		pYn2 <- pYn2 + ylab("") + xlab("Intercept")
		pYn2 <- pYn2  + ggtitle("Nation 2")
		pYn2 <- pYn2  +  theme(plot.title = element_text(hjust = 0.5))
		pYn2 <- pYn2  + guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                label.position = "bottom"))
		pYn2  <- pYn2  + geom_text(aes(x=dfY_plot$intercept_hcr[which.max(dfY_plot$Yield_nat2)],y=dfY_plot$slope_hcr[which.max(dfY_plot$Yield_nat2)],label="max"),fontface = "bold",colour="gray90")
		pYn2  <- pYn2  + geom_point(aes(x=dfY_plot$intercept_hcr[Yield_nat2_comparison],y=dfY_plot$slope_hcr[Yield_nat2_comparison]),colour="gray90", size=10)
		pYn2  <- pYn2  + geom_text(aes(x=dfY_plot$intercept_hcr[Yield_nat2_comparison],y=dfY_plot$slope_hcr[Yield_nat2_comparison],label="40:10"),fontface = "bold",colour="black")

		pYn2

		pYn<-plot_grid(pYn1,pYn2,pYt,ncol=3,rel_widths = c(1,1, 1.2))
		#pYn <- pYn+geom_contour()
		print(pYn)

		if(sv==TRUE){
			setwd("/Users/catarinawor/Documents/Lagrangian/report/HCR_coarse")
			ggsave(paste(nome,"_smYield.pdf",sep=""), plot=pYn, width = 15, height = 5)
		}

		return(pYn)

	}
	
}

	

	

