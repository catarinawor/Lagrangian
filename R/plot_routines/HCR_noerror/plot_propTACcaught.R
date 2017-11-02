#==============================================================
# Read in the results from the HCR search - no assessment error
#Author: Catarina Wor
#Date: Sep 14 2017

#==============================================================




#M<-SIMSdat

require(reshape2)
require(tidyr)
require(ggplot2)
library(scales)
library(RColorBrewer)


plot_pTAC <- function( M , Msq, sv=F, nome="",nations=F)
{
	cat("plot_propTACcaught")


	n <- length( M )

	slope_hcr<-NULL
	intercept_hcr<-NULL
	cTAC_total<-NULL
	cTAC_nat1<-NULL
	cTAC_nat2<-NULL

	for(i in 1:n){

		iniyr<-M[[i]][[1]]$"nyr"+1
		fimyr<-M[[i]][[1]]$"proj_yr"

		slope_hcr[i]<-M[[i]][[1]]$slope_hcr
		intercept_hcr[i]<-M[[i]][[1]]$intercept_hcr
		cTAC_total[i]<-mean((apply(M[[i]][[1]]$yYieldNat[iniyr:(fimyr),],1,sum)/
					apply(M[[i]][[1]]$histTAC[-1,],1,sum))[apply(M[[i]][[1]]$histTAC[-1,],1,sum)>0])


		cTAC_nat1[i]<-mean((M[[i]][[1]]$yYieldNat[iniyr:(fimyr),1]/M[[i]][[1]]$histTAC[-1,1])[M[[i]][[1]]$histTAC[-1,1]>0])


		cTAC_nat2[i]<-mean((M[[i]][[1]]$yYieldNat[iniyr:(fimyr),2]/M[[i]][[1]]$histTAC[-1,2])[M[[i]][[1]]$histTAC[-1,2]>0])

	}


	nsq <- length( Msq )

	slope_hcr_sq<-NULL
	intercept_hcr_sq<-NULL
	cTAC_total_sq<-NULL
	cTAC_nat1_sq<-NULL
	cTAC_nat2_sq<-NULL

	for(i in 1:nsq){

		iniyr<-Msq[[i]][[1]]$"nyr"+1
		fimyr<-Msq[[i]][[1]]$"proj_yr"

		slope_hcr_sq[i]<-Msq[[i]][[1]]$slope_hcr
		intercept_hcr_sq[i]<-Msq[[i]][[1]]$intercept_hcr
		cTAC_total_sq[i]<-mean((apply(Msq[[i]][[1]]$yYieldNat[iniyr:(fimyr),],1,sum)/
					apply(Msq[[i]][[1]]$histTAC[-1,],1,sum))[apply(Msq[[i]][[1]]$histTAC[-1,],1,sum)>0])


		cTAC_nat1_sq[i]<-mean((Msq[[i]][[1]]$yYieldNat[iniyr:(fimyr),1]/Msq[[i]][[1]]$histTAC[-1,1])[Msq[[i]][[1]]$histTAC[-1,1]>0])


		cTAC_nat2_sq[i]<-mean((Msq[[i]][[1]]$yYieldNat[iniyr:(fimyr),2]/Msq[[i]][[1]]$histTAC[-1,2])[Msq[[i]][[1]]$histTAC[-1,2]>0])

	}

	
	

	dfcTAC<-data.frame(slope_hcr=slope_hcr,intercept_hcr=intercept_hcr,
		cTAC_total=cTAC_total,cTAC_nat1=cTAC_nat1,cTAC_nat2=cTAC_nat2)
	
	dfcTACsq<-data.frame(slope_hcr=slope_hcr_sq,intercept_hcr=intercept_hcr_sq,
		cTAC_total=cTAC_total_sq,cTAC_nat1=cTAC_nat1_sq,cTAC_nat2=cTAC_nat2_sq)
	


	dfcTAC_plot<-aggregate(dfcTAC,list(dfcTAC$slope_hcr,dfcTAC$intercept_hcr), median)

	dfcTAC_plotsq<-aggregate(dfcTACsq,list(dfcTACsq$slope_hcr,dfcTACsq$intercept_hcr), median)
	
	cTAC_total_comparison<-which.min(abs(dfcTAC$cTAC_total-dfcTAC_plotsq$cTAC_total))
	cTAC_nat1_comparison<-which.min(abs(dfcTAC$cTAC_nat1-dfcTAC_plotsq$cTAC_nat1))
	cTAC_nat2_comparison<-which.min(abs(dfcTAC$cTAC_nat2-dfcTAC_plotsq$cTAC_nat2))


	summary(dfcTAC_plot)

	dfcTAC_plotn<-data.frame(slope_hcr=rep(dfcTAC_plot$slope_hcr,3),intercept_hcr=rep(dfcTAC_plot$intercept_hcr,3),
		cTAC=c(dfcTAC_plot$cTAC_total,dfcTAC_plot$cTAC_nat1,dfcTAC_plot$cTAC_nat2),nation=rep(c("Total", "Nation 1", "Nation 2"),each=length(dfcTAC_plot$intercept_hcr)))

	summary(dfcTAC_plotn)
	myPalette <- colorRampPalette(rev(brewer.pal(6, "Spectral")))
	#myPalette <- colorRampPalette((brewer.pal(8, "Greys")))
	

	if(nations==F){
		pct <- ggplot(dfcTAC_plot, aes(x=intercept_hcr,y=slope_hcr,z=cTAC_total))
		pct <- pct+geom_raster(aes(fill=cTAC_total))
		pct <- pct + scale_fill_gradientn(colours = myPalette(4))
		pct <- pct + theme_bw(16)
		print(pct)
	}else{
		#pctn <- ggplot(dfcTAC_plotn, aes(x=intercept_hcr,y=slope_hcr,z=cTAC))
		#pctn <- pctn +geom_raster(aes(fill=cTAC))
		#pctn <- pctn+facet_wrap(~nation)
		#pctn <- pctn + scale_fill_gradientn(colours = myPalette(4),name  ="% of TAC caught")
		#pctn <- pctn + theme_bw(16)
		#pctn <- pctn + ylab("Slope") + xlab("Intercept")
		#pctn

		pctt <- ggplot(dfcTAC_plot, aes(x=intercept_hcr,y=slope_hcr,z=cTAC_total))
		pctt <- pctt +geom_raster(aes(fill=cTAC_total))
		pctt <- pctt + scale_fill_gradientn(colours = myPalette(4),name  ="% of TAC caught")
		pctt <- pctt + theme_bw(16) + theme(legend.position="top")
		pctt <- pctt + ylab("") + xlab("Intercept")
		pctt <- pctt + ggtitle("Total")
		pctt <- pctt + theme(plot.title = element_text(hjust = 0.5))
		pctt <- pctt + guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                label.position = "bottom"))
		pctt <- pctt + geom_text(aes(x=dfcTAC_plot$intercept_hcr[which.max(dfcTAC_plot$cTAC_total)],y=dfcTAC_plot$slope_hcr[which.max(dfcTAC_plot$cTAC_total)],label="max"),fontface = "bold",colour="gray90")
		pctt <- pctt + geom_point(aes(x=dfcTAC$intercept_hcr[cTAC_total_comparison],y=dfcTAC$slope_hcr[cTAC_total_comparison]),colour="gray90", size=10)
		pctt <- pctt + geom_text(aes(x=dfcTAC$intercept_hcr[cTAC_total_comparison],y=dfcTAC$slope_hcr[cTAC_total_comparison],label="10:40"),fontface = "bold",colour="black")
		pctt

		pctn1 <- ggplot(dfcTAC_plot, aes(x=intercept_hcr,y=slope_hcr,z=cTAC_nat1))
		pctn1 <- pctn1 +geom_raster(aes(fill=cTAC_nat1))
		pctn1 <- pctn1 + scale_fill_gradientn(colours = myPalette(4),name  ="% of TAC caught")
		pctn1 <- pctn1 + theme_bw(16) + theme(legend.position="top")
		pctn1 <- pctn1 + ylab("Slope") + xlab("Intercept")
		pctn1 <- pctn1 + ggtitle("Nation 1 ")
		pctn1 <- pctn1 + theme(plot.title = element_text(hjust = 0.5))
		pctn1 <- pctn1 + guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                label.position = "bottom"))
		pctn1 <- pctn1 + geom_text(aes(x=dfcTAC_plot$intercept_hcr[which.max(dfcTAC_plot$cTAC_nat1)],y=dfcTAC_plot$slope_hcr[which.max(dfcTAC_plot$cTAC_nat1)],label="max"),fontface = "bold",colour="gray90")
		pctn1 <- pctn1 + geom_point(aes(x=dfcTAC$intercept_hcr[cTAC_nat1_comparison],y=dfcTAC$slope_hcr[cTAC_nat1_comparison]),colour="gray90", size=10)
		pctn1 <- pctn1 + geom_text(aes(x=dfcTAC$intercept_hcr[cTAC_nat1_comparison],y=dfcTAC$slope_hcr[cTAC_nat1_comparison],label="10:40"),fontface = "bold",colour="black")

		pctn1

		pctn2 <- ggplot(dfcTAC_plot, aes(x=intercept_hcr,y=slope_hcr,z=cTAC_nat2))
		pctn2 <- pctn2 +geom_raster(aes(fill=cTAC_nat2))
		pctn2 <- pctn2 + scale_fill_gradientn(colours = myPalette(4),name  ="% of TAC caught")
		pctn2 <- pctn2 + theme_bw(16) + theme(legend.position="top")
		pctn2 <- pctn2 + ylab("") + xlab("Intercept")
		pctn2 <- pctn2 + ggtitle("Nation 2 ")
		pctn2 <- pctn2 + theme(plot.title = element_text(hjust = 0.5))
		pctn2 <- pctn2 + guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                label.position = "bottom"))
		pctn2 <- pctn2 + geom_text(aes(x=dfcTAC_plot$intercept_hcr[which.max(dfcTAC_plot$cTAC_nat2)],y=dfcTAC_plot$slope_hcr[which.max(dfcTAC_plot$cTAC_nat2)],label="max"),fontface = "bold",colour="gray90")
		pctn2 <- pctn2 + geom_point(aes(x=dfcTAC$intercept_hcr[cTAC_nat2_comparison],y=dfcTAC$slope_hcr[cTAC_nat2_comparison]),colour="gray90", size=10)
		pctn2 <- pctn2 + geom_text(aes(x=dfcTAC$intercept_hcr[cTAC_nat2_comparison],y=dfcTAC$slope_hcr[cTAC_nat2_comparison],label="10:40"),fontface = "bold",colour="black")

		pctn2


		pctn<-plot_grid(pctn1,pctn2,pctt,ncol=3)


		#plun <- plun+geom_contour()
		print(pctn)

		if(sv==TRUE){
			setwd("/Users/catarinawor/Documents/Lagrangian/report/HCR_coarse")
			ggsave(paste(nome,"_pTAC.pdf",sep=""), plot=pctn, width = 15, height = 5)
		}

	}
	


		

	
}

