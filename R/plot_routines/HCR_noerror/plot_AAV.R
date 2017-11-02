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


plot_AAV <- function( M , Msq, sv=F, nome="",nations=F)
{
	cat("plot_AAV")


	n <- length( M )


	slope_hcr<-NULL
	intercept_hcr<-NULL
	AAV_total<-NULL
	AAV_nat1<-NULL
	AAV_nat2<-NULL


	for(i in 1:n){

		iniyr<-M[[i]][[1]]$"nyr"+1
		fimyr<-M[[i]][[1]]$"proj_yr"

		slope_hcr[i]<-M[[i]][[1]]$slope_hcr
		intercept_hcr[i]<-M[[i]][[1]]$intercept_hcr
		AAV_total[i]<-mean(abs(apply(M[[i]][[1]]$yYieldNat[iniyr:(fimyr-1),],1,sum)-
							apply(M[[i]][[1]]$yYieldNat[(iniyr+1):fimyr,],1,sum)))/mean(
							(apply(M[[i]][[1]]$yYieldNat[iniyr:(fimyr-1),],1,sum)+
							apply(M[[i]][[1]]$yYieldNat[(iniyr+1):fimyr,],1,sum)))
		
		AAV_nat1[i]<-mean(abs(M[[i]][[1]]$yYieldNat[iniyr:(fimyr-1),1]-
							M[[i]][[1]]$yYieldNat[(iniyr+1):fimyr,1]))/mean(
							(M[[i]][[1]]$yYieldNat[iniyr:(fimyr-1),1]+
							M[[i]][[1]]$yYieldNat[(iniyr+1):fimyr,1]))


		AAV_nat2[i]<-mean(abs(M[[i]][[1]]$yYieldNat[iniyr:(fimyr-1),2]-
						 M[[i]][[1]]$yYieldNat[(iniyr+1):fimyr,2]))/mean(
							(M[[i]][[1]]$yYieldNat[iniyr:(fimyr-1),2]+
							M[[i]][[1]]$yYieldNat[(iniyr+1):fimyr,2]))


	}


	nsq <- length( Msq )


	slope_hcr_sq<-NULL
	intercept_hcr_sq<-NULL
	AAV_total_sq<-NULL
	AAV_nat1_sq<-NULL
	AAV_nat2_sq<-NULL


	for(i in 1:nsq){

		iniyr<-Msq[[i]][[1]]$"nyr"+1
		fimyr<-Msq[[i]][[1]]$"proj_yr"

		slope_hcr_sq[i]<-Msq[[i]][[1]]$slope_hcr
		intercept_hcr_sq[i]<-Msq[[i]][[1]]$intercept_hcr
		AAV_total_sq[i]<-mean(abs(apply(Msq[[i]][[1]]$yYieldNat[iniyr:(fimyr-1),],1,sum)-
							apply(Msq[[i]][[1]]$yYieldNat[(iniyr+1):fimyr,],1,sum)))/mean(
							(apply(Msq[[i]][[1]]$yYieldNat[iniyr:(fimyr-1),],1,sum)+
							apply(Msq[[i]][[1]]$yYieldNat[(iniyr+1):fimyr,],1,sum)))
		
		AAV_nat1_sq[i]<-mean(abs(Msq[[i]][[1]]$yYieldNat[iniyr:(fimyr-1),1]-
							Msq[[i]][[1]]$yYieldNat[(iniyr+1):fimyr,1]))/mean(
							(Msq[[i]][[1]]$yYieldNat[iniyr:(fimyr-1),1]+
							Msq[[i]][[1]]$yYieldNat[(iniyr+1):fimyr,1]))


		AAV_nat2_sq[i]<-mean(abs(Msq[[i]][[1]]$yYieldNat[iniyr:(fimyr-1),2]-
						 Msq[[i]][[1]]$yYieldNat[(iniyr+1):fimyr,2]))/mean(
							(Msq[[i]][[1]]$yYieldNat[iniyr:(fimyr-1),2]+
							Msq[[i]][[1]]$yYieldNat[(iniyr+1):fimyr,2]))


	}	

	

	dfaav<-data.frame(slope_hcr=slope_hcr,intercept_hcr=intercept_hcr,
		AAV_total=AAV_total,AAV_nat1=AAV_nat1,AAV_nat2=AAV_nat2)

	dfaavsq<-data.frame(slope_hcr=slope_hcr_sq,intercept_hcr=intercept_hcr_sq,
		AAV_total=AAV_total_sq,AAV_nat1=AAV_nat1_sq,AAV_nat2=AAV_nat2_sq)
	


	dfaav_plot<-aggregate(dfaav,list(dfaav$slope_hcr,dfaav$intercept_hcr), mean)

	dfaav_plotsq<-aggregate(dfaavsq,list(dfaavsq$slope_hcr,dfaavsq$intercept_hcr), mean)

	AAV_total_comparison<-which.min(abs(dfaav_plot$AAV_total-dfaav_plotsq$AAV_total))
	AAV_nat1_comparison<-which.min(abs(dfaav_plot$AAV_nat1-dfaav_plotsq$AAV_nat1))
	AAV_nat2_comparison<-which.min(abs(dfaav_plot$AAV_nat2-dfaav_plotsq$AAV_nat2))


	summary(dfaav_plot)

	dfaav_plotn<-data.frame(slope_hcr=rep(dfaav_plot$slope_hcr,3),intercept_hcr=rep(dfaav_plot$intercept_hcr,3),
		AAV=c(dfaav_plot$AAV_total,dfaav_plot$AAV_nat1,dfaav_plot$AAV_nat2),nation=rep(c("Total", "Nation 1", "Nation 2"),each=length(dfaav_plot$intercept_hcr)))

	summary(dfaav_plotn)
	myPalette <-  colorRampPalette(rev(brewer.pal(8, "Spectral")))


	if(nations==F){
		pcl <- ggplot(dfaav_plot, aes(x=intercept_hcr,y=slope_hcr,z=AAV_total))
		pcl <- pcl+geom_raster(aes(fill=AAV_total))
		pcl <- pcl + scale_fill_gradientn(colours = myPalette(4))
		pcl <- pcl + theme_bw(16)
		print(pcl)
	}else{
		#paavn <- ggplot(dfaav_plotn, aes(x=intercept_hcr,y=slope_hcr,z=AAV))
		#paavn <- paavn +geom_raster(aes(fill=AAV))
		#paavn <- paavn+facet_wrap(~nation)
		#paavn <- paavn + scale_fill_gradientn(colours = myPalette(4),name  ="Median AAV")
		#paavn <- paavn + theme_bw(16)
		#paavn  <- paavn  + ylab("Slope") + xlab("Intercept")

		paavt <- ggplot(dfaav_plot, aes(x=intercept_hcr,y=slope_hcr,z=AAV_total))
		paavt <- paavt +geom_raster(aes(fill=AAV_total))
		paavt <- paavt + scale_fill_gradientn(colours = myPalette(4),name  ="Median AAV")
		paavt <- paavt + theme_bw(16)+ theme(legend.position="right")
		paavt <- paavt + ylab("") + xlab("Intercept")
		paavt <- paavt + ggtitle("Total")
		paavt <- paavt + theme(plot.title = element_text(hjust = 0.5))
		paavt <- paavt + guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5))
		paavt  <- paavt + geom_point(aes(x=dfaav_plot$intercept_hcr[AAV_total_comparison],y=dfaav_plot$slope_hcr[AAV_total_comparison]),colour="gray90", size=10)
		paavt  <- paavt + geom_text(aes(x=dfaav_plot$intercept_hcr[AAV_total_comparison],y=dfaav_plot$slope_hcr[AAV_total_comparison],label="10:40"),fontface = "bold",colour="black")
		paavt  <- paavt + geom_text(aes(x=dfaav_plot$intercept_hcr[which.min(dfaav_plot$AAV_total)],y=dfaav_plot$slope_hcr[which.min(dfaav_plot$AAV_total)],label="min"),fontface = "bold",colour="black", position = position_nudge(y = -0.02))
		paavt


		paavn1 <- ggplot(dfaav_plot, aes(x=intercept_hcr,y=slope_hcr,z=AAV_nat1))
		paavn1 <- paavn1 +geom_raster(aes(fill=AAV_nat1))
		paavn1 <- paavn1 + scale_fill_gradientn(colours = myPalette(4),name  ="Median AAV")
		paavn1 <- paavn1 + theme_bw(16)+ theme(legend.position="none")#egend.position = "top")e
		paavn1 <- paavn1 + ylab("Slope") + xlab("Intercept")
		paavn1 <- paavn1  + ggtitle("Nation 1")
		paavn1 <- paavn1  +  theme(plot.title = element_text(hjust = 0.5))
		paavn1 <- paavn1 + guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5))
                                #label.position = "bottom"))
		paavn1  <- paavn1 + geom_point(aes(x=dfaav_plot$intercept_hcr[AAV_nat1_comparison],y=dfaav_plot$slope_hcr[AAV_nat1_comparison]),colour="gray90", size=10)
		paavn1  <- paavn1 + geom_text(aes(x=dfaav_plot$intercept_hcr[AAV_nat1_comparison],y=dfaav_plot$slope_hcr[AAV_nat1_comparison],label="10:40"),fontface = "bold",colour="black")
		paavn1  <- paavn1 + geom_text(aes(x=dfaav_plot$intercept_hcr[which.min(dfaav_plot$AAV_nat1)],y=dfaav_plot$slope_hcr[which.min(dfaav_plot$AAV_nat1)],label="min"),fontface = "bold",colour="black", position = position_nudge(y = -0.02))
		paavn1

		paavn2 <- ggplot(dfaav_plot, aes(x=intercept_hcr,y=slope_hcr,z=AAV_nat2))
		paavn2 <- paavn2 +geom_raster(aes(fill=AAV_nat2))
		paavn2 <- paavn2 + scale_fill_gradientn(colours = myPalette(4),name  ="Median AAV")
		paavn2 <- paavn2 + theme_bw(16) +theme(legend.position="none")
		paavn2 <- paavn2 + ylab("") + xlab("Intercept")
		paavn2 <- paavn2  + ggtitle("Nation 2")
		paavn2 <- paavn2  +  theme(plot.title = element_text(hjust = 0.5))
		paavn2 <- paavn2 + guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5
                               ))
		paavn2  <- paavn2 + geom_point(aes(x=dfaav_plot$intercept_hcr[AAV_nat2_comparison],y=dfaav_plot$slope_hcr[AAV_nat2_comparison]),colour="gray90", size=10)
		paavn2  <- paavn2 + geom_text(aes(x=dfaav_plot$intercept_hcr[AAV_nat2_comparison],y=dfaav_plot$slope_hcr[AAV_nat2_comparison],label="10:40"),fontface = "bold",colour="black")
		paavn2  <- paavn2 + geom_text(aes(x=dfaav_plot$intercept_hcr[which.min(dfaav_plot$AAV_nat2)],y=dfaav_plot$slope_hcr[which.min(dfaav_plot$AAV_nat2)],label="min"),fontface = "bold",colour="black", position = position_nudge(y = -0.02))
	
		paavn2

		#plun <- plun+geom_contour()
		paavn<-plot_grid(paavn1,paavn2,paavt,ncol=3,rel_widths = c(1,1, 1.3))
		print(paavn)

		if(sv==TRUE){
			setwd("/Users/catarinawor/Documents/Lagrangian/report/HCR_coarse")
			ggsave(paste(nome,"AAV.pdf",sep=""), plot=paavn, width = 15, height = 5)
		}

		return(paavn)

	}
	
}

