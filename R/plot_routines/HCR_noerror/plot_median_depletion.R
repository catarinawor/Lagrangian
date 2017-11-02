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


plot_medianDepl <- function( M , sv=F, nome="",nations=F)
{
	cat("plot_logUtility")


	n <- length( M )


	slope_hcr<-NULL
	intercept_hcr<-NULL


	depletion_total<-NULL
	

	for(i in 1:n){

		iniyr<-M[[i]][[1]]$"nyr"+1
		fimyr<-M[[i]][[1]]$"proj_yr"

		slope_hcr[i]<-M[[i]][[1]]$slope_hcr
		intercept_hcr[i]<-M[[i]][[1]]$intercept_hcr
		depletion_total[i]<-median(M[[i]][[1]]$ytB[iniyr:fimyr]/M[[i]][[1]]$Bo)
		

	}


	dmd<-data.frame(slope_hcr=slope_hcr,intercept_hcr=intercept_hcr,
		depletion_total=depletion_total)


	dmd_plot<-aggregate(dmd,list(dmd$slope_hcr,dmd$intercept_hcr), mean)

	myPalette <- colorRampPalette(rev(brewer.pal(6, "Spectral")))
	
	pmd <- ggplot(dmd_plot, aes(x=intercept_hcr,y=slope_hcr,z=depletion_total))
	pmd <- pmd + geom_raster(aes(fill=depletion_total))
	pmd <- pmd + scale_fill_gradientn(colours = myPalette(4))
	pmd <- pmd + theme_bw(16)
	print(pmd)

	if(sv==TRUE){
		setwd("/Users/catarinawor/Documents/Lagrangian/report/HCR_coarse")
		ggsave(paste(nome,"logUtility.pdf",sep=""), plot=plun, width = 15, height = 5)
	}

	
}

