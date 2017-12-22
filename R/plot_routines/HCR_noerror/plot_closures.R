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


plot_closures <- function( M , Msq, sv=F, nome="",nations=F)
{
	cat("plot_closures")


	n <- length( M )


	slope_hcr<-NULL
	intercept_hcr<-NULL
	closure_total<-NULL
	closure_nat1<-NULL
	closure_nat2<-NULL


	for(i in 1:n){

		slope_hcr[i]<-M[[i]][[1]]$slope_hcr
		intercept_hcr[i]<-M[[i]][[1]]$intercept_hcr
		closure_total[i]<-sum(apply(M[[i]][[1]]$yYieldNat[101:160,],1,sum)==0)/length(101:160)
		closure_nat1[i]<-sum(M[[i]][[1]]$yYieldNat[101:160,1]==0)/length(101:160)
		closure_nat2[i]<-sum(M[[i]][[1]]$yYieldNat[101:160,2]==0)/length(101:160)


	}

	nsq <- length( Msq )

	slope_hcr_sq<-NULL
	intercept_hcr_sq<-NULL
	closure_total_sq<-NULL
	closure_nat1_sq<-NULL
	closure_nat2_sq<-NULL


	for(y in 1:nsq){

		slope_hcr_sq[y]<-Msq[[y]][[1]]$slope_hcr
		intercept_hcr_sq[y]<-Msq[[y]][[1]]$intercept_hcr
		closure_total_sq[y]<-sum(apply(Msq[[y]][[1]]$yYieldNat[101:160,],1,sum)==0)/length(101:160)
		closure_nat1_sq[y]<-sum(Msq[[y]][[1]]$yYieldNat[101:160,1]==0)/length(101:160)
		closure_nat2_sq[y]<-sum(Msq[[y]][[1]]$yYieldNat[101:160,2]==0)/length(101:160)


	}

	

	dfcl<-data.frame(slope_hcr=slope_hcr,intercept_hcr=intercept_hcr,
		closure_total=closure_total,closure_nat1=closure_nat1,closure_nat2=closure_nat2)

	dfclsq<-data.frame(slope_hcr=slope_hcr_sq,intercept_hcr=intercept_hcr_sq,
		closure_total=closure_total_sq,closure_nat1=closure_nat1_sq,closure_nat2=closure_nat2_sq)



	dfcl_plot<-aggregate(dfcl,list(dfcl$slope_hcr,dfcl$intercept_hcr), mean)
	dfcl_plotsq<-aggregate(dfclsq,list(dfclsq$slope_hcr,dfclsq$intercept_hcr), mean)
	
	closure_total_comparison<-which.min(abs(dfcl_plot$closure_total-dfcl_plotsq$closure_total))
	closure_nat1_comparison<-which.min(abs(dfcl_plot$closure_nat1-dfcl_plotsq$closure_nat1))
	closure_nat2_comparison<-which.min(abs(dfcl_plot$closure_nat2-dfcl_plotsq$closure_nat2))



	summary(dfcl_plot)

	dfcl_plotn<-data.frame(slope_hcr=rep(dfcl_plot$slope_hcr,3),intercept_hcr=rep(dfcl_plot$intercept_hcr,3),
		closure=c(dfcl_plot$closure_total,dfcl_plot$closure_nat1,dfcl_plot$closure_nat2),
		nation=rep(c("Total", "Nation 1", "Nation 2"),each=length(dfcl_plot$intercept_hcr)))

	summary(dfcl_plotn)
	myPalette <- colorRampPalette(rev(brewer.pal(8, "Spectral")))

	dfcl_plot2<-dfcl_plot

	if(nations==F){
		pcl <- ggplot(dfcl_plot, aes(x=intercept_hcr,y=slope_hcr,z=closure_total))
		pcl <- pcl+geom_raster(aes(fill=closure_total))
		pcl <- pcl + scale_fill_gradientn(colours = myPalette(4))
		pcl <- pcl + theme_bw(16)
		print(pcl)
	}else{
		#pcln <- ggplot(dfcl_plotn, aes(x=intercept_hcr,y=slope_hcr,z=closure))
		#pcln <- pcln +geom_raster(aes(fill=closure))
		#pcln <- pcln+facet_wrap(~nation)
		#pcln  <- pcln  + scale_fill_gradientn(colours = myPalette(8),name  ="Closure proportion")
		#pcln  <- pcln  + theme_bw(16)
		#pcln  <- pcln  + ylab("Slope") + xlab("Intercept")

		tlabs<-round(range(dfcl_plot2$closure_total)*0.9,4)
		pclt <- ggplot(dfcl_plot2, aes(x=intercept_hcr,y=slope_hcr,z=closure_total))
		pclt <- pclt + geom_raster(aes(fill=closure_total))
		pclt <- pclt + scale_fill_gradientn(colours = myPalette(8),name  ="Closure rate",breaks=tlabs,labels=tlabs)
		pclt <- pclt + theme_bw(16) + theme(legend.position="top")
		pclt <- pclt + ylab("Harvest rate") + xlab("Biomass threshold")
		pclt <- pclt  + ggtitle("Total")
		pclt <- pclt  +  theme(plot.title = element_text(hjust = 0.5))
		pclt <- pclt  + guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5))
		pclt <- pclt + geom_point(aes(x=dfcl_plot2$intercept_hcr[which.min(dfcl_plot2$closure_total)],y=dfcl_plot2$slope_hcr[which.min(dfcl_plot2$closure_total)]),colour="black", size=10,shape=15)
		pclt <- pclt + geom_point(aes(x=dfcl_plot2$intercept_hcr[closure_total_comparison],y=dfcl_plot2$slope_hcr[closure_total_comparison]),colour="gray90", size=10)
		pclt <- pclt + geom_text(aes(x=dfcl_plot2$intercept_hcr[closure_total_comparison],y=dfcl_plot2$slope_hcr[closure_total_comparison],label="40:10"),fontface = "bold",colour="black")
		#pclt <- pclt + geom_text(aes(x=dfcl_plot2$intercept_hcr[which.min(dfcl_plot2$closure_total)],y=dfcl_plot2$slope_hcr[which.min(dfcl_plot2$closure_total)],label="min"),fontface = "bold",colour="gray90",position = position_nudge(y = -0.02))
		
		pclt 

		pcln1 <- ggplot(dfcl_plot, aes(x=intercept_hcr,y=slope_hcr,z=closure_nat1))
		pcln1 <- pcln1+geom_raster(aes(fill=closure_nat1))
		pcln1 <- pcln1 + scale_fill_gradientn(colours = myPalette(8),name  ="Closure prop.")
		pcln1 <- pcln1 + theme_bw(16) + theme(legend.position = "none") 
		pcln1 <- pcln1 + ylab("Harvest rate") + xlab("Biomass threshold")
		pcln1 <- pcln1  + ggtitle("Nation 1")
		pcln1 <- pcln1  +  theme(plot.title = element_text(hjust = 0.5))
		pcln1 

		pcln2 <- ggplot(dfcl_plot, aes(x=intercept_hcr,y=slope_hcr,z=closure_nat2))
		pcln2 <- pcln2+geom_raster(aes(fill=closure_nat2))
		pcln2 <- pcln2 + scale_fill_gradientn(colours = myPalette(8),name  ="Closure prop.")
		pcln2 <- pcln2 + theme_bw(16) + theme(legend.position = "none") 
		pcln2 <- pcln2 + ylab("") + xlab("Biomass threshold")
		pcln2 <- pcln2  + ggtitle("Nation 2")
		pcln2 <- pcln2  +  theme(plot.title = element_text(hjust = 0.5))
		pcln2 


		pcln<-plot_grid(pcln1,pcln2,pclt,ncol=3)

		#plun <- plun+geom_contour()
		print(pclt)

		if(sv==TRUE){
			setwd("/Users/catarinawor/Documents/Lagrangian/report/HCR_coarse")
			ggsave(paste(nome,"closure.pdf",sep=""), plot=pcln, width = 5, height = 5)
		}
	return(pclt)

	}
	

}

