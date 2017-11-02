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


plot_logUtility <- function( M , Msq, sv=F, nome="",nations=F)
{
	cat("plot_logUtility")


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

	
	myPalette <- colorRampPalette(rev(brewer.pal(8, "Spectral")))

	

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
		
		#tlabs<-round(range(dflu_plot$utility_total)+c(2,-2))#+(c(tmp1,tmp2)*c(.1,-.1))
		#plut <- ggplot(dflu_plot, aes(x=intercept_hcr,y=slope_hcr,z=utility_total))
		#plut <- plut +geom_raster(aes(fill=utility_total))
		#plut  <- plut  + scale_fill_gradientn(colours = myPalette(6),name  ="log utility",breaks=tlabs,labels=(tlabs))
		#plut  <- plut  + theme_bw(16) + theme(legend.position="top")
		#plut  <- plut  + ylab("") + xlab("Intercept")
		#plut  <- plut  + ggtitle("Total sum")
		#plut  <- plut  +  theme(plot.title = element_text(hjust = 0.5))
		#plut  <- plut  + geom_text(aes(x=intercept_hcr[which.max(utility_total)],y=slope_hcr[which.max(utility_total)],label=round(max(dflu_plot$utility_total)/1000,2)),fontface = "bold",colour="gray90")
		#plut  <- plut  + geom_text(aes(x=intercept_hcr[which.max(utility_total)],y=slope_hcr[which.max(utility_total)],label="max"),fontface = "bold",colour="gray90")
		#plut  <- plut  + geom_point(aes(x=intercept_hcr[utility_total_comparison],y=slope_hcr[utility_total_comparison]),colour="gray90", size=10)
		#plut  <- plut  + geom_text(aes(x=intercept_hcr[utility_total_comparison],y=slope_hcr[utility_total_comparison],label=round(dflu_plotsq$utility_total/1000,2)),fontface = "bold",colour="black")		
		#plut  <- plut  + geom_text(aes(x=intercept_hcr[utility_total_comparison],y=slope_hcr[utility_total_comparison],label="10:40"),fontface = "bold",colour="black")		
		#plut  <- plut  + theme(legend.position = "none") 
		#plut


		tmp1<-range(dflu_plot$utility_total2)[2]-range(dflu_plot$utility_total2)[1]	

		tlabs<-round(range(dflu_plot$utility_total2))+(tmp1*c(.1,-.1))
		plut2 <- ggplot(dflu_plot, aes(x=intercept_hcr,y=slope_hcr,z=utility_total2))
		plut2 <- plut2 +geom_raster(aes(fill=utility_total2))
		plut2  <- plut2  + scale_fill_gradientn(colours = myPalette(6),name ="log utility",breaks=tlabs,labels=c("low","high"))
		plut2  <- plut2  + theme_bw(16) + theme(legend.position="right")
		plut2  <- plut2  + ylab("") + xlab("Intercept")
		plut2  <- plut2  + ggtitle("Total")
		plut2  <- plut2  + guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = 0.1,title.vjust = 10
                               ))
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
		plun1 <- plun1  + theme(legend.position = "none") 
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
		plun2  <- plun2  + theme(legend.position = "none") 
		plun2

		plun<-plot_grid(plun1,plun2,plut2,ncol=3,rel_widths = c(1,1, 1.2))
		print(plun)

		if(sv==TRUE){
			setwd("/Users/catarinawor/Documents/Lagrangian/report/HCR_coarse")
			ggsave(paste(nome,"logUtility.pdf",sep=""), plot=plun, width = 15, height = 5)
		}

		return(plun)
	}
	
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











	

