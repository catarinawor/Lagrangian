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


plot_B40 <- function( M,  Msq, sv=F, nome="",nations=F, limite=.4)
{
	cat("plot_Bnat_TAC")


	n <- length( M )

	#M<-SIMSdat[[1]]
	#Msq<-SIMSdat[[4]]

	slope_hcr<-NULL
	intercept_hcr<-NULL
	B40_total<-NULL



	for(i in 1:n){

		iniyr<-M[[i]][[1]]$"nyr"+1
		fimyr<-M[[i]][[1]]$"proj_yr"

		slope_hcr[i]<-M[[i]][[1]]$slope_hcr
		intercept_hcr[i]<-M[[i]][[1]]$intercept_hcr
		B40_total[i]<-sum((M[[i]][[1]]$ytB/M[[i]][[1]]$Bo)[iniyr:fimyr]<limite)/length(iniyr:fimyr)

	}


	slope_hcr_sq<-NULL
	intercept_hcr_sq<-NULL
	B40_total_sq<-NULL

	nsq<- length( Msq )

	

	for(y in 1:nsq){

		slope_hcr_sq[y]<-Msq[[y]][[1]]$slope_hcr
		intercept_hcr_sq[y]<-Msq[[y]][[1]]$intercept_hcr
		B40_total_sq[y]<-sum((Msq[[y]][[1]]$ytB/Msq[[y]][[1]]$Bo)[iniyr:fimyr]<limite)/length(iniyr:fimyr)

	}

	




	db40<-data.frame(slope_hcr=slope_hcr,intercept_hcr=intercept_hcr,
		b40=B40_total)

	db40sq<-data.frame(slope_hcr=slope_hcr_sq,intercept_hcr=intercept_hcr_sq,
		b40=B40_total_sq)


	db40_plot<-aggregate(db40,list(db40$slope_hcr,db40$intercept_hcr), mean)

	db40_plotsq<-aggregate(db40sq,list(db40sq$slope_hcr,db40sq$intercept_hcr), mean)


	b40_total_comparison<-which.min(abs(db40_plot$b40-db40_plotsq$b40))
	

	
	
	
	myPalette <- colorRampPalette(rev(brewer.pal(8, "Spectral")))
	#myPalette <- colorRampPalette((brewer.pal(8, "Greys")))
	

	if(nations==F){
		plu <- ggplot(dflu_plot, aes(x=intercept_hcr,y=slope_hcr,z=utility_total))
		plu <- plu+geom_raster(aes(fill=utility_total))
		plu <- plu + scale_fill_gradientn(colours = myPalette(4))
		plu <- plu + theme_bw(16)
		print(plu)
	}else{
		#
	
		tlabs<- round(range(db40_plot$b40)*c(1.2,0.95),2)
		pB40t <- ggplot(db40_plot, aes(x=intercept_hcr,y=slope_hcr,z=b40))
		pB40t <- pB40t +geom_raster(aes(fill=b40))
		pB40t <- pB40t  + scale_fill_gradientn(colours = myPalette(11),name  =paste("% time Bt <",limite,"Bo"),breaks=tlabs,labels=(tlabs))
		pB40t <- pB40t  + theme_bw(16) + theme(legend.position="top")
		pB40t <- pB40t  + ylab("Harvest rate") + xlab("Biomass threshold")
		pB40t <- pB40t  + ggtitle("Total")
		pB40t <- pB40t  +  theme(plot.title = element_text(hjust = 0.5))
		pB40t <- pB40t  + guides(fill = guide_colourbar(title.position = "top",
                                title.hjust = .5,
                                label.position = "bottom"))
		pB40t <- pB40t + geom_point(aes(x=db40_plot$intercept_hcr[b40_total_comparison],y=db40_plot$slope_hcr[b40_total_comparison]),colour="gray90", size=10)
		pB40t <- pB40t + geom_text(aes(x=db40_plot$intercept_hcr[b40_total_comparison],y=db40_plot$slope_hcr[b40_total_comparison],label="40:10"),fontface = "bold",colour="black")
		pB40t <- pB40t + geom_text(aes(x=db40_plot$intercept_hcr[which.min(db40_plot$b40)],y=db40_plot$slope_hcr[which.min(db40_plot$b40)],label="min"),fontface = "bold",colour="black", position = position_nudge(y = -0.02))
		#


		print(pB40t)

		if(sv==TRUE){
			setwd("/Users/catarinawor/Documents/Lagrangian/report/HCR_coarse")
			ggsave(paste(nome,"B40.pdf",sep=""), plot=pB40t, width = 5, height = 5)
		}

		return(pB40t)
	}
	
}





#plot_B40_allscn <- function( M, M2,M3,Msq,Msq2,Msq3, sv=F, nome="",nations=F)
#{
#	cat("plot_Bnat_TAC")
#
#	Ms<-list(M, M2,M3)
#
#	ns <- length( M )+length( M2 )+length( M3 )
#
#	#M<-SIMSdat[[1]]
#	#Msq<-SIMSdat[[4]]
#
#	slope_hcr<-NULL
#	intercept_hcr<-NULL
#	B40_total<-NULL
#	scenario<-c(rep(1,length( M )),rep(2,length( M2 )),rep(3,length( M3 )))
#
#	i=1
#	for(y in 1:length(Ms)){
#
#		M<-Ms[[y]]
#		n <- length( M )
#
#		for(a in 1:n){
#
#		
#
#		iniyr<-M[[a]][[1]]$"nyr"+1
#		fimyr<-M[[a]][[1]]$"proj_yr"
#
#		slope_hcr[i]<-M[[a]][[1]]$slope_hcr
#		intercept_hcr[i]<-M[[a]][[1]]$intercept_hcr
#		B40_total[i]<-sum((M[[a]][[1]]$ytB/M[[a]][[1]]$Bo)[iniyr:fimyr]>0.4)/length(iniyr:fimyr)
#		i=i+1
#	}
#}
#
#
#
#
#	slope_hcr_sq<-NULL
#	intercept_hcr_sq<-NULL
#	B40_total_sq<-NULL
#
#	nsqs<- length( Msq ) + length( Msq2 ) +length( Msq3 )
#
#	scenariosq<-c(rep(1,length( Msq )),rep(2,length( Msq2 )),rep(3,length( Msq3 )))
#	Msqs<-list(Msq, Msq2,Msq3)
#	ii=0
#	for(yy in 1:length(Msqs)){
#
#		Msq<-Msqs[[yy]]
#		nsq <- length( Msq )
#
#
#	for(aa in 1:nsq){
#		ii=ii+1
#
#		slope_hcr_sq[ii]<-Msq[[aa]][[1]]$slope_hcr
#		intercept_hcr_sq[ii]<-Msq[[aa]][[1]]$intercept_hcr
#		B40_total_sq[ii]<-sum((Msq[[aa]][[1]]$ytB/Msq[[aa]][[1]]$Bo)[iniyr:fimyr]>0.4)/length(iniyr:fimyr)
#
#	}
#}
#
#	length(slope_hcr)
#	length(scenario)
#
#	db40<-data.frame(slope_hcr=slope_hcr,intercept_hcr=intercept_hcr,
#		b40=B40_total, scenario=scenario)
#
#	db40sq<-data.frame(slope_hcr=slope_hcr_sq,intercept_hcr=intercept_hcr_sq,
#		b40=B40_total_sq,scenario=scenariosq)
#
#
#	db40_plot<-aggregate(db40,list(db40$slope_hcr,db40$intercept_hcr,db40$scenario), mean)
#
#	db40_plotsq<-aggregate(db40sq,list(db40sq$slope_hcr,db40sq$intercept_hcr,db40sq$scenario), mean)
#
#
#	b40_total_comparison<-NULL
#
#	for( sc in 1:length(unique(db40sq$scenario))){
#
#		b40_total_comparison[sc]<-which.min(abs(db40_plot$b40[db40_plot$scenario==sc]-db40_plotsq$b40[db40_plotsq$scenario==sc]))
#	}
#
#	df<-data.frame(x=db40_plot$intercept_hcr[b40_total_comparison],y=db40_plot$slope_hcr[b40_total_comparison],b40_total_comparison=b40_total_comparison, scenario=1:3,
#		xmax=c(db40_plot$intercept_hcr[which.max(db40_plot$b40[db40_plot$scenario==1])],
#			db40_plot$intercept_hcr[which.max(db40_plot$b40[db40_plot$scenario==2])],
#			db40_plot$intercept_hcr[which.max(db40_plot$b40[db40_plot$scenario==3])])
#		 ymax=c(db40_plot$slope_hcr[which.max(db40_plot$b40[db40_plot$scenario==1])],
#		 	db40_plot$slope_hcr[which.max(db40_plot$b40[db40_plot$scenario==2])],
#		 db40_plot$slope_hcr[which.max(db40_plot$b40[db40_plot$scenario==3])]))
#
#	
#	
#	
#	myPalette <- colorRampPalette(rev(brewer.pal(8, "Spectral")))
#	#myPalette <- colorRampPalette((brewer.pal(8, "Greys")))
#	
#
#	if(nations==F){
#		plu <- ggplot(dflu_plot, aes(x=intercept_hcr,y=slope_hcr,z=utility_total))
#		plu <- plu+geom_raster(aes(fill=utility_total))
#		plu <- plu + scale_fill_gradientn(colours = myPalette(4))
#		plu <- plu + theme_bw(16)
#		print(plu)
#	}else{
#		#
#	
#		tlabs<- round(range(db40_plot$b40),2)
#		pB40t <- ggplot(db40_plot, aes(x=intercept_hcr,y=slope_hcr,z=b40))
#		pB40t <- pB40t +geom_raster(aes(fill=b40))
#		pB40t <- pB40t + facet_wrap(~scenario)
#		pB40t <- pB40t  + scale_fill_gradientn(colours = myPalette(11),name  ="% time Bt > 0.4 Bo",breaks=tlabs,labels=(tlabs))
#		pB40t <- pB40t  + theme_bw(16) + theme(legend.position="top")
#		pB40t <- pB40t  + ylab("") + xlab("Intercept")
#		pB40t <- pB40t  + ggtitle("Total")
#		pB40t <- pB40t  +  theme(plot.title = element_text(hjust = 0.5))
#		pB40t <- pB40t  + guides(fill = guide_colourbar(title.position = "top",
#                                title.hjust = .5,
#                                label.position = "bottom"))
#		#pB40t <- pB40t + geom_point(data=df,aes(x=x,y=y),colour="gray90", size=10)
#		#pB40t <- pB40t + geom_text(aes(x=db40_plot$intercept_hcr[b40_total_comparison],y=db40_plot$slope_hcr[b40_total_comparison],label="10:40"),fontface = "bold",colour="black")
#		#pB40t <- pB40t + geom_text(aes(x=db40_plot$intercept_hcr[which.max(db40_plot$b40)],y=db40_plot$slope_hcr[which.max(db40_plot$b40)],label="max"),fontface = "bold",colour="black", position = position_nudge(y = -0.02))
#		#
#
#
#		print(pB40t)
#
#		if(sv==TRUE){
#			setwd("/Users/catarinawor/Documents/Lagrangian/report/HCR_coarse")
#			ggsave(paste(nome,"B40.pdf",sep=""), plot=pB40t, width = 5, height = 5)
#		}
#
#		return(pB40t)
#	}
#	
#}


	

