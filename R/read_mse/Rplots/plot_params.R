#==============================================
#Title:plot_params
#Author: Catarina Wor
#date: Oct 21th 2016
#Function to plot true biomass fro lagrangian model.extract
# code adapted from iSCAM stuff (Martell et al)
#==============================================
#for testing delete, when done
M<-SAdat

length(M)

source("/Users/catarinawor/Documents/Lagrangian/R/read_mse/Rplots/calc_quantile.R")


require(reshape2)
require(tidyr)

plot_h_Ro <- function( M , OM)
{
	cat("plot_h_Ro")

	n <- length( M )
	mdf <- NULL


	for(i in 1:n)
	{
		nn<- length( M[[i]] )
		for(y in 1:nn){
			df <- data.frame(Model=names(M)[i], Year=M[[i]][[y]]$nyr, h=M[[i]][[y]]$steepness,
			Ro=M[[i]][[y]]$ro)
			mdf <- rbind(mdf,df)
		}
		
	}
	
	hdf<-spread(mdf[,-4],key="Model", "h")

	qh<-apply(hdf,1,calc_quantile)
	medianh<-qh["50%",]
	lowh<-qh["2.5%",]
	highh<-qh["97.5%",]


	rodf<-spread(mdf[,-3],key="Model", "Ro")
	qro<-apply(rodf[,-1],1,calc_quantile)
	medianro<-qro["50%",]
	lowro<-qro["2.5%",]
	highro<-qro["97.5%",]
	
	

	adf<-rbind(data.frame(Median=medianro,low=lowro, high=highro,Year=rodf$"Year",parameter=rep("Ro",length(medianro))),
		data.frame(Median=medianh,low=lowh, high=highh,Year=hdf$"Year",parameter=rep("h",length(medianh))))

	p <- ggplot(adf,aes(x=Year,y=Median)) + geom_line()
	p <- p + geom_ribbon(aes(ymax=high, ymin=low),alpha=0.2)
	p <- p + facet_wrap(~parameter,scales="free")
	p <- p + labs(x="Year",y=parameter)
	p <- p + theme_bw(11)
	print(p)

	
}