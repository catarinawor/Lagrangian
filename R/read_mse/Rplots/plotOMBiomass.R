#==============================================
#Title:plotOMBiomass
#Author: Catarina Wor
#date: Oct 21th 2016
#Function to plot true biomass fro lagrangian model.extract
# code adapted from iSCAM stuff (Martell et al)
#==============================================

source("/Users/catarinawor/Documents/Lagrangian/R/Rplots/calc_quantile.R")
require(ggplot2)
plotOMBiomass <- function( M )
{
	cat("plotOMBiomass")

	n <- length( M )
	mdf <- NULL

	for(i in 1:n)
	{
		df <- data.frame(Model=as.factor(M[[i]]$seed),Biomass=M[[i]]$"ytB",Year=(M[[i]]$rep_yr+1):M[[i]]$proj_yr)
		mdf <- rbind(mdf,df)
	}

	p <- ggplot(mdf,aes(x=Year,y=Biomass, color=Model)) + geom_line()
	p <- p + labs(x="Year",y="Total Biomass")
	p <- p + ylim(0,max(df$Biomass))
	p <- p + theme_bw(11)
	print(p)



	cib <- data.frame(M[[1]]$"ytB")

	for(i in 2:n){

		bio <- data.frame(M[[i]]$"ytB")
		cib <- cbind(cib,bio)
	}

	qq<-apply(cib,1,calc_quantile)

	median<-qq["50%",]
	low<-qq["2.5%",]
	high<-qq["97.5%",]

	fdf<-data.frame (Median=median,Low=low, High=high,Year=(M[[i]]$rep_yr+1):M[[i]]$proj_yr)

	
	p <- ggplot(fdf,aes(x=Year,y=Median)) + geom_line()
	p <- p + geom_ribbon(aes(ymax=High, ymin=Low),alpha=0.2)
	p <- p + labs(x="Year",y="Total Biomass")
	p <- p + ylim(min(fdf$Low),max(fdf$High))
	p <- p + theme_bw(11)
	print(p)

	
}




