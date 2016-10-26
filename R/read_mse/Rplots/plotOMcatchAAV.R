#==============================================
#Title:plotOMcatch
#Author: Catarina Wor
#date: Oct 21th 2016
#Function to plot realized catches from OM runs
# code adapted from iSCAM stuff (Martell et al)
#==============================================



source("/Users/catarinawor/Documents/Lagrangian/R/Rplots/calc_quantile.R")
require(ggplot2)
require(reshape2)
library(tidyr)



plotOMCatchAAV <- function( M )
{
	cat("plotOMCatchAAV")

	n <- length( M )
	mdf <- NULL

	for(i in 1:n)
	{

		df <- data.frame(Model=as.factor(M[[i]]$seed),Nation=M[[i]]$"yCatchStateAge"[,2],Catch=apply(M[[i]]$"yCatchStateAge"[,-c(1,2)],1,sum),Year=M[[i]]$"yCatchStateAge"[,1])
		mdf <- rbind(mdf,df)
	}

	mdf<- mdf[mdf$Year>M[[1]]$rep_yr,]

	
	
	yrs<-unique(mdf$Year)

	mdf$AAV[mdf$Year==yrs[1]]<-0
	for(y in 2:length(yrs)){

		mdf$AAV[mdf$Year==yrs[y]]<-abs(mdf$Catch[mdf$Year==yrs[y]]-mdf$Catch[mdf$Year==yrs[y-1]])/(mdf$Catch[mdf$Year==yrs[y]]+mdf$Catch[mdf$Year==yrs[y-1]])
	}

		


	p <- ggplot(mdf,aes(x=Year,y=AAV, color=Model)) + geom_line()
	p <- p + labs(x="Year",y="Total Biomass")
	#p <- p + ylim(0,max(mdf$Catch))
	p <- p + facet_wrap(~Nation,scales="free")
	p <- p + theme_bw(11)
	print(p)

}




plotOMCICatchAAV <- function( M )
{
	cat("plotOMCICatchAAV")


	n <- length( M )
	mdf <- NULL

	for(i in 1:n)
	{

		df <- data.frame(Model=as.factor(M[[i]]$seed),Nation=M[[i]]$"yCatchStateAge"[,2],Catch=apply(M[[i]]$"yCatchStateAge"[,-c(1,2)],1,sum),Year=M[[i]]$"yCatchStateAge"[,1])
		mdf <- rbind(mdf,df)
	}

	mdf<- mdf[mdf$Year>M[[1]]$rep_yr,]

	
	
	yrs<-unique(mdf$Year)

	mdf$AAV[mdf$Year==yrs[1]]<-0
	for(y in 2:length(yrs)){

		mdf$AAV[mdf$Year==yrs[y]]<-(mdf$Catch[mdf$Year==yrs[y]]-mdf$Catch[mdf$Year==yrs[y-1]])/(mdf$Catch[mdf$Year==yrs[y]]+mdf$Catch[mdf$Year==yrs[y-1]])
	}

	cdf<-data.frame(Model=mdf$Model, Year=mdf$Year, AAV=mdf$AAV, Nation=mdf$Nation)

	cib<-spread(cdf, Model, AAV)

	

	qq<-apply(cib[,-c(1,2)],1,calc_quantile)

	median<-qq["50%",]
	low<-qq["2.5%",]
	high<-qq["97.5%",]

	fdf<-data.frame(Median=median,Low=low, High=high,Year=cib$Year,Nation=as.factor(cib$Nation))

	
	p <- ggplot(fdf,aes(x=Year,y=Median)) + geom_line(aes(color=Nation))
	p <- p + geom_ribbon(aes(ymax=High, ymin=Low,fill=Nation), alpha=0.3)
	p <- p + labs(x="Year",y="Catch")
	p <- p + ylim(min(fdf$Low),max(fdf$High))
	p <- p + facet_wrap(~Nation,scales="free")
	p <- p + theme_bw(11)
	print(p)

	
}
