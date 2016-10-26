#==============================================
#Title:plotOMcatch
#Author: Catarina Wor
#date: Oct 21th 2016
#Function to plot realized catches from OM runs
# code adapted from iSCAM stuff (Martell et al)
#==============================================



source("/Users/catarinawor/Documents/Lagrangian/R/Rplots/calc_quantile.R")
require(ggplot2)
plotOMCatch <- function( M )
{
	cat("plotOMCatch")

	n <- length( M )
	mdf <- NULL

	for(i in 1:n)
	{
		df <- data.frame(Model=as.factor(M[[i]]$seed),Nation=M[[i]]$"yCatchStateAge"[,2],Catch=apply(M[[i]]$"yCatchStateAge"[,-c(1,2)],1,sum),Year=M[[i]]$"yCatchStateAge"[,1])
		mdf <- rbind(mdf,df)
	}

	mdf<- mdf[mdf$Year>M[[1]]$rep_yr,]

	p <- ggplot(mdf,aes(x=Year,y=Catch, color=Model)) + geom_line()
	p <- p + labs(x="Year",y="Total Biomass")
	p <- p + ylim(0,max(mdf$Catch))
	p <- p + facet_wrap(~Nation,scales="free")
	p <- p + theme_bw(11)
	print(p)

}

plotOMCICatch <- function( M )
{
	cat("plotOMCatch")

	cib <- data.frame(apply(M[[1]]$"yCatchStateAge"[,-c(1,2)],1,sum))

	for(i in 2:n){

		ci <- data.frame(apply(M[[i]]$"yCatchStateAge"[,-c(1,2)],1,sum))
		cib <- cbind(cib,ci)
	}

	qq<-apply(cib,1,calc_quantile)

	median<-qq["50%",]
	low<-qq["2.5%",]
	high<-qq["97.5%",]

	fdf<-data.frame(Median=median,Low=low, High=high,Year=M[[1]]$"yCatchStateAge"[,1],Nation=as.factor(M[[1]]$"yCatchStateAge"[,2]))

	fdf<- fdf[fdf$Year>M[[1]]$rep_yr,]
	
	p <- ggplot(fdf,aes(x=Year,y=Median)) + geom_line(aes(color=Nation))
	p <- p + geom_ribbon(aes(ymax=High, ymin=Low,fill=Nation), alpha=0.3)
	p <- p + labs(x="Year",y="Catch")
	p <- p + ylim(min(fdf$Low),max(fdf$High))
	p <- p + theme_bw(11)
	print(p)

	
}
