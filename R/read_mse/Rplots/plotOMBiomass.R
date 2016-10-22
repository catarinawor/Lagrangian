#==============================================
#Title:plotOMBiomass
#Author: Catarina Wor
#date: Oct 21th 2016
#Function to plot true biomass fro lagrangian model.extract
# code adapted from iSCAM stuff (Martell et al)
#==============================================

#for testing - erase when done

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

	
}




