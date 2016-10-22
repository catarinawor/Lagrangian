#==============================================
#Title:plot_params
#Author: Catarina Wor
#date: Oct 21th 2016
#Function to plot true biomass fro lagrangian model.extract
# code adapted from iSCAM stuff (Martell et al)
#==============================================

M<-SAdat

require(reshape2)
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


	mdf<-melt(mdf, id=c("Year", "Model"),value.name="parameter")


	p <- ggplot(mdf,aes(x=Year,y=parameter)) + geom_line()
	p <- p + facet_wrap(~variable,scales="free")
	p <- p + labs(x="Year",y="steepness estimate")
	p <- p + theme_bw(11)
	p

	
}