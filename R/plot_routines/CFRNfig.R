# ------------------------------------------------------#
# Lagragian movement model for shiny application		#
# Author: Catarina Wor									#
# Date: Nov 15th 2015 									#
# Figures for 5 min presentation						#
# ------------------------------------------------------#

source("/Users/catarinawor/Documents/Lagrangian/shiny/data/lagr.R")

MP$effm<-6
MP$maxPos50<-8
highFcold<-lagr_func(MP)

#===================================================
MP$effm<-1
MP$maxPos50<-8
lowFcold<-lagr_func(MP)

#===================================================
MP$effm<-6
MP$maxPos50<-5
highFwarm<-lagr_func(MP)

#===================================================
MP$effm<-1
MP$maxPos50<-5
lowFwarm<-lagr_func(MP)


VBareasHC<-matrix(highFcold[["VBarea"]], ncol=length(highFcold[["areas"]]),dimnames=list(1:highFcold[["ntsp"]],sarea:narea))
VBareasLC<-matrix(lowFcold[["VBarea"]], ncol=length(lowFcold[["areas"]]),dimnames=list(1:lowFcold[["ntsp"]],sarea:narea))
VBareasHW<-matrix(highFwarm[["VBarea"]], ncol=length(highFwarm[["areas"]]),dimnames=list(1:highFwarm[["ntsp"]],sarea:narea))
VBareasLW<-matrix(lowFwarm[["VBarea"]], ncol=length(lowFwarm[["areas"]]),dimnames=list(1:lowFwarm[["ntsp"]],sarea:narea))

VBplotHC<-melt(VBareasHC)
VBplotLC<-melt(VBareasLC)
VBplotHW<-melt(VBareasHW)
VBplotLW<-melt(VBareasLW)

names(VBplotHC)<- c("time","area", "VB")
names(VBplotLC)<- c("time","area", "VB")
names(VBplotHW)<- c("time","area", "VB")
names(VBplotLW)<- c("time","area", "VB")	


VBplotHC<-VBplotHC[VBplotHC$time>468,]
VBplotLC<-VBplotLC[VBplotLC$time>468,]
VBplotHW<-VBplotHW[VBplotHW$time>468,]
VBplotLW<-VBplotLW[VBplotLW$time>468,]


regime<-rep("Cold",nrow(VBplotHC))
fishing<-rep("High",nrow(VBplotHC))
month<-rep(1:12,length(unique(VBplotHC$area)))
VBplotHC<-cbind(VBplotHC,regime,fishing,month)	

regime<-rep("Cold",nrow(VBplotLC))
fishing<-rep("Low",nrow(VBplotLC))
month<-rep(1:12,length(unique(VBplotLC$area)))
VBplotLC<-cbind(VBplotLC,regime,fishing,month)	

regime<-rep("Warm",nrow(VBplotHW))
fishing<-rep("High",nrow(VBplotHW))
month<-rep(1:12,length(unique(VBplotHW$area)))
VBplotHW<-cbind(VBplotHW,regime,fishing,month)	

regime<-rep("Warm",nrow(VBplotLW))
fishing<-rep("Low",nrow(VBplotLW))
month<-rep(1:12,length(unique(VBplotLW$area)))
VBplotLW<-cbind(VBplotLW,regime,fishing,month)	

		

		sumVBnat1HC<-NULL
		sumVBnat2HC<-NULL

		sumVBnat1LC<-NULL
		sumVBnat2LC<-NULL

		sumVBnat1HW<-NULL
		sumVBnat2HW<-NULL

		sumVBnat1LW<-NULL
		sumVBnat2LW<-NULL
		

		for(i in 1:12){
		
		  sumVBnat1HC[i]<-sum(VBplotHC$VB[VBplotHC$area<48.7&VBplotHC$month==i])
		  sumVBnat2HC[i]<-sum(VBplotHC$VB[VBplotHC$area>48.7&VBplotHC$month==i])

		  sumVBnat1LC[i]<-sum(VBplotLC$VB[VBplotLC$area<48.7&VBplotLC$month==i])
		  sumVBnat2LC[i]<-sum(VBplotLC$VB[VBplotLC$area>48.7&VBplotLC$month==i])
		  
		  sumVBnat1HW[i]<-sum(VBplotHW$VB[VBplotHW$area<48.7&VBplotHW$month==i])
		  sumVBnat2HW[i]<-sum(VBplotHW$VB[VBplotHW$area>48.7&VBplotHW$month==i])

		  sumVBnat1LW[i]<-sum(VBplotLW$VB[VBplotLW$area<48.7&VBplotLW$month==i])
		  sumVBnat2LW[i]<-sum(VBplotLW$VB[VBplotLW$area>48.7&VBplotLW$month==i])
		}
		
		#meanVBnat1p<-meanVBnat1/(meanVBnat1+meanVBnat2)
		#meanVBnat2p<-meanVBnat2/(meanVBnat1+meanVBnat2)
		
		
		VulB<-c(sumVBnat1HC,sumVBnat2HC,sumVBnat1LC,sumVBnat2LC,sumVBnat1HW,sumVBnat2HW,sumVBnat1LW,sumVBnat2LW)
		month<-rep(1:12,8)
		nation<-rep(c(rep("U.S.A",12),rep("Canada",12)),4)
		regime<-c(rep(rep("Cold",24),2),rep(rep("Warm",24),2))
		fishing<-rep(c(rep("High F",24),rep("Low F",24)),2)
			
		
		#measure<-c(rep("nominal",48),rep("proportion",48))
		
		
		df<- data.frame(VulB,month,nation,regime, fishing)
		
		setwd("/Users/catarinawor/Documents/hake/CFRN/AGM_5meeting/presentation")

		p <- ggplot(df, aes(x=as.factor(month), y=VulB, fill=as.factor(nation)))
		p <- p + geom_bar(stat = "identity")
		#p <- p + scale_y_continuous(limits=c(0, 1))
		p <- p + facet_grid(regime~fishing)
		p <- p + labs(x="month", y="vulnerable biomass",fill="nation")
		p <- p + theme_bw()
		print(p)

		ggsave(filename ="availEffect.png")


