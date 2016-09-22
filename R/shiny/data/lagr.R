# ------------------------------------------------------#
# Lagragian movement model for shiny application		#
# Author: Catarina Wor									#
# Date: Nov 10th 2015 									#
# ------------------------------------------------------#


#--------------------------------
#Model dimensions
#--------------------------------

syr<- 1
nyr<-40 
years<-1:40

sage <- 1
nage <- 20

smon <- 1
nmon <- 12
sarea <- 30
narea <- 60
nations <- 2
border <- 48.7 

MD<-list(syr=syr,nyr=nyr, years=years,sage=sage,nage=nage,smon=smon,nmon=nmon,sarea=sarea,
	narea=narea,nations=nations,border=border)

#--------------------------------
#Model Parameters
#--------------------------------

Ro 	<- 2470000000
h  	<- 0.862 
m  	<- 0.213
fe 	<- 0
q 	<- 1
sigR  <-1
tau_c <-0.1
mo  <-3
err <-0

effm<- 1

wa <- c(0.13,0.41,0.47,0.48,0.54,0.57,0.62,0.66,0.72,0.70,1.16,1.02,0.95,0.97,1.06,1.06,1.06,1.06,1.06,1.06)
fa <- c(1.714009e-05,1.760926e-03,5.438130e-02,2.865322e-01,4.646053e-01,5.085980e-01,6.264716e-01,7.148981e-01,7.309999e-01,8.338000e-01,9.989000e-01,1.075200e+00,1.230300e+00,1.118700e+00,1.068200e+00,1.054500e+00,1.054500e+00,1.054500e+00,1.054500e+00,1.054000e+00)
va <- c(0.0041512,0.084884,0.60249,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
minPos <- c(30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0,30.0)
maxPos50 <- 4					
maxPossd <- 3			 
cvPos <-0.05 				
	
TotEffyear <- matrix(c(rep(1,length(years)),rep(1,length(years))),nrow=nations,byrow=T)
TotEffmonth<- matrix(c(0.0,0.0,0.0,0.0,0.5,1.0,1.0,1.0,0.5,0.1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,0.5,0.3,0.0,0.0),nrow=nations,byrow=T)

effPwr<-c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.5,0.5,0.5,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.5,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1)

MP<-list(syr=syr,nyr=nyr, years=years,sage=sage,nage=nage,smon=smon,nmon=nmon,sarea=sarea,
	narea=narea,nations=nations,border=border,
	Ro=Ro,h=h,m=m,fe=fe,q=q,sigR=sigR,tau_c=tau_c,mo=mo,err=err,wa=wa,fa=fa,va=va,minPos=minPos,
	maxPos50=maxPos50,maxPossd=maxPossd,cvPos=cvPos,TotEffyear=TotEffyear,TotEffmonth=TotEffmonth,effPwr=effPwr,effm=effm)


lagr_func<- function(MP){

	with(MP,{

	#--------------------------------
	# accessory quantities and counters
	#--------------------------------
		
	ntstp <-length(smon:nmon)*length(years)
	age  <- sage:nage
	areas<- sarea:narea
	nationareas<-c(sum(areas<border),sum(areas>=border))
	wt<-rep(1,length(years))
	
	
	indyr <- rep(years,each=length(years))
	indmonth <- rep(smon:nmon,length(years))
	indnatarea <- rep(1,length(areas))
	indnatarea[areas>=border] <- 2

	#--------------------------------
	# other derived quantities
	#--------------------------------
	lxo<-vector(length=length(age))
	za <-vector(length=length(age))
	
	lxo <- exp(-m*(age-1))
	lxo[nage] <- lxo[nage]/(1.0 - exp(-m)) 
	
	kappa <- 4*h/(1-h)
	phie <- sum(lxo*fa)
	So 	<- kappa/phie
	Bo 	<- kappa/So*Ro
	beta <- (kappa-1)/Bo
	
	m_tsp <- m/nmon
	za <- m+va*fe;
		
	
	#--------------------------------
	# Storage objects
	#--------------------------------
	
	
	SB <-vector(length=ntstp)	
	maxPos<-vector(length=length(age))
	varPos<-vector(length=length(age))
	
	NationVulB<-matrix(NA,nrow=ntstp,ncol=nations)
	Nage<-matrix(NA,nrow=ntstp,ncol=length(age))
	VulB<-matrix(NA,nrow=ntstp,ncol=length(age))
	PosX<-matrix(NA,nrow=ntstp,ncol=length(age))
	Effage<-matrix(NA,nrow=ntstp,ncol=length(age))
	VBarea<-matrix(NA,nrow=ntstp,ncol=length(areas))	
	 	
	Effarea<-matrix(NA,nrow=ntstp,ncol=length(areas))
	
	propVBareaAge<-array(NA,dim=c(ntstp,length(areas),length(age))) 
	NAreaAge<-array(NA,dim=c(ntstp,length(areas),length(age)))
	CatchAreaAge<-array(NA,dim=c(ntstp,length(areas),length(age)))
	CatchNatAge<-array(0,dim=c(ntstp,nations,length(age)))
	CatchNat<-matrix(NA,nrow=ntstp,ncol=nations)
	EffNatAge<-array(NA,dim=c(nations,ntstp,length(age)+2))
	
	#--------------------------------
	# initialization
	#--------------------------------

	TotEffyear<-TotEffyear*effm

	Nage[1,1] = So*Bo/(1+beta*Bo)

	for(i in (sage+1):nage){
		Nage[1,i] = Nage[1,i-1] * exp(-za[i-1]/12)
	}

	VulB[1,]<- Nage[1,]*va*wa
	SB[1] <- sum(Nage[1,]*fa*wa)/2

	maxPos <- 1./(1.+exp(-(age-maxPos50)/maxPossd))
	maxPos <- maxPos*(narea-sarea)
	maxPos <- maxPos+sarea	
	
	varPos <- maxPos*cvPos

	PosX[1,] <- minPos + (maxPos - minPos) * (0.5+0.5*sin(indmonth[1]*pi/6 - mo*pi/6)) 


	for(r in 1:length(areas)){
		
		VBarea[1,r]<-sum(VulB[1,]*(pnorm(areas[r]+0.5,PosX[1,],varPos)-pnorm(areas[r]-0.5,PosX[1,],varPos)))
		NAreaAge[1,r,] <- Nage[1,]*(pnorm(areas[r]+0.5,PosX[1,],varPos)-pnorm(areas[r]-0.5,PosX[1,],varPos))
	}


	NationVulB[1,1] <- sum(VBarea[1,1:(nationareas[1]-1)]) 
	NationVulB[1,2] <- sum(VBarea[1,(nationareas[1]):(narea-sarea)])
	

	tmp1<-vector(length=length(areas))
	tmp2<-vector(length=length(areas))
	
	for(r in 1:length(areas)){
		tmp1[r] <- (VBarea[1,r]/ (sum(VulB[1,]) + 0.0001))* effPwr[r]
		tmp2[r] <- tmp1[r]*TotEffyear[indnatarea[r],indyr[1]]
		Effarea[1,r] <- tmp2[r]*TotEffmonth[indnatarea[r],indmonth[1]]
	}
	
	
	for(a in 1:length(age)){
		
		for(r in 1:length(areas)){
			
			propVBareaAge[1,r,a] <-(pnorm(areas[r]+0.5,PosX[1,],varPos)-pnorm(areas[r]-0.5,PosX[1,],varPos))[a]
			CatchAreaAge[1,r,a] <- q*Effarea[1,r]/12*va[a]/(q*Effarea[1,r]/12*va[a]+m_tsp)*(1-exp(-(q*Effarea[1,r]/12*va[a]+m_tsp)))*NAreaAge[1,r,a];
			CatchNatAge[1,indnatarea[r],a] <- CatchNatAge[1,indnatarea[r],a]+CatchAreaAge[1,r,a];


			EffNatAge[indnatarea[r],1,1] = 1;
			EffNatAge[indnatarea[r],1,2] = indnatarea[r];
			EffNatAge[indnatarea[r],1,a+2] = EffNatAge[indnatarea[r],1,a+2]+Effarea[1,r]*propVBareaAge[1,r,a];

		}
		
		Effage[1,a] <- sum(Effarea[1,]* propVBareaAge[1,,a])

	}	


	#--------------------------------
	# move grow and die
	#--------------------------------
	
	for(i in 2:ntstp){
		
		if(indmonth[i]==1){

			Nage[i,1] = (So*SB[i-nmon]/(1.0+beta*SB[i-nmon]))

            for(a in 2:(length(age))){
            	Nage[i,a] = Nage[i-1,a-1]*exp(-(m_tsp+q*Effage[i-1,a-1]/12*va[a-1]))
            }

		}else{
			Nage[i,] = Nage[i-1,]*exp(-(m_tsp+q*Effage[i-1,]/12))
		}
       
		
		VulB[i,] = Nage[i,]*va*wa
		SB[i] = sum(Nage[i,]*fa*wa)/2;
		
		#maxPos <- 1./(1.+exp(-(age-maxPos50)/maxPossd))
		#maxPos <- maxPos*(narea-sarea)
		#maxPos <- maxPos+sarea
		
		varPos <- maxPos*cvPos;

		PosX[i,] <- minPos + (maxPos - minPos) * (0.5+0.5*sin(indmonth[i]*pi/6 - mo*pi/6)); 


		for(r in 1:(length(areas))){
			VBarea[i,r]<-sum(VulB[i,] *(pnorm(areas[r]+0.5,PosX[i,],varPos)-pnorm(areas[r]-0.5,PosX[i,],varPos)))
			NAreaAge[i,r,] <- Nage[i,]*(pnorm(areas[r]+0.5,PosX[i,],varPos)-pnorm(areas[r]-0.5,PosX[i,],varPos))
		}	

		NationVulB[i,1] <- sum(VBarea[i,1:(nationareas[1]-1)]) 
		NationVulB[i,2] <- sum(VBarea[i,(nationareas[1]):(narea-sarea)])
	
		tmp1<-vector(length=length(areas))
		tmp2<-vector(length=length(areas))
	
		for(r in 1:(length(areas))){
			tmp1[r] <- (VBarea[i,r]/ (sum(VulB[i,]) + 0.0001))* effPwr[r]
			tmp2[r] <- tmp1[r]*TotEffyear[indnatarea[r],indyr[i]]
			Effarea[i,r] <- tmp2[r]*TotEffmonth[indnatarea[r],indmonth[i]]
		}
		
		for(a in 1:(length(age))){			
		
			for(r in 1:(length(areas))){
				propVBareaAge[i,r,a] <-(pnorm(areas[r]+0.5,PosX[i,],varPos)-pnorm(areas[r]-0.5,PosX[i,],varPos))[a]
				CatchAreaAge[i,r,a] <- q*Effarea[i,r]/12*va[a]/(q*Effarea[i,r]/12*va[a]+m_tsp)*(1-exp(-(q*Effarea[i,r]/12*va[a]+m_tsp)))*NAreaAge[i,r,a];
				CatchNatAge[i,indnatarea[r],a] <- CatchNatAge[i,indnatarea[r],a]+CatchAreaAge[i,r,a];

				EffNatAge[indnatarea[r],i,1] = i;
				EffNatAge[indnatarea[r],i,2] = indnatarea[r];
				EffNatAge[indnatarea[r],i,a+2] = EffNatAge[indnatarea[r],i,a+2]+Effarea[i,r]*propVBareaAge[i,r,a];
			}
			
			Effage[i,a] <- sum(Effarea[i,]* propVBareaAge[i,,a])
		}

	}

	for(n in 1:nations){
		CatchNat[,n]<-apply(CatchNatAge[,n,],1,sum)
	}
	

	Mres<- list(ntsp=ntstp,nage=nage,sage=sage,ages=age, indmonth=indmonth,indyr=indyr,maxPos50=maxPos50,areas=areas,
	maxPossd=maxPossd,VBarea=VBarea,CatchNat=CatchNat)
	return(Mres)	
	})
}
#=========================================================================================



#======================================================================== 
# Graphs for high and low abundance scenarios case scenario
#======================================================================== 


#Mres<- list(ntsp=ntstp,nage=nage,sage=sage,ages=age, indmonth=indmonth,indyr=ind,maxPos50=maxPos50,
#	maxPossd=maxPossd,VBarea=VBarea)

library(ggplot2)
library(reshape2)

pl_vb<-function(Mres){

		#Mres<-lagr_func(MP)
		
		VBareas<-matrix(Mres[["VBarea"]], ncol=length(Mres[["areas"]]),dimnames=list(1:Mres[["ntsp"]],sarea:narea))
	
		VBplot<-melt(VBareas)
		
		names(VBplot)<- c("time","area", "VB")
		
		#lat<-VBplot$area
		#lon<-rep(-131,length(lat))
		#nation<-rep("U.S.A",nrow(VBplot))
		#nation[VBplot$lat>48.1]<-"Canada"
		month <- rep(Mres[["indmonth"]],length(Mres[["areas"]]))
		#indYr <- indyr[VBplot$time]	
		
		VBplot<-cbind(VBplot,month)

		meanVBnat1<-NULL
		meanVBnat2<-NULL
		
		for(i in 1:12){
		
		  meanVBnat1[i]<-sum(VBplot$VB[VBplot$area<48.7&VBplot$month==i&VBplot$time>468])
		  meanVBnat2[i]<-sum(VBplot$VB[VBplot$area>48.7&VBplot$month==i&VBplot$time>468])
		  
		}
		
		meanVBnat1p<-meanVBnat1/(meanVBnat1+meanVBnat2)
		meanVBnat2p<-meanVBnat2/(meanVBnat1+meanVBnat2)
		
		
		VulB<-c(meanVBnat1,meanVBnat2)
		month<-rep(1:12,2)
		nation<-c(rep("U.S.A",12),rep("Canada",12))
		
		#measure<-c(rep("nominal",48),rep("proportion",48))
		
		
		df<- data.frame(VulB,month,nation)

		return(df)
		
		print("inplotvb")
		#print(meanVBnat1p+meanVBnat2p)
		
		#p <- ggplot(df, aes(x=as.factor(month), y=VulB, fill=as.factor(nation)))
		#p <- p + geom_bar(stat = "identity")
		#p <- p + scale_y_continuous(limits=c(0, 1))
		#p <- p + labs(x="month", y="vulnerable biomass",fill="nation")
		#p <- p + theme_bw()
		#print(p)
	

}

	

plotvb<-function(MresA,MresB){

	dfA<-pl_vb(MresA)
	scenario<-rep("A",nrow(dfA))
	dfA<-cbind(dfA,scenario)
	

	dfB<-pl_vb(MresB)
	scenario<-rep("B",nrow(dfB))
	dfB<-cbind(dfB,scenario)

	df<-rbind(dfA,dfB)

	p <- ggplot(df, aes(x=as.factor(month), y=VulB, fill=as.factor(nation)))
	p <- p + geom_bar(stat = "identity")
	#p <- p + scale_y_continuous(limits=c(0, 1))
	p <- p + labs(x="month", y="vulnerable biomass",fill="nation")
	p <- p + facet_grid(.~scenario)
	p <- p + theme_bw()
	print(p)



}

 
pl_catch<-function(Mres){

		Mres<-lagr_func(MP)
		names(Mres)
		
		CatchNat<-matrix(Mres[["CatchNat"]], ncol=2,dimnames=list(1:Mres[["ntsp"]],c("U.S.A","Canada")))
		CatchNat<-CatchNat[469:480,]	

		TotCatchNat<-apply(CatchNat,2,sum)


		VBareas<-matrix(Mres[["VBarea"]], ncol=length(Mres[["areas"]]),dimnames=list(1:Mres[["ntsp"]],sarea:narea))
	
		VBplot<-melt(VBareas)
		
		names(VBplot)<- c("time","area", "VB")
			
		
		sumVBnat1<-sum(VBplot$VB[VBplot$area<48.7&VBplot$time>468])
		sumVBnat2<-sum(VBplot$VB[VBplot$area>48.7&VBplot$time>468])

		value<-c(.7388*sum(TotCatchNat),TotCatchNat[1],.2612*sum(TotCatchNat),TotCatchNat[2])
		nation<-c(rep("U.S.A",2),rep("Canada",2))
		variable<-rep(c("TAC","Possible Catch"),2)
		
		print("inplotcatch")

		df<-data.frame(value,nation,variable)

		

}


plotcatch<-function(MresA,MresB){

	dfA<-pl_catch(MresA)
	scenario<-rep("A",nrow(dfA))
	dfA<-cbind(dfA,scenario)
	

	dfB<-pl_catch(MresB)
	scenario<-rep("B",nrow(dfB))
	dfB<-cbind(dfB,scenario)

	df<-rbind(dfA,dfB)
	
	p <- ggplot(df, aes(x=as.factor(variable), y=value, fill=as.factor(nation)))
	p <- p + geom_bar(stat = "identity")
	p <- p + labs(x="", y="vulnerable biomass",fill="nation")
	p <- p + facet_grid(scenario~nation,scales = "free_y")
	p <- p + theme_bw()
	print(p)	

}	





