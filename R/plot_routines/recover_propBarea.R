#================================================================
#routine ro recover propBarea from 5 area gtg scenarios

#================================================================
setwd("/Users/catarinawor/Documents/Lagrangian/")
source("R/read.admb.R")

.SIMDIRS   <- c(#"/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_3areas_tau04",
                #"/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_3areas_tau04_delta2",
                #"/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_3areas_tau1"#,
                #"/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_3areas_tau1_delta2",
                #"/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_5areas_tau04",
                #"/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_5areas_tau04_delta2"#,
                "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_5areas_tau1"#,
               #"/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_5areas_tau1_delta2"              
  )

ALLDF<-list(length(.SIMDIRS))

.SIMNAME<-list(length(.SIMDIRS))

for(i in 1:length(.SIMDIRS)){
.SIMNAME[[i]]   <- list.files(.SIMDIRS[i],pattern="\\.Rdata", full.name=TRUE)
  
  dfall<-NULL

  Fs<-NULL
  
  ctr<-0
  for( j in 1:length(.SIMNAME[[i]])){

		#read in parameter estimate
 			load(.SIMNAME[[i]][j])

 			names(sims[[3]])
 			sims[[3]]$maxgrad

 		if(sims[[3]]$maxgrad<0.01){
 			ctr<-ctr+1

 		setwd("/Users/catarinawor/Documents/RecLag/admb/OM/gtg")
 			sink("seed.txt")
 			cat(sims[[1]]$seed,"\n")
 			sink()

		#write a lagrangian_OM_gtg.dat

 		#setwd("/Users/catarinawor/Documents/RecLag/admb/OM/gtg")


 		sink("lagrangian_OM_gtg.dat")
 		cat("#syr\n")
		cat("1\n")
		cat("#nyr\n")
		cat("100\n")
		cat("#rep_yr\n")
		cat("70\n")
		cat("#proj_yr\n")
		cat("#100\n")
		cat("#sage\n")
		cat("1\n")
		cat("#nage\n")
		cat("20\n")
		cat("#smon\n")
		cat("1\n")
		cat("#nmon\n")
		cat("12\n")
		cat("#sarea\n")
		cat("30\n")
		cat("#narea\n")
		cat("60\n")
		cat("#ngroup\n")
		cat("20\n")
		cat("#nations\n")
		cat("2\n")
		cat("#border\n")
		cat("48.5\n")
		cat("#fisharea\n")
		cat("5\n")
		cat("#fishbound\n")
		cat("42 46 48.5 51 \n")
		cat("#Ro\n")
		cat("2.923\n")
		cat("#h\n")
		cat("0.814 \n")
		cat("#m\n")
		cat("0.223\n")
		cat("#fe\n")
		cat("0\n")
		cat("#q \n")
		cat("1\n")
		cat("#fbeta\n")
		cat("1.0\n")
		cat("#sigR\n")
		cat("0.6\n")
		cat("#tau_c\n")
		cat("1.0\n")
		cat("#err\n")
		cat("1.0\n")
		cat("#linf - cte\n")
		cat("53.2 53.2 53.2 53.2 53.2 53.2 53.2 53.2 53.2 53.2\n")
		cat("53.2 53.2 53.2 53.2 53.2 53.2 53.2 53.2 53.2 53.2\n")
		cat("#vbk\n")
		cat("0.3\n")
		cat("#to\n")
		cat("-0.5\n")
		cat("#mo\n")
		cat(exp(sims[[3]]$est[1]),"\n")
		cat("#maxPos50\n")
		cat(exp(sims[[3]]$est[3]),"\n")
		cat("#maxPossd\n")
		cat(exp(sims[[3]]$est[4]),"\n")
		cat("#cvpos\n")
		cat(exp(sims[[3]]$est[2]),"\n")	
		cat("#weight at age \n")
		cat("0.13 0.41 0.47 0.48 0.54 0.57 0.62 0.66 0.72 0.70 1.16 1.02 0.95 0.97 1.06 1.06 1.06 1.06 1.06 1.06 \n")
		cat("#fecundity at age\n")
		cat("1.393163e-05 1.431298e-03 4.420166e-02 2.328962e-01 3.776358e-01 4.133935e-01 5.092023e-01 5.810762e-01 5.941639e-01 6.777209e-01 8.119158e-01 8.739332e-01 1.000000e+00 9.092904e-01 8.682435e-01 8.571080e-01 8.571080e-01 8.571080e-01 8.571080e-01 8.567016e-01\n")
		cat("#vulnerability at age\n")
		cat("0.5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1\n")
		cat("#mininmum position by age\n")
		cat("32.0 32.0 32.0 32.0 32.0 32.0 32.0 32.0 32.0 32.0 32.0 32.0 32.0 32.0 32.0 32.0 32.0 32.0 32.0 32.0\n")
		cat("#Fmultiplier\n")
		cat(exp(sims[[3]]$est[5]),"\n")	
		cat("# Total effort by country and year\n")
		cat("#area1\n")
		cat("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n")
		cat("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n")
		cat("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n")
		cat("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n")
		cat("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n")
		cat("#area2\n")
		cat("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n")
		cat("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n")
		cat("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n")
		cat("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n")
		cat("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n")
		cat("#area3\n")
		cat("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n")
		cat("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n")
		cat("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n")
		cat("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n")
		cat("1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 \n")
		cat("#area4 \n")
		cat("0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 \n")
		cat("0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 \n")
		cat("0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 \n")
		cat("0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 \n")
		cat("0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2\n")
		cat("#area5\n")
		cat("0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 \n")
		cat("0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 \n")
		cat("0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 \n")
		cat("0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 \n")
		cat("0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2\n")
		cat("#Total effort by country and month\n")
		cat("0.0 0.0 0.0 0.0 0.5 1.0 1.0 1.0 0.5 0.1 0.0 0.0\n")  
		cat("0.0 0.0 0.0 0.0 0.5 1.0 1.0 1.0 0.5 0.1 0.0 0.0\n")
		cat("0.0 0.0 0.0 0.0 0.5 1.0 1.0 1.0 0.5 0.1 0.0 0.0\n")  
		cat("0.0 0.0 0.0 0.0 0.0 1.0 1.0 1.0 0.5 0.3 0.0 0.0\n")  
		cat("0.0 0.0 0.0 0.0 0.0 1.0 1.0 1.0 0.5 0.3 0.0 0.0\n") 
		cat("#Effort power by latitude deg\n")
		cat("0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.1 0.5 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 0.5 0.1 0.1 \n")
		cat("0.1 0.1 0.1 0.1 0.1 0.0 \n")
		cat("#satype\n")
		cat("2\n")
		cat("#nationTACprop\n")
		cat("73.88 26.12\n")
		cat("#eof\n")
		cat("999\n")
		sink()

		#run lagrangian_OM_gtg.dat
		system("./lagrangian_OM_gtg")


		#read in lagrangian_OM_gtg.dat.rep
		simrec <- read.rep("lagrangian_OM_gtg.rep")
		file.name <- paste("simrec",sims[[1]]$seed,".Rdata",sep="")
		#write recovered .Rdata
		setwd("/Users/catarinawor/Documents/Lagrangian/simeval/recovered_gtg_5areas_tau1_delta2")
		save(simrec,file=file.name)


		test<-aggregate(sims[[1]]$propVBarea[,4:23],by=list(sims[[1]]$propVBarea[,1],sims[[1]]$propVBarea[,3]),sum)

   #dim(test)
		
    #length(apply(test[,3:22],1,sum))

    dfsi<-data.frame(indmonth=sims[[1]]$indmonth[test[,1]],indyr=sims[[1]]$indyr[test[,1]],
      area=test[,2], VB=apply(test[,3:22],1,sum))
    summary(dfsi)

    #dfsi<-cbind(sims[[1]]$indmonth,sims[[1]]$indyr,melt(sims[[1]]$VBarea, value.name="vb"))
    #names(dfsi)<-c("indmonth", "indyr","tstp","area","vb")
    dfsi1<-dfsi[dfsi$indyr>71,]
   
    dfsi2<-aggregate(dfsi1$VB,by=list(dfsi1$area,dfsi1$indmonth), median)
    names(dfsi2)<-c("area", "month","medianvb")
    dfsi2$type<-"simulated"
    #summary(dfsi2)


    #names(sims[[2]])
    tempest<-aggregate(simrec$propVBarea[,4:23],by=list(simrec$propVBarea[,1],simrec$propVBarea[,3]),sum)
    # summary(tempest)
    
    dfe<-data.frame(indmonth=simrec$indmonth[tempest[,1]],indyr=simrec$indyr[tempest[,1]],
      area=tempest[,2], VB=apply(tempest[,3:22],1,sum))
    #summary(dfe)

    dfe1<-dfe[dfe$indyr>71,]

    dfe2<-aggregate(dfe1$VB,by=list(dfe1$area,dfe1$indmonth), median)
    names(dfe2)<-c("area", "month","medianvb")
    dfe2$type<-"estimated"
    #summary(dfe2)


    #cbind(indmonth,indyr,melt(sims[[2]]$VBarea, value.name="vb"))
    #names(dfes)<-c("indmonth", "indyr","tstp","area","vb")
    #dfes1<-dfes[dfes$indyr>71,]
    #dfes1$area<-as.factor(as.numeric(dfes1$area)+29)
   
    #dfes2<-aggregate(dfes1$vb,by=list(dfes1$area,dfes1$indmonth), median)
    #names(dfes2)<-c("area", "month","medianvb")
    #dfes2$type<-"estimated"

    dfall<-rbind(dfall,dfsi2,dfe2)



	}
	}
	ALLDF[[i]]<-dfall

}


names(ALLDF[[1]])
df<-ALLDF[[1]]
summary(df)

meses<-c("Jan","Feb","Mar", "Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
df$meses<-meses[df$month]
df$meses <- factor(df$meses, levels = meses)

df$type<-as.factor(dfall$type)

p <- ggplot(df, aes(x=as.factor(area), y=medianvb, color=type))
p <- p + geom_boxplot(outlier.shape = NA)
p <- p + geom_vline(xintercept = 48.5)
p <- p + facet_wrap(~meses,ncol=3)
p <- p + theme_bw(14)
p <- p + scale_colour_grey(start = 0.1, end = 0.6)
p <- p + labs(x="Latitude",y="Biomass")
p <- p + theme(legend.position="bottom",legend.title=element_blank(),
              axis.text.x = element_text( angle=45))
p <- p + scale_x_discrete(breaks=seq(30,60,by=3))
#p <- p + scale_x_discrete(limits=seq(30,60,by=2))
p


setwd("/Users/catarinawor/Documents/hake/Lag_Model_paper/")
 ggsave("movement_gtg_estimate5at1.pdf", plot=p,width =18, height = 24,units="cm")
  


for(i in 1:length(.SIMDIRS)){
.SIMNAME[[i]]   <- list.files(.SIMDIRS[i],pattern="\\.Rdata", full.name=TRUE)
  
  dfall<-NULL

  Fs<-NULL
  
  mean(Fs+(.1*.1/2))

  for( j in 1:length(.SIMNAME[[i]])){
  	load(.SIMNAME[[i]][j])
  	Fs[j]<-exp(sims[[3]]$est[5])

  }
}

