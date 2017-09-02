#==============================================
#Title: Read in simulation outcomes and plot the simulation results
#parameter estimates and biomass
#Author: Catarina Wor
#date: Oct 6 2015
#
#==============================================


library(ggplot2)
library(reshape2)




#=============================================================
#simple
.SIMDIRS   <- c("/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_3areas_tau04",
                "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_3areas_tau04_delta2",
                "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_3areas_tau1",
                "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_3areas_tau1_delta2",
                "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_5areas_tau04",
                "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_5areas_tau04_delta2",
                "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_5areas_tau1",
                "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_5areas_tau1_delta2"#,
 )



.SIMNAME<-list(length(.SIMDIRS))

estn<-list(length(.SIMDIRS))
pbias<-list(length(.SIMDIRS))
maxgrad<-list(length(.SIMDIRS))
initvals<-list(length(.SIMDIRS))
initvals_bad<-list(length(.SIMDIRS))



for( i in 1:length(.SIMDIRS)){
  .SIMNAME[[i]]   <- list.files(.SIMDIRS[i],pattern="\\.Rdata", full.name=TRUE)
  
  tmp_estn<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=5)
  tmp_pbias<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=5)
  tmp_maxgrad<-vector(length=length(.SIMNAME[[i]]))
  tmp_initvals<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=5)
  

  for( j in 1:length(.SIMNAME[[i]])){
    load(.SIMNAME[[i]][j])

    true_pars <- c(sims[[1]]$"mo",sims[[1]]$"cvPos",sims[[1]]$"maxPos50",sims[[1]]$"maxPossd",mean(sims[[1]]$Fmult))  


    #parameters
    tmp_estn[j,]<-exp(sims[[3]]$est[1:5])
    tmp_pbias[j,]<-((tmp_estn[j,]-true_pars)/true_pars)*100
    tmp_maxgrad[j]<-sims[[3]]$maxgrad
    tmp_initvals[j,]<-exp(unlist(sims[[5]][1:5]))
   }

  tmp_estn<- tmp_estn[tmp_maxgrad<=.1,]
  tmp_pbias<- tmp_pbias[tmp_maxgrad<=.1,]
  
  estn[[i]]<-tmp_estn
  pbias[[i]]<-tmp_pbias
  maxgrad[[i]]<-tmp_maxgrad
  initvals[[i]]<-tmp_initvals[tmp_maxgrad<=.1,]
  initvals_bad[[i]]<-tmp_initvals[tmp_maxgrad>.1,]


}


titulos<-c( 

            "3 areas, tau=0.4, B = 1.0",
            "3 areas, tau=0.4, B = 2.0",
            "3 areas, tau=1.0, B = 1.0",
            "3 areas, tau=1.0, B = 2.0",
            "5 areas, tau=0.4, B = 1.0",  
           "5 areas, tau=0.4, B = 2.0",
            "5 areas, tau=1.0, B = 1.0",
           "5 areas, tau=1.0, B = 2.0"
  )

#setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")
#setwd("/Users/catarinawor/Documents/hake/JTC_talk")
#pdf("single_version_simeval.pdf", width=14, height=7)
setwd("/Users/catarinawor/Documents/hake/Lag_Model_paper/")
pdf("Figure3.pdf", width=14, height=7)
par(mfcol=c(2,4))
for( i in 1:length(.SIMDIRS)){
  boxplot(pbias[[i]],names=c(expression("t"[0]),expression("CV"),expression("a"[50]),
    expression(sigma["X"["max"]]),expression("q")),ylim=c(-30,30),main=titulos[i],cex.axis=1.5,
    cex.lab=2,cex.main=2,cex=1.6)
  abline(h=0)
  #text(4.8, y = -65, labels = nrow(pbias[[i]]), cex=2)
  text(4.8, y = 25, labels = paste ("scenario",i), cex=1.2)
}
mtext("% relative error", 2, line = -2, outer = TRUE, font=2)
mtext("Parameter", 1, line = -2, outer = TRUE, font=2)
dev.off()



#=============================================================

#gtg

#parameter estimates
.SIMDIRS   <- c("/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_3areas_tau04",
                "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_3areas_tau04_delta2",
                "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_3areas_tau1",
                "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_3areas_tau1_delta2",
                "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_5areas_tau04",
                "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_5areas_tau04_delta2",
                 "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_5areas_tau1",
               "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_5areas_tau1_delta2"              
              
  )



.SIMNAME<-list(length(.SIMDIRS))

estn<-list(length(.SIMDIRS))
pbias<-list(length(.SIMDIRS))
maxgrad<-list(length(.SIMDIRS))
initvals<-list(length(.SIMDIRS))
initvals_bad<-list(length(.SIMDIRS))



for( i in 1:length(.SIMDIRS)){
  .SIMNAME[[i]]   <- list.files(.SIMDIRS[i],pattern="\\.Rdata", full.name=TRUE)
  
  tmp_estn<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=5)
  tmp_pbias<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=5)
  tmp_maxgrad<-vector(length=length(.SIMNAME[[i]]))
  tmp_initvals<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=5)
  

  for( j in 1:length(.SIMNAME[[i]])){
    load(.SIMNAME[[i]][j])

    true_pars <- c(sims[[1]]$"mo",sims[[1]]$"cvPos",sims[[1]]$"maxPos50",sims[[1]]$"maxPossd",mean(sims[[1]]$Fmult))  


    #parameters
    tmp_estn[j,]<-exp(sims[[3]]$est[1:5])
    tmp_pbias[j,]<-((tmp_estn[j,]-true_pars)/true_pars)*100
    tmp_maxgrad[j]<-sims[[3]]$maxgrad
    tmp_initvals[j,]<-exp(unlist(sims[[5]][1:5]))
   }

  tmp_estn<- tmp_estn[tmp_maxgrad<=.01,]
  tmp_pbias<- tmp_pbias[tmp_maxgrad<=.01,]
  
  estn[[i]]<-tmp_estn
  pbias[[i]]<-tmp_pbias
  maxgrad[[i]]<-tmp_maxgrad
  initvals[[i]]<-tmp_initvals[tmp_maxgrad<=.01,]
  initvals_bad[[i]]<-tmp_initvals[tmp_maxgrad>.01,]


}


titulos<-c( 
            "3 areas, tau=0.4, B = 1.0",
            "3 areas, tau=0.4, B = 2.0",
            "3 areas, tau=1.0, B = 1.0",
            "3 areas, tau=1.0, B = 2.0",
            "5 areas, tau=0.4, B = 1.0",  
           "5 areas, tau=0.4, B = 2.0",
            "5 areas, tau=1.0, B = 1.0",
           "5 areas, tau=1.0, B = 2.0"
  )

#setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")
#setwd("/Users/catarinawor/Documents/hake/JTC_talk")
setwd("/Users/catarinawor/Documents/hake/Lag_Model_paper/")
#pdf("Figure4.pdf", width=14, height=7)
par(mfcol=c(2,4))
for( i in 1:length(.SIMDIRS)){
  boxplot(pbias[[i]],names=c(expression("t"[0]),expression("CV"),expression("a"[50]),
    expression(sigma["X"["max"]]),expression("q")),ylim=c(-40,40),main=titulos[i],cex.axis=1.5,
    cex.lab=2,cex.main=2,cex=1.6)
  abline(h=0)
  text(4.8, y = -25, labels = nrow(pbias[[i]]), cex=2)
  text(4.8, y = 35, labels = paste ("scenario",i), cex=1.2)
}
mtext("% relative error", 2, line = -2, outer = TRUE, font=2)
mtext("Parameter", 1, line = -2, outer = TRUE, font=2)
#dev.off()


#===============================================================
#===============================================================
#===============================================================
#come up with an indicator of the impact of biased parameter estimates on movement

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

i=1
  .SIMNAME[[i]]   <- list.files(.SIMDIRS[i],pattern="\\.Rdata", full.name=TRUE)
  
  dfall<-NULL
  

  for( j in 1:length(.SIMNAME[[i]])){

    #j=1
    #i=1

    load(.SIMNAME[[i]][j])
    
    #names(sims[[1]])
    #names(sims[[2]])
    # dim(sims[[1]]$propVBarea)
    # sims[[2]]$propVBarea

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
    tempest<-aggregate(sims[[2]]$propVBarea[,4:23],by=list(sims[[2]]$propVBarea[,1],sims[[2]]$propVBarea[,3]),sum)
    # summary(tempest)
    
    dfe<-data.frame(indmonth=sims[[2]]$indmonth[tempest[,1]],indyr=sims[[2]]$indyr[tempest[,1]],
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
#ALLDF[[i]]<-dfall
#}


summary(dfall)

meses<-c("Jan","Feb","Mar", "Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
dfall$meses<-meses[dfall$month]
dfall$meses <- factor(dfall$meses, levels = meses)

dfall$type<-as.factor(dfall$type)

p <- ggplot(dfall, aes(x=as.factor(area), y=medianvb, color=type))
p <- p + geom_boxplot(outlier.shape = NA)
p <- p + geom_vline(xintercept = 48.5)
p <- p + facet_wrap(~meses,ncol=3)
p <- p + theme_bw(14)
p <- p + scale_colour_grey(start = 0.1, end = 0.6,labels = c("simulated", "estimated"))
p <- p + labs(x="Latitude",y="Biomass")
p <- p + theme(legend.position="bottom",legend.title=element_blank(),
              axis.text.x = element_text( angle=45))
p <- p + scale_x_discrete(breaks=seq(30,60,by=2))
#p <- p + scale_x_discrete(limits=seq(30,60,by=2))
p

setwd("/Users/catarinawor/Documents/hake/Lag_Model_paper/")
      ggsave("movement_gtg_estimate.pdf", plot=p)
  


axis.text.x = element_text(face="bold", 
                           size=14, angle=45),)


dfe<-cbind(indmonth,indyr,melt(sim$VBarea, value.name="vb"))
names(dfe)<-c("indmonth", "indyr","tstp","area","vb")
dfe1<-dfe[dfe$indyr>71,]
dfe1$area<-as.factor(as.numeric(dfe1$area)+29)
head(dfe1)


aggregate(dfe1$vb,by=list(dfe1$area,dfe1$indmonth), median)

p<-ggplot(dfe1, aes(x=(area), y=vb))
p<-p+geom_boxplot()
p<-p+geom_vline(xintercept = 49-29)
p<-p+facet_wrap(~indmonth)
p

boxplot(sim$VBarea[which(indmonth==7)[71:100],])

plot(apply(sim$VBarea[which(indmonth==7)[71:100],],2,median),type="l")



dfe<-melt(sim$VBarea[which(indmonth==7)[71:100],])
head(dfe)

p<-ggplot(dfe, aes(x=as.numeric(Var2), y=value))
p<-p+geom_line(aes(color=as.factor(Var1)))
p

#doing it for a series of simulations

.SIMNAME[[1]]   <- list.files(.SIMDIRS[i],pattern="\\.Rdata", full.name=TRUE)
  
  length(sims)
    length(sims[[1]])
    names(sims[[1]])
    length(sims[[2]])
    names(sims[[2]])

    sims[[1]]$phiE

finaldf<-NULL
for( j in 1:length(.SIMNAME[[i]])){
    load(.SIMNAME[[1]][j])

    
    indmonth<-rep(1:12,100)

    

    for(m in 1:12){
      tmpsim<-apply(sims[[1]]$"VBarea"[which(indmonth==m)[71:100],1:26],2,median)
      tmpest<-apply(sims[[2]]$"VBarea"[which(indmonth==m)[71:100],1:26],2,median)
      
      tsdf<-data.frame(value=tmpsim, month=m, src="simulated", area=30:55)
      tedf<-data.frame(value=tmpest, month=m, src="estimated", area=30:55)

      finaldf<-rbind(finaldf,tsdf,tedf)
    }
   }

dim(finaldf)
head(finaldf)
summary


p<-ggplot(finaldf, aes(x=as.factor(area), y=value))
p<-p+geom_boxplot( aes(color=src, fill=src))
p<-p+ facet_wrap(~month)
p<-p+ geom_vline(aes(xintercept=which(30:55==49)))
p




########

#old stuff 



.PWD        <- "/Users/catarinawor/Documents/Lagrangian/"
.THEME      <- theme_bw(16)

.SIMDIRS   <- list.dirs(paste(.PWD,"simeval",sep=""), full.name=TRUE)







########################################################################
#simeval one scenario


#=============================================================
#gtg

#parameter estimates
.SIMDIRS   <- c(#"/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_3areas_tau04",
                #"/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_3areas_tau04_delta2",
                #"/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_3areas_tau1",
                #"/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_3areas_tau1_delta2",
                #"/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_5areas_tau04",
                #"/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_5areas_tau04_delta2",
               #"/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_5areas_tau1",
               #"/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_5areas_tau1_delta2" ,  
                "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_3areas_tau04",
                "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_3areas_tau04_delta2",
               # "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_3areas_tau1",
                "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_3areas_tau1_delta2",
                "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_5areas_tau04",
                "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_5areas_tau04_delta2",
               "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_5areas_tau1",
               "/Volumes/3T_dom media/Catarina/new_simeval_lag/SimResult_gtg_5areas_tau1_delta2"              
  )



.SIMNAME<-list(length(.SIMDIRS))

estn<-list(length(.SIMDIRS))
pbias<-list(length(.SIMDIRS))
maxgrad<-list(length(.SIMDIRS))
initvals<-list(length(.SIMDIRS))
initvals_bad<-list(length(.SIMDIRS))



for( i in 1:length(.SIMDIRS)){
  .SIMNAME[[i]]   <- list.files(.SIMDIRS[i],pattern="\\.Rdata", full.name=TRUE)
  
  tmp_estn<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=5)
  tmp_pbias<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=5)
  tmp_maxgrad<-vector(length=length(.SIMNAME[[i]]))
  tmp_initvals<-matrix(NA,nrow=length(.SIMNAME[[i]]),ncol=5)
  

  for( j in 1:length(.SIMNAME[[i]])){
    load(.SIMNAME[[i]][j])

    true_pars <- c(sims[[1]]$"mo",sims[[1]]$"cvPos",sims[[1]]$"maxPos50",sims[[1]]$"maxPossd",mean(sims[[1]]$Fmult))  


    #parameters
    tmp_estn[j,]<-exp(sims[[3]]$est[1:5])
    tmp_pbias[j,]<-((tmp_estn[j,]-true_pars)/true_pars)
    tmp_maxgrad[j]<-sims[[3]]$maxgrad
    tmp_initvals[j,]<-exp(unlist(sims[[5]][1:5]))
   }

  tmp_estn<- tmp_estn[tmp_maxgrad<=1.,]
  tmp_pbias<- tmp_pbias[tmp_maxgrad<=1.,]
  
  estn[[i]]<-tmp_estn
  pbias[[i]]<-tmp_pbias
  maxgrad[[i]]<-tmp_maxgrad
  initvals[[i]]<-tmp_initvals[tmp_maxgrad<=1.,]
  initvals_bad[[i]]<-tmp_initvals[tmp_maxgrad>1.,]


}


titulos<-c("3 areas, tau=0.4, B = 1.0",
            "3 areas, tau=0.4, B = 2.0",
           # "3 areas, tau=1.0, B = 1.0",
            "3 areas, tau=1.0, B = 2.0",
          "5 areas, tau=0.4, B = 1.0",  
           "5 areas, tau=0.4, B = 2.0",
            "5 areas, tau=1.0, B = 1.0",
           "5 areas, tau=1.0, B = 2.0"
  )

#setwd("/Users/catarinawor/Documents/hake/Thesis/figs/chap2")
#setwd("/Users/catarinawor/Documents/hake/JTC_talk")
#pdf("single_version_simeval.pdf", width=14, height=7)
setwd("/Users/catarinawor/Documents/hake/Lag_Model_paper/")
pdf("Figure3.pdf", width=14, height=7)
simscn<-1:8

par(mfcol=c(2,4))
for( i in 1:length(.SIMDIRS)){
  boxplot(pbias[[i]],names=c(expression("t"[0]),expression("CV"),expression("a"[50]),
    expression(sigma["X"["max"]]),expression("q")),ylim=c(-50,50),main=titulos[i],cex.axis=1.5,
    cex.lab=2,cex.main=2,cex=1.6)
  abline(h=0)
 # text(4, y = 8, labels = nrow(pbias[[i]]), cex=2)
  text(4, y = 35, labels = gtgscn[i], cex=2)
}
mtext("% relative error", 2, line = -2, outer = TRUE, font=2)
mtext("Parameter", 1, line = -2, outer = TRUE, font=2)
#dev.off()



