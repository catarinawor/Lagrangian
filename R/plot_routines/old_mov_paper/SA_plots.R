
#=======================================================================
# Graphing routine for selectivity chapter - Is iSCAM recobering new parameters? 
#
#
#=======================================================================
rm(list=ls()); 


setwd("/Users/catarinawor/Documents/Lagrangian/R")
source("read.admb.R")

setwd("/Users/catarinawor/Documents/Lagrangian/admb/OM/simple")
sim <- read.rep("lagrangian_OM.rep")

#import data from iSCAM
setwd("/Users/catarinawor/Documents/iSCAM/examples/hakelag/DATA")
est <- read.rep("hakelag.rep")

names(sim)
names(est)

sim$totcatch
est$ct


mean(sim$Nage[sim$indmonth==1,1][71:100])


plot(est$yr,sim$itB, type="l", lty=2)
lines(est$yr,est$bt[1:30],lwd=2)



plot(est$yr,sim$SB[sim$indmonth==1][71:100], type="l", ylim=)
lines(est$yr,est$sbt[1:30],lwd=2, col="darkblue")

length(sim$SB)

Fi<-seq(0.000,4.00,by=0.001)
plot(Fi,est$allspr)
plot(Fi,est$diffspr)


est$allspr[1:30]
est$spr
est$fspr

sim$h
est$steepness

sim$Ro
est$ro

est$ro/sim$Ro



est$d3_survey_data[,2]
est$it_hat

plot(est$d3_survey_data[,1],est$d3_survey_data[,2])
lines(est$d3_survey_data[,1],est$it_hat,lwd=2)

#=======================================================================
#Plotting simulation results
#=======================================================================

setwd("/Users/catarinawor/Documents/Lagrangian/")
source("read.admb.R")


.SIMDIRS   <- c("/Users/catarinawor/Documents/Lagrangian/SAsim1")

.SIMNAME<-list(length(.SIMDIRS))

steepness_list<-list(length=length(.SIMDIRS))
ro_list<-list(length=length(.SIMDIRS))
ut_list<-list(length=length(.SIMDIRS))
ct_list<-list(length=length(.SIMDIRS))

for( i in 1:length(.SIMDIRS)){
  .SIMNAME[[i]] <- list.files(.SIMDIRS[i],pattern="\\.Rdata")


  steepness_mat<-matrix(NA, ncol=3, nrow=length(.SIMNAME[[i]]))
  ro_mat<-matrix(NA, ncol=3, nrow=length(.SIMNAME[[i]]))
  ut_array<-array(NA,c(length(.SIMNAME[[i]]),length(est$yr),3))
  ct_array<-array(NA,c(length(.SIMNAME[[i]]),length(est$yr),3))
  
  setwd(.SIMDIRS[i])
  for( j in 1:length(.SIMNAME[[i]])){

    load(.SIMNAME[[i]][j])

    steepness_mat[j,1]<-sims$sim$h
    steepness_mat[j,2]<-sims$est$steepness
    steepness_mat[j,3]<-((steepness_mat[j,2]-steepness_mat[j,1])/steepness_mat[j,1])*100

    ro_mat[j,1]<-sims$sim$Ro
    ro_mat[j,2]<-sims$est$ro
    ro_mat[j,3]<-((ro_mat[j,2]-ro_mat[j,1])/ro_mat[j,1])*100

    ut_array[j,,1]<-sims$sim$Ut
    ut_array[j,,2]<-sims$est$ut
    ut_array[j,,3]<-((ut_array[j,,2]-ut_array[j,,1])/ut_array[j,,1])*100

    ct_array[j,,1]<-sims$sim$totcatch
    ct_array[j,,2]<-sims$est$ct
    ct_array[j,,3]<-((ct_array[j,,2]-ct_array[j,,1])/ct_array[j,,1])*100


	}
	steepness_list[[i]]<-steepness_mat
	ro_list[[i]]<-ro_mat
	ut_list[[i]]<-ut_array
	ct_list[[i]]<-ct_array

}

steepness_list[[1]]

par(mfrow=c(2,2))
boxplot(steepness_list[[1]][,3], main="steepness")
abline(h=0,lwd=2,col="darkred")
boxplot(ro_list[[1]][,3], main="Ro")
abline(h=0,lwd=2,col="darkred")
boxplot(ut_list[[1]][,,3], main="Ut")
abline(h=0,lwd=2,col="darkred")
boxplot(ct_list[[1]][,,3], ylim=c(-0.01,0.8), main="Ct")
abline(h=0,lwd=2,col="darkred")


#==============================================================
#plots for the sim
#==============================================================

rm(list=ls()); 

setwd("/Users/catarinawor/Documents/Lagrangian/")
source("read.admb.R")

library(plyr)
library(data.table)
library(ggplot2)
library(reshape2)

sim <- read.rep("lagrangian_OM.rep")


names(sim)

head(sim$selnation)


tmpdf<-data.frame(nation=sim$selnation[,1],year=sim$selnation[,2],
  age1=sim$selnation[,3],age2=sim$selnation[,4],age3=sim$selnation[,5],
  age4=sim$selnation[,6],age5=sim$selnation[,7],age6=sim$selnation[,8],
  age7=sim$selnation[,9],age8=sim$selnation[,10],age9=sim$selnation[,11],
  age10=sim$selnation[,12],age11=sim$selnation[,13],age12=sim$selnation[,14],
  age13=sim$selnation[,15],age14=sim$selnation[,16],age15=sim$selnation[,17],
  age16=sim$selnation[,18],age17=sim$selnation[,19],age18=sim$selnation[,20],
  age19=sim$selnation[,21],age20=sim$selnation[,22])

tmpdf=tmpdf[tmpdf$year>70,]

head(tmpdf)

a<-melt(tmpdf,id.vars=c("nation","year"))
head(a)
?melt

p <- ggplot( a,aes(year,variable,size=value, col=factor(nation) ))
  p <- p + geom_point(alpha=0.75)
  # p <- p + scale_area(range = c(0,10))
  p 


 
seln1<-sim$selnation[sim$selnation[,1]==1,-c(1,2)]
seln1<-seln1[71:100,-c(11:20)]
seln1<-(t(seln1))
head(seln1)
dim(seln1)

seln2<-sim$selnation[sim$selnation[,1]==2,-c(1,2)]
seln2<-seln2[71:100,-c(11:20)]
seln2<-(t(seln2))

par(mfrow=c(1,2))
matplot(1:10,seln1*10,type="n",ylim=c(71,110))
matlines(1:10,(seln1*10+rep(seq(71,100, len=30), each=10)), lty=1,col=1,lwd=2)
matplot(1:10,seln2*10,type="n",ylim=c(71,110))
matlines(1:10,seln2*10+rep(seq(71,100, len=30), each=10), lty=1,col=1,lwd=2)

source("utilities.r")
install.packages.if.needed("r4ss", "r4ss/r4ss", github=TRUE)
library(r4ss)

par(mfrow=c(1,2))
mountains(zmat=t(seln1), yvec=71:100)
mountains(zmat=t(seln2), yvec=71:100)

#==============================================================
#SPR plotting -  plots for equilibrium conditions 
#==============================================================


rm(list=ls()); 

setwd("/Users/catarinawor/Documents/Lagrangian/R")
source("read.admb.R")

library(plyr)
library(data.table)
library(ggplot2)
library(reshape2)

setwd("/Users/catarinawor/Documents/Lagrangian/admb/OM/simple")
sim <- read.rep("lagrangian_OM.rep")

names(sim)

par(mfrow=c(1,2))
plot(apply(sim$yNage,1,sum), type="b",lwd=2)
plot(sim$spr, type="b",lwd=2)




setwd("/Users/catarinawor/Documents/Lagrangian/admb/OM/gtg")
sim_gtg<-read.rep("lagrangian_OM_gtg.rep")

names(sim_gtg)

par(mfrow=c(1,2))
plot(apply(sim_gtg$yNage,1,sum), type="b",lwd=2)
plot(sim_gtg$spr, type="b",lwd=2)


setwd("/Users/catarinawor/Documents/Lagrangian/admb/stock_assessment")
SA<-read.rep("lagrangian_SA.rep")

sprSA<-read.rep("spr.rep")

names(sprSA)
fs<-seq(length=4001, by=0.001)
plot(fs,sprSA$diffspr)


names(SA)

SA$fspr
SA$spr_opt


par(mfrow=c(1,2))
plot(apply(sim_gtg$yNage,1,sum), type="b",lwd=2)
plot(sim_gtg$spr, type="b",lwd=2)






