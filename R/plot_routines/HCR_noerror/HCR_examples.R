#========================================================================================
# Harvest control rules schematic plot
#========================================================================================

#need to run with no effort



SBo<-4421.47


ySB<-seq(0,4421.47,length=100)
yB<-ySB*1.239169


slope=.2
slope2=.2
slope3=.5

intercept=0
intercept2= 0.1
intercept3=.4

utarget=0.15

TAC<-NULL
TAC2<-NULL
TAC3<-NULL

TAC1040<-NULL

for(i in 1:length(yB)){

	if(ySB[i]/SBo>intercept){
		TAC[i]= max(min(slope*yB[i]*(ySB[i]-SBo*intercept)/SBo,600),0)
	}else{
		TAC[i]= 0
	}

	if(ySB[i]/SBo>intercept2){
		TAC2[i]= max(min(slope2*yB[i]*(ySB[i]-SBo*intercept2)/SBo,600),0)
	}else{
		TAC2[i]= 0
	}

	if(ySB[i]/SBo>intercept2){
		TAC3[i]= max(min(slope3*yB[i]*(ySB[i]-SBo*intercept3)/SBo,600),0)
	}else{
		TAC3[i]= 0
	}

	if(ySB[i]/SBo>.4){
		TAC1040[i]= min(utarget*yB[i],600)
	}else{
	 	if(ySB[i]/SBo>.1){
	 		TAC1040[i]=min(((ySB[i]-0.1*SBo)*((0.4/ySB[i])/(0.4-0.1)))* utarget*yB[i],600)
		}else{
			TAC1040[i]=0
		}
	}
}


#pdf("Figurec4_HCRex.pdf",width = 10, height = 8)

plot(ySB/SBo,TAC,type="l",lwd=2, xlab= expression(SB[t]/SB[0]))
text(x=0.32,75, labels ="slope = 0.2", adj=0.0,font=2)
text(x=0.32,100, labels ="intercept = 0.0", adj=0.0,font=2)

lines(ySB/SBo,TAC3, lwd=2, lty=2,col="gray40")
text(x=0.55,175, labels ="slope = 0.5", adj=0.0,col="gray40",font=2)
text(x=0.55,200, labels ="intercept = 0.4", adj=0.0,col="gray40",font=2)

lines(yB/SBo,TAC1040, lwd=2, lty=3,col="gray20")
text(x=0.15,175, labels ="40:10 rule", adj=0.0,col="gray20",font=2)

#dev.off()




par(las=1, mar=c(5,5,1,1))
plot(ySB/SBo,TAC,type="l",lwd=5, xlab= expression(SB[t]/SB[0],font=2),
 col="darkorange2", cex.lab=1.4,cex.axis=1.3,font.lab=2,font.axis=2)
lines(ySB/SBo,TAC3, lwd=5, lty=2,col="seagreen")
lines(yB/SBo,TAC1040, lwd=5, lty=1,col="dodgerblue4", cex=1.4)

text(x=0.32,75, labels ="ER = 0.2", adj=0.0,font=2, col="darkorange2", cex=1.4)
text(x=0.32,100, labels ="threshold = 0.0", adj=0.0,font=2, col="darkorange2", cex=1.4)

text(x=0.55,175, labels ="ER = 0.5", adj=0.0,col="seagreen",font=2, cex=1.4)
text(x=0.55,200, labels ="threshold = 0.4", adj=0.0,col="seagreen",font=2, cex=1.4)

text(x=0.18,175, labels ="40:10 rule", adj=0.0,col="dodgerblue4",font=2, cex=1.4)




#lines(ySB/SBo,TAC2, lwd=2,col="gray30")
yint=-0.2*Bo*0.4
TAC<- yint + 0.4*yB

yint2=-0.*Bo*1
TAC2<- yint2 + 1*yB

yint3=-0.4*Bo*1
TAC3<- yint3 + 1*yB


setwd("/Users/catarinawor/Documents/Lagrangian/report/HCR_coarse")
pdf("HCR_example.pdf", width=6,height=6)
par(xaxt="n",bty="l")
plot(yB,TAC, xlim=c(0,100), ylim=c(0,100),yaxs="i",xaxs="i", type="l", lwd=2, xlab=expression(B[t]/B[0]))
text(x=65,13, labels ="slope = 0.4", adj=0.0)
text(x=65,16, labels ="intercept = 0.2", adj=0.0)

lines(yB,TAC2,lwd=2)
lines(yB,TAC2,lwd=2,lty=2, col="gray70")
text(x=27,42, labels ="slope = 1.0", adj=0.0)
text(x=27,45, labels ="intercept = 0.0", adj=0.0)

lines(yB,TAC3,lwd=2)
lines(yB,TAC3,lwd=2,lty=3, col="gray90")
text(x=69,44, labels ="slope = 1.0", adj=0.0)
text(x=69,47, labels ="intercept = 0.4", adj=0.0)
axis(1, at = seq(0,100,by=20), labels = seq(0,1.0,by=.2),xaxt="s")
dev.off()




?par







