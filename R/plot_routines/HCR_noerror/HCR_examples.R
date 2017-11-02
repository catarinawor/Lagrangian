#========================================================================================
# Harvest control rules schematic plot
#========================================================================================

#need to run with no effort

Bo<-100


yB<-0:100

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







