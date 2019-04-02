


myfunc<-function(maxPos50,maxPossd){

	age<-0:20

	x<-(1./(1.+exp(-(age-maxPos50)/maxPossd)))

	x<-x* (60-32)+32

	return(x)
	
}


setwd("/Users/catarinawor/Documents/Thesis/")
			
pdf("Figurec4_movscn.pdf",width = 10, height = 8)
par(mar=c(4.1,6,2.1,2.1))
plot(myfunc(4,2), type="b", ylab= expression(paste(bar(X)[max*paste(",")~a]," (Latitude)")), 
	xlab="Age",ylim=c(32,60),lwd=2,cex.lab=1.5,cex.axis=1.3)
text(13,43,"Base",font=2,cex=1.2)
text(13,42,expression(paste(a[50]==4,", ", sigma[X[max]] == 2)),font=2,cex=1.2)
lines(myfunc(2.5,1), col="gray50", type="b",lwd=2)
text(13,40.5,"Early movement",font=2, col="gray50",cex=1.2)
text(13,39.5,expression(paste(a[50]==2.5,", ", sigma[X[max]] == 1)),font=2,col="gray50",cex=1.2)
lines(myfunc(5.5,3), col="gray30", type="b",lwd=2)
text(13,38,"Late movement",font=2, col="gray30",cex=1.2)
text(13,37,expression(paste(a[50]==5.5,", ", sigma[X[max]] == 3)),font=2,col="gray30",cex=1.2)
abline(h=49, col="gray70")
text(18,49.6,"Canada",font=2, col="gray70",cex=1.2)
text(18,48.4,"U.S.A",font=2, col="gray70",cex=1.2)
dev.off()



par(mar=c(4.1,6,2.1,2.1))
plot(myfunc(4,2), type="b", ylab= expression(paste(bar(X)[max*paste(",")~a]," (Latitude)")), 
	xlab="Age",ylim=c(32,60),lwd=5,cex.lab=1.5,cex.axis=1.3,col="darkorange2")
text(13,45,"Base",font=2,cex=1.4, col="darkorange2")
text(13,43.5,expression(paste(a[50]==4,", ", sigma[X[max]] == 2)),font=2,cex=1.6, col="darkorange2")
lines(myfunc(2.5,1), col="seagreen", type="b",lwd=5)
text(13,41,"Early movement",font=2, col="seagreen",cex=1.4)
text(13,39.5,expression(paste(a[50]==2.5,", ", sigma[X[max]] == 1)),font=2,col="seagreen",cex=1.6)
lines(myfunc(5.5,3), col="dodgerblue4", type="b",lwd=5)
text(13,37,"Late movement",font=2, col="dodgerblue4",cex=1.4)
text(13,35.5,expression(paste(a[50]==5.5,", ", sigma[X[max]] == 3)),font=2,col="dodgerblue4",cex=1.6)
abline(h=49, col="gray70")
text(18,49.6,"Canada",font=2, col="gray70",cex=1.6)
text(18,48.4,"U.S.A",font=2, col="gray70",cex=1.6)
