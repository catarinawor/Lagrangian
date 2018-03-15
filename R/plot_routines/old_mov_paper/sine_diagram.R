#=====================================================
#Produces a diagram of the dynamics of the sine curve 
#author: Catarina Wor
#Date June 2017
#=====================================================


#=====================================================



x <- seq(1,24,by=0.1)
tt<-1:12
y <- 30+(60-30)*(0.5+0.5*sin(x*2*pi/length(tt)-2*pi/length(tt)-pi/2))

setwd("/Users/catarinawor/Documents/hake/Lag_Model_paper")
pdf("sin_diagram.pdf", width=12,height=6)
par(bty="l",yaxs="i",las=1)
plot(x,y,type="l", ylab=" ", xlab="",lwd=3,
	ylim=c(28,62.5),axes=FALSE)
mtext("X spatial gradient", side=2,las=0, font=2,cex=1.5)
arrows(22.5,45,24,45,length=0.1,lwd=3)
text(x=23.5, y = 46, labels = "Time", font=2,cex=1.5)
arrows(x[length(x)-1],y[length(x)-1],x[length(x)]
     ,y[length(x)],length=0.1,lwd=3,xaxp=c(1,24,23))
axis(1,pos=45,xaxp=c(1,24,23),hadj=.1)
axis(2,pos=1,font=2)
arrows(1,45,1.1,45.5,length=0.0,lwd=3,col="gray50")
text(x=1.2, y = 46.1, labels = expression(paste("t"[0])),font=2,cex=1.5)
arrows(13,45,13,45.5,length=0.0,lwd=3,col="gray50")
text(x=13.0, y = 46.1, labels = expression(paste("t"[0])), font=2,cex=1.5)
arrows(7,60.5,7,45.,length=0.0,lwd=3,col="gray50")
text(x=7.1, y = 61.5, labels = expression(paste("X"[max])), font=2,cex=1.5)
arrows(10,45,10,45.5,length=0.0,lwd=3,col="gray50")
text(x=10.4, y = 46, labels = expression(paste("X"[mid])), font=2,cex=1.5)
arrows(13,30,13,45,length=0.0,lwd=3,col="gray50")
text(x=13.1, y = 29, labels = expression(paste("X"[min])), font=2,cex=1.5)
arrows(16,45,16,45.5,length=0.0,lwd=3,col="gray50")
text(x=15.5, y = 46, labels = expression(paste("X"[mid])), font=2,cex=1.5)
arrows(19,60.5,19,45,length=0.0,lwd=3,col="gray50")
text(x=19.1, y = 61.5, labels = expression(paste("X"[max])), font=2,cex=1.5)
dev.off()


meses<-c("Jan", "Feb", "Mar", "Apr","May", "Jun","Jul","Aug","Sep", "Oct","Nov", "Dec")



par(bty="l",yaxs="i",las=1, mar=c(7.1, 5.1, 3.1, 3.1))
plot(x,y,type="l", ylab=" ", xlab="",lwd=3,
	ylim=c(28,64),axes=FALSE)
mtext("X spatial gradient", side=2,las=0, font=2,cex=2, line=3)
axis(1,pos=30,xaxp=c(1,24,23),hadj=.1, labels=rep(meses,2), at=1:24, cex.axis=2, font=2)
mtext("Time", side = 1, line = 3, outer = FALSE, at = NA,
      adj = NA, padj = NA, cex = 2, font = 2)
#text(x=23.5, y = 46, labels = "Time", font=2,cex=1.5)
axis(2,pos=1,font=2,cex.axis=2,at=seq(30,60,5))
arrows(x[length(x)-1],y[length(x)-1],x[length(x)]
     ,y[length(x)],length=0.1,lwd=3,xaxp=c(1,24,23))
arrows(7,60.5,7,58.,length=0.0,lwd=3,col="gray50")
text(x=7.1, y = 61.9, labels = expression(paste("X"[max])), font=2,cex=2.5)
arrows(13,30,13,31,length=0.0,lwd=3,col="gray50")
text(x=13.1, y = 33, labels = expression(paste("X"[min])), font=2,cex=2.5)
arrows(19,60.5,19,58.,length=0.0,lwd=3,col="gray50")
text(x=19.1, y = 61.9, labels = expression(paste("X"[max])), font=2,cex=2.5)

