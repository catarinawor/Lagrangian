#=====================================================
#Produces a diagram of the dynamics of the sine curve 
#author: Catarina Wor
#Date June 2017
#=====================================================


#=====================================================



x <- seq(1,24,length.out=100)
tt<-1:12
y <- 30+(60-30)*(0.5+0.5*sin(x*2*pi/length(tt)-2*pi/length(tt)-pi/2))
par(bty="l",yaxs="i")
plot(x,y,type="l", ylab="X spatial gradient ", xlab="Time",lwd=3,
 xaxp=c(1,24,23) )
arrows(x[99],y[99],x[100]
	,y[100],length=0.1,lwd=3,xaxp=c(1,24,23))
arrows(x[1],y[1],x[1]
	,y[4],length=0.0,lwd=3,xaxp=c(1,24,23))
text(x=x[1], y = y[5], labels = expression(paste("t"[0])))


?text
?par
?arrows