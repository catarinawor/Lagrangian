#=========================================================================
# Area map
# Author: Catarina Wor
# Jul 18st 2017 
#=========================================================================


library(plyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(animation)
library(ggmap)
#if (Sys.info()["nodename"] =="sager")  setwd("~/Dropbox/LSRA/length_SRA/sim_est_lsra")
setwd("/Users/catarinawor/Documents/Lagrangian")
source("R/read.admb.R")


basemap<-get_map(location = c(lon = -130, lat = 45.5),
    zoom = 5, source = "stamen",maptype="toner-lite",color="bw") #)maptype="terrain", color="bw"

invert <- function(x) rgb(t(255-col2rgb(x))/255)    
m_inv <- as.raster(apply(basemap, 2, invert))

# copy attributes from original object
class(m_inv) <- class(basemap)
attr(m_inv, "bb") <- attr(basemap, "bb")
m_inv

basemap<-get_stamenmap(bbox = c(left = -138, bottom = 30, right = -115,
  top = 55), zoom = 5, maptype = "toner-lite", crop = TRUE, messaging = FALSE,
  color = "bw")


p2<-ggmap(basemap,extent = "panel", maprange=FALSE, crop=T) 
p2<- p2 + labs(x = 'Longitude', y = 'Latitude') 
#p2 <- p2 + geom_hline(yintercept=42, linetype=2, colour="grey60")
p2 <- p2 + geom_segment(aes(x=-138,xend=-124,y=42,yend=42),linetype=2, colour="black")
p2 <- p2 + geom_text(x=-134.5,y=41.5,label="fishing ground 1", colour="black")

#p2 <- p2 + geom_hline(yintercept=46, linetype=2, colour="grey60")
p2 <- p2 + geom_segment(aes(x=-138,xend=-124.5,y=46,yend=46),linetype=2, colour="black")
p2 <- p2 + geom_text(x=-134.5,y=45.5,label="fishing ground 2", colour="black")

#p2 <- p2 + geom_hline(yintercept=48.5, linetype=2, colour="grey60")
p2 <- p2 + geom_segment(aes(x=-138,xend=-124.8,y=48.5,yend=48.5),linetype=2, colour="black")

p2 <- p2 + geom_text(x=-134.5,y=48,label="fishing ground 3", colour="black")
p2 <- p2 + geom_segment(aes(x=-138,xend=-127.6,y=51,yend=51),linetype=2, colour="black")

#p2 <- p2 + geom_hline(yintercept=51, linetype=2, colour="grey60")
p2 <- p2 + geom_text(x=-134.5,y=50.5,label="fishing ground 4", colour="black")
p2 <- p2 + geom_text(x=-134.5,y=53.5,label="fishing ground 5", colour="black")
p2 <- p2 + theme_bw(16)
p2 <- p2 + geom_text(x=-120,y=51,label="Canada", colour="black", size=6)
p2 <- p2 + geom_text(x=-120,y=43,label="U.S.A", colour="black", size=6)

#-127.620
p2
setwd("/Users/catarinawor/Documents/hake/Lag_Model_paper")
ggsave("study_area.pdf", plot = p2)
