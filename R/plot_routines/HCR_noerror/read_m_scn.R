#==============================================================
# Read in the results from the HCR search - no assessment error
#Author: Catarina Wor
#Date: Sep 14 2017

#==============================================================


library(cowplot)



DIRS<-c("/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_linear_nostRec_100",
"/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_linear_nostRec_100_m2",
"/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_linear_nostRec_100_m3",
"/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_stquo_nostRec_100",
"/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_stquo_nostRec_100_m2",
"/Users/catarinawor/Documents/Lagrangian/HCRresult/HCR_stquo_nostRec_100_m3"
)




Rfiles<-list()
SIMSdat<-list()

for(d in 1:length(DIRS)){

	Rfiles[[d]]<-list.files(DIRS[[d]],pattern="\\.Rdata",full.name=TRUE)

	tmp_SIMS<-list()
	for(i in 1:length(Rfiles[[d]])){
	

		load(Rfiles[[d]][i])
	
		tmp_SIMS[[i]]<-sims


	}
	SIMSdat[[d]]<-tmp_SIMS

}


plotlib<-"/Users/catarinawor/Documents/Lagrangian/R/plot_routines/HCR_noerror"
plotfiles <- list.files(plotlib,pattern="plot",full.name=TRUE)


#load graphing routinesfor
for(p in 1:length(plotfiles)){
	source(plotfiles[p])
}

M<-SIMSdat[[1]]

M1<-SIMSdat[[1]]
M2<-SIMSdat[[2]]
M3<-SIMSdat[[3]]
Msq<-SIMSdat[[4]]
Msq1<-SIMSdat[[4]]
Msq2<-SIMSdat[[5]]
Msq3<-SIMSdat[[6]]


draft_dir<-"/Users/catarinawor/Documents/Lagrangian/report/manuscript/ICES"
#normal<-plot_logUtility( SIMSdat,SIMSsq , sv=TRUE, nome="nostrongrec",nations=TRUE)
no<-plot_logUtility( SIMSdat[[1]],SIMSdat[[4]] , sv=FALSE, nome="nostrongrec",nations=TRUE)
one<-plot_logUtility( SIMSdat[[2]],SIMSdat[[5]]  , sv=FALSE, nome="onestrongrec",nations=TRUE)
two<-plot_logUtility( SIMSdat[[3]],SIMSdat[[6]]  , sv=FALSE, nome="twostrongrec",nations=TRUE)
lu<-plot_grid(no,one,two,nrow=3, labels="AUTO")
lu
setwd("/Users/catarinawor/Documents/Lagrangian/figures/HCR")
ggsave("logUtility_allScn.pdf", plot=lu, width = 15, height = 13,dpi = 600)



noAAV<-plot_AAV( SIMSdat[[1]],SIMSdat[[4]] , sv=TRUE, nome="AAVnostrongrec",nations=TRUE)
oneAAV<-plot_AAV( SIMSdat[[2]],SIMSdat[[5]]  , sv=TRUE, nome="AAVonestrongrec",nations=TRUE)
twoAAV<-plot_AAV( SIMSdat[[3]],SIMSdat[[6]]  , sv=TRUE, nome="AAVtwostrongrec",nations=TRUE)
AAV<-plot_grid(noAAV,oneAAV,twoAAV,nrow=3, labels=c("no strong recruitment","one strong recruitments","two strong recruitments"),
	  hjust = 0)
AAV
setwd(draft_dir)
ggsave("AAV_allScn.pdf", plot=AAV, width = 15, height = 15,dpi = 600)



nocl<-plot_closures( SIMSdat[[1]],SIMSdat[[4]] , sv=TRUE, nome="clnostrongrec",nations=TRUE)
onecl<-plot_closures( SIMSdat[[2]],SIMSdat[[5]]  , sv=TRUE, nome="clonestrongrec",nations=TRUE)
twocl<-plot_closures( SIMSdat[[3]],SIMSdat[[6]]  , sv=TRUE, nome="cltwostrongrec",nations=TRUE)
cl<-plot_grid(nocl,onecl,twocl,ncol=3, labels="AUTO",rel_widths = c(1,1, 1))
cl
setwd(draft_dir)
ggsave("closure_allScn.pdf", plot=cl, width = 15, height = 6,dpi = 600)



noY<-plot_Yield( SIMSdat[[1]],SIMSdat[[4]] , sv=FALSE, nome="clnostrongrec",nations=TRUE)
oneY<-plot_Yield( SIMSdat[[2]],SIMSdat[[5]]  , sv=FALSE, nome="clonestrongrec",nations=TRUE)
twoY<-plot_Yield( SIMSdat[[3]],SIMSdat[[6]]  , sv=FALSE, nome="cltwostrongrec",nations=TRUE)
Yield<-plot_grid(noY,oneY,twoY,nrow=3, labels="AUTO",rel_widths = c(1,1, 1))
Yield
setwd("/Users/catarinawor/Documents/Lagrangian/figures/HCR")
ggsave("Yield_allScn.pdf", plot=Yield, width = 15, height = 13,dpi = 600)


Bno<-plot_B40( SIMSdat[[1]],SIMSdat[[4]], sv=FALSE, nome="nostrongrec",nations=TRUE)
Bone<-plot_B40( SIMSdat[[2]],SIMSdat[[5]], sv=FALSE, nome="nostrongrec",nations=TRUE)
Btwo<-plot_B40( SIMSdat[[3]],SIMSdat[[6]], sv=FALSE, nome="nostrongrec",nations=TRUE)

b40<-plot_grid(Bno,Bone,Btwo,ncol=3, labels="AUTO")
b40



setwd("/Users/catarinawor/Documents/Lagrangian/figures/HCR")
ggsave("B40_allScn.pdf", plot=b40, width = 15, height = 6,dpi = 600)

plot_B40( SIMSdat[[1]],SIMSdat[[4]], sv=FALSE, nome="nostrongrec",nations=TRUE, limite=.1)



Blimno<-plot_B40( SIMSdat[[1]],SIMSdat[[4]], sv=FALSE, nome="nostrongrec",nations=TRUE,limite=.2)
Blimone<-plot_B40( SIMSdat[[2]],SIMSdat[[5]], sv=FALSE, nome="nostrongrec",nations=TRUE,limite=.2)
Blimtwo<-plot_B40( SIMSdat[[3]],SIMSdat[[6]], sv=FALSE, nome="nostrongrec",nations=TRUE,limite=.2)

b20<-plot_grid(Blimno,Blimone,Blimtwo,ncol=3, labels="AUTO")
b20


plot_tradeoff(SIMSdat[[1]],SIMSdat[[2]],SIMSdat[[3]],SIMSdat[[4]],SIMSdat[[5]],SIMSdat[[6]],sv=TRUE)

