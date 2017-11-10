#========================================================================================
# Generate input data for Lagrangian_gtg_OM
#========================================================================================







writeWTdat<-function(DIR= "/Users/catarinawor/Documents/Lagrangian/admb/OM/gtg/", myseed){

	setwd(DIR)

	twt<-scan(file="fixed_wt_scn.dat",comment.char = "#")


	set.seed(myseed)
	perr<-twt[50:length(twt)]*exp(rnorm(length(twt[50:length(twt)]),0,0.1))

	rec_ts<-c(twt[1:49],perr)


	setwd(DIR)
	sink("wt.dat")
	cat("#recruitment deviations\n")
	cat(rec_ts,"\n")
	cat("#eoffw\n")
	cat(999, "\n")
	sink()



}










