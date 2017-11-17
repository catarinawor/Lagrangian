#========================================================================================
# Generate input data for Lagrangian_gtg_OM
#========================================================================================





writeHCRdat<-function(DIR= "../../admb/OM/gtg/", slope, intercept){


	setwd(DIR)
	sink("HCR.dat")
	#cat("#Harvest control rules\n")
	cat("#SAtype\n")
	cat("2 \n")
	cat("#nationTACprop\n")
	cat("0.7388 0.2612 \n")
	cat("#hcr\n")
	cat(2 ,"\n")
	cat("#hcr_slope\n")
	cat(slope, "\n")
	cat("#hcr_intercept\n")
	cat(intercept, "\n")
	cat("#eof\n")
	cat(999, "\n")
	sink()



}










