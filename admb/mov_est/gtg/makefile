## Makefile for cleaning up
## I don't think it works on pc! 
## Targets: 
##		all:   -copy executable and run the model with DAT & ARG
##		run:   -copy executable and force a run
##		mcmc:  -copy executable and run in mcmc mode and mceval
##		retro: -copy executable and run  retrospective analysis

## when @ is put at the beggining of a recipe, it suppressing the echoing of that recipe
## i.e. it does not print it on the terminal.
##variables

EXEC=lagrangian_est_gtg
DAT=lagrangian_est_gtg.dat
SEED = seed.txt


run: $(EXEC) $(DAT) $(ARG) 
	-./$(EXEC) -ind $(DAT) $(ARG)


all: 
	compi run


compi: $(EXEC)
	admb $(EXEC)






FORCE: 

cleanall: 
	-rm -f  admodel.* variance eigv.rpt fmin.log $(EXEC) variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor *.obj *.p02 *.b02 *.par *.r02 *.rep 

        
clean: 
	-rm -f  admodel.* variance eigv.rpt fmin.log $(EXEC) variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor *.obj *.p02 *.b02 *.par *.r02


cleanOM: 
	-rm -f  admodel.* variance eigv.rpt fmin.log $(OM) variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor *.par *.obj
