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

OM=lagrangian_OM
SIMDAT=lagrangian_OM.dat
SEED = ../../seed.txt
LAST := 100
SADIR= /Users/catarinawor/Documents/iSCAM/examples/hakelag/DATA


all: 
	compiOM rom 

compiOM: $(OM)
	admb $(OM)


rom: $(OM) $(SIMDAT) $(ARG) $(SEED)
	-./$(OM) -ind $(SIMDAT) $(ARG)

 

clean: 
	-rm -f  admodel.* variance eigv.rpt fmin.log variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor *.obj *.p02 *.b02 *.par *.r02 *.rep $(OM).cpp

        
