## Makefile for cleaning up
## I don't think it works on pc! 
## Targets: 
##		all:   -copy executable and run the model with DAT & ARG
##		run:   -copy executable and force a run
##		mcmc:  -copy executable and run in mcmc mode and mceval
##		retro: -copy executable and run  retrospective analysis



EXEC=lagrangian_est
DAT=lagrangian_est.dat
OM=lagrangian_OM
SIMDAT=lagrangian_OM.dat
LAST := 3
NUMBERS := $(shell seq 1 ${LAST})

RDATA ='source(file.path("../","Lagrangian","simRun.r"))'

 
all: 
	compiOM rom compi run

compiOM: $(OM)
	admb $(OM)

compi: $(EXEC)
	admb $(EXEC)

rom: $(OM)
	-./$(OM) -ind $(SIMDAT) $(ARG)


run: $(EXEC)
	./$(EXEC) -ind $(DAT) $(ARG)

readRdat: 
	-@echo $(RDATA) | R --vanilla --slave


doall: rom run readRdat

simeval: $(foreach var,$(NUMBERS),doall)
	
        
clean: 
	-rm -f  admodel.* variance eigv.rpt fmin.log $(EXEC) variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor  


cleanOM: 
	-rm -f  admodel.* variance eigv.rpt fmin.log $(OM) variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor *.par 
