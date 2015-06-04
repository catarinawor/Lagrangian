## Makefile for cleaning up
## I don't think it works on pc! 
## Targets: 
##		all:   -copy executable and run the model with DAT & ARG
##		run:   -copy executable and force a run
##		mcmc:  -copy executable and run in mcmc mode and mceval
##		retro: -copy executable and run  retrospective analysis



EXEC=lagrangian_est
DAT=lagrangian_est.dat
CTL=
ARG=
MCFLAG=-mcmc 10000 -mcsave 100 -nosdmcmc
NR=4
OM=lagrangian_OM
SIMDAT=lagrangian_OM.dat
LAST := 1
NUMBERS := $(shell seq 1 ${LAST})


 
all: $(EXEC)
	./$(EXEC) -ind $(DAT) $(ARG)

compi: admb $(EXEC) 


run: ./$(EXEC) -ind $(DAT) $(ARG)

mcmc: $(EXEC) $(EXEC).psv
	./$(EXEC) -ind $(DAT) -mceval


mceval: $(EXEC)
	cp $(CTL).psv $(EXEC).psv
	./$(EXEC) -ind $(DAT) -mceval

simeval: $(foreach var,$(NUMBERS),./$(OM);\
	./$(EXEC) -ind $(DAT) $(ARG);)
	
        
clean: 
	-rm -f  admodel.* variance eigv.rpt fmin.log $(EXEC) variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor  


cleanOM: 
	-rm -f  admodel.* variance eigv.rpt fmin.log $(OM) variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor *.par 
