## Makefile for cleaning up
## I don't think it works on pc! 
## Targets: 
##		all:   -copy executable and run the model with DAT & ARG
##		run:   -copy executable and force a run
##		mcmc:  -copy executable and run in mcmc mode and mceval
##		retro: -copy executable and run  retrospective analysis



EXEC=lagrangian_est
DAT=lagrangian_est.dat
OM=lagrangian_OM_gtg
SIMDAT=lagrangian_OM_gtg.dat
SEED = seed.txt
LAST := 100

NUMBERS := $(shell seq 1 ${LAST})

RDATA ='source(file.path("../","Lagrangian/SimResult_gtg_5areas_tau04","simRun.r"))'

 
all: 
	compiOM rom compi run

compiOM: $(OM)
	admb $(OM)

compi: $(EXEC)
	admb $(EXEC)

rom: $(OM) $(SIMDAT) $(ARG) $(SEED)
	-./$(OM) -ind $(SIMDAT) $(ARG)


run: $(EXEC) $(DAT) $(ARG)
	-./$(EXEC) -ind $(DAT) $(ARG)

readRdat: 
	@echo $(RDATA) | R --vanilla --slave


doall: rom run readRdat 



simeval: $(foreach var,$(NUMBERS), $(MAKE) doall ;)

qwert: 
	for F in $(NUMBERS) ; do \
    echo $$doall ; done


FORCE: 

        
clean: 
	-rm -f  admodel.* variance eigv.rpt fmin.log $(EXEC) variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor  


cleanOM: 
	-rm -f  admodel.* variance eigv.rpt fmin.log $(OM) variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor *.par *.obj
