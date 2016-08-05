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

EXEC=lagrangian_est
DAT=lagrangian_est.dat
OM=lagrangian_OM
SIMDAT=lagrangian_OM.dat
SEED = seed.txt
LAST := 100
SADIR= /Users/catarinawor/Documents/iSCAM/examples/hakelag/DATA


NUMBERS := $(shell seq 1 ${LAST})

RDATA ='source(file.path("../","Lagrangian/SimResult_5areas_tau04","simRun.r"))'

RDATASA ='source(file.path("../","Lagrangian/","RunSA.r"))' 

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


sarun:
	make clean --directory=$(SADIR)
	-$(MAKE) --directory=$(SADIR)


readRdat: 
	@echo $(RDATA) | R --vanilla --slave

readSA: 
	@echo $(RDATASA) | R --vanilla --slave



doall: rom run readRdat 

dosa: rom sarun readSA

simeval: $(foreach var,$(NUMBERS), $(MAKE) doall ;)

simsa: $(foreach var,$(NUMBERS), $(MAKE) dosa ;)


qwert: 
	for F in $(NUMBERS) ; do \
    echo $$doall ; done


FORCE: 

        
clean: 
	-rm -f  admodel.* variance eigv.rpt fmin.log $(EXEC) variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor  


cleanOM: 
	-rm -f  admodel.* variance eigv.rpt fmin.log $(OM) variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor *.par *.obj
