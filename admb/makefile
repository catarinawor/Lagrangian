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

DIROM=OM/gtg/
DIREXEC=mov_est/gtg/
EXEC=lagrangian_est_gtg
DAT=mov_est/gtg/lagrangian_est_gtg.dat
OM=lagrangian_OM_gtg
SIMDAT=OM/simple/lagrangian_OM_gtg.dat
SEED = seed.txt
LAST := 2
SADIR= /Users/catarinawor/Documents/iSCAM/examples/hakelag/DATA


NUMBERS := $(shell seq 1 ${LAST})

RDATA ='source(file.path("..","simeval/SimResult_gtg_5areas_tau04","simRun.r"))'

RDATASA ='source(file.path("../","Lagrangian/","RunSA.r"))' 

all: 
	omi exi

omi: 
	- $(MAKE) --directory $(DIROM) 

exi: 
	- $(MAKE) --directory $(DIREXEC) 



sarun:
	make clean --directory=$(SADIR)
	-$(MAKE) --directory=$(SADIR)


readRdat: 
	@echo $(RDATA) | R --vanilla --slave

readSA: 
	@echo $(RDATASA) | R --vanilla --slave



doall: omi exi readRdat  
	

dosa: rom sarun readSA

simeval: $(foreach var,$(NUMBERS), $(MAKE) doall ;)

simsa: $(foreach var,$(NUMBERS), $(MAKE) dosa ;)


qwert: 
	for F in $(NUMBERS) ; do \
    echo $$doall ; done


FORCE: 

cleanall: 
	-rm -f  admodel.* variance eigv.rpt fmin.log $(EXEC) variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor *.obj *.p02 *.b02 *.par *.r02 *.rep 

        
clean: 
	-rm -f  admodel.* variance eigv.rpt fmin.log $(EXEC) variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor *.obj *.p02 *.b02 *.par *.r02


cleanOM: 
	-rm -f  admodel.* variance eigv.rpt fmin.log $(OM) variance *.b01 *.p01 *.r01 *.eva *.bar *.log *.htp *.cor *.par *.obj
