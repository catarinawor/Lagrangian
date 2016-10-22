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


RDATA ='source(file.path("read_mse.r"))'
RDATATWO ='source(file.path("read_iscam.r"))'
RDATAOM ='source(file.path("read_OMrep.r"))'
NUMBERS := $(shell seq 1 ${LAST})
LAST := 10
OM=../../admb/OM/simple/lagrangian_OM


readRdat: 
	@echo $(RDATA) | R --vanilla --slave


readRdattwo: 
	@echo $(RDATATWO) | R --vanilla --slave

readROM: 
	@echo $(RDATAOM) | R --vanilla --slave

rom: $(OM) 
	-./$(OM)


simeval: $(foreach var,$(NUMBERS), $(MAKE) rom ;)
