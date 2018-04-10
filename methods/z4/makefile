MAIN = main

COEFFS = coeffs
RUNMETHODS = runmethods
SLAE = slae

INTEGRALMOD = integralmod
POLYNOMIALS = polynomials

INIT = init
#------------------------------------------------
MAIN_D = $(INIT) $(COEFFS) $(RUNMETHODS) $(SLAE)
COEFFS_D = $(INIT) $(POLYNOMIALS) $(INTEGRALMOD)
RUNMETHODS_D = $(INIT) $(POLYNOMIALS)
SLAE_D = $(INIT)
INTEGRALMOD_D = $(INIT)
POLYNOMIALS_D = $(INIT)
INIT_D = 
#------------------------------------------------
OBJECTFILES = $(patsubst %.f95,%.o,$(wildcard *.f95))
EXENAME = start
COMPOBJ = gfortran -c 
COMPEXE = gfortran -o $(EXENAME) -fcheck=all 

LIB = -llapack
#------------------------------------------------
aim = $($(1)).mod $($(1)).o 
depend = $($(1)).f95 $(addsuffix .mod,$($(1)_D))
whole_rule = $(call aim,$(1)) : $(call depend,$(1)) ; \
				$(COMPOBJ) $$<; touch $($(1)).mod
#------------------------------------------------

$(EXENAME) : $(OBJECTFILES)
	$(COMPEXE) $^

$(call whole_rule,MAIN)
$(call whole_rule,RUNMETHODS)
$(call whole_rule,COEFFS)
$(call whole_rule,SLAE)
$(call whole_rule,INTEGRALMOD)
$(call whole_rule,POLYNOMIALS)
$(call whole_rule,INIT)

.PHONY : clean run

run : $(EXENAME)
	./$(EXENAME)

clean : 
	rm -f *.mod
	rm -f *.o
