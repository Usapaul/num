MAIN = main

METHOD1 = runmethods
METHOD2 = coeffs

METHODS = 
FUNCTIONS = integralmod polynomials slae
INITIAL = init

OBJECTFILES = 

EXENAME = start
COMPOBJ = gfortran -c 
COMPEXE = gfortran -o $(EXENAME) -fcheck=all
#-fbounds-check
LIB = -llapack

run : $(EXENAME)
	./$(EXENAME)
	rm -f *.mod

$(EXENAME) : $(INITIAL) $(FUNCTIONS) $(METHODS) $(MAIN)
	$(COMP) $^ $(LIB)

$(EXENAME) : $(INITIAL) $(FUNCTIONS) $(METHODS) $(MAIN)
	$(COMPEXE) $(LIB) 
