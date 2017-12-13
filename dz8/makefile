results: start
	./start

start:  ffunction.f90 integral.f90
	gfortran -fbounds-check $^  -o start
	rm -f *.mod

tablichki:  bernul.f90 lead.f90 problems.f90 quad.f90
	gfortran -fbounds-check $^  -o tablichki
	rm -f *.mod
