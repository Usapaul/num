results: start data.dat
	./start 1
	rm -f *.mod

start:  funW.f90 preobrazovanie.f90 inandout.f90
	gfortran -fbounds-check $^  -o start
