results: start data.dat
	./start
	rm -f *.mod
	rm -f data.dat

start:  lead.f90 ffunction.f90 newtonmethod.f90 inandout.f90
	gfortran -fbounds-check $^  -o start

data.dat:
	touch data.dat
	echo '# 3' >> data.dat
	echo '1.2' >> data.dat
	echo '1.3' >> data.dat
	echo '1.0' >> data.dat