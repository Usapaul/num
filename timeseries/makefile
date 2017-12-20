#results: start data.dat
#	./start 1
#	rm -f *.mod


run : start
	./start

start :  functions.f95 funW.f95 preobrazovanie.f95 dpfurye.f95 functions2.f95 general.f95
	gfortran -fbounds-check $^  -o start
	rm -f *.mod
