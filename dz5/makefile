plot: results
	rm -f *.mod
	echo plot "'result.dat'" using 1:2 w l,\
	"'data.dat'" using 1:2 w p > plot.dat
	echo pause'(-1)' >> plot.dat
	gnuplot plot.dat
	rm plot.dat

results: start testfun
	./start

start: diag3comp.f90 diag3dim1.f90 progon5.f90 workmatrix.f90 makeABQ.f90 spline.f90 solve.f90 inandout.f90
	gfortran -fbounds-check $^  -o start

testfun: test
	./test

test: test.f90
	gfortran test.f90 -o test