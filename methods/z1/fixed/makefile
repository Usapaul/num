graph: results yyy.dat
	echo plot "'yyy.dat'" > plot.dat
	echo pause'(-1)' >> plot.dat
	gnuplot plot.dat
	rm plot.dat

results: start
	rm -f *.mod
	./start

start:  funs.f95 general.f95
	gfortran -fbounds-check $^ -o start
