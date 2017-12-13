graph: results yyy.dat
	echo plot "'ggg.dat'" using 1:2 w l,\
	"'yyy.dat'" using 1:2 w p > plot.dat
	echo pause'(-1)' >> plot.dat
	gnuplot plot.dat
	rm plot.dat

results: start
	rm -f *.mod
	./start

start:  funs.f95 general.f95
	gfortran -fbounds-check $^ -o start

