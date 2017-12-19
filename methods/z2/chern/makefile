run : start
	./start
	rm -f *.mod

start : functions.f95 general.f95
	gfortran $^ -o $@


plot : results.dat
	echo plot "'results.dat'" using 1:2 w l,\
	"'results.dat'" using 1:3 w l > plot.txt
	echo pause '(-1)' >> plot.txt
	gnuplot plot.txt
	rm plot.txt