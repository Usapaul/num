results: start
	rm -f *.mod
#!	echo '****       make rk - Метод Рунге-Кутты         ****'
#!	echo '**** make ae - Экстраполяционный метод Адамса  ****'
#!	echo '**** make ai - Интерполяционный метод Адамса   ****'


start:  mprecision.f90 lead.f90 newtonmethod.f90 ffunction.f90 ABfunctions.f90  rungekut.f90 adams.f90 general.f90
	gfortran -fbounds-check $^  -o start

rk: start
	./start rk

ae: start
	./start ae

ai: start
	./start ai
