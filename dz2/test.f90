implicit none

integer :: i, n=10
real :: a=-5.0, b=5.0
real :: h

h=(b-a)/2

open(100,file='sinus.dat')

    write(100,'("# ",i8)') n
    write(100,*) a, b
    do i=0,n
	write(100,*) sin(a+i*h)
    enddo

close(100)

end