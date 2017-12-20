program hard_systems

use functions

implicit none

! Параметр pr, задающий точность, описан в модуле
integer  :: Ntot ! Ntot(al) -- количество промежутков разбиения
real(pr), dimension(:), allocatable :: X, Y, Z
integer :: i
character(1) :: test

!------------------------------------------------
write(*,*) 'Input from keyboard? y/n'
read(*,*) test
if (test == 'y') then
	write(*,*) 'Input N...'
	read(*,*) Ntot
	write(*,*) 'Input A...'
	read(*,*) A
	! Процедура check_data_A_N не делает ничего больше, чем проверку
	! чисел A и N на >/<= 0
	call check_data_A_N()
else
	Ntot = 150
	A = -50
end if

allocate(X(0:Ntot*2), Y(0:Ntot*2), Z(0:Ntot*2))
X(0) = x0
X(1:) = 0
Y(0) = y0
Y(1:) = 0
Z(0) = z0
Z(1:) = 0

! xend является правой границей рассматриваемого промежутка
! и по условию задачи он задается как модуль 1/(A*d)
xend = abs(1/(d*A))
h = (xend - x0) / Ntot
forall (i=1:Ntot) X(i) = x0 + h*i

! Коэффициенты Czn1, Czn, Cyn, Csqrt описаны в модуле
call initC()
call solve(Y=Y,Z=Z,n1=1,n2=Ntot,h=h)

! Теперь вторая часть за точкой xend и до 1:
x0 = xend
xend = 1
h = (xend - x0) / Ntot
forall (i=Ntot+1:Ntot*2) X(i) = x0 + h*(i-Ntot)
call initC()
call solve(Y=Y,Z=Z,n1=Ntot+1,n2=2*Ntot,h=h)

open(100,file='results.dat',status='replace')
	do i = 0, 2*Ntot
		write(100,*) X(i), Y(i), Z(i)
	end do
close(100,status='keep')

!================================================

contains

subroutine check_data_A_N()
	do while (A>=0)
		write(*,*) 'Incorrect A. It should be < 0. Enter A again '
		read(*,*) A
	end do
	do while (Ntot<=0)
		write(*,*) 'Incorrect N. It should be > 0. Enter N again '
		read(*,*) Ntot
	end do	
end subroutine check_data_A_N

end program hard_systems