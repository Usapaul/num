program test
implicit none
integer :: i, n, nfun
real :: a, b, h, x

!---------------------------------------
write(*,*) 'Для тестирования программы следует выбрать известную функцию'
write(*,*) 'Например, sin(x), x**2, sqrt(x), затухающее колебание y=A*e**(-b*x)*sin(w*x)'
write(*,*) 'Если интересует одна из предложенных - введите ее номер: 1, 2, 3 или 4'
write(*,*) 'Или введите "0", если в директории уже есть файл data.dat с данными'
!---------------------------------------
read(*,*) nfun
do while ( nfun > 4 .or. nfun < 0 )
	write(*,*) 'Введено число, отличное от {0,1,2,3,4}. Введите заново'
	read(*,*) nfun
enddo
!---------------------------------------
if (nfun==0) go to 88
write(*,*) 'Введите  N - число интервалов аппроксимации'
read(*,*) n
a=0.0; b=10.0; h=(b-a)/n; x=a
open(100,file='data.dat')
    write(100,'("# ",i8)') n
    select case(nfun)
    case(1); do i=0,n; write(100,*) x+h*i, sin(a+i*h), 30.0; enddo
    case(2); do i=0,n; write(100,*) x+h*i-b/2, (a-b/2+i*h)**2, 30.0; enddo
    case(3); do i=0,n; write(100,*) x+h*i, sqrt(a+i*h), 30.0; enddo
    case(4); do i=0,n; write(100,*) x+h*i, 2.0*exp(-0.2*(a+i*h))*sin(1.5*(a+i*h)), 30.0; enddo
    end select
close(100)

88 end program test