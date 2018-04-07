program sobstvenno

use init
use runmethods

implicit none

! Параметр pr, задающий точность, описан в модуле
! Порядок матрицы также задается в модуле, это глобальная переменная s
! сама матрица -- динамический двумерный массив matr

integer :: num_iter ! количество итераций
real(pr) :: l_max ! наибольшее по модулю собственное число
real(pr) :: l_min ! а это наименьшее
real(pr), dimension(:), allocatable :: X0 ! нач. приближение

!------------------------------------------------

! matrix-initialization (согласно условию моей задачи):
call matrixinit()

allocate(X0(s))
X0 = 1.0e-1_pr ! Ничего больше, чем создание нач. приближения

! Решение задачи методом последовательных итераций:
call itermethod(matr,X0,l_max,num_iter)
write(*,*) 'l_max:', l_max
write(*,*) 'number of iterations:', num_iter

! Теперь методом итераций находим и наименьшее с.ч.
! при помощи сдвига:
call itermethod((matr-l_max*Ematr),X0,l_min,num_iter)
l_min = l_min + l_max ! потому что метод сдвига он такой!
write(*,*) 'l_min:', l_min 
write(*,*) 'number of iterations:', num_iter

! Теперь преобразуем матрицу к трехдиаг. виду по методу Гивенса:
diag3 = givens(matr)


contains

subroutine stdw(matr,nice)
	! Процедура просто в нормальном, понятном и управляемом
	! виде выводит на экран матрицу
	! Если nice=1, то кол-во знаков после запятой сокращено
	! Если nice=3, то будет вывод трехдиаг. матрицы, и
	! кол-во знаков после запятой у ненулевых эл-тов сокращено
	implicit none

	real(pr), dimension(1:,1:), intent(in) :: matr
	integer, optional :: nice
	integer :: i, j, n

	n = size(matr,1)
	!--------------------------------------------
	write(*,*) ' '
	if (present(nice)) then:
		if (nice == 3) then:
			do i=1,n
				write(*,1111) (/(0,j=1,i-2)/),mart(i-1:i+1,i),(/(0,j=i+2,n)/)
			end do
			1111 format fordiag3(*(i4,3x),*(x,f7.3),*(i4,3x))
		elseif (nice == 1) then:
			do i=1,n
				write(*,'(*(x,f7.3))') matr(:,i)
			end do			
		else
			write(*,*) 'indefined NICE parameter'
	else
		do i=1,n
			write(*,*) matr(:,i)
		end do		
	end if

end subroutine stdw

end program sobstvenno
