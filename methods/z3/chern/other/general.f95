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

call stdw(diag3,nice=3)

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
	if (present(nice)) then
		if (nice == 3) then
			write(*,'(2(x,f7.3),*(4x,i4))') matr(1:2,1), (/(0,j=3,n)/)
			write(*,'(3(x,f7.3),*(4x,i4))') matr(1:3,2), (/(0,j=4,n)/)
			do i=3,n-1
				write(*,1111) (/(0,j=1,i-2)/),matr(i,i-1:i+1),(/(0,j=i+2,n)/)
			end do
			write(*,1111) (/(0,j=1,n-2)/), matr(n-1:n,n)
			1111 format(*(4x,i4),*(x,f7.3),*(i4,3x))
		elseif (nice == 1) then
			do i=1,n
				write(*,'(*(x,f7.3))') matr(:,i)
			end do			
		else
			write(*,*) 'indefined NICE parameter'
		end if
	else
		do i=1,n
			write(*,*) matr(:,i)
		end do		
	end if
end subroutine stdw

end program sobstvenno
