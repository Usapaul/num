program sobstvenno

use init

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
call givens(matr,diag3)

write(*,*) ' '
do s=1,size(matr,1)
	write(*,*) diag3(:,s)
end do

contains

subroutine givens(matr,diag3)
	implicit none

	real(pr), dimension(1:,1:), intent(in) :: matr
	real(pr), dimension(1:size(matr,1),1:size(matr,1)), intent(out) :: diag3
	real(pr) :: theta, phi, c, s
	integer :: i, j, n

	!--------------------------------------------
	n = size(matr,1) 

	! Реализация метода вращений Гивенса,
	! приведение к трехдиагональному виду исходной матрицы:
	do i=1,n-1
		do j=i+1,n
			theta = 2 * matr(j,i)/(matr(i,i) - matr(j,j))
			phi = 0.5_pr * atan(theta)
			c = cos(phi)
			s = sin(phi)
			write(*,*) c, s, theta

			diag3(i,i) = c*c*matr(i,i) + 2*c*s*matr(j,i) + s*s*matr(j,j)	
			diag3(j,j) = s*s*matr(i,i) - 2*c*s*matr(j,i) + c*c*matr(j,j)
			diag3(j,i) = s*c*(matr(j,j) - matr(i,i)) + matr(j,i)*(c**2-s**2)
			diag3(i,j) = diag3(j,i)
		end do
	end do
end subroutine givens

end program sobstvenno
