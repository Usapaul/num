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
real(pr), dimension(:), allocatable :: L_all ! Все собственные значения

!------------------------------------------------

! matrix-initialization (согласно условию моей задачи):
call matrixinit()

allocate(X0(s),L_all(s))
X0 = 1.0_pr ! Ничего больше, чем создание нач. приближения

! Решение задачи методом последовательных итераций:
call itermethod(matr,X0,l_max,num_iter)
write(*,*) 'l_max:', l_max
write(*,*) 'number of iterations:', num_iter
just_l_max = abs(l_max) ! переменная just_l_max будет нужна для метода Гивенса

! Теперь методом итераций находим и наименьшее с.ч.
! при помощи сдвига:
call itermethod((matr-l_max*Ematr),X0,l_min,num_iter)
l_min = l_min + l_max ! потому что метод сдвига он такой!
write(*,*) 'l_min:', l_min 
write(*,*) 'number of iterations:', num_iter

! Теперь преобразуем матрицу к трехдиаг. виду по методу Гивенса:
diag3 = givens(matr)

!call stdw(diag3,nice=4)

! После того как получена подобная исходной трехдиагональная матрица
! функция minor выражает полином (в виде массива из коэффициентов),
! корни которого -- и есть искомые собственные числа.
! А функция find_roots уже говорит сама за себя, см. ниже вызов функции:
L_all = find_roots(minor(diag3))

L_all = sorted(L_all)
write(*,*) ' '
write(*,*) 'L_all:'
do j=1,size(L_all)
	write(*,'(x,i5,x,f15.10)') j, L_all(j)
end do

contains

pure function sorted(array)
	implicit none

	real(pr), dimension(1:), intent(in) :: array
	real(pr), dimension(1:size(array)) :: sorted
	integer :: i

	!--------------------------------------------
	sorted(1) = maxval(array)
	do i=2,size(array)
		sorted(i) = maxval(array, mask=(array < sorted(i-1)))
	end do
end function sorted

subroutine stdw(matr,nice)
	! Процедура просто в нормальном, понятном и управляемом
	! виде выводит на экран матрицу
	! Если nice=x, то кол-во знаков после запятой сокращено до x
	implicit none

	real(pr), dimension(1:,1:), intent(in) :: matr
	integer, optional :: nice
	integer :: i, j, n
	character(1) :: fornice

	n = size(matr,1)
	!--------------------------------------------
	write(*,*) ' '
	if (present(nice)) then
		write(fornice,'(i1)') nice
		do i=1,n
			write(*,'(*(x,f9.'//fornice//'))') matr(:,i)
		end do			
	else
		do i=1,n
			write(*,*) matr(:,i)
		end do		
	end if
end subroutine stdw

end program sobstvenno
