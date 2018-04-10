program collocation

use init
use runmethods
use coeffs
use slae

implicit none

real(pr), dimension(:,:), allocatable :: matr
! polynom_c -- искомые коэффициенты при разложении решения в виде
! линейной комбинации ортогональных многочленов
real(pr), dimension(:), allocatable :: polynom_c
integer :: i




integer :: iii, kkk
real(pr), dimension(10) :: xes 




!------------------------------------------------
call init_consts() ! Начнем с инициализации постоянных из условия
! В то время как K(x,t) и f(x) нужно явно описывать в модуле init

! n_grid -- количество узлов коллокации, определяется в модуле init
allocate(matr(n_grid,n_grid)) ! Матрица из коэфф., используемая в решении
allocate(polynom_c(n_grid)) ! массив из искомых коэффициентов в разложении
allocate(X(n_grid)) ! массив из узлов сетки

! Выбор полинома здесь, см. модуль init. 1 -- Лежандра, 2 -- Чебышева
what_the_polynomial = 1

! Процедура polynom просто создаст массив X из узлов сетки
call polynom()

! Для решения задачи нужно построить матрицу
! Процедура coeff() сохранит в matr значения этих коэффициентов
! для дальшейшей работы с ней (решения системы уравнений)
call coeff(matr)

! Для нахождения коэффициентов polynom_c нужно решить систему линейных
! уравнений, матрица которой -- matr, а столбец в правой части системы
! -- это вектор, состоящий из значений f(X(i)).



open(120,file='matrix.dat',status='replace')
do i=1,n_grid
	write(120,*) matr(i,:)
end do
close(120,status='keep')

open(130,file='fff.dat')
do i=1,n_grid
	write(130,*) f(X(i))
end do
close(130)



polynom_c = solve(matr,f(X))



!================================================

! Все, коэффициенты найдены, полиномы выбраны, замечательно!
! Теперь, чтобы убедиться в том, что программа работает правильно,
! нужно вывести таблицу значений функции и сравнить с тем, что
! как бы должно получиться (взяв модельную задачу с известным решением)

do i=1,n_grid
	write(*,*) polynom_c(i)
end do
call output_function(result_function,(/(-1 + 0.2_pr*i,i=0,10)/))



write(*,*) ' '
xes = (/(-1 + 0.2_pr*i,i=1,size(xes))/)
do i=1,size(xes)
	write(*,*) xes(i), result_function(xes(i))
end do
write(*,*) ' '
do i=1,size(xes)
	write(*,*) 'i ===', i, result_function(xes(i)) - f(xes(i))
	iii = i
	write(*,*) 'integral:', integral(KK_phi,a_left,b_right)
end do


!call la_syev(matr,polynom_c)
!do i=1,n_grid
!	write(*,*) polynom_c(i)
!end do



contains







pure real(pr) function KK_phi(t)
	implicit none

	real(pr), intent(in) :: t

	!--------------------------------------------
	KK_phi = K(xes(iii),t) * result_function(t)

end function KK_phi

















pure real(pr) function result_function(x)
	! Просто опишу функцию решения в том виде, в котором она
	! и есть: u(x) = сумма из произведений полиномов с коэфф. polynom_c
	implicit none

	real(pr), intent(in) :: x
	integer :: k

	!--------------------------------------------
	result_function = dot_product(polynom_c,(/(phi(k,x),k=1,size(polynom_c))/))
end function result_function

subroutine output_function(f,X,id_file)
	! Процедура выполняет запись табличных значений функции f в файл,
	! и ничего больше не делает. idfile можно задать как std_out
	implicit none

	real(pr), dimension(1:), intent(in) :: X
	integer, optional :: id_file
	integer :: aaa
	integer :: i

	interface
		pure real(pr) function f(t)
			use init, only: pr
			real(pr), intent(in) :: t
		end function f
	end interface

	!--------------------------------------------
	! if ( .not. present(id_file)) then
		! id_file = 1987 ! рандомное число
		! write(*,*) 'HERE OK'
	! end if
	aaa = 895
	open(aaa,file='fun_table.dat',status='replace')
	do i=1,size(X)
		write(aaa,*) X(i), f(X(i))
	end do
	close(aaa,status='keep')
end subroutine output_function

end program collocation