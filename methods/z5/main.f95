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

! b_m_vector -- это вектор правой части в получаемой СЛАУ 
real(pr), dimension(:), allocatable :: b_m_vector 

integer :: i

integer :: i_for_KK ! Нужно для вывода результатов
real(pr), dimension(:), allocatable :: X_results, Discrepancy ! тоже для вывода

!------------------------------------------------
call init_consts() ! Начнем с инициализации постоянных из условия
! В то время как K(x,t) и f(x) нужно явно описывать в модуле init

! n_grid -- количество узлов коллокации, определяется в модуле init
allocate(matr(n_grid,n_grid)) ! Матрица из коэфф., используемая в решении
allocate(polynom_c(n_grid)) ! массив из искомых коэффициентов в разложении
allocate(b_m_vector(n_grid)) ! вектор правой части в получаемой СЛАУ
allocate(X(n_grid)) ! массив из узлов сетки

! Выбор полинома здесь, см. модуль init. 1 -- Лежандра, 2 -- Чебышева
what_the_polynomial = 2

! Процедура polynom просто создаст массив X из узлов сетки
call polynom()

! Для решения задачи нужно построить матрицу
! Процедура coeff() сохранит в matr значения этих коэффициентов
! для дальшейшей работы с ней (решения системы уравнений)
call coeff(matr)

! Для нахождения коэффициентов polynom_c нужно решить систему линейных
! уравнений, матрица которой -- matr, а столбец в правой части системы
! -- это вектор, состоящий из значений интегралов от f(x)*phi(m,x), и
! эта подынтегральная функция описана в модуле coeffs

do i=1,n_grid
	m_from_main = i ! функции f_x_phi в coeffs нужно значение текущего m
	b_m_vector(i) = integral(f_x_phi,a_left,b_right)
end do
polynom_c = solve(matr,b_m_vector)

!================================================

! Все, коэффициенты найдены, полиномы выбраны, замечательно!
! Теперь, чтобы убедиться в том, что программа работает правильно,
! нужно вывести таблицу значений функции и сравнить с тем, что
! как бы должно получиться (взяв модельную задачу с известным решением)

write(*,*) ' '
write(*,*) 'polynom coefficients:'
call print_vector(polynom_c)

! В Discrepancy будет храниться вычисленная невязка. То есть то значение,
! которое получится, если в исходное интегральное уравнение вставить
! найденную мной функцию -- квадратурную, сумму из полиномов. Потом
! эту невязку я буду выводить вместе с результатами, чтобы было видно,
! что программа работает. А то сухие числа ничего и не дают ведь

! А X_results -- значения абсциссы, в которых я значение функции
! вычисляю и вывожу на экран (или в файл)

allocate(X_results(0:10),Discrepancy(0:10))
X_results = (/(0.0_pr + 0.1_pr*i,i=0,10)/)

do i=0,size(X_results)-1
	! i_for_KK говорит отчасти само за себя. Я там ниже определил функцию
	! KK_phi(t), которая представляет из себя произведение K*phi, но с 
	! нужными текущими значениями i для X(i). Вот эта переменная там нужна: 
	i_for_KK = i

	! В одну строку такое не поместилось, и в две тоже, поэтому их... три!
	! В этих трех строчках я действительно честно подставляю найденное
	! решение в исходное интегральное уравнение, да в итоге нахожу невязку
	! для каждой точки из X_results
	Discrepancy(i) = result_function(X_results(i)) * a_const
	Discrepancy(i) = Discrepancy(i) + integral(KK_phi,a_left,b_right)
	Discrepancy(i) = Discrepancy(i) - f(X_results(i))
end do

call output_function(result_function,X_results,Discrepancy)


contains


subroutine print_vector(X)
	! Тупо печать вектора в удобном виде
	real(pr), dimension(1:), intent(in) :: X
	integer :: i

	do i=1,size(X)
		write(*,'(f15.12)') X(i)
	end do
end subroutine print_vector

pure real(pr) function KK_phi(t)
	! Нужна только лишь для вывода невязки
	implicit none

	real(pr), intent(in) :: t

	!--------------------------------------------
	KK_phi = K(X_results(i_for_KK),t) * result_function(t)

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

subroutine output_function(f,X,Discrepancy,idfile)
	! Процедура выполняет запись табличных значений функции f в файл,
	! или на экран, и ничего больше не делает
	implicit none

	real(pr), dimension(1:), intent(in) :: X, Discrepancy
	integer, optional :: idfile
	integer :: id = 987 ! рандомное число
	integer :: i

	interface
		pure real(pr) function f(t)
			use init, only: pr
			real(pr), intent(in) :: t
		end function f
	end interface

	!--------------------------------------------
	if (present(idfile)) then
		id = idfile
		open(id,file='fun_table.dat',status='replace')
		write(id,*) X(i), f(X(i))
		close(id,status='keep')	
	end if

	write(*,*) ' '
	write(*,*) '        x        ', '     u(x)      ', '   Discrepancy      '
	do i=1,size(X)
		write(*,'(2(2x,f13.9),5x,e13.6)') X(i), f(X(i)), Discrepancy(i)
	end do

end subroutine output_function

end program collocation