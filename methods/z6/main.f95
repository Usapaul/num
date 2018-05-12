program method_setok

use init
use runmethods

implicit none

! pr -- целое число отвечает за параметр точности real

! Массивы X и t заданы в модуле init и называются Xarr и tarr, там же
! они становятся размещенными

! matr -- рабочая матрица, тоже в модуле init она описана и размещена

integer :: i, k

character(1) :: output_files = 'n' ! y/n -- нужно ли записывать в файл резул-ты

!------------------------------------------------
! Процедура init_const() задаст начальные значения количества узлов
! сетки по t и по x. Там же 
call init_const()

! Заполнение первого слоя по данной в модуле init функции phi(x) (=u(x,0))
forall (i=0:n_xgrid) matr(i,0) = phi(Xarr(i))
! И границ:
forall (i=1:m_tgrid)
	matr(0,i) = psi1(tarr(i))
	matr(n_xgrid,i) = psi2(tarr(i))
end forall

! процедура get_ABCD говорит сама за себя, а массивы A, B, C, D -- это
! коэффициенты, получающиеся в сеточном уравнении, где производные были уже
! заменены каким-нибудь выражением (сеточным) для приближения, и оно такое:
! u(i,k+1) = A(i,k)*u(i+1,k) - B(i,k)*u(i,k) + C(i,k)*u(i-1,k) + D(i,k)
! В данном случае это явный метод решения 6-й задачи по вычам:
call get_ABCD()

!------------------------------------------------
! В матрицу matr и будет записано решение задачи вообще:
call explicit_method(matr,A_c,B_c,C_c,D_c)

write(*,*) '    *** EXPLICIT METHOD ***   '
write(*,*) '      number  ', ' Discrepancy'
do i=0,m_tgrid,m_tgrid/20
	write(*,*) i, norma(matr(:,i)-exact_solution(Xarr,tarr(i)))
end do

if (output_files == 'y') then
	open(100,file='res_explicit.dat',status='replace')
	do i=0,m_tgrid
		write(100,*) matr(:,i)
	end do
	close(100)
end if
!------------------------------------------------
! Еще надо в файл записать известное решение функции в табличном виде:
if (output_files == 'y') then
	open(300,file='res_exact.dat',status='replace')
	do i=0,m_tgrid
		write(300,*) exact_solution(Xarr,tarr(i))
	end do
	close(300)
end if
! Вообще для gnuplot3D нужно записывать в файл значения следующим образом:
! open(300,file='res_exact.dat',status='replace')
	! do k=0,m_tgrid
		! do i=0,n_xgrid
			! write(300,*) Xarr(i), tarr(k), exact_solution(Xarr(i),tarr(k))
		! end do
	! end do
! close(300)

!------------------------------------------------
! Теперь отрабатывает неявный метод, подробности см. в 
! комментариях в процедурах

! Для неявного метода выражения для коэффициентов весьма похожи,
! точнее, они же, -- но со знаком минус, а для B_с еще и минус 2
A_c = -A_c
B_c = -B_c - 2
C_c = -C_c
D_c = -D_c 

! Предполагается, что, раз первый слой уже был заполнен еще до работы
! явного метода, а также были заполнены края, согласно краевым условиям, 
! то они и не изменились после отработки явного метода
call implicit_method(matr,A_c,B_c,C_c,D_c)

write(*,*) '  ===========================  '
write(*,*) '    *** IMPLICIT METHOD ***   '
write(*,*) '      number  ', ' Discrepancy'
do i=0,m_tgrid,m_tgrid/20
	write(*,*) i, norma(matr(:,i)-exact_solution(Xarr,tarr(i)))
end do

if (output_files == 'y') then
	open(300,file='res_implicit.dat',status='replace')
	do i=0,m_tgrid
		write(300,*) matr(:,i)
	end do
	close(300)
end if

!------------------------------------------------

! Самокиш и Сулягина просят показать, что при разных значениях n_xgrid
! значения найденного решения u в точках совпадают (с какой-то точностью).
! Для этого нужно проверить, что на сетке с любым разбиением значение
! u(x,t) получается (с какой-то точностью) одинаковым

! Это будет делать процедура find_value. Она !лишняя! в этой программе, нужна
! для проверки лишь Самокишу и Сулягиной, поэтому я ее выписываю отдельно
! и даже не пожалею времени работы программы -- повторю все вычисления.

call find_value()

!================================================

contains

elemental real(pr) function exact_solution(x,t)
	! К моей задаче известно точное решение, вот эта функция для него
	implicit none

	real(pr), intent(in) :: x, t

	!--------------------------------------------
	exact_solution = exp(-0.25_pr * t)*cos(0.5_pr * x) + exp(-t)*(1 - x)*x

end function exact_solution 

real(pr) function norma(A)
	! Функция возвращает норму невязки 
	implicit none

	real(pr), dimension(1:), intent(in) :: A

	!--------------------------------------------
	norma = sqrt(dot_product(A,A))

end function norma


!================================================


subroutine find_value()
	! Лишняя процедура в этой программе, см. комм. выше при вызове процедуры
	! С экрана считываются значения x и t, и находится ближайшие сеточные
	! значения x_i и t_k, и выводится на экран значение найденного 
	! решения в точке (x_i,t_k)
	implicit none

	real(pr) :: x_stdin, t_stdin ! Читаемые с экрана значения x и t
	character(1) :: checking_choice ! (y -- нужна проверка или n -- не нужна)
	integer :: num_x, num_t ! Индексы значений, ближайших к введенным с экрана

	!--------------------------------------------
	write(*,*) ' '
	write(*,*) 'Is the verification necessary? y/n'
	read(*,*) checking_choice
	if (checking_choice /= 'y') return
	write(*,*) 'Enter x and t: '
	read(*,*) x_stdin, t_stdin

	do while (x_stdin > Xarr(n_xgrid) .or. x_stdin < Xarr(0))
		write(*,*) ' Wrong input x. Enter the value between: '
		write(*,*) ' X(0) = ', Xarr(0)
		write(*,*) ' X(n) = ', Xarr(n_xgrid)
		read(*,*) x_stdin
	end do

	do while (t_stdin > tarr(m_tgrid) .or. t_stdin < tarr(0))
		write(*,*) ' Wrong input t. Enter the value between: '
		write(*,*) ' t(0) = ', tarr(0)
		write(*,*) ' t(n) = ', tarr(m_tgrid)
		read(*,*) t_stdin
	end do

	num_x = find_nearest(Xarr,x_stdin)
	num_t = find_nearest(tarr,t_stdin)

	write(*,*) ' '
	write(*,*) ' The nearest point: '
	write(*,*) ' x =', Xarr(num_x), '  t = ', tarr(num_t)
	write(*,*) ' exact_solution(x,t) = ',exact_solution(Xarr(num_x),tarr(num_t))

	!--------------------------------------------
	! Сначала проверка явного метода
	A_c = -A_c
	B_c = -B_c - 2
	C_c = -C_c
	D_c = -D_c 

	call explicit_method(matr,A_c,B_c,C_c,D_c)

	write(*,*) ' '
	write(*,*) ' Verification: explicit method '
	write(*,*) '              u(x,t) = ', matr(num_x,num_t)
	write(*,*) '         Discrepancy = ', & 
				& matr(num_x,num_t) - exact_solution(Xarr(num_x),tarr(num_t))
	!--------------------------------------------
	A_c = -A_c
	B_c = -B_c - 2
	C_c = -C_c
	D_c = -D_c 

	call implicit_method(matr,A_c,B_c,C_c,D_c)

	write(*,*) ' '
	write(*,*) ' Verification: implicit method '
	write(*,*) '              u(x,t) = ', matr(num_x,num_t)
	write(*,*) '         Discrepancy = ', & 
				& matr(num_x,num_t) - exact_solution(Xarr(num_x),tarr(num_t))

end subroutine find_value

integer function find_nearest(Array,x)
	! Возвращает индекс ближайшего к x элемента вещественного массива Array,
	! который имеет фактические индексы от нуля
	implicit none

	real(pr), dimension(0:), intent(in) :: Array
	real(pr), intent(in) :: x
	integer, dimension(1) :: index

	!--------------------------------------------
	index = minloc(abs(Array-x)) - 1 ! -1, т.к. Array имеет нулевой эл-т
	find_nearest = index(1)

end function find_nearest

end program method_setok