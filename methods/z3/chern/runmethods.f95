module runmethods

use init

implicit none

contains 

subroutine itermethod(matr,X0,l_1,count)
	! Реализация степенного метода (итер.) нахождения max по модулю с.ч.
	implicit none

	real(pr), dimension(1:), intent(in) :: X0
	real(pr), dimension(1:size(X0),1:size(X0)), intent(in) :: matr
	real(pr), parameter :: eps = 1e-3 * 0.1**pr ! точнее, чем ">" pr
	real(pr), dimension(1:size(X0)) :: Xnew, X 
	real(pr), intent(out) :: l_1 ! наибольшее собственное число 
	integer, parameter :: nummax=3e5 ! Максимально допустимое число итераций
	integer, intent(out) :: count ! счетчик числа итераций
	real(pr) :: nevyazka

	!--------------------------------------------
	! Начальное приближение -- это X0, с него и стартуем
	! X -- текущий вектор, а Xnew уже говорит сам за себя
	X = X0
	Xnew = X0 
	Xnew(1) = Xnew(1) + 1 ! чтобы в do while не было совпадения сразу
	count=0
	do while (norma(Xnew/Xnew(1) - X) > eps .and. count < nummax)
		X = Xnew / Xnew(1)
		Xnew = matmul(transpose(matr),X)
		count = count + 1
		if (mod(count,500)==0) write(*,*) count
	end do
	l_1 = Xnew(1)

	X = Xnew / Xnew(1)

	write(*,*) ' '

	if (count == nummax) then
		write(*,*) 'Number of iterations reached a value of NUMMAX!'
	end if

	nevyazka = norma(matmul(transpose(matr),X) - l_1*X)
	if (nevyazka < eps*10*s) then
		write(*,*) '||Au_1 - l_1*u_1|| =',  nevyazka, ' ALL RIGHT'
	else
		write(*,*) '||Au_1 - l_1*u_1|| =',  nevyazka, ' > eps, cry...'
	end if
end subroutine itermethod

pure real(pr) function norma(X)
	real(pr), dimension(1:), intent(in) :: X
	norma = sqrt(sum(X**2))
end function norma

!================================================

pure function givens(matr) result(diag3)
	implicit none

	real(pr), dimension(1:,1:), intent(in) :: matr
	real(pr), dimension(1:size(matr,1),1:size(matr,1)):: diag3
	! Следующие матрицы будут постоянно изменяемые, по сути, вспомогательные:
	real(pr), dimension(1:size(matr,1),1:size(matr,1)) :: matrix, rot_pq
	real(pr) :: theta, phi, c, s
	integer :: p, q, n

	!--------------------------------------------
	n = size(matr,1) 
	! с этого момента индексы строк и столбцов будут записаны в
	! привычном естественном порядке:
	matrix = transpose(matr)

	! Реализация метода вращений Гивенса,
	! приведение к трехдиагональному виду исходной матрицы:
	do p=2,n-1
		do q=p+1,n
			c = matrix(p-1,p) / sqrt(matrix(p-1,p)**2 + matrix(p-1,q)**2)
			s = matrix(p-1,q) / sqrt(matrix(p-1,p)**2 + matrix(p-1,q)**2)
			rot_pq = rotmatr(n,p,q,c,s)
			matrix = matmul(transpose(rot_pq),matrix)
			matrix = matmul(matrix,rot_pq)
		end do
	end do

	! Наконец, получили трехдиаг. матрицу в результате:
	diag3 = transpose(matrix)
end function givens

pure function rotmatr(n,p,q,c,s) result(rot)
	implicit none

	integer, intent(in) :: p, q, n ! номера строки, столбца; размер матрицы
	real(pr), intent(in) :: c, s
	real(pr), dimension(1:n,1:n) :: rot
	integer :: i

	!--------------------------------------------
	rot = 0
	forall (i=1:n) rot(i,i) = 1
	rot(p,p) = c
	rot(q,q) = c
	rot(p,q) = -s
	rot(q,p) = s
end function rotmatr

!================================================

pure recursive function minor(matr) result(D)
	implicit none

	! Функция minor работает с матрицей, индексы
	! строк и столбцов которой записываются
	! в естественном порядке (строка,столбец)
	real(pr), dimension(1:,1:), intent(in) :: matr
	! D -- это массив из коэфф. при соответствующих степенях (от 0 до n)
	! у полинома D_n(l), выражение которого и есть минор
	! D_n1 и D_n2 -- это те же полиномы степени n-1 и n-2
	real(pr), dimension(0:size(matr,1)) :: D 
	real(pr), dimension(0:size(matr,1)-1) :: D_n1
	real(pr), dimension(0:size(matr,1)-2) :: D_n2
	integer :: n

	n = size(matr,1)
	!--------------------------------------------
	if (n > 1) then
		D_n1 = minor(matr(1:n-1,1:n-1))
		D_n2 = minor(matr(1:n-2,1:n-2))

		! Пришлось потрудиться, чтобы додуматься до этих формул и
		! сложения массивов из коэффициентов. Все, что ниже, -- это
		! реализация рекурсивной формулы для многочлена
		! D_n(l) = (b_n - l)*D_n-1(l) - a_n*c_n-1*D_n-2(l)
		! где коэффициенты a_i, b_i, c_i берутся из трехдиагональной
		! матрицы, которую можно изобразить схематически так:
		! b_1 c_1 0 ................................ 0
		! a_2 b_2 c_2 0 ............................ 0
		!  0  a_3 b_3 c_3 0 ........................ 0
		!       ...............................
		!  0 ...................... 0  a_n-1 b_n-1 c_n-1
		!  0 ........................... 0    a_n   b_n

		D = 0
		D(0:n-2) = -matr(n,n-1) * matr(n-1,n) * D_n2
		D(0:n-1) = D(0:n-1) + matr(n,n) * D_n1
		D(1:n) = D(1:n) - D_n1
	else
		if (n == 1) then
			D(0) = matr(1,1)
			D(1) = -1
		else
			D(0) = 1
		end if
	end if
end function minor

!================================================

recursive function find_roots(P) result(X)
	! Ищем все корни полинома, коэфф. которого хранятся в массиве P
	implicit none

	real(pr), dimension(0:), intent(in) :: P
	! корней на 1 меньше чем размер P (т.к. P содержит все коэфф.,
	! в том числе и при нулевой степени. size(P_n) = n+1,
	! где n -- степень полинома. Кол-во корней = n.
	! Индексировать корни буду в обратном порядке:
	! первый найденный будет X_n, второй -- X_n-1
	! и так далее, пока не найду последний -- X_1
	real(pr), dimension(1:size(P)-1) :: X
	! X0 -- начальное приближение. Будет const, если не придумаю
	! что-то более изощренное
	real(pr) :: X0 = 1
	integer :: n

	n = size(P) - 1
	!--------------------------------------------
	if (n > 1) then
		X(n) = newton(P,X0) ! Нахождение корня методом Ньютона с нач.приб. X0
		X(n-1:1:-1) = find_roots(gorner(P,X(n)))
	else
		X(1) = -P(0) / P(1) ! Ну у линейной функции корень найти можно
	end if
end function find_roots

function newton(P,X0) result(root)
	implicit none
	! newton() возвращает корень полинома
	! используя начальное приближение X0

	real(pr), dimension(1:), intent(in) :: P
	real(pr), intent(in) :: X0
	real(pr) :: root
	real(pr) :: x, xnew, f_x, f_xeps
	real(pr), parameter :: eps = 1e-5 * 0.1**pr
	real(pr), parameter :: f_eps = 1e1 * 0.1**pr
	integer, parameter :: nummax = 1e4 ! Максимально допустимое число итераций
	integer :: count = 0 ! счетчик числа итераций
	integer :: i, n

	n = size(P) - 1
	!--------------------------------------------
	x = X0
	xnew = X0 + 1 ! для того, чтобы в do while не было совпадения сразу

	do while (abs(xnew - x) > eps .and. count < nummax)
		x = xnew
		! следующее выражение эквивалентно x-f(x)/f'(x):
		xnew = x - fp(P,x) * f_eps / (fp(P,x+f_eps) - fp(P,x))
		count = count + 1
	end do

	if (count == nummax) then
		write(*,*) 'In newton method function count reached NUMMAX'
		write(*,*) 'Polynomial degree = ', n
	end if

	root = xnew
end function newton

pure real(pr) function fp(P,x)
	implicit none

	real(pr), dimension(0:), intent(in) :: P
	real(pr), intent(in) :: x
	integer :: i, n

	n = size(P) - 1
	!--------------------------------------------
	fp = P(n)
	do i=n-1,0,-1
		fp = fp * x + P(i)
	end do
end function fp

pure function gorner(P,x) result(P1)
	implicit none
	! Просто делим многочлен P(y) степени n на бином (y-x),
	! где x -- заранее известный корень многочлена P(y)
	! Получится многочлен степени n-1, коэффициенты
	! которого будут записаны в P1, что и будет 
	! результатом применения функции gorner

	real(pr), dimension(0:), intent(in) :: P
	real(pr), dimension(0:size(P)-2) :: P1
	real(pr), intent(in) :: x
	integer :: n, i

	n = size(P) - 1
	!--------------------------------------------
	! Ниже просто реализована схема Горнера, ничего больше.
	P1(n-1) = P(n)
	do i=n-2,0,-1
		P1(i) = x * P1(i+1) + P(i+1)
	end do
end function gorner

end module runmethods