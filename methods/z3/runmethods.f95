module runmethods

use init

implicit none

contains 

subroutine itermethod(matr,X0,l_1,count)
	! Реализация степенного метода (итер.) нахождения max по модулю с.ч.
	implicit none

	real(pr), dimension(1:), intent(in) :: X0
	real(pr), dimension(1:size(X0),1:size(X0)), intent(in) :: matr
	real(pr), parameter :: eps = 0.1**(pr/2+2) ! точнее, чем ">" pr
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
		if (mod(count,500)==0) write(*,*) count ! Чтобы видеть, если долго
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
	! Построение трехдиаг. матрицы методом вращений
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
	! Построение матрицы вращений, которая нужна в методе Гивенса
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

function find_roots(P) result(X)
	! Ищем все корни полинома, коэфф. которого хранятся в массиве P
	implicit none

	real(pr), dimension(0:), intent(in) :: P
	! корней на 1 меньше чем размер P (т.к. P содержит все коэфф.,
	! в том числе и при нулевой степени. size(P_n) = n+1,
	real(pr), dimension(1:size(P)-1) :: X

	! X0 -- начальное приближение. Будет const, если не придумаю
	! что-то более изощренное... UPD: Придумал! Точнее, вынужден был
	! придумать, потому что из-за локальных экстремумов у полинома
	! у меня возникает зацикливание метода Ньютона и появление бесконечностей:
	real(pr), dimension(1:size(P)-1) :: X0 ! -- массив из начальных приближений

	integer :: n, i

	n = size(P) - 1
	!--------------------------------------------
	X0 = find_good_x0(P)
	do i=1,n
		X(i) = newton(P,X0(i))
	end do
end function find_roots

function find_good_x0(P) result(x0roots)
	! Функция возвращает x0 -- хорошее начальное приближение, которое
	! можно потом использовать для нахождения корня полинома методом Ньютона
	implicit none

	real(pr), dimension(0:), intent(in) :: P
	real(pr), dimension(1:size(P)-1) :: x0roots ! Массив из начальных приближений
	real(pr) :: l_max ! max по модулю собственное число

	! Буду работать с разбиением отрезка. Для простоты левую и правую
	! границы отрезка буду называть соответственно 'a' и 'b':
	real(pr) :: a, b	

	real(pr), parameter :: eps = 0.1**(pr/4*3)

	! число частей, на которые разбивается каждый отрезок -- m, и оно
	! выбрано практически случайным образом
	integer :: m ! выбрал n*8, считая оптимальным, ниже будет присваивание

	real(pr), dimension(:), allocatable :: t ! массив из точек
	real(pr), dimension(:), allocatable :: fun ! значения ф-ции в точках t

	integer :: count_roots ! Счетчик числа найденных интервалов с корнями
	integer :: n, i, k

	n = size(P) - 1
	m = n * 8
	!--------------------------------------------
	! Итак, мне надо найти хорошее начальное приближение к корню, зная, что
	! исследуемая функция -- полином, коэффициенты которого я знаю, и то, что
	! наибольший по модулю корень мне известен -- это наиб. по модулю с.ч., 
	! так как в этой задаче корни полинома -- это все собственные числа матрицы

	! Нужно построить сетку и записать значения полинома в ней, чтобы найти 
	! такие соседние точки сетки, где функция (полином) меняет знак. Затем 
	! следует разбить отрезок, соединяющий эти соседние две точки тоже на 
	! множество частей и в точках новой сетки искать две соседние, где 
	! функция меняет знак (уточняем месторасположение корня). Нужно когда-то
	! остановиться с таким разбиением -- когда уже можно быть уверенным, что
	! функция не успеет между двумя соседними точками образовать два локальных
	! экстремума. Так как при нахождении корня методом касательных (Ньютона)
	! в локальных точках экстремума метод зациклится, сходясь к точке
	! экстремума, но не к точке корня. Поэтому мне нужно придумать критерий
	! для остановки разбиения. 

	! Буду искать все корни, то есть сразу все интервалы, где находятся корни,
	! и останавливается разбиение именно тогда, когда число таких интервалов
	! становится равным числу возможных корней у полинома
	!--------------------------------------------
	count_roots = 0
	allocate(t(0:m),fun(0:m))

	l_max = just_l_max

	! Известно, что корни, которые я ищу у полинома, -- это
	! собственные числа матрицы, причем, я знаю наибольшее по
	! модулю. А значит, могу сузить границы поиска остальных 
	! корней до отрезка [-l_max,l_max], один из концов которого --
	! реальный корень. Чтобы по алгоритму был найден и этот корень,
	! отрезок следует расширить: расширю на 0.01eps на обоих краях:
	a = -l_max - eps*1e-2
	b = l_max + eps*1e-2

	! Строится сетка (разбиение на m частей):
	forall (i=0:m) t(i) = a + (b - a) * i / m
	forall (i=0:m) fun(i) = fp(P,t(i))	

	count_roots = count((/(fun(i)*fun(i-1) < 0,i=1,m)/))
	do while (count_roots < n)
		deallocate(t,fun)
		m = m * 2
		allocate(t(0:m),fun(0:m))

		forall (i=0:m) t(i) = a + (b - a) * i / m
		forall (i=0:m) fun(i) = fp(P,t(i))

		! Посчитаем, сколько в текущем разбиении есть интервалов с корнями:
		count_roots = count((/(fun(i)*fun(i-1) < 0,i=1,m)/))
	end do
	if (count_roots > n) stop 'COUNT_ROOTS > N' ! Ну вдруг

	! Ну всё, все интервалы с корнями найдены, возьму как начальное
	! приближение каждого корня середину отрезка, в котором, как
	! я к этому моменту уже вычислил, -- точно содержится корень
	k = 0
	do i=1,m
		if (fun(i)*fun(i-1) < 0) then
			x0roots(k+1) = (t(i) + t(i-1)) / 2
			k = k + 1
		end if
		if (k == count_roots) exit
	end do
end function find_good_x0

function newton(P,X0) result(root)
	implicit none
	! newton() возвращает корень полинома
	! используя начальное приближение X0

	real(pr), dimension(0:), intent(in) :: P
	real(pr), intent(in) :: X0
	real(pr) :: root
	real(pr) :: x, xnew, f_x, f_xeps
	real(pr), parameter :: eps = 0.1**(pr/3*4-1)
	real(pr), parameter :: h_eps = 0.1**(pr + 2)
	integer, parameter :: nummax = 1e4 ! Максимально допустимое число итераций
	integer :: count ! счетчик числа итераций
	integer :: i, n

	n = size(P) - 1
	!--------------------------------------------
	x = X0
	xnew = X0 + 1e3 * eps ! для того, чтобы в do while не было совпадения сразу

	count = 0
	do while (abs(xnew - x) > eps .and. count < nummax)
		x = xnew
		! следующее выражение эквивалентно x-f(x)/f'(x):
		xnew = x - fp(P,x) / (fp(P,x + h_eps) - fp(P,x - h_eps)) * 2 * h_eps
		count = count + 1
	end do

	if (count == nummax) then
		write(*,*) 'In newton method function count reached NUMMAX'
	end if

	root = xnew

end function newton

pure real(pr) function fp(P,x)
	! Функция-полином с коэффициентами a_n в массиве P
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

end module runmethods