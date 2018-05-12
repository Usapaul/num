module runmethods

use init

implicit none

! A_c, ..., D_c -- коэффициенты, см. комментарий в main перед вызовом get_ABCD
real(pr), dimension(:,:), allocatable :: A_c, B_c, C_c, D_c


contains

subroutine get_ABCD()
	! Процедура получает коэффициенты A_c, B_c, C_c, D_c
	! см. комментарий вверху в этом модуле при объявлении этих коэффициентов
	implicit none

	! коэффициенты a,b,c,f -- функции из самого условия задачи, которое
	! задается как:
	! du/dt = a(x,t)*d2u/dx2 + b(x,t)*du/dx+c(x,t)*u(x,t)+f(x,t)
	! и они заданы как (a/b/c/f)_init_fun в модуле init
	real(pr), dimension(:,:), allocatable :: a, b, c, f

	integer :: i, k
	integer :: n, m ! = n_xgrid и m_tgrid соответственно

	n = n_xgrid
	m = m_tgrid
	!--------------------------------------------
	! от 1 и до #-1 потому, что есть начальные и граничные условия
	! а последний слой по t высчитывается из предыдущего -- m_tgrid-1 и всё
	allocate(A_c(1:n-1,0:m))
	allocate(B_c(1:n-1,0:m))
	allocate(C_c(1:n-1,0:m))
	allocate(D_c(1:n-1,0:m))
	
	allocate(a(1:n-1,0:m))
	allocate(b(1:n-1,0:m))
	allocate(c(1:n-1,0:m))
	allocate(f(1:n-1,0:m))
	!--------------------------------------------
	forall (i=1:n-1, k=0:m)
		a(i,k) = a_init_fun(Xarr(i),tarr(k))
		b(i,k) = b_init_fun(Xarr(i),tarr(k))
		c(i,k) = c_init_fun(Xarr(i),tarr(k))
		f(i,k) = f_init_fun(Xarr(i),tarr(k))
	end forall

	forall (i=1:n-1, k=0:m)
		A_c(i,k) = tau * (a(i,k)/h**2 + b(i,k)/(2*h))
		B_c(i,k) = -1 + tau * (2*a(i,k)/h**2 - c(i,k))
		C_c(i,k) = tau * (a(i,k)/h**2 - b(i,k)/(2*h))
		D_c(i,k) = tau * f(i,k)
	end forall

end subroutine get_ABCD


subroutine explicit_method(u,a,b,c,d)
	! Процедура, согласно явному методу, заполняет
	! матрицу u (являющуюся сеточным решением) при известных
	! коэффициентах a, b, c, d, которые вместе составляют уравнение:
	! u(i,k+1) = a(i,k)*u(i+1,k) - b(i,k)*u(i,k) + c(i,k)*u(i-1,k) + d(i,k)
	implicit none

	real(pr), dimension(0:,0:), intent(inout) :: u
	real(pr), dimension(1:,0:), intent(in) :: a, b, c, d

	integer :: i, k, n, m

	n = size(u,dim=1) - 1
	m = size(u,dim=2) - 1
	!--------------------------------------------
	! Известно, что первая строчка матрицы u уже заполнена -- это начальная 
	! функция phi(x) = u(x,0), а также заполнены края

	do k=0,m-1
		forall (i=1:n-1)
			u(i,k+1) = a(i,k)*u(i+1,k)-b(i,k)*u(i,k)+c(i,k)*u(i-1,k)+d(i,k)
		end forall
	end do
end subroutine explicit_method

!================================================

subroutine implicit_method(u,a,b,c,d)
	! Процедура, согласно неявному методу, заполняет
	! матрицу u (являющуюся сеточным решением) при известных
	! коэффициентах a, b, c, d, которые вместе составляют уравнение:
	! a(i,k)*u(i+1,k) - b(i,k)*u(i,k) + c(i,k)*u(i-1,k) + d(i,k) = u(i,k-1)
	implicit none

	real(pr), dimension(0:,0:), intent(inout) :: u
	real(pr), dimension(1:,0:), intent(in) :: a, b, c, d

	integer :: i, k, n, m

	n = size(u,dim=1) - 1
	m = size(u,dim=2) - 1
	!--------------------------------------------
	! Известно, что первая строчка матрицы u уже заполнена -- это начальная 
	! функция phi(x) = u(x,0), а также заполнены края

	do k=1,m
		! Еще есть "свободный член", который не вошел в d, и равен
		! значению u в i-ой точке предыдущего (k-1-го) слоя по t, и для 
		! отправки массива d функции progonka еще его нужно взять со знаком 
		! минус, потому что в функции progonka d -- это массив -- столбец
		! свободных членов уравнений и находится в правой части каждого ур-я
		u(:,k) = progonka(a(:,k),b(:,k),c(:,k),-d(:,k)+u(1:n-1,k-1), &
							& 0._pr,u(0,k),0._pr,u(n,k))
	end do

end subroutine implicit_method

function progonka(a,b,c,d,kappa1,nu1,kappa2,nu2) result(Y)
	! Функция возвращает вектор Y решений следующей системы уравнений:
	! Y(0) = kappa1 * Y(1) + nu1
	! {a(i)*Y(i+1) - b(i)*Y(i) + c(i)*Y(i-1) = d(i)} {i=1:n-1}
	! Y(n) = kappa2 * Y(n-1) + nu2

	! Описанный здесь код с точностью до обозначений реализует алгоритм,
	! описанный в моем отчете по первой задаче по вычам

	implicit none

	real(pr), intent(in) :: kappa1, nu1
	real(pr), intent(in) :: kappa2, nu2

	! a, b, c, d -- массивы коэффициентов в системе уравнений
	! в методе прогонки, появляющиеся во всех строчках, кроме
	! нулевой и последней, там уравнения с двумя неизвестными,
	! и там коэффициенты -- каппа и ню
	real(pr), dimension(1:), intent(in) :: a, b, c, d

	real(pr), dimension(0:size(a)+1) :: Y ! Вектор решения системы
	real(pr), dimension(0:size(a)+1) :: U, V ! вспомогательные массивы
	integer :: i, n

	n = size(a)+1
	!--------------------------------------------
	! Ищем решение в виде Y(i) = U(i) * Y(i+1) + V(i),
	! где U(i), V(i) -- неизвестные пока коэффициенты
	
	U(0) = kappa1
	V(0) = nu1

	! Все согласно тому, что написано в моем отчете по первой задаче по вычам:
	do i=1,n-1
		U(i) = a(i) / (b(i) - c(i) * U(i-1))
		V(i) = (c(i) * V(i-1) - d(i)) / (b(i) - c(i) * U(i-1))
	end do

	Y(n) = (-kappa2 * V(n-1) - nu2) / (kappa2 * U(n-1) - 1)

	do i=n-1,0,-1
		Y(i) = U(i) * Y(i+1) + V(i)
	end do

end function progonka


end module runmethods