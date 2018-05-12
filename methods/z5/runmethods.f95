module runmethods

use init
use polynomials


implicit none

real(pr), parameter :: pi = 4.0_pr * atan(1.0_pr)

! n_for_legendre_fun нужно для превращения функции -- полином Лежандра
! в зависящую только от одного аргумента x, а "текущая" степень полинома
! записывается в эту переменную:
integer, private :: n_for_legendre_fun

contains

subroutine polynom()
	! В этой процедуре только лишь создается массив X из корней соотв. полинома
	implicit none

	integer :: i
	!--------------------------------------------
	select case(what_the_polynomial)
		case(1)
			! Для полинома Лежандра:
			X = legendre_roots(n_grid)*(b_right-a_left)/2 + (b_right+a_left)/2
		case(2)
			! Для полинома Чебышева:
			X = chebyshev_roots(n_grid)*(b_right-a_left)/2 + (b_right+a_left)/2
	end select	
end subroutine polynom

!================================================

function legendre_roots(n) result(X)
	! Получение массива из корней полинома Лежандра степени n
	use init, only: pr
	implicit none

	integer, intent(in) :: n
	real(pr), dimension(n) :: X
	real(pr) :: x0 ! x0 служит начальным приближением для корня
	integer :: i

	!--------------------------------------------
	n_for_legendre_fun = n
	do i=1,n
		x0 = cos(pi*(4*i-1)/(4*n+2)) ! приближение для корней п. Лежандра
		X(i) = root_by_newton(legendre_n,x0)
	end do
end function legendre_roots

pure real(pr) function legendre_n(x)
	! Чистая функция от x -- полином Лежандра заданной степени n.
	! Эта чистая функция нужна для нахождения корней полинома Лежандра
	! методом Ньютона
	use init, only: pr
	implicit none

	real(pr), intent(in) :: x

	!--------------------------------------------
	legendre_n = legendre(n_for_legendre_fun,x)

end function legendre_n

pure function chebyshev_roots(n) result(X)
	! Получение массива из корней полинома Чебышева степени n
	use init, only: pr
	implicit none

	integer, intent(in) :: n
	real(pr), dimension(n) :: X
	integer :: i

	!--------------------------------------------
	forall (i=1:n) X(i) = cos(pi * (2*i - 1) / (2*n))

end function chebyshev_roots

!================================================

function root_by_newton(f,x0) result(root)
	! Нахождение корня по методу Ньютона с нач. приближ. x0
	use init, only: pr
	implicit none

	interface
		pure real(pr) function f(t)
			use init, only: pr
			real(pr), intent(in) :: t
		end function f
	end interface

	real(pr), intent(in) :: x0
	real(pr) :: x, xnew, root
	real(pr), parameter :: eps = 0.1**(pr/2 + 3)
	! Производную буду считать по формуле (f(x+h) - f(x-h))/2h
	! для которой остаток <= f'''(c)*h^2/6
	! поэтому беру h (h_eps) такое, чтобы при итерациях
	! обеспечить сходимость приближений корня с точностью eps
	real(pr), parameter :: h_eps = 0.1**(pr + 2) 
	integer, parameter :: nummax = 1e4 ! Максимально допустимое число итераций	
	integer :: count ! счетчик числа итераций

	!--------------------------------------------
	count = 0
	x = x0
	xnew = x0 + 1e-2 ! для того, чтобы в do while не было совпадения сразу
	do while (abs(xnew - x) > eps)
		x = xnew
		! следующее выражение эквивалентно x-f(x)/f'(x):
		xnew = x - f(x) / (f(x + h_eps) - f(x - h_eps)) * 2 * h_eps
		count = count + 1
	end do

	if (count == nummax) then
		write(*,*) 'In newton method function count reached NUMMAX'
	end if	
	root = xnew
end function root_by_newton

end module runmethods