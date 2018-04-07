module fornewton

use init

implicit none

real(pr), dimension(:), allocatable :: Coef

contains

function newton(P,X0) result(root)
	implicit none
	! newton() возвращает корень полинома
	! используя начальное приближение X0

	real(pr), dimension(1:), intent(in) :: P
	real(pr), intent(in) :: X0
	real(pr) :: root
	real(pr) :: x, xnew
	real(pr), parameter :: eps = 1e-3 * 0.1**pr
	integer :: nummax = 1e4 ! Максимально допустимое число итераций
	integer :: count = 0 ! счетчик числа итераций

	!-------------------------------------------- 
	! Я тупо извращаюсь, усложняю алгоритм и увеличиваю число строк
	! кода, но мне очень хочется, чтобы я писал просто f(x) вместо
	! непосредственного вычисления значения полинома в точке x

	! В Coef будут храниться коэфф. "текущего" полинома
	allocate(Coef(0:size(P))) 
	Coef = P 
	!--------------------------------------------
	x = X0
	xnew = X0 + 1 ! для того, чтобы в do while не было совпадения сразу

	do while (abs(xnew - x) > eps .and. count < nummax)
		x = xnew
		xnew = x - 
	end do

end function newton


end module fornewton