module polynomials

use init

implicit none

contains

pure real(pr) function phi(k,x)
	implicit none
	integer, intent(in) :: k
	real(pr), intent(in) :: x
	!--------------------------------------------
	select case(what_the_polynomial)
		case(1)
			phi = legendre(k,x)
		case(2)
			phi = chebyshev(k,x)
	end select
end function phi

pure recursive real(pr) function legendre(n,x) result(P)
	implicit none

	real(pr), intent(in) :: x
	integer, intent(in) :: n

	!--------------------------------------------
	if (n > 1) then
		P = (2*n-1) * x / n * legendre(n-1,x) - legendre(n-2,x) * (n-1)/n
	else if (n == 1) then
		P = x
	else
		P = 1.0_pr
	end if
end function legendre

pure recursive real(pr) function chebyshev(n,x) result(T)
	implicit none

	real(pr), intent(in) :: x
	integer, intent(in) :: n

	!--------------------------------------------
	if (n > 1) then
		T = 2 * x * chebyshev(n-1,x) - chebyshev(n-2,x)
	else if (n == 1) then
		T = x
	else
		T = 1.0_pr
	end if
end function chebyshev

end module polynomials