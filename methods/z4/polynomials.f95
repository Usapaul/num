module polynomials

use init

implicit none

contains

pure real(pr) function phi(k,x)
	implicit none
	integer, intent(in) :: k
	real(pr), intent(in) :: x

	integer :: n 
	
	! Встал вопрос, каких степеней полиномы я должен брать. Точнее, какая
	! должна быть наименьшая -- 1 или 0 (соответственно наиб. k или k-1)
	! В отчете Ангелины от 1, но от нуля получается заметно точнее. Пока
	! не разберусь с этим вопросом наверняка, оставлю так:
	n = k - 1

	!--------------------------------------------
	select case(what_the_polynomial)
		case(1)
			phi = legendre(n,x)
		case(2)
			phi = chebyshev(n,x)
	end select
end function phi

pure recursive real(pr) function legendre(n,x) result(P)
	! Вычисление значения полинома Лежандра степени n в точке x
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
	! Вычисление значения полинома Чебышева степени n в точке x
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