module runmethods

use init
use polynomials


implicit none

real(pr), parameter :: pi = 4.0_pr * atan(1.0_pr)

contains

subroutine polynom()
	implicit none

	integer :: i
	!--------------------------------------------
	select case(what_the_polynomial)
		case(1)
			! Для полинома Лежандра:
			!forall (i=1:n_grid) X(i) = -1.0_pr + 2.0_pr * (i-1)/(n_grid - 1)
			forall (i=1:n_grid)
				X(i) = a_left + (b_right - a_left) * (i-1)/(n_grid - 1)
			end forall
		case(2)
			! Для полинома Чебышева:
			forall (i=1:n_grid) X(i) = 1.0_pr + cos(pi*(i-1)/(n_grid-1))/2.0_pr
	end select	
end subroutine polynom

end module runmethods