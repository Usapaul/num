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
	allocate(X(n_grid)) ! массив из узлов сетки
	select case(what_the_polynomial)
		case(1)
			forall (i=1:n_grid) X(i) = -1.0_pr + 2.0_pr * (i-1)/(n_grid - 1)
		case(2)
			forall (i=1:n_grid) X(i) = 1.0_pr + cos(pi*(i-1)/(n_grid-1))/2.0_pr
	end select	
end subroutine polynom

!================================================


end module runmethods