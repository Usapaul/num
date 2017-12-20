module funW

use functions

contains

pure function W1(q,p) result(W)
	implicit none

	integer, intent(in) :: q, p
	complex(pr) :: W
	!---------------------------------------
	W=exp(-2.0_pr*pi*cmplx(0,1)*q/p)

end function W1

pure function W2(q,p) result(W)
	implicit none

	integer, intent(in) :: q, p
	complex(pr) :: W
	!---------------------------------------
	W=exp(2.0_pr*pi*cmplx(0,1)*q/p)

end function W2

end module funW