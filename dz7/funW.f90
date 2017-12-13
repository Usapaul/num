module funW

contains

pure function W1(q,p) result(W)
implicit none
integer, intent(in) :: q, p
complex :: W
!---------------------------------------
W=exp(-2.0*(4*atan(1.0))*cmplx(0,1)*q/p)

end function W1

pure function W2(q,p) result(W)
implicit none
integer, intent(in) :: q, p
complex :: W
!---------------------------------------
W=exp(2.0*(4*atan(1.0))*cmplx(0,1)*q/p)

end function W2

end module funW