module funs
implicit none

real(8) :: a=0.0d0, b=1.0d0
real(8) :: alp1=1.0d0, alp2=0.0d0, Abig=1.3d0
real(8) :: bet1=0.76d0, bet2=1.0d0, Bbig=0.23d0

contains

pure function Ffun(x) result(f)
implicit none
real(8), intent(in) :: x
real(8) :: f
!---------------------------------------
f=1/sqrt(x**2+1.0d0)

end function Ffun

pure function Pfun(x) result(p)
implicit none
real(8), intent(in) :: x
real(8) :: p
!---------------------------------------
p=x+1.5d0

end function Pfun

pure function Qfun(x) result(q)
implicit none
real(8), intent(in) :: x
real(8) :: q
!---------------------------------------
q=(x+1.0d0)/(x+3.0d0)

end function Qfun


end module funs