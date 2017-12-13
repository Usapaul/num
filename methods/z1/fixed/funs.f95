module funs
implicit none

integer, parameter :: mp=8
real(mp) :: a=0.0_mp, b=1.0_mp
real(mp) :: alp1=1.0_mp, alp2=0.0_mp, Abig=1.3_mp
real(mp) :: bet1=0.76_mp, bet2=1.0_mp, Bbig=0.23_mp

contains

pure function Ffun(x) result(f)
implicit none
real(mp), intent(in) :: x
real(mp) :: f
!---------------------------------------
f=1/sqrt(x**2+1.0_mp)

end function Ffun

pure function Pfun(x) result(p)
implicit none
real(mp), intent(in) :: x
real(mp) :: p
!---------------------------------------
p=x+1.5_mp

end function Pfun

pure function Qfun(x) result(q)
implicit none
real(mp), intent(in) :: x
real(mp) :: q
!---------------------------------------
q=(x+1.0_mp)/(x+3.0_mp)

end function Qfun


end module funs