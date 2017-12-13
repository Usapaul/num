module funs
implicit none

real(8) :: a=0.0d0, b=1.0d0
real(8) :: alp1=1.0d0, alp2=0.0d0, Abig=3.0d0
real(8) :: bet1=0.76d0, bet2=1.0d0, Bbig=16.08d0


!решение должно быть x**3+3x**2+x+3


contains

pure function G(x)
implicit none
real(8), intent(in) :: x
real(8) :: g
!---------------------------------------
g=x**3+3.0d0*x**2+x+3.0d0

end function G

pure function Ffun(x) result(f)
implicit none
real(8), intent(in) :: x
real(8) :: f
!---------------------------------------
f=4.0d0*x**3+11.5d0*x**2+17.0d0*x+8.5d0

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