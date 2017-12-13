module ABfunctions
use mprecision

contains

function A(n,j)
implicit none
integer, intent(in) :: n, j
integer :: i
real(mp) :: A
real(mp) :: integral
!---------------------------------------
call integralf(0.0_mp,1.0_mp,Aintfun,integral)
A=(-1.0_mp)**j/product((/(i,i=1,j)/))/product((/(i,i=1,n-1-j)/))*integral

contains

function Aintfun(z)
implicit none
real(mp), intent(in) :: z
real(mp) :: Aintfun
integer :: i
!---------------------------------------
Aintfun=product((/(z+i,i=0,n-1)/))/(z+j)

end function Aintfun

end function A

!=======================================

function B(n,j)
implicit none
integer, intent(in) :: n, j
integer :: i
real(mp) :: B
real(mp) :: integral
!---------------------------------------
call integralf(0.0_mp,1.0_mp,Bintfun,integral)
B=(-1.0_mp)**(j+1)/product((/(i,i=1,j+1)/))/product((/(i,i=1,n-2-j)/))*integral

contains

function Bintfun(z)
implicit none
real(mp), intent(in) :: z
real(mp) :: Bintfun
integer :: i
!---------------------------------------
Bintfun=product((/(z+i,i=-1,n-2)/))/(z+j)

end function Bintfun

end function B

!=======================================

subroutine integralf(a0,b0,f,r)
! *** Процедура считает интеграл от f на [a0,b0] по методу Гаусса (n=5)
implicit none
integer, parameter :: n=5
real(mp), intent(in) :: a0, b0
real(mp), intent(out) :: r
real(mp), dimension(1:n) :: A, t
integer :: i
character(2) :: num
interface
function f(x)
use mprecision
real(mp), intent(in) :: x
real(mp) :: f
end function
end interface
!---------------------------------------
A(1)=0.236926794; t(1)=0.906179845
A(2)=0.236926794; t(2)=-0.906179845
A(3)=0.478628963; t(3)=0.538469315
A(4)=0.478628904; t(4)=-0.538469315
A(5)=0.568888545; t(5)=0.00000000
! Рассчет интеграла по формуле Гаусса. Функция f масштабирована под отрезок [-1,1]
r=0
do i=1,n
    r=r+A(i)*(b0-a0)/2*f(t(i)*(b0-a0)/2+(a0+b0)/2)
enddo

end subroutine integralf


end module ABfunctions
