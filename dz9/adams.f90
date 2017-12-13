module adams
use ABfunctions
use ffunction
use newtonmethod
use rungekut

contains

subroutine extradams(X0,t0,Xnew,X)
! *** Процедура находит вектор Xnew(t0+h) из вектора X0(t0)
implicit none
real(mp), dimension(1:), intent(in) :: X0
real(mp), intent(in) :: t0
real(mp), dimension(1:size(X0)), intent(out) :: Xnew
real(mp), dimension(-Nextradams+1:0,1:size(X0)) :: X
real(mp), dimension(1:size(X0)) :: r
integer :: j, n

n=Nextradams
!---------------------------------------
r=0
do j=-n+1,0
    r=r+A(n,-j)*f(t0+h*j,X(j,:))
enddo
Xnew=X0+h*r

end subroutine extradams

subroutine interadams(X0,t0,Xnew,X)
! *** Процедура находит вектор Xnew(t0+h) из вектора X0(t0)
implicit none
real(mp), dimension(1:), intent(in) :: X0
real(mp), intent(in) :: t0
real(mp), dimension(1:size(X0)), intent(out) :: Xnew
real(mp), dimension(-Ninteradams+2:1,1:size(X0)) :: X
real(mp), dimension(1:size(X0)) :: r
integer :: j, n

n=Ninteradams
!---------------------------------------
r=0
do j=-n+2,0
    r=r+B(n,-j)*f(t0+h*j,X(j,:))
enddo
!---------------------------------------
call newton(X0,10000,Xnew,ft)


contains

function ft(X) result(ff)
! *** Функция ft(x) равна функции f(t,x) при фиксированном t (t0+h)
implicit none
real(mp), dimension(1:), intent(in) :: X
real(mp), dimension(1:size(X)) :: ff

ff=r*h+h*B(Ninteradams,-1)*f(t0+h,X)+X0-X

end function ft

end subroutine interadams

end module adams
