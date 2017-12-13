module rungekut
use ffunction

contains

subroutine rk4(X0,t,X)
! *** Процедура из вектора X0(t) находит вектор X(t+h)
implicit none
real(mp), dimension(1:), intent(in) :: X0
real(mp), intent(in) :: t
real(mp), dimension(1:size(X0)), intent(out) :: X
real(mp), dimension(1:size(X0)) :: K1, K2, K3, K4
!---------------------------------------
K1=h*f(t,X0)
K2=h*f(t+h/2,X0+K1/2)
K3=h*f(t+h/2,X0+K2/2)
K4=h*f(t+h,X0+K3)

X=X0+1.0_mp/6*(K1+2.0_mp*K2+2.0_mp*K3+K4)

end subroutine rk4

end module rungekut