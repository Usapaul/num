module ffunction
use mprecision

implicit none
integer, parameter :: ndim=3
integer :: i
real(mp) :: tdata=4.0_mp
real(mp), dimension(ndim) :: Xdata=(/(1.0_mp,i=1,ndim)/)
real(mp) :: h=0.1_mp**(2)
integer :: Nextradams=4
integer :: Ninteradams=4

contains

function f(t,X) result(Y)
implicit none
real(mp), dimension(1:), intent(in) :: X
real(mp), intent(in) :: t
real(mp), dimension(1:size(X)) :: Y
integer :: i, n

n=size(X)
!---------------------------------------
Y=0
Y(1)=cos(X(1))*abs(X(1))-0.1_mp**(5)*t  ! Такая функция придумана для примера

end function f

end module ffunction