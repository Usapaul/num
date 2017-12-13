module ffunction

contains

function f(x)
implicit none
integer, parameter :: mp=4
real(mp), intent(in) :: x
real(mp) :: f
integer :: i
!---------------------------------------
f=2*x**2

end function f

end module ffunction