module preobrazovanie

contains

recursive subroutine preobraz(X,Y,Wx)
implicit none
complex, dimension(0:), intent(in) :: X
complex, dimension(0:size(X)-1), intent(out) :: Y
complex, dimension(0:), intent(in) :: Wx
integer :: N, i

N=size(X)
!---------------------------------------
if (N>2) then
    call preobraz(X(0:N/2-1)+X(N/2:N-1),Y(0:N-2:2),Wx(0:N-1:2))
    call preobraz(Wx(0:N/2-1)*(X(0:N/2-1)-X(N/2:N-1)),Y(1:N-1:2),Wx(0:N-1:2))
else
    forall (i=0:1) Y(i)=X(0)+(-1)**i*X(1)
endif

end subroutine preobraz

end module preobrazovanie