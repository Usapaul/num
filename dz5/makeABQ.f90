module makeABQ

contains

subroutine makeABQfun(X,P,A,B,Q)
! *** Процедура получает страшные матрицы A, B, Q, необходимые в решении
implicit none
real, dimension(0:), intent(in) :: X, P
real, dimension(-1:1,0:size(X)-1), intent(out) :: A, B, Q
integer :: i, n

n=size(X)-1
!---------------------------------------
A(0,0)=2*(X(1)-X(0)); A(-1,0)=0.0; A(1,0)=0.0
A(-1,1)=0; A(0,1)=2*(X(2)-X(0)); A(1,1)=X(2)-X(1)
forall (i=2:n-2)
    A(-1,i)=X(i)-X(i-1)
    A(0,i)=2*(X(i+1)-X(i-1))
    A(1,i)=X(i+1)-X(i)
end forall
A(-1,n-1)=X(n-1)-X(n-2); A(0,n-1)=2*(X(n)-X(n-2)); A(1,n-1)=0.0
A(0,n)=2*(X(n)-X(n-1)); A(-1,n)=0.0; A(1,n)=0.0
!---------------------------------------
B(-1:1,0)=0
forall (i=1:n-1)
    B(-1,i)=1.0/(X(i)-X(i-1))
    B(0,i)=-( 1.0/(X(i)-X(i-1))+1.0/(X(i+1)-X(i)) )
    B(1,i)=1.0/(X(i+1)-X(i))
end forall
B(-1:1,n)=0
!---------------------------------------
forall (i=0:n)
    Q(0,i)=1/P(i); Q(-1,i)=0; Q(1,i)=0
end forall

end subroutine makeABQfun

end module makeABQ