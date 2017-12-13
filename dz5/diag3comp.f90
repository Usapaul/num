module diag3comp

contains

subroutine diag3compfun(A,B,C)
! *** Процедура умножает друг на друга трехдиагональные матрицы
! *** и выдает пятидиагональную в формате {(-1:1,n)} ---> {(-2:2,n)}
implicit none
real, dimension(-1:,1:), intent(in) :: A, B
real, dimension(-2:2,1:size(A,dim=2)), intent(out) :: C
integer :: i, n

n=size(A,dim=2)
!---------------------------------------
! Рассчет элементов пятидиагональной матрицы ведется по пяти "полоскам":
C(-2,1)=0; C(-2,2)=0
forall (i=3:n) C(-2,i)=A(-1,i)*B(-1,i-1)

C(-1,1)=0; C(-1,2)=A(-1,2)*B(0,1)+A(0,2)*B(-1,2)
forall (i=3:n) C(-1,i)=A(-1,i)*B(0,i-1)+A(0,i)*B(-1,i)

C(0,1)=A(0,1)*B(0,1)+A(1,1)*B(-1,2)
forall (i=2:n-1) C(0,i)=A(-1,i)*B(1,i-1)+A(0,i)*B(0,i)+A(1,i)*B(-1,i+1)
C(0,n)=A(-1,n)*B(1,n-1)+A(0,n)*B(0,n)

forall (i=1:n-2) C(1,i)=A(0,i)*B(1,i)+A(1,i)*B(0,i+1)
C(1,n-1)=A(0,n-1)*B(1,n-1)+A(1,n-1)*B(0,n); C(1,n)=0

forall (i=1:n-2) C(2,i)=A(1,i)*B(1,i+1)
C(2,n-1)=0; C(2,n)=0

end subroutine diag3compfun

end module diag3comp
