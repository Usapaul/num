module diag3dim1

contains

subroutine diag3dim1fun(A,B,C)
! *** Процедура получает произведение трехдиагональной матрицы на вектор
implicit none
real, dimension(-1:,1:), intent(in) :: A
real, dimension(1:), intent(in) :: B
real, dimension(1:size(B)), intent(out) :: C
integer :: i, n

n=size(B)
!---------------------------------------
C(1)=A(0,1)*B(1)+A(1,1)*B(2)
forall (i=2:n-1) C(i)=A(-1,i)*B(i-1)+A(0,i)*B(i)+A(1,i)*B(i+1)
C(n)=A(-1,n)*B(n-1)+A(0,n)*B(n)

end subroutine diag3dim1fun

end module diag3dim1