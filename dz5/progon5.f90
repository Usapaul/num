module progon5

contains

subroutine progon5fun(M,D,X)
! *** Процедура находит по методу пятиточечной прогонки
! *** вектор решений системы MX=D, где М - пятидиагональная матрица
implicit none
real, dimension(-2:,1:), intent(in) :: M
real, dimension(1:), intent(in) :: D
real, dimension(1:size(D)), intent(out) :: X
real, dimension(1:size(D)) :: B, A, P, Q, R ! См. в задание, для чего нужны эти массивы
integer :: i, n

n=size(D)
!---------------------------------------
! Дальше происходят действия, списанные с условия:
B(1)=0
A(1)=M(0,1)
P(1)=M(1,1)/A(1)
Q(1)=M(2,1)/A(1)
R(1)=D(1)/A(1)

B(2)=M(-1,2)
A(2)=M(0,2)-P(1)*B(2)
P(2)=(M(1,2)-Q(1)*B(2))/A(2)
Q(2)=M(2,2)/A(2)
R(2)=(D(2)-R(1)*B(2))/A(2)
!---------------------------------------
do i=3,n-2
    B(i)=M(1,i-1)-P(i-2)*M(2,i-2)
    A(i)=M(0,i)-P(i-1)*B(i)-Q(i-2)*M(2,i-2)
    P(i)=(M(1,i)-Q(i-1)*B(i))/A(i)
    Q(i)=M(2,i)/A(i)
    R(i)=(D(i)-R(i-1)*B(i)-R(i-2)*M(2,i-2))/A(i)
enddo
!---------------------------------------
B(n-1)=M(1,n-2)-P(n-3)*M(2,n-3)
A(n-1)=M(0,n-1)-P(n-2)*B(n-1)-Q(n-3)*M(2,n-3)
P(n-1)=(M(1,n-1)-Q(n-2)*B(n-1))/A(n-1)
Q(n-1)=0
R(n-1)=(D(n-1)-R(n-2)*B(n-1)-R(n-3)*M(2,n-3))/A(n-1)

B(n)=M(1,n-1)-P(n-2)*M(2,n-2)
A(n)=M(0,n)-P(n-1)*B(n)-Q(n-2)*M(2,n-2)
P(n)=-Q(n-1)*B(n)/A(n)
Q(n)=0
R(n)=(D(n)-R(n-1)*B(n)-R(n-2)*M(2,n-2))/A(n)
!---------------------------------------
! По условию, именно так находится вектор решений X:
X(n)=R(n)
X(n-1)=R(n-1)-P(n-1)*X(n)
do i=n-2,1,-1
    X(i)=R(i)-P(i)*X(i+1)-Q(i)*X(i+2)
enddo

end subroutine progon5fun

end module progon5