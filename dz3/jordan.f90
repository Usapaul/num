module jordan
use solve

contains

subroutine jordanfun(M,X,n)
! *** Процедура, строящая матрицу для решения по методу Жордана
implicit none
character(6) :: choice='jordan'
real, dimension(1:n+1,n) :: M
real, dimension(1:n) :: X
integer :: i, j, k, n
!-------------------------------------------------
do k=1,n
    forall (j=k:n+1) M(j,k)=M(j,k)/M(k,k)                 ! Диагонализация включает алгоритм,
    forall (i=1:k-1, j=k:n+1) M(j,i)=M(j,i)-M(k,i)*M(j,k) ! подходящий к методу Гаусса, но
    forall (i=k+1:n, j=k:n+1) M(j,i)=M(j,i)-M(k,i)*M(j,k) ! счет ведется для всех i/=k
enddo
M(n:n+1,n)=M(n:n+1,n)/M(n,n) ! Элемент считается вне цикла из-за определенной конструкции forall
!-------------------------------------------------
call solution(choice,M,X,n)

end subroutine jordanfun

end module jordan