module gauss
use solve

contains

subroutine gaussfun(M,X,n)
! *** Процедура приводит расширенную матрицу к треугольному виду, затем отдает ее считающей процедуре
implicit none
character(6) :: choice='gauss'
real, dimension(1:n+1,n) :: M
real, dimension(1:n) :: X
integer :: i, j, k, n
real :: eps=1.0e-4 ! При появлении ведущего элемента, модуль которого < eps, будет выдано сообщение о его близости к нулю
!-------------------------------------------------
do k=1,n-1
    if (abs(M(k,k)) < eps) write(*,*) 'Ведущий элемент близок к нулю, номер строки:', k ! Делить на близкий к нулю элемент - плохо!
    forall (j=k:n+1) M(j,k)=M(j,k)/M(k,k)                 ! Делаем матрицу треугольной по известному
    forall (i=k+1:n, j=k:n+1) M(j,i)=M(j,i)-M(k,i)*M(j,k) ! алгоритму для метода Гаусса
enddo
M(n:n+1,n)=M(n:n+1,n)/M(n,n) ! Элемент (n,n) Матрицы считается вне цикла из-за определенной конструкции forall
!-------------------------------------------------
call solution(choice,M,X,n)

end subroutine gaussfun

end module gauss