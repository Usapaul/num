module zeidel

contains

subroutine zeidelfun(A,B,X,eps)
! *** Процедура находит вектор X по методу Зейделя
implicit none
real, dimension(1:,1:), intent(in) :: A
real, dimension(1:), intent(in) :: B
real, dimension(1:size(B)), intent(out) :: X
real, dimension(1:size(B),1:size(B)) :: P ! Массивы P и Q здесь описаны так же, как в условии
real, dimension(1:size(B)) :: Q, Help     ! задания. Help - вспомогательный массив
integer :: i, j, n
real, intent(in) :: eps
!-------------------------------------------------
n=size(B) ! Необходимо знать число n, поэтому узнаем его у данного массива B
forall (i=1:n, j=1:n) P(j,i)=-A(j,i)/A(i,i) ! См. в условие задания,
forall (i=1:n) Q(i)=B(i)/A(i,i)             ! откуда взялись эти формулы
!-------------------------------------------------
Help=1 ! Просто так приравняем к единице
X=1 ! И икс тоже
do while ( sqrt( dot_product(Help,Help) ) > eps ) ! Здесь Help является разностью текущего и предыдущего X
    Help=X ! Сначала в Help запишем сам X, а потом уже запишем в него и разность
    do i=1,n
	X(i)=sum(P(1:i-1,i)*X(1:i-1))+sum(P(i+1:n,i)*X(i+1:n))+Q(i) ! См. в условие задания
    enddo
    Help=X-Help ! А вот и разность
enddo

end subroutine zeidelfun

end module zeidel