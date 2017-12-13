module relax

contains

subroutine relaxfun(A,B,X,eps)
! *** Процедура находит вектор X по методу релаксации
implicit none
real, dimension(1:,1:), intent(in) :: A
real, dimension(1:), intent(in) :: B
real, intent(in) :: eps
real, dimension(1:size(B)), intent(out) :: X
real, dimension(1:size(B),1:size(B)) :: P ! Массивы P и Q здесь описаны так же,
real, dimension(1:size(B)) :: Q           ! как в условии задания
integer :: i, j, n
integer, dimension(1) :: nmax             ! В этот массив запишутся координаты наибольшего по модулю элемента текущего Q
!-------------------------------------------------
n=size(B) ! Необходимо знать число n, поэтому узнаем его у данного массива B
forall (i=1:n, j=1:n) P(j,i)=-A(j,i)/A(i,i) ! См. в условие задания,
forall (i=1:n) Q(i)=B(i)/A(i,i)             ! откуда взялись такие формулы
!-------------------------------------------------
X=0 ! X должен быть равен нулю изначально (по условию)
do while ( maxval(abs(Q))  > eps ) ! На это значение изменяется координата X. Если она < eps, то счет окончен
    nmax=maxloc(abs(Q)) ! Получаем координату наибольшего по модулю элемента Q
    j=nmax(1) ! Приписываем переменной j координату максимального по модулю значения Q
    X(j)=X(j)+Q(j)                       ! См. в условие задания,
    forall (i=1:n) Q(i)=Q(i)+P(j,i)*Q(j) ! откуда взялись эти формулы
enddo

end subroutine relaxfun

end module relax