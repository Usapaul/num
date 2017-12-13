module yakob

contains

subroutine yakobfun(A,B,X,eps)
! *** Процедура находит вектор X по методу Якоби
implicit none
real, dimension(1:,1:), intent(in) :: A
real, dimension(1:), intent(in) :: B
real, dimension(1:size(B)), intent(out) :: X
integer :: i, n
real, allocatable :: G(:), D(:,:), D0(:,:), Z(:,:) ! Массивы, которые описаны в самом условии задания, D0=D**(-1)
real, intent(in) :: eps
!-------------------------------------------------
n=size(B) ! Необходимо знать число n, поэтому узнаем его у данного массива B
allocate(G(1:n),D(1:n,1:n),D0(1:n,1:n),Z(1:n,1:n))
D=0
forall (i=1:n) D(i,i)=A(i,i) ! Нулевая матрица превратилась в диагональную
D0=0
forall (i=1:n) D0(i,i)=1/D(i,i) ! Нулевая матрица превратилась в диагональную, да обратную к D
!-------------------------------------------------
Z=matmul((D-A),D0) ! Высчитывается матрица Z(j,i), где i, j - номера строк и столбцов соответственно
G=matmul(B,D0) ! Высчитывается вектор G, как в условии
X=0 ! Приравниваем X к нулю, хотя это не требуется
do while ( sqrt( sum((matmul(X,Z)+G-X)**2) ) > eps ) ! Цикл работает, пока текущий X отличается от предыдущего больше чем на eps
    X=matmul(X,Z)+G ! А вот так считается текущий X по условию
enddo
!-------------------------------------------------
deallocate(Z,D,D0,G)

end subroutine yakobfun

end module yakob