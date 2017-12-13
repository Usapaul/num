module workmatrix
use diag3comp
use diag3dim1
use progon5

contains

subroutine workmatrixfun(A,B,Q,Y,S,R)
! *** Процедура выполняет всю работу с матрицами
! *** и в итоге получает векторы S и R
implicit none
real, dimension(-1:,0:), intent(in) :: A, B, Q ! Матрицы (трехдиагональные) названы так же, как в задании
real, dimension(0:), intent(in) :: Y
real, dimension(0:size(Y)-1), intent(out) :: S, R ! Векторы R и S названы так же, как в задании
real, dimension(0:size(Y)-1) :: Help0          ! Help0 и Help1 -
real, dimension(-2:2,0:size(Y)-1) :: Help1     ! вспомогательные массивы.
real, dimension(-1:1,0:size(Y)-1) :: BT        ! BT - транспонированная (трехдиагональная) матрица B
integer :: i, n

n=size(Y)-1
!---------------------------------------
! Так должна выглядеть транспонированная трехдиаг. матрица B, записанная как массив 3xn:
BT(-1,0)=0; BT(0,0)=0; BT(1,0)=B(-1,1)
forall (i=1:n-1)
    BT(-1,i)=B(1,i-1)
    BT(0,i)=B(0,i)
    BT(1,i)=B(-1,i+1)
end forall
BT(-1,n)=B(1,n-1); BT(0,n)=0; BT(1,n)=0
!---------------------------------------
call diag3compfun(Q,BT,Help1)  ! В алгоритме нужно перемножить матрицы Q и BT, они трехдиагональные, поэтому считает процедура

call diag3compfun(B,Q,Help1)


forall (i=0:n) BT(-1:1,i)=Help1(-1:1,i) ! Help1 - трехдиаг., потому что Q - диаг. Здесь и далее BT есть Q*BT
call diag3compfun(B,BT,Help1)

Help1=Help1*6 ! По алгоритму нужно сложить A+6*Help1, что и происходит здесь
Help1(-1:1,:)=Help1(-1:1,:)+A(-1:1,:)
call diag3dim1fun(B,Y,Help0)   ! Умножение трехдиагональной матрицы на вектор
Help0=Help0*6 ! По алгоритму: правая часть системы равна 6*BY
call progon5fun(Help1,Help0,S) ! Решение системы уравнений (относительно S) по методу пятиточечной прогонки
!---------------------------------------
call diag3dim1fun(BT,S,Help0)  ! Напоминание: К этому моменту BT является матрицей Q*BT
R=Y-Help0


end subroutine workmatrixfun

end module workmatrix