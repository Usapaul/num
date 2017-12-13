module spline

contains

subroutine splinefun(X,S,R,XYres)
! *** Процедура создает массив из точек {x,y}, по которым строится функция (сплайн)
implicit none
real, dimension(0:), intent(in):: X, S, R
real, dimension(0:size(X)-1) :: Xhelp ! Вспомогательный массив. В него будет записан X, и его придется изменять
real, dimension(1:2,0:100*(size(X)-1)), intent(out) :: XYres
real :: h, t
integer, dimension(1) :: nmax
integer :: i, j, k, n

n=size(R)-1
Xhelp=X
!---------------------------------------
! Для начала нужно построить сетку:
! Некоторые узлы в равномерной плотной сетке нам точно известны, это массив X для i=0:100*n:100
forall (i=0:100*n) XYres(1,i)=X(0)+i*(X(n)-X(0))/100/n
! Затем вычисляем в каждом узле этой сетки значение функции по определенному алгоритму:
do i=0,100*n-1
    where (Xhelp<=XYres(1,i)) Xhelp=Xhelp(n)+1.0
    nmax=minloc(Xhelp)-2 ! Здесь мы нашли наибольший X, который меньше текущего значения аргумента
    j=nmax(1)
    h=X(j+1)-X(j)
    t=(XYres(1,i)-X(j))/h
    XYres(2,i)=R(j)*(1-t)+R(j+1)*t-h**2*t*(1-t)/6.0*( (2.0-t)*S(j)+(1.0+t)*S(j+1)  )
enddo
h=X(n)-X(n-1)
XYres(2,100*n)=R(j+1)


end subroutine splinefun

end module spline