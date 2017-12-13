module solve
use makeABQ
use workmatrix
use spline

contains

subroutine solvefun(XYP,XYres)
! *** Процедура строит аппроксимирующую функцию по данным:
! *** массив узлов X, массив значений в узлах Y, массив значений весов P
implicit none
real, dimension(1:,0:), intent(in) :: XYP
real, dimension(0:size(XYP,dim=2)-1) :: X, Y, P, S, R ! Векторы S и R нужны для реализации алгоритма
real, dimension(1:2,0:100*(size(XYP,dim=2)-1)), intent(out) :: XYres
real, dimension(-1:1,0:size(XYP,dim=2)-1) :: A, B, Q ! Эти матрицы нужны для реализации алгоритма
integer :: i, n

n=size(XYP,dim=2)-1
!---------------------------------------
forall (i=0:n)
    X(i)=XYP(1,i)
    Y(i)=XYP(2,i)
    P(i)=XYP(3,i)
end forall
!---------------------------------------
call makeABQfun(X,P,A,B,Q) ! Как только получили матрицы A, B и Q, немедленно отправляем их на работу!
call workmatrixfun(A,B,Q,Y,S,R)
call splinefun(X,S,R,XYres) ! Получаем таблицу {x,y} для построения сплайна



end subroutine solvefun

end module solve