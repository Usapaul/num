module givexk
use funct

contains

subroutine funi(a,b,n,Xk,Yk,XkandYk)
! *** Процедура считает значения узлов равномерной сетки и передает управление считающей процедуре
implicit none
integer :: i, n
real :: a, b, h
real, dimension(0:n) :: Xk, Yk
real, dimension(2,0:100*n) :: XkandYk

h=(b-a)/n                  ! Шаг = расстояние между двумя соседними узлами
forall (i=0:n) Xk(i)=a+i*h ! Заполнение массива значениями узлов
h=h/100                    ! Считающей процедуре понадобится знать только шаг вместо значений на концах

call fun(Xk,Yk,XkandYk,h,n)

end subroutine funi


subroutine fche(a,b,n,Xk,Yk,XkandYk)
! *** Процедура считает значения узлов чебышевской сетки и передает управление считающей процедуре
implicit none
integer :: i, n
real :: a, b, h, seredka, amplituda ! См. ниже пояснение для странных переменных
real, dimension(0:n) :: Xk, Yk
real, dimension(2,0:100*n) :: XkandYk

h=4*atan(1.0)/(2*n+2)                                 ! Вспомогательная переменная. См. формулу задания чебышевской сетки
seredka=(a+b)/2                                       ! Середина отрезка, на котором интерполируем
amplituda=seredka-a                                   ! Амплитуда, нужная косинусу, чтобы он от середки дотягивал до краев [a,b]
forall(i=0:n) Xk(i)=-amplituda*cos((2*i+1)*h)+seredka ! Заполнение значениями в узлах чебышевской сетки массива Xk
h=(b-a)/n/100                                         ! Переменная h становится шагом, который нужен считающей процедуре

call fun(Xk,Yk,XkandYk,h,n)

end subroutine fche

end module givexk