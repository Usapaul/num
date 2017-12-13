module funct

contains

subroutine fun(Xk,Yk,XkandYk,h,n)
! *** Процедура строит полином Лагранжа и счтитает значения функции на [a,b], строит по 99 точек в промежутке между данными точками
implicit none
real :: h ! Шаг
real, dimension(0:n) :: Xk, Yk, Chis, Fi ! Chis - числитель каждого Фи k-то, Fi - вспомогательный массив, соответствующий Фи в формуле полинома
real, dimension(2,0:100*n) :: XkandYk
integer :: k, n, i ! i и k участвуют только как индексы в циклах

forall (k=0:n) Fi(k)=Yk(k)/product(Xk(k)-Xk,(Xk(k)-Xk)/=0) ! См. формулу полинома Лагранжа. Фи=Ук/знаменатель (/*числитель)
forall(i=0:100*n) XkandYk(1,i)=Xk(0)+h*i                   ! Заполнение конечного массива значениями в узлах
do i=0,100*n
    forall(k=0:n) Chis(k)=product(XkandYk(1,i)-Xk(0:k-1))*product(XkandYk(1,i)-Xk(k+1:n)) ! Расчет числителя (См. формулу Лагранжа)
    XkandYk(2,i)=dot_product(Fi,Chis)                      ! Значение полинома в точке равен сумме одночленов. Каждый одночлен - Фи*"числитель"
enddo

end subroutine fun


end module funct