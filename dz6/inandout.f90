program newtonslau
use newtonmethod
use ffunction

implicit none
real, allocatable :: X0(:), X(:) !X0 - вектор начальных приближений
integer :: n, nummax=1000 ! nummax - максимально возможное количество итераций
!---------------------------------------
open(100,file='data.dat') ! Вектор начальных приближений берется из файла (в первой строке #_n - число уравнений)
    read(100,'(2x,i8)') n
    allocate(X0(n),X(n))
    read(100,*) X0
close(100)
!---------------------------------------
call newton(X0,nummax,X,Ffun) ! Получаем решение системы, используя начальное приближение X0
write(*,*) '|f(X)|=', sqrt(dot_product(Ffun(X),Ffun(X))) ! Вывод модуля f(X), где X - найденное решение
deallocate(X,X0)

end program newtonslau