program inandout 
use yakob
use zeidel
use relax

! *** Программа решает систему уравнений AX=B методом итераций, где A - матрица nxn, B - массив из n элементов
! *** Это главная программа, которая читает данные, дает их процедурам и записывает решение в файл

implicit none
real, allocatable :: A(:,:), B(:), X(:)
integer :: i, n
real :: eps=10e-6 ! Вычисления будут проводиться до тех пор, пока изменение нормы вектора X не станет меньше eps
character(6) :: choice ! Название метода, использующегося при решении (call getarg(1,choice))
!-------------------------------------------------
open(100,file='data.dat')             ! Открываем файл, из которого берем матрицу А и вектор В

    read(100,'(2x,i8)') n ! Размер системы

    allocate(A(1:n,1:n),B(1:n),X(1:n))

    read(100,*) A ! Матрица А, где А(1:n,i) читается из i-ой строки (i:[2,n+1])
    do i=1,n
    if (abs(A(i,i)) < ( sum(abs(A(1:i-1,i))) + sum(abs(A(i+1:n,i))) ) ) then
    stop 'Диагонального преобладания нет'
    endif
    enddo
    read(100,*) B ! Вектор В, где B(i) читается из i-ой строки (i:[n+2,2n+1])

close(100)
!-------------------------------------------------
call getarg(1,choice) ! В зависимости от выбранного метода, вызываются различные процедуры
if (choice == 'yakob ') then
    call yakobfun(A,B,X,eps)
elseif (choice == 'zeidel') then
    call zeidelfun(A,B,X,eps)
elseif (choice == 'relax ') then
    call relaxfun(A,B,X,eps)
else
    stop '***** Используется GETARG. Выберите: yakob, zeidel, relax *****'
endif
!-------------------------------------------------
open(777,file='result.dat',status='replace') ! Открываем файл и записываем в него результаты

    write(777,'("# ",i8)') n
    do i=1,n
    write(777,*) X(i)
    enddo

close(777)
!-------------------------------------------------
write(*,*) '||AX-B||=', sqrt(sum( (matmul(transpose(A),X) -B)**2 )) ! Вывод вектора невязки (погрешность)
deallocate(A,B,X)

end program inandout
