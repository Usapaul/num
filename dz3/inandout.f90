program inandout 
use gauss
use jordan
use lead

! *** Программа решает систему уравнений AX=B, где A - матрица nxn, B - массив из n элементов
! *** Это главная программа, которая читает данные, дает их процедурам и записывает решение в файл

implicit none
real, allocatable :: A(:,:), B(:), X(:), M(:,:) ! M - Расширенная матрица
integer :: i, n
character(6) :: choice ! Название метода, использующегося при решении (call getarg(1,choice))
!-------------------------------------------------
open(100,file='data.dat')             ! Открываем файл, из которого берем матрицу А и вектор В

    read(100,'(2x,i8)') n ! Размер системы 

    allocate(A(1:n,1:n),B(1:n),X(1:n),M(1:n+1,1:n))

    read(100,*) A ! Матрица А, где А(1:n,i) читается из i-ой строки (i:[2,n+1])
    read(100,*) B ! Вектор В, где B(i) читается из i-ой строки (i:[n+2,2n+1])

    M(1:n,1:n)=A; M(n+1,1:n)=B ! Строим расширенную матрицу

close(100)
!-------------------------------------------------
call getarg(1,choice) ! В зависимости от выбранного метода, вызываются различные процедуры
if (choice == 'gauss ') then
    call gaussfun(M,X,n)
elseif (choice == 'jordan') then
    call jordanfun(M,X,n)
elseif (choice == 'lead  ') then
    call leadfun(M,X,n)
else
    stop 'using "getarg", choose: gauss, jordan, lead '
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
deallocate(A,B,M,X)

end program inandout
