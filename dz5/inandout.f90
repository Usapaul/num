program givespline
use solve

! *** Программа берет экспериментальные значения Y в узлах X и их веса P,
! *** строит по ним аппроксимирующую функцию, и записывает в файл
! *** ее значения на равномерной сетке в 100 раз гуще

implicit none
real, allocatable :: XYP(:,:), XYres(:,:) ! Пояснение к массивам дано в блоке их задания
integer :: i, n ! Количество промежутков (узлов - n+1 штук)
!---------------------------------------
open(100,file='data.dat')
    read(100,'(2x,i8)') n
    allocate(XYP(1:3,0:n),XYres(1:2,0:100*n)) ! Массив XYres состоит из пар {x,y} на равномерной сетке
    do i=0,n
	read(100,*) XYP(1:3,i) ! В файле с данными в каждой строке записаны X(i) (X(n+1)>X(n)), Y(i), P(i)
    enddo
close(100)
!---------------------------------------
call solvefun(XYP,XYres)

deallocate(XYP)
open(777,file='result.dat')
    do i=0,100*n
    write(777,*) XYres(1:2,i)
    enddo
close(777)
deallocate(XYres)

end program givespline