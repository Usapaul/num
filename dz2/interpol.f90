program interpolation
implicit none
real :: a, b, h ! a, b - Нижний и верхний пределы интерполирования. h - шаг, в чебышевской сетке - вспомогательная переменная
integer :: i, n, k
real, parameter :: pi=4*atan(1.0)
real, allocatable :: XkandYk(:,:) ! Массив, который следует вывести в файл в качестве результата решения задачи
real, allocatable :: C(:), Fi(:), Xk(:) ! Массив Xk состоит из значений узлов (от 0 до n), массивы C и Fi являются вспомогательными в вычислениях
real, allocatable :: Yk(:) ! Массив, состоящий из значений функции в узлах (начальные данные)
real, allocatable :: Znam(:) ! Массив, состоящий из значений знаменателя для каждого Фи k-того
character(9) :: setka ! В эту переменную запишется название выбранной сетки


call getarg(1,setka)

111 format(2x,i8)

if (setka .eq. 'uniform  ') then

    open(100,file='uniform.dat')

	read(100,111) n
	allocate(Yk(0:n),Xk(0:n))
	read(100,*) a, b
	read(100,*) Yk

    close(100)
    
    h=(b-a)/n           ! h является шагом

    do i=0,n            ! Заполнение массива Xk значениями узлов
        Xk(i)=a+i*h !
    enddo               !

    open(777,file='res_uniform.dat',status='replace')

    
elseif (setka .eq. 'chebyshev') then

    open(100,file='chebyshev.dat')

	read(100,111) n
	allocate(Yk(0:n),Xk(0:n))
	read(100,*) a, b
	read(100,*) Yk

    close(100)

    h=pi/(2*n+2)
    do i=0,n
	Xk(i)=-cos((2*i+1)*h)
    enddo

    Xk=Xk*(b-a)/2+(a+b)/2

    open(777,file='res_chebyshev.dat',status='replace')

else 
    stop 'Название сетки введено неверно'
endif


allocate(Znam(0:n),C(0:n),Fi(0:n))


do k=0,n
    C=Xk-Xk(k)
    C(k)=1
    Znam(k)=product(C) 
enddo

Fi=Yk/Znam

deallocate(C)


allocate (XkandYk(1:2,0:100*n))
h=(b-a)/n/100 ! Если используется чебышевская сетка, то тут переменная h изменяет свое предназначение
call poli(Xk,Fi,XkandYk,h,n)

do i=0,100*n
    write(777,*) XkandYk(1:2,i) ! Каждая строка этого массива содержит Xk и значение полинома в точке Xk
enddo

close(777)


end program interpolation


subroutine poli(Xk,Fi,XkandYk,h,n)
implicit none
real :: x, Lx, h
real, dimension(0:n) :: C, Xk, L, Chis, Fi ! Массив L состоит из значений одночленов (т.е. L(x)=sum(L)). Chis - числитель для каждого Фи k-то
real, dimension(2,0:100*n) :: XkandYk
integer :: k, n, i

x=Xk(0)

do i=0,100*n

    x=x+h ! Высчитываем значения аргумента в цикле 100*n
    C=x-Xk 

    do k=0,n
        Chis(k)=product(C(0:k-1))*product(C(k+1:n))
    enddo
    
    L=Fi*Chis
    Lx=sum(L)

    XkandYk(1,i)=x
    XkandYk(2,i)=Lx

enddo


end subroutine poli

