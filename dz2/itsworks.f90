program all
implicit none
real :: a, b ! a, b - Нижний и верхний пределы интерполирования
integer :: i, n ! i участвует как индекс в циклах. n - число интервалов между попарно соседними точками на заданном отрезке [a,b]
real, allocatable :: XkandYk(:,:) ! Массив, который следует вывести в файл в качестве результата решения задачи
real, allocatable :: Xk(:), Yk(:) ! Массив Yk состоит из значений функции в узлах Xk (начальные данные),
character(9) :: setka ! При запуске программы следует указать в качестве аргумента вид сетки. Ее название запишется в эту переменную 

! *** Программа занимается интерполяцией по данным {x,y}:
! ... {x,y} = таблица из значений узлов Xk(даны концы отрезка интерполяции и количество интервалов...
! ... и значений функции в них Yk)

call getarg(1,setka)
if (setka .eq. 'uniform  ') then            ! *** Открывается нужный файл в зависимости от выбора сетки
    open(100,file='uniform.dat')
    open(777,file='res_uniform.dat',status='replace')
elseif (setka .eq. 'chebyshev') then
    open(100,file='chebyshev.dat')
    open(777,file='res_chebyshev.dat',status='replace')
elseif (setka == 'sinus   ') then
    open(100,file='sinus.dat')
    open(777,file='res_sinus.dat')
else
    stop 'Неверно введено название сетки'
endif
!------------------------------------------------------------
111 format(2x,i8)               ! *** Чтение значений концов промежутка, количества интервалов и массива со значениями функции в узлах
read(100,111) n
allocate(Yk(0:n),Xk(0:n))
read(100,*) a, b
read(100,*) Yk

    close(100)                  ! Закрытие файла с данными
!------------------------------------------------------------
allocate (XkandYk(1:2,0:100*n))         ! *** Массив с конечными результатами заказываем у соответствующей процедуры

if (setka .ne. 'chebyshev') then
    call funi(a,b,n,Xk,Yk,XkandYk)
else
    call fche(a,b,n,Xk,Yk,XkandYk)
endif
!------------------------------------------------------------
do i=0,100*n                    ! *** Запись массива с результатами в файл
    write(777,*) XkandYk(1:2,i) ! Каждая строка этого массива содержит Xk и значение полинома в точке Xk
enddo
close(777)                      ! Закрытие файла с результатами
!------------------------------------------------------------
end program all


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
