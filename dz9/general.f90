program diffur
use rungekut
use adams
use ffunction

! *** Дана система: dX/dt=f(t,X). Производится интегрирование системы с постоянным шагом интегрирования

implicit none
real(mp) :: t
real(mp), dimension(1:size(Xdata)) :: X, Xnew
!real(mp), dimension(1:Nextradams,1:size(Xdata)) :: Xae
!real(mp), dimension(1:Ninteradams,1:size(Xdata)) :: Xai
real(mp), allocatable :: Xae(:,:), Xai(:,:)
character(2) :: choice ! Переменная служит для выбора метода (rk, ae, ai)
integer :: j

allocate(Xae(-Nextradams+1:0,1:size(Xdata)),Xai(-Ninteradams+2:1,1:size(Xdata)))


call getarg(1,choice)
!---------------------------------------
select case(choice)
case('rk')
open(100,file='rk.dat')
    write(100,*) 0.0_mp, Xdata
    t=0.0_mp; X=Xdata ! Это начальные данные
    do while (t<tdata+h) ! Итегрирование ведется на полуинтервале [0,tdata+h)
        call rk4(X,t,Xnew) ! Xnew - новое вычисленное значение X в точке t
        write(100,*) t, Xnew
        t=t+h ! Новый шаг
        X=Xnew
    enddo
close(100)
!---------------------------------------
case('ae')
open(200,file='ae.dat')
    write(200,*) 0.0_mp, Xdata
    t=(Nextradams-1)*h
    Xae(-Nextradams+1,:)=Xdata ! Для начального набора X самый первый член берется равным начальным данным
    do j=-Nextradams+1,-1
        call rk4(Xae(j,:),t+h*j,Xae(j+1,:)) ! Создается начальный набор X
        write(200,*) t, Xae(j+1,:)
    enddo
    X=Xae(0,:)
    do while (t<tdata+h)
        call extradams(X,t,Xnew,Xae)
        forall (j=-Nextradams+1:-1) Xae(j,:)=Xae(j+1,:)
        Xae(0,:)=Xnew
        write(200,*) t, Xnew
        t=t+h
        X=Xnew
    enddo
close(200)
!---------------------------------------
case('ai')
open(300,file='ai.dat')
    write(300,*) 0.0_mp, Xdata
    t=(Ninteradams-2)*h
    Xai(-Ninteradams+2,:)=Xdata ! Для начального набора X самый первый член берется равным начальным данным
    do j=-Ninteradams+2,-1
        call rk4(Xai(j,:),t+h*j,Xai(j+1,:)) ! Создается начальный набор X
        write(300,*) t, Xai(j+1,:)
    enddo
    X=Xai(0,:)
    do while (t<tdata+h)
        call interadams(X,t,Xnew,Xai)
        forall (j=-Ninteradams+2:0) Xai(j,:)=Xai(j+1,:)
        Xai(1,:)=Xnew
        write(300,*) t, Xnew
        t=t+h
        X=Xnew
    enddo
close(300)
!---------------------------------------
case default
    stop 'Используется GETARG. Выберите: rk, ae или ai'
end select

end program diffur