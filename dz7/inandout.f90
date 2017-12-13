program dpfurye
use preobrazovanie
use funW

implicit none
complex, allocatable :: X(:), Y(:), Wx(:)
real, allocatable :: ReImX(:,:), ReImY(:,:) ! Массивы (2,N), содержащие Re и Im чисел
integer :: N, N1, i
character(1) :: choice
!---------------------------------------
open(100,file='data.dat')
    read(100,'(2x,i30)') N
!------------------------
! Хотим, чтобы N было степенью двойки
! Если это изначально не было так, то дополняем данные нулями
    N1=1
    do while (N1<N)
        N1=N1*2
    enddo
    allocate(ReImX(2,N1))
    read(100,*) ReImX(:,1:N)
    ReImX(:,N+1:N1)=0
    N=N1
!------------------------
    allocate(X(N),Y(N),ReImY(2,N),Wx(0:N-1))
close(100)
X=cmplx(ReImX(1,:),ReImX(2,:))
deallocate(ReImX)
!---------------------------------------
call getarg(1,choice)
select case(choice)
case('1')
    forall (i=0:N-1) Wx(i)=W1(i,N)
    call preobraz(X,Y,Wx) ! Функции W1 и W2 описаны в модуле funW
case('2')
    forall (i=0:N-1) Wx(i)=W2(i,N)
    call preobraz(X,Y,Wx) ! В процедуре preobraz массивы будут индексироваться с 0
case default; stop 'Используется GETARG. Прямое ДПФ: 1; Обратное: 2'
end select
Y=Y/sqrt(N*1.0)
deallocate(X)
!---------------------------------------
open(778,file='abs.dat')
    write(778,*) '# ', N
    do i=1,N
        write(778,*) abs(Y(i))
    enddo
close(778)
ReImY(1,:)=real(Y)
ReImY(2,:)=aimag(Y)
open(777,file='result.dat')
    write(777,*) '# ', N
    do i=1,N
        write(777,*) ReImY(:,i)
    enddo
close(777)
deallocate(Y,ReImY)

end program dpfurye
