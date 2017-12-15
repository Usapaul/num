module dpfurye

use preobrazovanie
use funW


contains

subroutine do_bpf(choice)
implicit none

complex(pr), allocatable :: X(:), Y(:), Wx(:)
real(pr),    allocatable :: ReImX(:,:), ReImY(:,:) ! Массивы (2,N), содержащие Re и Im чисел
integer      :: N, N1, i
character(1), intent(in) :: choice
!---------------------------------------
open(100,file='data.dat',status='old')
    read(100,'(2x,i30)') N
! Хотим, чтобы N было степенью двойки
! Если это изначально не было так, то дополняем данные нулями
    N1=1
    do while (N1<N)
        N1=N1*2
    end do
    ! И даже еще в два раза увеличим, как делается в книжке Витязева
    N1=2*N1
    ! В файле в каждой строке записана пара значений -- 
    ! Re и Im числа
    allocate(ReImX(2,N1))
    ! Но данных все равно только N, так что по факту мы считываем
    ! N, и только потом дополняем все остальное нулями до N1
    read(100,*) ReImX(:,1:N)
    ReImX(:,N+1:N1)=0.0_pr
    N=N1
    allocate(X(N),Y(N),ReImY(2,N),Wx(0:N-1))
close(100,status='keep')
X=cmplx(ReImX(1,:),ReImX(2,:))
deallocate(ReImX)
!---------------------------------------
select case(choice)
    case('1')
        forall (i=0:N-1) Wx(i)=W1(i,N)
        call preobraz(X,Y,Wx) ! Функции W1 и W2 описаны в модуле funW
    case('2')
        forall (i=0:N-1) Wx(i)=W2(i,N)
        call preobraz(X,Y,Wx) ! В процедуре preobraz массивы будут индексироваться с 0
    case default; stop 'choice: Прямое ДПФ: 1; Обратное: 2'
end select
Y=Y/sqrt(N*1.0_pr)
deallocate(X)
!---------------------------------------
open(778,file='abs.dat')
    write(778,*) '# ', N
    do i=1,N
        write(778,*) abs(Y(i))
    end do
close(778)

! Сохраним в массив ReImY соответствующие значения вещественной и
! мнимой части каждого элемента массива Y, используя встроенные
! функции преобразования типа
ReImY(1,:)=real(Y)
ReImY(2,:)=aimag(Y)
open(777,file='result.dat',status='replace')
    write(777,*) '# ', N
    do i=1,N
        write(777,*) ReImY(:,i)
    end do
close(777,status='keep')

deallocate(Y,ReImY)

end subroutine do_bpf

end module dpfurye
