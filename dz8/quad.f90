program quad
use problems
use lead
! *** Программа вычисляет коэфф. квадратурной формулы Гаусса
implicit none
integer, parameter :: mp=4
integer :: n, i
real(mp), allocatable :: A(:), t(:), M(:,:)
character(2) :: num
!---------------------------------------
call getarg(1,num)
if (len_trim(num)<2) then   ! Запись двузначного числа n, записанного в текстовой переменной num
    n=ichar(trim(num))-48   ! в переменную типа integer
else
    n=(ichar(num(1:1))-48)*10+ichar(num(2:2))-48
endif
!---------------------------------------
allocate(A(1:n),t(n),M(1:n+1,0:n-1))
call lejandrX(n,t)

forall (i=0:n-1) M(1:n,i)=t**(i)
M(n+1,1:n-1:2)=0
forall (i=0:n-1:2) M(n+1,i)=2.0/(i+1)

call leadfun(M,A,n)
!---------------------------------------
open(777,file='quad'//trim(num)//'.dat')
    do i=1,n
        write(777,*) A(i), t(i)
    enddo
close(777)

end program quad