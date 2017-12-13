program intgauss
use ffunction
! *** Метод численного интегрирования Гаусса
implicit none
integer, parameter :: mp=4
real(mp) :: a=1.0, b=6.0
real(mp) :: result
integer :: n=5
!---------------------------------------
call integralf(a,b,n,f,result)

write(*,*) result

end program intgauss

subroutine integralf(a0,b0,n,f,r)
! *** Процедура считает интеграл от f на [a0,b0] по методу Гаусса
implicit none
integer, parameter :: mp=4
real(mp), intent(in) :: a0, b0
integer, intent(in) :: n
real(mp), intent(out) :: r
real(mp), dimension(1:n) :: A, t
integer :: i
character(2) :: num
interface
function f(x)
integer, parameter :: mp=4
real(mp), intent(in) :: x
real(mp) :: f
end function
end interface
!---------------------------------------
! В этом блоке происходит запись двузначного n (переменная integer) в текстовую переменную num
num='  '
if (n>9) then
    write(num(1:1),'(i1)') n/10; write(num(2:2),'(i1)') mod(n,10)
else
    write(num(1:1),'(i1)') n
endif
!---------------------------------------
open(100,file='quad'//trim(num)//'.dat')
    do i=1,n
        read(100,*) A(i), t(i)
    enddo
close(100)
!---------------------------------------
! Рассчет интеграла по формуле Гаусса. Функция f масштабирована под отрезок [-1,1]
r=0
do i=1,n
    r=r+A(i)*(b0-a0)/2*f(t(i)*(b0-a0)/2+(a0+b0)/2)
enddo

end subroutine integralf
