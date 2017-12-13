module problems
use bernul

contains

recursive subroutine lejandrA(n,A)
! *** Процедура вычисляет коэфф. полинома Лежандра степени n
implicit none
integer, parameter :: mp=4
integer, intent(in) :: n
real(mp), dimension(0:n), intent(out) :: A
real(mp), dimension(0:n-1) :: P1
real(mp), dimension(0:n-2) :: P2
integer :: i
!---------------------------------------
select case(n)
case(0); A(0)=1             ! P_1(x)=1
case(1); A(0:1)=(/1.0,0.0/) ! P_2(x)=x
case default                ! Полиномы Лежандра степени больше 2-х ищем по реккурентной формуле
    call lejandrA(n-1,P1) ! Полином Лежандра степени n-1
    call lejandrA(n-2,P2) ! Полином Лежандра степени n-2
    A(0)=(2*n-1.0)/n*P1(0)
    A(1)=(2*n-1.0)/n*P1(1)
    forall (i=2:n-1) A(i)=(2*n-1.0)/n*P1(i)-(n-1.0)/n*P2(i-2)
    A(n)=-(n-1.0)/n*P2(n-2)
end select

end subroutine lejandrA

subroutine lejandrX(n,X)
! *** Процедура вычисляет корни полинома Лежандра степени n
implicit none
integer, parameter :: mp=4
integer, intent(in) :: n
real(mp), dimension(1:n), intent(out) :: X
real(mp), dimension(0:n) :: A
!---------------------------------------
select case(n)      ! Корни полинома Лежандра степени меньше 4-х нельзя найти нашим исполнением метода Бернулли,
case(1); X=0.0      ! поэтому запишем их сами
case(2); X=(/0.5773503,-0.5773503/)
case(3); X=(/0.7745967,-0.7745967,0.0/)
case default
    call lejandrA(n,A) ! Сначала получим коэффициенты полинома Лежандра
    call bernulli(A,X) ! И затем найдем корни полученного мн-на по методу Бернулли
end select

end subroutine lejandrX


end module problems