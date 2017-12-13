module solve

contains

subroutine solution(choice,M,X,n)
! *** Процедура высчитывает массив решений, используя уже преобразованную расширенную матрицу
implicit none
integer :: i, n
real, dimension(1:n+1,n) :: M
real, dimension(1:n) :: X
character(6) :: choice ! Процедура, вызывающая solution, даст о себе знать, благодаря переменной "choice"
!-------------------------------------------------
if (choice /= 'jordan') then ! Только для метода Жордана счет иксов ведется по-другому
    do i=n,1,-1
	X(i)=M(n+1,i)-dot_product(M(i+1:n,i),X(i+1:n))
    enddo
else
    do i=1,n
	X(i)=M(n+1,i)
    enddo
endif

end subroutine solution

end module solve