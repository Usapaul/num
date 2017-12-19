module functions

implicit none

integer, parameter :: pr = 8

real(pr) :: A, d = 5

real(pr) :: x0 = 0, xend, h, y0 = 3, z0 = -1

! Коэффициенты, получаемые в выводе уравнения
! Czn1*z_n+1 = Czn*z_n + Cyn*y_n + Csqrt*V(x_n+1/2^2 + 1) - h/2*f(zn,z_n+1)
! где f(zn,z_n+1) условно обозначает следующую дробь:
! (z_n + z_{n+1}) / (z_n + z_{n+1} + 2)
real(pr) :: Czn1, Czn, Cyn, Csqrt


contains

subroutine solve(Y,Z,n1,n2,h)
	implicit none

	real(pr), dimension(0:), intent(inout) :: Y, Z
	integer, intent(in) :: n1, n2
	real(pr), intent(in) :: h
	real(pr) :: eps = 1e-9
	real(pr) :: znew
	real(pr) :: sl1, sl2, sl3, sl4 ! Просто слагаемые из разных выражений
	integer :: i

	!--------------------------------------------
	do i = n1, n2
		!Z(i) = 1 ! Взято случайное число для начального шага итераций
		znew = 1 ! --- // --- (аналогично)
		!write(*,*) 'Z(i-1) = ', Z(i-1)
		do while (abs(znew - Z(i)) > eps)
			Z(i) = znew
			sl1 = Z(i-1) * Czn
			sl2 = Y(i-1) * Cyn
			sl3 = sqrt((x0 + h*(i-0.5_pr))**2 + 1) * Csqrt
			sl4 = (Z(i-1) + Z(i)) / (Z(i-1) + Z(i) + 2) * (-h/2)
			znew = (sl1 + sl2 + sl3 + sl4) / Czn1
		!	write(*,*) znew - Z(i)
		end do
		Z(i) = znew
		!write(*,*) 'i = ', i, 'Z(i) =  ', Z(i)
		sl1 = Y(i-1) * (1 + A*h/8)
		sl2 = Z(i-1) * (-h/4)
		sl3 = Z(i) * (-h/4)
		sl4 = sqrt((x0 + h*(i-0.5_pr))**2 + 1) * h/2
		Y(i) = (sl1 + sl2 + sl3 + sl4) * 8/(8 - A*h)
	end do
end subroutine solve

end module functions