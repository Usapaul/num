module init

implicit none

integer, parameter :: pr = 8 ! параметр, отвечающий за точность real

real(pr) :: a_const ! const, которая мне попадется в условии (для K(x,t))
! левый и правый концы отрезка интегрирования a_left и b_right:
real(pr), parameter :: a_left = 0.0_pr, b_right = 1.0_pr
integer :: n_grid ! Количество узлов коллокации
! X -- массив, в котором будут храниться все узлы коллокации
real(pr), dimension(:), allocatable :: X 
! Переменной what_the_polynomial присваевается целое значение, являющееся
! индикатором для выбора полинома.
! 1 -- Лежандра
! 2 -- Чебышева
integer :: what_the_polynomial

contains

subroutine init_consts()
	! Просто инициализация постоянных из условия
	implicit none

	character(1) :: test

	!--------------------------------------------
	write(*,*) 'Input from keyboard? y/n/g (all/none/only n_grid)'
	read(*,*) test

	if (test == 'y') then
		write(*,*) 'Enter a_const...'
		read(*,*) a_const
	else
		a_const = 0.97_pr
		write(*,*) 'automatic initialization...'
		write(*,'(4x,a,f5.2)') 'const_a:   ', a_const
	end if

	if (test /= 'n') then
		write(*,*) 'Enter the number of grid nodes'
		read(*,*) n_grid
	else
		n_grid = 5
		write(*,*) 'automatic initialization...'
		write(*,'(5x,a,i5)') 'n_grid: ', n_grid
	end if
end subroutine init_consts

!================================================

elemental real(pr) function K(x,t)
	! Ядро интегрального оператора K(x,t)
	implicit none
	real(pr), intent(in) :: x, t

	!--------------------------------------------
	!K = 0.1_pr * sqrt(1.0_pr + a_const*sin(x - t))
	!K = 0.62_pr * log(1.0_pr + x*t/3)
	K = -0.72_pr * sin(x*(0.5_pr + t*t))
end function K

elemental real(pr) function f(x)
	! Функция f(x) из данного условия задачи
	implicit none
	real(pr), intent(in) :: x

	!--------------------------------------------
	!f = 1.0_pr + x
	!f = 0.62_pr + x
	f = x - 0.72_pr
end function f

end module init