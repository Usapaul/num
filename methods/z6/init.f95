module init

implicit none

integer, parameter :: pr = 4 ! задает точность real(pr)

real(pr) :: a_left, b_right ! левый и правый концы отрезка интегрирования
integer :: n_xgrid ! Число отрезков в сетке по Х
integer :: m_tgrid ! Число отрезков (слоёв) в сетке по t

real(pr) :: h, tau ! Шаг по Х и по t соответственно
real(pr) :: sigma ! sigma = tau/h^2

real(pr), dimension(:), allocatable :: Xarr ! Один слой вдоль оси X
real(pr), dimension(:), allocatable :: tarr ! Один слой вдоль оси t
real(pr), dimension(:,:), allocatable :: matr ! Рабочая матрица

contains

subroutine init_const()
	! Просто инициализация всех переменных, что у меня есть
	implicit none

	integer :: i

	!--------------------------------------------
	! количество отрезков сетки и слоев в программе вводятся вручную
	n_xgrid = 100
	m_tgrid = 1000
	! как и границы отрезка по x [a,b]:
	a_left = 0
	b_right = 1
	!--------------------------------------------
	allocate(Xarr(0:n_xgrid))
	h = (b_right - a_left) / n_xgrid
	forall (i=0:n_xgrid) Xarr(i) = a_left + i * h 

	allocate(tarr(0:m_tgrid))
	! тут должно быть tau <= 1/2A, где A = max|a(x,t)|, но у меня a(x,t)=1
	tau = h**2 / (2 * 1 + 2) ! + 2, чтобы тау было строго меньше чем h^2/2A
	sigma = tau / h**2
	forall (i=0:m_tgrid) tarr(i) = 0.0_pr + i * tau

	! matr -- рабочая матрица, хранящая в себе все значения в 
	! точках сетки по x и по t
	allocate(matr(0:n_xgrid,0:m_tgrid))

end subroutine init_const

!================================================
!================================================

! Дальше идут функции phi, psi1 и psi2, которые заданы в условии задачи
! и являются начальным, левым краевым и правым краевым условиями

elemental real(pr) function phi(x)
	! Функция phi(x) -- это начальное условие, phi(x) = u(x,0), где u -- реш.
	implicit none

	real(pr), intent(in) :: x

	!--------------------------------------------
	!phi = 1.0_pr / (1 + x**2)**2
	phi = cos(0.5_pr * x) + (1 - x) * x

end function phi


elemental real(pr) function psi1(t)
	! Функция psi(t) -- это краевое условие, psi1(t) = u(a_left,t),
	! где u -- сеточное решение задачи  
	implicit none

	real(pr), intent(in) :: t

	!--------------------------------------------
	psi1 = exp(-0.25_pr * t)

end function psi1

elemental real(pr) function psi2(t)
	! Функция psi(t) -- это краевое условие, psi2(t) = u(b_right,t),
	! где u -- сеточное решение задачи  
	implicit none

	real(pr), intent(in) :: t

	!--------------------------------------------
	psi2 = exp(-0.25_pr * t) * cos(0.5_pr)

end function psi2


!================================================

! Дальше a(x,t), b(x,t), c(x,t) и f(x,t)
! из условия задачи:
! du/dt = a(x,t)*d2u/dx2 + b(x,t)*du/dx+c(x,t)*u(x,t)+f(x,t)

! Для моего условия все очень просто будет: a=1, b=0, c=0 и f=exp(..)*(..)
! Но для универсальности и чисто для тренировки написания программ
! и тренировки навыков понимания, где лучше конкретные штуки вынести
! в отдельный модуль или в отдельные функции/процедуры, я эти начальные
! функции задаю явно как функции. Зато, если мне вдруг захочется
! похвастаться одногруппникам, что я могу решить и их шестую задачу,
! потому что мне для этого достаточно всего лишь поменять две строчки 
! в коде, то я это сделаю:D 

! Вообще мое условие такое:
! a_init_fun = 1
! b_init_fun = 0
! c_init_fun = 0
! f_init_fun = exp(-t) * (x**2 - x + 2)

! А дальше через a_ik, b_ik, c_ik и f_ik я буду вычислять
! коэффициенты A_c, B_c, C_c и D_c, см. модуль runmethods

elemental real(pr) function a_init_fun(x,t)
	implicit none

	real(pr), intent(in) :: x, t

	!--------------------------------------------
	a_init_fun = 1.0_pr

end function a_init_fun

elemental real(pr) function b_init_fun(x,t)
	implicit none

	real(pr), intent(in) :: x, t

	!--------------------------------------------
	b_init_fun = 0.0_pr

end function b_init_fun

elemental real(pr) function c_init_fun(x,t)
	implicit none

	real(pr), intent(in) :: x, t

	!--------------------------------------------
	c_init_fun = 0.0_pr

end function c_init_fun

elemental real(pr) function f_init_fun(x,t)
	implicit none

	real(pr), intent(in) :: x, t

	!--------------------------------------------
	f_init_fun = exp(-t) * (x**2 - x + 2)

end function f_init_fun


end module init