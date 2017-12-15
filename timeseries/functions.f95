module functions


implicit none

! Здесь N получает свое значение
integer,  parameter :: N = 230  
integer,  parameter :: pr = 8
real(pr), parameter :: pi = 4.0_pr*atan(1.0_pr)

contains

subroutine get_series(X,t)
! Процедура возвращает временной ряд некоторый длины и массив отсчетов времени. 
	implicit none

	real(pr), parameter :: dt = 1.0_pr ! Значение шага в секундах
	real(pr), dimension(0:N-1), intent(out) :: X, t
	integer :: k

	!--------------------------------------------
	! В функцию linear_with_noise заложены параметры, создающие ряд, описанный в книжке
	forall (k=0:N-1) t(k) = dt * k
	X = linear_with_noise(t)

end subroutine get_series


function linear_with_noise(t) result(X)
	implicit none

	real(pr), dimension(0:), intent(in) :: t
	real(pr), dimension(0:size(t)-1) :: X
	integer :: N, k

	! Параметры для задания линейного тренда
	real(pr), parameter :: alpha = 0.1_pr, beta = 0.05_pr
	! Параметр A_1, который упомянут в книжке -- для создания 
	! гармонического компонента, плюс фаза phase, частота nu1 --
	! все для того же гармонического компонента
	real(pr), parameter :: A1 = 1._pr, nu1 = 0.1_pr, phase = 0._pr
	! Параметр отношения S/N -- gamma, необходимый для нахождения сигмы
	real(pr), parameter :: gamma = 0.50_pr
	real(pr) :: sigma_x = sqrt(A1**2 / (2*gamma))

	N = size(t)
	!--------------------------------------------
	forall (k=0:N-1) X(k) = alpha + beta*t(k) + A1*cos(2*pi*nu1*t(k) - phase)
	! Добавим нормально распределенный шумовой компонент
	X = X + sigma_x * N_0_1(N)

end function linear_with_noise

function N_0_1(n) result(Z)
! Функция получает на вход значение N и возвращает массив
! из N элементов, составляющих выборку случайных величин, 
! распределенных по нормальному закону N(0,1) 
	integer, intent(in) :: n
	real(pr), dimension(:), allocatable :: Z
	real(pr), dimension(:), allocatable :: u, v

	allocate(Z(n),u(n),v(n))
	call random_seed()
	call random_number(u)
	call random_seed()
	call random_number(v)
	Z = sqrt(-2*log(u))*cos(2*pi*u)

end function N_0_1

!----------------------------------------------------------

subroutine delete_trend(X,t)
	implicit none

	real(pr), dimension(0:), intent(inout) :: X
	real(pr), dimension(0:), intent(in) :: t
	real(pr) :: a, b
	integer :: N, k

	N = size(X)
	!--------------------------------------------
	a = (N*sum(t*X) - sum(t)*sum(X)) / (N*sum(t**2) - (sum(t))**2)
	b = (sum(X) - a*sum(t)) / N

	X = X - a*t - b

end subroutine delete_trend

end module functions