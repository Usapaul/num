module coeffs

use init
use polynomials
use integralmod

implicit none

! Индексы i_pol и k_pol мне нужны для k-го полинома, чтобы 
! менять внутренние переменные для подинтегральной функции
integer, private :: k_for_funs, m_for_funs
real(pr), private :: x_for_funs

integer :: m_from_main

! Да, я решил дополнительно поизвращаться, хочу функции, вычисляющей
! интеграл, подать как параметр просто функцию fun(t), но ведь это
! так просто не сделать, поскольку для разных i -- индекса X
! и k -- индекса полинома это будут разные функции. Я вынес определение
! "чистой" (зависящей только от одного параметра -- переменной интегрирования)
! функции отдельно -- это K_phi(t) (и др.). Там я уже определяю саму функцию,
! используя индексы, которые храню вне процедуры -- k_for_funs, m_for_funs
! По аналогичной причине мне нужен и x_for_funs
! А также мне понадобилось менять значение m в главной программе при 
! построении вектора правой части системы уравнений. Текущее значение m
! будет храниться в переменной m_from_main

contains

subroutine coeff(matr)
	! Создаются коэффициенты для матрицы, которую в методе коллокаций
	! придется решать как систему с правой частью из f(X(i))
	implicit none

	real(pr), dimension(1:,1:), intent(inout) :: matr

	integer :: n, i, k, m

	n = size(matr,1)
	!--------------------------------------------
	! Так получается, что я просто программу из 4-ой задачи использую в
	! пятой, мне ведь нужно только матрицу коэффициентов перестроить
	! Тут я действительно меняю формулы для коэффициентов, и еще 
	! функцию K_phi тоже меняю, так как тут у меня вообще 
	! двойной интеграл присутствует

	do k=1,n
		k_for_funs = k
		do m=1,n
			m_for_funs = m
			matr(k,m) = a_const * integral(phi_x_phi,a_left,b_right)
			matr(k,m) = matr(k,m) + integral(phi_x_integral,a_left,b_right)
		end do
	end do
end subroutine coeff


pure real(pr) function phi_x_phi(x)
	! Эта функция мне нужна для того, чтобы с правильными параметрами 
	! m и k отдать ее на интегрирование:) 
	! А именно, произведение phi(k,x)*phi(m,x)
	implicit none

	real(pr), intent(in) :: x

	!--------------------------------------------

	phi_x_phi = phi(k_for_funs,x) * phi(m_for_funs,x)

end function phi_x_phi


real(pr) function phi_x_integral(x)
	! Эта функция мне нужна для того, чтобы с правильными параметрами 
	! m и k отдать ее на интегрирование:) 
	! А именно, перемножение phi(m,x) на интеграл от K(x,t)*phi(k,x)dt
	implicit none

	real(pr), intent(in) :: x

	!--------------------------------------------
	
	x_for_funs = x
	phi_x_integral = phi(m_for_funs,x) * integral(K_phi,a_left,b_right)

end function phi_x_integral

pure real(pr) function K_phi(t)
	! Эта функция мне нужна для того, чтобы с правильными параметрами 
	! m и k отдать ее на интегрирование:) 
	! А именно, произведение K(x,t)*phi(k,t)
	implicit none

	real(pr), intent(in) :: t

	!--------------------------------------------

	K_phi = K(x_for_funs,t) * phi(k_for_funs,t)

end function K_phi

!================================================

pure real(pr) function f_x_phi(x)
	! Эта функция мне нужна для того, чтобы с правильными параметрами 
	! m и k отдать ее на интегрирование:) 
	! А именно, произведение f(x)*phi(m,x)
	implicit none

	real(pr), intent(in) :: x

	!--------------------------------------------
	! значение m хранится в этом модуле, но меняется в главной программе
	f_x_phi = f(x) * phi(m_from_main,x)

end function f_x_phi

end module coeffs