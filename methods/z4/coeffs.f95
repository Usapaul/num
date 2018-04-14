module coeffs

use init
use polynomials
use integralmod

implicit none

! Индексы i_pol и k_pol мне нужны для k-го полинома, чтобы 
! менять внутренние переменные для подинтегральной функции
integer, private :: i_pol, k_pol 

! Да, я решил дополнительно поизвращаться, хочу функции, вычисляющей
! интеграл, подать как параметр просто функцию fun(t), но ведь это
! так просто не сделать, поскольку для разных i -- индекса X
! и k -- индекса полинома это будут разные функции. Я вынес определение
! "чистой" (зависящей только от одного параметра -- переменной интегрирования)
! функции отдельно -- это K_phi(t). Там я уже определяю саму функцию,
! используя индексы, которые храню вне процедуры -- i_pol и k_pol

contains

subroutine coeff(matr)
	! Создаются коэффициенты для матрицы, которую в методе коллокаций
	! придется решать как систему с правой частью из f(X(i))
	implicit none

	real(pr), dimension(1:,1:), intent(inout) :: matr

	integer :: n, i, k

	n = size(matr,1)
	!--------------------------------------------
	do i=1,n
		i_pol = i
		do k=1,n
			k_pol = k
			matr(i,k) = phi(k,X(i)) - integral(K_phi,a_left,b_right)
		end do
	end do
end subroutine coeff

pure real(pr) function K_phi(t)
	! Эта функция мне нужна для того, чтобы с правильными параметрами 
	! i и k отдать ее на интегрирование:) 
	implicit none

	real(pr), intent(in) :: t

	!--------------------------------------------
	K_phi = K(X(i_pol),t) * phi(k_pol,t)

end function K_phi

end module coeffs