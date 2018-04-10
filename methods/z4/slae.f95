module slae

use init, only: pr

implicit none

contains

function solve(A,B) result(X)
	implicit none

	real(pr), dimension(1:,1:), intent(in) :: A
	real(pr), dimension(1:), intent(in) :: B
	real(pr), dimension(1:size(B)) :: X

	real(pr), dimension(1:size(B),1:size(B)+1) :: M ! Расширенная матрица
	integer :: n

	n = size(B)
	!--------------------------------------------
	if (size(A,1) /= size(B)) then
		stop 'Function solve: A B have incompatible dimensions'
	end if

	M(1:n,1:n) = A
	M(:,n+1) = B
	! Решение системы будет выполняться по методу Гаусса.
	! Расширенная матрица M подается на вход процедуре, которая 
	! отдает обратно массив X решений СЛАУ:
	call gauss_slae(M,X)
end function solve

subroutine gauss_slae(M,X)
	implicit none

	real(pr), dimension(1:,1:), intent(inout) :: M
	real(pr), dimension(size(M,1)), intent(out) :: X

	! Так как в методе Гаусса с выбором ведущего элемента придется
	! переставлять столбцы местами, нужно записывать порядок переменных,
	! чтобы потом знать, в каком столбце какая переменная оказалась.
	! За хранение информации о перестановке переменных отвечает Order
	integer, dimension(size(M,1)) :: Order
	integer, dimension(2) :: index
	integer :: n, i, j, k, p, q

	n = size(M,1)
	Order = (/(i,i=1,n)/)
	!--------------------------------------------
	do k=1,n-1
		! В index хранятся два числа -- значения индексов строки и 
		! столбца. Так как функция maxloc выдает индексы, считая их
		! от единицы, то, понятное дело, я должен прибавить предыдущие
		! индексы, а это k-1, т.е. k-ый элемент (левый верхний угол
		! текущей матрицы) для функции maxloc имеет индексы (/1,1/).
		! Чтобы получить индексы в нормальном виде, надо сложить:
		! index(1:2) = (/k-1,k-1/) + maxloc(...)
		index = k - 1 + maxloc(abs(M(k:n,k:n)))
		!Просто для удобства введу обозначения p и q:
		p = index(1)
		q = index(2)

		! Перестановка строк:
		if (p /= k) then
			forall (i=k:p:p-k) M(i,:) = M(p+k-i,:)
		end if
		! Перестановка столбцов:
		if (q /= k) then
			forall (i=k:q:q-k)
				M(:,i) = M(:,q+k-i)
				Order(i) = Order(q+k-i) ! Следим за порядком переменных
			end forall
		end if

		! Собственно, сам метод Гаусса:
		forall (j=k:n+1) M(k,j) = M(k,j) / M(k,k)
		forall (i=k+1:n, j=k:n+1) M(i,j) = M(i,j) - M(k,j) * M(i,k)
	end do
	! (n,n)-ый элемент матрицы не нуждался в перестановке, поэтому не попал
	! в цикл, и, вообще говоря, с него можно уже начать решение СЛАУ AX = B
	! так как он выражает прямое равенство A(n,n) * X(n) = B(n). Но его
	! тоже нужно поделить на самого себя (вместе с B(n)=M(n,n+1)), вот:
	forall (j=n:n+1) M(n,j) = M(n,j) / M(n,n)

	! На этом метод Гаусса не закончился, только матрица была приведена
	! к треугольному виду, теперь же нужно получить решение системы:
	do i=n,1,-1
		X(i) = M(i,n+1) - dot_product(M(i,i+1:n),X(i+1:n))
	end do
	! Приведем массив в порядок
	X(Order) = X
end subroutine gauss_slae

end module slae