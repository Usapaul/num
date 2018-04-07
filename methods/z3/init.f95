module init

implicit none

integer, parameter :: pr = 16 ! параметр, отвечающий за точность real

integer :: s ! порядок матрицы
integer :: j

real(pr), dimension(:,:), allocatable :: matr  ! "рабочая" матрица
real(pr), dimension(:,:), allocatable :: Ematr ! единичная матрица размера s
real(pr), dimension(:,:), allocatable :: diag3 ! тут будет 3-диаг. матрица

contains

subroutine matrixinit()
	! Матрица может задаваться по-разному в зависимости от 
	! условия задачи. В этой процедуре я задаю свою.
	implicit none
	! const_a/b - произвольные константы из условия задачи
	real(pr) :: const_a, const_b
	integer :: i, k
	character(1) :: test

	!--------------------------------------------
	write(*,*) 'Input from keyboard? y/n/s (all/none/only s)'
	read(*,*) test
	if (test == 'y') then
		write(*,*) 'Input s...'
		read(*,*) s
		write(*,*) 'Input const a and b...'
		read(*,*) const_a, const_b
	else
		if (test == 's') then
			write(*,*) 'Input s...'
			read(*,*) s
		else
			s = 5
		end if

		const_a = 3.0_pr
		const_b = 1.0_pr
		write(*,*) 'automatic initialization...'
		write(*,'(10x,a,i5)') 's: ', s
		write(*,'(4x,a,f5.2)') 'const_a:   ', const_a
		write(*,'(4x,a,f5.2)') 'const_b:   ', const_b
	end if

	allocate(matr(s,s), Ematr(s,s), diag3(s,s))
	Ematr = 0.0_pr
	forall (i=1:s) Ematr(i,i) = 1.0_pr

	! Замечание: первый индекс массива -- номер столбца, 
	! а второй, соответственно, -- строки

	! Именно так у меня в условии задается матрица:
	forall (i=1:s) matr(i,i) = const_a
	forall (i=1:s-1) matr(i+1,i) = (i*1.0_pr/s - 0.5_pr) * const_b
	forall (i=1:s-1) matr(i,i+1) = matr(i+1,i)
	forall (i=1:s, k=1:s, abs(i-k) > 1) matr(k,i) = 4*1.0_pr/(i+k)**2
end subroutine matrixinit

end module init