program CKABP
! В программе все делается согласно алгоритмам, изложенным 
! в книжке Витязева

use dpfurye
use functions

implicit none


! X1 и X2 -- критические значения для разделения шумового и 
! детерминированного компонентов временного ряда, отвечающие 
! уровню значимости q
! X2=0.0, так как рассматриваем случай, когда частота заранее неизвестна
real(pr), parameter :: q = 0.01_pr, X1 = 9.0_pr, X2 = 0.0_pr
real(pr), dimension(0:N-1) :: X, t
real(pr) :: dt
integer :: i, j, k
complex(pr), dimension(:), allocatable :: X_FFT
integer :: N2
real(pr), dimension(:), allocatable :: D

!------------------------------------------------
! В процедуре get_series создается массив X(t_i). Там же задается
! временной шаг delta_t. После обращения к этой процедуре появляются
! массив со значениями X. известно N (оно в модуле functions),
! которое используется внутри этой процедуры и в главной программе.
call get_series(X,t)
call print(X,t,'series.dat')

call delete_trend(X,t)
call print(X,t,'centered.dat')


! Выполняется процедура БПФ
open(100,file='data.dat',status='replace')
	write(100,*) '# ', N
	do k = 0, N-1
		write(100,*) X(k), 0.0_pr
	end do
close(100,status='keep')

call do_bpf('1')

open(100,file='result.dat',status='old')
	read(100,'(2x,i10)') N2
	allocate(X_FFT(0:N2-1),D(0:N2-1))
	do j = 0,N2-1
		read(100,*) X_FFT(j)
	end do
close(100,status='keep')

forall (j=0:(N2+1)/2) D(j) = 1/N**2 * (real(X_FFT)**2 + aimag(X_FFT)**2)




!================================================

contains 

subroutine print(X,t,filename)
! Просто печать результатов в процедуре, чтобы не висело в главной 
! программе сто циклов с открытием файла
	implicit none

	character(*), intent(in) :: filename
	real(pr), dimension(0:), intent(in) :: X, t
	integer :: N

	N = size(X)
	!--------------------------------------------
	open(100,file=filename,status='replace')
		do k = 0,N-1
			write(100,*) t(k), X(k)
		end do
	close(100,status='keep')
	
end subroutine print

end program CKABP
