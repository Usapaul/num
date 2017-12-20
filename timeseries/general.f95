program CKABP
! В программе все делается согласно алгоритмам, изложенным 
! в книжке Витязева

use dpfurye
use functions
use functions2

implicit none

real(pr), dimension(0:N-1) :: X, t
real(pr) :: dt
integer :: i, j, k

!------------------------------------------------
! В процедуре get_series создается массив X(t_i). Там же задается
! временной шаг delta_t. После обращения к этой процедуре появляются
! массив со значениями X. известно N (оно в модуле functions),
! которое используется внутри этой процедуры и в главной программе.
call get_series(X,t)
call print_tX(t,X,'series.dat')

call delete_trend(X,t)
call print_tX(t,X,'centered.dat')

! Выполняется процедура БПФ и вычисление периодограммы
call print_tX(X,(/(0.0_pr,i=0,N-1)/),'data.dat')
call do_bpf('1','y')
call compute_Dj(dt=t(1)-t(0))

call estimate_dispersion(X,N)
call print_tX(nu,D,'period.dat')
call print_tX((/nu(0),nu(N2/4),nu((N2+1)/2)/),(/(disp_lvl,i=1,3)/),'disp.dat')

! Выполняется процедура обратного БПФ и вычисление коррелограммы
call corrgramm()
call print_tX(t,c_m,'corrgramm.dat')

call smooth()
call print_tX(nu,Dsm,'smoothed.dat')

!================================================

contains 

subroutine print_tX(t,X,filename)
! Просто печать результатов в процедуре, чтобы не висело в главной 
! программе сто циклов с открытием файла
	implicit none

	character(*), intent(in) :: filename
	real(pr), dimension(0:), intent(in) :: X, t
	integer :: N

	N = size(X)
	!--------------------------------------------
	open(100,file=filename,status='replace')
		write(100,*) '# ', N
		do k = 0,N-1
			write(100,*) t(k), X(k)
		end do
	close(100,status='keep')
	
end subroutine print_tX

end program CKABP
