program noise100

use functions
use functions2

implicit none

real(pr), dimension(0:N-1) :: X, t
real(pr) :: dt
integer :: i, ntimes = 200

!------------------------------------------------call get_series(X,t)
open(250,file='noisepeaks.dat',status='replace')
open(251,file='disp200.dat',status='replace')
do i = 1, ntimes
	call whitenoise(X,t)
	call print_tX(t,X,'noise.dat')
	
	call print_tX(X,(/(0.0_pr,i=0,N-1)/),'data.dat')
	call do_bpf('1','y')
	call compute_Dj(dt=t(1)-t(0))
	
	call estimate_dispersion(X,N)
	write(250,*) i, maxval(D)
	write(251,*) i, disp_lvl
	!call print_tX(nu,D,'period.dat')
	!call print_tX((/nu(0),nu(N2/4),nu((N2+1)/2)/),(/(disp_lvl,i=1,3)/),'disp.dat')
	deallocate(X_FFT,D,nu)
end do
close(251,status='keep')
close(250,status='keep')

! Теперь для q = 0.05
! разница в добавлении q в процедуре estimate_dispersion
ntimes = 100
open(250,file='noisepeaks5.dat',status='replace')
open(251,file='disp100.dat',status='replace')
do i = 1, ntimes
	call whitenoise(X,t)
	call print_tX(t,X,'noise.dat')
	
	call print_tX(X,(/(0.0_pr,i=0,N-1)/),'data.dat')
	call do_bpf('1','y')
	call compute_Dj(dt=t(1)-t(0))
	
	call estimate_dispersion(X,N,q=0.05_pr)
	write(250,*) i, maxval(D)
	write(251,*) i, disp_lvl
	! Просто для себя, чтобы не строить постоянно графики и смотреть, где 
	! 5 раз пики выше уровня, а видеть на экране их количество
	if (maxval(D) > disp_lvl) then
		write(*,*) maxval(D) - disp_lvl, i
	end if
	!call print_tX(nu,D,'period.dat')
	!call print_tX((/nu(0),nu(N2/4),nu((N2+1)/2)/),(/(disp_lvl,i=1,3)/),'disp.dat')
	deallocate(X_FFT,D,nu)
end do
close(251,status='keep')
close(250,status='keep')






contains

subroutine whitenoise(X,t)
	implicit none
	real(pr), dimension(0:N-1), intent(out) :: X, t
	real(pr), parameter :: dt = 1.0_pr ! Значение шага в секундах
	integer :: k
	!--------------------------------------------
	forall (k=0:N-1) t(k) = dt * k	
	X = N_0_1(N)


end subroutine whitenoise

subroutine print_tX(t,X,filename)
! Просто печать результатов в процедуре, чтобы не висело в главной 
! программе сто циклов с открытием файла
	implicit none

	character(*), intent(in) :: filename
	real(pr), dimension(0:), intent(in) :: X, t
	integer :: N, k

	N = size(X)
	!--------------------------------------------
	open(100,file=filename,status='replace')
		write(100,*) '# ', N
		do k = 0,N-1
			write(100,*) t(k), X(k)
		end do
	close(100,status='keep')
	
end subroutine print_tX

end program noise100