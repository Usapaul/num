module functions2

use functions
use dpfurye

implicit none

complex(pr), dimension(:), allocatable :: X_FFT
real(pr),    dimension(:), allocatable :: D
real(pr),    dimension(:), allocatable :: nu
real(pr) :: disp_lvl
integer :: N2

real(pr), dimension(:), allocatable :: c_m
integer, parameter :: NtoNx = 10
integer :: Nx = N/NtoNx ! то самое N* из книжки
real(pr), dimension(:), allocatable :: Dsm

contains

subroutine compute_Dj(dt)
	implicit none

	real(pr), dimension(:,:),  allocatable :: FFT_ReIm
	real(pr), intent(in) :: dt
	real(pr) :: delta_nu 
	integer :: j, k

	!--------------------------------------------
	open(100,file='result.dat',status='old')
		read(100,'(2x,i20)') N2
		allocate(X_FFT(0:N2-1),FFT_ReIm(2,0:N2-1),D(0:(N2+1)/2))
		do j = 0,N2-1
			read(100,*) FFT_ReIm(1:2,j)
		end do
	close(100,status='keep')

	X_FFT = cmplx(FFT_ReIm(1,:),FFT_ReIm(2,:))

	forall (j=0:(N2+1)/2)
		D(j) = 1._pr/N**2 * (real(X_FFT(j))**2 + aimag(X_FFT(j))**2)
	end forall

	delta_nu = 1._pr / (N2*dt)
	allocate(nu(0:(N2+1)/2))
	forall (j=0:(N2+1)/2) nu(j) = delta_nu * j

end subroutine compute_Dj

subroutine estimate_dispersion(X,N)
	implicit none

	real(pr), parameter :: X1 = 9.0_pr
	real(pr), dimension(0:), intent(in) :: X
	integer, intent(in) :: N
	real(pr) :: disp
	integer :: j

	!------------------------------------------------------
	disp = 1._pr/(N-1) * sum(X**2)
	disp_lvl = disp * X1 / N

end subroutine estimate_dispersion

subroutine corrgramm()
	implicit none

	real(pr), dimension(:), allocatable :: Xj2
	real(pr), dimension(:,:),  allocatable :: FFT2_ReIm
	complex(pr), dimension(:), allocatable :: X2_FFT2

	integer :: N22, i, j, m

	!--------------------------------------------
	open(100,file='abs.dat',status='old')
		read(100,'(2x,i20)') N22
		if (N22 /= N2) stop 'N22 in abs.dat /= N2' 
		allocate(Xj2(0:N2-1))
		do i = 0,N2-1
			read(100,*) Xj2(i)
		end do
	close(100,status='keep')

	open(200,file='data.dat',status='replace')
		write(200,*) '# ', N2
		do i = 0,N2-1
			write(200,*) Xj2(i)**2, 0.0_pr
		end do
	close(200,status='keep')

	call do_bpf('2','n')

	open(100,file='result.dat',status='old')
		read(100,'(2x,i20)') N22
		allocate(X2_FFT2(0:N22-1),FFT2_ReIm(2,0:N22-1))
		do j = 0,N22-1
			read(100,*) FFT2_ReIm(1:2,j)
		end do
	close(100,status='keep')
	X2_FFT2 = cmplx(FFT2_ReIm(1,:),FFT2_ReIm(2,:))

	allocate(c_m(0:N-1))

	forall (m=0:N-1)
		c_m(m) = 1._pr/N * real(X2_FFT2(m))
	end forall

end subroutine corrgramm

subroutine smooth()
	implicit none

	real(pr), parameter :: a = 0.25_pr
	real(pr), dimension(:), allocatable :: c_mW, Wm
	real(pr), dimension(:,:),  allocatable :: FFT_ReIm
	complex(pr), dimension(:), allocatable :: cm_FFT	
	integer :: m, i, j, N22

	!--------------------------------------------
	allocate(c_mW(0:N2-1),Wm(0:Nx-1))
	forall (m=0:Nx-1)
		Wm(m) = (1._pr - 2._pr*a) + 2._pr * a * cos(pi*m/Nx)
	end forall
	c_mW(0:Nx-1) = c_m(0:Nx-1) * Wm
	c_mW(Nx:N2-1) = 0._pr

	allocate(Dsm(0:(N2+1)/2))

	open(200,file='data.dat',status='replace')
		write(200,*) '# ', N2
		do i = 0,N2-1
			write(200,*) c_mW(i), 0.0_pr
		end do
	close(200,status='keep')	

	call do_bpf('1','n')

	open(100,file='result.dat',status='old')
		read(100,'(2x,i20)') N22
		allocate(cm_FFT(0:N22-1),FFT_ReIm(2,0:N22-1))
		do j = 0,N22-1
			read(100,*) FFT_ReIm(1:2,j)
		end do
	close(100,status='keep')
	cm_FFT = cmplx(FFT_ReIm(1,:),FFT_ReIm(2,:))

	forall (j=0:(N2+1)/2)
		Dsm(j) = 1._pr/Nx * (2._pr * real(cm_FFT(j)) - c_mW(0))
	end forall

end subroutine smooth

end module functions2