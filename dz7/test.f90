program testdpfurye
implicit none
integer :: i, N=2**12
real :: a=0.0, b=20.0
real :: A1=5.0, A2=2.0, A3=1.5
real :: w1=1.0, w2=1.5, w3=4.0
real :: h,qq
!---------------------------------------
h=(b-a)/(N-1)*100.0
open(100,file='data.dat')
    write(100,*) '# ', N
    do i=0,N-1
	call random_number(qq)
        write(100,*) A1*sin(w1*i*h)+A2*sin(w2*i*h)+A3*sin(w3*i*h)+(2*qq-1)*2.0, 0.0
    enddo
close(100)


end program testdpfurye