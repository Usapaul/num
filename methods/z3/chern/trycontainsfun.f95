program xxx

integer :: a = 10, n = 5


write(*,*) yyy(a,n)

contains

integer function yyy(a,n)

integer :: a, n

yyy = a*n
yyy = zzz(n)*a

contains

function zzz(n)

integer :: n

zzz = n*2

end function zzz

end function yyy


end program xxx