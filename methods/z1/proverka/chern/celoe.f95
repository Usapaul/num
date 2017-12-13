
real(8) :: h=3.34d0
real(16) :: q=3.34q0

write(*,*) h/2
write(*,*) 2+h
write(*,*) 2*h
write(*,*) 2.0+h
write(*,*) 5.2/2.0
write(*,*) '----------------'
write(*,*) 5.2/2.3d0
write(*,*) 5.2d0/2.3
write(*,*) 5.2/2.3*1.0d0
write(*,*) 5.2d0/(2.3*1.0d0)
write(*,*) 2.3*1.0d0
write(*,*) 2.3d0
write(*,*) 5.2d0/2.3d0
write(*,*) '---------------'
write(*,*) h**21
write(*,*) h**21.0d0
write(*,*) q**21
write(*,*) q**21.0q0

end