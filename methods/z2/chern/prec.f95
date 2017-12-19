

integer, parameter :: pr = 8

real(pr) :: a, one = 1.0_pr


write(*,*) 1, 1.0, 1*one, 1.0*one, 1.0_pr
write(*,*) 2, 2.0, 2*one, 2.0*one, 2.0_pr

write(*,*) 11/7, 11.0/7, 11.0/7.0, 11*one/7, 11.0*one/7, 11.0_pr/7
write(*,*) 11/7, 11/7.0, 11.0/7.0, 11/(7*one), 11/(7.0*one), 11/7.0_pr


write(*,*) 13.7, 13.7*one, 13.7_pr
write(*,*) 13.7_pr/7.1_pr
write(*,*) 13.7/7.1, 13.7_pr/7.1, 13.7/7.1_pr

write(*,*) 8 - 1.2345_pr
write(*,*) (8 - 1.2345_pr)/4
write(*,*) (8 - 1.2345_pr)/4._pr


end 