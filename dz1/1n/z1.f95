program matrix
implicit none
integer n, i, n1 !A и B - матрицы n*n и n1*n1 соответственно (n=n1)
real, allocatable :: A(:,:), B(:,:), C(:,:)

!программа получает произведение матриц AxB=C

111 format(2x,i8)

open(100,file='data1.dat') !файл с матрицей A (#_n\"матрица")
read(100,111) n

    allocate(A(n,n))

    do i=1,n
	read(100,*) A(i,:) !читаем матрицу A по строкам
    enddo

close(100)

open(200,file='data2.dat')  !файл с матрицей B (#_n\"матрица")
read(200,111) n1

    allocate(B(n,n))

    do i=1,n
	read(200,*) B(i,:) !читаем матрицу B по строкам
    enddo

close(200)


    allocate(C(n,n))
    call composition(A,B,C,n)

open(777,file='result.dat',status='replace') !запись матрицы C в файл (#_n\"матрица")


    write(777,'("# ",i8)') n
    do i=1,n
	write(777,*) C(i,:) !записываем матрицу C по строкам
    enddo

close(777)


end program matrix


subroutine composition(M1,M2,R,n) !процедура рассчитывает непосредственно произведение
implicit none
integer i, j, k, n
real, dimension(n,n) :: M1, M2, R

R=0
do i=1,n
do j=1,n
do k=1,n
    R(i,j)=R(i,j)+M1(i,k)*M2(k,j) !вычисление элемента i строки j столбца конечной матрицы
enddo
enddo
enddo

end subroutine composition
