program inverse_example
	use iso_fortran_env, only: dp => real64
	use shtawa
	implicit none
	real(dp), dimension(:,:), allocatable :: a,ainv,check
	integer :: i,n

	print*, "> Enter size of (square) matrix:"
	read(*,*) n
	
	allocate(a(n,n))
	allocate(check(n,n))
	print*, "> Enter matrix to be inverted:"
	do i=1,n
		read(*,*) a(i,:)
	end do

	ainv = inverse(a)

	print*, "Inverse matrix:"
	do i = 1,n
		write(*,"(*(f8.2))") ainv(i,:)
	end do

	print*, "Checking inversion: a @ ainv must be 1"
	check = matmul(a,ainv)
	do i = 1,n
		write(*,"(*(f8.2))") check(i,:)
	end do


end program