program tridiag_crout_example
	use iso_fortran_env, only: dp => real64
	use shtawa
	implicit none
	real(dp), dimension(:,:), allocatable :: a,acopy
	real(dp), dimension(:), allocatable :: x
	integer :: i,n

	print*, "> Enter number size of square matrix"
	read(*,*) n
	allocate(a(1:n,1:n+1))
	allocate(acopy(1:n,1:n+1))
	allocate(x(1:n))
	print* ,"> Enter augmented tridiagonal matrix:"
	do i=1,n
		read(*,*) a(i,:)
	end do

	acopy = a

	call tridiag_crout(a,x)

	print*, "Solution:"
	print*, x

	print*, "Verification of the solution: a @ x - b should be 0"
	print*, [(sum(acopy(i,:n)*x(:)) ,i=1,n)] -  acopy(:,n+1)


end program