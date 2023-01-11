program gauss_example
	use iso_fortran_env, only: dp => real64
	use shtawa
	implicit none
	real(dp), dimension(:,:), allocatable :: a,acopy
	real(dp), dimension(:,:), allocatable :: x
	integer :: istat, n,m,i

	print*, "> Enter shape of matrix:"
	read(*,*) n, m
	
	allocate(a(n,n+m))
	allocate(x(n,m))
	print*, "> Enter augmented matrix a(1:n,1:n+m):"
	!!  An augmented matrix is built such that
	!!  such that each eq. i is 
	!!  a(i,1)*x(1) + ... a(i,n)*x(n) = a(i,n+1)
	do i=1,n
		read(*,*) a(i,:)
	end do

	!! To verify the solution later, we need to make a copy of a. The subroutine
	!! linsys_gauss modifies the value of 'a'.
	! allocate(acopy(n,n+1))
	! acopy(:,:) = a(:,:)
	
	call linsys_gaussjordan(a,x,istat)

	if (istat /= 0) then
		print*, "Error code: ", istat
		error stop "(X) Error in the solution"
	end if

	print*, "Solution: x = "
	print*, x


	! print*, "Verification of the solution: a @ x - b should be 0"
	! print*, [(sum(acopy(i,:n)*x(:)) ,i=1,n)] -  acopy(:,n+1)



end program