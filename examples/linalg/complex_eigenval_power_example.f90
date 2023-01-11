program complex_eigenval_power_example
	use iso_fortran_env, only: dp => real64
	use shtawa
	implicit none
	complex(dp), dimension(:,:), allocatable :: a,x
	complex(dp) :: l
	real(dp) :: tol
	integer :: i,n,istat


	print*, "> Enter number size of square matrix:"
	read(*,*) n
	allocate(a(1:n,1:n))
	allocate(x(1:n,1))

	print*, "> Enter complex matrix, entries (re,im):"
	do i=1,n
		read(*,*) a(i,:)
	end do
	
	print*, "> Enter initial (ideally random) complex vector:"
	read(*,*) x

	call complex_eigenval_power(a=a,x=x,eigenval=l,tol=1d-6,niter=1000,istat=istat)
	if (istat /= 0) error stop istat

	print*, "Dominant eigenvalue:"
	print*, l
	print*, "Dominant eigenvector:"
	print*, x



end program