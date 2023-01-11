program determinant_example
	use iso_fortran_env, only: dp => real64
	use shtawa
	implicit none
	real(dp), dimension(:,:), allocatable :: a
	integer :: i,n

	print*, "> Enter size of matrix:"
	read(*,*) n

	allocate(a(1:n,1:n))

	print*, "> Enter matrix 'a':"
	do i=1,n
		read(*,*) a(i,:)
	end do

	print*, "determinant(a) = ", determinant(a)


end program