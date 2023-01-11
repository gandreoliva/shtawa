program permutations_example
	use iso_fortran_env, only: dp=> real64
	use shtawa
	implicit none
	integer, dimension(:,:), allocatable :: permutations
	integer :: i,n
	
	n = 3
	call get_permutations_slow(n,permutations)

	do i=1,size(permutations,2)
		print*, permutations(:,i)
	end do
	
	print*, "Computed", size(permutations,2), "permutations, which should be n!"
	print*, "n! = ", int(gamma(real(n+1)))


end program