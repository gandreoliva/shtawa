program quicksort_example
	use iso_fortran_env, only: dp => real64
	use shtawa
	implicit none

	real(dp), dimension(:), allocatable :: list
	list = [9,4,1,3,0]

	print*, list
	call quicksort(list)

	print*, list

end program