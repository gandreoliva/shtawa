recursive subroutine quicksort(arr)
	!! Sorts an array into increasing order with the quicksort algorithm
	!! Given an array, one element is chosen and the others are partitioned into two
	!! subsets: the subset of elements less than the partition element, and the subset
	!! of elements greater than or equal to it. The same process is applied recursively
	!! to the two subsets until each subset has only one element: the recursion stops.
	!! Reference: The C Programming Language, Sect. 4.10 (Recursion)
	real(wp), dimension(:), intent(inout) :: arr
	integer :: i, last, n

	n = size(arr)

	if (n <= 1) return

	!! Move partition element to the beginning
	call swap(arr(1),arr((1+n)/2))
	last = 1
	
	do i=2,n
		!! order the elements against the first element in the partition
		if (arr(i) < arr(1)) then
			last = last + 1
			call swap( arr(last), arr(i) )
		end if
	end do

	call swap( arr(1), arr(last) ) !! restore partition element
	call quicksort(arr(1:last-1))
	call quicksort(arr(last+1:n))

contains
	elemental subroutine swap(x,y)
		real(wp), intent(inout) :: x,y
		real(wp) :: temp
		temp = x
		x = y
		y = temp
	end subroutine

end subroutine