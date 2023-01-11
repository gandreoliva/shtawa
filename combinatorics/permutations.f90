recursive subroutine get_permutations_slow(n,permutations)
	!! Generates the permutations of all the elements of a given set by generating the
	!! Cartesian product and discarding repeated elements (this is why it's slow).

	!! Algorithm from rosettacode.org

	!! Theory
	!! ------
	!! Permutations are ways of arraging of all elements of a set, such that no element is
	!! repeated. The number of permutations is n! Note that this subroutine doesn't compute
	!! the so-called k-permutations (kPn, from elementary textbooks).
	!! Example: if the set is {x,y,z}, the permutations are
	!!  {[x,y,z],[x,z,y],[y,x,z],[y,z,x],[z,x,y],[z,y,x]}.
	!! This subroutine generates the indices of the permutations starting at 1,
	!! which can be used to get the permutations of the elements of any array that
	!! contains the set.
	integer, save :: ind_element = 1
		!! index of the currently examined element in the set
	integer, intent(in) :: n
		!! Total number of elements in the (original) set
	integer, dimension(:), allocatable, save :: permutation
		!! Work array (communicates information from outer calls into deeper calls)
	integer, dimension(:,:), intent(inout), allocatable :: permutations
		!! Output array
	integer :: i

	if (.not. allocated(permutations)) allocate(permutations(n,0))
	if (.not. allocated(permutation)) allocate(permutation(n))

	if (ind_element > n) then
		!! We are done, all elements have been examined
		permutations = reshape([permutations,permutation],&
			& shape=[size(permutations,1),size(permutations,2)+1])
	else
		!! One nested do loop per element in the set
		do i = 1,n
			!! If the element is not already in the permutation, continue
			if (.not. any(permutation(:ind_element-1) == i)) then
				permutation(ind_element) = i
				ind_element = ind_element + 1
				call get_permutations_slow(n,permutations)
			end if
		end do
	end if

	ind_element = ind_element - 1

end subroutine