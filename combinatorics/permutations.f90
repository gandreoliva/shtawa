subroutine get_permutations_slow(n,permutations)
    !! Gets the permutations of all the elements of a given set by generating the
    !! Cartesian product and discarding repeated elements (this is why it's slow).

    !! This is a wrapper to the function gen_permutations_slow, for convenience of call.
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
    
    integer, intent(in) :: n
        !! Number of elements of the set 
    integer, dimension(:,:), allocatable, intent(out) :: permutations
        !! Array shape(n,number_of_permutations) that contains one permutation of the set
        !! in each column.
    integer, dimension(:), allocatable :: permutation

    allocate(permutations(n,0))
    allocate(permutation(n))
    call gen_permutations_slow(1,n,permutation,permutations)

end subroutine

recursive subroutine gen_permutations_slow(ind_element,n,permutation,permutations)
    !! This function generates all the permutations recursively. 
    !! For an introduction to recursion, see the documentation of cartprod
    integer, intent(in) :: ind_element
        !! index of the currently examined element in the set
    integer, intent(in) :: n
        !! Total number of elements in the (original) set
    integer, dimension(:), intent(inout), allocatable :: permutation
        !! Work array (communicates information from outer calls into deeper calls)
    integer, dimension(:,:), intent(inout), allocatable :: permutations
        !! Output array
    integer :: i

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
                call gen_permutations_slow(ind_element+1,n,permutation,permutations)
            end if
        end do
    end if


end subroutine