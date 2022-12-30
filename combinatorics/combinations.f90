subroutine get_combinations_slow(n,k,combinations)
    !! Gets the combinations of all the elements of a given set by generating a
    !! Cartesian product and discarding elements (this is why it's slow).

    !! This is a wrapper to the function gen_combinations_slow, for convenience of call.
    !! Algorithm from rosettacode.org

    !! Theory
    !! ------
    !! Combinations are selections of k elements out of a set of n elements, in such
    !! a way that the order does not matter. For example, given 3 fruits for a smoothie,
    !! {banana, pineapple, papaya}, and 2 selections, one can choose in only 3 ways:
    !! {banana, pineapple}, {banana, papaya}, {pineapple, papaya} (once it's in the
    !! blender, who cares what the order was!).
    !! The number of combinations can be computed with the binomial coefficient C(n,k)
    !! C(n,k) = n!/(k!(n-k)!)
    
    integer, intent(in) :: n
        !! Number of elements of the set 
    integer, intent(in) :: k
        !! Number of choices
    integer, dimension(:,:), allocatable, intent(out) :: combinations
        !! Array shape(n,number_of_combinations) that contains one combination of the set
        !! in each column.
    integer, dimension(:), allocatable :: combination

    allocate(combinations(k,0))
    allocate(combination(k))
    call gen_combinations_slow(1,k,n,combination,combinations)

end subroutine

recursive subroutine gen_combinations_slow(ind_choice,k,n,combination,combinations)
    !! This function generates all the combinations recursively. 
    !! For an introduction to recursion, see the documentation of cartprod
    integer, intent(in) :: ind_choice
        !! index of the currently examined choice in the set
    integer, intent(in) :: n
        !! Total number of elements in the (original) set
    integer, intent(in) :: k
        !! Total number of choices
    integer, dimension(:), intent(inout), allocatable :: combination
        !! Work array (communicates information from outer calls into deeper calls)
    integer, dimension(:,:), intent(inout), allocatable :: combinations
        !! Output array
    integer :: i

    if (ind_choice > k) then
        !! We are done, all choices have been examined
        combinations = reshape([combinations,combination],&
            & shape=[size(combinations,1),size(combinations,2)+1])
    else
        !! One nested do loop per element in the set
        do i = 1,n
            !! At first, one can choose any element. Then,
            !! select only the remaining elements up for choice by choosing only
            !! the larger indices: the smaller, distinct indices are in already, and
            !! since the order doesn't matter, they represent the same combination.
            if ((ind_choice == 1) .or. (i > combination(ind_choice - 1)) ) then
                combination(ind_choice) = i
                call gen_combinations_slow(ind_choice+1,k,n,combination,combinations)
            end if
        end do
    end if


end subroutine