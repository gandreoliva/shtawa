recursive subroutine get_combinations_slow(n,k,combinations)
    !! Generates the combinations of all the elements of a given set by generating a
    !! Cartesian product and discarding elements (this is why it's slow).

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

    integer, save :: ind_choice = 1
        !! index of the currently examined choice in the set
    integer, intent(in) :: n
        !! Number of elements in the set
    integer, intent(in) :: k
        !! Number of choices
    integer, dimension(:), save, allocatable :: combination
        !! Work array (communicates information from outer calls into deeper calls)
    integer, dimension(:,:), intent(inout), allocatable :: combinations
        !! Array shape(n,number_of_combinations) that contains one combination of the set
        !! in each column (output array).
    integer :: i

    if (.not. allocated(combinations)) allocate(combinations(k,0))
    if (.not. allocated(combination)) allocate(combination(k))

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
                ind_choice = ind_choice + 1
                call get_combinations_slow(n,k,combinations)
            end if
        end do
    end if

    ind_choice = ind_choice - 1

end subroutine