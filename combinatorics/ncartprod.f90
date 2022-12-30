subroutine get_ncartprod(sets_nelements,cartprod)
    !! Gets the Cartesian product of a n sets. This is a wrapper function to
    !! gen_ncartprod, the recursive generator of the Cartesian product. It is here for
    !! convenience, in order to prepare call to the generation subroutine.

    !! Theory
    !! ------
    !! The Cartesian product of two sets, e.g., {x,y} and {1,2,3} can be easily obtained
    !! with a table:
    !!      1       2       3
    !!  x   x,1    x,2     x,3
    !!  y   y,1    y,2     y,3
    !! And it can programatically be generated with two nested do loops:
    !!      do i=1,3
    !!          do j=1,2
    !!              print*, i,j    ! (j=1 means x, j=2 means y)
    !!          end do
    !!      end do
    !! For n sets, however, we would need n nested loops. This is exactly the kind of problem
    !! better solved by recursion. More details in the documentation of gen_cartprod.
    !! In a Cartesian product of n sets, each arrangement is called an n-tuple.

    integer, dimension(:), intent(inout) :: sets_nelements
        !! Array with the number of elements of each set. 
        !! In terms of "do loops", there would be a
        !! "do i=1,[entry of sets_nelements]" for each entry of the array.
    integer, dimension(:,:), allocatable, intent(out) :: cartprod
        !! Cartesian product, which is a 2d allocatable array to output the results.
        !! Shape:  (size(sets_nelements), total_number_of_n_tuples), i.e., each column
        !! contains one arrangement of the Cartesian product.
    integer, dimension(:), allocatable :: ntuple
        !! Work array for the subroutine gen_ncartprod
    integer :: count,nsets

    !! Initial allocation of arrays
    nsets = size(sets_nelements)
    allocate(cartprod(nsets,0))
    allocate(ntuple(nsets))

    call gen_ncartprod(sets_nelements,ntuple,count,cartprod)
end subroutine

recursive subroutine gen_ncartprod(sets_nelements,ntuple,count,cartprod)
    !! Recursive generation of the Cartesian product of n sets

    !! Theory
    !! ------
    !! Recursive functions are functions that call themselves. This
    !! subroutine generates a "do loop" for each set (each entry of
    !! sets_nelements), with the lower boundary being 1 and the upper boundary, the
    !! entry of sets_nelements.

    !! Then the subroutine calls itself during each iteration of the do loop,
    !! but with sets_nelements having an entry fewer than the outer call, until
    !! the entries of sets_nelements have been exhausted. This is equivalent to nsets
    !! nested do loops, one per entry of sets_nelements. During each iteration, of each
    !! outer do loop, an element of the n-tuple is generated. However, this n-tuple is
    !! ready for output only after it has reached the deepest call (empty sets_nelements).
    !! Because of this, the information (partially built ntuple) from the outer call
    !! has to be passed to the deepest call (when we know ntuple is ready).

    !! When sets_nelements is an empty array, then a new column is appended to
    !! cartprod, the ntuple is outputted and it's reused for the next one.

    integer, dimension(:), intent(inout) :: sets_nelements
        !! Number of elements of each set
    integer, dimension(:), intent(inout) :: ntuple
        !! Work array, used to pass information during each recursive call of the subroutine
        !! It should have the same shape as the initial call to sets_nelements, i.e., =nsets
    integer, dimension(:,:), allocatable, intent(inout) :: cartprod
        !! Cartesian product, in the end, shape(nsets,total_number_of_ntuples)
    integer, intent(inout) :: count
        !! Counter of the total number of ntuples.
    integer :: remaining_nsets,i

    remaining_nsets = size(sets_nelements)

    if (remaining_nsets /= 0) then
        !! One do loop per entry of sets_nsets
        do i=1,sets_nelements(remaining_nsets)
            !! Collect the (outer) index for the deepest call
            ntuple(remaining_nsets) = i
            !! Now we dealt with one entry. Call the subroutine again for the remaining
            !! entries until we call it with an empty array.
            call gen_ncartprod(sets_nelements(:remaining_nsets-1),ntuple,count,cartprod)
        end do
    else
        !! If there are no remaining sets, we have reached the
        !! deepest call to the nested loops. The ntuple is finished.
        count = count + 1
        !! Append a column to cartprod with the finished ntuple
        cartprod = reshape([cartprod,ntuple],[size(cartprod,1),count])
    end if

end subroutine