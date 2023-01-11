recursive subroutine get_ncartprod(sets_nelements,cartprod)
	!! Generates the Cartesian product of n sets with different number of elements.
	!!
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
	!! better solved by recursion. In a Cartesian product of n sets, each arrangement
	!! is called an n-tuple.

	!! Algorithm
	!! ------
	!! Recursive functions are functions that call themselves. This
	!! subroutine generates a "do loop" for each set (each entry of
	!! sets_nelements), with the lower boundary being 1 and the upper boundary, the
	!! entry of sets_nelements.

	!! Then the subroutine calls itself during each iteration of the do loop,
	!! but with sets_nelements having an entry fewer than the outer call, until
	!! the entries of sets_nelements have been exhausted. This is equivalent to nsets
	!! nested do loops, one per entry of sets_nelements. During each iteration, of each
	!! outer do loop, an element of the n-tuple is generated.

	!! When sets_nelements is an empty array (deepest call), then a new column is 
	!! appended to cartprod (the finished n-tuple).


	integer, dimension(:), intent(inout) :: sets_nelements
		!! Number of elements of each set
	integer, dimension(:,:), allocatable, intent(inout) :: cartprod
		!! Cartesian product, in the end, shape(nsets,total_number_of_ntuples)
	integer, dimension(:), allocatable, save :: ntuple
		!! Work array, used to pass information during each recursive call of the subroutine.
		!! The "save" attribute is used to keep the array values during calls.
		!! Alternatively, one can pass this array as an argument, or use a global variable,
		!! but then the user has to handle an array intended for internal use.
	integer :: remaining_nsets,i


	remaining_nsets = size(sets_nelements)

	if (.not. allocated(ntuple)) allocate(ntuple(size(sets_nelements)))
	if (.not. allocated(cartprod)) allocate(cartprod(size(sets_nelements),0))

	if (remaining_nsets /= 0) then
		!! One do loop per entry of sets_nsets
		do i=1,sets_nelements(remaining_nsets)
			!! Collect the (outer) index. We fill ntuple backwards because we call
			!! the subroutine recursively backwards (see below)
			ntuple(size(ntuple)-remaining_nsets+1) = i
			!! Now we dealt with one entry (the last one of sets_nelements).
			!! Call the subroutine again for the remaining entries until we call
			!! it with an empty array.
			call get_ncartprod(sets_nelements(:remaining_nsets-1),cartprod)
		end do
	else
		!! If there are no remaining sets, we have reached the
		!! deepest call to the nested loops. The ntuple is finished.
		!! Append a column to cartprod with the finished ntuple
		cartprod = reshape([cartprod,ntuple],[size(cartprod,1),size(cartprod,2)+1])
	end if

end subroutine