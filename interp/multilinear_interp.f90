recursive subroutine multilinear_interpolate_point(point,grids,f_grid,f)
    !! Multilinear interpolation of one point

    !! Algorithm
    !! ---------
    !! This subroutine applies linear interpolation recursively to each dimension.
    !! For example, for 2 coordinates,
    !!  xL,yR       .   xR,yR
    !!              p               R: right, L: left, p: point
    !!  xL,yL       .   xR,yL
    !! the subroutine applies interpolation in this way:
    !!      f(x,y) = interp(y,yL,yR, f(x,yL), f(x,yR) )
    !!                    ______________|         |________________
    !!         interp(x,xL,xR,f(xL,yL),f(xR,yL))         interp(x,xL,xR,f(xL,yR),f(xR,yR))
    !! The subroutine sees that it can't compute f(x,xL) and f(x,xR), so it calls itself
    !! recursively until only scalars are involved. The subroutine locate_points_in_grid
    !! finds the left indices where the desired value (x,y) is located. The right indices
    !! are found by adding 1 to the left indices. The variable corner is used to build the
    !! evaluation corner of f, for example, for f(xR,yL), corner is (1,0).

    real(wp), dimension(:) :: point
        !! Array with the coordinates of the point to be interpolated: [x1,x2,x3,...]
    type(grid), dimension(:) :: grids
        !! The type(grid) holds one array that contains the grid points. Here, we need
        !! an array of type(grid). For example, for 3 coordinates,
        !! grids(1)%points = x1_grid; grids(2)%points = x2_grid; grids(3)%points = x3_grid
        !! This derived type is used because each grid might have a different number
        !! of cells, and an array of grids allows us to deal with all of them with one var.
    real(wp), dimension(*) :: f_grid
        !! Array rank(ncoords) with the values of the function at each grid point.
        !! For example, for 3 coordinates, f_grid(x1,x2,x3) is one value of the field f.
    real(wp) :: f
        !! Final interpolated value
    integer :: remaining_n, i
    integer, dimension(:), allocatable, save :: ind_left, f_grid_shape, corner
        !! Arrays to recursively build the indices of the left corner and the right corner
    real(wp) :: f_left, f_right, x_left, x_right

    remaining_n = size(point)


    if (.not. allocated(ind_left)) allocate(ind_left(remaining_n))
    if (.not. allocated(corner)) allocate(corner(remaining_n))
    if (.not. allocated(f_grid_shape)) then
        allocate(f_grid_shape(remaining_n))
        f_grid_shape = [(size(grids(i)%points),i=1,remaining_n)]
    end if

    if (remaining_n > 1) then
        !! Fill ind_left backwards (the recursion is called from last to first)
        !! This means that, e.g., in 3D interp., we first interpolate in z, then y, then x.
        call locate_points_in_grid([point(remaining_n)],grids(remaining_n)%points, &
            & ind_left(remaining_n:remaining_n))

        !! Recursive calls to find f_left and f_right for this dimension
        !! Left corner in this dimension (to find f_left)
        corner(remaining_n) = 0
        call multilinear_interpolate_point(point(:remaining_n-1),grids,f_grid,f_left)
        !! Right corner in this dimension (to find f_right)
        corner(remaining_n) = 1
        call multilinear_interpolate_point(point(:remaining_n-1),grids,f_grid,f_right)

        !! Interpolate in this dimension
        f = linear_interpolate(point(remaining_n),&
            &x_left=grids(remaining_n)%points(ind_left(remaining_n)),&
            &x_right=grids(remaining_n)%points(ind_left(remaining_n)+1),&
            &f_left=f_left,f_right=f_right)

    else
        !! Maximum depth: interpolate
        call locate_points_in_grid([point(1)],grids(1)%points,ind_left(1:1))
        !! Check that the desired point is actually inside of the grid
        if (any((ind_left == 0) .or. (ind_left == size(ind_left)))) error stop "Point outside of grid"
               
        x_left = grids(1)%points(ind_left(1))
        x_right = grids(1)%points(ind_left(1)+1)

        corner(1) = 0 
        f_left = narray_get_element(f_grid,f_grid_shape,corner+ind_left)

        corner(1) = 1
        f_right = narray_get_element(f_grid,f_grid_shape,corner+ind_left)

        !! Deepest interpolation (first dimension)
        f = linear_interpolate(point(1),x_left,x_right,f_left,f_right)
    end if
end subroutine


pure function narray_get_element(array,ashape,ind) result(element)
    !! Returns an element of an array of any rank.
    !! Implementation notes: there are 'assumed rank' arrays in Fortran, but they
    !! require a select rank construct that handles every rank individually.
    !! In this implementation, we use the older 'assumed size' with only one dimension.
    !! This has the effect of flattening the array.
    real(wp), dimension(*), intent(in) :: array
        !! Array of any rank
    integer, dimension(:), intent(in) :: ind
        !! Indices of the desired element
    integer, dimension(:),intent(in) :: ashape
        !! Shape of array
    real(wp) :: element
    integer :: loc,i

    !! 'array' is flattened to rank(1), in column order. The following
    !! computes the location of the desired element in array
    loc = ind(1)
    do i=2,size(ind)
        loc = loc + (ind(i) - 1)*product(ashape(1:i-1))
    end do
    element = array(loc)
end function



!!!
!! Assumed-rank version of get element
! pure function narray_get_element(narray,indices) result(element)
!     real(wp), dimension(..) :: narray
!     real(wp), dimension(:) :: ind
!     select rank(narray)
!     rank(1); element = narray(ind(1))
!     rank(2); element = narray(ind(1),ind(2))
!     rank(3); element = narray(ind(1),ind(2),ind(3))
!     rank(4); element = narray(ind(1),ind(2),ind(3),ind(4))
!     rank(5); element = narray(ind(1),ind(2),ind(3),ind(4),ind(5))
!     rank(6); element = narray(ind(1),ind(2),ind(3),ind(4),ind(5),ind(6))
!     rank(7); element = narray(ind(1),ind(2),ind(3),ind(4),ind(5),ind(6),ind(7))
!     end select
! end function