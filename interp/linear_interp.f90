elemental function linear_interpolate(x,x_left,x_right,f_left,f_right) result(f)
    !! Theory
    !! -----
    !! Consider a function fg(:) whose values are only known at given points in a grid xg(:).
    !! fg -|             *
    !!    -|     *
    !!    -| *
    !!     +-|---|---x---|-- xg
    !! To find the value of f at an arbitrary point x, we first locate the grid point where the
    !! desired value is contained (in the above example, it is located in the cell 2;
    !! the number is in reference to the left edge of the cell). Then, when we build the equation of 
    !! a straight line with the values of x and f at the left and right edges,
    !! and we finally use the straight line to find the desired value.
    !! Linear interpolation is a very simple form of interpolation and has the significant
    !! weakness that because the straight lines form sharp corners, the resulting interpolating
    !! function is not differentiable.
    real(wp), intent(in) :: x
        !! Desired value of the independent variable (not necessarily in the grid)
    real(wp), intent(in) :: x_left,x_right
        !! values of x from the grid, to the "left" and "right" of x,
        !! i.e., x_left <= x <= x_right
    real(wp), intent(in) :: f_left,f_right
        !! f evaluated at the left and right grid points
        !! (i.e., f_left=f(x_left), f_right=f(x_right))
    real(wp) :: f
        !! The interpolated value, f(x)
    f = (f_right-f_left)/(x_right-x_left)*(x-x_left) + f_left
end function


subroutine locate_points_in_grid(x,x_grid,ind_left)
    !! Locates a list of 1d points x in a grid using the bisection algorithm
    real(wp), dimension(:), intent(in) :: x
        !! List of x values (points) to be located
    real(wp), dimension(:), intent(in) :: x_grid
        !! Grid in the x direction
    integer, dimension(:), intent(out) :: ind_left
        !! Output: list of the indices of the cells where x is located.
        !! The "left" value is provided, i.e., x_grid(ind_left) <= x < x_grid(i_right)
        !! Warnings:
        !! * the indices of ind_left start at 1. If x_grid starts at another lower
        !! boundary, don't forget to do  ind_left = ind_left - 1 +lbound(x_grid,1) 
        !! before evaluating!
        !! * If ind_left = 0, it means that x is outside (to the left) of the grid.
        !! * If ind_left = n, it means that x is outside (to the right) of the grid.
    integer :: curr_i_left, curr_i_right, curr_i_midpoint
    integer :: ngrid, npoints, j
    logical :: is_ascending !! Whether x_grid is in ascending or descending order
    
    ngrid = size(x_grid)
    npoints = size(x)
    !! This handles both cases, of whether the grid is in ascending or descending order
    is_ascending =  ( x_grid(ngrid) >= x_grid(1) )
    
    do concurrent (j=1:npoints) !! for every point in x...
        curr_i_left = 0
        curr_i_right = ngrid + 1
        
        !! Bisection for integers
        do while (curr_i_right - curr_i_left > 1)
            curr_i_midpoint = (curr_i_left + curr_i_right)/2
            if ( (x(j) >= x_grid(curr_i_midpoint)) .eqv. is_ascending ) then
                curr_i_left = curr_i_midpoint
            else
                curr_i_right = curr_i_midpoint
            end if
        end do

        !! Output and handling of end points
        if (x(j) == x_grid(1)) then 
            ind_left(j) = 1
        else if (x(j) == x_grid(ngrid)) then
            ind_left(j) = ngrid - 1
        else
            ind_left(j) = curr_i_left
        end if
    end do

end subroutine