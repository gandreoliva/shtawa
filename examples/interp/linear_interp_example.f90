program linear_interp_example
    use iso_fortran_env, only: dp => real64
    use shtawa
    implicit none
    real(dp), dimension(1:10) :: x_grid, temperature
    integer, dimension(1:2) :: ind_left
    real(dp), dimension(1:2) :: x
    integer :: j

    x_grid = [(0.1d0*j, j=1, 10)]
    temperature = sin(x_grid) + 3d0

    x = [3.63d-1, 0.642d0]

    ! print*, "Grid"
    ! print*, x_grid

    call locate_points_in_grid(x,x_grid,ind_left)
    if(any(ind_left == 0) .or. any(ind_left == size(x_grid))) error stop 

    print*, "locating points in the grid"
    print*, "  ind_left        x_left                   x                   x_right"
    do j=1,size(ind_left)
        print*, ind_left(j), x_grid(ind_left(j)), x(j), x_grid(ind_left(j)+1)
    end do

    print*, "Interpolation at x =", x
    print*, linear_interpolate(x,x_grid(ind_left),x_grid(ind_left+1),temperature(ind_left),temperature(ind_left+1))
    print*, "True value"
    print*, sin(x) + 3d0
    ! print*, temperature

end program