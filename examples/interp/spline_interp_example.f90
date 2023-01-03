program spline_interp_example
    use iso_fortran_env, only: dp => real64
    use shtawa
    implicit none
    integer, parameter :: nx=50
    real(dp), dimension(nx) :: x_grid, f_grid, spline_coefs
    real(dp) :: f,x

    x = 3.2d0

    generate_data: block
        integer :: i
        x_grid = [(i*0.3d0,i=1,nx)]
        f_grid = sin(x_grid)-0.5*x_grid**0.2
    end block generate_data
    
    call build_cubic_spline(x_grid,f_grid,spline_coefs)
    
    print*, "x = ", x
    f = spline_interpolate(x,x_grid,f_grid,spline_coefs)
    print*, "Interpolated value:",f
    print*, "True value:",sin(x)-0.5*x**0.2

    x = 1.5d0
    print*, "x = ", x
    f = spline_interpolate(x,x_grid,f_grid,spline_coefs)
    print*, "Interpolated value:",f
    print*, "True value:",sin(x)-0.5*x**0.2

end program